/**
 * \copyright 2014 Danny Kaufman, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifdef WIN32
#define NOMINMAX
#include <windows.h>
#include <concrtrm.h>
#endif
#include "StrandImplicitManager.hh"
#include "MeshScriptingController.hh"
#include "FluidScriptingController.hh"
#include "ImplicitStepper.hh"
#include "StrandDynamicTraits.hh"
#include "DOFScriptingController.hh"
#include "../Collision/CollisionDetector.hh"
#include "../Collision/ElementProxy.hh"
#include "../Collision/EdgeFaceIntersection.hh"
#include "../Collision/ContinuousTimeCollision.hh"
#include "../Collision/EdgeEdgeCollision.hh"
#include "../Collision/VertexFaceCollision.hh"
#include "../Collision/EdgeFaceCollision.hh"
#include "../Collision/CollisionUtils.hh"
#include "../Core/ElasticStrand.hh"
#include "../Render/StrandRenderer.hh"
#include "../Utils/SpatialHashMap.hh"
#include "../Utils/LoggingTimer.hh"
#include "../Utils/Memory.hh"
#include "../Utils/MemUtilities.hh"
#include "../Utils/MathUtilities.hh"

#include "../../bogus/Interfaces/MecheEigenInterface.hpp"
#include "../../bogus/Interfaces/OneCollisionProblem.hpp"

#include <boost/lexical_cast.hpp>

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <mutex>
#include <condition_variable>

#if defined(_OPENMP)
#include <omp.h>
#endif

#ifndef M_SQRT3
#define M_SQRT3 1.73205080756887729352744634151
#endif

namespace strandsim
{

	StrandImplicitManager::StrandImplicitManager(const std::vector<ElasticStrand*>& strands,
		const std::map<std::pair<int, int>, std::set< std::pair<int, int> > >& collision_free,
		const std::vector< std::shared_ptr<MeshScriptingController> >& meshScripting_controllers,
		const std::vector< std::shared_ptr<FluidScriptingController> >& fluidScripting_controllers,
		const std::vector<ConstraintScriptingController*>& constraintScripting_controllers,
		Scalar startTime, Scalar dt, const SimulationParameters& params, SubStepCallback* sub_callback)
		: m_time(startTime)
		, m_dt(dt)
		, m_params(params)
		, m_strands(strands)
		, m_collision_free(collision_free)
		, m_steppers()
		, m_meshScriptingControllers(meshScripting_controllers)
		, m_fluidScriptingControllers(fluidScripting_controllers)
		, m_constraintScriptingControllers(constraintScripting_controllers)
		, m_collisionDetector(NULL)
		, m_collisionDatabase(strands.size())
		, m_hashMap(NULL)
		, m_mem_usage_accu(0)
		, m_mem_usage_divisor(0)
		, m_substep_callback(sub_callback)
	{
		// Start logging
		if (!g_log)
		{
			g_log = new TextLog(std::cout, static_cast<MsgInfo::Severity>(m_params.m_logLevel));
		}
		DebugStream(g_log, "") << "Creating new Strand Manager with dt = " << m_dt;

		// Enforce desired or maximum number of threads
		{
#ifdef WIN32
			SYSTEM_INFO sysinfo;
			GetSystemInfo(&sysinfo);
#endif
			const int numThreads =
				m_params.m_numberOfThreads > 0 ? m_params.m_numberOfThreads :
#ifdef WIN32
				sysinfo.dwNumberOfProcessors;
#else
				sysconf(_SC_NPROCESSORS_ONLN);
#endif
#if defined(_OPENMP)
			omp_set_num_threads(numThreads);
#endif

		}

		for (std::vector<ElasticStrand*>::const_iterator strand = m_strands.begin();
			strand != m_strands.end(); ++strand)
		{
			(*strand)->setParent(this);
		}

		// Store the stepper for each strand
		// Also store edge proxies for collision detection
		for (int i = 0; i < m_strands.size(); ++i)
		{
			ElasticStrand* strand = m_strands[i];
			strand->clearSharedRuntimeForces();

			ImplicitStepper* stepper = new ImplicitStepper(*strand, m_params);
			m_steppers.push_back(stepper);
			strand->setStepper(stepper);

			for (int vtx = 0; vtx < strand->getNumEdges(); ++vtx)
			{
				const CollisionParameters& cp = strand->collisionParameters();
				m_elementProxies.push_back(
					new CylinderProxy(*strand, vtx, cp.externalCollisionsRadius(vtx)));
			}
		}

		for (auto controller = m_meshScriptingControllers.begin();
			controller != m_meshScriptingControllers.end(); ++controller)
		{
			auto mesh = (*controller)->getCurrentMesh();

			for (unsigned f = 0; f < mesh->nf(); ++f)
			{
				m_elementProxies.push_back(new FaceProxy(f, *controller));
			}
		}

		for (auto controller = m_fluidScriptingControllers.begin();
			controller != m_fluidScriptingControllers.end(); ++controller)
		{
			(*controller)->initialize();
		}

		// Now that we have all the external forces in place, do reverse hairdo if needed

		m_externalContacts.resize(m_strands.size());

		m_collisionDetector = new CollisionDetector(m_elementProxies);
	}

	StrandImplicitManager::~StrandImplicitManager()
	{
		for (auto stepper = m_steppers.begin(); stepper != m_steppers.end(); ++stepper)
		{
			delete* stepper;
		}

		for (auto elem = m_elementProxies.begin(); elem != m_elementProxies.end(); ++elem)
		{
			delete* elem;
		}

		delete m_hashMap;
		m_hashMap = NULL;

		delete m_collisionDetector;
	}

	void StrandImplicitManager::execute(int total_num_substeps, int total_substep_id, const Scalar total_substep_dt)
	{
		//    InfoStream( g_log, "" ) << "Simulation time: " << getTime();
		m_dt = total_substep_dt;

		m_timings.push_back(StepTimings());

		// Clean up debug drawing data
		for (auto strand = m_strands.begin(); strand != m_strands.end(); ++strand)
			(*strand)->dynamics().clearDebugDrawing();

		//step(total_num_substeps, total_substep_id);
		step();
		if (m_substep_callback) m_substep_callback->executeCallback();	// update drawing data

		/*m_statTotalCollisions += m_statExternalContacts + m_statMutualCollisions;
		InfoStream(g_log, "Contact")
			<< "external: " << m_statExternalContacts
			<< " mutual: " << m_statMutualCollisions
			<< " total: " << m_statExternalContacts + m_statMutualCollisions
			<< " overall: " << m_statTotalCollisions;*/

		InfoStream(g_log, "Timing") << "Last frame (ms):";
		print(m_timings.back());
		InfoStream(g_log, "Timing") << "Average after " << m_timings.size() << " frames (ms):";
		print(m_timings);

		if (m_params.m_statGathering)
		{
			//printNewtonSolverBreakdownTiming<InfoStream>();
			print<InfoStream>(m_cdTimings);
			print<InfoStream>(m_solverStat);
		}

		m_cdTimings.reset();
		m_solverStat.reset();
	}


	Scalar StrandImplicitManager::getCFL() const
	{
		if (!m_fluidScriptingControllers.size() || !m_fluidScriptingControllers[0])
			return 1e+20;

		return m_fluidScriptingControllers[0]->cfl();
	}

	//void StrandImplicitManager::step(int total_num_substeps, int total_substep_id)
	//{
	//	MemoryDiff<DebugStream> mem("step");
	//	SubStepTimings timings;

	//	std::cout << "[Prepare Simulation Step]" << std::endl;
	//	Timer timer("step", false);
	//	step_prepare(m_dt);
	//	timings.prepare = timer.elapsed();

	//	m_collidingGroups.clear();
	//	m_collidingGroupsIdx.assign(m_strands.size(), -1);

	//	if (m_params.m_solveCollision && !m_params.m_useCTRodRodCollisions) {
	//		std::cout << "[Setup Hair-Hair Collisions]" << std::endl;
	//		timer.restart();
	//		setupHairHairCollisions(m_dt);
	//		timings.hairHairCollisions = timer.elapsed();
	//	}

	//	std::cout << "[Step Dynamics]" << std::endl;
	//	timer.restart();
	//	switch (m_params.m_solverType)
	//	{
	//	case 0:
	//		step_dynamics(total_num_substeps, total_substep_id, m_dt); break;
	//	case 1:
	//		step_dynamics_Jacobi(m_dt); break;
	//	case 2:
	//		step_dynamics_Gauss(m_dt); break;
	//	case 3:
	//		step_dynamics_CG(m_dt); break;
	//	}
	//	timings.dynamics = timer.elapsed();

	//	if (m_params.m_solveCollision) {
	//		std::cout << "[Setup Continous-Time Collisions]" << std::endl;
	//		timer.restart();
	//		setupMeshHairCollisions(m_dt);
	//		timings.meshHairCollisions = timer.elapsed();

	//		std::cout << "[Process Collisions]" << std::endl;
	//		timer.restart();
	//		step_processCollisions(m_dt);
	//		timings.processCollisions = timer.elapsed();

	//		std::cout << "[Solve Collisions]" << std::endl;
	//		timer.restart();
	//		step_solveCollisions(total_num_substeps, total_substep_id);
	//		timings.solve += timer.elapsed();
	//	}

	//	assert(m_collisionDetector->empty());

	//	m_mutualContacts.clear();
	//	m_elasticMutualContacts.clear();

	//	m_time += m_dt;

	//	print<CopiousStream>(timings);

	//	InfoStream(g_log, "") << "Database size: " << m_collisionDatabase.computeSizeInBytes();
	//	m_timings.back().push_back(timings);

	//	printMemStats();
	//}

	void StrandImplicitManager::printMemStats()
	{
		int peak_idx = 0;
		int cur_idx = 0;

		const char* mem_units[] = {
			"B", "KB", "MB", "GB", "TB", "PB"
		};

		size_t cur_usage = memutils::getCurrentRSS();
		m_mem_usage_accu += (Scalar)cur_usage;
		m_mem_usage_divisor += 1.0;

		size_t peak_usage = memutils::getPeakRSS();
		Scalar peak_mem = (Scalar)peak_usage;
		while (peak_mem > 1024.0 && peak_idx < (int)(sizeof(mem_units) / sizeof(char*))) { peak_mem /= 1024.0; peak_idx++; }

		Scalar avg_mem = m_mem_usage_accu / m_mem_usage_divisor;
		while (avg_mem > 1024.0 && cur_idx < (int)(sizeof(mem_units) / sizeof(char*))) { avg_mem /= 1024.0; cur_idx++; }

		std::cout << "Peak Mem Usage, " << peak_mem << mem_units[peak_idx] << ", Avg Mem Usage, " << avg_mem << mem_units[cur_idx] << std::endl;
	}

	//void StrandImplicitManager::setupMeshHairCollisions(Scalar dt)
	//{
	//	if (m_params.m_skipRodMeshCollisions || !m_strands.size())
	//		return;

	//	Timer tt("MeshHair", false, "controllers");
	//	DebugStream(g_log, "") << "StrandImplicitManager::setup continuous collision detection";

	//	tt.restart("buildBVH");
	//	m_collisionDetector->buildBVH(false);
	//	m_cdTimings.buildBVH = tt.elapsed();

	//	tt.restart("findCollisions");
	//	EdgeFaceIntersection::s_doProximityDetection = true;
	//	// TODO:
	//	// split hair/hair and hair/mesh CCD - build different BVH (avoid compare between them)

	//	// ignoreCTRodRod = false : CCD of edges among hairs
	//	//     -> add EdgeEdgeCollision to m_continuousTimeCollisions
	//	// ignoreContinuousTime = false : CCD of vertex-face and edge-face (CCD of edge and 3 edges of the face)
	//	//     -> add VertexFaceCollision and EdgeFaceCollision to m_continuousTimeCollisions
	//	// ignoreProximity = false : edge-face intersection test (s_doProximityDetection = true) using position before unconstraint update
	//	//     -> add EdgeFaceIntersection to m_proximityCollisions
	//	m_collisionDetector->findCollisions(!m_params.m_useCTRodRodCollisions, false, false); // set whether cd should do ctc for hairs, etc
	//	m_cdTimings.findCollisionsBVH = tt.elapsed();
	//	DebugStream(g_log, "") << "CollisionDetector Stat: " << *m_collisionDetector;

	//	// We do CT collisions first in order to guess meshes normals signs as soon as possible
	//	//tt.restart( "continuous" );
	//	tt.restart("narrowPhase");
	//	// compact m_continuousTimeCollisions in ProximityCollision, and add these collisions to m_externalContacts
	//	doContinuousTimeDetection(dt);
	//	//tt.restart( "proximity" );
	//	// compact m_proximityCollisions in ProximityCollision, and add these collisions to m_externalContacts
	//	doProximityMeshHairDetection(dt);
	//	m_cdTimings.narrowPhase = tt.elapsed();

	//	tt.restart("buildBVH");
	//	m_collisionDetector->clear();
	//	m_cdTimings.buildBVH += tt.elapsed();
	//}

	void StrandImplicitManager::doProximityMeshHairDetection(Scalar dt)
	{
		std::list<CollisionBase*>& collisionsList = m_collisionDetector->getProximityCollisions();

		TraceStream(g_log, "") << "doProximityMeshHairDetection: before pruning: " << collisionsList.size();

		unsigned nInt = 0;
		for (auto intersection = collisionsList.begin(); intersection != collisionsList.end(); ++intersection)
		{
			if (EdgeFaceIntersection * efi = dynamic_cast<EdgeFaceIntersection*>(*intersection))
			{
				EdgeProxy* edge = efi->getEdge();

				const unsigned strIdx = edge->getStrand().getGlobalIndex();
				const unsigned edgeIdx = edge->getVertexIndex();

				ProximityCollision edgeFaceCollision;

				if ((m_strands[strIdx]->isVertexFreezed(edgeIdx) && m_strands[strIdx]->isVertexFreezed(edgeIdx + 1)) || (m_strands[strIdx]->isVertexGoaled(edgeIdx) && m_strands[strIdx]->isVertexGoaled(edgeIdx + 1)))
					continue;

				edgeFaceCollision.normal = efi->getFaceNormal();
				edgeFaceCollision.mu = sqrt(
					efi->faceFrictionCoefficient()
					* edge->getStrand().collisionParameters().frictionCoefficient(edgeIdx));

				edgeFaceCollision.objects.second.globalIndex = -1;
				edgeFaceCollision.objects.second.vertex = efi->getFaceId();
				edgeFaceCollision.objects.second.freeVel = efi->getFaceVelocity(dt); //+ ContinuousTimeCollision::s_extraRadius * collision.normal ;

				if (addExternalContact(strIdx, edgeIdx, efi->getEdgeAbscissa(), edgeFaceCollision)) 
				{
					++nInt;
				}
			}
		}

		DebugStream(g_log, "") << "We found " << nInt << " edge/face intersection";
	}

	void StrandImplicitManager::doContinuousTimeDetection(Scalar dt)
	{
		EdgeFaceIntersection::s_doProximityDetection = false;

		std::list<CollisionBase*>& collisionsList = m_collisionDetector->getContinuousTimeCollisions();

		TraceStream(g_log, "") << "doContinuousTimeDetection: before pruning: " << collisionsList.size() << " hits";

		if (collisionsList.empty())
		{
			TraceStream(g_log, "") << "No more collisions for this time step";
			return;
		}

		int nExt = 0;

		for (auto collIt = collisionsList.begin(); collIt != collisionsList.end(); ++collIt)
		{
			ContinuousTimeCollision* const ctCollision = dynamic_cast<ContinuousTimeCollision*>(*collIt);
			if (!ctCollision)
				continue;

			ProximityCollision collision;

			collision.m_originalCTCollision = ctCollision;
			collision.normal = ctCollision->normal();
			collision.distance = ctCollision->distance();

			EdgeEdgeCollision* eeCollision = dynamic_cast<EdgeEdgeCollision*>(ctCollision);

			if (eeCollision)
			{
				ElasticStrand* sP = eeCollision->getFirstStrand();
				ElasticStrand* sQ = eeCollision->getSecondStrand();
				int iP = eeCollision->getFirstVertex();
				int iQ = eeCollision->getSecondVertex();
				const CollisionParameters& cpP = sP->collisionParameters();
				const CollisionParameters& cpQ = sQ->collisionParameters();

				auto itr_cf = m_collision_free.find(std::pair<int, int>(sP->getGlobalIndex(), iP));
				if (itr_cf != m_collision_free.end()) {
					const std::set<std::pair<int, int> >& free_set = itr_cf->second;
					if (free_set.find(std::pair<int, int>(sQ->getGlobalIndex(), iQ)) != free_set.end()) {
						continue;
					}
				}

				bool acceptFirst = cpP.reactsToSelfCollisions() && (!(sP->isVertexFreezed(iP) || sP->isVertexGoaled(iP)) || (iP < sP->getNumVertices() - 2 && !(sP->isVertexFreezed(iP + 1) || sP->isVertexGoaled(iP + 1))));
				bool acceptSecond = cpQ.reactsToSelfCollisions() && (!(sQ->isVertexFreezed(iQ) || sQ->isVertexGoaled(iQ)) || (iQ < sP->getNumVertices() - 2 && !(sQ->isVertexFreezed(iQ + 1) || sQ->isVertexGoaled(iQ + 1))));
				
				// && -> Discard collisions on root immunity length
				// || -> make them as external objects
				if (!acceptFirst && !acceptSecond)
					continue;
				if (cpP.usesFakeLayering() && cpQ.usesFakeLayering())
				{
					if (!(acceptFirst && acceptSecond))
						continue;
				}

				collision.mu = sqrt(cpP.frictionCoefficient(iP) * cpQ.frictionCoefficient(iQ));

				collision.objects.first.globalIndex = sP->getGlobalIndex();
				collision.objects.first.vertex = iP;
				collision.objects.first.abscissa = eeCollision->getFirstAbscissa();
				collision.objects.first.freeVel.setZero();

				collision.objects.second.globalIndex = sQ->getGlobalIndex();
				collision.objects.second.vertex = iQ;
				collision.objects.second.abscissa = eeCollision->getSecondAbscissa();
				collision.objects.second.freeVel.setZero();

				// if accept both: add to m_mutualContacts
				// if only accept one: add to m_externalContacts

				if (acceptFirst && acceptSecond) {
					collision.swapIfNecessary();
#pragma omp critical(mutualContact)
					m_mutualContacts.push_back(collision);
				}
				else {
					makeExternalContact(collision, acceptFirst);
				}

			}
			else {

				ElasticStrand* const strand = ctCollision->getFirstStrand();
				const unsigned edgeIdx = ctCollision->getFirstVertex();

				if (FaceCollision* fCollision = dynamic_cast<FaceCollision*>(ctCollision)) // Assignment =
				{
					collision.mu = sqrt(
						fCollision->faceFrictionCoefficient()
						* strand->collisionParameters().frictionCoefficient(edgeIdx));
				}

				collision.objects.second.globalIndex = -1;

				const unsigned strIdx = strand->getGlobalIndex();
				const Vec3x offset = ctCollision->offset() / dt;

				if (VertexFaceCollision* vfCollision = dynamic_cast<VertexFaceCollision*>(ctCollision)) // Assignment =
				{
					if (strand->isVertexFreezed(edgeIdx) || strand->isVertexGoaled(edgeIdx))
						continue;

					collision.objects.second.vertex = vfCollision->face()->uniqueId();
					collision.objects.second.freeVel = vfCollision->meshVelocity(dt) + offset;

					if (addExternalContact(strIdx, edgeIdx, 0, collision))
					{
						++nExt;
					}
				}
				else if (EdgeFaceCollision* efCollision = dynamic_cast<EdgeFaceCollision*>(ctCollision)) // Assignment =
				{
					if ((strand->isVertexFreezed(edgeIdx) && strand->isVertexFreezed(edgeIdx + 1)) || (strand->isVertexGoaled(edgeIdx) && strand->isVertexGoaled(edgeIdx + 1)))
						continue;

					collision.objects.second.vertex = efCollision->faceEdgeId();
					collision.objects.second.freeVel = efCollision->meshVelocity(dt) + offset;

					// add collision to m_externalContacts, set object.first(freeVel = 0)
					if (addExternalContact(strIdx, edgeIdx, efCollision->abscissa(), collision))
					{
						++nExt;
					}
				}
			}
		}

		DebugStream(g_log, "") << "We found " << m_mutualContacts.size() << " hair/hair contacts, "
			<< nExt << " hair/external contacts";
	}

	bool StrandImplicitManager::addExternalContact(const unsigned strIdx, const unsigned edgeIdx,
		const Scalar abscissa, const ProximityCollision& externalContact)
	{
		auto acceptCollision = [&](int sIdxP, int iP) {
			auto itr_cf = m_collision_free.find(std::pair<int, int>(sIdxP, iP));

			if (itr_cf != m_collision_free.end()) {
				const std::set<std::pair<int, int> >& free_set = itr_cf->second;
				if (free_set.find(std::pair<int, int>(-1, -1)) != free_set.end()) {
					return false;
				}
			}

			return true;
		};

		if (!acceptCollision(strIdx, edgeIdx) || !acceptCollision(strIdx, edgeIdx + 1))
			return false;

#pragma omp critical(externalContact)
		m_externalContacts[strIdx].push_back(externalContact);

		ProximityCollision& collision = m_externalContacts[strIdx].back();

		collision.objects.first.globalIndex = strIdx;
		collision.objects.first.vertex = edgeIdx;
		collision.objects.first.abscissa = abscissa;
		collision.objects.first.freeVel.setZero();

		return true;
	}

	void StrandImplicitManager::setupHairHairCollisions(Scalar dt)
	{
		// hair proximity test based on spacial hash map
		// movable hair pair stored in m_mutualContacts
		// movable hair and freezed hair pair stored in m_externalContacts

		// do hair-hair setup basics and proceed to proximity collisions if requested
		Timer tt("HashMap", false, "update");
		//LoggingTimer<InfoStream> tt( "HashMap", "update" );

		if (m_params.m_skipRodRodCollisions)
			return;

		if (!m_hashMap)
		{
			m_hashMap = new SpatialHashMapT();
		}

		// Compute maximum number of vertices and min colliding radius
		unsigned maxNumVert = 0;
		Scalar meanRadius = 0.;
		Scalar maxEdgeLen = 0.;
		for (int i = 0; i < m_strands.size(); ++i)
		{
			const unsigned nv = m_strands[i]->getNumVertices();
			if (nv > maxNumVert)
			{
				maxNumVert = nv;
			}

			const Scalar rootRad = m_strands[i]->collisionParameters().selfCollisionsRadius(0);
			meanRadius += rootRad;
			maxEdgeLen = std::max(maxEdgeLen, m_strands[i]->getTotalRestLength() / nv);
		}
		m_hashMap->setCellSize(std::max(maxEdgeLen, meanRadius / m_strands.size()));

		m_hashMap->batchUpdate(m_strands, maxNumVert);

		m_cdTimings.updateHashMap = tt.elapsed();

		tt.restart("process");
		//tt.restart( "compute" );

		// compute will actually do nothing, except initializing result
		// The collisions will actually be retrived by batches of 10 objects at each call to result.next
		SpatialHashMapT::Result result(true, 10);
		m_hashMap->compute(result);

		unsigned nRough = 0, nExt = 0;

		//tt.restart( "analyze" );

		// Processes batches of collision in parallel
		SpatialHashMapT::Result::Collisions collisions;

#pragma omp parallel private( collisions ) reduction( + : nRough, nExt )
		while (result.next(collisions))
		{

			const unsigned nColObj = collisions.size();

			std::vector<SpatialHashMapT::Result::Collisions::const_iterator> iters(nColObj);

			int i = 0;
			for (SpatialHashMapT::Result::Collisions::const_iterator first = collisions.begin();
				first != collisions.end(); ++first)
			{
				iters[i++] = first;
			}

			for (int itIdx = 0; itIdx < nColObj; ++itIdx)
			{
				SpatialHashMapT::Result::Collisions::const_iterator& first = iters[itIdx];
				ElasticStrand* sP = first->first;
				const CollisionParameters& cpP = sP->collisionParameters();

				if (sP->getParameters().collisionFree())
					continue;

				if (!cpP.createsSelfCollisions())
					continue;

				const int sIdxP = sP->getGlobalIndex();

				for (auto second = first->second.begin(); second != first->second.end(); ++second)
				{
					ElasticStrand* sQ = second->first;

					const CollisionParameters& cpQ = sQ->collisionParameters();
					if (!cpQ.createsSelfCollisions())
						continue;

					if (sQ->getParameters().collisionFree())
						continue;

					const int sIdxQ = sQ->getGlobalIndex();

					for (auto collision = second->second.begin(); collision != second->second.end();
						++collision)
					{
						int iP = collision->first;
						int iQ = collision->second;

						if ((sP == sQ && std::abs(iP - iQ) < 4) || !(iP && iQ))
							continue;

						auto itr_cf = m_collision_free.find(std::pair<int, int>(sIdxP, iP));
						if (itr_cf != m_collision_free.end()) {
							const std::set<std::pair<int, int> >& free_set = itr_cf->second;
							if (free_set.find(std::pair<int, int>(sIdxQ, iQ)) != free_set.end()) {
								continue;
							}
						}

						++nRough;

						Vec3x normal;
						Scalar s, t, d;
						Scalar rel_vel(0.);
						
						if (!analyseRoughRodRodCollision(sP, sQ, iP, iQ, normal, s, t, d, rel_vel))
						{
							continue;
						}

						bool acceptFirst = cpP.reactsToSelfCollisions() && (!(sP->isVertexFreezed(iP) || sP->isVertexGoaled(iP)) || (iP < sP->getNumVertices() - 2 && !(sP->isVertexFreezed(iP + 1) || sP->isVertexGoaled(iP + 1))));
						bool acceptSecond = cpQ.reactsToSelfCollisions() && (!(sQ->isVertexFreezed(iQ) || sQ->isVertexGoaled(iQ)) || (iQ < sP->getNumVertices() - 2 && !(sQ->isVertexFreezed(iQ + 1) || sQ->isVertexGoaled(iQ + 1))));

						// && -> Discard collisions on root immunity length
						// || -> make them as external objects
						if (!(acceptFirst || acceptSecond))
							continue;
						if (cpP.usesFakeLayering() && cpQ.usesFakeLayering())
						{
							if (!(acceptFirst && acceptSecond))
								continue;
						}

						ProximityCollision mutualContact;

						mutualContact.normal = normal;
						mutualContact.distance = d;
						mutualContact.relative_vel = rel_vel;
						mutualContact.mu = sqrt(cpP.frictionCoefficient(iP) * cpQ.frictionCoefficient(iQ));

						mutualContact.objects.first.globalIndex = sP->getGlobalIndex();
						mutualContact.objects.first.vertex = iP;
						mutualContact.objects.first.abscissa = s;
						mutualContact.objects.first.freeVel.setZero();

						mutualContact.objects.second.globalIndex = sQ->getGlobalIndex();
						mutualContact.objects.second.vertex = iQ;
						mutualContact.objects.second.abscissa = t;
						mutualContact.objects.second.freeVel.setZero();

						// if accept both and approching: add to m_mutualContacts
						// if only accept one and approching: add to m_externalContacts

						if (acceptFirst && acceptSecond) {
							mutualContact.swapIfNecessary();
#pragma omp critical(mutualContact)
							m_mutualContacts.push_back(mutualContact);
						}
						else {
							makeExternalContact(mutualContact, acceptFirst);
							++nExt;
						}
					}
				}
			}
		}


		DebugStream(g_log, "") << "Rough possible rod/rod collisions: " << nRough << ", real "
			<< (nExt + m_mutualContacts.size());

		if (m_params.m_simulationManager_limitedMemory)
		{
			delete m_hashMap;
			m_hashMap = NULL;
		}
		else
		{
			m_hashMap->clear();
		}

		m_cdTimings.processHashMap = tt.elapsed();
	}

	void StrandImplicitManager::makeExternalContact(ProximityCollision& externalContact, bool onFirstObject)
	{
		int extObjId;
		if (onFirstObject)
		{
			extObjId = externalContact.objects.second.globalIndex;
			externalContact.objects.second.globalIndex = -1;
		}
		else
		{
			extObjId = externalContact.objects.first.globalIndex;
			externalContact.objects.first.globalIndex = -1;
		}

		// This will put the external object on c.objects.second
		externalContact.swapIfNecessary();

		externalContact.objects.second.globalIndex = extObjId;
		const VecXx& velocities = m_steppers[extObjId]->velocities();
		const int iQ = externalContact.objects.second.vertex;
		const Scalar t = externalContact.objects.second.abscissa;

		externalContact.objects.second.freeVel = (1 - t) * velocities.segment<3>(4 * iQ)
			+ t * velocities.segment<3>(4 * iQ + 4);

#pragma omp critical(externalContact)
		m_externalContacts[externalContact.objects.first.globalIndex].push_back(externalContact);
	}

	//void StrandImplicitManager::computeDeformationGradient(ProximityCollision::Object& object) const
	//{
	//	m_strands[object.globalIndex]->getFutureState().computeDeformationGradient(object.vertex,
	//		object.abscissa,
	//		m_steppers[object.globalIndex]->velocities(),
	//		object.defGrad);
	//}

	void StrandImplicitManager::computeDeformationGradient(ProximityCollision::Object& object) const
	{
		m_strands[object.globalIndex]->getFutureState().computeDeformationGradient(
			object.vertex, object.abscissa, object.defGrad
		);
	}

	void StrandImplicitManager::setupDeformationBasis(ProximityCollision& collision) const
	{
		const ProximityCollision* oldCollision = m_collisionDatabase.find(collision);
		if (oldCollision)
		{
			collision.force = oldCollision->force; //Warm start
			collision.updateTransformationMatrix(oldCollision->transformationMatrix);
		}
		else
		{
			collision.generateTransformationMatrix();
		}

		computeDeformationGradient(collision.objects.first);
		if (collision.objects.second.globalIndex != -1)
		{
			computeDeformationGradient(collision.objects.second);
		}
	}
	
	void StrandImplicitManager::step_prepare(Scalar dt)
	{
		// Execute the controllers to script meshes to the future time of this step
#pragma omp parallel for
		for (int i = 0; i < m_meshScriptingControllers.size(); ++i)
		{
			DebugStream(g_log, "") << "Executing controller number " << i;
			m_meshScriptingControllers[i]->setTime(m_time + dt);
			m_meshScriptingControllers[i]->execute(true); // Advance the mesh; compute level set if we do rod/mesh penalties
		}
	}
	/*
	void StrandImplicitManager::redo_step_dynamics(int total_num_substeps, int total_substep_id, Scalar substepDt)
	{
		DebugStream(g_log, "") << "Dynamics";

		std::vector< Scalar > residuals(m_strands.size(), 1e+99);

#pragma omp parallel for
		for (int i = 0; i < m_strands.size(); i++)
		{
			m_steppers[i]->rewind();
			m_steppers[i]->stepFlowData();
			m_steppers[i]->updateAdditionalInertia();
		}

		int hair_subSteps = std::max(1, (int)ceil((Scalar)m_params.m_rodSubSteps / (Scalar)total_num_substeps));

		for (int k = 0; k < hair_subSteps; ++k) {
			std::cout << "start rod substep " << k << std::endl;

			std::vector< unsigned char > passed(m_strands.size(), false);

			while (true) {
#pragma omp parallel for
				for (int i = 0; i < m_strands.size(); i++) {
					if (!passed[i]) {
						m_steppers[i]->startSubstep(k, substepDt / hair_subSteps, 1.0 / hair_subSteps);
						m_steppers[i]->prepareDynamics();
						m_steppers[i]->prepareSolveNonlinear();
					}
				}

				std::vector< int > all_done(m_strands.size(), false);

				int iter = 0;
				while (true) {
#pragma omp parallel for
					for (int i = 0; i < m_strands.size(); i++)
					{
						if (!passed[i] && !all_done[i]) {
							m_steppers[i]->prepareNewtonIteration();
							all_done[i] = m_steppers[i]->performNewtonIteration();
						}
					}

					bool done = true;
					for (int i = 0; i < m_strands.size(); i++)
						if (!passed[i])
							done = done && all_done[i];

					// std::cout << "step dynamics iter: " << iter << std::endl;
					++iter;

					if (done) break;
				}

				std::cout << "redo step dynamics total iter: " << iter << std::endl;

#pragma omp parallel for
				for (int i = 0; i < m_strands.size(); i++)
				{
					if (!passed[i])
						passed[i] = m_steppers[i]->postSolveNonlinear();
				}

				// check if all passed
				bool all_passed = true;

				for (int i = 0; i < (int)m_strands.size(); i++) {
					all_passed = all_passed && passed[i];

					if (!passed[i]) {
						std::cout << "Strand " << i << " failed! Redo with " << m_steppers[i]->getStretchMultiplier() << std::endl;
					}
				}

				if (all_passed)
					break;
			}

			if (k == hair_subSteps - 1) {
				// recover at the last substep
#pragma omp parallel for
				for (int i = 0; i < (int)m_strands.size(); i++)
				{
					m_steppers[i]->update();
					m_steppers[i]->recoverFromBeginning();
					residuals[i] = m_steppers[i]->getNewtonResidual();
				}
			}
			else {
#pragma omp parallel for
				for (int i = 0; i < (int)m_strands.size(); i++)
				{
					m_steppers[i]->update();
					m_steppers[i]->updateRodAccelerationAndVelocity();
					residuals[i] = m_steppers[i]->getNewtonResidual();
				}
			}

			if (m_strands.size())
				residualStats("redo_step_dynamics", residuals);
		}

		// recover stepper Dt
#pragma omp parallel for
		for (int i = 0; i < (int)m_strands.size(); i++)
		{
			m_steppers[i]->setDt(substepDt);
			m_steppers[i]->setFraction(1.0);
			m_steppers[i]->updateRodAcceleration();
		}
	}


	void StrandImplicitManager::step_dynamics_CG(Scalar dt) 
	{
		std::vector<Scalar> residuals(m_strands.size(), 1e+99);
		std::vector<int> iters(m_strands.size());

		//#pragma omp parallel for
		for (int s = 0; s < (int)m_strands.size(); ++s) {
			StrandDynamicTraits& dynamics = m_strands[s]->dynamics();
			VecXx& vel = m_steppers[s]->velocities();
			VecXx& newVel = m_steppers[s]->newVelocities();

			// initialize
			m_strands[s]->setFutureDegreesOfFreedom(m_strands[s]->getCurrentDegreesOfFreedom());
			dynamics.computeViscousForceCoefficients(m_dt);
			dynamics.computeDOFMasses();
			dynamics.getDisplacements().setZero();
			dynamics.getAccelerations().setZero();

			// compute RHS
			VecXx b = vel;
			dynamics.multiplyByMassMatrix(b);
			dynamics.computeFutureForces(!m_params.m_usePreFilterGeometry, false, false, false, std::cout);
			b += m_strands[s]->getFutureTotalForces() * m_dt;

			// compute LHS
			dynamics.computeFutureJacobian(!m_params.m_usePreFilterGeometry, false, false, false, std::cout);
			JacobianMatrixType& A = m_strands[s]->getTotalJacobian();
			A *= m_dt * m_dt;
			dynamics.addMassMatrixTo(A);

			A.multiply(b, 1., vel);
			dynamics.getScriptingController()->fixLHSAndRHS(A, b, m_dt);

			SparseMatx sparseA(vel.size(), vel.size());
			A.convertToSparseMat(sparseA);

			Eigen::ConjugateGradient< SparseMatx > iterativeSolver;
			iterativeSolver.compute(sparseA);
			iterativeSolver.setTolerance(m_params.m_gaussSeidelTolerance);
			iterativeSolver.setMaxIterations(m_params.m_maxNewtonIterations);
			newVel = iterativeSolver.solve(b);
			std::cout << "[cg total iter: " << iterativeSolver.iterations() << ", res: " << iterativeSolver.error() << "]" << std::endl;
			
			dynamics.getAccelerations() = (newVel - vel) / m_dt;
			residuals[s] = iterativeSolver.error() * b.norm();

			// update x_t+1
			VecXx displacements = newVel * m_dt;
			dynamics.getScriptingController()->enforceDisplacements(displacements);
			m_strands[s]->setCurrentDegreesOfFreedom(m_strands[s]->getCurrentDegreesOfFreedom() + displacements);
			m_strands[s]->getFutureState().freeCachedQuantities();
			vel = newVel;

			//residuals[s] = iterativeSolver.error() * b.norm();
		}

		// output residual
		residualStats("step_dynamics_CG", residuals);
	}

	void StrandImplicitManager::step_dynamics(int total_num_substeps, int total_substep_id, Scalar substepDt)
	{
		DebugStream(g_log, "") << "Dynamics";

		std::vector< Scalar > residuals(m_strands.size(), 1e+99);

		// Dynamics system assembly
#pragma omp parallel for
		for (int i = 0; i < (int)m_strands.size(); i++)
		{
			m_steppers[i]->initStepping(substepDt);
		}

		int hair_subSteps = std::max(1, (int)ceil((Scalar)m_params.m_rodSubSteps / (Scalar)total_num_substeps));

		for (int k = 0; k < hair_subSteps; ++k) {
			std::cout << "start rod substep " << k << " with " << (substepDt / hair_subSteps) << std::endl;

			std::vector< unsigned char > passed(m_strands.size(), false);
			while (true) {
#pragma omp parallel for
				for (int i = 0; i < (int)m_strands.size(); i++)
				{
					if (!passed[i]) {
						m_steppers[i]->startSubstep(k, substepDt / hair_subSteps, 1.0 / hair_subSteps);
						m_steppers[i]->prepareDynamics();
						m_steppers[i]->prepareSolveNonlinear();
					}
				}

				std::vector< int > all_done(m_strands.size(), false);

				int iter = 0;
				while (true) {
//#pragma omp parallel for
					for (int i = 0; i < (int)m_strands.size(); i++)
					{
						if (!passed[i] && !all_done[i]) {
							m_steppers[i]->prepareNewtonIteration();
							all_done[i] = m_steppers[i]->performNewtonIteration();
						}
					}

					bool done = true;
					for (int i = 0; i < (int)m_strands.size(); i++) {
						if (!passed[i])
							done = done && all_done[i];
					}

					// std::cout << "step dynamics iter: " << iter << std::endl;
					++iter;

					if (done) break;
				}

				std::cout << "step dynamics total iter: " << iter << std::endl;

#pragma omp parallel for
				for (int i = 0; i < (int)m_strands.size(); i++)
				{
					if (!passed[i])
						passed[i] = m_steppers[i]->postSolveNonlinear();	// fail-safe
				}

				// check if all passed
				bool all_passed = true;

				for (int i = 0; i < (int)m_strands.size(); i++) {
					all_passed = all_passed && passed[i];

					if (!passed[i]) {
						std::cout << "Strand " << i << " failed! Redo with " << m_steppers[i]->getStretchMultiplier() << std::endl;
					}
				}

				if (all_passed)
					break;
			}


			if (k == hair_subSteps - 1) {
				// recover at the last substep
#pragma omp parallel for
				for (int i = 0; i < (int)m_strands.size(); i++)
				{
					m_steppers[i]->update();
					m_steppers[i]->recoverFromBeginning();
					residuals[i] = m_steppers[i]->getNewtonResidual();
				}
			}
			else {
#pragma omp parallel for
				for (int i = 0; i < (int)m_strands.size(); i++)
				{
					m_steppers[i]->update();
					m_steppers[i]->updateRodAccelerationAndVelocity();
					residuals[i] = m_steppers[i]->getNewtonResidual();
				}
			}

			if (m_strands.size())
				residualStats("step_dynamics", residuals);
		}

		// recover stepper Dt
#pragma omp parallel for
		for (int i = 0; i < (int)m_strands.size(); i++)
		{
			m_steppers[i]->setDt(substepDt);
			m_steppers[i]->setFraction(1.0);
			m_steppers[i]->updateRodAcceleration();

			if (m_params.m_solveLiquids) {
				m_steppers[i]->stepFlowAdvForce();
			}
		}
	}
	*/
	static const unsigned maxObjForOuterParallelism = 8 * omp_get_max_threads();
	/*
	void StrandImplicitManager::step_processCollisions(Scalar dt)
	{
		DebugStream(g_log, "") << " Postprocessing collisions ";

		// here we create the "penalty" collisions
		// 1. prune mutual contact
		//     a. random prune based on distance gauss and sort with distance
		//     b. accept the first maxNumCollisionsPerEdge collision 
		// 2. convert some contact to external contact (m_nonSPD ??)
		// 3. construct collision group by bfs
		// -> result stored in m_collidingGroups(m_elasticCollidingGroups) and m_collidingGroupsIdx(m_elasticCollidingGroupsIdx)
		computeCollidingGroups(m_mutualContacts, dt, false);
		computeCollidingGroups(m_elasticMutualContacts, dt, true);

		// construct 18 cells for every edge and extract one contact per cell based on their reletive velocity
		if (m_params.m_pruneExternalCollisions)
		{
			pruneExternalCollisions(m_externalContacts);
			pruneExternalCollisions(m_elasticExternalContacts);
		}

		// Deformation gradients at constraints
		unsigned nExternalContacts = 0;
		unsigned nElasticExternalContacts = 0;
#pragma omp parallel for reduction ( + : nExternalContacts,nElasticExternalContacts )
		for (int i = 0; i < (int)m_strands.size(); i++)
		{
			if (m_params.m_useDeterministicSolver && !m_params.m_pruneExternalCollisions)
			{
				std::sort(m_externalContacts[i].begin(), m_externalContacts[i].end());
				std::sort(m_elasticExternalContacts[i].begin(), m_elasticExternalContacts[i].end());
			}
			for (unsigned k = 0; k < m_externalContacts[i].size(); ++k)
			{
				setupDeformationBasis(m_externalContacts[i][k]);
			}
			for (unsigned k = 0; k < m_elasticExternalContacts[i].size(); ++k)
			{
				setupDeformationBasis(m_elasticExternalContacts[i][k]);
			}

			nExternalContacts += m_externalContacts[i].size();
			nElasticExternalContacts += m_elasticExternalContacts[i].size();
		}
		CopiousStream(g_log, "") << "Number of external contacts: " << nExternalContacts;
		CopiousStream(g_log, "") << "Number of elastic external contacts: " << nElasticExternalContacts;
		m_statExternalContacts = nExternalContacts + nElasticExternalContacts;

#pragma omp parallel for
		for (int i = 0; i < (int)m_collidingGroups.size(); i++)
		{
			CollidingGroup& cg = m_collidingGroups[i];

			if (cg.first.size() < maxObjForOuterParallelism)
			{
				for (unsigned k = 0; k < cg.second.size(); ++k)
				{
					setupDeformationBasis(cg.second[k]);
				}
			}
		}

		for (std::vector<CollidingGroup>::size_type i = 0; i < m_collidingGroups.size(); i++)
		{
			CollidingGroup& cg = m_collidingGroups[i];
			m_statMutualCollisions += cg.second.size();

			if (cg.first.size() >= maxObjForOuterParallelism)
			{
#pragma omp parallel for
				for (int k = 0; k < (int)cg.second.size(); ++k)
				{
					setupDeformationBasis(cg.second[k]);
				}
			}
		}

#pragma omp parallel for
		for (int i = 0; i < (int)m_elasticCollidingGroups.size(); i++)
		{
			CollidingGroup& cg = m_elasticCollidingGroups[i];

			if (cg.first.size() < maxObjForOuterParallelism)
			{
				for (unsigned k = 0; k < (int)cg.second.size(); ++k)
				{
					setupDeformationBasis(cg.second[k]);
				}
			}
		}

		for (std::vector<CollidingGroup>::size_type i = 0; i < m_elasticCollidingGroups.size(); i++)
		{
			CollidingGroup& cg = m_elasticCollidingGroups[i];
			m_statMutualCollisions += cg.second.size();

			if (cg.first.size() >= maxObjForOuterParallelism)
			{
#pragma omp parallel for
				for (int k = 0; k < (int)cg.second.size(); ++k)
				{
					setupDeformationBasis(cg.second[k]);
				}
			}
		}
	}
	*/
	void StrandImplicitManager::residualStats(const std::string& name, const std::vector< Scalar >& residuals)
	{
		using namespace boost::accumulators;
		typedef accumulator_set<Scalar, stats<tag::mean, tag::variance> > AccType;

		AccType acc;

		Scalar res_min = 1e+99;
		Scalar res_max = 0.0;

		for (const Scalar& s : residuals)
		{
			if (s < 0.0) continue;

			acc(s);
			res_min = std::min(s, res_min);
			res_max = std::max(s, res_max);
		}

		Scalar m = mean(acc);
		Scalar var = variance(acc);

		std::cout << "[" << name << " statistics: mean " << m << ", variance " << var << ", min " << res_min << ", max " << res_max << "]" << std::endl;
	}

//	void StrandImplicitManager::step_solveFlowFrictions(int total_num_substeps, int total_substep_id)
//	{
//		if (m_params.m_skipFlowFrictions)
//			return;
//
//		DebugStream(g_log, "") << " Solving ";
//		// m_contactProblemsStats.clear(); //DK: [stats]
//
//		// Dynamics solve
//
//		std::vector< Scalar > residuals(m_elasticCollidingGroups.size(), 1e+99);
//		std::vector< Scalar > newton_iters(m_elasticCollidingGroups.size(), 0.);
//		std::vector< Scalar > newton_residuals(m_strands.size(), -1.0);
//
//#pragma omp parallel for
//		for (int i = 0; i < (int)m_elasticCollidingGroups.size(); ++i)
//		{
//			if (m_elasticCollidingGroups[i].first.size() <= maxObjForOuterParallelism)
//			{
//				int iters = 0;
//				residuals[i] = solveCollidingGroup(m_elasticCollidingGroups[i], m_elasticExternalContacts, false, true, false, false, newton_residuals, iters, total_num_substeps, total_substep_id);
//				newton_iters[i] = iters;
//			}
//		}
//
//		for (int i = 0; i < (int)m_elasticCollidingGroups.size(); ++i)
//		{
//			if (m_elasticCollidingGroups[i].first.size() > maxObjForOuterParallelism)
//			{
//				int iters = 0;
//				residuals[i] = solveCollidingGroup(m_elasticCollidingGroups[i], m_elasticExternalContacts, false, true, false, false, newton_residuals, iters, total_num_substeps, total_substep_id);
//				newton_iters[i] = iters;
//			}
//		}
//
//		if (m_elasticCollidingGroups.size()) {
//			residualStats(std::string("Flow Colliding Group So-Bogus"), residuals);
//		}
//
//		std::vector< int > all_done(m_strands.size(), false);
//
//		std::vector< Scalar > residuals_single(m_strands.size(), -1.);
//		std::vector< Scalar > newton_iters_single(m_strands.size(), -1.);
//		std::vector< Scalar > newton_residuals_single(m_strands.size(), -1.0);
//
//		int count = 0;
//#pragma omp parallel for reduction(+:count)
//		for (int i = 0; i < (int)m_strands.size(); i++)
//		{
//			if (m_elasticCollidingGroupsIdx[i] == -1 && needsElasticExternalSolve(i))
//			{
//				int iters = 0;
//				residuals_single[i] = solveSingleObject(m_elasticExternalContacts, i, false, true, false, false, newton_residuals_single, iters, total_num_substeps, total_substep_id);
//				newton_iters_single[i] = iters;
//
//				count++;
//			}
//		}
//
//		if (count) {
//			residualStats(std::string("Single So-Bogus"), residuals_single);
//		}
//
//		Scalar maxWI = 0.;
//		int maxWI_idx = -1;
//		int maxWI_subidx = -1;
//
//		for (int i = 0; i < (int)m_steppers.size(); i++)
//		{
//			int subidx = -1;
//			// const Scalar Wi = m_steppers[i]->maxAdditionalImpulseNorm(subidx);
//			const Scalar Wi = m_steppers[i]->maxCollisionImpulseNorm(subidx);
//			//            maxWI = std::max(maxWI, Wi);
//			if (Wi > maxWI) {
//				maxWI = Wi;
//				maxWI_idx = i;
//				maxWI_subidx = subidx;
//			}
//		}
//
//		std::cout << "Max WI (Elastic): " << (maxWI / m_dt) << " @ " << maxWI_idx << ", " << maxWI_subidx << std::endl;
//	}

	// DK: solve all collisions and then finalize (relax theta's and accept state update).
//	void StrandImplicitManager::step_solveCollisions(int total_num_substeps, int total_substep_id)
//	{
//		DebugStream(g_log, "") << " Solving ";
//		// m_contactProblemsStats.clear(); //DK: [stats]
//
//		// Dynamics solve
//
//		std::vector< Scalar > residuals(m_collidingGroups.size(), 1e+99);
//		std::vector< Scalar > newton_iters(m_collidingGroups.size(), 0.);
//		std::vector< Scalar > newton_residuals(m_strands.size(), -1.0);
//
//#pragma omp parallel for
//		for (int i = 0; i < (int)m_collidingGroups.size(); ++i)
//		{
//			if (m_collidingGroups[i].first.size() <= maxObjForOuterParallelism)
//			{
//				int iters = 0;
//				residuals[i] = solveCollidingGroup(m_collidingGroups[i], m_externalContacts, false, false, true, false, newton_residuals, iters, total_num_substeps, total_substep_id);
//				newton_iters[i] = iters;
//			}
//		}
//
//		for (int i = 0; i < (int)m_collidingGroups.size(); ++i)
//		{
//			if (m_collidingGroups[i].first.size() > maxObjForOuterParallelism)
//			{
//				int iters = 0;
//				residuals[i] = solveCollidingGroup(m_collidingGroups[i], m_externalContacts, false, false, true, false, newton_residuals, iters, total_num_substeps, total_substep_id);
//				newton_iters[i] = iters;
//			}
//		}
//
//		if (m_params.m_useAdditionalExternalFailSafe) {
//			// redo it again without mutual collisions
//#pragma omp parallel for
//			for (int i = 0; i < (int)m_collidingGroups.size(); ++i)
//			{
//				if (m_collidingGroups[i].first.size() <= maxObjForOuterParallelism)
//				{
//					int iters = 0;
//					residuals[i] = solveCollidingGroup(m_collidingGroups[i], m_externalContacts, false, false, true, true, newton_residuals, iters, total_num_substeps, total_substep_id);
//					newton_iters[i] = iters;
//				}
//			}
//
//			for (int i = 0; i < (int)m_collidingGroups.size(); ++i)
//			{
//				if (m_collidingGroups[i].first.size() > maxObjForOuterParallelism)
//				{
//					int iters = 0;
//					residuals[i] = solveCollidingGroup(m_collidingGroups[i], m_externalContacts, false, false, true, true, newton_residuals, iters, total_num_substeps, total_substep_id);
//					newton_iters[i] = iters;
//				}
//			}
//		}
//
//		if (m_collidingGroups.size()) {
//			residualStats(std::string("Colliding Group So-Bogus"), residuals);
//			if (m_params.m_useNonlinearContacts) {
//				residualStats(std::string("Colliding Group Newton Iters"), newton_iters);
//				residualStats(std::string("Colliding Group Newton Residual"), newton_residuals);
//			}
//		}
//
//		std::vector< int > all_done(m_strands.size(), false);
//
//		std::vector< Scalar > residuals_single(m_strands.size(), -1.);
//		std::vector< Scalar > newton_iters_single(m_strands.size(), -1.);
//		std::vector< Scalar > newton_residuals_single(m_strands.size(), -1.0);
//
//		int count = 0;
//#pragma omp parallel for reduction(+:count)
//		for (int i = 0; i < (int)m_strands.size(); i++)
//		{
//			if (m_collidingGroupsIdx[i] == -1 && needsExternalSolve(i))
//			{
//				int iters = 0;
//				residuals_single[i] = solveSingleObject(m_externalContacts, i, false, false, true, false, newton_residuals_single, iters, total_num_substeps, total_substep_id);
//				newton_iters_single[i] = iters;
//
//				count++;
//			}
//		}
//
//		if (count) {
//			residualStats(std::string("Single So-Bogus"), residuals_single);
//			if (m_params.m_useNonlinearContacts) {
//				residualStats(std::string("Single Newton Iters"), newton_iters_single);
//				residualStats(std::string("Single Newton Residual"), newton_residuals_single);
//			}
//		}
//
//		Scalar maxWI = 0.;
//		int maxWI_idx = -1;
//		int maxWI_subidx = -1;
//
//		for (int i = 0; i < (int)m_steppers.size(); i++)
//		{
//			int subidx = -1;
//			// const Scalar Wi = m_steppers[i]->maxAdditionalImpulseNorm(subidx);
//			const Scalar Wi = m_steppers[i]->maxCollisionImpulseNorm(subidx);
//			//            maxWI = std::max(maxWI, Wi);
//			if (Wi > maxWI) {
//				maxWI = Wi;
//				maxWI_idx = i;
//				maxWI_subidx = subidx;
//			}
//		}
//
//		std::cout << "Max WI: " << (maxWI / m_dt) << " @ " << maxWI_idx << ", " << maxWI_subidx << std::endl;
//
//#pragma omp parallel for
//		for (int i = 0; i < (int)m_strands.size(); i++)
//		{
//			m_steppers[i]->finalize();
//		}
//	}

	//bool StrandImplicitManager::needsExternalSolve(unsigned strandIdx) const
	//{
	//	return !m_externalContacts[strandIdx].empty();
	//}

	//bool StrandImplicitManager::needsElasticExternalSolve(unsigned strandIdx) const
	//{
	//	return !m_elasticExternalContacts[strandIdx].empty();
	//}
	/*
	bool StrandImplicitManager::assembleBogusFrictionProblem(CollidingGroup& collisionGroup,
		bogus::MecheFrictionProblem& mecheProblem,
		std::vector<ProximityCollisions>& externalContacts,
		std::vector<unsigned>& globalIds,
		std::vector<ProximityCollision*>& colPointers, VecXx& vels, VecXx& worldImpulses, VecXx& impulses, VecXx& adhesions, VecXx& filters,
		VecXu& startDofs, VecXu& nDofs, int& numSubSys, bool herschelBulkleyProblem, bool ignoreMutualCollision)
	{
		Timer tt("solver", false, "assemble");

		unsigned nContacts = ignoreMutualCollision ? 0 : collisionGroup.second.size();

		unsigned nMutualContacts = nContacts;

		std::vector<unsigned> nndofs;
		unsigned nSubsystems = collisionGroup.first.size();
		numSubSys = nSubsystems;
		globalIds.reserve(nSubsystems);

		nndofs.resize(nSubsystems);
		nDofs.resize(nSubsystems);
		startDofs.resize(nSubsystems);

		unsigned dofCount = 0;
		unsigned subCount = 0;

		for (IndicesMap::const_iterator it = collisionGroup.first.begin(); it != collisionGroup.first.end(); ++it)
		{ // Computes the number of DOFs per strand and the total number of contacts
			startDofs[subCount] = dofCount;
			const unsigned sIdx = it->first;
			globalIds.push_back(sIdx);
			nDofs[subCount] = m_steppers[sIdx]->velocities().rows();
			nndofs[subCount] = m_steppers[sIdx]->velocities().rows();
			nContacts += externalContacts[sIdx].size();
			dofCount += m_steppers[sIdx]->velocities().rows();
			++subCount;
		}
		unsigned nExternalContact = nContacts - nMutualContacts;
		colPointers.resize(nContacts, NULL);

#pragma omp parallel for
		for (int i = 0; i < (int)globalIds.size(); ++i)
		{ // Prepare (rewind) the steppers
			ImplicitStepper* stepper = m_steppers[globalIds[i]];
			stepper->rewind();
		}

		std::vector < strandsim::SymmetricBandMatrixSolver<double, 10>* > MassMat;
		MassMat.resize(nSubsystems);

		VecXx forces(dofCount);
		VecXx dof(dofCount);

		vels.resize(dofCount);
		worldImpulses.resize(dofCount);
		impulses.resize(nContacts * 3);
		
		worldImpulses.setZero();
		vels.setZero();
		impulses.setZero();

		// we don't use adhesion for the frictional directions,
		// but we allocate their space here for computational easiness
		adhesions.resize(nContacts * 3);

		VecXx mu(nContacts);
		VecXx yields(nContacts);
		VecXx etas(nContacts);
		VecXx powers(nContacts);

		filters.resize(nContacts);

		std::vector < Eigen::Matrix< double, 3, 3 > > E;
		E.resize(nContacts);

		VecXx u_frees(nContacts * 3);

		std::vector<int> ObjA(nContacts);
		std::vector<int> ObjB(nContacts);

		std::vector < SparseRowMatx* > H_0;
		std::vector < SparseRowMatx* > H_1;
		H_0.resize(nContacts);
		H_1.resize(nContacts);

		unsigned objectId = 0, collisionId = 0, currDof = 0;
		bool oneNonSPD = false;

		// Setting up objects and adding external constraints
		for (IndicesMap::const_iterator it = collisionGroup.first.begin(); it != collisionGroup.first.end(); ++it)
		{
			const unsigned sIdx = it->first;
			ImplicitStepper& stepper = *m_steppers[sIdx];

			if (stepper.notSPD())
				oneNonSPD = true;

			dof.segment(currDof, stepper.m_strand.getCurrentDegreesOfFreedom().size()) = stepper.m_strand.getCurrentDegreesOfFreedom();

			if (!m_params.m_useImpulseMethod)
			{
				strandsim::JacobianSolver* M = &stepper.linearSolver();
				MassMat[objectId++] = M;

				forces.segment(currDof, stepper.rhs().size()) = -stepper.rhs();
				currDof += stepper.rhs().size();
			}
			else
			{
				// Pass a "pure" mass matrix - impulse problem
				strandsim::JacobianSolver* M = &stepper.massMatrixLinearSolver();
				MassMat[objectId++] = M;

				forces.segment(currDof, stepper.impulse_rhs().size()) = -stepper.impulse_rhs();
				currDof += stepper.impulse_rhs().size();
			}

			ProximityCollisions& externalCollisions = externalContacts[sIdx];

			for (int i = 0; i < (int)externalCollisions.size(); ++i)
			{
				ProximityCollision& c = externalCollisions[i];
				const Scalar frictionCoeff = c.mu;

				mu[collisionId] = frictionCoeff;
				E[collisionId] = c.transformationMatrix;
				u_frees.segment<3>((collisionId) * 3) = c.objects.first.freeVel - c.objects.second.freeVel;
				ObjA[collisionId] = (int)it->second;
				ObjB[collisionId] = -1;
				filters[collisionId] = m_steppers[c.objects.first.globalIndex]->getStrand().collisionParameters().m_impulseMaxNorm * m_dt;
				H_0[collisionId] = c.objects.first.defGrad;
				H_1[collisionId] = NULL;
				impulses.segment<3>((collisionId) * 3) = c.force;
				colPointers[collisionId++] = &c;
				// std::cout << "Col passed into SoBogus: " << std::endl;
				// c.print( std::cout );
			}
		}

		if (!ignoreMutualCollision) {
#pragma omp parallel for
			for (int i = 0; i < (int)collisionGroup.second.size(); ++i)
			{
				// Setting up mutual constraints
				// Solid touching
				ProximityCollision& collision = collisionGroup.second[i];
				const int oId1 = collisionGroup.first.find(collision.objects.first.globalIndex)->second;
				const int oId2 = collisionGroup.first.find(collision.objects.second.globalIndex)->second;
				const Scalar frictionCoeff = oneNonSPD ? 0. : collision.mu;
				mu[collisionId + i] = frictionCoeff;
				E[collisionId + i] = collision.transformationMatrix;
				u_frees.segment<3>((collisionId + i) * 3) = collision.objects.first.freeVel - collision.objects.second.freeVel;
				ObjA[collisionId + i] = oId1;
				ObjB[collisionId + i] = oId2;
				filters[collisionId + i] = m_steppers[collision.objects.first.globalIndex]->getStrand().collisionParameters().m_impulseMaxNorm * m_dt;

				H_0[collisionId + i] = collision.objects.first.defGrad;
				H_1[collisionId + i] = collision.objects.second.defGrad;
				impulses.segment<3>((collisionId + i) * 3) = collision.force;
				colPointers[collisionId + i] = &collision;
			}
			assert(collisionId + collisionGroup.second.size() == nContacts);
		}

		m_solverStat.m_assembleTime += tt.elapsed();
		mecheProblem.fromPrimal(nSubsystems, nndofs, dof / m_dt, MassMat, forces, adhesions, filters,
			nContacts, mu, yields, etas, powers, E, u_frees, &ObjA[0], &ObjB[0], H_0, H_1, m_params.m_useImpulseMethod);

		return true;
	}

	int StrandImplicitManager::postProcessBogusFrictionProblem(bool solveDryFrictions, CollidingGroup& collisionGroup,
		const bogus::MecheFrictionProblem& mecheProblem, const std::vector<unsigned>& globalIds,
		const std::vector<ProximityCollision*>& colPointers, VecXx& vels, VecXx& worldImpulses, VecXx& impulses,
		VecXu& startDofs, VecXu& nDofs, std::vector< Scalar >& newton_residuals, int total_num_substeps, int total_substep_id)
	{
		int iter = 0;
		assert(globalIds.size() == startDofs.size());
		assert(globalIds.size() == nDofs.size());

		m_globalIds.clear();
		m_globalIds = globalIds;

		// here we need complete Newton solve to ensure stability
#pragma omp parallel for
		for (int i = 0; i < globalIds.size(); ++i)
		{
			const unsigned sIdx = globalIds[i];
			const unsigned subSystem = collisionGroup.first.find(sIdx)->second;
			m_steppers[sIdx]->additionalImpulses() += worldImpulses.segment(startDofs[subSystem], nDofs[subSystem]);
		}

		//        std::cout << "[PPBFP 0]" << std::endl;

		if (solveDryFrictions)
		{
			if (m_params.m_useNonlinearContacts) {
				//                std::cout << "[PPBFP 1]" << std::endl;
#pragma omp parallel for
				for (int i = 0; i < (int)globalIds.size(); ++i)
				{
					const unsigned sIdx = globalIds[i];
					m_steppers[sIdx]->updateAdditionalInertia();
				}

				int hair_subSteps = std::max(1, (int)ceil((Scalar)m_params.m_rodSubSteps / (Scalar)total_num_substeps));

				//                std::cout << "[PPBFP 2]" << std::endl;
				for (int k = 0; k < hair_subSteps; ++k) {
					//                    std::cout << "start rod substep " << k << std::endl;

					std::vector< unsigned char > passed(globalIds.size(), false);

					//                    std::cout << "[PPBFP 21]" << std::endl;

					while (true) {
#pragma omp parallel for
						for (int i = 0; i < (int)globalIds.size(); ++i) {
							if (!passed[i]) {
								const unsigned sIdx = globalIds[i];
								const unsigned subSystem = collisionGroup.first.find(sIdx)->second;
								m_steppers[sIdx]->startSubstep(k, m_dt / hair_subSteps, 1.0 / hair_subSteps);
								m_steppers[sIdx]->prepareDynamics();
								if (k == 0)
									m_steppers[sIdx]->newVelocities() = vels.segment(startDofs[subSystem], nDofs[subSystem]);
								m_steppers[sIdx]->prepareSolveNonlinear();
							}
						}

						std::vector< int > all_done(globalIds.size(), false);

						while (true) {
#pragma omp parallel for
							for (int i = 0; i < (int)globalIds.size(); i++)
							{
								if (!passed[i] && !all_done[i]) {
									const unsigned sIdx = globalIds[i];

									m_steppers[sIdx]->prepareNewtonIteration();
									all_done[i] = m_steppers[sIdx]->performNewtonIteration();
								}
							}

							bool done = true;
							for (int i = 0; i < (int)globalIds.size(); i++)
								if (!passed[i])
									done = done && all_done[i];

							// std::cout << "step dynamics iter: " << iter << std::endl;
							++iter;

							if (done) break;
						}

#pragma omp parallel for
						for (int i = 0; i < (int)globalIds.size(); i++)
						{
							if (!passed[i]) {
								const unsigned sIdx = globalIds[i];
								passed[i] = m_steppers[sIdx]->postSolveNonlinear();
							}
						}

						// check if all passed
						bool all_passed = true;

						for (int i = 0; i < (int)globalIds.size(); i++) {
							all_passed = all_passed && passed[i];

							if (!passed[i]) {
								std::cout << "Strand " << globalIds[i] << " failed! Redo with " << m_steppers[globalIds[i]]->getStretchMultiplier() << std::endl;
							}
						}

						if (all_passed)
							break;
					}

					//                    std::cout << "[PPBFP 22]" << std::endl;

					if (k == hair_subSteps - 1) {
						// recover at the last substep
#pragma omp parallel for
						for (int i = 0; i < (int)globalIds.size(); i++)
						{
							const unsigned sIdx = globalIds[i];
							m_steppers[sIdx]->update();
							m_steppers[sIdx]->recoverFromBeginning();
							newton_residuals[sIdx] = m_steppers[sIdx]->getNewtonResidual();
						}
					}
					else {
#pragma omp parallel for
						for (int i = 0; i < (int)globalIds.size(); i++)
						{
							const unsigned sIdx = globalIds[i];
							m_steppers[sIdx]->update();
							m_steppers[sIdx]->updateRodAccelerationAndVelocity();
							newton_residuals[sIdx] = m_steppers[sIdx]->getNewtonResidual();
						}
					}

					//                    std::cout << "[PPBFP 22]" << std::endl;
				}
				//                std::cout << "[PPBFP 3]" << std::endl;
				// recover stepper Dt
#pragma omp parallel for
				for (int i = 0; i < (int)globalIds.size(); i++)
				{
					const unsigned sIdx = globalIds[i];
					m_steppers[sIdx]->setDt(m_dt);
					m_steppers[sIdx]->setFraction(1.0);
					m_steppers[sIdx]->updateRodAcceleration();

					// must perform this to clean up RHS and LHS for possibly future use
					m_steppers[sIdx]->solveLinear();
				}

				//                std::cout << "[PPBFP 4]" << std::endl;
			}
			else {
				for (int i = 0; i < (int)globalIds.size(); ++i)
				{
					const unsigned sIdx = globalIds[i];
					const unsigned subSystem = collisionGroup.first.find(sIdx)->second;
					m_steppers[sIdx]->newVelocities() = vels.segment(startDofs[subSystem], nDofs[subSystem]);
					m_steppers[sIdx]->updateRHSwithImpulse(m_steppers[sIdx]->additionalImpulses());
					m_steppers[sIdx]->update(true);
				}
			}
		}
		else {
			// only update rhs
#pragma omp parallel for
			for (int i = 0; i < (int)globalIds.size(); ++i)
			{
				const unsigned sIdx = globalIds[i];
				const unsigned subSystem = collisionGroup.first.find(sIdx)->second;
				m_steppers[sIdx]->newVelocities() = vels.segment(startDofs[subSystem], nDofs[subSystem]);
				m_steppers[sIdx]->updateRHSwithImpulse(m_steppers[sIdx]->additionalImpulses());
				m_steppers[sIdx]->update(true);
			}
		}

		//        std::cout << "[PPBFP 5]" << std::endl;

		if (!m_params.m_simulationManager_limitedMemory)
		{ // If we are low on memory, do not bother storing collisions
			for (int i = 0; i < colPointers.size(); ++i)
			{
				if (!colPointers[i]) continue;

				ProximityCollision& col = *colPointers[i];
				col.force = impulses.segment<3>(i * 3); // ??? Missing Adhesion forces
				// std::cout << "impulse[" << i << "]: " << impulses.segment<3>(i) << std::endl;
				m_collisionDatabase.insert(col);
			}
		}

		//        std::cout << "[PPBFP 6]" << std::endl;

		return iter;
	}
	

	Scalar StrandImplicitManager::solveBogusFrictionProblem(bogus::MecheFrictionProblem& mecheProblem,
		const std::vector<unsigned>& globalIds, bool asFailSafe, bool herschelBulkleyProblem, bool doFrictionShrinking, VecXx& vels, VecXx& worldImpulses, VecXx& impulses, int& numSubSys)
	{
		bogus::MecheFrictionProblem::Options options;
		options.maxThreads = m_params.m_numberOfThreads;
		options.maxIters = m_params.m_gaussSeidelIterations;
		options.cadouxIters = 0;
		options.tolerance = m_params.m_gaussSeidelTolerance;
		options.useInfinityNorm = false;	// if ture, use max res, else use average res
		options.algorithm = m_params.m_bogusAlgorithm;
		options.ignoreVelocity = false;		// calculate reletive velocity (u)

		options.gsRegularization = 0.0;

		//        std::cout << "[solveBogus 0]" << std::endl;

		const Scalar residual = mecheProblem.solve(
			impulses,   // impulse guess and returned impulse
			vels,       // returned velocities
			worldImpulses,
			options,
			false, // static problem
			herschelBulkleyProblem,
			doFrictionShrinking,
			0.0,
			m_params.m_useImpulseMethod); // cadoux iters

		//        std::cout << "[solveBogus 1]" << std::endl;

		return residual;
	}
	

	Scalar StrandImplicitManager::solveSingleObject(std::vector<ProximityCollisions>& externalContacts, unsigned objectIdx, bool asFailSafe, bool herschelBulkleyProblem, bool updateVelocity, bool ignoreMutualCollision, std::vector< Scalar >& newtonResiduals, int& numNewtonIters, int total_num_substeps, int total_substep_id)
	{
		CollidingGroup collisionGroup;
		collisionGroup.first[objectIdx] = 0;

		return solveCollidingGroup(collisionGroup, externalContacts, asFailSafe, herschelBulkleyProblem, updateVelocity, ignoreMutualCollision, newtonResiduals, numNewtonIters, total_num_substeps, total_substep_id);
	}

	// DK: another place the failsafes kick in -- here we have it outside the GS solve *LOOK HERE*
	Scalar StrandImplicitManager::solveCollidingGroup(CollidingGroup& collisionGroup, std::vector<ProximityCollisions>& externalContacts, bool asFailSafe,
		bool herschelBulkleyProblem, bool updateVelocity, bool ignoreMutualCollision, std::vector< Scalar >& newtonResiduals, int& numNewtonIters, int total_num_substeps, int total_substep_id)
	{
		if (collisionGroup.first.empty())
			return 0.0;

		std::vector<unsigned> globalIds;
		std::vector<ProximityCollision*> colPointers;

		Scalar res = 1e+99;

		numNewtonIters = 0;

		{
			VecXx vels;
			VecXx impulses;
			VecXx worldImpulses;
			VecXx adhesions;
			VecXx filters;
			bogus::MecheFrictionProblem mecheProblem;
			VecXu startDofs;
			VecXu nDofs;
			int numSubSystems;

			if (assembleBogusFrictionProblem(collisionGroup, mecheProblem, externalContacts, globalIds, colPointers, vels, worldImpulses, impulses, adhesions, filters, startDofs, nDofs, numSubSystems, herschelBulkleyProblem, ignoreMutualCollision))
			{
				res = solveBogusFrictionProblem(mecheProblem, globalIds, asFailSafe, herschelBulkleyProblem, !herschelBulkleyProblem && !m_params.m_useApproxRodElasticFriction, vels, worldImpulses, impulses, numSubSystems);
				// accumulate worldImpulse to rhs
				// set new velocity
				// update currentStates with FutureStates
				// use FTL or if stretch energy is larger than threshold, set m_lastStepWasRejected = true
				// update collision database
				Timer tt("solver", false, "pose precess");
				numNewtonIters = postProcessBogusFrictionProblem(updateVelocity, collisionGroup, mecheProblem, globalIds, colPointers, vels, worldImpulses, impulses, startDofs, nDofs, newtonResiduals, total_num_substeps, total_substep_id);
				m_solverStat.m_poseProcessTime += tt.elapsed();
			}

			m_solverStat.addStat(mecheProblem, collisionGroup);
		}

		/* If either the solver result was really bad or one strand was stretching,
		 we trigger a hierarchy of failsafes.

		 - First retry with the non-linear solver if this was not already the case
		 - Then discard all mutual contacts, solve each strand with non-linear solver
		 - Finally solve each strand with linear solver but with length constraints
		 - If that is still failing, trust the ImplicitStepper's reprojection
		 

		 // DK: failsafes here:
		bool mustRetry = false;
		for (int i = 0; i < globalIds.size(); ++i)
		{
			const unsigned sIdx = globalIds[i];
			if (m_steppers[sIdx]->lastStepWasRejected())
			{
				mustRetry = true;
				break;
			}
		}

		if (mustRetry)
		{

			// DK: if not failsafe & not nonlinear always try nonlinear as failsafe (possibly not rewinded?)
			// DK: only if nonlinear as failsafe and not nonlinear always ... so need to fix this to to do the right thing
			// DK: for now always using "nonlinear" due to changes just made in friction solver so this is unused:
			if (globalIds.size() > 1)
			{
				//ContactStream( g_log, "" ) << " Dropping mutual constraints ";
				InfoStream(g_log, "") << "\033[31;1m Dropping mutual constraints \033[m\n";

				// Failed again, drop mutual collisions
				std::vector< Scalar > single_iters(globalIds.size(), 0.0);

#pragma omp parallel for
				for (int i = 0; i < globalIds.size(); ++i)
				{
					const unsigned sIdx = globalIds[i];
					if (m_steppers[sIdx]->lastStepWasRejected())
					{
						int iters = 0;
						res = solveSingleObject(externalContacts, globalIds[i], true, herschelBulkleyProblem, updateVelocity, ignoreMutualCollision, newtonResiduals, iters, total_num_substeps, total_substep_id);
						single_iters[i] = iters;
					}
				}

				residualStats(std::string("retry colliding group iters"), single_iters);
			}
		}

		if (!asFailSafe)
		{
			for (int i = 0; i < colPointers.size(); ++i)
			{
				if (!colPointers[i]) continue;

				ProximityCollision& col = *colPointers[i];
				delete col.objects.first.defGrad;
				col.objects.first.defGrad = NULL;
				if (col.objects.second.globalIndex != -1)
				{
					delete col.objects.second.defGrad;
					col.objects.second.defGrad = NULL;
				}
			}

		}

		return res;
	}
	*/
	void StrandImplicitManager::updateParameters(const SimulationParameters& params)
	{

		g_log->SetMinSeverity((MsgInfo::Severity) params.m_logLevel);

		// Enforce desired or maximum number of threads
		{
#ifdef WIN32
			SYSTEM_INFO sysinfo;
			GetSystemInfo(&sysinfo);
#endif
			const int numThreads =
				m_params.m_numberOfThreads > 0 ? m_params.m_numberOfThreads :
#ifdef WIN32
				sysinfo.dwNumberOfProcessors;
#else
				sysconf(_SC_NPROCESSORS_ONLN);
#endif
#if defined(_OPENMP)
			omp_set_num_threads(numThreads);
#endif
		}

		m_params = params;
	}

	struct ProxColPointerCompare
	{
		bool operator()(const ProximityCollision* lhs, const ProximityCollision* rhs) const
		{
			return *lhs < *rhs;
		}
	};

	struct PairHash
	{
		template <class T1, class T2>
		std::size_t operator()(std::pair<T1, T2> const& pair) const
		{
			std::size_t h1 = std::hash<T1>()(pair.first);
			std::size_t h2 = std::hash<T2>()(pair.second);

			return h1 ^ h2;
		}
	};

	void StrandImplicitManager::pruneCollisions(const ProximityCollisions& origMutualCollisions,
		ProximityCollisions& mutualCollisions, const Scalar stochasticPruning)
	{
		// Transform unacceptable self collision info external contacts
		std::vector<std::map<unsigned, std::set<const ProximityCollision*, ProxColPointerCompare> > > tentativeCols(m_strands.size());
		
		for (int i = 0; i < origMutualCollisions.size(); ++i)
		{
			const ProximityCollision& proxyCol = origMutualCollisions[i];
			const unsigned s1 = proxyCol.objects.first.globalIndex;
			const unsigned s2 = proxyCol.objects.second.globalIndex;

			if (tentativeCols[s1][s2].insert(&proxyCol).second) {
				// std::cout << "JOIN tentativeCols:: "<< &proxyCol <<" ::["<< s1 << "]["<< s2 <<"]:" << proxyCol.normal << std::endl;
			}
			else {
				std::cout << "taken by: " << *(tentativeCols[s1][s2].find(&proxyCol)) << std::endl;
				std::cout << "INSERT FAILED on tentativeCols || " << &proxyCol << " ::[" << s1 << "][" << s2 << "]:" << proxyCol.normal << std::endl;
			}
		}

		// Stochastic pruning:
		// On each edge, compute the distribituion of the distances of all contacts
		// Then if there are more than one collision on this edge, we will prune
		// them with a probabily based on their position in this distribution

		using namespace boost::accumulators;
		typedef accumulator_set<double, stats<tag::mean, tag::variance> > AccType;

		std::vector<std::map<unsigned, AccType> > distanceDistribution(m_strands.size());

		for (unsigned s1 = 0; s1 < m_strands.size(); ++s1)
		{
			for (auto s2It = tentativeCols[s1].begin(); s2It != tentativeCols[s1].end(); ++s2It)
			{
				unsigned s2 = s2It->first;
				auto& cols = s2It->second;
				// std::cout << " distanceDistribution Cols size vs origMutualCollisions size:: " << cols.size() << " / " << origMutualCollisions.size() << std::endl;

				unsigned c = 0;
				for (auto cIt = cols.begin(); cIt != cols.end(); ++cIt)
				{
					const unsigned e1 = (*cIt)->objects.first.vertex;
					const unsigned e2 = (*cIt)->objects.second.vertex;

					distanceDistribution[s1][e1]((*cIt)->distance);
					distanceDistribution[s2][e2]((*cIt)->distance);
				}
			}
		}

		boost::mt19937 generator;
		//    generator.seed( 1.e6 * getTime() ) ;
		boost::uniform_real<> distribution(0, 1);

		for (unsigned s1 = 0; s1 < m_strands.size(); ++s1)
		{
			std::vector< std::map< Scalar, const ProximityCollision* > > closest_per_edges(m_strands[s1]->getNumVertices());

			for (auto s2It = tentativeCols[s1].begin(); s2It != tentativeCols[s1].end(); ++s2It)
			{
				unsigned s2 = s2It->first;
				auto& cols = s2It->second;

				for (const auto& col : cols)
				{
					closest_per_edges[col->objects.first.vertex][col->distance] = col;
				}

				/*std::map< Scalar, const ProximityCollision* > sorted_coll;
				std::unordered_set< std::pair< unsigned, unsigned >, PairHash > vertex_pair_set;

				for (const auto& col : cols)
				{
					sorted_coll[col->distance] = col;
				}

				std::vector<const ProximityCollision*> new_collisions;
				new_collisions.reserve(cols.size());

				for (const auto& item : sorted_coll)
				{
					unsigned v1 = item.second->objects.first.vertex;
					unsigned v2 = item.second->objects.second.vertex;

					if (vertex_pair_set.find(std::pair< unsigned, unsigned >(v1 + 1, v2)) == vertex_pair_set.end()
						&& vertex_pair_set.find(std::pair< unsigned, unsigned >(v1 - 1, v2)) == vertex_pair_set.end()
						&& vertex_pair_set.find(std::pair< unsigned, unsigned >(v1, v2 + 1)) == vertex_pair_set.end()
						&& vertex_pair_set.find(std::pair< unsigned, unsigned >(v1, v2 - 1)) == vertex_pair_set.end()
						&& vertex_pair_set.find(std::pair< unsigned, unsigned >(v1 - 1, v2 - 1)) == vertex_pair_set.end()
						&& vertex_pair_set.find(std::pair< unsigned, unsigned >(v1 + 1, v2 + 1)) == vertex_pair_set.end()
						&& vertex_pair_set.find(std::pair< unsigned, unsigned >(v1 - 1, v2 + 1)) == vertex_pair_set.end()
						&& vertex_pair_set.find(std::pair< unsigned, unsigned >(v1 + 1, v2 - 1)) == vertex_pair_set.end())
					{
						new_collisions.push_back(item.second);
						vertex_pair_set.insert(std::pair< unsigned, unsigned >(v1, v2));
					}
				}

				for (auto& col : new_collisions)
				{
					closest_per_edges[col->objects.first.vertex][col->distance] = col;
				}*/

				//std::vector<bool> accept(cols.size());
				//std::vector<bool> taken(m_strands[s2]->getNumVertices());
				//std::vector<bool> seen(m_strands[s1]->getNumVertices());

				//Vec3x prevNormal;

				//// First accept a collision for each of the first strand edges
				//unsigned c = 0;
				//prevNormal.setZero();
				//for (auto cIt = cols.begin(); cIt != cols.end(); ++cIt)
				//{
				//	const unsigned e1 = (*cIt)->objects.first.vertex;
				//	const unsigned e2 = (*cIt)->objects.second.vertex;

				//	if (!seen[e1])
				//	{

				//		const Scalar m1 = mean(distanceDistribution[s1][e1]);
				//		const Scalar v1 = variance(distanceDistribution[s1][e1]);
				//		boost::normal_distribution<> gaussian(m1, std::sqrt(v1));
				//		boost::variate_generator<decltype(generator), decltype(gaussian)> vg(
				//			generator, gaussian);

				//		if (isSmall(v1) || stochasticPruning * (*cIt)->distance <= vg())
				//		{

				//			seen[e1] = true;
				//			taken[e2] = true;
				//			accept[c] = true;
				//		}
				//	}
				//	++c;
				//}

				//// Then accept a collision for each of the not-yet-seen second strand edges
				//c = cols.size() - 1;
				////            prevNormal.setZero() ; // keep las normal
				//for (auto cIt = cols.rbegin(); cIt != cols.rend(); ++cIt)
				//{
				//	const unsigned e2 = (*cIt)->objects.second.vertex;
				//	if (!taken[e2])
				//	{
				//		const Scalar m2 = mean(distanceDistribution[s2][e2]);
				//		const Scalar v2 = variance(distanceDistribution[s2][e2]);

				//		boost::normal_distribution<> gaussian(m2, std::sqrt(v2));
				//		boost::variate_generator<decltype(generator), decltype(gaussian)> vg(
				//			generator, gaussian);

				//		if (isSmall(v2) || stochasticPruning * (*cIt)->distance <= vg())
				//		{
				//			taken[e2] = true;
				//			accept[c] = true;
				//		}
				//	}
				//	--c;
				//}

				//c = 0;
				//for (auto cIt = cols.begin(); cIt != cols.end(); ++cIt)
				//{
				//	if (accept[c])
				//	{
				//		const ProximityCollision& col = **cIt;

				//		closest_per_edges[(*cIt)->objects.first.vertex][(*cIt)->distance] = *cIt;
				//	}
				//	++c;
				//}
			}

			// final accept according to distance order
			const int maxNumCollisionsPerEdge = m_strands[s1]->collisionParameters().m_maxNumCollisionsPerEdge;
			for (auto& sorter : closest_per_edges)
			{
				int c = 0;
				for (auto& pair : sorter) {
					if (c >= maxNumCollisionsPerEdge)
						break;

					const ProximityCollision& col = *(pair.second);
					mutualCollisions.push_back(col);
					++c;
				}
			}
		}

		std::cout << "After Pruning: " << mutualCollisions.size() << " self collisions ( down from " << origMutualCollisions.size() << " ) \n";
	}

//	void StrandImplicitManager::computeCollidingGroups(const ProximityCollisions& origMutualCollisions,
//		Scalar dt, bool elastic)
//	{
//		ProximityCollisions mutualCollisions;
//		mutualCollisions.reserve(origMutualCollisions.size());
//
//		if (m_params.m_pruneSelfCollisions)
//		{
//			// here we put penalty collisions aside
//			pruneCollisions(origMutualCollisions, mutualCollisions, m_params.m_stochasticPruning, elastic);
//		}
//		else
//		{
//			for (int i = 0; i < (int)origMutualCollisions.size(); ++i)
//			{
//				const ProximityCollision& proxyCol = origMutualCollisions[i];
//				const unsigned s1 = proxyCol.objects.first.globalIndex;
//				const unsigned s2 = proxyCol.objects.second.globalIndex;
//
//				if (m_steppers[s1]->refusesMutualContacts()
//					|| m_steppers[s2]->refusesMutualContacts())
//				{
//					if (elastic) {
//						ProximityCollision copy(proxyCol);
//						makeElasticExternalContact(copy, m_steppers[s2]->refusesMutualContacts());
//					}
//					else {
//						ProximityCollision copy(proxyCol);
//						makeExternalContact(copy, m_steppers[s2]->refusesMutualContacts());
//					}
//
//				}
//				else
//				{
//					mutualCollisions.push_back(proxyCol);
//				}
//			}
//			if (m_params.m_useDeterministicSolver)
//			{
//				std::sort(mutualCollisions.begin(), mutualCollisions.end());
//			}
//		}
//
//		CopiousStream(g_log, "") << "Number of mutual collisions: " << mutualCollisions.size();
//
//		// For each strand, list all other contacting ones
//
//		std::vector<std::deque<unsigned> > objsGroups(m_strands.size());
//		for (int i = 0; i < (int)mutualCollisions.size(); ++i)
//		{
//			const ProximityCollision& proxyCol = mutualCollisions[i];
//			const unsigned s1 = proxyCol.objects.first.globalIndex;
//			const unsigned s2 = proxyCol.objects.second.globalIndex;
//
//			objsGroups[s1].push_back(s2);
//			objsGroups[s2].push_back(s1);
//		}
//
//		// Extract connected subgraphs from the global constraint graph
//		// std::vector<int> collidingGroupsIdx: Index of colliding group in which each strand should be. Can be -1.
//		auto bfsGraph = [&](std::vector<int>& collidingGroupsIdx, std::vector<CollidingGroup>& collidingGroups) {
//			for (unsigned s1 = 0; s1 < objsGroups.size(); ++s1)
//			{
//				if (collidingGroupsIdx[s1] != -1 || objsGroups[s1].empty())
//					continue;
//
//				const unsigned groupIdx = collidingGroups.size();
//
//				collidingGroups.push_back(CollidingGroup());
//				CollidingGroup& cg = collidingGroups.back();
//
//				collidingGroupsIdx[s1] = groupIdx;
//				cg.first[s1] = 0;
//
//				std::deque<unsigned> toVisit = objsGroups[s1];
//
//				while (toVisit.size())
//				{
//					const unsigned s2 = toVisit.front();
//					toVisit.pop_front();
//
//					if (collidingGroupsIdx[s2] != -1)
//						continue;
//
//					collidingGroupsIdx[s2] = groupIdx;
//					cg.first[s2] = 0;
//
//					toVisit.insert(toVisit.end(), objsGroups[s2].begin(), objsGroups[s2].end());
//				}
//			}
//
//			for (int i = 0; i < (int)mutualCollisions.size(); ++i)
//			{
//				const ProximityCollision& mutualCollision = mutualCollisions[i];
//				const unsigned s1 = mutualCollision.objects.first.globalIndex;
//
//				collidingGroups[collidingGroupsIdx[s1]].second.push_back(mutualCollision);
//			}
//
//#pragma omp parallel for
//			for (int i = 0; i < (int)collidingGroups.size(); ++i)
//			{
//				unsigned k = 0;
//				IndicesMap& indices = collidingGroups[i].first;
//				for (IndicesMap::iterator it = indices.begin(); it != indices.end(); ++it)
//				{
//					it->second = k++;
//				}
//			}
//		};
//
//		if (elastic) {
//			bfsGraph(m_elasticCollidingGroupsIdx, m_elasticCollidingGroups);
//		}
//		else {
//			bfsGraph(m_collidingGroupsIdx, m_collidingGroups);
//		}
//
//		DebugStream(g_log, "") << "Number of colliding groups = " << m_collidingGroups.size();
//	}

	void StrandImplicitManager::pruneExternalCollisions(std::vector<ProximityCollisions>& externalCollisions)
	{

		// Create 2 x 6 buckets on each edge:
		// 2 buckets corresponding to voronoi cells on local abscissa
		// 6 corresponding to collision's normnal angle w.r.t major radius
		// For each bucket, keep only the first occuring collision
		// ( i.e. the one with the highest freeVel on the normal direction )

		const Scalar minRhs = -1.e6;

		typedef Eigen::Matrix<double, 3, 6> Buckets;
		typedef Eigen::Matrix<int, 3, 6> AcceptedCols;

		Buckets empty;
		empty.setConstant(minRhs);
		AcceptedCols none;
		none.setConstant(-1);

		const Scalar invBucketWidth = Buckets::ColsAtCompileTime / (2. * M_PI);

		unsigned nPruned = 0;

#pragma omp parallel for reduction ( + : nPruned )
		for (int i = 0; i < (int)externalCollisions.size(); ++i)
		{
			ProximityCollisions& externalContacts = externalCollisions[i];
			if (externalContacts.empty())
				continue;

			std::vector<Buckets> buckets(m_strands[i]->getNumVertices(), empty);
			std::vector<AcceptedCols> accepted(m_strands[i]->getNumVertices(), none);

			for (unsigned j = 0; j < externalContacts.size(); ++j)
			{
				const ProximityCollision& externalContact = externalContacts[j];
				const int vtx = externalContact.objects.first.vertex;

				Buckets& vtxBucket = buckets[vtx];

				const Scalar angle = m_strands[i]->getSignedAngleToMajorRadius(vtx,
					externalContact.normal) + M_PI;

				const int angleBucket = clamp((int)std::floor(angle * invBucketWidth), 0,
					Buckets::ColsAtCompileTime - 1);
				assert(angleBucket >= 0);

				const int absBucket = clamp(
					(int)std::ceil(externalContact.objects.first.abscissa * Buckets::RowsAtCompileTime),
					0, Buckets::RowsAtCompileTime - 1);

				const Scalar wn = externalContact.objects.second.freeVel.dot(externalContact.normal);
				if (wn > vtxBucket(absBucket, angleBucket))
				{
					vtxBucket(absBucket, angleBucket) = wn;
					accepted[vtx](absBucket, angleBucket) = j;
				}
			}

			ProximityCollisions prunedExternalContacts;

			for (unsigned vtx = 0; vtx < accepted.size(); ++vtx)
			{
				AcceptedCols& acc = accepted[vtx];

				for (unsigned r = 0; r < acc.rows(); ++r)
				{
					for (unsigned c = 0; c < acc.cols(); ++c)
					{
						if (acc(r, c) != -1)
						{
							prunedExternalContacts.push_back(externalContacts[acc(r, c)]);
						}
					}
				}
			}

			const unsigned np = externalContacts.size() - prunedExternalContacts.size();

			nPruned += np;

			externalContacts.swap(prunedExternalContacts);
		}

		std::cout << "Pruned " << nPruned << " external collisions" << std::endl;

	}

	void StrandImplicitManager::print(const Timings& timings) const
	{
		SubStepTimings tot;

		for (auto stepT = timings.begin(); stepT != timings.end(); ++stepT)
		{
			for (auto substepT = stepT->begin(); substepT != stepT->end(); ++substepT)
			{
				tot = tot + *substepT;
			}
		}

		print<InfoStream>(tot / timings.size());
	}

	void StrandImplicitManager::print(const StrandImplicitManager::StepTimings& timings) const
	{
		SubStepTimings tot;
		for (int i = 0; i < timings.size(); ++i)
		{
			if (!i)
				tot = timings[i];
			else
				tot = tot + timings[i];
		}

		print<InfoStream>(tot);
	}

	template<typename StreamT>
	void StrandImplicitManager::print(const StrandImplicitManager::SubStepTimings& timings) const
	{
		StreamT(g_log, "Timing") 
			<< "PR: " << timings.prepare 
			<< " PCD: " << timings.proximityCollisions
			<< " DY: " << timings.dynamics << " LI: " << timings.linearIteration 
			<< " CCD: " << timings.continousTimeCollisions 
			<< " PC: " << timings.processCollisions 
			<< " SC: " << timings.solve;

		StreamT(g_log, "Timing") << "Total: " << timings.sum();
	}

	template<typename StreamT>
	void StrandImplicitManager::print(const StrandImplicitManager::CDTimings& timings) const
	{
		StreamT(g_log, "CD breakdown timing")
			<< "updateHashMap: " << timings.updateHashMap
			<< " HH Proxmity: " << timings.processHashMap
			<< " buildBVH: " << timings.buildBVH
			<< " CCD: " << timings.findCollisionsBVH
			<< " assemble: " << timings.narrowPhase
			<< " total: " << timings.sum();
	}

	template<typename StreamT>
	void StrandImplicitManager::print(const StrandImplicitManager::SolverStat& stat) const
	{
		StreamT(g_log, "Solver breakdown timing") << " Ass: " << stat.m_assembleTime << " PrimalCopy: " << stat.m_primalCopyTime
			<< " Minv: " << stat.m_MinvTime << " ComputeDual: " << stat.m_computeDualTime
			<< " Solve: " << stat.m_solveTime << " PostProc: " << stat.m_poseProcessTime
			<< " Total: " << stat.sum();

		for (int i = 0; i < stat.m_collisionSize.size(); ++i) {
			std::cout << "Collision group " << i << " (strands: " << stat.m_collisionSize[i].first 
				<< " contacts: " << stat.m_collisionSize[i].second << "):\n";
			int iter = 0;
			for (auto it = stat.m_solverStat[i].begin(); it != stat.m_solverStat[i].end(); ++it) {
				std::cout << "\tIter: " << iter++ << " err: " << it->first << " time: " << it->second << '\n';
			}
		}
	}
	/*
	template<typename StreamT>
	void StrandImplicitManager::printNewtonSolverBreakdownTiming() const
	{
		double pre = 0, LHS = 0, RHS = 0, store = 0, solve = 0, post = 0;
		double timesM = 0, F = 0, composeRhs = 0, addM = 0, J = 0;

		for (int i = 0; i < (int)m_strands.size(); i++) {
			m_steppers[i]->getTimings(pre, RHS, timesM, F, composeRhs, LHS, J, addM, store, solve, post);
		}

		StreamT(g_log, "Newton breakdown timing") << "pre: " << pre 
			<< " RHS: " << RHS << " (timesM: " << timesM << " F: " << F << " add together: " << composeRhs
			<< ") LHS: " << LHS << " (J: " << J << " addM: " << addM 
			<< ") store&fab: " << store << " solve: " << solve << " post: " << post 
			<< " total: " << pre + RHS + LHS + store + solve + post;
	}
	*/

	StrandImplicitManager::SubStepTimings operator+(const StrandImplicitManager::SubStepTimings& lhs,
		const StrandImplicitManager::SubStepTimings& rhs)
	{
		StrandImplicitManager::SubStepTimings sum;
		sum.prepare = lhs.prepare + rhs.prepare;
		sum.proximityCollisions = lhs.proximityCollisions + rhs.proximityCollisions;
		sum.dynamics = lhs.dynamics + rhs.dynamics;
		sum.continousTimeCollisions = lhs.continousTimeCollisions + rhs.continousTimeCollisions;
		sum.linearIteration = lhs.linearIteration + rhs.linearIteration;
		sum.processCollisions = lhs.processCollisions + rhs.processCollisions;
		sum.solve = lhs.solve + rhs.solve;

		return sum;
	}

	StrandImplicitManager::SubStepTimings operator/(const StrandImplicitManager::SubStepTimings& lhs,
		const Scalar rhs)
	{
		StrandImplicitManager::SubStepTimings avg = lhs;
		avg.prepare /= rhs;
		avg.proximityCollisions /= rhs;
		avg.dynamics /= rhs;
		avg.linearIteration /= rhs;
		avg.continousTimeCollisions /= rhs;
		avg.processCollisions /= rhs;
		avg.solve /= rhs;

		return avg;
	}

	void StrandImplicitManager::drawContacts() const
	{
		m_collisionDatabase.draw(m_strands);
	}

	void StrandImplicitManager::exportStrandRestShapes(const std::string& fileName) const
	{
		std::ofstream out(fileName.c_str());
		out.setf(std::ios::fixed);

		for (auto strand = m_strands.begin(); strand != m_strands.end(); ++strand)
		{
			std::copy((*strand)->getRestLengths().begin(), (*strand)->getRestLengths().end(),
				std::ostream_iterator<Scalar>(out, " "));
			out << '\n';
			std::copy((*strand)->getRestKappas().begin(), (*strand)->getRestKappas().end(),
				std::ostream_iterator<Vec4x>(out, " "));
			out << '\n';
			std::copy((*strand)->getRestTwists().begin(), (*strand)->getRestTwists().end(),
				std::ostream_iterator<Scalar>(out, " "));
			out << '\n';
		}
	}

	/************************************** My Step *****************************************/

	bool g_one_iter = false;
	std::mutex g_iter_mutex;
	std::condition_variable g_cv;

	void StrandImplicitManager::step()
	{
		SubStepTimings timings;

		std::cout << "[Prepare Simulation Step]" << std::endl;
		Timer timer("step", false);
		step_prepare(m_dt);
#pragma omp parallel for
		for (int i = 0; i < m_steppers.size();++i) {
			m_steppers[i]->initSolver(m_dt);
		}
		timings.prepare += timer.elapsed();

		if (m_params.m_solveCollision && !m_params.m_useCTRodRodCollisions) {
			std::cout << "[Proximity Collision Detection]" << std::endl;
			timer.restart();
			step_prepareCollision();
			/* TODO: add hair/mesh proximity test
			 * EdgeFaceIntercestion in m_collisionDetector->findCollision()
			 */
			setupHairHairCollisions(m_dt);
			// doProximityMeshHairDetection(m_dt);
			timings.proximityCollisions += timer.elapsed();

			timer.restart();
			step_processCollisions();
			timings.processCollisions += timer.elapsed();
		}

		std::cout << "[Step Dynamics]" << std::endl;
		for (int k = 0; k < m_params.m_nonlinearIterations; ++k) {
			timer.restart();
#pragma omp parallel for
			for (int i = 0; i < m_steppers.size(); ++i) {
				m_steppers[i]->linearize();
			}
			timings.dynamics += timer.elapsed();

			std::vector<bool> passed(m_strands.size(), false);
			bool all_done = true;

			for (int i = 0; i < m_params.m_linearIterations; ++i) {
				timer.restart();
#pragma omp parallel for
				for (int s = 0; s < m_steppers.size(); ++s) {
					passed[s] = m_steppers[s]->solveLinear();
					m_steppers[s]->postSolveLinear();
				}
				all_done = true;
				for (bool p : passed) {
					all_done = all_done && p;
				}
				if (all_done) break;

				timings.linearIteration += timer.elapsed();

				if (g_one_iter) {
					if (m_substep_callback) m_substep_callback->executeCallback();
					std::vector<Scalar> delta_v(m_strands.size());
					for (int s = 0; s < m_strands.size(); ++s) {
						delta_v[s] = m_steppers[s]->getVelocityDiff();
					}
					residualStats("Iter", delta_v);
					{
						std::unique_lock<std::mutex> lk(g_iter_mutex);
						g_one_iter = false;
					}
					std::unique_lock<std::mutex> lk(g_iter_mutex);
					g_cv.wait(lk, [] {return g_one_iter; });
				}

				if (m_params.m_solveCollision) {
					if (m_params.m_useCTRodRodCollisions && i == 0 && k == 0) {
						timer.restart();
						step_prepareCollision();

						step_continousCollisionDetection();
						timings.continousTimeCollisions += timer.elapsed();

						timer.restart();
						step_processCollisions();
						timings.processCollisions += timer.elapsed();
					}

					if (m_statTotalContacts > 0) {
						timer.restart();
						step_solveCollisions();
						timings.solve += timer.elapsed();
					}
				}
			}
			if (all_done) break;
		}

		m_timings.back().push_back(timings);

		step_postCollision();

		std::vector<Scalar> residuals(m_strands.size());
		std::vector<Scalar> delta_v(m_strands.size());
		int numUnsolved = 0;
		for (int i = 0; i < m_strands.size(); ++i) {
			residuals[i] = m_steppers[i]->getBestResidual();
			delta_v[i] = m_steppers[i]->getVelocityDiff();
			if (delta_v[i] > m_params.m_velocityDiffTolerance) ++numUnsolved;
		}
		std::cout << "Unsolved Strand: " << numUnsolved << std::endl;
		residualStats("Nonlinear Solver Residual", residuals);
		residualStats("Nonlinear Solver Delta v", delta_v);

		printMemStats();
	}

	void StrandImplicitManager::step_prepareCollision()
	{
		m_collisionDatabase.ageAll();

		m_mutualContacts.clear();
		for (int i = 0; i < m_strands.size(); ++i)
		{
			m_externalContacts[i].clear();
		}

		m_statTotalContacts = 0;
	}
	
	void StrandImplicitManager::step_continousCollisionDetection()
	{
		// now currentState = x_{t+1}, futureState = x_t
		// must input false to add AABB and futureAABB
		m_collisionDetector->buildBVH(false);

		EdgeFaceIntersection::s_doProximityDetection = true;

		// findCollision assume currentState = x_{t+1}

		// ignoreCTRodRod = false : CCD of edges among hairs
		//     -> add EdgeEdgeCollision to m_continuousTimeCollisions
		// ignoreContinuousTime = false : CCD of vertex-face and edge-face (CCD of edge and 3 edges of the face)
		//     -> add VertexFaceCollision and EdgeFaceCollision to m_continuousTimeCollisions
		// ignoreProximity = false : edge-face intersection test (s_doProximityDetection = true) using position before unconstraint update
		//     -> add EdgeFaceIntersection to m_proximityCollisions
		m_collisionDetector->findCollisions(false, false, false);

		// compact m_continuousTimeCollisions in ProximityCollision, and add these collisions to m_externalContacts
		doContinuousTimeDetection(m_dt);
		// compact m_proximityCollisions in ProximityCollision, and add these collisions to m_externalContacts
		doProximityMeshHairDetection(m_dt);

		m_collisionDetector->clear();
	}

	void StrandImplicitManager::step_processCollisions()
	{
		/* prune and sort m_mutualContacts and m_externalContacts
		 * compute deformation basis
		 */

		 // ---------------- mutual ---------------------
		ProximityCollisions mutualCollisions;
		mutualCollisions.reserve(m_mutualContacts.size());

		//for (int i = 0; i < (int)m_mutualContacts.size(); ++i)
		//{
		//	const ProximityCollision& proxyCol = m_mutualContacts[i];
		//	const unsigned s1 = proxyCol.objects.first.globalIndex;
		//	const unsigned s2 = proxyCol.objects.second.globalIndex;

		//	if (m_steppers[s1]->refusesMutualContacts() || m_steppers[s2]->refusesMutualContacts()) {
		//		ProximityCollision copy(proxyCol);
		//		makeExternalContact(copy, m_steppers[s2]->refusesMutualContacts());
		//	}
		//	else {
		//		mutualCollisions.push_back(proxyCol);
		//	}
		//}
		
		if (m_params.m_pruneSelfCollisions) {
			pruneCollisions(m_mutualContacts, mutualCollisions, m_params.m_stochasticPruning);
			m_mutualContacts.swap(mutualCollisions);
		}
		if (m_params.m_useDeterministicSolver) {
			std::sort(m_mutualContacts.begin(), m_mutualContacts.end());
		}

#pragma omp parallel for
		for (int i = 0; i < m_mutualContacts.size(); ++i) {
			setupDeformationBasis(m_mutualContacts[i]);
		}

		// --------------- external ------------------
		if (m_params.m_pruneExternalCollisions) 
		{
			pruneExternalCollisions(m_externalContacts);
		}
		int nExt = 0;
#pragma omp parallel for reduction(+: nExt)
		for (int i = 0; i < (int)m_strands.size(); ++i) 
		{
			if (m_params.m_useDeterministicSolver && !m_params.m_pruneExternalCollisions) {
				std::sort(m_externalContacts[i].begin(), m_externalContacts[i].end());
			}
			for (int j = 0; j < m_externalContacts[i].size(); ++j) {
				setupDeformationBasis(m_externalContacts[i][j]);
				++nExt;
			}
		}

		m_statTotalContacts = nExt + m_mutualContacts.size();

		InfoStream(g_log, "Contact") << "external: " << nExt << " mutual: " << m_mutualContacts.size();
	}

	void StrandImplicitManager::step_solveCollisions()
	{
		std::vector<ProximityCollision*> colPointers(m_statTotalContacts);
		std::vector<Scalar> residuals(m_statTotalContacts);

		/*auto rng = std::default_random_engine{};
		std::shuffle(m_mutualContacts.begin(), m_mutualContacts.end(), rng);*/

//#pragma omp parallel for
		for (int i = 0; i < m_mutualContacts.size(); ++i) {
			ProximityCollision& collision = m_mutualContacts[i];
			ImplicitStepper& first_stepper = *m_steppers[collision.objects.first.globalIndex];
			ImplicitStepper& second_stepper = *m_steppers[collision.objects.second.globalIndex];
			int first_vert = collision.objects.first.vertex;
			int second_vert = collision.objects.second.vertex;

			// Assemble
			Mat14x M = Mat14x::Zero();
			M.block<7, 7>(0, 0) = first_stepper.getMass().segment<7>(4 * first_vert).asDiagonal();
			M.block<7, 7>(7, 7) = second_stepper.getMass().segment<7>(4 * second_vert).asDiagonal();

			Mat3x14x H = Mat3x14x::Zero();
			H.block<3, 7>(0, 0) = MatXx(*collision.objects.first.defGrad);
			H.block<3, 7>(0, 7) = -MatXx(*collision.objects.second.defGrad);

			Vec14x f;
			f.segment<7>(0) = first_stepper.getb_hat().segment<7>(4 * first_vert)
				+ first_stepper.getCollisionImpulse().segment<7>(4 * first_vert);
			f.segment<7>(7) = second_stepper.getb_hat().segment<7>(4 * second_vert)
				+ second_stepper.getCollisionImpulse().segment<7>(4 * second_vert);

			Vec3x uf = Vec3x::Zero();

			Mat3x E = collision.transformationMatrix;

			Scalar mu = collision.mu;

			bogus::MutualContactSolver collision_solver(M, H, f, uf, E, mu);

			// Solve
			Vec14x vel, impulse;
			Scalar res = collision_solver.solve(collision.force, vel, impulse);

			// Store velocity and impulse
			first_stepper.accumulateCollision(first_vert, impulse.segment<7>(0));
			second_stepper.accumulateCollision(second_vert,  impulse.segment<7>(7));

			residuals[i] = res;
			colPointers[i] = &collision;
		}

		int col_idx = m_mutualContacts.size();
#pragma omp parallel for
		for (int i = 0; i < m_strands.size(); ++i) {
			for (int j = 0; j < m_externalContacts[i].size(); ++j) {
				ProximityCollision& collision = m_externalContacts[i][j];
				ImplicitStepper& stepper = *m_steppers[collision.objects.first.globalIndex];
				int vert = collision.objects.first.vertex;

				Mat7x M = stepper.getMass().segment<7>(4 * vert).asDiagonal();
				Mat3x7x H = MatXx(*collision.objects.first.defGrad);
				Vec3x uf = -collision.objects.second.freeVel;
				Mat3x E = collision.transformationMatrix;
				Scalar mu = collision.mu;

				Vec7x f = stepper.getb_hat().segment<7>(4 * vert) 
					+ stepper.getCollisionImpulse().segment<7>(4 * vert);

				bogus::ExternalContactSolver solver(M, H, f, uf, E, mu);

				Vec7x vel, impulse;
				residuals[col_idx] = solver.solve(collision.force, vel, impulse);
				stepper.accumulateCollision(vert, impulse);

				colPointers[col_idx++] = &collision;
			}
		}
		
		if (g_one_iter) {
			for (int i = 0; i < residuals.size(); ++i) 
			{
				ContactStream(g_log, "") << colPointers[i]->force << " @ (" << colPointers[i]->objects.first.globalIndex << ", "
					<< colPointers[i]->objects.first.vertex << ") (" << colPointers[i]->objects.second.globalIndex
					<< ", " << colPointers[i]->objects.second.vertex << ") res = " << residuals[i];
			}

			{
				std::unique_lock<std::mutex> lk(g_iter_mutex);
				g_one_iter = false;
			}
			std::unique_lock<std::mutex> lk(g_iter_mutex);
			g_cv.wait(lk, [] {return g_one_iter; });
		}

		// delete collision with small impulse
		m_mutualContacts.erase(std::remove_if(m_mutualContacts.begin(), m_mutualContacts.end(), [&](const ProximityCollision& col) {
			return isSmall(col.force.norm());
		}), m_mutualContacts.end());

		int nExt = 0;
		for (int i = 0; i < m_strands.size(); ++i)
		{
			m_externalContacts[i].erase(std::remove_if(m_externalContacts[i].begin(), m_externalContacts[i].end(),
				[&](const ProximityCollision& col) { return isSmall(col.force.norm()); }),
				m_externalContacts[i].end());
			nExt += m_externalContacts[i].size();
		}

		m_statTotalContacts = m_mutualContacts.size() + nExt;
	}

	void StrandImplicitManager::step_postCollision()
	{
		std::vector<ProximityCollision*> colPointers(m_statTotalContacts);

		for (int i = 0; i < m_mutualContacts.size(); ++i) {
			colPointers[i] = &m_mutualContacts[i];
		}

		int idx = m_mutualContacts.size();
		for (int i = 0; i < m_strands.size(); ++i) {
			for (auto& col : m_externalContacts[i]) {
				colPointers[idx++] = &col;
			}
		}

		// update collision database
		if (!m_params.m_simulationManager_limitedMemory) {
			for (int i = 0; i < colPointers.size(); ++i) {
				m_collisionDatabase.insert(*colPointers[i]);
			}
		}

		// print impulse
		for (int i = 0; i < colPointers.size(); ++i) {
			ContactStream(g_log, "") << colPointers[i]->force << " @ (" 
				<< colPointers[i]->objects.first.globalIndex << ", " << colPointers[i]->objects.first.vertex << ") (" 
				<< colPointers[i]->objects.second.globalIndex << ", " << colPointers[i]->objects.second.vertex << ")";
		}
	}

} // namespace strandsim
