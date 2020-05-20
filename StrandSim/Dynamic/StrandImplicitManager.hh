/**
 * \copyright 2014 Danny Kaufman, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef STRANDSIM_STRANDIMPLICITMANAGER_HH
#define STRANDSIM_STRANDIMPLICITMANAGER_HH

#include "SimulationParameters.hh"
#include "CohesionTable.hh"

#include "../Core/Definitions.hh"
#include "../Utils/SpatialHashMapFwd.hh"
#include "../Collision/ProximityCollision.hh"

#include "../../bogus/Interfaces/MecheEigenInterface.hpp"

#include <vector>
#include <list>
#include <set>
#include <memory>
#include <mutex>
#include <condition_variable>

namespace strandsim
{

class ElasticStrand;
class ImplicitStepper;
class DOFScriptingController;
class CollisionDetector;
class ElementProxy;
struct DistanceFieldObject;
class MeshScriptingController;
class FluidScriptingController;
class ConstraintScriptingController;
class CollisionBase;
class LevelSetForce;

extern bool g_one_iter;
extern std::mutex g_iter_mutex;
extern std::condition_variable g_cv;

class StrandImplicitManager
{
public:
    //! Map between a index in the simulation to an index in a colliding group
    typedef std::map<unsigned, unsigned> IndicesMap;
    //! Colliding group: set of strands and contacts that should be solved together
    typedef std::pair<IndicesMap, ProximityCollisions> CollidingGroup;

    //! Simulation step timings -- in milliseconds
	struct SubStepCallback
	{
		virtual void executeCallback() = 0;
        virtual void projectConstraint(const std::vector<ImplicitStepper*>&) = 0;
	};

    struct SubStepTimings
    {
        SubStepTimings() :
            prepare( 0 ), proximityCollisions( 0 ), dynamics( 0 ), linearIteration( 0 ), continousTimeCollisions( 0 ),
            processCollisions( 0 ), solve( 0 )
        {
        }

        double prepare;                     // Executing controllers, etc
        double proximityCollisions;         // Hair hair collision detection
        double dynamics;                    // Assembling (+ pre-solving linear systems)
        double linearIteration;             // Iterative solver
        double continousTimeCollisions;     // Mesh hair collision detection
        double processCollisions;           // Updating collision gradients
        double solve;                       // Proper solve

        double sum() const
        {
            return prepare + proximityCollisions + dynamics + linearIteration + continousTimeCollisions 
                + processCollisions + solve;
        }
    };

    //! collision detection timings -- in milliseconds
    struct CDTimings
    {
        CDTimings() :
        buildBVH( 0. ), findCollisionsBVH( 0. ), narrowPhase( 0. ), updateHashMap ( 0. ), processHashMap( 0. )
        {
        }
        
        double buildBVH;            
        double findCollisionsBVH;  
        double narrowPhase;              
        double updateHashMap;
        double processHashMap;                 
        
        void reset()
        {
            buildBVH = 0.;
            findCollisionsBVH = 0.;
            narrowPhase = 0. ;
            updateHashMap = 0. ;
            processHashMap = 0. ;
        }

        double sum() const {
            return buildBVH + findCollisionsBVH + narrowPhase + updateHashMap + processHashMap;
        }
    };

    //! collision solver timings & residuals -- in milliseconds
    struct SolverStat
    {
        SolverStat() :
            m_assembleTime(0.), m_primalCopyTime(0.), m_MinvTime(0.), m_computeDualTime(0.), m_solveTime(0.),
            m_poseProcessTime(0.)
        {}

        double sum() const
        {
            return m_assembleTime + m_primalCopyTime + m_MinvTime + m_computeDualTime + m_solveTime + m_poseProcessTime;
        }

        void reset()
        {
            m_assembleTime = 0;
            m_primalCopyTime = 0;
            m_MinvTime = 0;
            m_computeDualTime = 0;
            m_solveTime = 0;
            m_poseProcessTime = 0;

            m_collisionSize.clear();
            m_solverStat.clear();
        }

        void addStat(const bogus::MecheFrictionProblem& meche, const CollidingGroup& cg)
        {
            m_primalCopyTime += meche.primalCopyTime();
            m_MinvTime += meche.MinvTime();
            m_computeDualTime += meche.computeDualTime();
            m_solveTime += meche.solveTime();

            m_collisionSize.push_back(std::make_pair(cg.first.size(), cg.second.size()));
            m_solverStat.push_back(meche.iterateStat());
        }

        double m_assembleTime;
        double m_primalCopyTime;
        double m_MinvTime;
        double m_computeDualTime;
        double m_solveTime;
        double m_poseProcessTime;

        std::vector<std::pair<int, int> > m_collisionSize;
        std::vector<bogus::MecheFrictionProblem::IterateStat> m_solverStat;
    };

    typedef std::vector<SubStepTimings> StepTimings;
    typedef std::list<StepTimings> Timings;

    StrandImplicitManager( const std::vector<ElasticStrand*>& strands,
            const std::map<std::pair<int, int>, std::set< std::pair<int, int> > >& collision_free,
            const std::vector<DistanceFieldObject>& fields,
            const std::vector< std::shared_ptr<MeshScriptingController> >& meshScripting_controllers,
            const std::vector< std::shared_ptr<FluidScriptingController> >& fluidScripting_controllers,
            const std::vector<ConstraintScriptingController*>& constraintScripting_controllers,
            Scalar startTime, Scalar dt, const SimulationParameters& params, SubStepCallback* sub_callback );

    virtual ~StrandImplicitManager();

    //! Step the simulation forward by m_dt
    /*! \param substepsInDt the number of simulation substeps to perform */
    void execute(int total_num_substeps, int total_substep_id, const Scalar total_substep_dt);

    //! Performs a simulation substep
    void step();
    
    Scalar getCFL() const;

    void setDt( Scalar dt )
    {
        m_dt = dt;
    }
    Scalar getDt() const
    {
        return m_dt;
    }
    Scalar getTime() const
    {
        return m_time;
    }
    const std::string& getOutputDirectory() const
    {
        return m_output_directory;
    }
    
    void setOutputDirectory( const std::string& dir )
    {
        m_output_directory = dir;
    }
    
    void setTime( Scalar m_time )
    {
        this->m_time = m_time;
    }
    void updateParameters( const SimulationParameters& params );

    void print( const StepTimings& timings ) const;
    void print( const Timings& timings ) const;

    const Timings& timings() const
    {
        return m_timings;
    }

    void drawContacts() const;

    const std::vector<ElasticStrand*>& getStrands() const
    {
        return m_strands;
    }

    std::shared_ptr<FluidScriptingController> getFluidScriptingController(int idx = 0) const
	{
		if(idx >= m_fluidScriptingControllers.size()) {
			return NULL;
		} else {
			return m_fluidScriptingControllers[idx];
		}
	}

private:
    void printMemStats();
    template<typename StreamT>
    void print( const SubStepTimings& timings ) const;
    template<typename StreamT>
    void print(const CDTimings& timings) const;
    template<typename StreamT>
    void print( const SolverStat& stat ) const;

    //! Prepares the substep, executes the externals objects controllers
    void step_prepare( Scalar dt );
    //! Clears collision set, ages collision database
    void step_prepareCollision();
    //! Does continous collision detection
    void step_continousCollisionDetection();
    //! Vertex-face collision detection of distance field objects
    void step_vertexFaceCollisionDetection();
    //! Edge-face collision detection
    void step_edgeFaceCollisionDetection();
    //! Sorts and prunes collisions set
    void step_processCollisions();
    //! Solves collisions and updates collision impulses
    void step_solveCollisions();
    //! Updates collision database and prints collision info
    void step_postCollision();

    void residualStats( const std::string& name, const std::vector< Scalar >& residuals );

    //! Hair/hair proximity collision detection
    void setupHairHairCollisions( Scalar dt );

    //! Proximity mesh/hair collision detection
    void doProximityMeshHairDetection( Scalar dt );
    //! Continuous-time mesh/hair collisision detection
    void doContinuousTimeDetection( Scalar dt );

    //! Adds an external contact on strand \p strIdx, edge \p edgeIdx, abscissa \p abscissa
    /*! \return whether this collision has been accepted */
    bool addExternalContact( const unsigned strIdx, const unsigned edgeIdx, const Scalar abscissa,
            const ProximityCollision& collision );
    //! Transform a mutual collision into an external contact on the (onFirstObject ? first : second) object
    void makeExternalContact(ProximityCollision& c, bool onFirstObject);

    //! Computes the deformation gradient of a strand at one contact point, ie dq/dx
    //void computeDeformationGradient( ProximityCollision::Object &object ) const;
    //! Setup the local frame for one contact and calls computeDeformationGradient() for each object
    void setupDeformationBasis( ProximityCollision &collision ) const;

    //! Discards mesh/hair collisions that are unlikely to be activated
    void pruneExternalCollisions( std::vector<ProximityCollisions>& contacts );

    //! Discards hair/hair collisions that are unlikely to be activated
    void pruneCollisions( const ProximityCollisions &origMutualCollisions,
            ProximityCollisions &mutualCollisions, const Scalar stochasticPruning );

    //! Updates collision times
    void updateCollisionTimes(const ProximityCollisions& mutualContacts);

    //! solve individual collision in local coordinate
    Vec3x solveOneCollision(const Vec3x& vel, Scalar mass, const Mat3x& R, Scalar mu);
    
    void exportStrandRestShapes( const std::string& fileName ) const;

    Scalar m_time; //!< Current time
    Scalar m_dt;   //!< Time per "frame"
    SimulationParameters m_params;
    std::string m_output_directory;

    const std::vector<ElasticStrand*>& m_strands;
    std::vector<ImplicitStepper*> m_steppers;
    const std::vector< DistanceFieldObject >& m_distanceFields;
    const std::vector< std::shared_ptr<MeshScriptingController> >& m_meshScriptingControllers;
    const std::vector< std::shared_ptr<FluidScriptingController> >& m_fluidScriptingControllers;
    const std::vector<ConstraintScriptingController*>& m_constraintScriptingControllers;

    std::vector<ElementProxy*> m_elementProxies; //!< List of all proxies that should be inserted in the BVH
    CollisionDetector* m_collisionDetector;        //!< BVH-based collision detector

    ProximityCollisions m_mutualContacts;           //!< List of all mutual contacts (edge/edge)
    std::vector<ProximityCollisions> m_externalContacts;         //!< List of all contacts between movable hair and fixed edges
    //! Structure for storing collisions and forces, useful for drawaing and warm-starting solver
    ProximityCollisionDatabase m_collisionDatabase;
    const std::map<std::pair<int, int>, std::set< std::pair<int, int> > >& m_collision_free;

    //!< Spatial Hash Map for hair/hair proximity collision detetection
    typedef SpatialHashMap<ElasticStrand, unsigned, true> SpatialHashMapT;
    SpatialHashMapT * m_hashMap;

	SubStepCallback* m_substep_callback;
    // Stat gathering 
    SolverStat m_solverStat;
    Timings m_timings;
    unsigned m_statTotalContacts;
    Scalar m_mem_usage_accu;
    Scalar m_mem_usage_divisor;
    CDTimings m_cdTimings;
};

StrandImplicitManager::SubStepTimings operator+( const StrandImplicitManager::SubStepTimings& lhs,
        const StrandImplicitManager::SubStepTimings& rhs );
StrandImplicitManager::SubStepTimings operator/( const StrandImplicitManager::SubStepTimings& lhs,
        const Scalar rhs );

} // namespace strandsim

#endif // STRANDSIM_STRANDIMPLICITMANAGER_HH
