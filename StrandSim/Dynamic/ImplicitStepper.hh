#ifndef STRANDSIM_ImplicitStepper_HH
#define STRANDSIM_ImplicitStepper_HH

#if defined(_OPENMP)
#include <omp.h>
#endif
#include "../Core/Definitions.hh"
#include "../Utils/SymmetricBandMatrixSolver.hh"

namespace strandsim
{
	class ElasticStrand;
	class StrandDynamicTraits;
	struct SimulationParameters;

	class ImplicitStepper
	{
	public:
		enum LinearSolverType { DIRECT, JACOBI, GAUSS_SEIDEL, CONJ_GRAD };

		ImplicitStepper(ElasticStrand& strand, const SimulationParameters& params);
		virtual ~ImplicitStepper();

		virtual void prepareStep(Scalar dt) = 0;
		virtual bool performOneIteration() = 0;
		virtual void postStep() = 0;

		void clearCollisionImpulse() { m_collisionImpulse.setZero(); }
		Scalar maxCollisionImpulseNorm(int& idx) const;

		void finalize() {}

		bool refusesMutualContacts() const { return m_notSPD; }
		Scalar getStretchMultiplier() const { return 1.0; }
		VecXx flowComponents() const { return VecXx::Ones(m_velocities.size()); }	// for StrandRender::computeFlowQuads and ProblemStepper::dumpRods

		VecXx& velocities() { return m_velocities; }
		const VecXx& velocities() const { return m_velocities; }

		const VecXx& getCollisionImpulse() const { return m_collisionImpulse; }

	protected:
		omp_lock_t m_lock;

		ElasticStrand& m_strand;
		const SimulationParameters& m_params;
		
		VecXx m_velocities;
		VecXx m_collisionImpulse;

		bool m_notSPD;
	};
}

#endif // !STRANDSIM_ImplicitStepper_HH
