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
		~ImplicitStepper();

		// Save states
		void initSolver(Scalar dt);
		// Linearize at v_{t+1}^{k} (compute A and b)
		void linearize();
		// Update v (Direct solver or one iteration in iterative solver)
		bool solveLinear();
		// update current state and compute \hat{b}
		void postSolveLinear();
		// Register velocities and impulses outputed by collision solver
		void accumulateCollision(int vid, const Vec7x& impulse);

		void clearCollisionImpulse() { m_collisionImpulse.setZero(); }
		Scalar maxCollisionImpulseNorm(int& idx) const;

		void finalize() {}

		bool refusesMutualContacts() const { return m_notSPD; }
		Scalar getStretchMultiplier() const { return 1.0; }
		VecXx flowComponents() const { return VecXx::Ones(m_velocities.size()); }	// for StrandRender::computeFlowQuads and ProblemStepper::dumpRods

		VecXx& velocities() { return m_velocities; }
		const VecXx& velocities() const { return m_velocities; }

		VecXx& newVelocities() { return m_newVelocities; }
		const VecXx& newVelocities() const { return m_newVelocities; }

		ElasticStrand& getStrand() { return m_strand; }
		const ElasticStrand& getStrand() const { return m_strand; }

		const JacobianMatrixType& getA() const { return m_A; }
		const VecXx& getb() const { return m_b; }
		const VecXx& getMass() const { return m_mass; }
		const VecXx& getb_hat() const { return m_b_hat; }
		const VecXx& getCollisionImpulse() const { return m_collisionImpulse; }
		Scalar getBestResidual() const { return m_bestError; }
		Scalar getVelocityDiff() const { return m_velocityDiff; }

	private:
		void directSolver(const VecXx& b);
		void JacobiStep(const VecXx& b);
		void GaussSeidelStep(const VecXx& b);
		void ConjgradStep(const VecXx& b);

		Scalar m_dt;

		omp_lock_t m_lock;

		ElasticStrand& m_strand;
		const SimulationParameters& m_params;
		StrandDynamicTraits& m_dynamics;

		LinearSolverType m_linearSolverType;
		int m_iteration;
		Scalar m_tolerance;

		JacobianSolver m_directSolver;
		bool m_notSPD;

		VecXx m_mass;
		JacobianMatrixType m_A;
		VecXx m_b;
		VecXx m_b_hat;

		VecXx m_velocities;
		VecXx m_newVelocities;
		VecXx m_savedVelocities;
		VecXx m_substepVelocities;
		
		VecXx m_collisionImpulse;

		VecXx m_residual;
		Scalar m_bestError;
		Scalar m_velocityDiff;
	};
}

#endif // !STRANDSIM_ImplicitStepper_HH
