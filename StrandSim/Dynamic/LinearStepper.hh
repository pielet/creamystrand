#ifndef STRANDSIM_LINEARSTEPPER_HH
#define STRANDSIM_LINEARSTEPPER_HH

#include "../Core/Definitions.hh"
#include "../Utils/SymmetricBandMatrixSolver.hh"

namespace strandsim
{
	class ElasticStrand;
	class StrandDynamicTraits;
	struct SimulationParameters;

	class LinearStepper
	{
		/* solve Ax = b
		 * A = M + H(x_t) * h^2, b = Mv_{t+1} - f(x_t) + r
		 * optional solver type: direct solver (LLT decomposition), Jacobi, Gauss-Seidel, Conjugate Gradient
		 */
	public:
		enum SolverType { DIRECT, JACOBI, GAUSS_SEIDEL, CONJ_GRAD };

		LinearStepper(ElasticStrand& strand, const SimulationParameters& params);
		~LinearStepper();

		// Linearize at v_{t+1}=0 (compute A and b)
		void initSolver(Scalar dt);
		// Direct solver or one iteration in iterative solver
		bool solveLinear();
		// Register velocities and impulses outputed by collision solver
		void accumulateCollision(int vid, const Vec7x& vel, const Vec7x& impulse);
		// Use solved velocities to update current state
		void addCollisionVelocity();

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
		const VecXx& getResidual() const { return m_residual; }
		const VecXx& getExactResidual() const { return m_resWithoutImpulse; }
		Scalar getBestResidual() const { return m_bestError; }

	private:
		void directSolver(const VecXx& b);
		void JacobiStep(const VecXx& b);
		void GaussSeidelStep(const VecXx& b);
		void ConjgradStep(const VecXx& b);

		void updateCurrentState();

		Scalar m_dt;

		ElasticStrand& m_strand;
		const SimulationParameters& m_params;
		StrandDynamicTraits& m_dynamics;

		SolverType m_solverType;
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
		
		VecXx m_collisionImpulse;
		VecXx m_collisionVelocity;

		VecXx m_resWithoutImpulse;
		VecXx m_residual;
		Scalar m_bestError;
	};
}

#endif // !STRANDSIM_LINEARSTEPPER_HH
