#include "ImplicitStepper.hh"
#include "../Core/ElasticStrand.hh"
#include "SimulationParameters.hh"
#include "DOFScriptingController.hh"
#include "StrandDynamicTraits.hh"

namespace strandsim
{
	ImplicitStepper::ImplicitStepper(ElasticStrand& strand, const SimulationParameters& params) :
		m_dt(0.),
		m_strand(strand),
		m_params(params),
		m_dynamics(strand.dynamics()),
		m_linearSolverType(params.m_linearSolverType),
		m_iteration(0),
		m_tolerance(params.m_velocityDiffTolerance),
		m_notSPD(false),
		m_velocities(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
		m_newVelocities(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
		m_savedVelocities(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
		m_substepVelocities(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
		m_collisionImpulse(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
		m_residual(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
		m_bestError(1e+99)
	{
		m_dynamics.computeDOFMasses();
		omp_init_lock(&m_lock);
	}


	ImplicitStepper::~ImplicitStepper()
	{

	}

	void ImplicitStepper::initSolver(Scalar dt)
	{
		m_dt = dt;

		m_mass = m_dynamics.getDOFMasses();

		m_savedVelocities = m_velocities;
		m_velocities.setZero();
		m_newVelocities.setZero();

		m_strand.setSavedDegreesOfFreedom(m_strand.getCurrentDegreesOfFreedom());

		m_dynamics.getDisplacements().setZero();
		m_dynamics.getAccelerations().setZero();

		m_iteration = 0;
		m_bestError = 1e99;
		m_collisionImpulse.setZero();
	}

	void ImplicitStepper::linearize()
	{
		m_strand.setFutureDegreesOfFreedom(m_strand.getCurrentDegreesOfFreedom()); //< duplicate
		m_dynamics.computeViscousForceCoefficients(m_dt);

		// compute A-M
		m_dynamics.computeFutureJacobian(true, m_params.m_energyWithTwist, m_params.m_energyWithBend, false, false);
		m_A = m_strand.getTotalJacobian();
		m_A *= m_dt * m_dt;

		//// compute b
		//if (!m_iteration) {	//< first iteration
		//	m_b = m_savedVelocities;
		//	m_dynamics.multiplyByMassMatrix(m_b);
		//	m_dynamics.computeFutureForces(true, true, true, false, false);
		//	m_b += m_strand.getFutureTotalForces() * m_dt;
		//}
		//else {
		//	m_b = m_b_hat;
		//	m_A.multiply(m_b, 1., m_velocities);
		//}
		m_b = m_savedVelocities - m_velocities;
		m_dynamics.multiplyByMassMatrix(m_b);
		m_dynamics.computeFutureForces(true, m_params.m_energyWithTwist, m_params.m_energyWithBend, false, false);
		m_b += m_strand.getFutureTotalForces() * m_dt;

		m_dynamics.addMassMatrixTo(m_A);
		m_A.multiply(m_b, 1., m_velocities);

		m_dynamics.getScriptingController()->fixRHS(m_mass, 1.0);
		m_dynamics.getScriptingController()->fixLHSAndRHS(m_A, m_b, m_dt);

		if (m_linearSolverType == LinearSolverType::DIRECT) {
			m_directSolver.store(m_A);
			m_notSPD = m_directSolver.notSPD();
		}

		m_strand.getFutureState().freeCachedQuantities();

		// reset linear solver
		m_substepVelocities = m_iteration ? m_velocities : m_savedVelocities;
		++m_iteration;
	}


	bool ImplicitStepper::solveLinear()
	{
		VecXx b = m_b;

		b += m_collisionImpulse;
		m_dynamics.getScriptingController()->fixRHS(m_A, b, m_dt);

		switch (m_linearSolverType)
		{
		case strandsim::ImplicitStepper::DIRECT:
			directSolver(b);
			break;
		case strandsim::ImplicitStepper::JACOBI:
			JacobiStep(b);
			break;
		case strandsim::ImplicitStepper::GAUSS_SEIDEL:
			GaussSeidelStep(b);
			break;
		case strandsim::ImplicitStepper::CONJ_GRAD:
			ConjgradStep(b);
			break;
		default:
			break;
		}

		VecXx Av = VecXx::Zero(m_velocities.size());
		m_A.multiply(Av, 1., m_newVelocities);
		m_residual = b - Av;
		m_bestError = std::min(m_bestError, m_residual.squaredNorm());

		m_substepVelocities = m_newVelocities;

		m_velocityDiff = (m_newVelocities - m_velocities).norm();
		return m_velocityDiff < m_tolerance;
	}

	void ImplicitStepper::postSolveLinear()
	{
		VecXx displacements = m_newVelocities * m_dt;
		m_dynamics.getScriptingController()->enforceDisplacements(displacements);
		m_strand.setCurrentDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom() + displacements);

		m_dynamics.getDisplacements() = displacements;
		m_dynamics.getAccelerations() = (m_newVelocities - m_velocities) / m_dt;

		if (m_params.m_linearizebHat)
		{
			m_b_hat = m_b;
			JacobianMatrixType A = m_strand.getTotalJacobian();
			A.multiply(m_b_hat, -m_dt * m_dt, m_newVelocities);
		}
		else
		{
			m_strand.setFutureDegreesOfFreedom(m_strand.getCurrentDegreesOfFreedom());
			m_b_hat = m_savedVelocities;
			m_dynamics.multiplyByMassMatrix(m_b_hat);
			m_dynamics.computeFutureForces(true, m_params.m_energyWithTwist, m_params.m_energyWithBend, false, false);
			m_b_hat += m_strand.getFutureTotalForces() * m_dt;
		}

		m_velocities = m_newVelocities;
	}

	void ImplicitStepper::accumulateCollision(int vid, const Vec7x& impulse)
	{
		omp_set_lock(&m_lock);
		m_collisionImpulse.segment<7>(4 * vid) += impulse;
		omp_unset_lock(&m_lock);
	}

	void ImplicitStepper::directSolver(const VecXx& b)
	{
		m_directSolver.solve(m_newVelocities, b);
	}

	void ImplicitStepper::JacobiStep(const VecXx& b)
	{
		Scalar omega = m_params.m_relaxationFactor;

		VecXx res = VecXx::Zero(m_velocities.size());
		m_A.multiply(res, 1., m_substepVelocities);
		res = b - res;

		m_newVelocities = (res.array() / m_A.diagonal().array() / omega).matrix() + m_substepVelocities;
		//m_newVelocities = (res.array() / m_mass.array() / omega).matrix() + m_velocities;

	}

	void ImplicitStepper::GaussSeidelStep(const VecXx& b)
	{
		const JacobianMatrixType A = m_A;

		m_newVelocities = m_substepVelocities;
		for (int i = 0; i < m_newVelocities.size(); ++i) {
			Scalar sum = - A(i, i) * m_newVelocities(i);
			for (int j = 0; j < m_newVelocities.size(); ++j)
				sum += A(i, j) * m_newVelocities(j);
			m_newVelocities(i) = (b(i) - sum) / A(i, i);
		}
	}

	void ImplicitStepper::ConjgradStep(const VecXx &b)
	{

	}

	Scalar ImplicitStepper::maxCollisionImpulseNorm(int& idx) const
	{
		Scalar maxlen = 0.;
		idx = -1;
		const int num_verts = m_strand.getNumVertices();
		for (int i = 0; i < num_verts; ++i)
		{
			const Scalar len = m_collisionImpulse.segment<3>(i * 4).norm();
			if (len > maxlen) {
				maxlen = len;
				idx = i;
			}
		}
		return maxlen;
	}

}