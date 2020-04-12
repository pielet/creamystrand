#include "LinearStepper.hh"
#include "../Core/ElasticStrand.hh"
#include "SimulationParameters.hh"
#include "DOFScriptingController.hh"
#include "StrandDynamicTraits.hh"

namespace strandsim
{
	LinearStepper::LinearStepper(ElasticStrand& strand, const SimulationParameters& params) :
		m_dt(0.),
		m_strand(strand),
		m_params(params),
		m_dynamics(strand.dynamics()),
		m_solverType(params.m_solverType),
		m_iteration(0),
		m_tolerance(params.m_linearSolverTolerance),
		m_notSPD(false),
		m_velocities(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
		m_newVelocities(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
		m_collisionImpulse(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
		m_collisionVelocity(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
		m_resWithoutImpulse(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
		m_residual(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
		m_bestError(1e+99)
	{
		m_dynamics.computeDOFMasses();
		m_mass = m_dynamics.getDOFMasses();
	}


	LinearStepper::~LinearStepper()
	{

	}


	void LinearStepper::initSolver(Scalar dt)
	{
		m_dt = dt;

		m_strand.setFutureDegreesOfFreedom(m_strand.getCurrentDegreesOfFreedom());
		m_dynamics.computeViscousForceCoefficients(dt);
		m_dynamics.getDisplacements().setZero();
		m_dynamics.getAccelerations().setZero();

		// compute b
		m_b = m_velocities;
		m_dynamics.multiplyByMassMatrix(m_b);
		m_dynamics.computeFutureForces(true, true, true, false, false);
		m_b += m_strand.getFutureTotalForces() * dt;

		// compute A
		m_dynamics.computeFutureJacobian(true, true, true, false, false);
		m_A = m_strand.getTotalJacobian();
		m_A *= dt * dt;
		m_dynamics.addMassMatrixTo(m_A);

		m_dynamics.getScriptingController()->fixRHS(m_mass, 1.0);
		m_dynamics.getScriptingController()->fixLHSAndRHS(m_A, m_b, dt);

		if (m_solverType == SolverType::DIRECT) {
			m_directSolver.store(m_A);
			m_notSPD = m_directSolver.notSPD();
		}

		m_strand.getFutureState().freeCachedQuantities();

		// reset solver
		m_iteration = 0;
		m_bestError = 1e+99;
		m_collisionImpulse.setZero();
		m_collisionVelocity.setZero();
	}


	bool LinearStepper::solveLinear()
	{
		VecXx b = m_b;
		if (m_params.m_registerImpulse) {
			b += m_collisionImpulse;
			m_dynamics.getScriptingController()->fixRHS(m_A, b, m_dt);
			m_collisionImpulse.setZero();
		}

		switch (m_solverType)
		{
		case strandsim::LinearStepper::DIRECT:
			directSolver(b);
			break;
		case strandsim::LinearStepper::JACOBI:
			JacobiStep(b);
			break;
		case strandsim::LinearStepper::GAUSS_SEIDEL:
			GaussSeidelStep(b);
			break;
		case strandsim::LinearStepper::CONJ_GRAD:
			ConjgradStep(b);
			break;
		default:
			break;
		}

		VecXx Av = VecXx::Zero(m_velocities.size());
		m_A.multiply(Av, 1., m_newVelocities);
		m_residual = b - Av;
		m_resWithoutImpulse = m_b - Av;
		m_bestError = std::min(m_bestError, m_residual.squaredNorm());

		Av.setZero();
		m_strand.getTotalJacobian().multiply(Av, 1., m_velocities);
		m_b_hat = m_b - m_dt * m_dt * Av;

		updateCurrentState();

		++m_iteration;

		return m_bestError < m_tolerance;
	}

	void LinearStepper::directSolver(const VecXx& b)
	{
		m_directSolver.solve(m_newVelocities, b);
	}

	void LinearStepper::JacobiStep(const VecXx& b)
	{
		Scalar omega = m_params.m_relaxationFactor;

		VecXx res = VecXx::Zero(m_velocities.size());
		m_A.multiply(res, 1., m_velocities);
		res = b - res;

		m_newVelocities = (res.array() / m_A.diagonal().array() / omega).matrix() + m_velocities;
		//m_newVelocities = (res.array() / m_mass.array() / omega).matrix() + m_velocities;

	}

	void LinearStepper::GaussSeidelStep(const VecXx& b)
	{
		const JacobianMatrixType A = m_A;

		m_newVelocities = m_velocities;
		for (int i = 0; i < m_newVelocities.size(); ++i) {
			Scalar sum = - A(i, i) * m_newVelocities(i);
			for (int j = 0; j < m_newVelocities.size(); ++j)
				sum += A(i, j) * m_newVelocities(j);
			m_newVelocities(i) = (b(i) - sum) / A(i, i);
		}
	}

	void LinearStepper::ConjgradStep(const VecXx &b)
	{

	}

	void LinearStepper::accumulateCollision(int vid, const Vec7x& vel, const Vec7x& impulse)
	{
		m_collisionImpulse.segment<7>(4 * vid) += impulse;
		m_collisionVelocity.segment<7>(4 * vid) += vel;
	}

	void LinearStepper::addCollisionVelocity()
	{
		m_newVelocities = (m_collisionVelocity.array().abs() < 1e-12)
			.select(m_newVelocities, m_collisionVelocity);
		m_collisionVelocity.setZero();

		updateCurrentState();
	}

	void LinearStepper::updateCurrentState()
	{
		VecXx displacements = m_newVelocities * m_dt;
		m_dynamics.getScriptingController()->enforceDisplacements(displacements);
		m_strand.setCurrentDegreesOfFreedom(m_strand.getFutureDegreesOfFreedom() + displacements);

		m_dynamics.getDisplacements() = displacements;
		m_dynamics.getAccelerations() = (m_newVelocities - m_velocities) / m_dt;
		m_velocities = m_newVelocities;
	}

	Scalar LinearStepper::maxCollisionImpulseNorm(int& idx) const
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