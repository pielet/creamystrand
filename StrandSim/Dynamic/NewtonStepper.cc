#include "NewtonStepper.hh"
#include "../Core/ElasticStrand.hh"
#include "SimulationParameters.hh"
#include "DOFScriptingController.hh"
#include "StrandDynamicTraits.hh"
#include "../Utils/SymmetricBandMatrixSolver.hh"

namespace strandsim
{

	NewtonStepper::NewtonStepper(ElasticStrand& strand, const SimulationParameters& params) :
		m_strand(strand), m_dynamics(strand.dynamics()), m_params(params)
	{
		int ndof = strand.getCurrentDegreesOfFreedom().rows();
		m_velocities = VecXx::Zero(ndof);
		m_savedVelocities = VecXx::Zero(ndof);
	}

	NewtonStepper::~NewtonStepper()
	{

	}

	void NewtonStepper::prepareStep(Scalar dt)
	{
		m_dt = dt;

		m_savedVelocities = m_velocities;
		// m_velocities.setZero(); // v_{t+1}^0 = v_t

		m_strand.setSavedDegreesOfFreedom(m_strand.getCurrentDegreesOfFreedom());
	}

	bool NewtonStepper::performOneIteration()
	{
		m_strand.setFutureDegreesOfFreedom(m_strand.getCurrentDegreesOfFreedom() + m_velocities * m_dt);

		// gradient
		VecXx gradient = m_velocities - m_savedVelocities;
		m_dynamics.multiplyByMassMatrix(gradient);
		m_dynamics.computeFutureForces();
		gradient -= m_strand.getFutureTotalForces() * m_dt;

		// hessian
		m_dynamics.computeFutureJacobian();
		JacobianMatrixType hessian = m_strand.getFutureTotalJacobian();
		hessian *= m_dt * m_dt;
		m_dynamics.addMassMatrixTo(hessian);

		// fix point
		m_dynamics.getScriptingController()->fixLHSAndRHS(hessian, gradient, m_dt);

		// solve linear equation
		VecXx descent_dir = VecXx::Zero(m_velocities.size());
		JacobianSolver directSolver;
		directSolver.store(hessian);
		directSolver.solve(descent_dir, gradient);
		descent_dir = -descent_dir;

		// line search
		Scalar step_size = lineSearch(m_velocities, gradient, descent_dir);
		m_velocities += step_size * descent_dir;

		m_strand.getFutureState().freeCachedQuantities();

		if (descent_dir.dot(-gradient) < square(SMALL_NUMBER<Scalar>()))
			return true;
		else
			return false;
	}

	void NewtonStepper::postStep()
	{
		VecXx displacements = m_velocities * m_dt;
		//m_dynamics.getScriptingController()->enforceDisplacements(displacements);
		m_strand.setCurrentDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom() + displacements);

		m_dynamics.getDisplacements() = displacements;
	}
	
	Scalar NewtonStepper::lineSearch(const VecXx& current_v, const VecXx& gradient_dir, const VecXx& descent_dir)
	{
		if (m_params.m_useLineSearch) {
			Scalar current_obj_value = evaluateObjectValue(current_v), next_obj_value;
			Scalar alpha = 1. / m_params.m_ls_beta;
			Scalar rhs;

			do {
				alpha *= m_params.m_ls_beta;
				if (alpha < 1e-5) break;

				next_obj_value = evaluateObjectValue(current_v + alpha * descent_dir);
				rhs = current_obj_value + m_params.m_ls_alpha * alpha * gradient_dir.dot(descent_dir);
			
			} while (next_obj_value > rhs);

			if (alpha < 1e-5)
				alpha = 0;

			return alpha;
		}
		else {
			return 1.0;
		}
	}

	Scalar NewtonStepper::evaluateObjectValue(const VecXx& v)
	{
		// internal energy
		m_strand.setFutureDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom() + v * m_dt);
		m_dynamics.computeFutureStrandEnergy();
		Scalar g = m_strand.getFutureTotalEnergy();

		// inertia term
		VecXx u_t = m_dynamics.getExternalForce();
		m_dynamics.multiplyByMassMatrixInverse(u_t);
		u_t = v - (m_savedVelocities + m_dt * u_t);
		VecXx delta_v = u_t;
		m_dynamics.multiplyByMassMatrix(u_t);

		g += 0.5 * delta_v.transpose() * u_t;

		return g;
	}
}
