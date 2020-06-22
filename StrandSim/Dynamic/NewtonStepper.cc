#include "NewtonStepper.hh"
#include "../Core/ElasticStrand.hh"
#include "SimulationParameters.hh"
#include "DOFScriptingController.hh"
#include "StrandDynamicTraits.hh"
#include "../Utils/SymmetricBandMatrixSolver.hh"
#include "../Utils/MathUtilities.hh"

namespace strandsim
{

	NewtonStepper::NewtonStepper(ElasticStrand& strand, const SimulationParameters& params) :
		ImplicitStepper(strand, params), m_dynamics(strand.dynamics())
	{
		int ndof = strand.getCurrentDegreesOfFreedom().rows();
		m_velocities = VecXx::Zero(ndof);
		m_savedVelocities = VecXx::Zero(ndof);
		m_prevVelocities = VecXx::Zero(ndof);

		m_dynamics.computeDOFMasses();
	}

	NewtonStepper::~NewtonStepper()
	{

	}

	void NewtonStepper::prepareStep(Scalar dt)
	{
		m_dt = dt;

		m_iteration = 0;
		m_alpha = 1.0;

		m_timing.reset();

		clearCollisionImpulse();

		m_savedVelocities = m_velocities;
		m_prevVelocities = m_velocities;
		m_strand.setSavedDegreesOfFreedom(m_strand.getCurrentDegreesOfFreedom());

		m_dynamics.getScriptingController()->enforceVelocities(m_velocities, m_dt);

		m_dynamics.computeViscousForceCoefficients(dt);  // used for updating dt (air draging force);
	}

	bool NewtonStepper::performOneIteration()
	{
		m_strand.setFutureDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom() + m_velocities * m_dt);

		// gradient
		m_timer.restart();
		VecXx gradient = m_velocities - m_savedVelocities;
		m_dynamics.multiplyByMassMatrix(gradient);
		m_dynamics.computeFutureForces(m_params.m_energyWithStretch, m_params.m_energyWithTwist, m_params.m_energyWithBend);
		gradient -= m_strand.getFutureTotalForces() * m_dt + m_totalCollisionImpulse;
		m_dynamics.getScriptingController()->fixRHS(gradient);  // fix points
		m_timing.gradient += m_timer.elapsed();

		// hessian
		m_timer.restart();
		m_dynamics.computeFutureJacobian(m_params.m_energyWithStretch, m_params.m_energyWithTwist, m_params.m_energyWithBend);
		JacobianMatrixType hessian = m_strand.getFutureTotalJacobian();
		hessian *= m_dt * m_dt;
		m_dynamics.addMassMatrixTo(hessian);
		m_dynamics.getScriptingController()->fixLHS(hessian);  // fix points
		m_timing.hessian += m_timer.elapsed();

		// solve linear equation
		VecXx descent_dir = VecXx::Zero(m_velocities.size());
		m_timer.restart();
		m_directSolver.store(hessian);
		m_notSPD = m_directSolver.notSPD();
		m_timing.factorize += m_timer.elapsed();
		m_timer.restart();
		m_directSolver.solve(descent_dir, gradient);
		m_timing.solveLinear += m_timer.elapsed();
		descent_dir = -descent_dir;

		// line search
		m_timer.restart();
		Scalar step_size = lineSearch(m_velocities, gradient, descent_dir);
		m_timing.lineSearch += m_timer.elapsed();
		m_velocities += step_size * descent_dir;
		m_dynamics.getScriptingController()->enforceVelocities(m_velocities, m_dt);

		m_strand.setCurrentDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom() + m_velocities * m_dt);
		
		m_strand.getFutureState().freeCachedQuantities();

		++m_iteration;

		VecXx velDiff = m_velocities - m_prevVelocities;
		m_prevVelocities = m_velocities;
		m_dynamics.getScriptingController()->fixRHS(velDiff);
		
		if (velDiff.norm() / velDiff.size() < m_params.m_velocityDiffTolerance) return true;

		return false;
	}

	void NewtonStepper::postStep()
	{
		VecXx displacements = m_velocities * m_dt;
		m_dynamics.getScriptingController()->enforceDisplacements(displacements);
		m_strand.setCurrentDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom() + displacements);

		m_dynamics.getDisplacements() = displacements;
	}

	void NewtonStepper::rewind()
	{
		VecXx displacements = m_velocities * m_dt;
		m_dynamics.getScriptingController()->enforceDisplacements(displacements);
		m_dynamics.getDisplacements() = displacements;
		m_strand.setCurrentDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom() + displacements);
		m_strand.setFutureDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom());
	}

	void NewtonStepper::updateVelocities()
	{
		VecXx delta_v = VecXx::Zero(m_velocities.size());
		m_dynamics.getScriptingController()->fixRHS(m_deltaCollisionImpulse);
		m_directSolver.solve(delta_v, m_deltaCollisionImpulse);
		m_velocities += delta_v;
	}
	
	Scalar NewtonStepper::lineSearch(const VecXx& current_v, const VecXx& gradient_dir, const VecXx& descent_dir)
	{
		if (m_params.m_useLineSearch) {
			Scalar current_obj_value = evaluateObjectValue(current_v), next_obj_value;
			m_alpha = std::min(1., 2 * m_alpha) / m_params.m_ls_beta;
			Scalar rhs;

			do {
				m_alpha *= m_params.m_ls_beta;
				if (m_alpha < 1e-5) break;

				next_obj_value = evaluateObjectValue(current_v + m_alpha * descent_dir);
				rhs = current_obj_value + m_params.m_ls_alpha * m_alpha * gradient_dir.dot(descent_dir);
			} while (next_obj_value > rhs);

			if (m_alpha < 1e-5)
				m_alpha = 0;

			return m_alpha;
		}
		else {
			return 1.0;
		}
	}

	Scalar NewtonStepper::evaluateObjectValue(const VecXx& vel)
	{
		VecXx v = vel;
		m_dynamics.getScriptingController()->enforceVelocities(v, m_dt);

		// internal energy
		m_strand.setFutureDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom() + v * m_dt);
		m_dynamics.computeFutureStrandEnergy(m_params.m_energyWithStretch, m_params.m_energyWithTwist, m_params.m_energyWithBend);;
		Scalar g = m_strand.getFutureTotalEnergy();

		// inertia term
		VecXx u_t = m_dynamics.getExternalForce() * m_dt + m_totalCollisionImpulse;
		m_dynamics.multiplyByMassMatrixInverse(u_t);
		u_t = v - (m_savedVelocities + u_t);
		VecXx delta_v = u_t;
		m_dynamics.multiplyByMassMatrix(u_t);

		g += 0.5 * delta_v.transpose() * u_t;

		return g;
	}
}
