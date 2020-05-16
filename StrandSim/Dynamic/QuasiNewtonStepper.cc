#include "QuasiNewtonStepper.hh"
#include "../Core/ElasticStrand.hh"
#include "SimulationParameters.hh"
#include "StrandDynamicTraits.hh"
#include "DOFScriptingController.hh"
#include "../Utils/SymmetricBandMatrixSolver.hh"

namespace strandsim
{
	QuasiNewtonStepper::QuasiNewtonStepper(ElasticStrand& strand, const SimulationParameters& params):
		NewtonStepper(strand, params), m_windowSize(params.m_windowSize)
	{
		int ndof = strand.getCurrentDegreesOfFreedom().rows();

		m_last_v = VecXx::Zero(ndof);
		m_last_gradient = VecXx::Zero(ndof);
	}

	QuasiNewtonStepper::~QuasiNewtonStepper()
	{

	}

	void QuasiNewtonStepper::prepareStep(Scalar dt)
	{
		m_dt = dt;

		m_v_queue.clear();
		m_g_queue.clear();

		m_last_v.setZero();
		m_last_gradient.setZero();

		m_savedVelocities = m_velocities;
		m_dynamics.getScriptingController()->enforceVelocities(m_velocities, dt);

		m_isFirstIter = true;

		// compute initial hessian
		m_strand.setSavedDegreesOfFreedom(m_strand.getCurrentDegreesOfFreedom());
		m_strand.setFutureDegreesOfFreedom(m_strand.getCurrentDegreesOfFreedom() + m_velocities * m_dt);
		m_dynamics.computeFutureJacobian();
		JacobianMatrixType hessian = m_strand.getFutureTotalJacobian();
		hessian *= m_dt * m_dt;
		m_dynamics.addMassMatrixTo(hessian);
		m_dynamics.getScriptingController()->fixLHS(hessian);

		m_directSolver.store(hessian);	// pre-fab
	}

	bool QuasiNewtonStepper::performOneIteration()
	{
		m_strand.setFutureDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom() + m_velocities * m_dt);

		VecXx gradient = m_velocities - m_savedVelocities;
		m_dynamics.multiplyByMassMatrix(gradient);
		m_dynamics.computeFutureForces();
		gradient -= m_dt * m_strand.getFutureTotalForces();
		m_dynamics.getScriptingController()->fixRHS(gradient);

		if (gradient.squaredNorm() < square(SMALL_NUMBER<Scalar>())) return true;

		// debug
		//m_dynamics.computeFutureJacobian();
		//JacobianMatrixType hessian = m_strand.getFutureTotalJacobian();
		//hessian *= m_dt * m_dt;
		//m_dynamics.addMassMatrixTo(hessian);
		//m_dynamics.getScriptingController()->fixLHS(hessian);
		//JacobianSolver solver(hessian);
		//VecXx real_descent_dir = VecXx::Zero(m_velocities.size());
		//solver.solve(real_descent_dir, gradient);
		//real_descent_dir = -real_descent_dir;

		// update queue
		if (!m_isFirstIter)
		{
			m_v_queue.push_back(m_velocities - m_last_v);
			m_g_queue.push_back(gradient - m_last_gradient);

			if (m_v_queue.size() > m_windowSize)
			{
				m_v_queue.pop_front();
				m_g_queue.pop_front();
			}
		}

		m_isFirstIter = false;

		m_last_v = m_velocities;
		m_last_gradient = gradient;

		int w = std::min((int)m_v_queue.size(), m_windowSize);
		std::vector<Scalar> pho(w), alpha(w);
		// loop 1
		for (int i = 0; i < w; ++i)
		{
			VecXx v_i = m_v_queue[m_v_queue.size() - i - 1];
			VecXx g_i = m_g_queue[m_g_queue.size() - i - 1];
			Scalar vi_dot_gi = v_i.dot(g_i);
			if (vi_dot_gi < square(SMALL_NUMBER<Scalar>()))
			{
				return false;
			}
			pho[i] = 1. / vi_dot_gi;
			alpha[i] = v_i.dot(gradient) * pho[i];
			gradient -= alpha[i] * g_i;
		}
		// solve linear equation
		VecXx descent_dir = VecXx::Zero(m_velocities.size());
		m_directSolver.solve(descent_dir, gradient);
		// loop 2
		for (int i = w - 1; i >= 0; --i)
		{
			VecXx v_i = m_v_queue[m_v_queue.size() - i - 1];
			VecXx g_i = m_g_queue[m_g_queue.size() - i - 1];
			descent_dir += v_i * (alpha[i] - g_i.dot(descent_dir) * pho[i]);
		}
		descent_dir = -descent_dir;

		// line search
		Scalar step_size = lineSearch(m_velocities, gradient, descent_dir);
		m_velocities += step_size * descent_dir;
		m_dynamics.getScriptingController()->enforceVelocities(m_velocities, m_dt);

		// debug
		//std::cout << descent_dir.dot(real_descent_dir) / descent_dir.norm() / real_descent_dir.norm() << std::endl;
		
		m_strand.setCurrentDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom() + m_velocities * m_dt);

		m_strand.getFutureState().freeCachedQuantities();

		if (descent_dir.dot(-gradient) < square(SMALL_NUMBER<Scalar>()))
			return true;
		else
			return false;
	}

}