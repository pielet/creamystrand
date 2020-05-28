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
		NewtonStepper::prepareStep(dt);

		m_v_queue.clear();
		m_g_queue.clear();

		m_last_v.setZero();
		m_last_gradient.setZero();

		// compute initial hessian
		m_dynamics.getScriptingController()->enforceVelocities(m_velocities, dt);
		m_strand.setFutureDegreesOfFreedom(m_strand.getCurrentDegreesOfFreedom() + m_velocities * m_dt);
		m_timer.restart();
		m_dynamics.computeFutureJacobian();
		JacobianMatrixType hessian = m_strand.getFutureTotalJacobian();
		hessian *= m_dt * m_dt;
		m_dynamics.addMassMatrixTo(hessian);
		m_dynamics.getScriptingController()->fixLHS(hessian);
		m_timing.hessian += m_timer.elapsed();

		m_timer.restart();
		m_directSolver.store(hessian);	// pre-fab
		m_timing.factorize += m_timer.elapsed();
	}

	bool QuasiNewtonStepper::performOneIteration()
	{
		m_strand.setFutureDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom() + m_velocities * m_dt);

		// gradiant
		m_timer.restart();
		VecXx gradient = m_velocities - m_savedVelocities;
		m_dynamics.multiplyByMassMatrix(gradient);
		m_dynamics.computeFutureForces();
		gradient -= m_dt * m_strand.getFutureTotalForces() + m_collisionImpulse;
		m_dynamics.getScriptingController()->fixRHS(gradient);
		m_timing.gradient += m_timer.elapsed();

		// check convergence
		Scalar err = gradient.squaredNorm() / gradient.size();
		//if (isSmall(err) || (m_iteration > 3 && err < 1e-6))
		//	return true;

		// update saved info
		if (m_iteration)
		{
			// update queue
			m_v_queue.push_back(m_velocities - m_last_v);
			m_g_queue.push_back(gradient - m_last_gradient);

			if (m_v_queue.size() > m_windowSize)
			{
				m_v_queue.pop_front();
				m_g_queue.pop_front();
			}
		}

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

			if (isSmallSquare(vi_dot_gi)) continue;

			pho[i] = 1. / vi_dot_gi;
			alpha[i] = v_i.dot(gradient) * pho[i];
			gradient -= alpha[i] * g_i;
		}
		// solve linear equation
		VecXx descent_dir = VecXx::Zero(m_velocities.size());
		m_timer.restart();
		m_directSolver.solve(descent_dir, gradient);
		m_timing.solveLinear += m_timer.elapsed();
		// loop 2
		for (int i = w - 1; i >= 0; --i)
		{
			VecXx v_i = m_v_queue[m_v_queue.size() - i - 1];
			VecXx g_i = m_g_queue[m_g_queue.size() - i - 1];

			if (isSmallSquare(v_i.dot(g_i))) continue;

			descent_dir += v_i * (alpha[i] - g_i.dot(descent_dir) * pho[i]);
		}
		descent_dir = -descent_dir;

		// line search
		m_timer.restart();
		Scalar step_size = lineSearch(m_velocities, gradient, descent_dir);
		m_timing.lineSearch += m_timer.elapsed();
		m_velocities += step_size * descent_dir;
		m_dynamics.getScriptingController()->enforceVelocities(m_velocities, m_dt);
		resetCollisionVelocities();

		m_strand.setCurrentDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom() + m_velocities * m_dt);

		m_strand.getFutureState().freeCachedQuantities();

		++m_iteration;

		VecXx velDiff = m_velocities - m_prevVelocities;
		m_prevVelocities = m_velocities;
		m_dynamics.getScriptingController()->fixRHS(velDiff);

		if (velDiff.norm() / velDiff.size() < m_params.m_velocityDiffTolerance) return true;

		return false;
	}

}