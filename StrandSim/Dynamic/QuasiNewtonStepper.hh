#ifndef _QUASI_NEWTON_STEPPER_HH_
#define _QUASI_NEWTON_STEPPER_HH_

#include "../Core/Definitions.hh"
#include "NewtonStepper.hh"
#include "../Utils/SymmetricBandMatrixSolver.hh"

#include <deque>

namespace strandsim
{
	class ElasticStrand;
	class StrandDynamicTraits;
	struct SimulationParameters;

	class QuasiNewtonStepper : public NewtonStepper
	{
	public:
		QuasiNewtonStepper(ElasticStrand& strand, const SimulationParameters& params);
		~QuasiNewtonStepper();

		virtual void prepareStep(Scalar dt);
		virtual bool performOneIteration();

	protected:
		// parameters
		int m_windowSize;
		std::deque<VecXx> m_v_queue;
		std::deque<VecXx> m_g_queue;

		JacobianSolver m_directSolver;

		VecXx m_last_v;
		VecXx m_last_gradient;
		bool m_isFirstIter;
	};
}

#endif // !_QUASI_NEWTON_STEPPER_HH_