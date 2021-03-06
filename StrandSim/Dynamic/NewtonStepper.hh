#ifndef _NEWTON_STEPPER_HH_
#define _NEWTON_STEPPER_HH_

#include "../Core/Definitions.hh"
#include "ImplicitStepper.hh"

namespace strandsim
{
	class ElasticStrand;
	class StrandDynamicTraits;
	struct SimulationParameters;

	class NewtonStepper :public ImplicitStepper
	{
	public:
		NewtonStepper(ElasticStrand& strand, const SimulationParameters& params);
		virtual ~NewtonStepper();

		virtual void prepareStep(Scalar dt);
		virtual bool performOneIteration();
		virtual void postStep();

		virtual void rewind();
		virtual void updateVelocities();

		virtual Scalar stepSize() { return m_alpha; }

	protected:
		Scalar evaluateObjectValue(const VecXx& v);
		Scalar lineSearch(const VecXx& current_v, const VecXx& gradient_dir, const VecXx& descent_dir);
		
		StrandDynamicTraits& m_dynamics;

		VecXx m_savedVelocities;
		VecXx m_prevVelocities;

		int m_iteration;
		Scalar m_alpha;
		Scalar m_stepSize;
	};
}

#endif // !_NEWTON_STEPPER_HH_
