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

	protected:
		Scalar evaluateObjectValue(const VecXx& v);
		Scalar lineSearch(const VecXx& current_v, const VecXx& gradient_dir, const VecXx& descent_dir);
		
		//ElasticStrand& m_strand;
		StrandDynamicTraits& m_dynamics;
		//const SimulationParameters& m_params;

		Scalar m_dt;
		VecXx m_velocities;
		VecXx m_savedVelocities;

		int m_iteration;
		VecXx m_bestVelocities;
		Scalar m_bestErr;
		Scalar m_prevErr;
		Scalar m_alpha;
	};
}

#endif // !_NEWTON_STEPPER_HH_
