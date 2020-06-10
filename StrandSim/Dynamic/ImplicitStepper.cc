#include "ImplicitStepper.hh"
#include "../Core/ElasticStrand.hh"
#include "SimulationParameters.hh"
#include "DOFScriptingController.hh"
#include "StrandDynamicTraits.hh"

namespace strandsim
{
	ImplicitStepper::ImplicitStepper(ElasticStrand& strand, const SimulationParameters& params) :
		m_notSPD(false),
		m_strand(strand),
		m_params(params),
		m_velocities(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
		m_collisionImpulse(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows())),
		m_timer(Timer("stepper", false))
	{

	}


	ImplicitStepper::~ImplicitStepper()
	{

	}

	void ImplicitStepper::accumulateCollisionImpulse(int vid, const Vec3x& r)
	{
		//std::lock_guard<std::mutex> gaurd(m_lock);
		m_collisionImpulse.segment<3>(4 * vid) += r;
	}

	
	void ImplicitStepper::commitVelocity() {
		m_strand.setCurrentDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom() + m_velocities * m_dt); 
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

	void ImplicitStepper::outputTiming() const
	{
		InfoStream(g_log, "Stepper breakdown timing")
			<< "hessian: " << m_timing.hessian
			<< "  fac: " << m_timing.factorize
			<< "  graident: " << m_timing.gradient
			<< "  solve: " << m_timing.solveLinear
			<< "  ls: " << m_timing.lineSearch
			<< "  total: " << m_timing.sum();
	}

}