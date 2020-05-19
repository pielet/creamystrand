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
		m_collisionImpulse(VecXx::Zero(strand.getCurrentDegreesOfFreedom().rows()))
	{
		omp_init_lock(&m_lock);
	}


	ImplicitStepper::~ImplicitStepper()
	{

	}

	void ImplicitStepper::rewind()
	{
		m_strand.setCurrentDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom() + m_velocities * m_dt);
		m_strand.setFutureDegreesOfFreedom(m_strand.getSavedDegreesOfFreedom());
	}

	void ImplicitStepper::accumulateCollisionImpulse(int vid, const Vec3x& r)
	{
		m_collisionImpulse.segment<3>(4 * vid) += r;
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