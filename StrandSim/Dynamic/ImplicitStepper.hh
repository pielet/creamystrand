#ifndef STRANDSIM_ImplicitStepper_HH
#define STRANDSIM_ImplicitStepper_HH

#if defined(_OPENMP)
#include <omp.h>
#endif
#include "../Core/Definitions.hh"

namespace strandsim
{
	class ElasticStrand;
	class StrandDynamicTraits;
	struct SimulationParameters;

	class ImplicitStepper
	{
	public:
		enum LinearSolverType { DIRECT, JACOBI, GAUSS_SEIDEL, CONJ_GRAD };

		ImplicitStepper(ElasticStrand& strand, const SimulationParameters& params);
		virtual ~ImplicitStepper();

		virtual void prepareStep(Scalar dt) = 0;
		virtual bool performOneIteration() = 0;
		virtual void postStep() = 0;

		virtual void rewind() = 0;

		void accumulateCollisionImpulse(int vid, const Vec3x& r, Scalar decay);
		void clearCollisionImpulse() { m_collisionImpulse.setZero(); }
		Scalar maxCollisionImpulseNorm(int& idx) const;

		void finalize() {}

		bool refusesMutualContacts() const { return m_notSPD; }
		Scalar getStretchMultiplier() const { return 1.0; }
		VecXx flowComponents() const { return VecXx::Ones(m_velocities.size()); }	// for StrandRender::computeFlowQuads and ProblemStepper::dumpRods

		Vec3x getVelocity(int vid) { return m_velocities.segment<3>(4 * vid); }
		void setVelocity(int vid, const Vec3x& vel) { m_velocities.segment<3>(4 * vid) = vel; }
		VecXx& velocities() { return m_velocities; }
		const VecXx& velocities() const { return m_velocities; }

		const VecXx& getCollisionImpulse() const { return m_collisionImpulse; }

	protected:
		omp_lock_t m_lock;

		ElasticStrand& m_strand;
		const SimulationParameters& m_params;
		
		Scalar m_dt;

		VecXx m_velocities;
		VecXx m_collisionImpulse;

		bool m_notSPD;
	};
}

#endif // !STRANDSIM_ImplicitStepper_HH
