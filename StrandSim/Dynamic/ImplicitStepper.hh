#ifndef STRANDSIM_ImplicitStepper_HH
#define STRANDSIM_ImplicitStepper_HH

#if defined(_OPENMP)
#include <omp.h>
#endif
#include "../Core/Definitions.hh"
#include "../Utils/LoggingTimer.hh"
#include "../Utils/SymmetricBandMatrixSolver.hh"

namespace strandsim
{
	class ElasticStrand;
	class StrandDynamicTraits;
	struct SimulationParameters;

	struct StepperTiming
	{
		double hessian;
		double gradient;
		double factorize;
		double solveLinear;
		double lineSearch;

		StepperTiming() : hessian(0), gradient(0), factorize(0), solveLinear(0), lineSearch(0) {}
		StepperTiming(double h, double g, double f, double s, double l) :
			hessian(h), gradient(g), factorize(f), solveLinear(s), lineSearch(l) {}

		void reset()
		{
			hessian = 0;
			gradient = 0;
			factorize = 0;
			solveLinear = 0;
			lineSearch = 0;
		}
		double sum() const
		{
			return hessian + gradient + factorize + solveLinear + lineSearch;
		}
		const StepperTiming& operator+=(const StepperTiming& other)
		{
			hessian += other.hessian;
			gradient += other.gradient;
			factorize += other.factorize;
			solveLinear += other.factorize;
			lineSearch += other.lineSearch;
			return *this;
		}
		const StepperTiming& operator/=(int time)
		{
			hessian /= time;
			gradient /= time;
			factorize /= time;
			solveLinear /= time;
			lineSearch /= time;
			return *this;
		}
	};

	class ImplicitStepper
	{
	public:
		enum LinearSolverType {DIRECT, JACOBI, GAUSS_SEIDEL};

		ImplicitStepper(ElasticStrand& strand, const SimulationParameters& params);
		virtual ~ImplicitStepper();

		virtual void prepareStep(Scalar dt) = 0;
		virtual bool performOneIteration() = 0;
		virtual void postStep() = 0;

		virtual void rewind() = 0;

		virtual Scalar stepSize() = 0;

		const StepperTiming& getTiming() { return m_timing; }

		void outputTiming() const;

		//void accumulateCollisionImpulse(int vid, const Vec3x& r);
		void clearCollisionImpulse() { m_totalCollisionImpulse.setZero(); }
		Scalar maxCollisionImpulseNorm(int& idx) const;
		VecXx& totalCollisionImpulses() { return m_totalCollisionImpulse; }
		const VecXx& totalCollisionImpulses() const { return m_totalCollisionImpulse; }

		void commitVelocity();

		void accumulateDeltaCollisionImpulse(int vid, const Vec3x& vel);
		void clearDeltaCollisionImpulse() { m_deltaCollisionImpulse.setZero(); }
		const VecXx& deltaCollisionImpulse() const { return m_deltaCollisionImpulse; }

		void updateVelocities();

		Vec3x getVelocity(int vid) { return m_velocities.segment<3>(4 * vid); }
		void setVelocity(int vid, const Vec3x& vel) { m_velocities.segment<3>(4 * vid) = vel; }
		VecXx& velocities() { return m_velocities; }
		const VecXx& velocities() const { return m_velocities; }

		void finalize() {}

		ElasticStrand& getStrand() { return m_strand; }

		bool refusesMutualContacts() const { return m_notSPD; }
		Scalar getStretchMultiplier() const { return 1.0; }
		VecXx flowComponents() const { return VecXx::Ones(m_velocities.size()); }	// for StrandRender::computeFlowQuads and ProblemStepper::dumpRods

	protected:
#if defined(_OPENMP)
		omp_lock_t m_lock;
#endif

		ElasticStrand& m_strand;
		const SimulationParameters& m_params;
		
		Scalar m_dt;

		VecXx m_velocities;
		VecXx m_totalCollisionImpulse;
		VecXx m_deltaCollisionImpulse;

		JacobianSolver m_directSolver;

		StepperTiming m_timing;
		Timer m_timer;

		bool m_notSPD;
	};
}

#endif // !STRANDSIM_ImplicitStepper_HH
