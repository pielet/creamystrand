#ifndef BOGUS_ONE_COLLISION_PROBLEM_H
#define BOGUS_ONE_COLLISION_PROBLEM_H

#include "../Core/Definitions.hh"
#include "ProximityCollision.hh"

namespace strandsim
{
	class CollisionSolver
	{
	public:
		CollisionSolver(ProximityCollision& col);

		~CollisionSolver();

		//! Solve the dual friction problem and update impulse
		virtual Scalar solve(Vec3x& r, Vec7x& r_world) { return 0.0; };
		virtual Scalar solve(Vec3x& r, Vec14x& r_world) { return 0.0; };

		//! update individual collision impulse
		virtual void local2World(Vec7x& r_world) {};
		virtual void local2World(Vec14x& r_world) {};

	protected:
		ProximityCollision& m_collision;
		double m_mu;
	};

	class ExternalContactSolver : public CollisionSolver
	{
	public:
		//! Input Paramsters of the primal problem
		/* Mv - f = H^T r
		 * u = Hv + u_f
		 * M: general mass matrix
		 * H: deformation gradient
		 * f: general force
		 * u_f: free velocity for external objects
		 * E: rotation matrix
		 */
		ExternalContactSolver( Mat7x M, Mat3x7x H, Vec3x uf, Mat3x E, ProximityCollision& col);

		virtual Scalar solve(Vec3x& r, Vec7x& r_world);
		virtual void local2World(Vec7x& r_world) { r_world = m_Hinv * m_collision.force; }

	protected:
		Vec7x m_M;
		Mat3x7x m_H;
		Mat7x3x m_Hinv;

		// Dual parameters
		Eigen::Vector3d m_b;
		Eigen::Matrix3d m_W;
	};

	class MutualContactSolver :public CollisionSolver
	{
	public:
		MutualContactSolver( Mat14x M, Mat3x14x H, Vec3x uf, Mat3x E, ProximityCollision& col );

		virtual Scalar solve(Vec3x& r, Vec14x& r_world);
		virtual void local2World(Vec14x& r_world) { r_world = m_Hinv * m_collision.force; }

	protected:
		Vec14x m_M;
		Mat3x14x m_H;
		Mat14x3x m_Hinv;

		// Dual parameters
		Eigen::Vector3d m_b;
		Eigen::Matrix3d m_W;
	};
}

#endif // !BOGUS_ONE_COLLISION_PROBLEM_H
