#ifndef BOGUS_ONE_COLLISION_PROBLEM_H
#define BOGUS_ONE_COLLISION_PROBLEM_H

#include <Eigen/Core>
#include "../Core/Eigen/EigenLinearSolvers.hpp"
#include "../Extra/SecondOrder.fwd.hpp"

namespace bogus
{
	template <int ndof>
	class OneCollisionProblem
	{
	public:
		typedef Eigen::Matrix<double, ndof, ndof> MassMat;
		typedef LU<Eigen::MatrixBase<MassMat> > MInvType;
		typedef Eigen::Matrix<double, 3, ndof> DefGradMat;
		typedef Eigen::Matrix<double, ndof, 1> ForceVec;

		typedef SOCLaw<3, double, true> CoulombLaw;

		//! Input Paramsters of the primal problem
		/* Mv - f = H^T r
		 * u = Hv + u_f
		 * M: general mass matrix
		 * H: deformation gradient
		 * f: general force
		 * u_f: free velocity for external objects
		 * E: rotation matrix
		 * mu: friction coeff
		 */
		OneCollisionProblem(
			const MassMat& M,
			const DefGradMat& H,
			const ForceVec& f,
			const Eigen::Vector3d& uf,
			const Eigen::Matrix3d& E,
			double mu
		);

		~OneCollisionProblem();

		//! Solve the dual friction problem
		double solve(Eigen::Vector3d& r, ForceVec& v, ForceVec& r_world);

	protected:
		double m_mu;

		// Primal parameters
		MInvType m_MInv;
		DefGradMat m_defGrad;
		ForceVec m_force;

		// Dual parameters
		Eigen::Vector3d m_b;
		Eigen::Matrix3d m_W;
	};

	typedef OneCollisionProblem<7> ExternalContactSolver;
	typedef OneCollisionProblem<14> MutualContactSolver;
}

#endif // !BOGUS_ONE_COLLISION_PROBLEM_H
