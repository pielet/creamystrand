#include "OneCollisionProblem.hpp"
#include "../Extra/SecondOrder.impl.hpp"
#include "../Core/Utils/NumTraits.hpp"
#include <Eigen/LU>

namespace bogus
{
	template <int ndof> 
	OneCollisionProblem<ndof>::OneCollisionProblem( const MassMat& M, const DefGradMat& H,
		const ForceVec& f, const Eigen::Vector3d& uf, const Eigen::Matrix3d& E, double mu ):
		m_mu(mu), m_defGrad(E.transpose() * H)
	{
		m_MInv.compute(M);
		m_defGradInv.transpose() = H.transpose() * (H * H.transpose()).lu().solve(E);

		m_W = m_defGrad * (m_MInv * m_defGradInv.transpose());
		m_b = E.transpose() * uf + m_defGrad * (m_MInv * f);
	}

	template <int ndof>
	OneCollisionProblem<ndof>::~OneCollisionProblem()
	{

	}

	template <int ndof>
	double OneCollisionProblem<ndof>::solve(Eigen::Vector3d& r, ForceVec& v_world, ForceVec& delta_r_world)
	{
		CoulombLaw law(1, &m_mu);

		//if (!law.solveLocal(0, m_W, m_b, r, 1.0)) {
		//	r *= 0.5;
		//}

		Eigen::Vector3d old_r = r;
		m_b -= m_W * old_r;

		if (m_b(0) > 0) {
			r.setZero();
		}
		else {
			double f_n = -m_b(0);
			double f_t = sqrt(m_b(1) * m_b(1) + m_b(2) * m_b(2));

			r = -m_b;
			
			if (f_t > m_mu * f_n) {
				r(1) *= m_mu * f_n / f_t;
				r(2) *= m_mu * f_n / f_t;
			}
		}
		r = m_W.lu().solve(r);

		delta_r_world = m_defGradInv.transpose() * (r - old_r);
		//delta_r_world = m_defGradInv.transpose() * r;
		
		Eigen::Vector3d v = m_W * r + m_b;
		Eigen::Vector3d s;
		law.dualityCOV(0, v, s);

		return law.eval(0, r, v + s, 1.0);
	}

	template class OneCollisionProblem<7>;
	template class OneCollisionProblem<14>;
}