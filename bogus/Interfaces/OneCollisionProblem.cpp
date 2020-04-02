#include "OneCollisionProblem.hpp"
#include "../Extra/SecondOrder.impl.hpp"

namespace bogus
{
	template <int ndof> 
	OneCollisionProblem<ndof>::OneCollisionProblem( const MassMat& M, const DefGradMat& H,
		const ForceVec& f, const Eigen::Vector3d& uf, const Eigen::Matrix3d& E, double mu ):
		m_mu(mu), m_defGrad(E.transpose() * H), m_force(f)
	{
		m_MInv.compute(M);

		m_W = m_defGrad * (m_MInv * m_defGrad.transpose());
		m_b = E.transpose() * uf + m_defGrad * (m_MInv * f);
	}

	template <int ndof>
	OneCollisionProblem<ndof>::~OneCollisionProblem()
	{

	}

	template <int ndof>
	double OneCollisionProblem<ndof>::solve(Eigen::Vector3d& r, ForceVec& v_world, ForceVec& r_world)
	{
		CoulombLaw law(1, &m_mu);
		
		if (!law.solveLocal(0, m_W, m_b, r, 1.0)) {
			r *= 0.5;
		}
		
		r_world = m_defGrad.transpose() * r;
		v_world = m_MInv * (m_force + r_world);

		Eigen::Vector3d v = m_W * r + m_b;
		Eigen::Vector3d s;
		law.dualityCOV(0, v, s);

		return law.eval(0, r, v + s, 1.0);
	}

	template class OneCollisionProblem<7>;
	template class OneCollisionProblem<14>;
}