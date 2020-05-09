#include "CollisionSolver.hh"
#include "ProximityCollision.hh"
#include "../Dynamic/ImplicitStepper.hh"
#include <Eigen/LU>

namespace strandsim
{

	CollisionSolver::CollisionSolver(ProximityCollision& col)
		:m_collision(col), m_mu(col.mu)
	{

	}

	CollisionSolver::~CollisionSolver()
	{

	}

	ExternalContactSolver::ExternalContactSolver(Mat7x M, Mat3x7x H, Vec3x uf, Mat3x E, ProximityCollision& col):
		CollisionSolver(col), m_M(M.diagonal())
	{
		auto Minv_lu = M.lu();
		m_H = E.transpose() * H;
		m_Hinv = H.transpose() * (H * H.transpose()).lu().solve(E);
		m_W = m_H * Minv_lu.solve(m_Hinv);
		m_b = E.transpose() * uf;
	}

	Scalar ExternalContactSolver::solve(Vec3x& r, Vec7x& world_r)
	{
		Vec3x old_r = r;
		Vec3x b = m_b + m_H * world_r.cwiseQuotient(m_M) - m_W * old_r;

		if (b(0) > 0) {
			r.setZero();
		}
		else {
			double f_n = -b(0);
			double f_t = sqrt(b(1) * b(1) + b(2) * b(2));

			r = -b;

			if (f_t > m_mu* f_n) {
				r(1) *= m_mu * f_n / f_t;
				r(2) *= m_mu * f_n / f_t;
			}
		}
		r = m_W.lu().solve(r);

		world_r = m_Hinv * (r - old_r);

		return (r - old_r).norm();
	}


	MutualContactSolver::MutualContactSolver(Mat14x M, Mat3x14x H, Vec3x uf, Mat3x E, ProximityCollision& col) :
		CollisionSolver(col),
		m_M(M.diagonal())
	{
		auto Minv_lu = M.lu();
		m_H = E.transpose() * H;
		m_Hinv = H.transpose() * (H * H.transpose()).lu().solve(E);
		m_W = m_H * Minv_lu.solve(m_Hinv);
		m_b = E.transpose() * uf;
	}

	Scalar MutualContactSolver::solve(Vec3x& r, Vec14x& r_world)
	{
		Vec3x old_r = r;
		Vec3x b = m_b + m_H * r_world.cwiseQuotient(m_M) - m_W * old_r;

		if (b(0) > 0) {
			r.setZero();
		}
		else {
			double f_n = -b(0);
			double f_t = sqrt(b(1) * b(1) + b(2) * b(2));

			r = -b;

			if (f_t > m_mu* f_n) {
				r(1) *= m_mu * f_n / f_t;
				r(2) *= m_mu * f_n / f_t;
			}
		}
		r = m_W.lu().solve(r);

		r_world = m_Hinv * (r - old_r);

		return (r - old_r).norm();
	}

}