#include "GridTransfer.hh"
#include "../Core/ElasticStrand.hh"
#include "../Dynamic/ImplicitStepper.hh"

namespace strandsim
{
	GridTransfer::GridTransfer(Scalar gridSize, int gridNum) : m_gridSize(gridSize), m_invGridSize(1. / gridSize)
	{
		m_grids.reserve(gridNum);
	}

	GridTransfer::~GridTransfer()
	{

	}

	void GridTransfer::insertParticle(const Vec3x& p, Scalar mass, const Vec3x& vel)
	{
		const Vec3x pos = p * m_invGridSize;

		const int min_i = std::floor(pos(0)) - 1;
		const int min_j = std::floor(pos(1)) - 1;
		const int min_k = std::floor(pos(2)) - 1;

		for (int i = min_i; i < min_i + 4; ++i) {
			for (int j = min_j; j < min_j + 4; ++j) {
				for (int k = min_k; k < min_k + 4; ++k) {
					Scalar weight = getWeight(i, j, k, pos);
					if (isSmall(weight)) continue;

					auto grid_iter = m_grids.find(Grid(i, j, k));
					if (grid_iter != m_grids.end()) {
						grid_iter->second.first += weight * mass;
						grid_iter->second.second += weight * mass * vel;
					}
					else {
						GridMap::mapped_type g_value;
						g_value.first = weight * mass;
						g_value.second = weight * mass * vel;
						m_grids.insert(GridMap::value_type(Grid(i, j, k), g_value));
					}
				}
			}
		}
	}

	Vec3x GridTransfer::getValue(const Vec3x& pos)
	{
		Vec3x vel = Vec3x::Zero();

		Vec3x p = pos * m_invGridSize;
		const int min_i = std::floor(p(0)) - 1;
		const int min_j = std::floor(p(1)) - 1;
		const int min_k = std::floor(p(2)) - 1;

		for (int i = min_i; i < min_i + 4; ++i) {
			for (int j = min_j; j < min_j + 4; ++j) {
				for (int k = min_k; k < min_k + 4; ++k) {
					auto grid_iter = m_grids.find(Grid(i, j, k));
					if (grid_iter != m_grids.end()) {
						vel += getWeight(i, j, k, p) * (grid_iter->second.second / grid_iter->second.first);
					}
				}
			}
		}

		return vel;
	}

	Scalar GridTransfer::getWeight(const Vec3x& p1, const Vec3x& p2) const
	{
		Scalar weight = 1.;
		for (int i = 0; i < 3; ++i) {
			weight *= quadraticKernel((p1(i) - p2(i)) * m_invGridSize);
		}
		return weight;
	}

	Scalar GridTransfer::getWeight(Scalar dis) const
	{
		return quadraticKernel(dis * m_invGridSize);
	}

	Scalar GridTransfer::getWeight(int i, int j, int k, const Vec3x& pos)
	{
		Vec3x grid_pos = Vec3x(i, j, k);

		Scalar weight = 1.;
		for (int i = 0; i < 3; ++i)
		{
			weight *= quadraticKernel(pos(i) - grid_pos(i));	// pos is normalized
		}

		return weight;
	}

	Scalar GridTransfer::quadraticKernel(Scalar x_in) const
	{
		Scalar x = fabs(x_in);

		if (x < 0.5)
			return 0.75 - x * x;
		else if (x < 1.5)
			return 0.5 * square(1.5 - x);
		else
			return 0;
	}
}