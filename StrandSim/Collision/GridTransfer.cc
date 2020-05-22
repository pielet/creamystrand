#include "GridTransfer.hh"
#include "../Core/ElasticStrand.hh"
#include "../Dynamic/ImplicitStepper.hh"

namespace strandsim
{
	GridTransfer::GridTransfer(Scalar gridSize, int gridNum = 10) : m_gridSize(gridSize), m_invGridSize(1. / gridSize)
	{
		m_grids.reserve(gridNum);
	}

	GridTransfer::GridTransfer(const std::vector<ElasticStrand*>& strand, const std::vector<ImplicitStepper*>& stepper, Scalar gridSize, int gridNum = 10) :
		m_gridSize(gridSize), m_invGridSize(1. / gridSize)
	{
		m_grids.reserve(gridNum);
		buildGrid(strand, stepper);
	}

	GridTransfer::~GridTransfer()
	{

	}

	void GridTransfer::buildGrid(const std::vector<ElasticStrand*>& strand, const std::vector<ImplicitStepper*>& stepper)
	{
		for (int sidx = 0; sidx < strand.size(); ++sidx)
		{
			for (int vidx = 0; vidx < strand[sidx]->getNumVertices(); ++vidx)
			{
				const Vec3x pos = strand[sidx]->getVertex(vidx);
				const Vec3x vel = stepper[sidx]->getVelocity(vidx);

				const int min_i = std::floor(pos(0)) - 1;
				const int min_j = std::floor(pos(1)) - 1;
				const int min_k = std::floor(pos(2)) - 1;

				for (int i = min_i; i < min_i + 4; ++i) {
					for (int j = min_j; j < min_j + 4; ++j) {
						for (int k = min_k; k < min_k + 4; ++k) {
							auto grid_iter = m_grids.find(Grid(i, j, k));
							if (grid_iter != m_grids.end()) {
								grid_iter->second.second += getWeight(i, j, k, pos);
								grid_iter->second.first += getWeight(i, j, k, pos) * vel;
							}
							else {
								GridMap::mapped_type g_value;
								g_value.second = getWeight(i, j, k, pos);
								g_value.first = g_value.second * vel;
								m_grids.insert(GridMap::value_type(Grid(i, j, k), g_value));
							}
						}
					}
				}
			}
		}
		
		for (auto iter = m_grids.begin(); iter != m_grids.end(); ++iter)
		{
			iter->second.first /= iter->second.second;
		}
	}

	Vec3x GridTransfer::getValue(const Vec3x& pos)
	{
		Vec3x vel = Vec3x::Zero();

		const int min_i = std::floor(pos(0)) - 1;
		const int min_j = std::floor(pos(1)) - 1;
		const int min_k = std::floor(pos(2)) - 1;

		for (int i = min_i; i < min_i + 4; ++i) {
			for (int j = min_j; j < min_j + 4; ++j) {
				for (int k = min_k; k < min_k + 4; ++k) {
					auto grid_iter = m_grids.find(Grid(i, j, k));
					if (grid_iter != m_grids.end()) {
						vel += getWeight(i, j, k, pos) * grid_iter->second.first;
					}
				}
			}
		}

		return vel;
	}

	Scalar GridTransfer::getWeight(int i, int j, int k, const Vec3x& pos)
	{
		Vec3x grid_pos = Vec3x(i, j, k);
		grid_pos *= m_gridSize;

		Scalar weight = 1.;
		for (int i = 0; i < 3; ++i)
		{
			weight *= quadraticKernel((pos(i) - grid_pos(i)) * m_invGridSize);
		}

		return weight;
	}

	Scalar GridTransfer::quadraticKernel(Scalar x_in)
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