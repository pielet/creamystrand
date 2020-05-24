#ifndef _GRID_TRANSFER_HH_
#define _GRID_TRANSFER_HH_

#include "../Core/Definitions.hh"
#include <unordered_map>

namespace strandsim
{
	class ElasticStrand;
	class ImplicitStepper;

	class GridTransfer
	{
	public:
		GridTransfer(Scalar girdSize, int gridNum = 10);
		GridTransfer(const std::vector<ElasticStrand*>& strand, const std::vector<ImplicitStepper*>& stepper, Scalar gridSize, int gridNum = 10);
		~GridTransfer();

		void buildGrid(const std::vector<ElasticStrand*>& strand, const std::vector<ImplicitStepper*>& stepper);
		Vec3x getValue(const Vec3x& pos);

	private:
		struct Grid
		{
			int m_i, m_j, m_k;
			long m_hash;

			Grid(int i, int j, int k) :m_i(i), m_j(j), m_k(k),
				m_hash(std::hash<int>{}(m_i) ^ std::hash<int>{}(m_j) ^ std::hash<int>{}(m_k)) {};

			bool operator==(const Grid& other) const
			{
				return m_i == other.m_i && m_j == other.m_j && m_k == other.m_k;
			}

			struct Hasher
			{
				long operator()(const Grid& g) const
				{
					return g.m_hash;
				}
			};

			struct Comparer
			{
				bool operator()(const Grid& lhs, const Grid& rhs) const
				{
					return lhs == rhs;
				}
			};
		};

		Scalar getWeight(int i, int j, int k, const Vec3x& pos);
		Scalar quadraticKernel(Scalar x);

		typedef std::unordered_map<Grid, std::pair<Vec3x, Scalar>, typename Grid::Hasher, typename Grid::Comparer> GridMap;
		GridMap m_grids;

		Scalar m_gridSize;
		Scalar m_invGridSize;
	};
}

#endif // _GRID_TRANSFER_HH_
