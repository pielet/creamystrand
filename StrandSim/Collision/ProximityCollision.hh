/**
 * \copyright 2013 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef STRANDSIM_PROXIMITYCOLLISION_HH
#define STRANDSIM_PROXIMITYCOLLISION_HH

#include <Eigen/Sparse>
#include "../Core/BandMatrixFwd.hh"
#include <map>
#include <unordered_map>


namespace strandsim
{
    class ElasticStrand;
    class ContinuousTimeCollision;
    
    struct ProximityCollision
    {
    public:
        struct Object
        {
            int globalIndex;
            int vertex;
            Scalar abscissa;

            Vec3x freeVel;
            Vec3x direction;
            Mat3x invInertia;	// compute by generateInverseInertia()

            Object() : globalIndex(0), vertex(0), abscissa(0) {}

            void generateInverseInertia(const Vec3x& normal);

            unsigned long computeSizeInBytes() const
            {
                unsigned long total = 0;
                total += sizeof(int) + sizeof(int) + sizeof(Scalar);
                total += freeVel.size() * sizeof(Vec3x::Scalar) * 2;
                total += invInertia.size() * sizeof(Mat3x::Scalar);
                return total;
            }
        };
        
        ProximityCollision(): 
            m_originalCTCollision( NULL ), mu(0.), force(Vec3x::Zero()), solved(false), valid(true)
        {
        }
        
        ContinuousTimeCollision* m_originalCTCollision; // If the proximity collision was actually built from a continuous time collision, this keeps a pointer to the original in case we need it for impulse resolution
        
		Scalar mu;		// need external initialization
		Scalar distance;	// need external initialization
        Vec3x normal;		// need external initialization
        Vec3x force;
		Vec3x offset;		// need external initialization
        Mat3x transformationMatrix ;	// compute by generateTransformationMatrix()

		bool solved;
		bool valid;
        
        std::pair<Object, Object> objects;
        
        void generateInverseInertia();
        void generateTransformationMatrix() ;
        void updateTransformationMatrix( const Mat3x& previous ) ;
        
        unsigned long computeSizeInBytes() const;
        
        void swapIfNecessary()
        {
            if (objects.first.globalIndex == -1 || (
                   objects.second.globalIndex != -1
                    && objects.second.globalIndex < objects.first.globalIndex ) )
            {
                normal = -normal;
                force = -force;
				offset = -offset;
                std::swap( objects.first, objects.second );
            }
        }
        
        void print( std::ostream& os ) const
        {
            os << "Collision: strand edge " << objects.first.globalIndex << ' ' << objects.first.vertex
            << " vs. strand edge " << objects.second.globalIndex << ' ' << objects.second.vertex
            << " N: " << normal << "\n";
        }
        
        bool operator<( const ProximityCollision &rhs ) const
        {
			if (objects.first.globalIndex == rhs.objects.first.globalIndex) {
				if (objects.second.globalIndex == rhs.objects.second.globalIndex) {
					if (objects.first.vertex == rhs.objects.first.vertex) {
						return objects.second.vertex < rhs.objects.second.vertex;
					}
					else
						return objects.first.vertex < rhs.objects.first.vertex;
				}
				else
					return objects.second.globalIndex < rhs.objects.second.globalIndex;
			}
			else
				return objects.first.globalIndex < rhs.objects.first.globalIndex;
        }

		bool operator== (const ProximityCollision & other) const
		{
			if (objects.first.globalIndex == other.objects.first.globalIndex
				&& objects.first.vertex == other.objects.first.vertex
				&& objects.second.globalIndex == other.objects.second.globalIndex
				&& objects.second.vertex == other.objects.second.vertex)
				return true;
			else
				return false;
		}

		struct Hasher
		{
			long operator()(const ProximityCollision& col) const
			{
				return std::hash<int>{}(col.objects.first.globalIndex) 
					^ std::hash<int>{}(col.objects.first.vertex) 
					^ std::hash<int>{}(col.objects.second.globalIndex) 
					^ std::hash<int>{}(col.objects.second.vertex);
			}
		};
    };
    
    class ProximityCollisionDatabase
    {
        struct Record
        {
            ProximityCollision collision ;
            unsigned short age ;
            bool firstTime ;
        } ;
        typedef std::pair< int, int > RecordKey ;
        typedef std::map< RecordKey, Record > Records;
        typedef std::unordered_map<int, Records> Table;
        typedef std::vector<Table> Base;
        typedef std::map< RecordKey, int > Counter;
        
        static const unsigned s_maxAge ;
        
    public:
        
        unsigned long computeSizeInBytes() const;
        
        explicit ProximityCollisionDatabase( const unsigned nObjs ) :
        m_nQueries( 0 ), m_nFound( 0 )
        {
            m_base.resize( nObjs );
        }
        
        const ProximityCollision* find( const ProximityCollision &needle ) const;
        
        void insert( const ProximityCollision& collision );
        
        void ageAll();
        
        void draw( const std::vector<ElasticStrand*>& strands ) const ;
        
        int numCollisions( int sIdx, int eIdx ) const;
        
        bool connectedLastStep(  int sP, int iP, int sQ, int iQ ) const;
        
    private:
        
        mutable unsigned m_nQueries;
        mutable unsigned m_nFound;
        
        Base m_base;
        Counter m_counter;
    };
    
    typedef std::vector<ProximityCollision> ProximityCollisions;
    
}

#endif // PROXIMITYCOLLISION_HH

