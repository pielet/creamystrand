/**
 * \copyright 2012 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "VertexFaceCollision.hh"
#include "../Core/ElasticStrand.hh"
#include "ElementProxy.hh"
#include "CollisionUtils.hh"
#include "../Utils/Distances.hh"
#include "../Utils/TextLog.hh"
#include "../Dynamic/StrandDynamicTraits.hh"
#include "../Utils/MathUtilities.hh"

namespace strandsim
{

static const double SQ_TOLERANCE = 1e-12;
static const double EXTRA_RADIUS = 1e-3;

void VertexFaceCollision::print( std::ostream& os ) const
{
    os << "VertexFaceCollision: strand vertex " << m_firstStrand->getGlobalIndex() << ' '
            << m_firstVertex << " vs. " << m_faceProxy << " by " << -m_normalRelativeDisplacement
            << " at time = " << m_time;

    os << '\n';
    os << "Normal displacement = " << m_normal << '\n';
    os << "Vertex moved from: " << m_firstStrand->getFutureState().getVertex( m_firstVertex )
            << " to " << m_firstStrand->getVertex( m_firstVertex );
    os << '\n';
    os << "Face moved: " << *m_faceProxy << '\n';
    os << "Normal relative displacement: " << m_normalRelativeDisplacement;
    os << '\n';

}

bool VertexFaceCollision::analyse()
{
    static double times[4];

    const int num_verts = m_firstStrand->getNumVertices();
    const Vec3x p_off = m_firstStrand->getVertex( m_firstVertex );
    // NB after taking off the offset p = 0
    Vec3x pf0 = m_faceProxy->getVertex( 0 ) - p_off;
    Vec3x pf1 = m_faceProxy->getVertex( 1 ) - p_off;
    Vec3x pf2 = m_faceProxy->getVertex( 2 ) - p_off;

    const Vec3x dp = m_firstStrand->dynamics().getDisplacement( m_firstVertex );
    Vec3x df0 = m_faceProxy->getDisplacement( 0 );
    Vec3x df1 = m_faceProxy->getDisplacement( 1 );
    Vec3x df2 = m_faceProxy->getDisplacement( 2 );
	
	const Vec3x dfc = (df0 + df1 + df2) / 3.;

    int fnsign = m_faceProxy->knowsNormalSign( true, *m_firstStrand, m_firstVertex ) ;
    m_normal = m_faceProxy->getNormal() ;

    const Scalar colRadius = m_firstStrand->collisionParameters().externalCollisionsRadius( m_firstVertex, M_PI / 2 ) ;
    const Vec3x dpc = dp - dfc;
	const Scalar disp_rel = dpc.dot(m_normal);
    const Scalar tang_rel = sqrt(std::max(0., dpc.squaredNorm() - disp_rel * disp_rel));

    unsigned num_times ;
    getCoplanarityTimes( -dp, pf0 - df0, pf1 - df1, pf2 - df2, Vec3x(), pf0, pf1, pf2, times,
            NULL, num_times );

    for ( size_t j = 0; j < num_times; ++j )
    {
        // TODO: Use barycentric coordinates or point-triangle closest point < epsilon here? closest point < epsilon really just extends the triangle a little bit.
        // Determine if the collision actually happens
        const Scalar dtime = times[j] - 1.0;
        const Vec3x pcol = dtime * ( dp );
        const Vec3x f0col = pf0 + dtime * df0;
        const Vec3x f1col = pf1 + dtime * df1;
        const Vec3x f2col = pf2 + dtime * df2;

        const Vec3x cp = ClosestPtPointTriangle( pcol, f0col, f1col, f2col );

        const Scalar pqdist2 = ( pcol - cp ).squaredNorm();
        // If, when they are coplanar, the objects are sufficiently close, register a collision
        if ( pqdist2 < colRadius * colRadius )
        {
            m_time = times[j];
        
            computeBarycentricCoordinates( f0col, f1col, f2col, pcol, m_u, m_v, m_w );
            // computeBarycentricCoordinates coords could be outside of [0,1] right now because we've extended the triangles a little bit
            assert( isSmall(m_u + m_v + m_w - 1.0) );

            m_meshDisplacement = m_u * df0 + m_v * df1 + m_w * df2;
            const Vec3x relativeDisplacement = ( 1 - m_time ) * ( dp - m_meshDisplacement );

			m_offset = -dp - ((m_u * pf0 + m_v * pf1 + m_w * pf2) - m_meshDisplacement); // orig point on mesh

            const Scalar nDepl = relativeDisplacement.dot( m_normal ) ;
            if( !fnsign )
            {
                // Normal sign was unknown, now we know that it should be opposed to relativeDisplacement
                fnsign = ( nDepl > 0. ? -1 : 1 ) ;
                m_faceProxy->setNormalSign( fnsign, m_time, *m_firstStrand, m_firstVertex );
            }
            else {
                if ( fnsign * nDepl > 0. )
                {
                    return false;
                }
            }
            m_normal = fnsign * m_normal;
            
            //postAnalyse( relativeDisplacement );

            return true;
        }
    }
    return false;

}

Vec3x VertexFaceCollision::offset() const
{
	const Scalar colRadius = m_firstStrand->collisionParameters().externalResponseRadius();
	return (m_offset.dot(m_normal) - colRadius) * m_normal;
}

bool compare( const VertexFaceCollision* vf1, const VertexFaceCollision* vf2 )
{
    if ( vf1->m_firstStrand == vf2->m_firstStrand )
        if ( vf1->m_firstVertex == vf2->m_firstVertex )
            if ( vf1->m_faceProxy == vf2->m_faceProxy )
                return vf1->m_time < vf2->m_time;
            else
                return vf1->m_faceProxy < vf2->m_faceProxy;
        else
            return vf1->m_firstVertex < vf2->m_firstVertex;
    else
        return vf1->m_firstStrand < vf2->m_firstStrand;
}

}
