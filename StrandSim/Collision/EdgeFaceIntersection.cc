/**
 * \copyright 2012 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "../Core/ElasticStrand.hh"
#include "EdgeFaceIntersection.hh"
#include "ElementProxy.hh"
#include "CollisionUtils.hh"
#include "../Utils/Distances.hh"
#include "../Dynamic/StrandDynamicTraits.hh"
#include "../Utils/MathUtilities.hh"

namespace strandsim
{
    
    static const double SQ_TOLERANCE = 1e-12;
    double EdgeFaceIntersection::s_doProximityDetection = false;
    
    EdgeFaceIntersection::EdgeFaceIntersection( const ElasticStrand* strand, int edge, EdgeProxy* edgeP,
                                               const strandsim::FaceProxy* triangle ) :
    m_strand( strand ), m_edge( edge ), m_edgeP( edgeP ), m_triangle( triangle )
    {
    }
    
    EdgeFaceIntersection::~EdgeFaceIntersection()
    {
    }
    
    void EdgeFaceIntersection::print( std::ostream& os ) const
    {
        os << "EdgeFaceIntersection: strand vertex " << m_strand << ' ' << m_edge << " vs. "
        << m_triangle;
    }
    
    bool EdgeFaceIntersection::barycentricCoordinatesIfCloser( const Scalar radius, Scalar s0,
                                                              const Vec3x& p0, const Vec3x &p1, const Vec3x &q0, const Vec3x &q1, const Vec3x &q2,
                                                              Scalar &dist, Scalar &s, Scalar &u, Scalar &v, Scalar &w )
    {
        const Vec3x pcol = p0 + s0 * ( p1 - p0 );
        const Vec3x cp = ClosestPtPointTriangle( pcol, q0, q1, q2 );
        
        const Scalar d = ( pcol - cp ).dot( m_normal );

        
        if ( d > -radius && d < dist )
        {
            Scalar lu, lv, lw;
            const Vec3x projOnTri = pcol - d * m_normal;
            computeBarycentricCoordinates( q0, q1, q2, projOnTri, lu, lv, lw );
            
            if ( lu > 0 && lv > 0 && lw > 0 )
            {
                s = s0;
                u = lu;
                v = lv;
                w = lw;
                dist = d;
                
                // calculate adhesive length
                const Vec3x panother = p0 + (1. - s0) * (p1 - p0);
                const Vec3x cpanother = ClosestPtPointTriangle(panother, q0, q1, q2);
                const Scalar dcol = std::max(0., d);
                const Scalar danother = std::max(0., ( panother - cpanother ).dot( m_normal ));
                
                return true;
            }
            
        }
        
        return false;
    }
    
    bool EdgeFaceIntersection::analyseProximity()
    {
        static double times[4];
        
        // Tolerance for the interior edge/edge collisions that point in a direction perpendicular to the face normal
        const Scalar perpEdgeTol = .1 ;
        
        const Vec3x offset = m_strand->getVertex( m_edge );
        
        Vec3x p0( Vec3x::Zero() );
        Vec3x p1 = m_strand->getVertex( m_edge + 1 ) - offset;

		const Vec3x dp0 = m_strand->dynamics().getDisplacement( m_edge );
		const Vec3x dp1 = m_strand->dynamics().getDisplacement( m_edge + 1 );
		const Vec3x dpc = (dp0 + dp1) * 0.5;
		
        p0 -= dp0;
        p1 -= dp1;
        
        Vec3x q0 = m_triangle->getVertex( 0 ) - offset;
        Vec3x q1 = m_triangle->getVertex( 1 ) - offset;
        Vec3x q2 = m_triangle->getVertex( 2 ) - offset;
        
        m_normal = m_triangle->getNormal();
		
		const Vec3x dq0 = m_triangle->getDisplacement( 0 );
		const Vec3x dq1 = m_triangle->getDisplacement( 1 );
		const Vec3x dq2 = m_triangle->getDisplacement( 2 );
		
		q0 -= dq0;
		q1 -= dq1;
		q2 -= dq2;

		const Vec3x dqc = (dq0 + dq1 + dq2) / 3.;
		
        int nsign = m_triangle->knowsNormalSign( false, *m_strand, m_edge ) ;
        if( nsign )
        {
            m_normal *= nsign ;
        } else if( !m_triangle->getFace().hasBoundaryEdges() ) {
            return false ;  // Do not try to make educated guesses, trust the CT collisions
        }
		
        Vec3x dpq = dpc - dqc;
		const Scalar disp_rel = dpq.dot(m_normal);
        const Scalar tang_rel = sqrt(std::max(dpq.squaredNorm() - disp_rel * disp_rel, 0.));
		
		Scalar radius = getEllipticExternalCollisionOffset( *m_strand, m_edge, m_normal );
		
        bool found = false;
        
        m_distance = radius;
        m_isBoundaryedgeCollision = false ;
        bool intersecting = false ;
        
        // Perform the face proximity tests only if we have a valid sign
        
        //   // I - Vertex face proximity ( each edge checks only its second vertex )
        if( nsign )
        {
            found = barycentricCoordinatesIfCloser( radius, 1, p0, p1, q0, q1, q2, m_distance, m_s, m_u, m_v, m_w ) || found;
        }
        
        if ( !(m_strand->isVertexFreezed(m_edge) && m_strand->isVertexFreezed(m_edge + 1)) )
        {
            
            // II - Edge Triangle intersection
            
            unsigned num_times ;
            getIntersectionPoint( p0, p1, q0, q1, q2, times, NULL, num_times );
            
            
            // If nsign is not set, we cannot accept the collision, but we still want to know if there is an intersection
            for ( unsigned j = 0; j < num_times ; ++j )
            {
                const Scalar s0 = times[j] ;
                const Vec3x pcol = p0 + s0 * ( p1 - p0 );
                
                Scalar lu, lv, lw;
                computeBarycentricCoordinates( q0, q1, q2, pcol, lu, lv, lw );
                
                if ( lu > 0 && lv > 0 && lw > 0 )
                {
                    intersecting = true ;
                    
                    if( nsign && m_distance > 0 )
                    {
                        m_s = s0;
                        m_u = lu;
                        m_v = lv;
                        m_w = lw;
                        
                        found = true ;
                        m_distance = 0.;
                        
                        // find the adhesive length
                        const Vec3x cp0 = ClosestPtPointTriangle(p0, q0, q1, q2);
                        const Vec3x cp1 = ClosestPtPointTriangle(p1, q0, q1, q2);
                        
                        // there must be one of them equals zero, but doesn't matter since wedge_wet_abscissa would naturally handle this case
                        const Scalar d0 = std::max(0., ( p0 - cp0 ).dot( m_normal ));
                        const Scalar d1 = std::max(0., ( p1 - cp1 ).dot( m_normal ));
                    }
                    
                    break;
                }
            }
            
        }
        
        // III - Edge Edge proximity
        
        Vec3x* tri[3] = { &q0, &q1, &q2 };
        
        const Vec3x edge = ( p1 - p0 ).normalized();
        const Vec3x faceNormal = m_normal;
        for ( unsigned k = 0; k < 3; ++k )
        {
            bool isBoundaryEdge = m_triangle->getFace().edgeOnBoundary( k ) ;
            if( !isBoundaryEdge && !nsign ) continue ;
            
            const Vec3x& e0 = *tri[k];
            const Vec3x& e1 = *tri[( k + 1 ) % 3];
            const Vec3x meshEdge = ( e1 - e0 ).normalized() ;
            
            Vec3x edgeNormal = meshEdge.cross( edge );
            const Scalar enn = edgeNormal.norm();
            
            if ( isSmall( enn ) )
                continue;
            edgeNormal /= enn;
            
            // We need a reference normal to deduce to sign of the edge/edge collision normal
            Vec3x origNormal ;
            if( isBoundaryEdge )
            {
                const Vec3x& e2 = *tri[( k + 2 ) % 3];
                
                origNormal = ( e0 - e2 ) - ( e0 - e2 ).dot( meshEdge ) * meshEdge ;
                assert( !isSmall( origNormal.norm() ) ) ;
                origNormal = origNormal.normalized() ;
            } else {
                origNormal = faceNormal ;
            }
            
            const Scalar proj = edgeNormal.dot( origNormal );
            
            // Ignore the collision perpendicular to a face normal on interior edges
            if ( !isBoundaryEdge && std::fabs( proj ) < perpEdgeTol )
                continue;
            
            if ( proj < 0. )
                edgeNormal = -edgeNormal;
            
            
            Scalar es, et;
            Vec3x pcols ;
            Vec3x pcolt ;
            ClosestPtSegmentSegment( p0, p1, e0, e1, es, et, pcols, pcolt );
            
            // Check if both closest points are inside the segments
            if ( es <= 0 || es == 1 || et <= 0 || et == 1 )
            {
                // If not, but we know that the segemtn is intersecting the triangle, we can use
                // the projection of the intersection point instead
                if( isBoundaryEdge && intersecting )
                {
                    es = times[0] ;
                    pcols = p0 + es * ( p1 - p0 ) ;
                    pcolt = ClosestPtPointSegment( pcols, e0, e1 ) ;
                } else
                    continue;
            }
            
            const Vec3x depl = ( pcols - pcolt ) ;
            const Scalar dist = depl.dot( edgeNormal ) ;
            
            // For boundary edges, we accept the collision that has the least penetration distance ;
            // for inner edges, we accept the maximum one
            bool accept = false ;
            
            if( isBoundaryEdge ) {
                
                accept =  dist < radius && ( !m_isBoundaryedgeCollision || dist > m_distance )
                && ( std::fabs( proj ) > perpEdgeTol || dist > -radius )
                && ( dist > 0 || intersecting ) ;
                
                // Discard edge collisions opposed to face normal when the strand is far from the boundary edge
                if( nsign && accept )
                {
                    accept = dist > -radius || edgeNormal.dot( faceNormal ) > perpEdgeTol ;
                }
                
                m_isBoundaryedgeCollision = m_isBoundaryedgeCollision || accept ;
                
            } else {
                
                accept = dist > -radius && !m_isBoundaryedgeCollision && dist < m_distance ;

            }
            
            if ( accept )
            {
                
                computeBarycentricCoordinates( q0, q1, q2, pcolt, m_u, m_v, m_w );
                
                m_distance = dist;
                found = true;
                m_s = (float) es;
                m_normal = edgeNormal;
               
            }
            
        }
        
        return found ;
    }
    
    bool EdgeFaceIntersection::analyse()
    {
        
        if ( s_doProximityDetection )
        {
            return analyseProximity();
        }
        m_distance = 0.;
        
        const Vec3x offset = m_strand->getVertex( m_edge );
        
        const Vec3x p0( Vec3x::Zero() );
        const Vec3x p1 = m_strand->getVertex( m_edge + 1 ) - offset;
        
        const Vec3x q0 = m_triangle->getVertex( 0 ) - offset;
        const Vec3x q1 = m_triangle->getVertex( 1 ) - offset;
        const Vec3x q2 = m_triangle->getVertex( 2 ) - offset;
        
        const Vec3x dp0 = m_strand->dynamics().getDisplacement( m_edge );
        const Vec3x dp1 = m_strand->dynamics().getDisplacement( m_edge + 1 );
        
        const Vec3x dq0 = m_triangle->getDisplacement( 0 );
        const Vec3x dq1 = m_triangle->getDisplacement( 1 );
        const Vec3x dq2 = m_triangle->getDisplacement( 2 );
        
        m_normal = getFaceNormal();
        
        double times[4];
        double errors[4];
        unsigned num_times ;
        
        getIntersectionPoint( p0, p1, q0, q1, q2, times, errors, num_times );
        
        for ( size_t j = 0; j < num_times ; ++j )
        {
            double dtime = times[j] * 1.0;
            
            // Determine if the collision actually happens
            const Vec3x pcol = p0 + dtime * ( p1 - p0 );
            const Vec3x f0col = q0;
            const Vec3x f1col = q1;
            const Vec3x f2col = q2;
            
            Vec3x cp = ClosestPtPointTriangle( pcol, f0col, f1col, f2col );
            
            // If the intersection point and the close point on the triangle are close and that point is interior, register a collision
            if ( ( pcol - cp ).squaredNorm() <= SQ_TOLERANCE )
            {
                m_s = dtime;
                
                computeBarycentricCoordinates( f0col, f1col, f2col, pcol, m_u, m_v, m_w );
                // computeBarycentricCoordinates coords could be outside of [0,1] right now because we've extended the triangles a little bit
                assert( isSmall(m_u + m_v + m_w- 1.0) );
                
                if ( m_u > 0 && m_v > 0 && m_w > 0 ) {
                    return true;
                }
            }
        }
        
        return false;
        
    }
    
    Vec3x EdgeFaceIntersection::getFaceNormal() const
    {
        return m_normal;
    }
    
    Vec3x EdgeFaceIntersection::getFaceVelocity( const Scalar dt ) const
    {
        const Vec3x d0 = m_triangle->getDisplacement( 0 );
        const Vec3x d1 = m_triangle->getDisplacement( 1 );
        const Vec3x d2 = m_triangle->getDisplacement( 2 );
        
        Vec3x depl = ( m_u * d0 + m_v * d1 + m_w * d2 );
        
        ////        std::cout << m_distance << std::endl ;
        
        if( m_isBoundaryedgeCollision )
        {
            //        const Scalar radius = 1.e-3 ;
            const Scalar radius = m_strand->collisionParameters().externalCollisionsRadius( m_edge ) ;
            if( radius > m_distance )
            {
                depl += ( 1.e-3 - std::max( m_distance, -radius ) ) * m_normal  ;
            }
        }
        
        return depl / dt;
    }
    
    uint32_t EdgeFaceIntersection::getFaceId() const
    {
        return m_triangle->uniqueId() + 0x80000000u;
    }
    
} /* namespace strandsim */

