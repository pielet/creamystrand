/**
 * \copyright 2012 Jean-Marie Aubry, 2020 Shiyang Jia
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "EdgeEdgeCollision.hh"
#include "../Core/ElasticStrand.hh"
#include "../Utils/Distances.hh"
#include "CollisionUtils.hh"

#include "../Utils/TextLog.hh"

namespace strandsim
{
    static const double PARALLEL_TOL = 1e-12;

    void EdgeEdgeCollision::print(std::ostream& os) const
    {
        os << "EdgeEdgeCollision: strand " << m_firstStrand->getGlobalIndex() << " edge " << m_firstVertex
            << " s = " << m_firstAbscissa;
        os << " vs. strand " << m_secondStrand->getGlobalIndex() << " edge " << m_secondVertex
            << " t = " << m_secondAbscissa;
        os << " by " << -m_normalRelativeDisplacement << " at time = " << m_time << '\n';

        // more information
        os << "normal: " << m_normal << '\n';
        os << "strand edge1 move from " << m_firstStrand->getFutureState().getVertex(m_firstVertex)
            << " --- " << m_firstStrand->getFutureState().getVertex(m_firstVertex + 1)
            << " to " << m_firstStrand->getVertex(m_firstVertex) << " --- " << m_firstStrand->getVertex(m_firstVertex + 1) << '\n';
        os << "strand edge2 move from " << m_secondStrand->getFutureState().getVertex(m_secondVertex)
            << " --- " << m_secondStrand->getFutureState().getVertex(m_secondVertex + 1)
            << " to " << m_secondStrand->getVertex(m_secondVertex) << " --- " << m_secondStrand->getVertex(m_secondVertex + 1) << '\n';
    }

    bool EdgeEdgeCollision::analyse()
    {
        const Vec3x pp0 = m_firstStrand->getVertex(m_firstVertex);
        const Vec3x pp1 = m_firstStrand->getVertex(m_firstVertex + 1);
        const Vec3x pq0 = m_secondStrand->getVertex(m_secondVertex);
        const Vec3x pq1 = m_secondStrand->getVertex(m_secondVertex + 1);

        Vec3x dp0 = m_firstStrand->dynamics().getDisplacement(m_firstVertex);
        Vec3x dp1 = m_firstStrand->dynamics().getDisplacement(m_firstVertex + 1);
        Vec3x dq0 = m_secondStrand->dynamics().getDisplacement(m_secondVertex);
        Vec3x dq1 = m_secondStrand->dynamics().getDisplacement(m_secondVertex + 1);

        static double times[4];
        unsigned num_times;
        getCoplanarityTimes(pp0 - dp0, pp1 - dp1, pq0 - dq0, pq1 - dq1, pp0, pp1, pq0, pq1, times, NULL, num_times);

        for (unsigned i = 0; i < num_times; ++i) { 
            const Scalar dtime = times[i] - 1.0;

            const Vec3x p0col = pp0 + dtime * dp0;
            const Vec3x p1col = pp1 + dtime * dp1;
            const Vec3x q0col = pq0 + dtime * dq0;
            const Vec3x q1col = pq1 + dtime * dq1;

            Vec3x cp, cq;
            Scalar t;
            // if p//q, s = t = 0.0
            const double sqrt_dist = ClosestPtSegmentSegment(p0col, p1col, q0col, q1col, m_firstAbscissa, m_secondAbscissa, cp, cq);
            
            Scalar colRadius = m_firstStrand->collisionParameters().selfCollisionsRadius(m_firstVertex)
                                + m_secondStrand->collisionParameters().selfCollisionsRadius(m_secondVertex);
            
            if (sqrt_dist < colRadius * colRadius) {    // collision happens
                m_normal = ((pp1 - dp1) - (pp0 - dp0)).cross((pq1 - dq1) - (pq0 - -dq0));
                double nnorm = m_normal.norm();

                if (nnorm * nnorm <= PARALLEL_TOL) {    // parallel -> pass
                    continue;
                }

                m_normal /= nnorm;
                m_time = times[i];

                m_offset = (1. - m_firstAbscissa) * (pp0 - dp0) + m_firstAbscissa * (pp1 - dp1)   // p orig
                    - ((1. - m_secondAbscissa) * (pq0 - dq0) + m_secondAbscissa * (pq1 - dq1));   // q orig

                const Vec3x relativeDisplacement = m_time * 
                    (((1.0 - m_firstAbscissa) * dp0 + m_firstAbscissa * dp1)    // edge1 displacement
                    - ((1.0 - m_secondAbscissa) * dq0 + m_secondAbscissa * dq1));   // edge2 displacement
                postAnalyse(relativeDisplacement);  // adjust normal direction
                return true;
            }
        }

        return false;
    }

    bool compare(const EdgeEdgeCollision* ee1, const EdgeEdgeCollision* ee2)
    {
        if (ee1->m_firstStrand == ee2->m_firstStrand)
            if (ee1->m_firstVertex == ee2->m_firstVertex)
                if (ee1->m_secondStrand == ee2->m_secondStrand)
                    if (ee1->m_secondVertex == ee2->m_secondVertex)
                        return ee1->m_time < ee2->m_time;
                    else
                        return ee1->m_secondVertex < ee2->m_secondVertex;
                else
                    return ee1->m_secondStrand < ee2->m_secondStrand;
            else
                return ee1->m_firstVertex < ee2->m_firstVertex;
        else
            return ee1->m_firstStrand < ee2->m_firstStrand;
    }
}