/**
 * \copyright 2012 Jean-Marie Aubry, 2020 Shiyang Jia
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EDGEEDGECOLLISION_HH_
#define EDGEEDGECOLLISION_HH_

#include "ContinuousTimeCollision.hh"
#include "ElementProxy.hh"

namespace strandsim
{
    class EdgeEdgeCollision : public ContinuousTimeCollision
    {
    public:
        EdgeEdgeCollision(ElasticStrand* firstStrand, int firstVertex, ElasticStrand* secondStrand, int secondVertex) :
            ContinuousTimeCollision(firstStrand, firstVertex), m_secondStrand(secondStrand), m_secondVertex(secondVertex),
            m_firstAbscissa(0.), m_secondAbscissa(0.), m_doSOCSolve(false)
        {
        }

        ElasticStrand* getSecondStrand() const { return m_secondStrand; }

        int getSecondVertex() const { return m_secondVertex; }

        Scalar getFirstAbscissa() const {  return m_firstAbscissa; }

        Scalar getSecondAbscissa() const { return m_secondAbscissa; }

        bool doSOCSolve() const { return m_doSOCSolve; }

        virtual bool analyse();

        friend bool compare(const EdgeEdgeCollision* ee1, const EdgeEdgeCollision* ee2);

    private:
        void print(std::ostream& os) const;

        ElasticStrand* m_secondStrand;
        int m_secondVertex;

        Scalar m_firstAbscissa;
        Scalar m_secondAbscissa;

        bool m_doSOCSolve;
    };
}

#endif // !EDGEEDGECOLLISION_HH_
