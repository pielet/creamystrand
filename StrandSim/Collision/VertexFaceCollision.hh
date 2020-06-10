/**
 * \copyright 2012 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef VERTEXFACECOLLISION_HH_
#define VERTEXFACECOLLISION_HH_

#include "ContinuousTimeCollision.hh"
#include "ElementProxy.hh"

namespace strandsim
{
    
    class VertexFaceCollision: public ContinuousTimeCollision, public FaceCollision
    {
    public:
        VertexFaceCollision( ElasticStrand* firstStrand, int firstVertex, const FaceProxy* faceProxy ) :
        ContinuousTimeCollision( firstStrand, firstVertex ), m_faceProxy( faceProxy )
        {
        }
        
        virtual ~VertexFaceCollision()
        {
        }
        
        virtual bool analyse();

		virtual Vec3x offset() const;
        
        friend bool compare( const VertexFaceCollision* vf1, const VertexFaceCollision* vf2 );
        
        Vec3x meshVelocity( Scalar dt ) const
        {
            return m_meshDisplacement / dt ;
        }
        
        const FaceProxy* face() const
        {
            return m_faceProxy ;
        }

    protected:
        void print( std::ostream& os ) const;
        
        const FaceProxy* const m_faceProxy;
        Vec3x m_meshDisplacement;
        Vec3x m_collisionOffset;
    };
    
}

#endif

