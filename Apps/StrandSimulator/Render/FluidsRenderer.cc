/**
 * \copyright 2010 Breannan Smith
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#define GL_SILENCE_DEPRECATION
#include "FluidsRenderer.hh"

#include "../../StrandSim/Render/Color.hh"
#include "../../StrandSim/Render/OpenGLDecl.hh"

#include "RenderingUtilities.hh"

namespace strandsim
{
    
    FluidsRenderer::FluidsRenderer( const std::shared_ptr<FluidScriptingController>& pcontrol, MutexType& particle_mutex, MutexType& grid_mutex, FluidsRenderer::DrawMode dm )
    : m_pcontroller(pcontrol), m_particle_mutex(particle_mutex), m_grid_mutex(grid_mutex), m_mode(dm)
    {
        const int num_components = m_pcontroller->getNumComponents();
        
        auto idx2rgb = [&] (int i) -> Vec3x {
            Scalar hh = (Scalar) i / (Scalar) max(1, num_components - 1) * 5.0;
            if(hh == 6.0)
                hh = 0.0;
            
            int fi = (int) floor(hh);
            Scalar ff = hh - floor(hh);
            
            switch(fi) {
                case 0:
                    return Vec3x(1.0, ff, 0.0);
                case 1:
                    return Vec3x(1.0 - ff, 1.0, 0.0);
                case 2:
                    return Vec3x(0.0, 1.0, ff);
                case 3:
                    return Vec3x(0.0, 1.0 - ff, 1.0);
                case 4:
                    return Vec3x(ff, 0.0, 1.0);
                case 5:
                default:
                    return Vec3x(1.0, 0.0, 1.0 - ff);
            }
        };

        m_component_colors.resize(num_components);
        for(int i = 0; i < num_components; ++i)
        {
            m_component_colors[i] = idx2rgb(i);
        }
    }
    
    void FluidsRenderer::render()
    {
        assert(m_pcontroller);
        
        int ni, nj, nk;
        m_pcontroller->getDimension(ni, nj, nk);
        Vec3x origin = m_pcontroller->getOrigin();
        Scalar dx = m_pcontroller->getDx();
        
        
        // Draw Grid
        glDisable(GL_DEPTH_TEST);
        glEnable (GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        
        // Draw Particle
        glPointSize(4.0);
        {
            LockGuard lock( m_particle_mutex );
            
            glBegin(GL_POINTS);
            int N = m_pcontroller->numParticles();
            for(int i = 0; i < N; ++i)
            {
                const Vec3x& pos = m_pcontroller->getParticles(i);
//                if(m_pcontroller->isWeakened(i))
//                    glColor4d(1.0, 0.0, 0.0, 0.1);
//                else
//                    glColor4d(0.0, 0.0, 1.0, 0.1);
                
                const VecXx& c = m_pcontroller->getComponents(i);
                Vec3x color = Vec3x::Zero();
                for(int j = 0; j < c.size(); ++j)
                {
                    color += m_component_colors[j] * c(j);
                }
                glColor4d(color(0), color(1), color(2), 0.1);
                glVertex3d(pos[0], pos[1], pos[2]);
            }
            glEnd();
        }
        
        // Draw Distance Grid
//        glPointSize(2.0);
//        {
//            LockGuard lock( m_grid_mutex );
//            const std::vector< VecXx >& distance = m_pcontroller->getNodeSignedDistance();
//            const Scalar r = m_pcontroller->getCellSize() * 0.5;
//            const Scalar dx = m_pcontroller->getCellSize();
//            const int num_nodes = m_pcontroller->getDefaultNumNodes() * 2;
//            
//            glBegin(GL_POINTS);
//            for(int bucket_idx = 0; bucket_idx < (int) distance.size(); ++bucket_idx)
//            {
//                const VecXx& vec = distance[ bucket_idx ];
//                if(vec.size() == 0) continue;
//                
//                Vec3x base_pos = m_pcontroller->getNodePos(bucket_idx, 0, 4);
//                for(int i = 0; i < num_nodes; ++i) for(int j = 0; j < num_nodes; ++j) for(int k = 0; k < num_nodes; ++k)
//                {
//                    Vec3x pos = base_pos + Vec3x(i, j, k) * r;
//                    int node_idx = k * num_nodes * num_nodes + j * num_nodes + i;
//                    const Scalar phi = vec(node_idx);
//                    if(phi < -0.5 * dx) continue;
//                    
////                    const Scalar c = clamp((phi - (-3.0 * dx)) / (6.0 * dx), 0.0, 1.0);
//                    
//                    const Vec3x c = renderingutils::interpolateColor(phi, -3.0 * dx, 3.0 * dx);
//                    
//                    glColor4d(c(0), c(1), c(2), 0.25);
//                    glVertex3d(pos(0), pos(1), pos(2));
//                }
//            }
//            glEnd();
//        }
        
        glDisable (GL_BLEND);
    }
    
    
    
}

