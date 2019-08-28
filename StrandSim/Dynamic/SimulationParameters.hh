/**
 * \copyright 2012 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef SIMULATIONPARAMETERS_HH_
#define SIMULATIONPARAMETERS_HH_

#include "../Utils/MsgInfo.hh"
#include "../Core/Definitions.hh"

namespace strandsim
{

struct SimulationParameters
{
    MsgInfo::Severity m_logLevel;

    int m_numberOfThreads;
	int m_rodSubSteps;

    // Number of substeps
//    int m_solver_steps_per_frame;

    int m_statGathering;
    
    bool m_simulationManager_limitedMemory;

    unsigned m_maxNewtonIterations ;
    //    int m_maxLineSearchIterations;
    /**
     * Rod-rod collisions
     */
    bool m_skipRodMeshCollisions;
    bool m_skipRodRodCollisions; // whether we should use rod-rod proximity collisions (TODO: should rename; keeping name for now to maintain backwards comaptibility of older examples)
    bool m_skipFlowFrictions;
    bool m_useCTRodRodCollisions; // whether we should use rod-rod ctc collisions
    bool m_useApproxRodElasticFriction;
    Scalar m_percentCTRodRodCollisionsAccept;

    bool m_useImpulseMethod;

    double m_hairMeshFrictionCoefficient;
    double m_hairHairFrictionCoefficient;

    bool m_pruneExternalCollisions ;
    bool m_pruneSelfCollisions ;
    double m_stochasticPruning ;

    double m_airDrag ;
    /**
     * Inextensibility 
     */
    bool m_useLengthProjection;
    bool m_usePreFilterGeometry;

    bool m_useDeterministicSolver;

    unsigned m_gaussSeidelIterations;
    double m_gaussSeidelTolerance;
    bool m_useNonlinearContacts;
    bool m_useAdditionalExternalFailSafe;

    bool m_useSoftAttachConstraints;
    bool m_solveLiquids;
    
    unsigned m_subSteps;
};

}

#endif /* SIMULATIONPARAMETERS_HH_ */
