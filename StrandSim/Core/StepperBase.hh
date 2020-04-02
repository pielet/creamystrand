/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef STEPPERBASE_HH_
#define STEPPERBASE_HH_

namespace strandsim
{

class StepperBase
{
public:
    StepperBase(ElasticStrand& strand): m_strand(strand) {}
    virtual ~StepperBase() {}

    VecXx& velocities() { return m_velocities; }
    const VecXx& velocities() const { return m_velocities; }

    ElasticStrand& getStrand() { return m_strand; }
    const ElasticStrand& getStrand() const { return m_strand; }

protected:
    VecXx m_velocities;
    ElasticStrand& m_strand;
};

}

#endif /* STEPPERBASE_HH_ */
