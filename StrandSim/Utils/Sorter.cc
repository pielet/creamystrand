/**
 * \copyright 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Sorter.hh"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

namespace strandsim {
Sorter::Sorter()
: ni(0), nj(0), nk(0)
{
	array_sup.resize(0);
}

Sorter::Sorter( int ni_, int nj_, int nk_ )
: ni(ni_), nj(nj_), nk(nk_)
{
	resize(ni, nj, nk);
}

Sorter::~Sorter() {
}

void Sorter::resize( int ni_, int nj_, int nk_ )
{
	array_sup.resize(ni_ * nj_ * nk_);
	ni = ni_; nj = nj_; nk = nk_;
}
};
