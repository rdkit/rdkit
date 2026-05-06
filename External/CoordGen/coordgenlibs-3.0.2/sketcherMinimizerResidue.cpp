/*
 *  sketcherMinimizerResidue.cpp
 *
 *  Created by Nicola Zonta on 13/05/2011.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#include "sketcherMinimizerResidue.h"

using namespace std;

sketcherMinimizerResidue::sketcherMinimizerResidue()
{
    m_closestLigandAtom = nullptr;
}

sketcherMinimizerResidue::~sketcherMinimizerResidue() = default;

bool sketcherMinimizerResidue::isResidue() const
{
    return true;
}
