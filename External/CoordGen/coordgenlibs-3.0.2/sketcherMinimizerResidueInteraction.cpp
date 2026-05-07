/*
 *  sketcherMinimizerResidueInteraction.cpp
 *
 *  Created by Nicola Zonta on 13/05/2011.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#include "sketcherMinimizerResidueInteraction.h"

using namespace std;

sketcherMinimizerResidueInteraction::sketcherMinimizerResidueInteraction()
    : sketcherMinimizerBond()
{
}

sketcherMinimizerResidueInteraction::~sketcherMinimizerResidueInteraction() =
    default;

bool sketcherMinimizerResidueInteraction::isResidueInteraction()
{
    return true;
}

vector<sketcherMinimizerAtom*>
sketcherMinimizerResidueInteraction::getAllEndAtoms()
{
    vector<sketcherMinimizerAtom*> out = m_otherEndAtoms;
    out.insert(out.begin(), endAtom);
    return out;
}

vector<sketcherMinimizerAtom*>
sketcherMinimizerResidueInteraction::getAllStartAtoms()
{
    vector<sketcherMinimizerAtom*> out = m_otherStartAtoms;
    out.insert(out.begin(), startAtom);
    return out;
}