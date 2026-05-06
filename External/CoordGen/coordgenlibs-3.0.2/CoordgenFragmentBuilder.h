/*
 Contributors: Nicola Zonta
 Copyright Schrodinger, LLC. All rights reserved
 */

#ifndef COORDGEN_FRAGMENT_BUILDER_H
#define COORDGEN_FRAGMENT_BUILDER_H

#include <queue>
#include <set>
#include <stack>
#include <vector>

#include "CoordgenConfig.hpp"
#include "CoordgenMacrocycleBuilder.h"

class sketcherMinimizerAtom;
class sketcherMinimizerBond;
class sketcherMinimizerRing;
class sketcherMinimizerFragment;
class sketcherMinimizerPointF;

/*
 class that handles the creation of 2d coordinates for a molecular fragment
 */
class EXPORT_COORDGEN CoordgenFragmentBuilder
{
  public:
    /*
     create coordinates for a molecular fragment
     */
    void initializeCoordinates(sketcherMinimizerFragment* fragment) const;

    /*
     return a vector of ring atoms so that bound atoms are placed next to each
     other
     */
    static std::vector<sketcherMinimizerAtom*>
    orderRingAtoms(const sketcherMinimizerRing* r);

    /*
     return a vector of atoms so that bound atoms are placed next to each other
     */
    static std::vector<sketcherMinimizerAtom*>
    orderChainOfAtoms(const std::vector<sketcherMinimizerAtom*>& atoms,
                      sketcherMinimizerAtom* startAtom);

    /*
     return a list of coordinates representing a regular polygon for the given
     atoms in a ring
     */
    static std::vector<sketcherMinimizerPointF>
    listOfCoordinatesFromListofRingAtoms(
        const std::vector<sketcherMinimizerAtom*>& atoms);

    /*
     set a flag that forces the macrocycle builder to skip expensive polyomino
     matching routines and go straight to the breaking a bond approach
     */
    void setForceOpenMacrocycles(bool b)
    {
        m_macrocycleBuilder.m_forceOpenMacrocycles = b;
    }

    /* set precision of the calculations. Higher precisions settings result
     better
     quality but slower
     calculations */
    void setPrecision(float f) { m_macrocycleBuilder.setPrecision(f); }

    /*
     all bonds are placed at even intervals around the atom, as opposed for
     instance to the 90°-90°-120°-60° around tetracoordinated centers
     */
    bool m_evenAngles;

  private:
    /*
     find if the present ring is fused with another than has already gotten
     coordinates for
     */
    sketcherMinimizerRing* getSharedAtomsWithAlreadyDrawnRing(
        const sketcherMinimizerRing* ring,
        std::vector<sketcherMinimizerAtom*>& fusionAtoms,
        sketcherMinimizerBond*& fusionBond) const;

    /*
     assign coordinates to a ring
     */
    void buildRing(sketcherMinimizerRing* ring) const;

    /*
     generate coordinates for a group of fused rings that share more than two
     atoms with each other
     */
    void generateCoordinatesCentralRings(
        std::vector<sketcherMinimizerRing*> centralRings) const;
    sketcherMinimizerRing* findCentralRingOfSystem(
        const std::vector<sketcherMinimizerRing*>& rings) const;

    /*
     find a template to generate coordinates for a ring system
     */
    bool findTemplate(const std::vector<sketcherMinimizerRing*>& rings) const;

    /*
     generate coordinates for rings that have been stripped away from the core
     (see buildRings)
     */
    void generateCoordinatesSideRings(
        std::stack<sketcherMinimizerRing*> sideRings) const;

    /*
     after coordinates are generated, find an orientation for the main fragment
     */
    void rotateMainFragment(sketcherMinimizerFragment* fragment) const;

    /*
     assign coordinates to a fragment
     */
    void buildFragment(sketcherMinimizerFragment* fragment) const;

    /*
     assign coordinates to all ring atoms. Start by stripping out side rings
     that only share two atoms with other rings to find a core of central rings
     */
    void buildRings(sketcherMinimizerFragment* fragment) const;

    /*
     assign coordinates to atoms that are not in rings
     */
    void buildNonRingAtoms(sketcherMinimizerFragment* fragment) const;

    /*
     initialize information about connectivity of rings
     */
    void
    initializeFusedRingInformation(sketcherMinimizerFragment* fragment) const;

    /*
     split ring system into side rings and central rings by stripping away
     recursively rings that only
     share two atoms with other rings
     */
    void
    simplifyRingSystem(const std::vector<sketcherMinimizerRing*>& allRings,
                       std::stack<sketcherMinimizerRing*>& sideRings,
                       std::vector<sketcherMinimizerRing*>& centralRings) const;

    /* if the fragment contains any NaN coordinates and 3d coords are available,
     * use thouse instead */
    void fallbackIfNanCoordinates(sketcherMinimizerFragment* fragment) const;

    /* generate the coordinates of atoms bound to the first atom in the queue */
    void generateCoordinatesNeighborsOfFirstAtomInQueue(
        std::queue<sketcherMinimizerAtom*>& atomQueue,
        std::set<sketcherMinimizerAtom*>& isAtomVisited,
        const sketcherMinimizerFragment* fragment) const;

    /* return a list of angles that bonds from the given atoms should form */
    std::vector<float>
    neighborsAnglesAtCenter(const sketcherMinimizerAtom* atom) const;

    /* initialize data to generate coordinates of atoms bound to a non-ring atom
     */
    void initializeVariablesForNeighboursCoordinates(
        sketcherMinimizerAtom* atom,
        std::set<sketcherMinimizerAtom*>& isAtomVisited,
        sketcherMinimizerPointF& startCoordinates,
        std::vector<sketcherMinimizerAtom*>& orderedNeighbors,
        std::vector<float>& angles) const;

    /* initialize data to generate coordinates of atoms bound to a ring atom */
    void initializeVariablesForNeighboursCoordinatesRingAtom(
        const sketcherMinimizerAtom* atom,
        std::set<sketcherMinimizerAtom*>& isAtomVisited,
        sketcherMinimizerPointF& startCoordinates,
        std::vector<sketcherMinimizerAtom*>& orderedNeighbors,
        std::vector<float>& angles) const;

    /* check if the atom is part of a macrocycle and has some degrees of freedom
     that can be added to be used in the minimizer */
    void maybeAddMacrocycleDOF(sketcherMinimizerAtom* atom) const;

    /* make sure ZE chirality is maintained */
    void
    avoidZEInversions(const sketcherMinimizerAtom* at,
                      std::set<sketcherMinimizerAtom*>& isAtomVisited) const;

    /* assign a score to the possibility of rings to be drawn on a plane */
    float
    newScorePlanarity(const std::vector<sketcherMinimizerRing*>& rings) const;

    /* the macrocycle builder */
    CoordgenMacrocycleBuilder m_macrocycleBuilder;
};

#endif /* defined(COORDGEN_FRAGMENT_BUILDER_H) */
