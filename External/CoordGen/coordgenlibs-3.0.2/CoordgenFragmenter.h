/*
   Contributors: Nicola Zonta
   Copyright Schrodinger, LLC. All rights reserved
 */
#pragma once
#include "CoordgenConfig.hpp"

#include <cstddef>
#include <vector>

class sketcherMinimizerFragment;
class sketcherMinimizerMolecule;
class sketcherMinimizerBond;
class sketcherMinimizerRing;
class sketcherMinimizerAtom;

/*
 class to divide a molecule into rigid fragments
 */
class EXPORT_COORDGEN CoordgenFragmenter
{
  public:
    /*
     divide the molecule into rigid fragments
     */
    static void splitIntoFragments(sketcherMinimizerMolecule* molecule);

  private:
    /*
     put all atoms and bonds of fragment2 into fragment1 and delete it
     */
    static void
    joinFragments(sketcherMinimizerFragment* fragment1,
                  sketcherMinimizerFragment* fragment2,
                  std::vector<sketcherMinimizerFragment*>& fragments);

    /*
     process bond between two fragments
     */
    static void processInterFragmentBond(
        sketcherMinimizerBond* bond,
        std::vector<sketcherMinimizerFragment*>& fragments);

    /*
     process bond internal a fragment
     */
    static void processBondInternalToFragment(
        sketcherMinimizerBond* bond,
        std::vector<sketcherMinimizerFragment*>& fragments);

    /*
     initialize info from bond
     */
    static void addBondInformation(sketcherMinimizerBond* bond);

    /*
     initialize info from ring
     */
    static void addRingInformation(sketcherMinimizerRing* ring);

    /*
     initialize info of molecule
     */
    static void
    initializeInformation(std::vector<sketcherMinimizerFragment*> fragments,
                          sketcherMinimizerMolecule* molecule);

    /*
      compare two fragments for priority (e.g. to find the main fragment)
     */
    static bool hasPriority(const sketcherMinimizerFragment* fragment1,
                            const sketcherMinimizerFragment* fragment2);

    /*
     get the score of a particular descriptor for the given fragment (used to
     assign priorities)
     */
    static size_t getValueOfCheck(const sketcherMinimizerFragment* fragment,
                                  int checkN, bool& checkNoMore);

    /*
     find the main fragment
     */
    static sketcherMinimizerFragment*
    findMainFragment(const std::vector<sketcherMinimizerFragment*>& fragments);

    /*
     if molecule has a long enough chain return the start of the chain as main
     fragment instead
     */
    static sketcherMinimizerFragment*
    considerChains(const std::vector<sketcherMinimizerFragment*>& fragments,
                   sketcherMinimizerFragment* mainFragment);

    /* empirical minimum length of zigzag chain of fragments
     that makes the chain get priority over the main fragment
     in determining the molecule's layout
     */
    static unsigned int
    acceptableChainLength(sketcherMinimizerFragment* mainFragment);

    /*
     find the longest chain of connected fragments
     */
    static std::vector<sketcherMinimizerFragment*>
    findLongestChain(const std::vector<sketcherMinimizerFragment*>& fragments);

    /*
     initialize parent-child info of bound fragments
     */
    static void addParentRelationsToFragments(
        sketcherMinimizerFragment* mainFragment,
        const std::vector<sketcherMinimizerFragment*>& fragments);

    /*
     order the vector of fragments so that bound fragments are consecutive,
     starting from the main fragment
     */

    static void
    orderFragments(std::vector<sketcherMinimizerFragment*>& fragments,
                   sketcherMinimizerFragment* mainFragment);
};
