/*
 Contributors: Nicola Zonta
 Copyright Schrodinger, LLC. All rights reserved
 */

#include "CoordgenFragmenter.h"
#include "sketcherMinimizerAtom.h"
#include "sketcherMinimizerBond.h"
#include "sketcherMinimizerFragment.h"
#include "sketcherMinimizerMolecule.h"
#include "sketcherMinimizerRing.h"

#include "sketcherMinimizerMaths.h"
#include <algorithm>
#include <queue>
using namespace std;

void CoordgenFragmenter::splitIntoFragments(sketcherMinimizerMolecule* molecule)
{
    /*
     split input molecule into rigid fragments
     */
    vector<sketcherMinimizerFragment*> fragments;
    for (sketcherMinimizerAtom* atom : molecule->getAtoms()) {
        atom->setFragment(nullptr);
    }

    /* molecule with a single atom and no bonds */
    if (molecule->getAtoms().size() == 1) {
        auto* fragment = new sketcherMinimizerFragment();
        sketcherMinimizerAtom* onlyAtom = molecule->getAtoms().at(0);
        fragment->addAtom(onlyAtom);
        fragments.push_back(fragment);
    }
    for (sketcherMinimizerBond* bond : molecule->_bonds) {
        if (bond->isResidueInteraction()) {
            continue;
        }
        if (bond->isInterFragment()) {
            processInterFragmentBond(bond, fragments);
        } else {
            processBondInternalToFragment(bond, fragments);
        }
    }
#ifndef NDEBUG
    for (auto atom : molecule->getAtoms()) {
        assert(atom->getFragment() != nullptr);
    }
#endif
    if (!fragments.empty()) {
        initializeInformation(fragments, molecule);
    }
}

void CoordgenFragmenter::addBondInformation(sketcherMinimizerBond* bond)
{
    if (bond->isResidueInteraction()) {
        return;
    }
    if (bond->getStartAtom()->getFragment() ==
        bond->getEndAtom()->getFragment()) {
        bond->getStartAtom()->getFragment()->addBond(bond);

    } else {
        bond->getStartAtom()->getFragment()->addInterFragmentBond(bond);
        bond->getEndAtom()->getFragment()->addInterFragmentBond(bond);
    }
}

void CoordgenFragmenter::addRingInformation(sketcherMinimizerRing* ring)
{
    ring->_atoms.at(0)->getFragment()->addRing(ring);
}

void CoordgenFragmenter::processInterFragmentBond(
    sketcherMinimizerBond* bond, vector<sketcherMinimizerFragment*>& fragments)
{
    if (bond->getStartAtom()->getFragment() == nullptr) {
        auto* fragment = new sketcherMinimizerFragment();
        fragment->addAtom(bond->getStartAtom());
        fragments.push_back(fragment);
    }
    if (bond->getEndAtom()->getFragment() == nullptr) {
        auto* fragment = new sketcherMinimizerFragment();
        fragment->addAtom(bond->getEndAtom());
        fragments.push_back(fragment);
    }
}

void CoordgenFragmenter::processBondInternalToFragment(
    sketcherMinimizerBond* bond, vector<sketcherMinimizerFragment*>& fragments)
{
    if (bond->getStartAtom()->getFragment() == nullptr &&
        bond->getEndAtom()->getFragment() == nullptr) {
        /* add the two atoms to a new fragment */
        auto* fragment = new sketcherMinimizerFragment();
        fragment->addAtom(bond->getStartAtom());
        fragment->addAtom(bond->getEndAtom());
        fragments.push_back(fragment);
    } else if (bond->getEndAtom()->getFragment() == nullptr) {
        /* extend fragment of start atom */
        bond->getStartAtom()->getFragment()->addAtom(bond->getEndAtom());
    } else if (bond->getStartAtom()->getFragment() == nullptr) {
        /* extend fragment of end atom */
        bond->getEndAtom()->getFragment()->addAtom(bond->getStartAtom());
    } else if (bond->getStartAtom()->getFragment() !=
               bond->getEndAtom()->getFragment()) {
        /* join the two fragments */
        joinFragments(bond->getStartAtom()->getFragment(),
                      bond->getEndAtom()->getFragment(), fragments);
    }
}

void CoordgenFragmenter::joinFragments(
    sketcherMinimizerFragment* fragment1, sketcherMinimizerFragment* fragment2,
    vector<sketcherMinimizerFragment*>& fragments)
{
    for (sketcherMinimizerAtom* atom : fragment2->atoms()) {
        fragment1->addAtom(atom);
    }

    // remove fragment2 from fragments and free memory
    auto positionOfFragment2 =
        remove(fragments.begin(), fragments.end(), fragment2);
    fragments.erase(positionOfFragment2, fragments.end());
    delete fragment2;
}

// Check if atom has fixed cooordinates
static bool isAtomFixed(const sketcherMinimizerAtom* atom)
{
    // needed for the find_if algorithm
    return atom->fixed;
}

// Check if atom has constrained cooordinates
static bool isAtomConstrained(const sketcherMinimizerAtom* atom)
{
    // needed for the find_if algorithm
    return atom->constrained;
}

// Check if fragment is part of an aliphatic chain
static bool isChain(const sketcherMinimizerFragment* fragment)
{
    /* can this fragment be part of a zig-zag chain, e.g. exane? */
    vector<sketcherMinimizerAtom*> fragmentAtoms = fragment->getAtoms();
    if (fragmentAtoms.size() > 3) {
        return false;
    }
    for (sketcherMinimizerAtom* atom : fragmentAtoms) {
        if (atom->getBonds().size() > 3) {
            return false;
        }
        if (!atom->getRings().empty()) {
            return false;
        }
    }
    vector<sketcherMinimizerBond*> fragmentBonds = fragment->getBonds();
    for (sketcherMinimizerBond* bond : fragmentBonds) {
        if (bond->getBondOrder() > 2) {
            return false;
        }
    }
    return true;
}

// Check if fragment can be flipped and sets the appropriate flag
static void setFlipInfo(sketcherMinimizerFragment* fragment)
{
    fragment->constrainedFlip = false;
    if (fragment->getAtoms().size() == 1 && fragment->countConstrainedAtoms() == 1) {
        // fragment has exactly 1 atom. should only be constrained if a child is constrained
        for (auto child : fragment->_children) {
            if(child->constrained) {
                fragment->constrainedFlip = true;
            }
        }
    } else {
        fragment->constrainedFlip  = (fragment->countConstrainedAtoms() > 1);
    }
}

// Check if fragment has fixed coordinates and sets the appropriate flag
static bool setFixedInfo(sketcherMinimizerFragment* fragment)
{
    fragment->fixed =
        (find_if(fragment->atoms().begin(), fragment->atoms().end(),
                 isAtomFixed) != fragment->atoms().end());
    return fragment->fixed;
}

// Check if fragment has constrained coordinates and set the appropriate flag
static bool setConstrainedInfo(sketcherMinimizerFragment* fragment)
{
    fragment->constrained =
        (find_if(fragment->atoms().begin(), fragment->atoms().end(),
                 isAtomConstrained) != fragment->atoms().end());
    return fragment->constrained;
}

// Check if fragment is part of an aliphatic chain and sets the appropriate
// flag
static void setChainInfo(sketcherMinimizerFragment* fragment)
{
    fragment->isChain = isChain(fragment);
}

void CoordgenFragmenter::initializeInformation(
    vector<sketcherMinimizerFragment*> fragments,
    sketcherMinimizerMolecule* molecule)
{
    for_each(molecule->getBonds().begin(), molecule->getBonds().end(),
             addBondInformation);
    for_each(molecule->getRings().begin(), molecule->getRings().end(),
             addRingInformation);
    for_each(fragments.begin(), fragments.end(), setChainInfo);
    molecule->hasConstrainedFragments =
        count_if(fragments.begin(), fragments.end(), setConstrainedInfo) > 0;
    molecule->hasFixedFragments =
        count_if(fragments.begin(), fragments.end(), setFixedInfo) > 0;
    sketcherMinimizerFragment* mainFragment = findMainFragment(fragments);
    addParentRelationsToFragments(mainFragment, fragments);
    orderFragments(fragments, mainFragment);
    molecule->setMainFragment(mainFragment);
    molecule->setFragments(fragments);
    for_each(fragments.begin(), fragments.end(), setFlipInfo);
}

bool CoordgenFragmenter::hasPriority(const sketcherMinimizerFragment* fragment1,
                                     const sketcherMinimizerFragment* fragment2)
{
    bool checkNoMore = false;
    int checkN = 0;
    while (!checkNoMore) {
        size_t leftValue = getValueOfCheck(fragment1, checkN, checkNoMore);
        size_t rightValue = getValueOfCheck(fragment2, checkN, checkNoMore);
        if (leftValue > rightValue) {
            return true;
        }
        if (leftValue < rightValue) {
            return false;
        }
        ++checkN;
    }
    return false;
}

size_t
CoordgenFragmenter::getValueOfCheck(const sketcherMinimizerFragment* fragment,
                                    int checkN, bool& checkNoMore)
{
    switch (checkN) {
    case 0:
        return fragment->countFixedAtoms();
    case 1:
        return fragment->countConstrainedAtoms();
    case 2:
        return fragment->getRings().size();
    case 3:
        return fragment->getAtoms().size();
    case 4:
        return fragment->_interFragmentBonds.size();
    case 5:
        return fragment->countHeavyAtoms();
    case 6:
        return fragment->totalWeight();
    case 7:
        return fragment->countDoubleBonds();
    default:
        checkNoMore = true;
        return 0;
    }
}

sketcherMinimizerFragment* CoordgenFragmenter::findMainFragment(
    const vector<sketcherMinimizerFragment*>& fragments)
{
    sketcherMinimizerFragment* mainFragment =
        *(min_element(fragments.begin(), fragments.end(), hasPriority));
    mainFragment = considerChains(fragments, mainFragment);
    return mainFragment;
}

sketcherMinimizerFragment* CoordgenFragmenter::considerChains(
    const vector<sketcherMinimizerFragment*>& fragments,
    sketcherMinimizerFragment* mainFragment)
{

    for (sketcherMinimizerFragment* fragment : fragments) {
        if (fragment->fixed || fragment->constrained) {
            return mainFragment;
        }
    }
    vector<sketcherMinimizerFragment*> longestChain =
        findLongestChain(fragments);
    if (longestChain.size() >= acceptableChainLength(mainFragment)) {
        mainFragment = longestChain.at(0);
    }
    return mainFragment;
}

unsigned int CoordgenFragmenter::acceptableChainLength(
    sketcherMinimizerFragment* mainFragment)
{

    switch (mainFragment->getRings().size()) {
    case 0:
        return 1;
    case 1:
        return 5;
    case 2:
        return 8;
    case 3:
        return 10;
    default:
        return 12;
    }
}

vector<sketcherMinimizerFragment*> CoordgenFragmenter::findLongestChain(
    const vector<sketcherMinimizerFragment*>& fragments)
{
    vector<sketcherMinimizerFragment*> longestChain;
    for (sketcherMinimizerFragment* fragment : fragments) {
        if (!fragment->isChain) {
            continue;
        }
        int chainN = 0;
        for (sketcherMinimizerBond* b : fragment->_interFragmentBonds) {
            sketcherMinimizerFragment* childFragment =
                (b->getStartAtom()->getFragment() != fragment
                     ? b->getStartAtom()->getFragment()
                     : b->getEndAtom()->getFragment());
            if (childFragment->isChain) {
                ++chainN;
            }
        }
        if (chainN > 1) {
            continue; // it's in the middle of a chain
        }
        queue<sketcherMinimizerFragment*> q;
        map<sketcherMinimizerFragment*, sketcherMinimizerFragment*> parentMap;
        parentMap[fragment] = fragment; // mark first fragment as parent of
                                        // itself so it results visited
        q.push(fragment);
        sketcherMinimizerFragment* lastFragment = nullptr;
        while (!q.empty()) {
            lastFragment = q.front();
            q.pop();
            for (sketcherMinimizerBond* b : lastFragment->_interFragmentBonds) {
                sketcherMinimizerFragment* childFragment =
                    (b->getStartAtom()->getFragment() != lastFragment
                         ? b->getStartAtom()->getFragment()
                         : b->getEndAtom()->getFragment());
                if (parentMap[childFragment] != nullptr) {
                    continue;
                }
                if (!childFragment->isChain) {
                    continue;
                }
                parentMap[childFragment] = lastFragment;
                q.push(childFragment);
            }
        }
        vector<sketcherMinimizerFragment*> chain;
        sketcherMinimizerFragment* recursiveFragment = lastFragment;
        if (recursiveFragment) {
            while (recursiveFragment != fragment) {
                chain.insert(chain.begin(), recursiveFragment);
                recursiveFragment = parentMap[recursiveFragment];
            }
            chain.insert(chain.begin(), recursiveFragment);
        }
        if (chain.size() > longestChain.size()) {
            longestChain = chain;
        }
    }
    return longestChain;
}

void CoordgenFragmenter::addParentRelationsToFragments(
    sketcherMinimizerFragment* mainFragment,
    const vector<sketcherMinimizerFragment*>& fragments)
{
    queue<sketcherMinimizerFragment*> fragmentsQueue;
    fragmentsQueue.push(mainFragment);
    while (!fragmentsQueue.empty()) {
        sketcherMinimizerFragment* fragment = fragmentsQueue.front();
        fragmentsQueue.pop();
        for (sketcherMinimizerBond* bond : fragment->_interFragmentBonds) {
            sketcherMinimizerFragment* childFragment =
                (bond->getStartAtom()->getFragment() == fragment
                     ? bond->getEndAtom()->getFragment()
                     : bond->getStartAtom()->getFragment());
            if (childFragment == fragment->getParent()) {
                continue;
            }
            fragment->_children.push_back(childFragment);
            childFragment->setParent(fragment);
            childFragment->_bondToParent = bond;
            fragmentsQueue.push(childFragment);
        }
    }

    for (sketcherMinimizerFragment* fragment : fragments) {
        /* swap bonds to parent so that startAtom is always the parent's and
         * endAtom always the child's */
        if (fragment->_bondToParent) {
            if (fragment->_bondToParent->getEndAtom()->getFragment() !=
                fragment) {
                sketcherMinimizerAtom* swap = fragment->_bondToParent->endAtom;
                fragment->_bondToParent->endAtom =
                    fragment->_bondToParent->startAtom;
                fragment->_bondToParent->startAtom = swap;
                fragment->_bondToParent->isReversed =
                    !fragment->_bondToParent
                         ->isReversed; // bond stereochemistry
            }
            assert(fragment->_bondToParent->getEndAtom()->getFragment() ==
                   fragment);
        }
    }
}

void CoordgenFragmenter::orderFragments(
    vector<sketcherMinimizerFragment*>& fragments,
    sketcherMinimizerFragment* mainFragment)
{
    queue<sketcherMinimizerFragment*> q;
    vector<sketcherMinimizerFragment*> new_fragments;
    q.push(mainFragment);
    while (!q.empty()) {
        sketcherMinimizerFragment* fragment = q.front();
        q.pop();
        new_fragments.push_back(fragment);
        for (sketcherMinimizerFragment* child : fragment->_children)
            q.push(child);
    }
    assert(fragments.size() == new_fragments.size());
    fragments = new_fragments;
}
