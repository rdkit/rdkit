/*
 Contributors: Nicola Zonta
 Copyright Schrodinger, LLC. All rights reserved
 */

#include "sketcherMinimizerBond.h"
#include "sketcherMinimizerAtom.h"
#include "sketcherMinimizerMaths.h"
#include "sketcherMinimizerMolecule.h"
#include "sketcherMinimizerRing.h"
#include <algorithm>

using namespace std;

sketcherMinimizerAtom* sketcherMinimizerBond::startAtomCIPFirstNeighbor() const
{
    if (bondOrder != 2) {
        return nullptr;
    }
    sketcherMinimizerAtom* a = startAtom;
    if (a->neighbors.size() == 2) {
        if (a->neighbors[0] == endAtom) {
            return a->neighbors[1];
        } else {
            return a->neighbors[0];
        }
    } else if (a->neighbors.size() == 3) {
        std::vector<sketcherMinimizerAtom*> ats;
        for (sketcherMinimizerAtom* n : a->neighbors) {
            if (n != endAtom) {
                ats.push_back(n);
            }
        }
        if (ats.size() == 2) {
            return sketcherMinimizerAtom::CIPPriority(ats[0], ats[1],
                                                      startAtom);
        }
        return nullptr;
    } else {
        return nullptr;
    }
}

sketcherMinimizerAtom* sketcherMinimizerBond::endAtomCIPFirstNeighbor() const
{
    if (bondOrder != 2) {
        return nullptr;
    }
    sketcherMinimizerAtom* a = endAtom;
    if (a->neighbors.size() == 2) {
        if (a->neighbors[0] == startAtom) {
            return a->neighbors[1];
        } else {
            return a->neighbors[0];
        }
    } else if (a->neighbors.size() == 3) {
        std::vector<sketcherMinimizerAtom*> ats;
        for (sketcherMinimizerAtom* n : a->neighbors) {
            if (n != startAtom) {
                ats.push_back(n);
            }
        }
        if (ats.size() == 2) {
            return sketcherMinimizerAtom::CIPPriority(ats[0], ats[1], endAtom);
        }
        return nullptr;
    } else {
        return nullptr;
    }
}

void sketcherMinimizerBond::setAbsoluteStereoFromStereoInfo()
{
    if (isStereo() && m_stereo.atom1 != nullptr && m_stereo.atom2 != nullptr) {
        auto firstCIPNeighborStart = startAtomCIPFirstNeighbor();
        auto firstCIPNeighborEnd = endAtomCIPFirstNeighbor();
        if (firstCIPNeighborStart != nullptr && firstCIPNeighborEnd) {
            bool invert = false;
            if (m_stereo.atom1 != firstCIPNeighborStart &&
                m_stereo.atom1 != firstCIPNeighborEnd) {
                invert = !invert;
            }
            if (m_stereo.atom2 != firstCIPNeighborStart &&
                m_stereo.atom2 != firstCIPNeighborEnd) {
                invert = !invert;
            }
            bool settingIsZ =
                (m_stereo.stereo == sketcherMinimizerBondStereoInfo::cis);
            if (invert) {
                settingIsZ = !settingIsZ;
            }
            isZ = settingIsZ;
        }
    }
    if (m_stereo.stereo == sketcherMinimizerBondStereoInfo::unspecified) {
        m_ignoreZE = true;
    }
}

bool sketcherMinimizerBond::checkStereoChemistry() const
{
    if (!isStereo()) {
        return true;
    }
    if (isInSmallRing()) {
        return true;
    }
    sketcherMinimizerAtom* firstCIPNeighborStart = startAtomCIPFirstNeighbor();
    if (firstCIPNeighborStart == nullptr) {
        return true;
    }
    sketcherMinimizerAtom* firstCIPNeighborEnd = endAtomCIPFirstNeighbor();
    if (firstCIPNeighborEnd == nullptr) {
        return true;
    }
    return (sketcherMinimizerMaths::sameSide(
                firstCIPNeighborStart->getCoordinates(),
                firstCIPNeighborEnd->getCoordinates(),
                getStartAtom()->getCoordinates(),
                getEndAtom()->getCoordinates()) == isZ);
}

bool sketcherMinimizerBond::isInSmallRing() const
{
    for (auto ring : rings) {
        if (!ring->isMacrocycle()) {
            return true;
        }
    }
    return false;
}

bool sketcherMinimizerBond::isInMacrocycle() const
{
    for (auto ring : rings) {
        if (ring->isMacrocycle()) {
            return true;
        }
    }
    return false;
}

bool sketcherMinimizerBond::isTerminal() const
{
    return (getStartAtom()->getBonds().size() == 1 ||
            getEndAtom()->getBonds().size() == 1);
}

bool sketcherMinimizerBond::isInterFragment() const
{
    if (getStartAtom()->getBonds().size() == 1) {
        return false;
    }
    if (getEndAtom()->getBonds().size() == 1) {
        return false;
    }
    if (sketcherMinimizerAtom::shareARing(getStartAtom(), getEndAtom())) {
        return false;
    }
    if (isStereo()) {
        return false;
    }
    return true;
}

bool sketcherMinimizerBond::markedAsCis(sketcherMinimizerAtom* atom1,
                                        sketcherMinimizerAtom* atom2) const
{
    sketcherMinimizerAtom* firstCIPNeighbor1 = startAtomCIPFirstNeighbor();
    sketcherMinimizerAtom* firstCIPNeighbor2 = endAtomCIPFirstNeighbor();
    bool cis = isZ;
    if (atom1 != firstCIPNeighbor1 && atom1 != firstCIPNeighbor2) {
        cis = !cis;
    }
    if (atom2 != firstCIPNeighbor1 && atom2 != firstCIPNeighbor2) {
        cis = !cis;
    }
    return cis;
}

bool sketcherMinimizerBond::isStereo() const
{
    if (bondOrder != 2) {
        return false;
    }
    if (m_ignoreZE) {
        return false;
    }
    sketcherMinimizerRing* ring =
        sketcherMinimizerAtom::shareARing(getStartAtom(), getEndAtom());
    return !(ring && !ring->isMacrocycle());
}

void sketcherMinimizerBond::flip()
{
    size_t totalAtomsNumber = getStartAtom()->getMolecule()->getAtoms().size();
    vector<sketcherMinimizerAtom*> atoms =
        getStartAtom()->getSubmolecule(getEndAtom());
    if (atoms.size() > totalAtomsNumber / 2) {
        atoms = getEndAtom()->getSubmolecule(getStartAtom());
    }
    vector<sketcherMinimizerBond*> allBonds =
        getStartAtom()->getMolecule()->getBonds();

    for (sketcherMinimizerAtom* atom : atoms) {
        sketcherMinimizerAtom::mirrorCoordinates(atom, this);
    }
    for (sketcherMinimizerBond* bond : allBonds) {
        if (find(atoms.begin(), atoms.end(), bond->getStartAtom()) !=
                atoms.end() &&
            find(atoms.begin(), atoms.end(), bond->getEndAtom()) !=
                atoms.end()) {
            bond->isWedge = !bond->isWedge;
        }
    }
}
