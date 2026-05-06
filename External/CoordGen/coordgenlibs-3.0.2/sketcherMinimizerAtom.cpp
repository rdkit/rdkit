/*
 *  sketcherMinimizerAtom.h
 *
 *  Created by Nicola Zonta on 19/10/2010.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#include <algorithm>
#include <numeric>
#include <queue>

#include "sketcherMinimizer.h"
#include "sketcherMinimizerAtom.h"
#include "sketcherMinimizerBond.h"
#include "sketcherMinimizerMaths.h"
#include "sketcherMinimizerRing.h"

using namespace std;

bool CIPAtom::operator<(const CIPAtom& rhs) const
{
    /* check if this has priority over rhs. An atom is better than another if it
       has a parent in the chain that gets priority. Parents are evaluated
       starting from the furthest to the closest.
       Priority is assigned to the atom with highest atomic number, or has been
       found to have priority due to its children in previous iterations of the
       algorithm (scores) or in the present iteration (medals)
    */
    assert(allParents.size() == rhs.allParents.size());
    for (size_t i = 0; i < allParents.size(); i++) {

        if (allParents[i]->atomicNumber > rhs.allParents[i]->atomicNumber) {
            return true;
        }
        if (allParents[i]->atomicNumber < rhs.allParents[i]->atomicNumber) {
            return false;
        }

        if ((*scores)[allParents[i]] < (*rhs.scores)[rhs.allParents[i]]) {
            return true;
        }
        if ((*scores)[allParents[i]] > (*rhs.scores)[rhs.allParents[i]]) {
            return false;
        }

        vector<int> meds = (*medals)[allParents[i]];
        vector<int> meds2 = (*rhs.medals)[rhs.allParents[i]];
        size_t s = (meds.size() < meds2.size()) ? meds.size() : meds2.size();

        for (size_t mm = 0; mm < s; mm++) {
            if (meds[mm] > meds2[mm]) {
                return true;
            }
            if (meds[mm] < meds2[mm]) {
                return false;
            }
        }
        if (meds.size() > meds2.size()) {
            return true;
        }
        if (meds2.size() > meds.size()) {
            return false;
        }
    }
    size_t siz = theseAtoms.size();
    if (rhs.theseAtoms.size() < siz) {
        siz = rhs.theseAtoms.size();
    }
    for (size_t i = 0; i < siz; i++) {
        if (theseAtoms[i].first > rhs.theseAtoms[i].first) {
            return true;
        }
        if (theseAtoms[i].first < rhs.theseAtoms[i].first) {
            return false;
        }
    }
    if (theseAtoms.size() > rhs.theseAtoms.size()) {
        return true;
    }
    if (theseAtoms.size() < rhs.theseAtoms.size()) {
        return false;
    }

    return false;
}

bool CIPAtom::operator==(const CIPAtom& rhs) const
{
    assert(allParents.size() == rhs.allParents.size());
    for (size_t i = 0; i < allParents.size(); i++) {
        if (allParents[i]->atomicNumber != rhs.allParents[i]->atomicNumber) {
            return false;
        }
        if ((*scores)[allParents[i]] != (*rhs.scores)[rhs.allParents[i]]) {
            return false;
        }
    }
    if (theseAtoms.size() != rhs.theseAtoms.size()) {
        return false;
    }
    for (size_t i = 0; i < theseAtoms.size(); i++) {
        if (theseAtoms[i].first != rhs.theseAtoms[i].first) {
            return false;
        }
    }
    return true;
}

std::ostream& operator<<(std::ostream& os, const CIPAtom& a)
{

    for (size_t i = 0; i < a.allParents.size(); i++) {
        os << a.allParents[i]->atomicNumber << "("
           << (*a.scores)[a.allParents[i]] << ")";
        if (!(*a.medals)[a.allParents[i]].empty()) {
            cerr << "<";
            for (int ii : (*a.medals)[a.allParents[i]]) {
                cerr << ii << "|";
            }
            cerr << ">";
        }
        cerr << "   ";
    }
    os << "-";
    for (const auto& theseAtom : a.theseAtoms) {
        os << "    " << theseAtom.first;
    }
    return os;
}

bool CIPAtom::isBetter(CIPAtom& rhs,
                       map<sketcherMinimizerAtom*, unsigned int>* m) const
{
    /*
     Similar to the < operator, but in estabilishing priorities also considers
     scores stored in m
     */
    assert(allParents.size() == rhs.allParents.size());
    for (size_t i = 0; i < allParents.size(); i++) {

        if ((*m)[allParents[i]] > (*m)[rhs.allParents[i]]) {
            return true;
        }
        if ((*m)[allParents[i]] < (*m)[rhs.allParents[i]]) {
            return false;
        }

        if (allParents[i]->atomicNumber > rhs.allParents[i]->atomicNumber) {
            return true;
        }
        if (allParents[i]->atomicNumber < rhs.allParents[i]->atomicNumber) {
            return false;
        }

        if ((*scores)[allParents[i]] < (*rhs.scores)[rhs.allParents[i]]) {
            return true;
        }
        if ((*scores)[allParents[i]] > (*rhs.scores)[rhs.allParents[i]]) {
            return false;
        }

        vector<int> meds = (*medals)[allParents[i]];
        vector<int> meds2 = (*rhs.medals)[rhs.allParents[i]];
        size_t s = (meds.size() < meds2.size()) ? meds.size() : meds2.size();

        for (size_t mm = 0; mm < s; mm++) {
            if (meds[mm] > meds2[mm]) {
                return true;
            }
            if (meds[mm] < meds2[mm]) {
                return false;
            }
        }
        if (meds.size() > meds2.size()) {
            return true;
        }
        if (meds2.size() > meds.size()) {
            return false;
        }
    }
    size_t siz = theseAtoms.size();
    if (rhs.theseAtoms.size() < siz) {
        siz = rhs.theseAtoms.size();
    }
    for (size_t i = 0; i < siz; i++) {
        if (theseAtoms[i].first > rhs.theseAtoms[i].first) {
            return true;
        }
        if (theseAtoms[i].first < rhs.theseAtoms[i].first) {
            return false;
        }
    }
    if (theseAtoms.size() > rhs.theseAtoms.size()) {
        return true;
    }
    if (theseAtoms.size() < rhs.theseAtoms.size()) {
        return false;
    }

    return false;
}

sketcherMinimizerAtom::~sketcherMinimizerAtom() = default;

sketcherMinimizerAtom::sketcherMinimizerAtom()
    : crossLayout(false), fixed(false), constrained(false), rigid(false),
      isSharedAndInner(false), atomicNumber(6), charge(0), _valence(-10),
      _generalUseN(-1), _generalUseN2(-1), m_chmN(-1),
      _generalUseVisited(false), _generalUseVisited2(false), fragment(nullptr),
      needsCheckForClashes(false), visited(false), coordinatesSet(false),
      isR(true), hasStereochemistrySet(false), _hasRingChirality(false)
{
    hidden = false;
    m_pseudoZ = 0.f;
    m_pocketDistance = 0.f;
    m_x3D = INVALID_COORDINATES;
    m_y3D = INVALID_COORDINATES;
    m_z3D = INVALID_COORDINATES;
    m_isClashing = false;
    m_isLigand = false;
    m_isWaterMap = false;
    m_clockwiseInvert = false;
    m_isStereogenic = false;
    m_ignoreRingChirality = false;
}

sketcherMinimizerRing*
sketcherMinimizerAtom::shareARing(const sketcherMinimizerAtom* atom1,
                                  const sketcherMinimizerAtom* atom2)
{
    /* return a ring shared by the two atoms. return a non-macrocycle if
     * possible */
    if (atom1->rings.empty()) {
        return nullptr;
    }
    if (atom2->rings.empty()) {
        return nullptr;
    }

    for (sketcherMinimizerRing* ring : atom1->rings) {
        if (ring->isMacrocycle()) {
            continue;
        }
        for (sketcherMinimizerRing* ring2 : atom2->rings) {
            if (ring == ring2) {
                return ring;
            }
        }
    }
    for (sketcherMinimizerRing* ring : atom1->rings) {
        for (sketcherMinimizerRing* ring2 : atom2->rings) {
            if (ring == ring2) {
                return ring;
            }
        }
    }
    return nullptr;
}

int sketcherMinimizerAtom::findHsNumber() const
{
    int valence = _valence;
    if (valence == -10) {
        valence = expectedValence(atomicNumber); // valence is not yet set
    }
    int nBondOrders = 0;
    for (auto bond : bonds) {
        nBondOrders += bond->bondOrder;
    }
    if (atomicNumber == 16) { // sulfite & sulfate
        int nOs = 0;
        for (size_t i = 0; i < neighbors.size(); ++i) {
            if (neighbors[i]->atomicNumber == 8 && bonds[i]->bondOrder == 2) {
                ++nOs;
            }
        }
        if ((nOs) < 3) {
            valence += nOs * 2;
        }
    }
    if (atomicNumber == 15) { // P
        int nOs = 0;
        for (size_t i = 0; i < neighbors.size(); ++i) {
            if (neighbors[i]->atomicNumber == 8 && bonds[i]->bondOrder == 2) {
                ++nOs;
            }
        }
        if (nOs < 2) {
            valence += nOs * 2;
        }
    }
    int out = valence - nBondOrders + charge;
    if (out < 0) {
        out = 0;
    } else if (out > 4) {
        out = 4;
    }
    return out;
}

unsigned int
sketcherMinimizerAtom::expectedValence(unsigned int atomicNumber) const
{
    switch (atomicNumber) {
    case 1:
        return 1;

    case 5:
        return 3;

    case 6:
        return 4;

    case 7:
        return 3;

    case 8:
        return 2;

    case 9:
        return 1;

    case 14:
        return 4;

    case 15:
        return 3;

    case 16:
        return 2;

    case 17:
        return 1;

    case 34:
        return 2;
    case 35:
        return 1;

    case 53:
        return 1;

    default:
        return 4;
    }
    return 4;
}

vector<sketcherMinimizerAtom*>
sketcherMinimizerAtom::clockwiseOrderedNeighbors() const
{
    vector<pair<float, sketcherMinimizerAtom*>> rankedNeighbors;
    rankedNeighbors.reserve(neighbors.size());
    for (auto&& neighbor : neighbors) {
        float newAngle = sketcherMinimizerMaths::signedAngle(
            neighbors[0]->coordinates, coordinates, neighbor->coordinates);
        if (std::isnan(newAngle)) {
            newAngle = 361;
        } else if (newAngle < 0) {
            newAngle += 360;
        }
        rankedNeighbors.emplace_back(newAngle, neighbor);
    }
    std::sort(rankedNeighbors.begin(), rankedNeighbors.end());
    vector<sketcherMinimizerAtom*> orderedNeighs;
    orderedNeighs.reserve(neighbors.size());
    for (const auto& rankedNeighbor : rankedNeighbors) {
        orderedNeighs.push_back(rankedNeighbor.second);
    }
    return orderedNeighs;
}

void sketcherMinimizerAtom::writeStereoChemistry() // sets stereochemistry for
                                                   // this atom and from
                                                   // hasStereochemistrySet and
                                                   // isR
{

    assert(neighbors.size() == bonds.size());

    if (!hasStereochemistrySet) {
        return;
    }
    if (!canBeChiral()) {

        hasStereochemistrySet = false;
        return;
    } else {

        sketcherMinimizerAtom dummyH;
        dummyH.atomicNumber = 1;
        dummyH.molecule = molecule;

        sketcherMinimizerAtom dummyLP;
        dummyLP.atomicNumber = 0;
        dummyLP.molecule = molecule;

        sketcherMinimizerAtom* dummy = &dummyH;

        vector<sketcherMinimizerAtom*> neighs = neighbors;
        vector<sketcherMinimizerAtom*> orderedNeighs;
        vector<sketcherMinimizerBond*> bbonds = bonds;
        vector<sketcherMinimizerBond*> orderedBonds;
        vector<float> angles;

        int lastPoppedIndex = 0;
        sketcherMinimizerAtom* lastPoppedAtom = neighs[lastPoppedIndex];
        orderedNeighs.push_back(lastPoppedAtom);
        neighs.erase(neighs.begin() + lastPoppedIndex);
        orderedBonds.push_back(bbonds[lastPoppedIndex]);
        bbonds.erase(bbonds.begin() + lastPoppedIndex);

        // TODO: consider using sketcherMinimizerAtom::clockwiseOrderedNeighbors
        while (!neighs.empty()) { // order atoms
            float smallestAngle = 361;
            for (unsigned int i = 0; i < neighs.size(); i++) {
                float newAngle = sketcherMinimizerMaths::signedAngle(
                    sketcherMinimizerPointF(lastPoppedAtom->coordinates.x(),
                                            lastPoppedAtom->coordinates.y()),
                    sketcherMinimizerPointF(coordinates.x(), coordinates.y()),
                    sketcherMinimizerPointF(neighs[i]->coordinates.x(),
                                            neighs[i]->coordinates.y()));
                if (newAngle < 0) {
                    newAngle += 360;
                }
                if (newAngle < smallestAngle) {
                    smallestAngle = newAngle;
                    lastPoppedIndex = i;
                }
            }
            angles.push_back(smallestAngle);
            lastPoppedAtom = neighs[lastPoppedIndex];
            orderedNeighs.push_back(lastPoppedAtom);
            neighs.erase(neighs.begin() + lastPoppedIndex);
            orderedBonds.push_back(bbonds[lastPoppedIndex]);
            bbonds.erase(bbonds.begin() + lastPoppedIndex);
        }
        if ((atomicNumber == 7 || atomicNumber == 16) && _implicitHs == 0 &&
            orderedBonds.size() == 3) {
            dummy = &dummyLP;
        }

        bool four = true;
        float totalAngle = std::accumulate(angles.begin(), angles.end(), 0.f);
        angles.push_back(360.f - totalAngle);

        vector<sketcherMinimizerAtomPriority> atomPriorities,
            orderedAtomPriorities;
        for (auto& orderedNeigh : orderedNeighs) {
            sketcherMinimizerAtomPriority p;
            p.a = orderedNeigh;
            atomPriorities.push_back(p);
        }
        if (atomPriorities.size() == 3) {
            four = false;
            sketcherMinimizerAtomPriority p;
            p.a = dummy;
            atomPriorities.push_back(p);
        }

        bool isStereocenter = setCIPPriorities(atomPriorities, this);

        if (!isStereocenter) {
            if (!m_ignoreRingChirality) {
                _hasRingChirality = true;
                isStereocenter = setCIPPriorities(atomPriorities, this);
            }
        }
        if (!isStereocenter) {
            _hasRingChirality = false;
        }

        orderedAtomPriorities = atomPriorities;
        orderAtomPriorities(orderedAtomPriorities, this);

        if (isStereocenter) {
            if (!four) {
                if (orderedAtomPriorities[0].a == dummy) {
                    orderedAtomPriorities.push_back(orderedAtomPriorities[0]);
                    orderedAtomPriorities.insert(orderedAtomPriorities.begin(),
                                                 orderedAtomPriorities[3]);
                    orderedAtomPriorities.erase(orderedAtomPriorities.begin() +
                                                1);
                    orderedAtomPriorities.erase(orderedAtomPriorities.begin() +
                                                3);
                }
            }

            unsigned int startIndex = 0;
            for (unsigned int i = 0; i < atomPriorities.size(); i++) {
                if (atomPriorities[i].a == orderedAtomPriorities[0].a) {
                    startIndex = i;
                    break;
                }
            }
            for (unsigned int i = 0; i < startIndex; i++) {
                atomPriorities.push_back(atomPriorities[0]);
                atomPriorities.erase(atomPriorities.begin());
            }

            if (four) {
                if (atomPriorities[1].a == orderedAtomPriorities[3].a) {
                    atomPriorities.push_back(atomPriorities[0]);
                    atomPriorities.erase(atomPriorities.begin());
                }
            }
            sketcherMinimizerAtom* mainAtom = nullptr;
            if (four) {
                if (atomPriorities[0].a == orderedAtomPriorities[0].a) {
                    mainAtom = atomPriorities[0].a;
                }
                if (atomPriorities[3].a == orderedAtomPriorities[0].a) {
                    mainAtom = atomPriorities[3].a;
                }
            }
            bool invert = false;
            sketcherMinimizerBond* b1 = bondTo(atomPriorities[0].a);
            if (b1) {
                if (!four || mainAtom == atomPriorities[0].a ||
                    (!sketcherMinimizer::sameRing(this, atomPriorities[0].a) &&
                     !atomPriorities[0].a->hasStereochemistrySet)) {
                    bool reverse = false;
                    if (b1->startAtom != this) {
                        reverse = true;
                    }
                    b1->isWedge = !invert;
                    b1->hasStereochemistryDisplay = true;
                    b1->isReversed = reverse;
                } else {
                    b1->hasStereochemistryDisplay = false;
                }
            }
            if (four) {
                sketcherMinimizerBond* b2 = bondTo(atomPriorities[3].a);
                if (b2) {
                    if (mainAtom == atomPriorities[3].a ||
                        (!sketcherMinimizer::sameRing(this,
                                                      atomPriorities[3].a) &&
                         !atomPriorities[3].a->hasStereochemistrySet)) {
                        bool reverse = false;
                        if (b2->startAtom != this) {
                            reverse = true;
                        }
                        b2->isWedge = invert;
                        b2->hasStereochemistryDisplay = true;
                        b2->isReversed = reverse;
                    } else {
                        b2->hasStereochemistryDisplay = false;
                    }
                }
            }

            sketcherMinimizerBond* b3 = bondTo(atomPriorities[1].a);
            if (b3) {

                b3->hasStereochemistryDisplay = false;
            }
            sketcherMinimizerBond* b4 = bondTo(atomPriorities[2].a);
            if (b4) {

                b4->hasStereochemistryDisplay = false;
            }

            int readS = readStereochemistry(true);
            readS = -readS; // inverting stereochemistry cause isR is in
                            // sketcher coords and readStereo in coordgen coords
                            // (flipped y)

            if ((readS == -1 && isR) || (readS == 1 && !isR)) { // inverting
                invert = true;

                if (b1) {
                    if (!four || mainAtom == atomPriorities[0].a ||
                        (!sketcherMinimizer::sameRing(this,
                                                      atomPriorities[0].a) &&
                         !atomPriorities[0].a->hasStereochemistrySet)) {
                        //          cerr << "and setting it "<<endl;
                        bool reverse = false;
                        if (b1->startAtom != this) {
                            reverse = true;
                        }
                        b1->isWedge = !invert;
                        b1->hasStereochemistryDisplay = true;
                        b1->isReversed = reverse;
                    } else {
                        b1->hasStereochemistryDisplay = false;
                    }
                }

                if (four) {
                    sketcherMinimizerBond* b2 = bondTo(atomPriorities[3].a);

                    if (b2) {
                        if (mainAtom == atomPriorities[3].a ||
                            (!sketcherMinimizer::sameRing(
                                 this, atomPriorities[3].a) &&
                             !atomPriorities[3].a->hasStereochemistrySet)) {
                            bool reverse = false;
                            if (b2->startAtom != this) {
                                reverse = true;
                            }
                            b2->isWedge = invert;
                            b2->hasStereochemistryDisplay = true;
                            b2->isReversed = reverse;
                        } else {
                            b2->hasStereochemistryDisplay = false;
                        }
                    }
                }
            }

        }

        else {
            for (auto b : bonds) {
                b->hasStereochemistryDisplay = false;
            }
        }
    }
}

sketcherMinimizerAtomChiralityInfo::sketcherMinimizerChirality
sketcherMinimizerAtom::getRelativeStereo(sketcherMinimizerAtom* lookingFrom,
                                         sketcherMinimizerAtom* atom1,
                                         sketcherMinimizerAtom* atom2)
{
    readStereochemistry(); // to set m_RSPriorities
    auto RSpriorities = m_RSPriorities;
    if (RSpriorities.size() < 3) {
        return sketcherMinimizerAtomChiralityInfo::unspecified;
    }
    vector<int> priorities(4, 3);

    /* order the CIP priority of the atoms in the following order
     atom1 - atom2 - atom3 - atomLookingFrom
     */
    for (unsigned int nn = 0; nn < neighbors.size(); nn++) {
        sketcherMinimizerAtom* n = neighbors[nn];
        if (n == atom1) {
            priorities[0] = RSpriorities[nn];

        } else if (n == atom2) {
            priorities[1] = RSpriorities[nn];
        } else if (n == lookingFrom) {
            priorities[3] = RSpriorities[nn];
        } else {
            priorities[2] = RSpriorities[nn];
        }
    }
    vector<int> can(4);
    for (unsigned int i = 0; i < 4; i++) {
        can[i] = i;
    }
    /*
     this represents a molecule with
     atom1 (priority 0 - highest)
     atom2 (priority 1)
     atom3 (priority 2)
     atomLookingFrom (priority 3 -lowest)
     which is the opposite of the CIP rules, where the the lowest priority atom
     is AWAY from the observer. This is the reason why we return CCW for R and
     CW for S.
     */
    bool match = sketcherMinimizerAtom::matchCIPSequence(priorities, can);
    bool isClockWise = (match ? !isR : isR);
    if (isClockWise) {
        return sketcherMinimizerAtomChiralityInfo::clockwise;
    }
    return sketcherMinimizerAtomChiralityInfo::counterClockwise;
}

bool sketcherMinimizerAtom::setAbsoluteStereoFromChiralityInfo()
{
    auto info = m_chiralityInfo;
    if (info.direction == sketcherMinimizerAtomChiralityInfo::unspecified) {
        return true;
    }
    readStereochemistry(); // to set m_RSPriorities
    auto RSpriorities = m_RSPriorities;
    if (RSpriorities.size() < 3) {
        cerr << "CHMMol-> sketcher stereo error: wrong number for RSpriorities"
             << endl;
        return false;
    }

    vector<int> priorities(4, 5);

    bool at3 = false;
    for (unsigned int nn = 0; nn < neighbors.size(); nn++) {
        sketcherMinimizerAtom* n = neighbors[nn];
        if (n == info.atom1) {
            priorities[0] = RSpriorities[nn];

        } else if (n == info.atom2) {

            priorities[1] = RSpriorities[nn];
        } else if (n == info.lookingFrom) {
            priorities[3] = RSpriorities[nn];
        } else {
            if (at3) {
                cerr << "CHMMol-> sketcher stereo error: more than 1 atom not "
                        "matching"
                     << endl;
                return false;

            } else {
                at3 = true;
                priorities[2] = RSpriorities[nn];
            }
        }
    }
    int addingHN = 0;
    if (priorities[0] == 5) {
        priorities[0] = 3;
        addingHN++;
    }
    if (priorities[1] == 5) {
        priorities[1] = 3;
        addingHN++;
    }
    if (priorities[2] == 5) {
        priorities[2] = 3;
        addingHN++;
    }
    if (priorities[3] == 5) {
        priorities[3] = 3;
        addingHN++;
    }
    if (addingHN > 1) {
        cerr << "CHMMol-> sketcher stereo error: more than 1 H on chiral center"
             << endl;
        return false;
    }

    bool invert = false;

    vector<int> can(4);
    for (unsigned int i = 0; i < 4; i++) {
        can[i] = i;
    }
    if (!sketcherMinimizerAtom::matchCIPSequence(priorities, can)) {
        invert = !invert;
    }
    bool isRBool = true;
    if (info.direction == sketcherMinimizerAtomChiralityInfo::clockwise) {
        isRBool = false;
    }
    if (invert) {
        isRBool = !isRBool;
    }
    isR = isRBool;
    hasStereochemistrySet = true;
    return true;
}

bool sketcherMinimizerAtom::matchCIPSequence(vector<int>& v1, vector<int>& v2)

{
    if (v1.size() < v2.size()) {
        v1.push_back(3);
    } else if (v2.size() < v1.size()) {
        v2.push_back(3);
    }

    int outofPlaceNs = 0;
    for (unsigned int i = 0; i < v1.size(); i++) {
        if (v1[i] != v2[i]) {
            outofPlaceNs++;
        }
    }
    if (outofPlaceNs == 2) {
        return false;
    } else if (outofPlaceNs == 4) {
        int n1 = v1[0];
        int index2 = 0;
        for (unsigned int j = 0; j < v2.size(); j++) {
            if (v2[j] == n1) {
                index2 = j;
                break;
            }
        }
        if (v1[index2] != v2[0]) {
            return false;
        }
    }
    return true;
}

void sketcherMinimizerAtom::setCoordinates(sketcherMinimizerPointF coords)
{
    coordinates = std::move(coords);
    coordinates.round();
    coordinatesSet = true;
}

bool sketcherMinimizerAtom::hasNoStereoActiveBonds() const
{
    for (auto bond : bonds) {
        if (bond->isStereo()) {
            return false;
        }
    }
    return true;
}

void sketcherMinimizerAtom::orderAtomPriorities(
    vector<sketcherMinimizerAtomPriority>& atomPriorities,
    sketcherMinimizerAtom* center) // orders trying to keep long chains in
                                   // position 2 and 3 and side substituents in
                                   // 1 and 4
{
    assert(atomPriorities.size() == 4);
    vector<float> weights(4);
    for (unsigned int i = 0; i < 4; i++) {
        queue<sketcherMinimizerAtom*> q;
        for (auto& _atom : center->molecule->_atoms) {
            _atom->_generalUseVisited = false;
        }

        q.push(atomPriorities[i].a);

        center->_generalUseVisited = true;
        atomPriorities[i].a->_generalUseVisited = true;
        int counter = 0;
        while (!q.empty()) {
            counter++;
            sketcherMinimizerAtom* at = q.front();
            q.pop();
            for (auto n : at->neighbors) {
                if (!n->_generalUseVisited) {
                    q.push(n);
                    n->_generalUseVisited = true;
                }
            }
        }
        weights[i] = static_cast<float>(counter);
        sketcherMinimizerBond* b = center->bondTo(atomPriorities[i].a);
        if (b) {
            if (b->bondOrder == 2) {
                weights[i] -=
                    0.25; // so that =O get lower priority than -OH in phosphate
            }
            if (center->atomicNumber == 16 && b->bondOrder == 2) {
                weights[i] += 2000; // forcing the wedge away from double bond
            }
            // in sulphoxide

            if (sketcherMinimizer::sameRing(b->startAtom, b->endAtom)) {
                weights[i] +=
                    500; // force same ring atoms to be in position 3 and 4
            }
        }
        if (atomPriorities[i].a->atomicNumber == 6) {
            weights[i] += 0.5; // carbons get priority over other heavy atoms
        }
        if (atomPriorities[i].a->atomicNumber == 1) {
            weights[i] -= 0.5;
        }
        if (atomPriorities[i].a->isSharedAndInner &&
            !center->isSharedAndInner) {
            weights[i] -= 2000; // forced bond to shared and inner
        }

        if (center->crossLayout) {
            if (atomPriorities[i].a->neighbors.size() > 1) {
                weights[i] += 200;
            }
        }
        if (/* atomPriorities[i].a->isStereogenic && */ atomPriorities[i]
                .a->hasStereochemistrySet) {
            weights[i] += 10000; // to avoid problems with wedges when 2
        }
        // stereocenters are near
        for (unsigned int j = 0; j < atomPriorities[i].a->bonds.size(); j++) {
            if (atomPriorities[i].a->bonds[j]->bondOrder == 2) {
                weights[i] += 100;
                break;
            }
        }
    }

    float lowestWeight = weights[0];
    int index = 0;
    sketcherMinimizerAtomPriority firstAtom = atomPriorities[0];
    for (unsigned int i = 1; i < 4; i++) {
        if (weights[i] < lowestWeight) {
            lowestWeight = weights[i];
            index = i;
        }
    }
    firstAtom = atomPriorities[index];
    atomPriorities.erase(atomPriorities.begin() + index);
    weights.erase(weights.begin() + index);

    index = 0;
    lowestWeight = weights[0];
    sketcherMinimizerAtomPriority secondAtom = atomPriorities[0];
    for (unsigned int i = 1; i < 3; i++) {
        if (weights[i] < lowestWeight) {
            lowestWeight = weights[i];
            index = i;
        }
    }

    secondAtom = atomPriorities[index];
    atomPriorities.erase(atomPriorities.begin() + index);

    if ((center->atomicNumber != 16 && center->atomicNumber != 15) ||
        center->neighbors.size() != 4) {
        atomPriorities.push_back(secondAtom);
        atomPriorities.insert(atomPriorities.begin(), firstAtom);
    } else {
        atomPriorities.insert(atomPriorities.begin() + 1, secondAtom);
        atomPriorities.insert(atomPriorities.begin(), firstAtom);
    }
}

bool sketcherMinimizerAtom::setCIPPriorities(
    vector<sketcherMinimizerAtomPriority>& atomPriorities,
    sketcherMinimizerAtom* center)
{
    for (auto& atomPrioritie : atomPriorities) {
        atomPrioritie.priority = 3;
    }
    if (atomPriorities.size() != 4) {
        //    cerr << "coordgen: stereo error. (wrong number of atom priorities:
        //    "<< atomPriorities.size () << ")"<<endl; // commented for
        //    Ev:134037
        return false;
    }
    for (unsigned int i = 0; i < atomPriorities.size() - 1; i++) {
        for (unsigned int j = i + 1; j < atomPriorities.size(); j++) {

            sketcherMinimizerAtom* res =
                CIPPriority(atomPriorities[i].a, atomPriorities[j].a, center);

            if (res == atomPriorities[i].a) {
                atomPriorities[i].priority--;
            } else if (res == atomPriorities[j].a) {
                atomPriorities[j].priority--;
            }
        }
    }
    vector<bool> found(4, false);

    for (auto& atomPrioritie : atomPriorities) {
        if (found[atomPrioritie.priority]) {
            return false; // two atoms have the same priority
        }
        found[atomPrioritie.priority] = true;
    }
    return true;
}

sketcherMinimizerAtom*
sketcherMinimizerAtom::CIPPriority(sketcherMinimizerAtom* at1,
                                   sketcherMinimizerAtom* at2,
                                   sketcherMinimizerAtom* center)
{
    assert(center->molecule);
    assert(at1->molecule == center->molecule);
    assert(at2->molecule == center->molecule);
    assert(!center->molecule->_atoms.empty());
    assert(at1);
    assert(at2);

    if (at1->atomicNumber > at2->atomicNumber) {
        return at1;
    } else if (at2->atomicNumber > at1->atomicNumber) {
        return at2;
    }
    if (center->_hasRingChirality && !center->m_ignoreRingChirality) {
        if (sketcherMinimizer::sameRing(center, at1, at2)) {
            if (at1->_generalUseN > at2->_generalUseN) {
                return at1;
            } else {
                return at2;
            }
        }
    }

    vector<CIPAtom> AN1, AN2;

    map<sketcherMinimizerAtom*, int> score1,
        score2; // used to keep track if a parent atom has been found to have
                // priority over another
    map<sketcherMinimizerAtom*, vector<int>> medals1,
        medals2; // marks if an atom is a parent of the atoms being evaluated in
                 // the current iteration
    map<sketcherMinimizerAtom*, int> visited1,
        visited2; // marks at which iteration this atom was evaluated

    visited1[center] = 1;
    visited2[center] = 1;
    visited1[at1] = 2;
    visited2[at2] = 2;

    vector<pair<int, sketcherMinimizerAtom*>> v1, v2;
    v1.emplace_back(at1->atomicNumber, at1);
    v2.emplace_back(at2->atomicNumber, at2);

    vector<sketcherMinimizerAtom*> parents1, parents2;
    parents1.push_back(center);
    parents1.push_back(at1);
    parents2.push_back(center);
    parents2.push_back(at2);

    AN1.emplace_back(v1, center, parents1, &score1, &medals1, &visited1);
    AN2.emplace_back(v2, center, parents2, &score2, &medals2, &visited2);

    while (!AN1.empty() || !AN2.empty()) {
        stable_sort(AN1.begin(), AN1.end());
        stable_sort(AN2.begin(), AN2.end());

        sketcherMinimizerAtom::assignMedals(AN1);
        sketcherMinimizerAtom::assignMedals(AN2);

        sketcherMinimizerAtom::chooseFirstAndSortAccordingly(AN1);
        sketcherMinimizerAtom::chooseFirstAndSortAccordingly(AN2);

        sketcherMinimizerAtom::finalizeScores(AN1);
        sketcherMinimizerAtom::finalizeScores(AN2);

        size_t nn = AN1.size();
        if (AN2.size() < nn) {
            nn = AN2.size();
        }

        for (size_t i = 0; i < nn; ++i) {
            if (AN1[i] < AN2[i]) {
                return at1;
            }
            if (AN2[i] < AN1[i]) {
                return at2;
            }
        }

        if (AN1.size() > AN2.size()) {
            return at1;
        }
        if (AN2.size() > AN1.size()) {
            return at2;
        }

        AN1 = sketcherMinimizerAtom::expandOneLevel(AN1);
        AN2 = sketcherMinimizerAtom::expandOneLevel(AN2);
    }
    return nullptr;
}

void sketcherMinimizerAtom::chooseFirstAndSortAccordingly(vector<CIPAtom>& V)
{
    if (V.size() < 2) {
        return;
    }
    vector<CIPAtom> copyV = V;
    V.clear();
    map<sketcherMinimizerAtom*, unsigned int> friendsMask;
    while (!copyV.empty()) {
        int bestI = 0;
        for (unsigned int i = 1; i < copyV.size(); i++) {
            if (copyV[i].isBetter(copyV[bestI], &friendsMask)) {
                bestI = i;
            }
        }
        CIPAtom newBest = copyV[bestI];
        copyV.erase(copyV.begin() + bestI);
        V.push_back(newBest);

        for (auto allParent : newBest.allParents) {
            friendsMask[allParent] |= (1 << copyV.size());
        }
    }
}

vector<CIPAtom> sketcherMinimizerAtom::expandOneLevel(vector<CIPAtom>& oldV)
{
    // we need to keep the bound atoms together, one by one will not work.

    map<sketcherMinimizerAtom*, bool> visitedThisRound;
    vector<CIPAtom> newV;
    for (auto& an : oldV) {
        for (unsigned int aa = 0; aa < an.theseAtoms.size(); aa++) {
            sketcherMinimizerAtom* a = an.theseAtoms[aa].second;
            if (a == nullptr) {
                continue; // dummy atom
            }
            //    if (visitedThisRound[a]) continue; // a is present twice
            //    because closing a ring and has already been dealt with
            visitedThisRound[a] = true;
            map<sketcherMinimizerAtom*, int>* visited = an.visited;
            map<sketcherMinimizerAtom*, vector<int>>* medals = an.medals;
            map<sketcherMinimizerAtom*, int>* scores = an.scores;

            vector<sketcherMinimizerAtom*> allParents = an.allParents;
            allParents.push_back(a);
            vector<pair<int, sketcherMinimizerAtom*>> theseAts;

            for (unsigned int n = 0; n < a->bonds.size(); n++) {
                sketcherMinimizerAtom* neigh = a->neighbors[n];
                if (neigh != an.parent) { // if n is not the direct parent of a
                    bool ghost =
                        (neigh->atomicNumber == 1 ||
                         ((*visited)[neigh] &&
                          (*visited)[neigh] !=
                              (*visited)[a] + 1) // closing a ring to an atom
                                                 // already visited in a
                                                 // previous cycle
                        );
                    theseAts.emplace_back(
                        neigh->atomicNumber,
                        ghost ? ((sketcherMinimizerAtom*) nullptr)
                              : neigh); // put a ghost for hydrogens and atoms
                                        // closing a ring
                    if (!ghost) {
                        (*visited)[neigh] = (*visited)[a] + 1;
                    }
                }
                if (a->bonds[n]->bondOrder == 2) { // put ghosts for multiple
                                                   // order bonds, even to the
                                                   // parent
                    theseAts.emplace_back(neigh->atomicNumber,
                                          (sketcherMinimizerAtom*) nullptr);
                }

                else if (a->bonds[n]->bondOrder == 3) {
                    theseAts.emplace_back(neigh->atomicNumber,
                                          (sketcherMinimizerAtom*) nullptr);
                    theseAts.emplace_back(neigh->atomicNumber,
                                          (sketcherMinimizerAtom*) nullptr);
                }
            }

            for (int counter = 0; counter < a->_implicitHs;
                 counter++) { // put ghosts for implicit Hs
                theseAts.emplace_back(1, (sketcherMinimizerAtom*) nullptr);
            }
            stable_sort(theseAts.begin(), theseAts.end());
            reverse(theseAts.begin(), theseAts.end());
            newV.emplace_back(theseAts, a, allParents, scores, medals, visited);
        }
    }
    return newV;
}

void sketcherMinimizerAtom::assignMedals(vector<CIPAtom>& v)
{

    if (v.empty()) {
        return;
    }
    map<sketcherMinimizerAtom*, vector<int>>* medals = v[0].medals;

    vector<bool> isEqualToPrevious(v.size());
    for (unsigned int i = 1; i < v.size();
         i++) { // need to be done before assigning the medals because they are
                // considered in the == operator

        isEqualToPrevious[i] = (v[i] == v[i - 1]);
    }
    unsigned int medalLvl = 0;
    for (unsigned int i = 0; i < v.size(); i++) {
        if (i > 0) {
            if (isEqualToPrevious[i]) {
                assert(medalLvl > 0);
                medalLvl--;
            }
        }
        for (auto allParent : v[i].allParents) {
            vector<int> medalsV = (*medals)[allParent];
            while (medalsV.size() < medalLvl) {
                medalsV.push_back(0);
            }
            if (medalsV.size() > medalLvl) {
                medalsV[medalLvl] = medalsV[medalLvl] + 1;
            } else {
                medalsV.push_back(1);
            }
            (*medals)[allParent] = medalsV;
        }
        medalLvl++;
    }
}

void sketcherMinimizerAtom::finalizeScores(vector<CIPAtom>& v)
{

    if (v.empty()) {
        return;
    }
    vector<bool> isEqualToPrevious(v.size());
    for (unsigned int i = 1; i < v.size();
         i++) { // need to be done before assigning the scores because they are
                // considered in the == operator

        isEqualToPrevious[i] = (v[i] == v[i - 1]);
    }
    map<sketcherMinimizerAtom*, vector<int>>* medals = v[0].medals;
    map<sketcherMinimizerAtom*, int>* scores = v[0].scores;
    scores->clear();
    int score = 1;
    for (unsigned int i = 0; i < v.size(); i++) {
        if (i > 0) {
            if (isEqualToPrevious[i]) {
                score--;
            }
        }
        /* write the score */
        for (auto allParent : v[i].allParents) {
            if ((*scores)[allParent] == 0) {
                (*scores)[allParent] = score;
            }
        }
        score++;
    }
    medals->clear();
}

sketcherMinimizerBond*
sketcherMinimizerAtom::bondTo(sketcherMinimizerAtom* at) const
{
    for (unsigned int i = 0; i < neighbors.size(); i++) {
        if (neighbors[i] == at) {
            return bonds[i];
        }
    }
    return nullptr;
}

bool sketcherMinimizerAtom::isResidue() const
{
    return false;
}

int sketcherMinimizerAtom::readStereochemistry(
    bool readOnly) // 0 if not assigned, 1 if R, -1 if S
{
    if (!readOnly) {
        _hasRingChirality = false;
        m_isStereogenic = false;
    }
    if (!canBeChiral()) {
        return 0;
    }
    //    if (neighbors.size () != 4 && neighbors.size () != 3) return 0;

    sketcherMinimizerAtom dummyH;
    dummyH.atomicNumber = 1;
    dummyH.molecule = molecule;
    sketcherMinimizerAtom dummyLP;
    dummyLP.atomicNumber = 0;
    dummyLP.molecule = molecule;

    vector<sketcherMinimizerAtom*> neighs = neighbors;
    vector<sketcherMinimizerAtom*> orderedNeighs;
    vector<sketcherMinimizerBond*> bnds = bonds;
    vector<sketcherMinimizerBond*> orderedBonds;
    vector<float> angles;
    int lastPoppedIndex = 0;

    sketcherMinimizerAtom* lastPoppedAtom = neighs[lastPoppedIndex];
    orderedNeighs.push_back(lastPoppedAtom);
    neighs.erase(neighs.begin() + lastPoppedIndex);
    orderedBonds.push_back(bnds[lastPoppedIndex]);
    bnds.erase(bnds.begin() + lastPoppedIndex);

    // TODO: consider using sketcherMinimizerAtom::clockwiseOrderedNeighbors
    while (!neighs.empty()) { // order atoms
        float smallestAngle = 361;
        for (unsigned int i = 0; i < neighs.size(); i++) {
            float newAngle = sketcherMinimizerMaths::signedAngle(
                lastPoppedAtom->coordinates, coordinates,
                neighs[i]->coordinates);
            if (newAngle < 0) {
                newAngle += 360;
            }
            if (newAngle < smallestAngle) {
                smallestAngle = newAngle;
                lastPoppedIndex = i;
            }
        }
        angles.push_back(smallestAngle);
        lastPoppedAtom = neighs[lastPoppedIndex];
        orderedNeighs.push_back(lastPoppedAtom);
        neighs.erase(neighs.begin() + lastPoppedIndex);
        orderedBonds.push_back(bnds[lastPoppedIndex]);
        bnds.erase(bnds.begin() + lastPoppedIndex);
    }

    float totalAngle = 0;
    for (float angle : angles) {
        totalAngle += angle;
    }
    angles.push_back(360.f - totalAngle);

    bool semiplane = false;
    sketcherMinimizerAtom* centralAtom = nullptr;
    if (angles.size() == 3) {
        for (unsigned int i = 0; i < angles.size(); i++) {
            if (angles[i] > 180.f) {
                semiplane = true;
                size_t precI;
                if (i == 0) {
                    precI = angles.size() - 1;
                } else {
                    precI = i - 1;
                }
                centralAtom = orderedNeighs[precI];
            }
        }
    }

    int wedgeN = -1;
    int dashedN = -1;
    bool giveUp = false;
    assert(orderedBonds.size() == 4 || orderedBonds.size() == 3);
    for (unsigned int i = 0; i < orderedBonds.size(); i++) {
        if (orderedBonds[i]->hasStereochemistryDisplay) {
            if (!orderedBonds[i]->isReversed &&
                orderedBonds[i]->startAtom != this) {
                continue;
            }
            if (orderedBonds[i]->isReversed &&
                orderedBonds[i]->endAtom != this) {
                continue;
            }

            bool wedge = orderedBonds[i]->isWedge;

            if (orderedBonds[i]->isReversed) {
                wedge = !wedge;
            }

            if (orderedBonds[i]->startAtom != this) {
                wedge = !wedge;
            }

            if (wedge) {
                if (wedgeN == -1) {
                    wedgeN = i;
                } else {
                    giveUp = true;
                }
            } else {
                if (dashedN == -1) {
                    dashedN = i;
                } else {
                    giveUp = true;
                }
            }
        }
    }

    int startIndex = 0;
    bool invert = false;
    if (dashedN == -1 && wedgeN == -1) {
        giveUp = true;
    } else if (dashedN == -1) {
        startIndex = wedgeN;

    } else if (wedgeN == -1) {
        startIndex = dashedN;
        invert = true;

    } else {
        if (orderedBonds.size() == 3) {
            return 0;
        }
        if (wedgeN - dashedN == 1 || wedgeN - dashedN == -3) {
            startIndex = wedgeN;

        } else if (wedgeN - dashedN == -1 || wedgeN - dashedN == 3) {
            startIndex = dashedN;
            invert = true;
        }

        else {
            giveUp = true;
        }
    }

    int nOfHs = _implicitHs;
    size_t totalSubstituentsN = neighbors.size() + nOfHs;

    if (orderedNeighs.size() == 3 && nOfHs == 1) {
        if (semiplane) {
            if (centralAtom == orderedNeighs[startIndex]) {
                invert = !invert;
            }
        }

        orderedNeighs.insert(orderedNeighs.begin() + startIndex, &dummyH);
        invert = !invert;
    }
    if (orderedNeighs.size() == 3 &&
        (atomicNumber == 7 ||
         atomicNumber == 16)) { // nitrogen or sulphoxide // should we check for
                                // the =O on the sulphoxide?

        orderedNeighs.insert(orderedNeighs.begin() + startIndex, &dummyLP);
        invert = !invert;
    }

    vector<sketcherMinimizerAtomPriority> atomPriorities;
    for (auto& orderedNeigh : orderedNeighs) {
        sketcherMinimizerAtomPriority p;
        p.a = orderedNeigh;
        atomPriorities.push_back(p);
    }
    if (atomPriorities.size() != 4) {
        return 0;
    }

    vector<int> canonical;
    canonical.push_back(0);
    canonical.push_back(1);
    canonical.push_back(2);
    canonical.push_back(3);

    m_RSPriorities.clear();
    bool isStereocenter = setCIPPriorities(atomPriorities, this);

    if (!isStereocenter) {
        if (!readOnly) {
            if (!m_ignoreRingChirality) {
                _hasRingChirality = true;
                isStereocenter = setCIPPriorities(atomPriorities, this);
            }
        }
    }

    if (!isStereocenter) {
        if (!readOnly) {
            _hasRingChirality = false;
        }
    }

    if (!isStereocenter) {
        giveUp = true;
    } else {
        if (!readOnly) {
            m_isStereogenic = true;
        }
    }
    if (totalSubstituentsN < 4 && (atomicNumber != 7 && atomicNumber != 16)) {
        if (!readOnly) {
            m_isStereogenic = false;
        }
    }

    if (!m_isStereogenic) {
        if (!readOnly) {
            _hasRingChirality = false;
        }
    }
    for (auto n : neighbors) {
        for (auto& atomPrioritie : atomPriorities) {
            if (atomPrioritie.a == n) {
                m_RSPriorities.push_back(atomPrioritie.priority);
                break;
            }
        }
    }

    if (!matchCIPSequence(canonical, m_RSPriorities)) {
        m_clockwiseInvert = true;
    } else {
        m_clockwiseInvert = false;
    }

    if (!giveUp) {
        int outofPlaceAtoms = 0;
        for (unsigned int i = 0; i < atomPriorities.size(); i++) {
            //  cerr <<i<<"    "<<atomPriorities[i].a->atomicNumber<<endl;
            int n = startIndex + i;
            if (n > 3) {
                n -= 4;
            }
            if (atomPriorities[n].priority != i) {
                outofPlaceAtoms++;
            }
        }
        if (outofPlaceAtoms == 2) {
            invert = !invert;
        } else if (outofPlaceAtoms == 4) {
            int n2 = atomPriorities[startIndex].priority + startIndex;
            if (n2 > 3) {
                n2 -= 4;
            }
            if (atomPriorities[n2].priority != 0) {
                invert = !invert; // if I swapped position 0 with another will
            }
            // that settle both atoms?
        }

        if (invert) {
            return -1;
        } else {
            return 1;
        }
    }
    return 0;
}

bool sketcherMinimizerAtom::canBeChiral() const
{
    if (atomicNumber == 16) {
        if (neighbors.size() == 3) {
            return true;
        }
    }
    if (atomicNumber == 7) {
        if (neighbors.size() == 3 || neighbors.size() == 4) {
            return true;
        }
    }
    if (neighbors.size() != 3 && neighbors.size() != 4) {
        return false;
    }
    if ((neighbors.size() + _implicitHs != 4)) {
        return false;
    }

    return true;
}

/* get a vector that points out towards the most free region of space around the
 * atom. Used to determined where an arrow to the atom will point from */
sketcherMinimizerPointF sketcherMinimizerAtom::getSingleAdditionVector(
    const vector<sketcherMinimizerAtom*>& ats)
{
    sketcherMinimizerPointF out(0.f, 0.f);
    int n = 0;
    for (unsigned int i = 0; i < ats.size(); i++) {
        sketcherMinimizerAtom* at = ats[i];
        for (unsigned int j = 0; j < at->neighbors.size(); j++) {
            sketcherMinimizerAtom* ne = at->neighbors[j];
            if (ne->neighbors.size() > 1 &&
                find(ats.begin(), ats.end(), ne) == ats.end()) {
                out += ne->coordinates - at->coordinates;
                n++;
            }
        }
    }
    if (n > 0) {
        out /= n;
    } else {
        return sketcherMinimizerPointF(50, 0);
    }
    out *= -1;
    return out;
}

sketcherMinimizerPointF sketcherMinimizerAtom::getSingleAdditionVector() const
{
    sketcherMinimizerPointF out(0.f, 0.f);
    float totalf = 0.f;
    if (!neighbors.empty()) {
        for (auto neighbor : neighbors) {
            float f = 1.f;
            if (sketcherMinimizer::sameRing(this, neighbor)) {
                f = 4.f;
            }
            out += f * (neighbor->coordinates - coordinates);
            totalf += f;
        }
        out /= totalf;
    }
    out *= -1;
    return out;
}

bool sketcherMinimizerAtom::isMetal(const unsigned int atomicNumber)
{
    if (atomicNumber >= 3 && atomicNumber <= 4) {
        return true;
    }
    if (atomicNumber >= 11 && atomicNumber <= 12) {
        return true;
    }
    if (atomicNumber >= 19 && atomicNumber <= 20) {
        return true;
    }
    if (atomicNumber >= 37 && atomicNumber <= 38) {
        return true;
    }
    if (atomicNumber >= 55 && atomicNumber <= 56) {
        return true;
    }
    if (atomicNumber >= 87 && atomicNumber <= 88) {
        return true;
    }

    if (atomicNumber >= 21 && atomicNumber <= 30) {
        return true;
    }
    if (atomicNumber >= 39 && atomicNumber <= 48) {
        return true;
    }
    if (atomicNumber >= 72 && atomicNumber <= 80) {
        return true;
    }
    if (atomicNumber >= 104 && atomicNumber <= 112) {
        return true;
    }
    if (atomicNumber >= 57 && atomicNumber <= 71) {
        return true; // lanthanoids
    }
    if (atomicNumber >= 89 && atomicNumber <= 103) {
        return true; // actinoids
    }
    if (atomicNumber == 13) {
        return true;
    }
    if (atomicNumber == 31) {
        return true;
    }
    if (atomicNumber == 49) {
        return true;
    }
    if (atomicNumber == 81) {
        return true;
    }
    if (atomicNumber == 32) {
        return true;
    }
    if (atomicNumber == 50) {
        return true;
    }
    if (atomicNumber == 82) {
        return true;
    }
    if (atomicNumber == 51) {
        return true;
    }
    if (atomicNumber == 83) {
        return true;
    }
    if (atomicNumber == 84) {
        return true;
    }
    return false;
}

void sketcherMinimizerAtom::mirrorCoordinates(sketcherMinimizerAtom* at,
                                              const sketcherMinimizerBond* bond)
{
    at->setCoordinates(sketcherMinimizerMaths::mirrorPoint(
        at->getCoordinates(), bond->getStartAtom()->getCoordinates(),
        bond->getEndAtom()->getCoordinates()));
}

bool sketcherMinimizerAtom::hasValid3DCoordinates() const
{
    return (m_x3D < INVALID_COORDINATES && m_y3D < INVALID_COORDINATES &&
            m_z3D < INVALID_COORDINATES);
}

vector<sketcherMinimizerAtom*>
sketcherMinimizerAtom::getSubmolecule(sketcherMinimizerAtom* excludedAtom)
{
    vector<sketcherMinimizerAtom*> subMolecule;
    queue<sketcherMinimizerAtom*> q;
    map<sketcherMinimizerAtom*, bool> isVisited;
    isVisited[excludedAtom] = true;
    q.push(this);
    isVisited[this] = true;
    while (!q.empty()) {
        sketcherMinimizerAtom* atom = q.front();
        subMolecule.push_back(atom);
        q.pop();
        for (sketcherMinimizerAtom* neighbor : atom->neighbors) {
            if (!isVisited[neighbor]) {
                q.push(neighbor);
                isVisited[neighbor] = true;
            }
        }
    }
    return subMolecule;
}
