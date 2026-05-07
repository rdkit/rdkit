/*
 Contributors: Nicola Zonta
 Copyright Schrodinger, LLC. All rights reserved
 */

#include <algorithm>
#include <map>
#include <numeric>

#include "CoordgenFragmentBuilder.h"
#include "CoordgenMinimizer.h"
#include "sketcherMinimizer.h"
#include "sketcherMinimizerFragment.h"
#include "sketcherMinimizerRing.h"
#include "sketcherMinimizerStretchInteraction.h"

using namespace std;

const int bondLength = BONDLENGTH;
const int PERFECTLY_PLANAR_SYSTEM_SCORE = 50;
const int NON_PLANAR_SYSTEM_SCORE = 1000;
const int UNTREATABLE_SYSTEM_PLANARITY_SCORE = 200000;

const int MACROCYCLE_CENTRAL_RING_SCORE = 1000;
const int NUMBER_OF_FUSED_RINGS_CENTRAL_RING_SCORE = 40;
const int NUMBER_OF_FUSION_ATOMS_CENTRAL_RING_SCORE = 15;
const int NEIGHBOR_ALREADY_BUILT_RING_SCORE = 100000;

void CoordgenFragmentBuilder::initializeCoordinates(
    sketcherMinimizerFragment* fragment) const
{
    assert(!fragment->getAtoms().empty());
    buildFragment(fragment);
    fragment->storeCoordinateInformation();
}

void CoordgenFragmentBuilder::rotateMainFragment(
    sketcherMinimizerFragment* f) const
{
    if (f->fixed) {
        return;
    }
    if (f->isTemplated) {
        return;
    }
    if (!f->constrained) {
        return;
    }

    sketcherMinimizerPointF constrainOldCenter(0.f, 0.f);
    sketcherMinimizerPointF constrainNewCenter(0.f, 0.f);
    vector<sketcherMinimizerAtom*> constrainedAtoms;
    for (const auto& atom : f->getAtoms()) {
        if (atom->constrained) {
            constrainedAtoms.push_back(atom);
        }
    }
    for (const auto& child : f->_children) {
        sketcherMinimizerAtom* atom = child->_bondToParent->endAtom;
        if (atom->constrained) {
            constrainedAtoms.push_back(atom);
        }
    }
    for (const auto& a : constrainedAtoms) {
        constrainOldCenter += a->templateCoordinates;
        constrainNewCenter += a->coordinates;
    }
    if (!constrainedAtoms.empty()) {
        constrainOldCenter /= constrainedAtoms.size();
        constrainNewCenter /= constrainedAtoms.size();
    }
    vector<sketcherMinimizerPointF> v1, v2;
    for (const auto& a : constrainedAtoms) {
        v2.push_back(a->coordinates - constrainNewCenter);
        v1.push_back(a->templateCoordinates - constrainOldCenter);
    }
    float rotMat[4];
    sketcherMinimizer::alignmentMatrix(v1, v2, rotMat);
    vector<sketcherMinimizerPointF> rotatedV2;
    for (const auto& p : v2) {
        auto rotatedPoint =
            sketcherMinimizerPointF(p.x() * rotMat[0] + p.y() * rotMat[1],
                                    p.x() * rotMat[2] + p.y() * rotMat[3]);

        rotatedV2.push_back(rotatedPoint);
    }
    for (sketcherMinimizerAtom* a : f->getAtoms()) {
        sketcherMinimizerPointF v = a->getCoordinates() - constrainNewCenter;
        v = sketcherMinimizerPointF(v.x() * rotMat[0] + v.y() * rotMat[1],
                                    v.x() * rotMat[2] + v.y() * rotMat[3]);
        a->setCoordinates(v + constrainOldCenter);
    }
    for (sketcherMinimizerFragment* child : f->_children) {
        sketcherMinimizerAtom* at = child->_bondToParent->endAtom;
        sketcherMinimizerPointF v = at->getCoordinates() - constrainNewCenter;
        v = sketcherMinimizerPointF(v.x() * rotMat[0] + v.y() * rotMat[1],
                                    v.x() * rotMat[2] + v.y() * rotMat[3]);
        at->setCoordinates(v + constrainOldCenter);
        at->coordinatesSet = false;
    }
}

bool CoordgenFragmentBuilder::findTemplate(
    const vector<sketcherMinimizerRing*>& rings) const
{
    vector<sketcherMinimizerAtom*> allAtoms =
        rings[0]->_atoms[0]->fragment->getAtoms();

    bool foundTemplate = false;
    sketcherMinimizer::loadTemplates();

    vector<sketcherMinimizerAtom*> atoms;
    vector<sketcherMinimizerBond*> bonds;
    map<sketcherMinimizerAtom*, bool> isVisited;

    vector<int> oldIndices;
    for (sketcherMinimizerRing* r : rings) {
        for (sketcherMinimizerAtom* a : r->_atoms) {
            if (!isVisited[a]) {
                isVisited[a] = true;
                oldIndices.push_back(a->_generalUseN);
                a->_generalUseN = static_cast<int>(atoms.size());

                atoms.push_back(a);
            }
        }
    }
    for (sketcherMinimizerRing* r : rings) {
        for (sketcherMinimizerBond* b : r->_bonds) {
            bool found = false;
            for (sketcherMinimizerBond* b2 : bonds) {
                if (b2 == b) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                bonds.push_back(b);
            }
        }
    }
    vector<unsigned int> mapping;
    for (auto& temp : sketcherMinimizer::m_templates.getTemplates()) {
        foundTemplate = sketcherMinimizer::compare(atoms, bonds, temp, mapping);
        if (foundTemplate) {
            if (!atoms.empty()) {
                atoms[0]->fragment->isTemplated = true;
            }
            for (unsigned int i = 0; i < atoms.size(); i++) {
                atoms[i]->setCoordinates(temp->_atoms[mapping[i]]->coordinates *
                                         bondLength);
                atoms[i]->rigid = true;
            }
            for (sketcherMinimizerRing* r : rings) {
                r->coordinatesGenerated = true;
            }
            break;
        }
    }
    for (unsigned int i = 0; i < oldIndices.size(); i++) {
        atoms[i]->_generalUseN = oldIndices[i];
    }
    return foundTemplate;
}

sketcherMinimizerRing* CoordgenFragmentBuilder::findCentralRingOfSystem(
    const vector<sketcherMinimizerRing*>& rings) const
{
    sketcherMinimizerRing* highest = nullptr;
    size_t high_score = 0;
    for (sketcherMinimizerRing* r : rings) {
        size_t priority = 0;
        /*keep growing the system building rings neighboring already built
         * rings*/
        for (auto neighborRing : r->fusedWith) {
            if (neighborRing->coordinatesGenerated) {
                priority += NEIGHBOR_ALREADY_BUILT_RING_SCORE;
                break;
            }
        }
        if (r->isMacrocycle()) {
            priority += MACROCYCLE_CENTRAL_RING_SCORE;
        }
        if (r->_atoms.size() == 6) {
            priority += 10;
        }
        priority += r->_atoms.size();
        priority +=
            NUMBER_OF_FUSED_RINGS_CENTRAL_RING_SCORE * (r->fusedWith.size());
        for (auto fusionAtoms : r->fusionAtoms) {
            priority +=
                NUMBER_OF_FUSION_ATOMS_CENTRAL_RING_SCORE * fusionAtoms.size();
        }
        if (priority > high_score || highest == nullptr) {
            highest = r;
            high_score = priority;
        }
    }
    if (highest == nullptr) {
        return rings.at(0);
    }
    return highest;
}

void CoordgenFragmentBuilder::generateCoordinatesCentralRings(
    vector<sketcherMinimizerRing*> rings) const
{
    if (rings.size() == 1) { // linearfusedRings
        buildRing(rings[0]);
    } else {
        bool foundTemplate = findTemplate(rings);
        if (!foundTemplate) {
            float planarityScore = newScorePlanarity(rings);
            if (planarityScore < NON_PLANAR_SYSTEM_SCORE) {
                bool needsTemplate =
                    planarityScore > PERFECTLY_PLANAR_SYSTEM_SCORE;
                if (needsTemplate) {
                    findTemplate(rings);
                }
                while (!rings.empty()) {
                    auto centralRing = findCentralRingOfSystem(rings);
                    buildRing(centralRing);
                    rings.erase(
                        std::remove(rings.begin(), rings.end(), centralRing),
                        rings.end());
                }
                CoordgenMinimizer::maybeMinimizeRings(rings);
            } else if (planarityScore > UNTREATABLE_SYSTEM_PLANARITY_SCORE) {
                return;
            } else {
                sketcherMinimizerRing* firstRing =
                    findCentralRingOfSystem(rings);
                m_macrocycleBuilder.openCycleAndGenerateCoords(firstRing);
                firstRing->getAtoms()
                    .at(0)
                    ->getMolecule()
                    ->requireMinimization();
            }
        }
    }
}

float CoordgenFragmentBuilder::newScorePlanarity(
    const vector<sketcherMinimizerRing*>& rings)
    const // if score > 1000 then it is not planar
{
    float score = 0.f;
    for (const auto& ring : rings) {
        if (ring->isMacrocycle() &&
            m_macrocycleBuilder.findBondToOpen(ring) == nullptr) {
            continue;
        }
        if (ring->isMacrocycle()) {
            for (const auto& otherRing : ring->fusedWith) {
                if (otherRing->isMacrocycle()) {
                    score += NON_PLANAR_SYSTEM_SCORE;
                }
            }
        }
        for (const auto& bond : ring->_bonds) {
            if (bond->rings.size() > 2) {
                score += NON_PLANAR_SYSTEM_SCORE * (bond->rings.size() - 2);
            }
        }
        for (const auto& atom : ring->getAtoms()) {
            if (atom->neighbors.size() > 3) {
                float angle = 0;
                for (sketcherMinimizerRing* r : atom->rings) {
                    angle += static_cast<float>(M_PI -
                                                (2 * M_PI / r->_atoms.size()));
                }
                if (angle >= 1.99 * M_PI) {
                    score += NON_PLANAR_SYSTEM_SCORE;
                }
            }
        }
    }
    return score;
}

void CoordgenFragmentBuilder::generateCoordinatesSideRings(
    stack<sketcherMinimizerRing*> sideRings) const
{
    while (!sideRings.empty()) {
        sketcherMinimizerRing* ring = sideRings.top();
        sideRings.pop();
        buildRing(ring);
    }
}

sketcherMinimizerRing*
CoordgenFragmentBuilder::getSharedAtomsWithAlreadyDrawnRing(
    const sketcherMinimizerRing* ring,
    vector<sketcherMinimizerAtom*>& fusionAtoms,
    sketcherMinimizerBond*& fusionBond) const
{
    sketcherMinimizerRing* parent = nullptr;
    for (auto i : ring->fusedWith) {
        if (i->coordinatesGenerated) {
            if (parent != nullptr) {
                if (i->getFusionAtomsWith(ring).size() <
                        parent->getFusionAtomsWith(ring).size() ||
                    (i->size() < parent->size()))
                    continue;
            }
            parent = i;
        }
    }
    if (parent) {
        for (unsigned int i = 0; i < parent->fusedWith.size(); i++) {
            if (parent->fusedWith[i] == ring) {
                fusionAtoms = parent->fusionAtoms[i];
            }
        }
        for (sketcherMinimizerBond* b : parent->fusionBonds) {
            if (ring->containsAtom(b->startAtom) ||
                ring->containsAtom(b->endAtom)) {
                fusionBond = b;
                break;
            }
        }
    }

    return parent;
}

vector<sketcherMinimizerAtom*> CoordgenFragmentBuilder::orderChainOfAtoms(
    const vector<sketcherMinimizerAtom*>& atoms,
    sketcherMinimizerAtom* startAtom)
{
    vector<sketcherMinimizerAtom*> orderedAtoms;
    map<sketcherMinimizerAtom*, bool> stillToAdd;
    for (sketcherMinimizerAtom* a : atoms) {
        stillToAdd[a] = true;
    }
    sketcherMinimizerAtom* atomToAdd = startAtom;
    while (orderedAtoms.size() < atoms.size()) {
        orderedAtoms.push_back(atomToAdd);
        stillToAdd[atomToAdd] = false;
        if (orderedAtoms.size() >= atoms.size()) {
            break;
        }

        for (sketcherMinimizerAtom* neighbor : atomToAdd->neighbors) {
            if (stillToAdd[neighbor]) {
                atomToAdd = neighbor;
                break;
            }
        }
    }
    return orderedAtoms;
}

vector<sketcherMinimizerAtom*>
CoordgenFragmentBuilder::orderRingAtoms(const sketcherMinimizerRing* ring)
{
    vector<sketcherMinimizerAtom*> ringAtoms = ring->getAtoms();
    assert(!ringAtoms.empty());
    return orderChainOfAtoms(ringAtoms, ringAtoms.at(0));
}

void CoordgenFragmentBuilder::buildRing(sketcherMinimizerRing* ring) const
{
    if (ring->coordinatesGenerated) {
        return;
    }

    vector<sketcherMinimizerAtom*> atoms;

    sketcherMinimizerRing* parent = nullptr;
    sketcherMinimizerBond* fusionBond = nullptr;
    vector<sketcherMinimizerAtom*> fusionAtoms;
    parent = getSharedAtomsWithAlreadyDrawnRing(ring, fusionAtoms, fusionBond);
    atoms = orderRingAtoms(ring);

    vector<sketcherMinimizerPointF> coords;
    if (atoms.size() >= MACROCYCLE) {
        coords = m_macrocycleBuilder.newMacrocycle(ring, atoms);
    } else {
        coords = listOfCoordinatesFromListofRingAtoms(atoms);
    }
    if (parent) {
        if (fusionAtoms.size() < 2) {
            sketcherMinimizerPointF coordinatesOfPivotAtom;

            sketcherMinimizerAtom* pivotAtom = nullptr;
            sketcherMinimizerAtom* pivotAtomOnParent = nullptr;
            if (!fusionAtoms.empty()) {
                pivotAtom = fusionAtoms.at(0);
                pivotAtomOnParent = fusionAtoms.at(0);
            } else {
                assert(fusionBond != nullptr);
                if (find(atoms.begin(), atoms.end(),
                         fusionBond->getStartAtom()) != atoms.end()) {
                    pivotAtom = fusionBond->getStartAtom();
                    pivotAtomOnParent = fusionBond->getEndAtom();
                } else {
                    pivotAtom = fusionBond->getEndAtom();
                    pivotAtomOnParent = fusionBond->getStartAtom();
                }
            }
            assert(find(atoms.begin(), atoms.end(), pivotAtom) != atoms.end());
            assert(parent->containsAtom(pivotAtomOnParent));

            for (unsigned int i = 0; i < atoms.size(); i++) {
                if (atoms[i] == pivotAtom) {
                    coordinatesOfPivotAtom = coords[i];
                    break;
                }
            }

            for (auto& coord : coords) {
                coord -= coordinatesOfPivotAtom;
            }
            sketcherMinimizerPointF center =
                accumulate(coords.begin(), coords.end(),
                           sketcherMinimizerPointF(0, 0)) /
                coords.size();
            center.normalize();

            sketcherMinimizerPointF neighborsMean(0.f, 0.f);
            for (sketcherMinimizerAtom* neighbor :
                 pivotAtomOnParent->neighbors) {
                if (parent->containsAtom(neighbor)) {
                    neighborsMean += neighbor->getCoordinates();
                }
            }
            neighborsMean *= 0.5;
            sketcherMinimizerPointF pivotAtomCoordinates =
                pivotAtomOnParent->getCoordinates();
            sketcherMinimizerPointF direction =
                pivotAtomCoordinates - neighborsMean;
            direction.normalize();
            sketcherMinimizerPointF bondVector(0.f, 0.f);
            if (fusionBond != nullptr) {
                bondVector = direction * bondLength;
            }
            for (auto& coord : coords) {
                sketcherMinimizerPointF coordinates = coord;
                coordinates.rotate(center.y(), center.x());
                coordinates.rotate(-direction.y(), direction.x());
                coord = coordinates + pivotAtomCoordinates + bondVector;
            }
            for (unsigned int i = 0; i < coords.size(); i++) {
                atoms[i]->setCoordinates(coords[i]);
            }
            if (fusionBond) {
                // check ZE inversions
                if (!fusionBond->checkStereoChemistry()) {
                    for (sketcherMinimizerAtom* atom : atoms) {
                        sketcherMinimizerAtom::mirrorCoordinates(atom,
                                                                 fusionBond);
                    }
                }
            }

        } else {
            // build 2 lists of coordinates, mirror of each other.

            vector<sketcherMinimizerPointF> coords2 = coords;
            for (auto& i : coords2) {
                i = sketcherMinimizerPointF(i.x(), -i.y());
            }
            map<sketcherMinimizerAtom*, sketcherMinimizerPointF> map1, map2;
            for (unsigned int i = 0; i < atoms.size(); i++) {
                map1[atoms[i]] = coords[i];
                map2[atoms[i]] = coords2[i];
            }
            assert(!fusionAtoms.empty());
            size_t lastI = fusionAtoms.size() - 1;
            sketcherMinimizerPointF mean = (fusionAtoms[0]->coordinates +
                                            fusionAtoms[lastI]->coordinates) *
                                           0.5;
            sketcherMinimizerPointF newMean =
                (map1[fusionAtoms[0]] + map1[fusionAtoms[lastI]]) * 0.5;
            sketcherMinimizerPointF v =
                fusionAtoms[lastI]->coordinates - fusionAtoms[0]->coordinates;
            float distance = v.length();
            if (distance < SKETCHER_EPSILON) {
                distance = SKETCHER_EPSILON;
            }
            float sine = v.y() / distance;
            float cosine = v.x() / distance;
            sketcherMinimizerPointF newV =
                (map1[fusionAtoms[lastI]] - map1[fusionAtoms[0]]);
            float newDistance = newV.length();
            float newSine = newV.y() / newDistance;
            float newCosine = newV.x() / newDistance;

            for (unsigned int i = 0; i < atoms.size();
                 i++) { // rotate and translate coords
                sketcherMinimizerPointF d = coords[i] - newMean;
                d.rotate(newSine, newCosine);
                d.rotate(-sine, cosine);
                coords[i] = mean + d;
            }
            newV = (map2[fusionAtoms[lastI]] - map2[fusionAtoms[0]]);
            newMean = (map2[fusionAtoms[0]] + map2[fusionAtoms[lastI]]) * 0.5;

            newDistance = newV.length();
            newSine = newV.y() / newDistance;
            newCosine = newV.x() / newDistance;

            for (unsigned int i = 0; i < atoms.size();
                 i++) { // rotate and translate coords2
                sketcherMinimizerPointF d = coords2[i] - newMean;
                d.rotate(newSine, newCosine);
                d.rotate(-sine, cosine);
                coords2[i] = mean + d;
            }
            // pick the list of coordinates that will place the ring outside the
            // parent i.e further away from its center.
            sketcherMinimizerPointF center = parent->findCenter();
            map<sketcherMinimizerAtom*, bool> fusionMap;
            for (auto fusionAtom : fusionAtoms) {
                fusionMap[fusionAtom] = true;
            }
            float distance1 = 0.f, distance2 = 0.f;
            for (unsigned int i = 0; i < atoms.size(); i++) {
                if (fusionMap[atoms[i]]) {
                    continue;
                }
                sketcherMinimizerPointF v1 = coords[i] - center;
                sketcherMinimizerPointF v2 = coords2[i] - center;
                distance1 += v1.length();
                distance2 += v2.length();
            }
            vector<sketcherMinimizerPointF>* targetCoords = &coords;
            if (distance1 < distance2) {
                targetCoords = &coords2;
            }
            for (unsigned int i = 0; i < atoms.size(); i++) {
                atoms[i]->setCoordinates((*targetCoords)[i]);
            }
        }
    } else {
        for (unsigned int i = 0; i < coords.size(); i++) {
            atoms[i]->setCoordinates(coords[i]);
        }
    }
    ring->coordinatesGenerated = true;
}

vector<sketcherMinimizerPointF>
CoordgenFragmentBuilder::listOfCoordinatesFromListofRingAtoms(
    const vector<sketcherMinimizerAtom*>& atoms)
{
    vector<sketcherMinimizerPointF> out;
    assert(!atoms.empty());
    auto a = static_cast<float>(2 * M_PI / atoms.size());
    sketcherMinimizerPointF coords(0.f, 0.f);
    float angle = 0;
    for (unsigned int n = 0; n < atoms.size(); n++) {
        out.push_back(coords);
        angle = a * n;
        coords += sketcherMinimizerPointF(cos(angle) * bondLength,
                                          -sin(angle) * bondLength);
    }
    return out;
}

void CoordgenFragmentBuilder::buildNonRingAtoms(
    sketcherMinimizerFragment* fragment) const
{
    set<sketcherMinimizerAtom*> isAtomVisited;
    queue<sketcherMinimizerAtom*> q;
    for (sketcherMinimizerAtom* atom : fragment->getAtoms()) {
        if (!atom->getRings().empty()) {
            q.push(atom);
            isAtomVisited.insert(atom);
        }
    }
    if (q.empty()) {
        if (fragment->getParent()) {
            fragment->_bondToParent->getStartAtom()->setCoordinates(
                sketcherMinimizerPointF(-bondLength, 0));
        }
        sketcherMinimizerAtom* atomToStartFrom =
            (fragment->getParent() ? fragment->_bondToParent->getEndAtom()
                                   : fragment->getAtoms()[0]);
        atomToStartFrom->setCoordinates(sketcherMinimizerPointF(0, 0));
        q.push(atomToStartFrom);
        isAtomVisited.insert(atomToStartFrom);
    }

    while (!q.empty()) {
        generateCoordinatesNeighborsOfFirstAtomInQueue(q, isAtomVisited,
                                                       fragment);
    }
}

vector<float> CoordgenFragmentBuilder::neighborsAnglesAtCenter(
    const sketcherMinimizerAtom* atom) const
{

    size_t angleDivision = atom->neighbors.size();
    vector<float> angles;
    if (!m_evenAngles) {
        if (atom->neighbors.size() == 2) {
            if (atom->getAtomicNumber() == 6 ||
                !atom->neighbors.at(0)->crossLayout ||
                !atom->neighbors.at(1)->crossLayout) {
                angleDivision = 3;
            }
            if (atom->bonds[0]->getBondOrder() +
                    atom->bonds[1]->getBondOrder() >=
                4) {
                angleDivision = 2;
            }

        } else if (atom->neighbors.size() == 4 && !atom->crossLayout) {
            angles.push_back(static_cast<float>(M_PI / 3));
            angles.push_back(static_cast<float>(M_PI * 0.5));
            angles.push_back(static_cast<float>(2 * M_PI / 3));
            angles.push_back(static_cast<float>(M_PI * 0.5));
        }
    }
    if (angles.empty()) {
        for (unsigned int i = 0; i < atom->neighbors.size(); i++) {
            angles.push_back(static_cast<float>(2 * M_PI / angleDivision));
        }
    }
    return angles;
}

void CoordgenFragmentBuilder::initializeVariablesForNeighboursCoordinates(
    sketcherMinimizerAtom* atom, set<sketcherMinimizerAtom*>& isAtomVisited,
    sketcherMinimizerPointF& startCoordinates,
    vector<sketcherMinimizerAtom*>& orderedNeighbors,
    vector<float>& angles) const
{
    if (!atom->rings.empty()) {
        return initializeVariablesForNeighboursCoordinatesRingAtom(
            atom, isAtomVisited, startCoordinates, orderedNeighbors, angles);
    }
    orderedNeighbors.clear();
    if (atom->neighbors.size() != 4) {
        orderedNeighbors = atom->neighbors;
    } else {
        vector<sketcherMinimizerAtomPriority> atomPriorities;
        for (auto& neighbor : atom->neighbors) {
            sketcherMinimizerAtomPriority p;
            p.a = neighbor;
            atomPriorities.push_back(p);
        }
        sketcherMinimizerAtom::orderAtomPriorities(atomPriorities, atom);
        for (sketcherMinimizerAtomPriority ap : atomPriorities) {
            orderedNeighbors.push_back(ap.a);
        }
    }
    unsigned int startN = 0;
    for (unsigned int i = 0; i < orderedNeighbors.size(); i++) {
        if (isAtomVisited.find(orderedNeighbors[i]) != isAtomVisited.end()) {
            startCoordinates =
                orderedNeighbors[i]->getCoordinates() - atom->getCoordinates();
            startN = i;
            break;
        }
    }
    for (unsigned int i = 0; i < startN; i++) {
        orderedNeighbors.push_back(orderedNeighbors.at(0));
        orderedNeighbors.erase(orderedNeighbors.begin());
    }
    angles = neighborsAnglesAtCenter(atom);
}

void CoordgenFragmentBuilder::
    initializeVariablesForNeighboursCoordinatesRingAtom(
        const sketcherMinimizerAtom* atom, set<sketcherMinimizerAtom*>&,
        sketcherMinimizerPointF& startCoordinates,
        vector<sketcherMinimizerAtom*>& orderedNeighbors,
        vector<float>& angles) const
{

    vector<pair<float, sketcherMinimizerAtom*>> ringNeighboursAndAngles;
    orderedNeighbors.clear();
    for (sketcherMinimizerAtom* neigh : atom->neighbors) {
        if (sketcherMinimizer::sameRing(neigh, atom)) {
            float ang = atan2(neigh->coordinates.y() - atom->coordinates.y(),
                              neigh->coordinates.x() - atom->coordinates.x());
            if (ang < 0) {
                ang += static_cast<float>(2 * M_PI);
            }
            pair<float, sketcherMinimizerAtom*> pairToAdd(ang, neigh);
            ringNeighboursAndAngles.push_back(pairToAdd);
        } else {
            orderedNeighbors.push_back(neigh);
        }
    }
    stable_sort(ringNeighboursAndAngles.begin(), ringNeighboursAndAngles.end());
    vector<float> gaps;
    vector<float> scaledGaps;
    for (unsigned int i = 0; i < ringNeighboursAndAngles.size(); i++) {
        int next = (i + 1) % ringNeighboursAndAngles.size();
        float gap = ringNeighboursAndAngles[next].first -
                    ringNeighboursAndAngles[i].first;
        if (gap < 0) {
            gap += static_cast<float>(2 * M_PI);
        }
        gaps.push_back(gap);
        bool rin = false;
        float mid = ringNeighboursAndAngles[i].first + gap * 0.5f;
        float sine = sin(-mid);
        float cosine = cos(-mid);
        sketcherMinimizerPointF p(bondLength * 0.1f, 0);
        p.rotate(sine, cosine);
        sketcherMinimizerPointF point = atom->coordinates + p;
        vector<sketcherMinimizerRing*> rings = atom->getFragment()->getRings();
        for (sketcherMinimizerRing* r : rings) {
            if (r->isMacrocycle()) {
                continue;
            }
            if (r->contains(point)) {
                rin = true;
                break;
            }
        }
        float scaledGap = gap;
        if (gap > M_PI) {
            scaledGap *= 10.f;
        } else if (rin) {
            scaledGap *= 0.2f;
        }
        scaledGaps.push_back(scaledGap);
    }
    int bestI = 0;
    for (unsigned int i = 0; i < gaps.size(); i++) {
        if (scaledGaps[i] > scaledGaps[bestI]) {
            bestI = i;
        }
    }
    angles.clear();
    float biggestGap = gaps[bestI];
    startCoordinates = ringNeighboursAndAngles[bestI].second->getCoordinates() -
                       atom->getCoordinates();
    for (unsigned int i = 0; i < orderedNeighbors.size(); i++) {
        angles.push_back(-biggestGap / (orderedNeighbors.size() + 1));
    }
}

void CoordgenFragmentBuilder::generateCoordinatesNeighborsOfFirstAtomInQueue(
    queue<sketcherMinimizerAtom*>& atomQueue,
    set<sketcherMinimizerAtom*>& isAtomVisited,
    const sketcherMinimizerFragment* fragment) const
{

    sketcherMinimizerAtom* at = atomQueue.front();
    atomQueue.pop();
    sketcherMinimizerPointF initialCoords(-bondLength, 0.f);
    vector<sketcherMinimizerAtom*> orderedNeighbors;
    vector<float> angles;

    initializeVariablesForNeighboursCoordinates(
        at, isAtomVisited, initialCoords, orderedNeighbors, angles);

    for (unsigned int i = 0; i < orderedNeighbors.size(); i++) {
        int n = i;
        sketcherMinimizerAtom* neigh = orderedNeighbors[n];
        if (isAtomVisited.find(neigh) != isAtomVisited.end()) {
            continue;
        }
        float s = sin(angles[n]);
        float c = cos(angles[n]);
        initialCoords.rotate(s, c);
        neigh->setCoordinates(at->coordinates + initialCoords);
        if (neigh->fragment == fragment) {
            isAtomVisited.insert(neigh);
            atomQueue.push(neigh);
        } else {
            neigh->coordinatesSet = false;
        }
        if (at->needsCheckForClashes) {
            neigh->needsCheckForClashes = true;
            CoordgenMinimizer::checkForClashes(neigh);
        }
        for (auto& dof : at->fragment->getDofsOfAtom(at)) {
            if (dof->getFragment() == fragment) {
                dof->addAtom(neigh);
            }
        }
    }
    avoidZEInversions(at, isAtomVisited);
    maybeAddMacrocycleDOF(at);
    for (auto& neighbor : at->neighbors) {
        if (!sketcherMinimizerAtom::shareARing(at, neighbor) &&
            at->getFragment() == neighbor->getFragment()) {
            auto* dof = new CoordgenScaleAtomsDOF(at);
            dof->addAtom(neighbor);
            at->fragment->addDof(dof);
        }
    }
}

void CoordgenFragmentBuilder::maybeAddMacrocycleDOF(
    sketcherMinimizerAtom* atom) const
{
    if (atom->getRings().size() == 1 &&
        atom->getRings().at(0)->isMacrocycle() && atom->neighbors.size() == 3) {
        for (auto& bond : atom->getBonds()) {
            if (bond->isStereo() && !bond->isTerminal()) {
                return;
            }
        }
        for (auto& neighbor : atom->neighbors) {
            if (!sketcherMinimizerAtom::shareARing(atom, neighbor)) {
                auto* dof = new CoordgenInvertBondDOF(atom, neighbor);
                atom->fragment->addDof(dof);
            }
        }
    }
}

void CoordgenFragmentBuilder::avoidZEInversions(
    const sketcherMinimizerAtom* at,
    set<sketcherMinimizerAtom*>& isAtomVisited) const
{
    if (!at->getRings().empty()) {
        return;
    }
    sketcherMinimizerBond* doubleBond = nullptr;
    vector<sketcherMinimizerAtom*> atomsToMirror;
    for (unsigned int i = 0; i < at->bonds.size(); i++) {
        if (at->bonds[i]->isStereo() &&
            isAtomVisited.find(at->neighbors[i]) != isAtomVisited.end()) {
            doubleBond = at->bonds[i];
        } else {
            atomsToMirror.push_back(at->neighbors[i]);
        }
    }
    if (doubleBond == nullptr) {
        return;
    }
    if (!atomsToMirror.empty() && doubleBond) {
        sketcherMinimizerAtom* firstCIPNeighborStart =
            doubleBond->startAtomCIPFirstNeighbor();
        if (firstCIPNeighborStart == nullptr) {
            return;
        }
        sketcherMinimizerAtom* firstCIPNeighborEnd =
            doubleBond->endAtomCIPFirstNeighbor();
        if (firstCIPNeighborEnd == nullptr) {
            return;
        }
        if (!doubleBond->checkStereoChemistry()) {
            for (sketcherMinimizerAtom* a : atomsToMirror) {
                sketcherMinimizerAtom::mirrorCoordinates(a, doubleBond);
            }
        }
    }
}

void CoordgenFragmentBuilder::buildFragment(
    sketcherMinimizerFragment* fragment) const
{
    buildRings(fragment);
    buildNonRingAtoms(fragment);
    CoordgenMinimizer::avoidInternalClashes(fragment);
    fallbackIfNanCoordinates(fragment);
    if (!fragment->getParent() && fragment->constrained) {
        rotateMainFragment(fragment);
    }

    // reset coordinates for fixed fragments:
    if (fragment->fixed) {
        fragment->setAllCoordinatesToTemplate();
    }
}

void CoordgenFragmentBuilder::fallbackIfNanCoordinates(
    sketcherMinimizerFragment* fragment) const
{
    vector<sketcherMinimizerAtom*> fragmentAtoms = fragment->getAtoms();
    if (CoordgenMinimizer::hasNaNCoordinates(fragmentAtoms) &&
        CoordgenMinimizer::hasValid3DCoordinates(fragmentAtoms)) {
        CoordgenMinimizer::fallbackOn3DCoordinates(fragmentAtoms);
    }
}

void CoordgenFragmentBuilder::initializeFusedRingInformation(
    sketcherMinimizerFragment* fragment) const
{
    /*
     for each ring in the fragment find out which other rings it shares atoms
     with,
     or if they are connected in other ways (i.e. via a double bond).
     */
    if (fragment->getRings().size() < 2) {
        return;
    }
    vector<sketcherMinimizerAtom*> fragmentAtoms = fragment->getAtoms();
    for (sketcherMinimizerAtom* atom : fragmentAtoms) {
        if (atom->rings.size() > 1) {
            for (unsigned int i = 0; i < atom->rings.size(); i++) {
                for (unsigned int j = i + 1; j < atom->rings.size(); j++) {

                    bool found = false;
                    for (unsigned int r = 0;
                         r < atom->rings[i]->fusedWith.size(); r++) {
                        if (atom->rings[i]->fusedWith[r] == atom->rings[j]) {
                            atom->rings[i]->fusionAtoms[r].push_back(atom);
                            for (unsigned int rr = 0;
                                 rr < atom->rings[j]->fusedWith.size(); rr++) {

                                if (atom->rings[j]->fusedWith[rr] ==
                                    atom->rings[i]) {
                                    atom->rings[j]->fusionAtoms[rr].push_back(
                                        atom);
                                    break;
                                }
                            }
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        atom->rings[i]->fusedWith.push_back(atom->rings[j]);
                        atom->rings[j]->fusedWith.push_back(atom->rings[i]);
                        vector<sketcherMinimizerAtom*> ats;
                        ats.push_back(atom);
                        atom->rings[i]->fusionAtoms.push_back(ats);
                        atom->rings[j]->fusionAtoms.push_back(ats);
                    }
                }
            }
        }
    }
    vector<sketcherMinimizerBond*> bonds = fragment->getBonds();
    for (sketcherMinimizerBond* b : bonds) {
        if ((b->bondOrder != 1) && (b->startAtom->rings.size()) &&
            (b->endAtom->rings.size()) &&
            (!sketcherMinimizer::sameRing(b->startAtom, b->endAtom))) {
            for (sketcherMinimizerRing* r : b->startAtom->rings) {
                for (sketcherMinimizerRing* r2 : b->endAtom->rings) {
                    r->fusedWith.push_back(r2);
                    r2->fusedWith.push_back(r);
                    vector<sketcherMinimizerAtom*> ats;
                    r->fusionAtoms.push_back(ats);
                    r2->fusionAtoms.push_back(ats);
                    r->fusionBonds.push_back(b);
                    r2->fusionBonds.push_back(b);
                }
            }
        }
    }
    vector<sketcherMinimizerRing*> fragmentRings = fragment->getRings();
    for (sketcherMinimizerRing* r : fragmentRings) {
        for (unsigned int i = 0; i < r->fusedWith.size(); i++) {
            vector<sketcherMinimizerAtom*> fusionAtoms = r->fusionAtoms[i];
            if (fusionAtoms.size() <= 2) {
                continue;
            }
            // find an atom with only one neighbor in the vector
            sketcherMinimizerAtom* startAtom = nullptr;
            for (sketcherMinimizerAtom* a : fusionAtoms) {
                int counter = 0;
                for (sketcherMinimizerAtom* n : a->neighbors) {
                    if (find(fusionAtoms.begin(), fusionAtoms.end(), n) !=
                        fusionAtoms.end()) {
                        counter++;
                    }
                }
                if (counter == 1) {
                    startAtom = a;
                    break;
                }
            }
            assert(startAtom != nullptr);
            r->fusionAtoms[i] = orderChainOfAtoms(fusionAtoms, startAtom);
        }
    }

    for (sketcherMinimizerRing* r : fragmentRings) {
        for (unsigned int ii = 0; ii < r->fusedWith.size(); ii++) {
            if (r->fusionAtoms[ii].size() > 2) {
                for (sketcherMinimizerAtom* a :
                     r->fusionAtoms[ii]) { // looking for atoms with no
                                           // neighbors in r but not in
                                           // r->fusedWith[ii]
                    bool found = false;
                    for (sketcherMinimizerAtom* n : a->neighbors) {
                        int nn = 0;
                        for (sketcherMinimizerRing* nr : n->rings) {
                            if (nr == r) {
                                nn++;
                            } else if (nr == r->fusedWith[ii]) {
                                nn++;
                            }
                        }
                        if (nn == 1) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        a->isSharedAndInner = true;
                    }
                }
            }
        }
    }
}

void CoordgenFragmentBuilder::simplifyRingSystem(
    const vector<sketcherMinimizerRing*>& allRings,
    stack<sketcherMinimizerRing*>& sideRings,
    vector<sketcherMinimizerRing*>& centralRings) const
{
    map<sketcherMinimizerRing*, bool> alreadyChosenAsSide;

    bool found = true;
    /* simplify the structure by iteratively removing rings that have only one
     fusion partner. These side rings will have coordinates generated
     separately,
     one at a time, starting from the closest to the core.
     If we think of rings as atoms and fusion between two rings as bonds, this
     is the same as stripping away non-ring atoms from a fragment to generate
     coordinate
     for the ring system first
     */
    while (found) {
        found = false;
        for (sketcherMinimizerRing* r : allRings) {
            if (alreadyChosenAsSide[r]) {
                continue;
            }
            if (r->isMacrocycle()) {
                continue;
            }
            int n = 0;
            for (auto ringCounter = 0u; ringCounter < r->fusedWith.size();
                 ++ringCounter) {
                auto fusedRing = r->fusedWith.at(ringCounter);
                if (!alreadyChosenAsSide[fusedRing]) {
                    n++;
                    if (fusedRing->isMacrocycle()) {
                        /* disable macrocycles */
                        n++;
                    }
                    if (r->fusionAtoms.at(ringCounter).size() > 3) {
                        /* disable rings that share too many atoms */
                        n++;
                    }
                    if (r->fusionAtoms.at(ringCounter).size() == 3 &&
                        r->size() == 4 && fusedRing->size() == 4) {
                        /* don't separate rings of bicyclo (1,1,1) pentane so we
                         * can use a template instead */
                        n++;
                    }
                }
            }
            if (n == 1) /* ring is fused with only one ring,
                         it will be removed from the core of the system
                         and treated as a side ring
                         */
            {
                sideRings.push(r);
                found = true;
                alreadyChosenAsSide[r] = true;
                r->side = true;
            }
        }
    }
    /*
     everything that is left after removing the side rings is the core of the
     ring system
     */
    for (sketcherMinimizerRing* ring : allRings) {
        if (!alreadyChosenAsSide[ring]) {
            centralRings.push_back(ring);
        }
    }
}

void CoordgenFragmentBuilder::buildRings(
    sketcherMinimizerFragment* fragment) const
{
    if (fragment->getRings().empty()) {
        return;
    }
    initializeFusedRingInformation(fragment);
    vector<sketcherMinimizerRing*> fragmentRings = fragment->getRings();
    if (fragmentRings.size() == 1) {
        buildRing(fragmentRings[0]);
    } else { // more than 1 ring
        stack<sketcherMinimizerRing*> sideRings;
        vector<sketcherMinimizerRing*> centralRings;
        simplifyRingSystem(fragmentRings, sideRings, centralRings);
        generateCoordinatesCentralRings(centralRings);
        generateCoordinatesSideRings(sideRings);
    }
}
