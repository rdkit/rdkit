/*
 Contributors: Nicola Zonta
 Copyright Schrodinger, LLC. All rights reserved
 */

#include "CoordgenMinimizer.h"
#include "sketcherMinimizer.h" // should be removed at the end of refactoring
#include "sketcherMinimizerAtom.h"
#include "sketcherMinimizerBendInteraction.h"
#include "sketcherMinimizerBond.h"
#include "sketcherMinimizerClashInteraction.h"
#include "sketcherMinimizerConstraintInteraction.h"
#include "sketcherMinimizerEZConstrainInteraction.h"
#include "sketcherMinimizerFragment.h"
#include "sketcherMinimizerMaths.h"
#include "sketcherMinimizerResidue.h"
#include "sketcherMinimizerResidueInteraction.h"
#include "sketcherMinimizerRing.h"
#include "sketcherMinimizerStretchInteraction.h"
#include <algorithm>
#include <limits>
#include <queue>

using namespace std;
static const float bondLength = BONDLENGTH;
static const float clashEnergyThreshold = 10;

#define SAME_SIDE_DPR_PENALTY 100
#define SAME_SIDE_DPR_PENALTY_2 50
static const float FORCE_MULTIPLIER = 0.3f;

static const float STANDARD_CROSSING_BOND_PENALTY = 2500.f;
static const float TERMINAL_BOND_CROSSING_MULTIPLIER = 0.5f;
static const float MACROCYCLE_BOND_CROSSING_MULTIPLIER = 8.f;
static const float RING_BOND_CROSSING_MULTIPLIER = 2.f;

static const unsigned int MAXIMUM_NUMBER_OF_SCORED_SOLUTIONS = 10000;
static const float REJECTED_SOLUTION_SCORE = 99999999.f;

static const unsigned int ITERATION_HISTORY_SIZE = 100;
static const float MAX_NET_ENERGY_CHANGE = 20.f;

CoordgenMinimizer::CoordgenMinimizer()
{
    m_maxIterations = 1000;
    skipMinimization = false;
    skipFlipFragments = false;
    skipAvoidClashes = false;
    m_scoreResidueInteractions = true;
    m_precision = 1.f;
    energy_list = {};
    all_coordinates = {};
}

CoordgenMinimizer::~CoordgenMinimizer()
{
    clearInteractions();
}

void CoordgenMinimizer::clearInteractions()
{
    for (auto& _interaction : _interactions) {
        delete _interaction;
    }
    _interactions.clear();
    _intramolecularClashInteractions.clear();
    _extraInteractions.clear();
    _stretchInteractions.clear();
    _bendInteractions.clear();
}

void CoordgenMinimizer::run()
{
    if (skipMinimization) {
        return;
    }
    if (_interactions.empty()) {
        setupInteractions();
    }

#ifdef DEBUG_MINIMIZATION_COORDINATES
    // to seperate energy and DOF minimization
    energy_list.push_back(-1.f);
    all_coordinates.push_back({sketcherMinimizerPointF(0,0)});
#endif

    std::vector<float> local_energy_list(m_maxIterations);
    std::vector<sketcherMinimizerPointF> lowest_energy_coords(m_atoms.size());
    float min_energy = std::numeric_limits<float>::max();
    for (unsigned int iterations = 0; iterations < m_maxIterations; ++iterations) {
        local_energy_list[iterations] = scoreInteractions();
        // track coordinates with lowest energy
        if (local_energy_list[iterations] < min_energy) {
            for (size_t i = 0; i < m_atoms.size(); ++i) {
                lowest_energy_coords[i] = m_atoms[i]->coordinates;
            }
        }
#ifdef DEBUG_MINIMIZATION_COORDINATES
        // store data from this minimization step to be written to a file later
        energy_list.push_back(local_energy_list[iterations]);
        std::vector<sketcherMinimizerPointF> these_coordinates;
        for (auto atom : m_atoms) {
            these_coordinates.push_back(atom->coordinates);
        }
        all_coordinates.push_back(these_coordinates);
#endif
        if (!applyForces(0.1f)) {
            break;
        }
        if (iterations < 2 * ITERATION_HISTORY_SIZE) {
            continue;
        }
        if (local_energy_list[iterations - ITERATION_HISTORY_SIZE] - local_energy_list[iterations] < MAX_NET_ENERGY_CHANGE) {
            break;
        }
    }
    // set coordinates back to lowest energy state
    if (min_energy < std::numeric_limits<float>::max()) {
        for (size_t i = 0; i < m_atoms.size(); ++i) {
            m_atoms[i]->coordinates = lowest_energy_coords[i];
        }
    }
}

bool CoordgenMinimizer::applyForces(float maxd)
{
    float delta = 0.001f; // minimum squared displacement
    float distance = 0.f;
    for (auto atom : m_atoms) {
        if (atom->fixed) {
            continue;
        }
        sketcherMinimizerPointF displacement = atom->force * FORCE_MULTIPLIER;
        if (displacement.x() != displacement.x() ||
            displacement.y() != displacement.y()) {
            displacement = sketcherMinimizerPointF(0.f, 0.f);
        }
        float dsquare = displacement.x() * displacement.x() +
                        displacement.y() * displacement.y();
        if (dsquare < SKETCHER_EPSILON) {
            dsquare = SKETCHER_EPSILON;
        }
        if (dsquare > maxd * maxd) {
            displacement *= maxd / sqrt(dsquare);
        }
        atom->coordinates += displacement;
        distance += displacement.squareLength();
        atom->force = sketcherMinimizerPointF(0, 0);
    }
    return distance >= delta;
}

/* store extra interaction to be used when minimizing molecule.
 cis amides constraints are an example as they need 3d coordinates
 to be detected */
void CoordgenMinimizer::addExtraInteraction(
    sketcherMinimizerMolecule* molecule,
    sketcherMinimizerInteraction* interaction)
{
    _extraInteractionsOfMolecule[molecule].push_back(interaction);
}

void CoordgenMinimizer::addClashInteractionsOfMolecule(
    sketcherMinimizerMolecule* molecule, bool intrafragmentClashes)
{
    vector<sketcherMinimizerAtom*> atoms = molecule->getAtoms();
    vector<sketcherMinimizerBond*> bonds = molecule->getBonds();

    if (atoms.size() > 1) {
        for (sketcherMinimizerAtom* atom : atoms) {
            if (atom->isResidue()) {
                continue;
            }

            for (sketcherMinimizerBond* bond : bonds) {
                if (bond->isResidueInteraction()) {
                    continue;
                }
                sketcherMinimizerAtom* at2 = atom;
                sketcherMinimizerAtom* at1 = bond->startAtom;
                sketcherMinimizerAtom* at3 = bond->endAtom;
                if (at1 == at2 || at1 == at3 || at2 == at3) {
                    continue;
                }
                if (at1->fragment->getDofsOfAtom(at1).empty() &&
                    at2->fragment->getDofsOfAtom(at2).empty() &&
                    at3->fragment->getDofsOfAtom(at3).empty() &&
                    !intrafragmentClashes) {
                    if (at1->fragment == at2->fragment) {
                        continue;
                    }
                    if (at3->fragment == at2->fragment) {
                        continue;
                    }
                }
                if (at2->fixed && at1->fixed && at3->fixed) {
                    continue;
                }

                if (at1->isNeighborOf(at2)) {
                    continue;
                }
                for (sketcherMinimizerAtom* n : at1->neighbors) {
                    if (n->isNeighborOf(at2)) {
                        continue;
                    }
                }
                if (at3->isNeighborOf(at2)) {
                    continue;
                }
                for (sketcherMinimizerAtom* n : at3->neighbors) {
                    if (n->isNeighborOf(at2)) {
                        continue;
                    }
                }

                if (!(at1->rigid && at2->rigid && at3->rigid)) {

                    //                }
                    auto* interaction =
                        new sketcherMinimizerClashInteraction(at1, at2, at3);
                    float restVK = 0.8f;
                    if (at2->atomicNumber == 6 && at2->charge == 0) {
                        restVK -= 0.1f;
                    }
                    if (at1->atomicNumber == 6 && at1->charge == 0 &&
                        at3->atomicNumber == 6 && at3->charge == 0) {
                        restVK -= 0.1f;
                    }

                    interaction->restV =
                        (bondLength * restVK) * (bondLength * restVK);
                    _intramolecularClashInteractions.push_back(interaction);
                    _interactions.push_back(interaction);
                }
            }
        }
    }
}

void CoordgenMinimizer::addStretchInteractionsOfMolecule(
                                sketcherMinimizerMolecule* molecule)
{
    vector<sketcherMinimizerBond*> bonds = molecule->getBonds();
    for (sketcherMinimizerBond* bo : bonds) {
        if (bo->isResidueInteraction()) {
            continue;
        }
        sketcherMinimizerAtom* at1 = bo->startAtom;
        sketcherMinimizerAtom* at2 = bo->endAtom;
        auto* interaction = new sketcherMinimizerStretchInteraction(at1, at2);
        interaction->k *= 0.1f;
        interaction->restV = bondLength;
        if (at1->rigid && at2->rigid) {
            sketcherMinimizerPointF v = at2->coordinates - at1->coordinates;
            interaction->restV = v.length();
        }
        auto sharedRing = sketcherMinimizer::sameRing(at1, at2);
        if (sharedRing && !sharedRing->isMacrocycle()) {
            interaction->k *= 50;
        }
        _interactions.push_back(interaction);
        _stretchInteractions.push_back(interaction);
    }
}

/* return a set of all carbons part of a carbonyl group, i.e. doubly bonded to
 * an oxygen. */
std::set<sketcherMinimizerAtom*> CoordgenMinimizer::getChetoCs(
    const std::vector<sketcherMinimizerAtom*>& allAtoms) const
{
    std::set<sketcherMinimizerAtom*> chetoCs;
    for (auto atom : allAtoms) {
        if (atom->atomicNumber != 6) {
            continue;
        }
        for (auto bondedAtom : atom->neighbors) {
            if (bondedAtom->atomicNumber == 8) {
                auto bond = sketcherMinimizer::getBond(atom, bondedAtom);
                if (bond && bond->bondOrder == 2) {
                    chetoCs.insert(atom);
                    continue;
                }
            }
        }
    }
    return chetoCs;
}

/* return a set of all amino nitrogens. Not chemically accurate, doesn't filter
 * out nitro Ns for instance. */
std::set<sketcherMinimizerAtom*> CoordgenMinimizer::getAminoNs(
    const std::vector<sketcherMinimizerAtom*>& allAtoms) const
{
    std::set<sketcherMinimizerAtom*> aminoNs;
    for (auto atom : allAtoms) {
        if (atom->atomicNumber == 7) {
            aminoNs.insert(atom);
        }
    }
    return aminoNs;
}

/* return a set of all aminoacid alpha carbon, i.e. a carbon that is bound to a
 nitrogen
 and a cheto carbon. */
std::set<sketcherMinimizerAtom*> CoordgenMinimizer::getAlphaCs(
    const std::vector<sketcherMinimizerAtom*>& allAtoms,
    const std::set<sketcherMinimizerAtom*>& chetoCs,
    const std::set<sketcherMinimizerAtom*>& aminoNs) const
{
    std::set<sketcherMinimizerAtom*> alphaCs;
    for (auto atom : allAtoms) {
        bool bondedToCheto = false;
        bool bondedToAminoN = false;
        if (atom->atomicNumber != 6) {
            continue;
        }
        if (chetoCs.find(atom) != chetoCs.end()) {
            continue;
        }
        for (auto bondedAtom : atom->neighbors) {
            if (chetoCs.find(bondedAtom) != chetoCs.end()) {
                bondedToCheto = true;
            }
            if (aminoNs.find(bondedAtom) != aminoNs.end()) {
                bondedToAminoN = true;
            }
        }
        if (bondedToCheto && bondedToAminoN) {
            alphaCs.insert(atom);
        }
    }
    return alphaCs;
}

/* add constraints to keep all backbone atoms of a peptide in a straight trans
 * line. */
void CoordgenMinimizer::addPeptideBondInversionConstraintsOfMolecule(
    sketcherMinimizerMolecule* molecule)
{
    auto atoms = molecule->getAtoms();
    auto chetoCs = getChetoCs(atoms);
    if (chetoCs.size() < 2) {
        return;
    }
    auto aminoNs = getAminoNs(atoms);
    if (aminoNs.size() < 2) {
        return;
    }
    auto alphaCs = getAlphaCs(atoms, chetoCs, aminoNs);
    if (alphaCs.size() < 2) {
        return;
    }
    std::vector<std::vector<sketcherMinimizerAtom*>> consecutiveAtomsGroups;
    getFourConsecutiveAtomsThatMatchSequence(consecutiveAtomsGroups, chetoCs,
                                             aminoNs, alphaCs, chetoCs);
    getFourConsecutiveAtomsThatMatchSequence(consecutiveAtomsGroups, aminoNs,
                                             alphaCs, chetoCs, aminoNs);
    getFourConsecutiveAtomsThatMatchSequence(consecutiveAtomsGroups, alphaCs,
                                             chetoCs, aminoNs, alphaCs);
    for (auto torsionAtoms : consecutiveAtomsGroups) {
        bool cis = false;
        auto interaction = new sketcherMinimizerEZConstrainInteraction(
            torsionAtoms[0], torsionAtoms[1], torsionAtoms[2], torsionAtoms[3],
            cis);
        _extraInteractions.push_back(interaction);
        _interactions.push_back(interaction);
    }
}

/* find chains of four bound atoms that are part of the four provided sets.
 Useful to detect portions of a peptide backbone for instance. */
void CoordgenMinimizer::getFourConsecutiveAtomsThatMatchSequence(
    std::vector<std::vector<sketcherMinimizerAtom*>>& consecutiveAtomsGroups,
    const std::set<sketcherMinimizerAtom*>& firstSet,
    const std::set<sketcherMinimizerAtom*>& secondSet,
    const std::set<sketcherMinimizerAtom*>& thirdSet,
    const std::set<sketcherMinimizerAtom*>& fourthSet) const
{
    for (auto firstAtom : firstSet) {
        for (auto secondAtom : firstAtom->neighbors) {
            if (secondSet.find(secondAtom) == secondSet.end()) {
                continue;
            }
            for (auto thirdAtom : secondAtom->neighbors) {
                if (thirdSet.find(thirdAtom) == thirdSet.end()) {
                    continue;
                }
                for (auto fourthAtom : thirdAtom->neighbors) {
                    if (fourthSet.find(fourthAtom) == fourthSet.end()) {
                        continue;
                    }
                    std::vector<sketcherMinimizerAtom*> fourMatchingAtoms(4);
                    fourMatchingAtoms.at(0) = firstAtom;
                    fourMatchingAtoms.at(1) = secondAtom;
                    fourMatchingAtoms.at(2) = thirdAtom;
                    fourMatchingAtoms.at(3) = fourthAtom;
                    consecutiveAtomsGroups.push_back(fourMatchingAtoms);
                }
            }
        }
    }
}

void CoordgenMinimizer::addConstrainedInteractionsOfMolecule(
    sketcherMinimizerMolecule* molecule)
{
    for (auto atom : molecule->getAtoms()) {
        if (atom->constrained) {
            auto interaction = new sketcherMinimizerConstraintInteraction(
                atom, atom->templateCoordinates);
            _intramolecularClashInteractions.push_back(interaction);
            _interactions.push_back(interaction);
        }
    }
}

void CoordgenMinimizer::addChiralInversionConstraintsOfMolecule(
    sketcherMinimizerMolecule* molecule)
{
    for (auto ring : molecule->getRings()) {
        if (ring->isMacrocycle()) {
            vector<sketcherMinimizerAtom*> atoms =
                CoordgenFragmentBuilder::orderRingAtoms(ring);
            for (unsigned int i = 0; i < atoms.size(); i++) {
                int size = static_cast<int>(atoms.size());
                int a1 = (i - 1 + size) % size;
                int a11 = (i - 2 + size) % size;
                int a2 = (i + 1) % size;

                sketcherMinimizerBond* bond =
                    sketcherMinimizer::getBond(atoms[a1], atoms[i]);
                if (bond->isStereo()) {
                    bool cis = bond->markedAsCis(atoms[a11], atoms[a2]);
                    auto* ezint = new sketcherMinimizerEZConstrainInteraction(
                        atoms[a11], atoms[a1], atoms[i], atoms[a2], cis);
                    _interactions.push_back(ezint);
                }
            }
        }
    }
}

void CoordgenMinimizer::addBendInteractionsOfMolecule(
    sketcherMinimizerMolecule* molecule)
{
    vector<sketcherMinimizerAtom*> atoms = molecule->getAtoms();
    vector<sketcherMinimizerBond*> bonds = molecule->getBonds();
    for (sketcherMinimizerAtom* at : atoms) {
        vector<sketcherMinimizerBendInteraction*> interactions;
        vector<sketcherMinimizerBendInteraction*> ringInteractions;
        vector<sketcherMinimizerBendInteraction*> nonRingInteractions;
        int nbonds = static_cast<int>(at->neighbors.size());
        bool invertedMacrocycleBond = false;
        if (nbonds > 1) {
            // order bonds so that they appear in clockwise order.
            vector<sketcherMinimizerAtom*> orderedNeighs =
                at->clockwiseOrderedNeighbors();
            float angle = 120.f;
            if (nbonds == 2) {
                if (at->bonds[0]->bondOrder + at->bonds[1]->bondOrder > 3) {
                    angle = 180.f;
                }
            }
            if (nbonds > 2) {
                angle = 360.f / nbonds;
            }
            for (int i = 0; i < nbonds; i++) {
                int j = (i - 1 + nbonds) % nbonds;
                if (nbonds == 2 && i == 1) {
                    continue; // first and last interaction are the same if
                }
                // there is just one interaction
                sketcherMinimizerAtom* at1 = orderedNeighs[i];
                sketcherMinimizerAtom* at2 = at;
                sketcherMinimizerAtom* at3 = orderedNeighs[j];
                auto* interaction =
                    new sketcherMinimizerBendInteraction(at1, at2, at3);
                interactions.push_back(interaction);
                interaction->restV = angle;
                sketcherMinimizerRing* r = sketcherMinimizer::sameRing(
                    at, orderedNeighs[i], orderedNeighs[j]);
                if (r) {
                    if (!r->isMacrocycle()) {
                        int extraAtoms = 0;
                        /* if the rings are to be drawn as fused, they will
                         * result in a bigger ring */
                        for (unsigned int i = 0; i < r->fusedWith.size(); i++) {
                            if (r->fusedWith[i]->isMacrocycle()) {
                                continue;
                            }
                            if (r->fusionAtoms[i].size() > 2) {
                                extraAtoms += static_cast<int>(
                                    r->fusedWith[i]->_atoms.size() -
                                    r->fusionAtoms[i].size());
                            }
                        }
                        interaction->isRing = true;
                        interaction->k *= 100;
                        interaction->restV = static_cast<float>(
                            180. - (360. / (r->size() + extraAtoms)));
                        ringInteractions.push_back(interaction);
                    } else {
                        if (nbonds == 3) {
                            sketcherMinimizerAtom* otherAtom = nullptr;
                            for (auto atom : orderedNeighs) {
                                if (atom != at1 && atom != at3) {
                                    otherAtom = atom;
                                    break;
                                }
                            }
                            if (otherAtom) {
                                if (sketcherMinimizerMaths::sameSide(
                                        at3->coordinates,
                                        otherAtom->coordinates,
                                        at1->coordinates, at2->coordinates)) {
                                    invertedMacrocycleBond = true;
                                }
                            }
                        }
                        bool fusedToRing = false;
                        if (orderedNeighs.size() > 2) {
                            fusedToRing = true;
                        }
                        for (auto atom : orderedNeighs) {
                            if (!sketcherMinimizer::sameRing(at, atom)) {
                                fusedToRing = false;
                                break;
                            }
                        }
                        if (fusedToRing || invertedMacrocycleBond) {
                            ringInteractions.push_back(interaction);
                        } else {
                            /* macrocycles are treated as non rings */
                            nonRingInteractions.push_back(interaction);
                        }
                    }
                } else {
                    nonRingInteractions.push_back(interaction);
                }
                if (interaction->atom1->rigid && interaction->atom2->rigid &&
                    interaction->atom3->rigid) {
                    interaction->restV = sketcherMinimizerMaths::unsignedAngle(
                        interaction->atom1->coordinates,
                        interaction->atom2->coordinates,
                        interaction->atom3->coordinates);
                }
            }
        }
        if (ringInteractions.size() != 1 || nonRingInteractions.size() != 2) {
            invertedMacrocycleBond = false;
        }
        if (!ringInteractions.empty()) { // subtract all the ring angles from 360
                                       // and divide the remaining equally
                                       // between the other interactions
            float totalAngleInRings = 0;
            for (sketcherMinimizerBendInteraction* i : ringInteractions) {
                totalAngleInRings += i->restV;
            }
            if (invertedMacrocycleBond) {
                totalAngleInRings = 360 - totalAngleInRings;
            }
            for (sketcherMinimizerBendInteraction* i : nonRingInteractions) {
                i->restV =
                    (360 - totalAngleInRings) / nonRingInteractions.size();
            }
        } else { // do nothing if 1 or 3 interactions (defaults to 120 degrees)
                 // if 4 or more set angles accordingly

            if (nonRingInteractions.size() == 4) {
                if (at->crossLayout || m_evenAngles) {
                    nonRingInteractions[0]->restV = 90;
                    nonRingInteractions[1]->restV = 90;
                    nonRingInteractions[2]->restV = 90;
                    nonRingInteractions[3]->restV = 90;
                } else {
                    int indexOfBiggestAngle = 0;
                    float biggestAngle = 0;
                    int counter = 0;
                    for (auto interaction : nonRingInteractions) {
                        float angle = sketcherMinimizerMaths::unsignedAngle(
                            interaction->atom1->coordinates,
                            interaction->atom2->coordinates,
                            interaction->atom3->coordinates);
                        if (angle > biggestAngle) {
                            biggestAngle = angle;
                            indexOfBiggestAngle = counter;
                        }
                        counter++;
                    }
                    nonRingInteractions[indexOfBiggestAngle]->restV = 120;
                    nonRingInteractions[(indexOfBiggestAngle + 1) % 4]->restV =
                        90;
                    nonRingInteractions[(indexOfBiggestAngle + 2) % 4]->restV =
                        60;
                    nonRingInteractions[(indexOfBiggestAngle + 3) % 4]->restV =
                        90;
                }
            } else if (nonRingInteractions.size() > 4) {
                for (sketcherMinimizerBendInteraction* i :
                     nonRingInteractions) {
                    i->restV =
                        static_cast<float>(360. / nonRingInteractions.size());
                }
            }
        }
        for (auto interaction : interactions) {
            if (!(interaction->atom1->fixed && interaction->atom2->fixed &&
                  interaction->atom3->fixed)) {
                _interactions.push_back(interaction);
                _bendInteractions.push_back(interaction);
            } else {
                delete interaction;
            }
        }
    }
}

void CoordgenMinimizer::minimizeMolecule(sketcherMinimizerMolecule* molecule)
{
    std::map<sketcherMinimizerAtom*, sketcherMinimizerPointF>
        previousCoordinates;
    for (auto atom : molecule->getAtoms()) {
        previousCoordinates[atom] = atom->getCoordinates();
    }
    clearInteractions();
    addInteractionsOfMolecule(molecule, true);
    run();
    for (auto bond : molecule->getBonds()) {
        if (!bond->checkStereoChemistry()) {
            for (auto atom : molecule->getAtoms()) {
                atom->setCoordinates(previousCoordinates[atom]);
            }
            break;
        }
    }
}

void CoordgenMinimizer::minimizeResidues()
{
    setupInteractionsOnlyResidues();
    run();
}

void CoordgenMinimizer::minimizeProteinOnlyLID(
    const std::map<std::string, std::vector<sketcherMinimizerResidue*>>& chains)
{
    setupInteractionsProteinOnly(chains);
    run();
}

void CoordgenMinimizer::minimizeAll()
{
    setupInteractions(true);
    run();
}

void CoordgenMinimizer::addInteractionsOfMolecule(
                                                  sketcherMinimizerMolecule* molecule, bool intrafragmentClashes)
{
    addClashInteractionsOfMolecule(molecule, intrafragmentClashes);
    addStretchInteractionsOfMolecule(molecule);
    addBendInteractionsOfMolecule(molecule);
    addChiralInversionConstraintsOfMolecule(molecule);
}

void CoordgenMinimizer::setupInteractionsOnlyResidues()
{
    const float CLASH_DISTANCE = bondLength * 1.5f;
    for (auto res : m_residues) {
        for (auto res2 : m_residues) {
            if (res2 >= res) {
                continue;
            }
            auto* minimizerInteraction =
                new sketcherMinimizerClashInteraction(res, res2, res);
            minimizerInteraction->restV = CLASH_DISTANCE * CLASH_DISTANCE;
            _interactions.push_back(minimizerInteraction);
        }
    }
}

void CoordgenMinimizer::setupInteractionsProteinOnly(
    const std::map<std::string, std::vector<sketcherMinimizerResidue*>>& chains)
{
    clearInteractions();
    std::set<sketcherMinimizerBond*> interactions;
    std::set<sketcherMinimizerResidue*> residues;
    for (const auto& chain : chains) {
        for (auto res : chain.second) {
            residues.insert(res);
            for (auto interaction : res->residueInteractions) {
                interactions.insert(interaction);
            }
        }
    }
    for (auto res : residues) {
        for (auto interaction : interactions) {
            if (res == interaction->startAtom || res == interaction->endAtom) {
                continue;
            }
            auto* minimizerInteraction = new sketcherMinimizerClashInteraction(
                interaction->startAtom, res, interaction->endAtom);
            minimizerInteraction->restV = bondLength * bondLength;
            _interactions.push_back(minimizerInteraction);
        }
    }
}

void CoordgenMinimizer::setupInteractions(bool intrafragmentClashes)
{
    clearInteractions();
    for (sketcherMinimizerMolecule* molecule : m_molecules) {
        addInteractionsOfMolecule(molecule, intrafragmentClashes);
    }
}

float CoordgenMinimizer::scoreInteractions()
{
    float totalEnergy = 0.f;
    for (auto interaction : _interactions) {
        interaction->score(totalEnergy);
    }
    return totalEnergy;
}

// returns true if the two molecules have a atom-atoms or atom-bond pair with
// distance < threshold or crossing bonds
bool CoordgenMinimizer::findIntermolecularClashes(
    sketcherMinimizerMolecule* mol1, sketcherMinimizerMolecule* mol2,
    float threshold)
{
    // could be made faster for instance checking the molecules bounding boxes
    // first
    if (mol1 == mol2) {
        return false;
    }
    float threshold2 = threshold * threshold;
    for (sketcherMinimizerAtom* a : mol1->_atoms) {
        for (sketcherMinimizerAtom* a2 : mol2->_atoms) {
            if (sketcherMinimizerMaths::squaredDistance(
                    a->coordinates, a2->coordinates) < threshold2) {
                return true;
            }
        }
    }

    for (sketcherMinimizerAtom* a : mol1->_atoms) {
        for (sketcherMinimizerBond* b : mol2->_bonds) {
            if (sketcherMinimizerMaths::squaredDistancePointSegment(
                    a->coordinates, b->startAtom->coordinates,
                    b->endAtom->coordinates) < threshold2) {
                return true;
            }
        }
    }
    for (sketcherMinimizerAtom* a : mol2->_atoms) {
        for (sketcherMinimizerBond* b : mol1->_bonds) {
            if (sketcherMinimizerMaths::squaredDistancePointSegment(
                    a->coordinates, b->startAtom->coordinates,
                    b->endAtom->coordinates) < threshold2) {
                return true;
            }
        }
    }
    for (sketcherMinimizerBond* b : mol1->_bonds) {
        for (sketcherMinimizerBond* b2 : mol2->_bonds) {
            if (sketcherMinimizerMaths::intersectionOfSegments(
                    b->startAtom->coordinates, b->endAtom->coordinates,
                    b2->startAtom->coordinates, b2->endAtom->coordinates)) {
                return true;
            }
        }
    }

    return false;
}

bool CoordgenMinimizer::findIntermolecularClashes(
    const vector<sketcherMinimizerMolecule*>& mols, float threshold)
{
    for (unsigned int i = 0; i < mols.size(); i++) {
        for (unsigned int j = i + 1; j < mols.size(); j++) {
            if (findIntermolecularClashes(mols[i], mols[j], threshold)) {
                return true;
            }
        }
    }
    return false;
}

float CoordgenMinimizer::scoreClashes(
    sketcherMinimizerMolecule* molecule, bool residueInteractions,
    bool scoreProximityRelationsOnOppositeSid) const
{
    float E = 0.f;
    for (auto i : _intramolecularClashInteractions) {
        i->score(E, true);
    }
    for (sketcherMinimizerInteraction* i : _extraInteractions) {
        i->score(E, true);
    }
    E += scoreDofs(molecule);
    E += scoreCrossBonds(molecule, residueInteractions);
    E += scoreAtomsInsideRings();
    if (scoreProximityRelationsOnOppositeSid) {
        E += scoreProximityRelationsOnOppositeSides();
    }
    return E;
}

float CoordgenMinimizer::scoreDofs(sketcherMinimizerMolecule* molecule) const
{
    float E = 0.f;
    for (const auto& fragment : molecule->getFragments()) {
        for (const auto& dof : fragment->getDofs()) {
            E += dof->getCurrentPenalty();
        }
    }
    return E;
}

float CoordgenMinimizer::scoreCrossBonds(sketcherMinimizerMolecule* molecule,
                                         bool residueInteractions) const
{
    if (!m_scoreResidueInteractions) {
        residueInteractions = false;
    }

    float out = 0.f;
    vector<sketcherMinimizerBond*> bonds = molecule->getBonds();
    if (molecule->getBonds().size() > 2) {

        for (unsigned int b = 0; b < bonds.size() - 1; b++) {
            sketcherMinimizerBond* b1 = bonds[b];
            if (b1->isResidueInteraction()) {
                continue;
            }
            for (unsigned int bb = b + 1; bb < bonds.size(); bb++) {
                sketcherMinimizerBond* b2 = bonds[bb];
                if (b2->isResidueInteraction()) {
                    continue;
                }
                if (b2->startAtom->molecule != b1->startAtom->molecule) {
                    continue;
                }
                if (bondsClash(b1, b2)) {
                    float penalty = STANDARD_CROSSING_BOND_PENALTY *
                        b1->crossingBondPenaltyMultiplier *
                        b2->crossingBondPenaltyMultiplier;
                    if (b1->isTerminal() || b2->isTerminal()) {
                        penalty *= TERMINAL_BOND_CROSSING_MULTIPLIER;
                    }
                    if (b1->isInMacrocycle() || b2->isInMacrocycle()) {
                        penalty *= MACROCYCLE_BOND_CROSSING_MULTIPLIER;
                    }
                    if (b1->isInSmallRing() || b2->isInSmallRing()) {
                        penalty *= RING_BOND_CROSSING_MULTIPLIER;
                    }
                    out += penalty;
                }
            }
        }
    }
    if (!m_residueInteractions.empty() && residueInteractions) {
        for (auto r : m_residues) {
            if (r->residueInteractions.size() > 1) {
                for (unsigned int ri1 = 0;
                     ri1 < r->residueInteractions.size() - 1; ri1++) {
                    for (unsigned int ri2 = 1;
                         ri2 < r->residueInteractions.size(); ri2++) {

                        sketcherMinimizerAtom* a1 =
                            r->residueInteractions[ri1]->endAtom;
                        sketcherMinimizerAtom* a2 =
                            r->residueInteractions[ri2]->endAtom;
                        if (sketcherMinimizerMaths::intersectionOfSegments(
                                a1->coordinates +
                                    a1->getSingleAdditionVector() * 0.2f,
                                a2->coordinates +
                                    a2->getSingleAdditionVector() * 0.2f,
                                a1->coordinates, a2->coordinates)) {
                            out += 15.f;
                        }

                        for (auto b2 : m_bonds) {
                            if (b2->startAtom ==
                                r->residueInteractions[ri1]->endAtom) {
                                continue;
                            }
                            if (b2->endAtom ==
                                r->residueInteractions[ri1]->endAtom) {
                                continue;
                            }
                            if (b2->startAtom ==
                                r->residueInteractions[ri2]->endAtom) {
                                continue;
                            }
                            if (b2->endAtom ==
                                r->residueInteractions[ri2]->endAtom) {
                                continue;
                            }

                            if (sketcherMinimizerMaths::intersectionOfSegments(
                                    a1->coordinates, a2->coordinates,
                                    b2->startAtom->coordinates,
                                    b2->endAtom->coordinates)) {
                                out += 10.f;
                            }
                        }
                    }
                }
            }
        }
    }
    return out;
}

float CoordgenMinimizer::scoreAtomsInsideRings() const
{
    float out = 0.f;
    float cutOff = bondLength;
    for (sketcherMinimizerMolecule* m : m_molecules) {
        for (sketcherMinimizerRing* r : m->_rings) {
            if (r->_atoms.size() > MACROCYCLE) {
                continue;
            }
            if (r->_atoms.size() < 3) {
                continue;
            }
            sketcherMinimizerPointF c = r->findCenter();
            for (sketcherMinimizerAtom* a : m->_atoms) {
                if (a->fragment == r->_atoms[0]->fragment) {
                    continue;
                }
                sketcherMinimizerPointF d = c - a->coordinates;
                if (d.x() > cutOff) {
                    continue;
                }
                if (d.y() > cutOff) {
                    continue;
                }
                if (d.x() < -cutOff) {
                    continue;
                }
                if (d.y() < -cutOff) {
                    continue;
                }
                float sq = d.squareLength();
                if (sq > cutOff * cutOff) {
                    continue;
                }
                float dist = d.length();
                if (dist < cutOff) {
                    out += 50 + 100 * (1 - (dist / cutOff));
                }
            }
        }
    }
    return out;
}

float CoordgenMinimizer::scoreProximityRelationsOnOppositeSides() const
{
    float out = 0.f;
    for (sketcherMinimizerMolecule* m : m_molecules) {
        if (m->_atoms.size() < 2) {
            continue;
        }
        for (unsigned int i = 0; i < m->m_proximityRelations.size(); i++) {
            sketcherMinimizerPointF v1, v2;
            sketcherMinimizerMolecule* otherMol1 = nullptr;
            sketcherMinimizerBond* pr1 = m->m_proximityRelations[i];
            sketcherMinimizerFragment* f1 = nullptr;
            if (pr1->startAtom->molecule == m) {
                f1 = pr1->startAtom->fragment;
                v1 = pr1->startAtom->getSingleAdditionVector();
                otherMol1 = pr1->endAtom->molecule;
            } else {
                f1 = pr1->endAtom->fragment;
                v1 = pr1->endAtom->getSingleAdditionVector();
                otherMol1 = pr1->startAtom->molecule;
            }
            if (otherMol1 == m) {
                continue;
            }

            for (unsigned int j = i + 1; j < m->m_proximityRelations.size();
                 j++) {
                sketcherMinimizerMolecule* otherMol2 = nullptr;

                sketcherMinimizerBond* pr2 = m->m_proximityRelations[j];
                if (pr2->startAtom->molecule == m) {
                    if (pr2->startAtom->fragment == f1) {
                        continue;
                    }
                    v2 = pr2->startAtom->getSingleAdditionVector();
                    otherMol2 = pr2->endAtom->molecule;

                } else {
                    if (pr2->endAtom->fragment == f1) {
                        continue;
                    }
                    v2 = pr2->endAtom->getSingleAdditionVector();
                    otherMol2 = pr2->startAtom->molecule;
                }
                if (otherMol2 == m) {
                    continue;
                }
                if (otherMol1 != otherMol2) {
                    continue;
                }
                float angle = sketcherMinimizerMaths::unsignedAngle(
                    v1, sketcherMinimizerPointF(0.f, 0.f), v2);
                if (angle > 90) {
                    out += SAME_SIDE_DPR_PENALTY +
                           SAME_SIDE_DPR_PENALTY_2 * (angle - 90);
                }
            }
        }
    }
    return out;
}

bool CoordgenMinimizer::runExhaustiveSearch(sketcherMinimizerMolecule* molecule,
                                            vector<CoordgenFragmentDOF*> dofs,
                                            float& clashE,
                                            CoordgenDOFSolutions& solutions)
{
    float bestResult = clashE;
    bool abort = false;
    runExhaustiveSearchLevel(molecule, dofs.begin(), dofs, bestResult, abort,
                             solutions);
    for (auto dof : dofs) {
        dof->setToOptimalValue();
    }
    clashE = bestResult;
    return (bestResult < clashEnergyThreshold);
}

void CoordgenMinimizer::runExhaustiveSearchLevel(
    sketcherMinimizerMolecule* molecule,
    const vector<CoordgenFragmentDOF*>::iterator& iterator,
    vector<CoordgenFragmentDOF*>& dofs, float& bestResult, bool& abort,
    CoordgenDOFSolutions& solutions)
{
    if (abort) {
        return;
    }
    if (iterator == dofs.end()) {
        float result = solutions.scoreCurrentSolution();
        if (result < clashEnergyThreshold) {
            for (auto dof : dofs) {
                dof->storeCurrentValueAsOptimal();
            }
            abort = true;
        } else if (result < bestResult - SKETCHER_EPSILON) {
            bestResult = result;
            for (auto dof : dofs) {
                dof->storeCurrentValueAsOptimal();
            }
        }
    } else {
        vector<CoordgenFragmentDOF*>::iterator nextIter = iterator;
        ++nextIter;
        for (int i = 0; i < (*iterator)->numberOfStates(); ++i) {
            runExhaustiveSearchLevel(molecule, nextIter, dofs, bestResult,
                                     abort, solutions);
            (*iterator)->changeState();
        }
    }
}

std::vector<std::vector<CoordgenFragmentDOF*>>
CoordgenMinimizer::buildTuplesOfDofs(const vector<CoordgenFragmentDOF*>& dofs,
                                     unsigned int order) const
{
    std::vector<std::vector<CoordgenFragmentDOF *>> growingVector,
        lastOrderVector;
    for (auto dof : dofs) {
        std::vector<CoordgenFragmentDOF*> tuple;
        tuple.push_back(dof);
        growingVector.push_back(tuple);
    }
    for (unsigned int i = 1; i < order; ++i) {
        lastOrderVector = growingVector;
        growingVector.clear();
        for (auto lastOrderTuple : lastOrderVector) {
            bool copy = false;
            for (auto dof : dofs) {
                if (copy) {
                    auto newTuple = lastOrderTuple;
                    newTuple.push_back(dof);
                    growingVector.push_back(newTuple);
                } else if (dof == *(lastOrderTuple.rbegin())) {
                    copy = true;
                }
            }
        }
    }
    return growingVector;
}

bool CoordgenMinimizer::growSolutions(
    std::set<std::vector<short unsigned int>>& allScoredSolutions,
    int& currentTier,
    std::map<std::vector<short unsigned int>, float>& growingSolutions,
    CoordgenDOFSolutions& solutions, float& bestScore)
{
    std::map<std::vector<short unsigned int>, float> oldGrowingSolutions =
        growingSolutions;
    float bestScoreForRun = bestScore;
    std::vector<std::pair<float, std::vector<short unsigned int>>>
        bestSolutions;
    bestSolutions.reserve(growingSolutions.size());
    for (const auto& solution : growingSolutions) {
        bestSolutions.emplace_back(solution.second, solution.first);
    }
    sort(bestSolutions.begin(), bestSolutions.end());
    growingSolutions.clear();
    int maxN = static_cast<int>(6 * getPrecision());
    if (maxN < 1) {
        maxN = 1;
    }
    int n = 0;

    for (const auto& solution : bestSolutions) {
        if (n > maxN) {
            break;
        }
        for (auto dof : solutions.getAllDofs()) {
            if (dof->tier() > currentTier) {
                continue;
            }
            solutions.loadSolution(solution.second);
            for (int i = 1; i < dof->numberOfStates(); ++i) {
                dof->changeState();

                auto newSolution = solutions.getCurrentSolution();
                if (allScoredSolutions.find(newSolution) ==
                    allScoredSolutions.end()) {
                    float score = solutions.scoreCurrentSolution();
                    if (score == REJECTED_SOLUTION_SCORE) {
                        return false;
                    }
                    allScoredSolutions.insert(newSolution);
                    if (score < bestScore) {
                        bestScore = score;
                    }
                    if (score < bestScoreForRun &&
                        score < REJECTED_SOLUTION_SCORE) {
                        growingSolutions[newSolution] = score;
                    }
                }
            }
        }
        n++;
    }
    if (growingSolutions.empty() && currentTier < 5) {
        currentTier += 3;
        growingSolutions = oldGrowingSolutions;
    }
    return true;
}

bool CoordgenMinimizer::runSearch(int tier, CoordgenDOFSolutions& solutions)
{
    std::map<std::vector<short unsigned int>, float> growingSolutions;
    std::set<std::vector<short unsigned int>> allScoredSolutions;
    float bestScore = solutions.scoreCurrentSolution();
    growingSolutions[solutions.getCurrentSolution()] = bestScore;
    int i = 0;
    bool hasValidSolution = true;
    do {
#ifdef DEBUG_MINIMIZATION_COORDINATES
        // store data from this minimization step to be written to a file later
        energy_list.push_back(solutions.scoreCurrentSolution());
        std::vector<sketcherMinimizerPointF> these_coordinates;
        for (auto atom : _atoms) {
            these_coordinates.push_back(atom->coordinates);
        }
        all_coordinates.push_back(these_coordinates);
#endif
        ++i;
        hasValidSolution = growSolutions(
            allScoredSolutions, tier, growingSolutions, solutions, bestScore);
    } while ((hasValidSolution && !growingSolutions.empty()) && i < 100);
    std::pair<std::vector<short unsigned int>, float> bestSolution =
        solutions.findBestSolution();
    solutions.loadSolution(bestSolution.first);
    return bestSolution.second < clashEnergyThreshold;
}

bool CoordgenMinimizer::runLocalSearch(sketcherMinimizerMolecule* molecule,
                                       const vector<CoordgenFragmentDOF*>& dofs,
                                       int levels, float& clashE,
                                       CoordgenDOFSolutions& solutions)
{
    bool downhill = false;
    auto combinationsOfDofs = buildTuplesOfDofs(dofs, levels);
    do {
        downhill = false;
        for (const auto& combinationOfDofs : combinationsOfDofs) {
            float lastResult = clashE;
            bool foundOptimalPosition = runExhaustiveSearch(
                molecule, combinationOfDofs, clashE, solutions);
            if (foundOptimalPosition) {
                return true;
            } else if (clashE < lastResult - SKETCHER_EPSILON) {
                downhill = true;
            }
        }
    } while (downhill);
    return false;
}

bool CoordgenMinimizer::flipFragments(sketcherMinimizerMolecule* molecule,
                                      float& clashE)
{
    float bestResult = clashE;
    if (skipFlipFragments) {
        return true;
    }
    if (bestResult < clashEnergyThreshold) {
        return true;
    }
    vector<CoordgenFragmentDOF *> dofs, onlyFlipDofs;
    vector<sketcherMinimizerFragment*> fragments = molecule->getFragments();
    reverse(fragments.begin(), fragments.end());
    for (auto fragment : fragments) {
        if (!fragment->fixed) {
            for (auto dof : fragment->getDofs()) {
                if (dof->numberOfStates() > 1) {
                    dofs.push_back(dof);
                    if (dof == *(fragment->getDofs().begin())) {
                        onlyFlipDofs.push_back(dof);
                    }
                }
            }
        }
    }
    CoordgenDOFSolutions solutions(this, molecule, dofs);
    bool cleanPose = runSearch(0, solutions);
    //  if (!cleanPose) cleanPose = runSearch(6, solutions);
    buildMoleculeFromFragments(molecule, false);
    return cleanPose;
}

bool CoordgenMinimizer::avoidClashesOfMolecule(
    sketcherMinimizerMolecule* molecule,
    const std::vector<sketcherMinimizerInteraction*>& extraInteractions)
{
    clearInteractions();
    addClashInteractionsOfMolecule(molecule, false);
    addPeptideBondInversionConstraintsOfMolecule(molecule);
    for (sketcherMinimizerInteraction* interaction : extraInteractions) {
        _interactions.push_back(interaction);
        _extraInteractions.push_back(interaction);
    }
    for (auto interaction : _extraInteractionsOfMolecule[molecule]) {
        _extraInteractions.push_back(interaction);
        _interactions.push_back(interaction);
    }
    bool scoreResidueInteractions = true;
    bool doNotComputeForces = true;
    float clashE =
        scoreClashes(molecule, scoreResidueInteractions, doNotComputeForces);
    bool cleanPose = flipFragments(molecule, clashE);
    if (!cleanPose) {
        avoidTerminalClashes(molecule, clashE);
        molecule->requireMinimization();
    }
    if (molecule->minimizationIsRequired()) {
        minimizeMolecule(molecule);
    }
    return cleanPose;
}

bool CoordgenMinimizer::avoidClashes()
{
    bool allCleanPoses = true;
    if (skipAvoidClashes) {
        return true;
    }
    for (sketcherMinimizerMolecule* molecule : m_molecules) {
        auto cleanPose = avoidClashesOfMolecule(molecule);
        allCleanPoses = allCleanPoses && cleanPose;
    }
    return allCleanPoses;
}

void CoordgenMinimizer::avoidInternalClashes(
    sketcherMinimizerFragment* fragment)
{
    // avoid intraFragmentClashes
    vector<sketcherMinimizerAtom*> fragmentAtoms = fragment->getAtoms();
    for (sketcherMinimizerAtom* a : fragmentAtoms) {
        if (a->neighbors.size() != 1) {
            continue;
        }
        if (a->needsCheckForClashes) {
            continue;
        }
        if (a->fixed) {
            continue;
        }
        if (!fragment->getDofsOfAtom(a).empty()) {
            continue;
        }
        for (sketcherMinimizerAtom* a2 : fragmentAtoms) {
            if (a == a2) {
                continue;
            }
            if (!fragment->getDofsOfAtom(a2).empty()) {
                continue;
            }
            if (sketcherMinimizer::getBond(a, a2)) {
                continue;
            }
            float dx = a2->coordinates.x() - a->coordinates.x();
            if (dx > bondLength * 0.5f) {
                continue;
            }
            if (dx < -bondLength * 0.5f) {
                continue;
            }
            float dy = a2->coordinates.y() - a->coordinates.y();
            if (dy > bondLength * 0.5f) {
                continue;
            }
            if (dy < -bondLength * 0.5f) {
                continue;
            }
            float squareD = dx * dx + dy * dy;
            if (squareD > bondLength * 0.5f * bondLength * 0.5f) {
                continue;
            }

            sketcherMinimizerPointF vec =
                a->coordinates - a->neighbors[0]->coordinates;
            vec *= 0.3f;
            a->coordinates -= vec;
            if (a2->neighbors.size() == 1) {
                a2->coordinates += vec;
                a2->coordinates.round();
            }
        }
    }
}

bool CoordgenMinimizer::bondsClash(sketcherMinimizerBond* bond,
                                   sketcherMinimizerBond* bond2) const
{
    if (bond == bond2) {
        return false;
    }
    if (bond->getStartAtom() == bond2->getStartAtom() ||
        bond->getStartAtom() == bond2->getEndAtom() ||
        bond->getEndAtom() == bond2->getStartAtom() ||
        bond->getEndAtom() == bond2->getEndAtom()) {
        return false;
    }
    auto& start1 = bond->getStartAtom()->coordinates;
    auto& start2 = bond2->getStartAtom()->coordinates;
    auto& end1 = bond->getEndAtom()->coordinates;
    auto& end2 = bond2->getEndAtom()->coordinates;
    // coincidence and intersection calculations are expensive. Often bonds
    // are nowhere near each other, so skip the remaining work if a bond is
    // strictly to the left or right of another bond.
    if (max(start1.x(), end1.x()) < min(start2.x(), end2.x()) ||
        max(start1.y(), end1.y()) < min(start2.y(), end2.y()) ||
        min(start1.x(), end1.x()) > max(start2.x(), end2.x()) ||
        min(start1.y(), end1.y()) > max(start2.y(), end2.y())) {
        return false;
    }

    if (sketcherMinimizerMaths::pointsCoincide(
            bond->getStartAtom()->coordinates,
            bond2->getStartAtom()->coordinates) ||
        sketcherMinimizerMaths::pointsCoincide(
            bond->getStartAtom()->coordinates,
            bond2->getEndAtom()->coordinates) ||
        sketcherMinimizerMaths::pointsCoincide(
            bond->getEndAtom()->coordinates,
            bond2->getStartAtom()->coordinates) ||
        sketcherMinimizerMaths::pointsCoincide(
            bond->getEndAtom()->coordinates,
            bond2->getEndAtom()->coordinates)) {
        return true;
    }
    return (sketcherMinimizerMaths::intersectionOfSegments(
        bond->startAtom->coordinates, bond->endAtom->coordinates,
        bond2->startAtom->coordinates, bond2->endAtom->coordinates));
}

void CoordgenMinimizer::avoidTerminalClashes(
    sketcherMinimizerMolecule* molecule, float& clashE)
{
    if (clashE < 0.1) {
        return;
    }
    for (auto bond : molecule->getBonds()) {
        if (bond->isResidueInteraction()) {
            continue;
        }
        if (!bond->isTerminal()) {
            continue;
        }
        sketcherMinimizerAtom* terminalAtom = bond->getEndAtom();
        sketcherMinimizerAtom* rootAtom = bond->getStartAtom();
        if (terminalAtom->getBonds().size() != 1) {
            terminalAtom = bond->getStartAtom();
            rootAtom = bond->getEndAtom();
        }
        if (terminalAtom->fixed) {
            continue;
        }
        for (auto bond2 : molecule->getBonds()) {
            if (bond2->isResidueInteraction()) {
                continue;
            }
            if (bondsClash(bond, bond2)) {
                terminalAtom->setCoordinates(rootAtom->getCoordinates() +
                                             (terminalAtom->getCoordinates() -
                                              rootAtom->getCoordinates()) *
                                                 0.1);
            }
        }
    }
    clashE = scoreClashes(molecule);
}

void CoordgenMinimizer::maybeMinimizeRings(
    const vector<sketcherMinimizerRing*>& rings)
{
    bool found = false;
    for (auto r : rings) {
        if (r->_atoms.size() == 5) {
            for (auto& _atom : r->_atoms) {
                if (_atom->rings.size() > 2) {
                    found = true;
                }
            }
        }
        if (r->isMacrocycle() && r->_atoms.size() % 2 != 0) {
            for (auto& _atom : r->_atoms) {
                if (_atom->rings.size() > 2) {
                    found = true;
                }
            }
        }
    }

    if (!found) {
        return;
    }
    rings.at(0)->getAtoms().at(0)->molecule->requireMinimization();
}

void CoordgenMinimizer::buildMoleculeFromFragments(
    sketcherMinimizerMolecule* molecule, bool firstTime) const
{
    for (auto fragment : molecule->getFragments()) {
        float angle = 0;
        sketcherMinimizerPointF position(0.f, 0.f);
        if (fragment->getParent()) {
            sketcherMinimizerPointF p1 =
                fragment->_bondToParent->startAtom->coordinates;
            sketcherMinimizerPointF p2 =
                fragment->_bondToParent->endAtom->coordinates;
            sketcherMinimizerPointF p = p2 - p1;
            angle = atan2(-p.y(), p.x());
            position = fragment->_bondToParent->endAtom->coordinates;
            if (firstTime) {
                sketcherMinimizer::alignWithParentDirection(fragment, position,
                                                            angle);
            }
        }
        fragment->setCoordinates(position, angle);
    }
}

void CoordgenMinimizer::buildFromFragments(bool firstTime) const
{
    for (sketcherMinimizerMolecule* molecule : m_molecules) {
        buildMoleculeFromFragments(molecule, firstTime);
    }
}

bool CoordgenMinimizer::hasValid3DCoordinates(
    const vector<sketcherMinimizerAtom*>& atoms)
{
    for (sketcherMinimizerAtom* atom : atoms) {
        if (!atom->hasValid3DCoordinates()) {
            return false;
        }
    }
    return true;
}

void CoordgenMinimizer::fallbackOn3DCoordinates(
    const vector<sketcherMinimizerAtom*>& atoms)
{
    float scale = 35.f; // ratio between average bond length and 2d bond length
    /* TODO find best projection */
    for (sketcherMinimizerAtom* atom : atoms) {
        atom->setCoordinates(
            sketcherMinimizerPointF(atom->m_x3D * scale, -atom->m_y3D * scale));
    }
}

bool CoordgenMinimizer::hasNaNCoordinates(
    const std::vector<sketcherMinimizerAtom*>& atoms)
{
    for (sketcherMinimizerAtom* a : atoms) {
        if (std::isnan(a->coordinates.x()) || std::isnan(a->coordinates.y())) {
            return true;
        }
    }
    return false;
}

bool CoordgenMinimizer::hasNaNCoordinates()
{
    return hasNaNCoordinates(m_atoms);
}

void CoordgenMinimizer::checkForClashes(sketcherMinimizerAtom* a)
{

    if (a->fixed) {
        return;
    }

    sketcherMinimizerPointF oldCoordinates = a->coordinates;
    vector<sketcherMinimizerPointF> coordsVect;
    coordsVect.push_back(oldCoordinates);
    coordsVect.push_back(oldCoordinates +
                         sketcherMinimizerPointF(bondLength * 0.25f, 0.f));
    coordsVect.push_back(oldCoordinates +
                         sketcherMinimizerPointF(-bondLength * 0.25f, 0.f));
    coordsVect.push_back(oldCoordinates +
                         sketcherMinimizerPointF(0.f, bondLength * 0.25f));
    coordsVect.push_back(oldCoordinates +
                         sketcherMinimizerPointF(0.f, -bondLength * 0.25f));
    coordsVect.push_back(oldCoordinates + sketcherMinimizerPointF(
                                              bondLength * 0.25f * 0.7071f,
                                              -bondLength * 0.25f * 0.7071f));
    coordsVect.push_back(oldCoordinates + sketcherMinimizerPointF(
                                              -bondLength * 0.25f * 0.7071f,
                                              -bondLength * 0.25f * 0.7071f));
    coordsVect.push_back(oldCoordinates +
                         sketcherMinimizerPointF(-bondLength * 0.25f * 0.7071f,
                                                 bondLength * 0.25f * 0.7071f));
    coordsVect.push_back(oldCoordinates +
                         sketcherMinimizerPointF(bondLength * 0.25f * 0.7071f,
                                                 bondLength * 0.25f * 0.7071f));
    float bestE = 999999.f;
    int bestI = 0;
    for (unsigned int i = 0; i < coordsVect.size(); i++) {
        a->coordinates = coordsVect[i];

        // solves intrafragment clashes by shifting the atomic coordinates up
        // down left right or diagonally
        sketcherMinimizerClashInteraction clashI(a, a, a);
        clashI.restV = 300;
        float clashE = 0;
        vector<sketcherMinimizerBond*> bonds = a->getFragment()->getBonds();
        for (sketcherMinimizerBond* b : bonds) {
            if (!b->startAtom->coordinatesSet) {
                continue;
            }
            if (!b->endAtom->coordinatesSet) {
                continue;
            }
            if (b->startAtom == a) {
                continue;
            }
            if (b->endAtom == a) {
                continue;
            }
            clashI.atom1 = b->startAtom;
            clashI.atom2 = a;
            clashI.atom3 = b->endAtom;
            clashI.energy(clashE);
        }

        for (sketcherMinimizerBond* b : a->bonds) {
            vector<sketcherMinimizerAtom*> atoms = a->getFragment()->getAtoms();
            for (sketcherMinimizerAtom* atom : atoms) {
                if (atom == a) {
                    continue;
                }
                if (!b->startAtom->coordinatesSet) {
                    continue;
                }
                if (!b->endAtom->coordinatesSet) {
                    continue;
                }
                if (b->startAtom == atom) {
                    continue;
                }
                if (b->endAtom == atom) {
                    continue;
                }
                clashI.atom1 = b->startAtom;
                clashI.atom2 = atom;
                clashI.atom3 = b->endAtom;
                clashI.energy(clashE);
            }
        }

        if (clashE < SKETCHER_EPSILON) {
            return;
        }
        if (i == 0) {
            bestE = clashE;
        }
        if (clashE < bestE) {
            bestE = clashE;
            bestI = i;
        }
    }
    a->setCoordinates(coordsVect[bestI]);
}

float CoordgenMinimizer::getPrecision() const
{
    return m_precision;
}

void CoordgenMinimizer::setPrecision(float f)
{
    m_precision = f;
}

std::pair<std::vector<short unsigned int>, float>
CoordgenDOFSolutions::findBestSolution() const
{
    std::pair<std::vector<short unsigned int>, float> bestSolution =
        *m_solutions.begin();
    for (auto solution : m_solutions) {
        if (solution.second < bestSolution.second) {
            bestSolution = solution;
        }
    }
    return bestSolution;
}

std::vector<short unsigned int> CoordgenDOFSolutions::getCurrentSolution()
{
    std::vector<short unsigned int> solution;
    for (auto dof : m_allDofs) {
        solution.push_back(dof->getCurrentState());
    }
    return solution;
}

void CoordgenDOFSolutions::loadSolution(
    const std::vector<short unsigned int>& solution)
{
    assert(solution.size() == m_allDofs.size());
    for (unsigned int i = 0; i < solution.size(); ++i) {
        m_allDofs.at(i)->setState(solution.at(i));
    }
}

bool CoordgenDOFSolutions::hasSolution(
    const std::vector<short unsigned int>& solution)
{
    return m_solutions.find(solution) != m_solutions.end();
}

float CoordgenDOFSolutions::scoreCurrentSolution()
{
    std::vector<short unsigned int> solution;
    for (auto dof : m_allDofs) {
        solution.push_back(dof->getCurrentState());
    }
    //   for (auto dof : solution) cerr <<dof;
    //   cerr << endl;
    auto position = m_solutions.find(solution);
    if (position != m_solutions.end()) {
        return position->second;
    } else {
        if (m_solutions.size() >
            MAXIMUM_NUMBER_OF_SCORED_SOLUTIONS * m_minimizer->getPrecision()) {
            return REJECTED_SOLUTION_SCORE;
        }
        m_minimizer->buildMoleculeFromFragments(m_molecule, false);
        float result = m_minimizer->scoreClashes(m_molecule, true);
        m_solutions[solution] = result;
        return result;
    }
}
