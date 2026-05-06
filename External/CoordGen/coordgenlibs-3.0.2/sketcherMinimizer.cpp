/*
 *  sketcherMinimizer.cpp
 *  2d_sketcher
 *
 *  Created by Nicola Zonta on 13/04/2010.
 *  Copyright Schrodinger, LLC. All rights reserved
 *
 */

#include "sketcherMinimizer.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <queue>
#include <stack>

#include "CoordgenTemplates.h"
#ifdef USE_MAEPARSER
#include "sketcherMaeReading.h"
#endif

#include "CoordgenFragmenter.h"
#include "CoordgenMacrocycleBuilder.h"
#include "sketcherMinimizerBendInteraction.h"
#include "sketcherMinimizerClashInteraction.h"
#include "sketcherMinimizerMaths.h"
#include "sketcherMinimizerRing.h"
#include "sketcherMinimizerStretchInteraction.h"

using namespace std;
using namespace schrodinger;

#define RESIDUE_CLASH_DISTANCE_SQUARED 2.0 * 2.0

#ifndef M_PI
#define M_PI 3.1415926535897931
#endif

const int bondLength = BONDLENGTH;

static const unsigned int MINIMUM_LIGAND_ATOMS = 8;
static const float SCORE_MULTIPLIER_FOR_DOUBLE_BONDS = 0.82f;
static const float SCORE_MULTIPLIER_FOR_SINGLE_BONDED_HETEROATOMS = 0.9f;
static const float SCORE_MULTIPLIER_FOR_FRAGMENTS = 0.1f;
const int MAX_NUMBER_OF_RINGS = 40;

sketcherMinimizer::sketcherMinimizer(float precision)
{
    setEvenAngles(false);
    m_minimizer.setPrecision(precision);
    m_fragmentBuilder.setPrecision(precision);
}

sketcherMinimizer::~sketcherMinimizer()
{
    clear();
}

void sketcherMinimizer::setScoreResidueInteractions(bool b)
{
    m_minimizer.m_scoreResidueInteractions = b;
}

void sketcherMinimizer::canonicalOrdering(sketcherMinimizerMolecule* minMol)
{
    vector<int> scores;
    for (unsigned int i = 0; i < minMol->_atoms.size(); i++) {
        minMol->_atoms[i]->_generalUseN = i;
    }

    sketcherMinimizer::morganScores(minMol->_atoms, minMol->_bonds, scores);

    if (scores.size() == minMol->_atoms.size()) {

        for (unsigned int i = 0; i < scores.size(); i++) {
            scores[i] *= 100;
            scores[i] += minMol->_atoms[i]->atomicNumber;
        }
        for (sketcherMinimizerAtom* at : minMol->_atoms) {
            at->neighbors.clear();
            at->bonds.clear();
        }
        for (sketcherMinimizerBond* bo : minMol->_bonds) {
            bo->startAtom->neighbors.push_back(bo->endAtom);
            bo->endAtom->neighbors.push_back(bo->startAtom);
            bo->startAtom->bonds.push_back(bo);
            bo->endAtom->bonds.push_back(bo);
        }

        vector<sketcherMinimizerAtom*> newAtoms;
        vector<sketcherMinimizerBond*> newBonds;

        for (unsigned int i = 0; i < minMol->_atoms.size(); i++) {
            minMol->_atoms[i]->_generalUseN = i;
            minMol->_atoms[i]->_generalUseVisited = false;
        }
        for (auto& _bond : minMol->_bonds) {
            _bond->_SSSRVisited = false;
        }
        bool found = true;
        do {
            int scoreMaxI = -1;
            for (unsigned int i = 0; i < scores.size(); i++) {
                if (minMol->_atoms[i]->_generalUseVisited) {
                    continue;
                }
                if (scoreMaxI == -1) {
                    scoreMaxI = i;
                } else if (scores[i] > scores[scoreMaxI]) {
                    scoreMaxI = i;
                }
            }
            if (scoreMaxI > -1) {
                queue<sketcherMinimizerAtom*> q;
                q.push(minMol->_atoms[scoreMaxI]);
                minMol->_atoms[scoreMaxI]->_generalUseVisited = true;

                while (!q.empty()) {
                    sketcherMinimizerAtom* at = q.front();
                    newAtoms.push_back(at);
                    q.pop();
                    int neighI = -1;
                    do {
                        neighI = -1;
                        for (unsigned int i = 0; i < at->neighbors.size();
                             i++) {
                            if (at->bonds[i]->_SSSRVisited) {
                                continue;
                            } else {
                                if (neighI == -1) {
                                    neighI = i;
                                } else {
                                    if (scores[at->neighbors[neighI]
                                                   ->_generalUseN] <
                                        scores[at->neighbors[i]
                                                   ->_generalUseN]) {
                                        neighI = i;
                                    }
                                }
                            }
                        }
                        if (neighI > -1) {
                            if (!at->neighbors[neighI]->_generalUseVisited) {
                                at->neighbors[neighI]->_generalUseVisited =
                                    true;
                                q.push(at->neighbors[neighI]);
                            }
                            at->bonds[neighI]->_SSSRVisited = true;
                            newBonds.push_back(at->bonds[neighI]);
                        }
                    } while (neighI > -1);
                }
            } else {
                found = false;
            }
        } while (found);

        minMol->_atoms = newAtoms;
        minMol->_bonds = newBonds;
    }
}

void sketcherMinimizer::initialize(
    sketcherMinimizerMolecule* minMol) // min mol is split into molecules if
                                       // needed and then added to the minimizer
{
    clear();
    m_referenceAtoms = minMol->_atoms;
    m_referenceBonds = minMol->_bonds;

    std::map<sketcherMinimizerAtom*, int> bondsToAtom;
    for (auto& _bond : minMol->_bonds) {
        bondsToAtom[_bond->startAtom]++;
        bondsToAtom[_bond->endAtom]++;
    }
    for (auto& _bond : minMol->_bonds) {
        if (_bond->skip) {
            continue;
        }
        if (getTreatNonterminalBondsToMetalAsZOBs()) {
            if (_bond->bondOrder == 1 || _bond->bondOrder == 2) {
                bool terminalBond = bondsToAtom[_bond->startAtom] == 1 ||
                                    bondsToAtom[_bond->endAtom] == 1;
                if (!terminalBond && (sketcherMinimizerAtom::isMetal(
                                          _bond->startAtom->atomicNumber) ||
                                      sketcherMinimizerAtom::isMetal(
                                          _bond->endAtom->atomicNumber))) {
                    _bond->bondOrder = 0;
                }
            }
        }
    }

    for (auto& _bond : minMol->_bonds) {
        if (_bond->skip) {
            continue;
        }
        if (_bond->bondOrder == 0) {
            m_proximityRelations.push_back(_bond);
        } else if (_bond->isResidueInteraction()) {
            if (!_bond->startAtom->isResidue() &&
                !_bond->endAtom->isResidue()) {
                m_proximityRelations.push_back(_bond);
            }
        }
    }
    for (auto& m_extraBond : m_extraBonds) {
        if (m_extraBond->skip) {
            continue;
        }
        if (m_extraBond->bondOrder == 0) {
            m_proximityRelations.push_back(m_extraBond);
        } else if (m_extraBond->isResidueInteraction()) {
            if (!m_extraBond->startAtom->isResidue() &&
                !m_extraBond->endAtom->isResidue()) {
                m_proximityRelations.push_back(m_extraBond);
            }
        }
    }

    {
        // remove skipped and hidden bonds
        auto new_end = std::remove_if(
            minMol->_bonds.begin(), minMol->_bonds.end(),
            [](sketcherMinimizerBond* b) {
                return (b->skip || b->bondOrder == 0 || b->startAtom->hidden ||
                        b->endAtom->hidden);
            });
        minMol->_bonds.erase(new_end, minMol->_bonds.end());
    }
    {
        // remove hidden atoms atoms
        auto new_end = std::remove_if(
            minMol->_atoms.begin(), minMol->_atoms.end(),
            [](sketcherMinimizerAtom* a) { return (a->hidden); });
        minMol->_atoms.erase(new_end, minMol->_atoms.end());
    }

    // order atoms and bonds using morgan indices to make the result input
    // order independent
    canonicalOrdering(minMol);

    for (auto a : minMol->_atoms) {
        if (!a->hidden) {
            m_atoms.push_back(a);
        }
        if (a->isResidue()) {
            m_residues.push_back(static_cast<sketcherMinimizerResidue*>(a));
        }
    }

    for (auto b : minMol->_bonds) {
        if (!b->startAtom->hidden && !b->endAtom->hidden) {
            m_bonds.push_back(b);
        }
        if (b->isResidueInteraction()) {
            m_residueInteractions.push_back(
                static_cast<sketcherMinimizerResidueInteraction*>(b));
        }
    }

    minMol->forceUpdateStruct(minMol->_atoms, minMol->_bonds, minMol->_rings);
    splitIntoMolecules(minMol, m_molecules);

    for (auto& b : m_proximityRelations) {
        b->startAtom->molecule->m_proximityRelations.push_back(b);
        if (b->endAtom != b->startAtom) {
            b->endAtom->molecule->m_proximityRelations.push_back(b);
        }
    }

    flagCrossAtoms();

    m_minimizer.m_atoms = m_atoms;
    m_minimizer.m_bonds = m_bonds;
    m_minimizer.m_molecules = m_molecules;
    m_minimizer.m_residues = m_residues;
    m_minimizer.m_residueInteractions = m_residueInteractions;
}

bool sketcherMinimizer::structurePassSanityCheck() const
{
    if (m_atoms.empty()) {
        return false;
    }
    for (auto molecule : m_molecules) {
        if (molecule->_rings.size() > MAX_NUMBER_OF_RINGS) {
            return false;
        }
    }
    return true;
}

bool sketcherMinimizer::runGenerateCoordinates()
{
    bool cleanPose = true;
    if (structurePassSanityCheck()) {
        findFragments();
        buildFromFragments(true);
        cleanPose = m_minimizer.avoidClashes();
        bestRotation();
        maybeFlip();
        arrangeMultipleMolecules();
        writeStereoChemistry();
#ifdef DEBUG_MINIMIZATION_COORDINATES
        // write minimization data and atom mapping to file
        writeMinimizationData();
#endif
    }
    return cleanPose;
}

void sketcherMinimizer::flagCrossAtoms()
{
    for (auto at : m_atoms) {
        if (at->atomicNumber == 16 || at->atomicNumber == 15) {
            at->crossLayout = true;
        }
    }

    for (auto at : m_atoms) {
        if (at->crossLayout) {
            continue;
        }
        int cross = 0;
        for (auto n : at->neighbors) {
            if (n->neighbors.size() > 3) {
                cross++;
            }
        }
        if (cross > 2) {
            at->crossLayout = true;
        }
    }
}

void sketcherMinimizer::clear()
{
    for (auto& _referenceAtom : m_referenceAtoms) {
        delete _referenceAtom;
    }
    m_referenceAtoms.clear();
    m_residues.clear();

    for (auto& _referenceBond : m_referenceBonds) {
        delete _referenceBond;
    }

    m_referenceBonds.clear();

    for (auto& m_extraBond : m_extraBonds) {
        delete m_extraBond;
    }

    m_extraBonds.clear();

    for (auto& _fragment : _fragments) {
        delete _fragment;
    }

    _fragments.clear();

    for (auto& _molecule : m_molecules) {

        delete _molecule;
    }

    m_molecules.clear();
}

void sketcherMinimizer::splitIntoMolecules(
    sketcherMinimizerMolecule* mol, vector<sketcherMinimizerMolecule*>& mols)
{

    if (mol->_atoms.empty()) {
        mols.push_back(mol);
        return;
    }
    for (sketcherMinimizerAtom* a : mol->_atoms) {
        a->_generalUseVisited = false;
    }
    queue<sketcherMinimizerAtom*> q;
    q.push(mol->_atoms[0]);
    for (sketcherMinimizerAtom* a : mol->_atoms) {
        if (!a->hidden) {
            q.push(a);
            break;
        }
    }
    while (!q.empty()) {
        sketcherMinimizerAtom* a = q.front();
        q.pop();
        a->_generalUseVisited = true;
        for (sketcherMinimizerAtom* n : a->neighbors) {
            if (!n->_generalUseVisited && !n->hidden) {
                q.push(n);
            }
        }
    }
    vector<sketcherMinimizerAtom*> newAtoms;

    for (sketcherMinimizerAtom* a : mol->_atoms) {
        if (!a->_generalUseVisited && !a->hidden) {
            newAtoms.push_back(a);
        }
    }
    if (newAtoms.empty()) {
        mols.push_back(mol);
        for (sketcherMinimizerMolecule* m : mols) {
            for (sketcherMinimizerAtom* a : m->_atoms) {
                a->_generalUseVisited = false;
            }
        }

    } else {
        auto* newMol = new sketcherMinimizerMolecule;
        for (unsigned int i = 0; i < mol->_rings.size(); i++) {

            if (!mol->_rings[i]->_atoms[0]->_generalUseVisited) {
                newMol->_rings.push_back(mol->_rings[i]);
                mol->_rings.erase(mol->_rings.begin() + i);
                i--;
            }
        }
        for (unsigned int i = 0; i < mol->_bonds.size(); i++) {
            if (!mol->_bonds[i]->startAtom->_generalUseVisited) {
                newMol->_bonds.push_back(mol->_bonds[i]);
                mol->_bonds.erase(mol->_bonds.begin() + i);
                i--;
            }
        }

        for (unsigned int i = 0; i < mol->_atoms.size(); i++) {
            if (!mol->_atoms[i]->_generalUseVisited) {
                mol->_atoms[i]->molecule = newMol;
                newMol->_atoms.push_back(mol->_atoms[i]);
                mol->_atoms.erase(mol->_atoms.begin() + i);
                i--;
            }
        }
        mols.push_back(mol);
        splitIntoMolecules(newMol, mols);
    }
}

sketcherMinimizerRing*
sketcherMinimizer::sameRing(const sketcherMinimizerAtom* at1,
                            const sketcherMinimizerAtom* at2)
{
    return sketcherMinimizerAtom::shareARing(at1, at2);
}

sketcherMinimizerRing*
sketcherMinimizer::sameRing(const sketcherMinimizerAtom* at1,
                            const sketcherMinimizerAtom* at2,
                            const sketcherMinimizerAtom* at3)
{
    if (at1->rings.empty()) {
        return nullptr;
    }
    if (at2->rings.empty()) {
        return nullptr;
    }
    if (at3->rings.empty()) {
        return nullptr;
    }
    sketcherMinimizerRing* r = nullptr;
    for (sketcherMinimizerRing* ring : at1->rings) {
        if (ring->isMacrocycle()) {
            continue;
        }
        for (sketcherMinimizerRing* ring2 : at2->rings) {
            if (ring != ring2) {
                continue;
            }
            for (sketcherMinimizerRing* ring3 : at3->rings) {
                if (ring3 == ring2) {
                    if (!r) {
                        r = ring2;
                    } else if (ring2->_atoms.size() < r->_atoms.size()) {
                        r = ring2;
                    }
                }
            }
        }
    }
    for (sketcherMinimizerRing* ring : at1->rings) {
        for (sketcherMinimizerRing* ring2 : at2->rings) {
            if (ring != ring2) {
                continue;
            }
            for (sketcherMinimizerRing* ring3 : at3->rings) {
                if (ring3 == ring2) {
                    if (!r) {
                        r = ring2;
                    } else if (ring2->_atoms.size() < r->_atoms.size()) {
                        r = ring2;
                    }
                }
            }
        }
    }
    return r;
}

void sketcherMinimizer::writeStereoChemistry()
{
    for (sketcherMinimizerAtom* a : m_atoms) {
        if (a->hasStereochemistrySet) {
            a->writeStereoChemistry();
        }
    }
    assignPseudoZ();
}

void sketcherMinimizer::assignPseudoZ()
{
    for (sketcherMinimizerMolecule* mol : m_molecules) {
        for (sketcherMinimizerAtom* a : mol->_atoms) {
            a->_generalUseVisited = false;
        }
        sketcherMinimizerAtom* lastAtom = nullptr;
        bool finished = false;
        while (!finished) {
            lastAtom = nullptr;
            for (sketcherMinimizerAtom* a : mol->_atoms) {
                if (!a->_generalUseVisited) {
                    lastAtom = a;
                    break;
                }
            }
            if (lastAtom) {
                queue<sketcherMinimizerAtom*> q;
                q.push(lastAtom);
                while (!q.empty()) {
                    lastAtom = q.front();
                    q.pop();
                    lastAtom->_generalUseVisited = true;
                    for (unsigned int i = 0; i < lastAtom->neighbors.size();
                         i++) {
                        if (lastAtom->neighbors[i]->_generalUseVisited) {
                            continue;
                        }
                        float Z = lastAtom->m_pseudoZ;
                        sketcherMinimizerBond* b = lastAtom->bonds[i];
                        if (b->hasStereochemistryDisplay) {
                            if (b->isWedge) {
                                if ((b->startAtom == lastAtom &&
                                     !b->isReversed) ||
                                    (b->endAtom == lastAtom && b->isReversed)) {
                                    Z += 1.f;
                                } else if ((b->startAtom == lastAtom &&
                                            b->isReversed) ||
                                           (b->endAtom == lastAtom &&
                                            !b->isReversed)) {
                                    Z -= 1.f;
                                }

                            } else {
                                if ((b->startAtom == lastAtom &&
                                     !b->isReversed) ||
                                    (b->endAtom == lastAtom && b->isReversed)) {
                                    Z -= 1.f;
                                } else if ((b->startAtom == lastAtom &&
                                            b->isReversed) ||
                                           (b->endAtom == lastAtom &&
                                            !b->isReversed)) {
                                    Z += 1.f;
                                }
                            }
                        }
                        lastAtom->neighbors[i]->m_pseudoZ = Z;
                        q.push(lastAtom->neighbors[i]);
                    }
                }
            } else {
                finished = true;
            }
        }
    }
}

void sketcherMinimizer::maybeFlipPeptides(
    const std::vector<sketcherMinimizerAtom*>& atoms, float& scoreX)
{
    auto chetoCs = m_minimizer.getChetoCs(atoms);
    auto aminoNs = m_minimizer.getAminoNs(atoms);
    auto alphaCs = m_minimizer.getAlphaCs(atoms, chetoCs, aminoNs);
    for (auto alphaC : alphaCs) {
        sketcherMinimizerAtom* aminoN = nullptr;
        sketcherMinimizerAtom* chetoC = nullptr;
        for (auto neighbor : alphaC->neighbors) {
            if (aminoNs.find(neighbor) != aminoNs.end()) {
                aminoN = neighbor;
            } else if (chetoCs.find(neighbor) != chetoCs.end()) {
                chetoC = neighbor;
            }
        }
        if (aminoN && chetoC) {
            auto direction = aminoN->coordinates - chetoC->coordinates;
            const float PEPTIDE_SCORE = 100.f;
            if (direction.x() > 0) {
                scoreX -= PEPTIDE_SCORE;
            } else {
                scoreX += PEPTIDE_SCORE;
            }
        }
    }
}

void sketcherMinimizer::maybeFlip()
{

    for (sketcherMinimizerMolecule* mol : m_molecules) {
        if (mol->hasFixedFragments) {
            continue;
        }
        if (mol->hasConstrainedFragments) {
            continue;
        }
        if (mol->_atoms.size() < 2) {
            continue;
        }

        float scoreY = 0.f, scoreX = 0.f;
        maybeFlipPeptides(mol->getAtoms(), scoreX);
        sketcherMinimizerPointF cent(0.f, 0.f);
        for (sketcherMinimizerAtom* a : mol->_atoms)
            cent += a->coordinates;
        if (!mol->_atoms.empty()) {
            cent /= mol->_atoms.size();
        }

        for (sketcherMinimizerFragment* f : mol->_fragments) {
            vector<sketcherMinimizerRing*> rings = f->getRings();
            if (rings.size() == 2) {
                sketcherMinimizerRing* r1 = rings[0];
                sketcherMinimizerRing* r2 = rings[1];
                sketcherMinimizerPointF c1 = r1->findCenter();
                sketcherMinimizerPointF c2 = r2->findCenter();
                if (c1.x() - c2.x() > SKETCHER_EPSILON) {
                    sketcherMinimizerPointF swapP = c2;
                    c2 = c1;
                    c1 = swapP;
                    sketcherMinimizerRing* swapR;
                    swapR = r2;
                    r2 = r1;
                    r1 = swapR;
                }
                if (r2->isBenzene() && !r1->isBenzene()) {
                    scoreX -= 20;
                } else if (r1->isBenzene() && !r2->isBenzene()) {
                    scoreX += 20;
                }
            }

            if (rings.size() > 3) {

                sketcherMinimizerPointF center(0.f, 0.f);
                sketcherMinimizerPointF weightedCenter(0.f, 0.f);

                size_t totalN = 0;
                int totalRings = 0;
                for (sketcherMinimizerRing* r : rings) {
                    if (r->_atoms.size() < 4) {
                        continue;
                    }

                    sketcherMinimizerPointF c = r->findCenter();
                    center += c;
                    weightedCenter += c * r->_atoms.size();
                    totalN += r->_atoms.size();
                    totalRings++;
                }
                if (totalRings && totalN) {
                    center /= totalRings;
                    weightedCenter /= totalN;
                }
                if (weightedCenter.y() - center.y() < -SKETCHER_EPSILON) {
                    scoreY += 50.f;
                } else if (weightedCenter.y() - center.y() > SKETCHER_EPSILON) {
                    scoreY -= 50.f;
                }

                if (weightedCenter.x() - center.x() < -SKETCHER_EPSILON) {
                    scoreX += 50.f;
                } else if (weightedCenter.x() - center.x() > SKETCHER_EPSILON) {
                    scoreX -= 50.f;
                }
            }
        }

        float minx = 9999.f, miny = 9999.f, maxx = -9999.f, maxy = -9999.f;

        for (sketcherMinimizerAtom* a : mol->_atoms) {
            float x = a->coordinates.x();
            float y = a->coordinates.y();
            if (x < minx) {
                minx = x;
            } else if (x > maxx) {
                maxx = x;
            }
            if (y < miny) {
                miny = y;
            } else if (y > maxy) {
                maxy = y;
            }
        }

        float meanx = (maxx + minx) * 0.5f;
        float meany = (maxy + miny) * 0.5f;

        if (meanx - cent.x() > SKETCHER_EPSILON) {
            scoreX -= 0.5f;
        } else if (meanx - cent.x() < -SKETCHER_EPSILON) {
            scoreX += 0.5f;
        }
        if (meany - cent.y() > SKETCHER_EPSILON) {
            scoreY += 0.5f;
        } else if (meany - cent.y() < -SKETCHER_EPSILON) {
            scoreY -= 0.5f;
        }

        for (sketcherMinimizerBond* b : mol->_bonds) {
            if (b->bondOrder == 2) {

                if (b->startAtom->neighbors.size() == 1 &&
                    b->endAtom->neighbors.size() > 1) {

                    float diff = b->startAtom->coordinates.y() -
                                 b->endAtom->coordinates.y();
                    if (diff > SKETCHER_EPSILON) {
                        scoreY += 1;
                    } else if (diff < -SKETCHER_EPSILON) {
                        scoreY -= 1;
                    }
                } else if (b->endAtom->neighbors.size() == 1 &&
                           b->startAtom->neighbors.size() > 1) {

                    float diff = b->endAtom->coordinates.y() -
                                 b->startAtom->coordinates.y();
                    if (diff > SKETCHER_EPSILON) {
                        scoreY += 1;
                    } else if (diff < -SKETCHER_EPSILON) {
                        scoreY -= 1;
                    }
                }
            }
        }

        if (0.f > scoreY) {
            flipY = -1;
            for (sketcherMinimizerAtom* a : mol->_atoms) {
                a->coordinates.setY(-a->coordinates.y());
            }
        }
        if (0.f > scoreX) {
            flipX = -1;
            for (sketcherMinimizerAtom* a : mol->_atoms) {
                a->coordinates.setX(-a->coordinates.x());
            }
        }
    }
}

void sketcherMinimizer::addToVector(float weight, float angle,
                                    vector<pair<float, float>>& angles)
{
    angle = roundToTwoDecimalDigits(angle);
    while (angle <= 0) {
        angle += static_cast<float>(M_PI);
    }
    for (unsigned int i = 0; i < angles.size(); i++) {
        if (angles[i].second < angle - SKETCHER_EPSILON) {
            if (i == angles.size() - 1) {
                angles.emplace_back(weight, angle);
                break;
            }
        } else if (angles[i].second - angle < SKETCHER_EPSILON &&
                   angles[i].second - angle > -SKETCHER_EPSILON) {
            angles[i].first += weight;
            break;
        } else {
            angles.insert(angles.begin() + i,
                          pair<float, float>(weight, angle));
            break;
        }
    }
    if (angles.empty()) {
        angles.emplace_back(weight, angle);
    }
}

/* if a peptide chain is present rotate the molecule so it's horizontal */
void sketcherMinimizer::addBestRotationInfoForPeptides(
    vector<pair<float, float>>& angles,
    const std::vector<sketcherMinimizerAtom*>& atoms)
{
    auto chetoCs = m_minimizer.getChetoCs(atoms);
    auto aminoNs = m_minimizer.getAminoNs(atoms);
    auto alphaCs = m_minimizer.getAlphaCs(atoms, chetoCs, aminoNs);
    for (auto alphaC : alphaCs) {
        sketcherMinimizerAtom* aminoN = nullptr;
        sketcherMinimizerAtom* chetoC = nullptr;
        for (auto neighbor : alphaC->neighbors) {
            if (aminoNs.find(neighbor) != aminoNs.end()) {
                aminoN = neighbor;
            } else if (chetoCs.find(neighbor) != chetoCs.end()) {
                chetoC = neighbor;
            }
        }
        if (aminoN && chetoC) {
            auto direction = aminoN->coordinates - chetoC->coordinates;
            float weight = 1000.f;
            float angle = atan2(-direction.y(), direction.x());
            addToVector(weight, angle, angles);
        }
    }
}

void sketcherMinimizer::bestRotation()
{
    for (sketcherMinimizerMolecule* mol : m_molecules) {
        vector<pair<float, float>> angles;
        if (mol->hasFixedFragments || mol->hasConstrainedFragments) {
            continue;
        }
        addBestRotationInfoForPeptides(angles, mol->getAtoms());
        float angle = 0.f;
        float lastAngle;
        unsigned int i = 0, j = 0;
        float weight = 1.f;
        auto increment = static_cast<float>(M_PI / 6);
        for (sketcherMinimizerAtom* a : mol->_atoms) {
            if (!a->rings.empty()) {
                continue;
            }
            if (a->neighbors.size() > 1) {
                for (i = 0; i < a->neighbors.size() - 1; i++) {
                    for (j = i + 1; j < a->neighbors.size(); j++) {

                        weight = 6;

                        if (a->neighbors[i]->neighbors.size() != 1) {
                            weight += 2;
                        }
                        if (a->neighbors[j]->neighbors.size() != 1) {
                            weight += 2;
                        }
                        if (a->neighbors[j]->atomicNumber == 6) {
                            weight += 1;
                        }
                        if (a->neighbors[j]->atomicNumber == 6) {
                            weight += 1;
                        }
                        if (a->neighbors[i]->charge == 0) {
                            weight += 1;
                        }
                        if (a->neighbors[j]->charge == 0) {
                            weight += 1;
                        }

                        sketcherMinimizerPointF p =
                            a->neighbors[i]->coordinates -
                            a->neighbors[j]->coordinates;
                        angle = atan2(-p.y(), p.x());
                        addToVector(weight, angle, angles);
                    }
                }
            }
        }
        for (sketcherMinimizerBond* b : mol->_bonds) {
            sketcherMinimizerPointF p =
                b->endAtom->coordinates - b->startAtom->coordinates;
            weight = 1;
            angle = atan2(-p.y(), p.x());
            angle = roundToTwoDecimalDigits(angle);

            while (angle <= 0) {
                angle += static_cast<float>(M_PI);
            }
            lastAngle = angle;
            for (unsigned int i = 0; i < 6; i++) {
                if (i == 1 || i == 5) {
                    weight = 5.f;
                } else if (i == 0 || i == 3) {
                    weight = 1.5f;
                } else {
                    weight = 1.f;
                }
                if (b->bondOrder == 2 && i == 3 &&
                    (b->startAtom->neighbors.size() == 1 ||
                     b->endAtom->neighbors.size() == 1)) {
                    weight += 1.5;
                }

                if (b->startAtom->neighbors.size() == 1 &&
                    b->endAtom->neighbors.size() == 1 && i == 0) {
                    weight += 10;
                }
                addToVector(weight, lastAngle, angles);
                lastAngle += increment;
                if (lastAngle > M_PI) {
                    lastAngle -= static_cast<float>(M_PI);
                }
            }
        }

        for (sketcherMinimizerFragment* f : mol->_fragments) {

            vector<sketcherMinimizerRing*> rings = f->getRings();
            vector<sketcherMinimizerRing*> inPlaneRings = rings;

            size_t ringsN = inPlaneRings.size();
            if (ringsN == 2) {

                sketcherMinimizerRing* r1 = inPlaneRings[0];
                sketcherMinimizerRing* r2 = inPlaneRings[1];
                sketcherMinimizerPointF p = r2->findCenter() - r1->findCenter();
                p.normalize();
                angle = atan2(-p.y(), p.x());
                weight = 25.f;

                addToVector(weight, angle, angles);

            } else if (ringsN == 3) {

                sketcherMinimizerPointF r1 = inPlaneRings[0]->findCenter();
                sketcherMinimizerPointF r2 = inPlaneRings[1]->findCenter();

                for (sketcherMinimizerRing* r : inPlaneRings) {
                    vector<sketcherMinimizerRing*> fusedWith;
                    vector<vector<sketcherMinimizerAtom*>> fusionAtoms;

                    for (unsigned int fw = 0; fw < r->fusedWith.size(); fw++) {
                        fusedWith.push_back(r->fusedWith[fw]);
                        fusionAtoms.push_back(r->fusionAtoms[fw]);
                    }
                    if (fusedWith.size() == 2) {
                        if (fusionAtoms[0].size() == 2 &&
                            fusionAtoms[1].size() == 2) {
                            r1 = (fusionAtoms[0][0]->coordinates +
                                  fusionAtoms[0][1]->coordinates) *
                                 0.5;
                            r2 = (fusionAtoms[1][0]->coordinates +
                                  fusionAtoms[1][1]->coordinates) *
                                 0.5;
                            break;
                        }
                    }
                }

                sketcherMinimizerPointF p = r2 - r1;
                angle = atan2(-p.y(), p.x());
                weight = 50.f;
                addToVector(weight, angle, angles);

            } else {
                vector<sketcherMinimizerRing*> rings;
                for (sketcherMinimizerRing* r : inPlaneRings) {
                    if (r->_atoms.size() == 6) {
                        rings.push_back(r);
                    }
                }
                for (sketcherMinimizerRing* r : rings) {
                    for (auto fusionAts : r->fusionAtoms) {
                        if (fusionAts.size() == 2) {
                            sketcherMinimizerPointF p =
                                fusionAts[0]->coordinates -
                                fusionAts[1]->coordinates;
                            //  if (p.x () != p.x () || p.y () != p.y ()) p =
                            //  sketcherMinimizerPointF (50.f, 0.f);
                            sketcherMinimizerPointF rotatedP(p.y(), p.x());
                            angle = static_cast<float>(atan2(-p.y(), p.x()) -
                                                       M_PI * 0.5);
                            weight = 25.f;
                            addToVector(weight, angle, angles);
                        }
                    }
                }
            }
        }
        if (angles.size() > 1) {
            if (angles[angles.size() - 1].second - angles[0].second >=
                M_PI - 2 * SKETCHER_EPSILON) {
                angles[0].first += angles[angles.size() - 1].first;
                angles.erase(angles.begin() + angles.size() - 1);
            }
        }

        if (!angles.empty()) {
            int bestI = 0;
            for (i = 0; i < angles.size(); i++) {
                if (angles[i].first > angles[bestI].first) {
                    bestI = i;
                }
            }
            float s = -sin(angles[bestI].second);
            float c = cos(angles[bestI].second);
            sketcherMinimizerPointF center(0.f, 0.f);
            for (sketcherMinimizerAtom* at : mol->_atoms)
                center += at->coordinates;
            if (!mol->_atoms.empty()) {
                center /= mol->_atoms.size();
            }

            for (sketcherMinimizerAtom* at : mol->_atoms) {
                sketcherMinimizerPointF v = at->coordinates - center;
                v.rotate(s, c);
                at->setCoordinates(center + v);
            }
            sin_flip = s;
            cos_flip = c;
            centerX = center.x();
            centerY = center.y();
        }
    }
}

void sketcherMinimizer::writeMinimizationData()
{
    // print all coordinates to output file
    sketcherMinimizerPointF center(centerX, centerY);
    std::ofstream energy_file("minimization_data.txt");
    for (size_t i = 0; i < m_minimizer.energy_list.size(); ++i) {
        energy_file << m_minimizer.energy_list[i] << ";";
        for (auto coord : m_minimizer.all_coordinates[i]) {
            sketcherMinimizerPointF v = coord - center;
            v.rotate(sin_flip, cos_flip);
            sketcherMinimizerPointF new_coord = center + v;
            energy_file << new_coord.x() * flipX << "," << new_coord.y() * flipY
                        << ";";
        }
        energy_file << "\n";
    }
    energy_file.close();

    // print atom mapping to output file
    std::ofstream atom_mapping_file("atom_mapping.txt");
    for (size_t i = 0; i < m_referenceAtoms.size(); ++i) {
        size_t actual_idx = -1;
        for (size_t j = 0; j < m_atoms.size(); ++j) {
            if (m_referenceAtoms[i] == m_atoms[j]) {
                actual_idx = j;
                break;
            }
        }
        atom_mapping_file << i << "," << actual_idx << ";";
    }
    atom_mapping_file.close();
}

void sketcherMinimizer::findFragments()
{

    assert(!m_molecules.empty());
    for (sketcherMinimizerMolecule* mol : m_molecules) {
        CoordgenFragmenter::splitIntoFragments(mol);
        if (mol->_fragments.empty()) {
            continue;
        }
        vector<sketcherMinimizerFragment*> fragments = mol->_fragments;
        _fragments.reserve(_fragments.size() + fragments.size());
        _fragments.insert(_fragments.end(), fragments.begin(), fragments.end());
        m_independentFragments.push_back(mol->getMainFragment());
    }

    m_minimizer.m_fragments = _fragments;

    initializeFragments();
}

void sketcherMinimizer::placeResiduesProteinOnlyModeCircleStyle(
    const std::map<std::string, std::vector<sketcherMinimizerResidue*>>& chains)
{
    size_t totalResiduesNumber = m_residues.size() + chains.size();

    auto angle = static_cast<float>(2.f * M_PI / totalResiduesNumber);
    const float residueRadius = 30.f;
    const auto circumference =
        static_cast<float>(totalResiduesNumber * residueRadius * 2);
    const auto radius = static_cast<float>(circumference * 0.5 / M_PI);
    int i = 0;
    for (const auto& chain : chains) {
        ++i; // gap between chains
        auto residues = chain.second;
        sort(residues.begin(), residues.end(),
             [](const sketcherMinimizerResidue* firstRes,
                const sketcherMinimizerResidue* secondRes) {
                 int firstN = firstRes->resnum;
                 int secondN = secondRes->resnum;
                 return firstN < secondN;
             });

        for (auto res : residues) {
            sketcherMinimizerPointF p(radius, 0);
            // place residues in a circle
            p.rotate(sin(angle * i), cos(angle * i));
            res->coordinates = p;
            res->coordinatesSet = true;
            res->molecule->isPlaced = true;
            ++i;
        }
    }
}

std::map<std::string, sketcherMinimizerPointF>
sketcherMinimizer::computeChainsStartingPositionsMetaMol(
    const std::map<std::string, std::vector<sketcherMinimizerResidue*>>& chains)
{
    map<std::string, sketcherMinimizerAtom*> molMap;

    auto* metaMol = new sketcherMinimizerMolecule;
    for (const auto& chainPair : chains) {
        auto* a = new sketcherMinimizerAtom;
        a->molecule = metaMol;
        metaMol->_atoms.push_back(a);
        molMap[chainPair.first] = a;
    }

    for (const auto& chainPair : chains) {
        for (auto residue : chainPair.second) {
            for (auto interaction : residue->residueInteractions) {
                if (interaction->startAtom->isResidue() &&
                    interaction->endAtom->isResidue()) {
                    auto* r1 = static_cast<sketcherMinimizerResidue*>(
                        interaction->startAtom);
                    auto* r2 = static_cast<sketcherMinimizerResidue*>(
                        interaction->endAtom);
                    if (r1->chain != r2->chain) {
                        // add a bond to the metaMol if it doesn't exist already
                        sketcherMinimizerAtom* at1 = molMap[r1->chain];
                        sketcherMinimizerAtom* at2 = molMap[r2->chain];
                        bool found = false;
                        for (sketcherMinimizerBond* b : metaMol->_bonds) {
                            if ((b->startAtom == at1 && b->endAtom == at2) ||
                                (b->startAtom == at2 && b->endAtom == at1)) {
                                found = true;
                                break;
                            }
                        }
                        if (!found) {
                            auto* newBond = new sketcherMinimizerBond;
                            newBond->startAtom = at1;
                            newBond->endAtom = at2;
                            metaMol->_bonds.push_back(newBond);
                        }
                    }
                }
            }
        }
    }

    sketcherMinimizer min;
    if (!metaMol->_atoms.empty()) {
        min.setEvenAngles(true);
        min.initialize(metaMol);
        findFragments();
        min.buildFromFragments(true);
        min.m_minimizer.avoidClashes();
        min.bestRotation();
        min.maybeFlip();
        min.arrangeMultipleMolecules();
    }
    std::map<std::string, sketcherMinimizerPointF> positions;
    for (const auto& iter : molMap) {
        positions[iter.first] = iter.second->coordinates * 10.;
    }
    return positions;
}

void sketcherMinimizer::shortenInteractions(
    const std::map<std::string, std::vector<sketcherMinimizerResidue*>>& chains)
{
    for (const auto& chain : chains) {
        for (auto res : chain.second) {
            for (auto interaction : res->residueInteractions) {
                sketcherMinimizerPointF midPoint =
                    0.5 * (interaction->startAtom->coordinates +
                           interaction->endAtom->coordinates);
                res->coordinates += (midPoint - res->coordinates) * 0.1;
            }
        }
    }
}

std::vector<sketcherMinimizerResidue*> sketcherMinimizer::orderResiduesOfChains(
    const std::map<std::string, std::vector<sketcherMinimizerResidue*>>& chains)
{
    std::vector<sketcherMinimizerResidue*> vec;
    for (const auto& chain : chains) {
        for (auto res : chain.second) {
            vec.push_back(res);
        }
    }
    sort(vec.begin(), vec.end(),
         [](const sketcherMinimizerResidue* firstRes,
            const sketcherMinimizerResidue* secondRes) {
             return firstRes->residueInteractions.size() >
                    secondRes->residueInteractions.size();
         });
    std::set<sketcherMinimizerResidue*> visitedResidues;
    std::queue<sketcherMinimizerResidue*> residueQueue;
    std::vector<sketcherMinimizerResidue*> finalVec;
    for (auto residue : vec) {
        if (visitedResidues.find(residue) != visitedResidues.end()) {
            continue;
        }
        residueQueue.push(residue);
        visitedResidues.insert(residue);
        while (!residueQueue.empty()) {
            auto topResidue = residueQueue.front();
            finalVec.push_back(topResidue);
            residueQueue.pop();
            for (auto partner : topResidue->residueInteractionPartners) {
                auto* partnerRes =
                    static_cast<sketcherMinimizerResidue*>(partner);
                if (visitedResidues.find(partnerRes) == visitedResidues.end()) {
                    residueQueue.push(partnerRes);
                    visitedResidues.insert(partnerRes);
                }
            }
        }
    }
    return finalVec;
}

void sketcherMinimizer::placeResiduesProteinOnlyModeLIDStyle(
    const std::map<std::string, std::vector<sketcherMinimizerResidue*>>& chains)
{
    auto positions = computeChainsStartingPositionsMetaMol(chains);
    sketcherMinimizerPointF p;
    for (const auto& chain : chains) {
        p = positions[chain.first];
        for (auto res : chain.second) {
            res->coordinates = p;
        }
    }
    shortenInteractions(chains);
    auto residues = orderResiduesOfChains(chains);
    for (auto res : residues) {
        sketcherMinimizerResidue* firstPartner = nullptr;
        for (auto partner : res->residueInteractionPartners) {
            if (partner->coordinatesSet) {
                auto partnerResidue =
                    static_cast<sketcherMinimizerResidue*>(partner);
                if (!firstPartner && partnerResidue->chain != res->chain) {
                    firstPartner = partnerResidue;
                }
            }
        }

        /* when searching for a position for res prefer a direction
         * perpendicular to the direction of interactions to optimize use of
         * space */
        sketcherMinimizerPointF direction(0, 1);
        if (firstPartner) {
            sketcherMinimizerPointF chainCenterDirections =
                positions[firstPartner->chain] - positions[res->chain];
            direction = sketcherMinimizerPointF(-chainCenterDirections.y(),
                                                chainCenterDirections.x());
            direction.normalize();
            direction *= 4.;
        }
        res->coordinates = exploreGridAround(res->coordinates, 10, 5, 0, 0,
                                             -1.f, false, res, direction);
        res->coordinatesSet = true;
        res->molecule->isPlaced = true;
    }
}

void sketcherMinimizer::placeResiduesProteinOnlyMode()
{

    std::map<std::string, std::vector<sketcherMinimizerResidue*>> chains;
    for (auto residue : m_residues) {
        string chainOfResidue = residue->chain;
        chains[chainOfResidue].push_back(residue);
    }
    placeResiduesProteinOnlyModeLIDStyle(chains);
    m_minimizer.minimizeProteinOnlyLID(chains);
}

void sketcherMinimizer::placeResiduesInCrowns()
{
    auto SSEs = groupResiduesInSSEs(m_residues);
    /* sort secondary structure elements so that the most importants are placed
     * first. prefer longer SSEs and ones that make more interactions */
    sort(SSEs.begin(), SSEs.end(),
         [](const vector<sketcherMinimizerResidue*>& firstSSE,
            const vector<sketcherMinimizerResidue*>& secondSSE) {
             float interactionsOfFirst = 0, interactionsOfSecond = 0;
             for (auto res : firstSSE) {
                 interactionsOfFirst += res->residueInteractions.size();
             }
             for (auto res : secondSSE) {
                 interactionsOfSecond += res->residueInteractions.size();
             }
             float interactionScaling = 3.f;
             float score1 = firstSSE.size() + interactionScaling *
                                                  interactionsOfFirst /
                                                  firstSSE.size();
             float score2 = secondSSE.size() + interactionScaling *
                                                   interactionsOfSecond /
                                                   secondSSE.size();
             return score1 > score2;
         });
    bool needOtherShape = true;
    int shapeCounter = 0;

    // place residues in a crowns around the ligand. Keep expanding to further
    // away crowns until all residues are placed
    while (needOtherShape) {
        vector<sketcherMinimizerPointF> shape =
            shapeAroundLigand(shapeCounter++);
        needOtherShape = fillShape(SSEs, shape, shapeCounter);
    }
}

/* place residues in SSEs in the current shape. Return false if all residues are
 * place, true otherwise */
bool sketcherMinimizer::fillShape(
    vector<vector<sketcherMinimizerResidue*>>& SSEs,
    const vector<sketcherMinimizerPointF>& shape, int shapeN)
{
    vector<bool> penalties(shape.size(), false);
    std::set<sketcherMinimizerResidue*> outliers;
    for (const auto& SSE : SSEs) {
        placeSSE(SSE, shape, shapeN, penalties, outliers);
    }
    return !outliers.empty();
}

/* assign a penalty for the stretching of bonds between residues in the same
 * SSE */
float sketcherMinimizer::scoreSSEBondStretch(
    const sketcherMinimizerPointF& coordinates1,
    const sketcherMinimizerPointF& coordinates2)
{
    const float stretchPenalty = 400.f;
    auto squaredLength = (coordinates2 - coordinates1).squareLength();
    return squaredLength * stretchPenalty;
}

float sketcherMinimizer::getResidueDistance(
    float startF, float increment, sketcherMinimizerResidue* resToConsider,
    const vector<sketcherMinimizerResidue*>& SSE) const
{
    float totalF = startF;
    sketcherMinimizerResidue* lastRes = nullptr;
    for (auto res : SSE) {
        if (lastRes) {
            auto result = static_cast<float>(res->resnum - lastRes->resnum);
            /* if the gap is more than 1, make the distance a bit smaller for
             * aesthetic reasons */
            result = static_cast<float>(1 + (result - 1) * 0.8);
            if (result < 1.f) {
                result = 1.f;
            }
            totalF += increment * result;
        }
        if (res == resToConsider) {
            break;
        }
        lastRes = res;
    }
    return totalF;
}

/* return a score for the placing on the SSE starting at startingPosition and
 * separated by increment */
float sketcherMinimizer::scoreSSEPosition(
    const vector<sketcherMinimizerResidue*>& SSE,
    const vector<sketcherMinimizerPointF>& shape, int shapeN,
    vector<bool>& penalties, float startingPosition, float increment)
{
    float score = 0.f;
    sketcherMinimizerResidue* lastResidue = nullptr;
    int lastResiduePosition = 0;
    sketcherMinimizerPointF lastResidueCoordinates;
    for (auto res : SSE) {
        int index = getShapeIndex(
            shape, getResidueDistance(startingPosition, increment, res, SSE));
        auto residueCoordinates = shape.at(index);
        int residuePosition = 0;
        if (res->coordinatesSet) {
            residuePosition = -1;
            residueCoordinates = res->coordinates;
        } else {
            if (!penalties[index]) {
                residuePosition = 0;
            } else {
                residuePosition = 1;
            }
        }
        if (residuePosition != -1) {
            score += scoreResiduePosition(index, shape, shapeN, penalties, res);
        }

        if (lastResidue && (residuePosition != lastResiduePosition))

        {

            score +=
                scoreSSEBondStretch(residueCoordinates, lastResidueCoordinates);
        }
        lastResiduePosition = residuePosition;
        lastResidueCoordinates = residueCoordinates;
        lastResidue = res;
    }
    return score;
}

void sketcherMinimizer::placeSSE(const vector<sketcherMinimizerResidue*>& SSE,
                                 const vector<sketcherMinimizerPointF>& shape,
                                 int shapeN, vector<bool>& penalties,
                                 set<sketcherMinimizerResidue*>& outliers,
                                 bool placeOnlyInteracting)
{
    int residuesToPlace = 0;
    for (auto res : SSE) {
        if (!res->coordinatesSet) {
            residuesToPlace++;
        }
    }
    if (residuesToPlace == 0) {
        return;
    }
    typedef pair<float, float> Solution;
    vector<pair<float, Solution>> scoredSolutions;
    /* move around the crown scoring possible solutions, varying the starting
     * position and the separation between consecutive residues */
    for (float f = 0.f; f < 1.f; f += 0.004f) {
        float distance = 5.f / shape.size();
        for (float increment = -1 * distance; increment <= 1 * distance;
             increment += distance) {
            if (increment == 0) {
                continue;
            }
            float score =
                scoreSSEPosition(SSE, shape, shapeN, penalties, f, increment);
            scoredSolutions.emplace_back(score, Solution(f, increment));
        }
    }
    auto bestResult =
        min_element(scoredSolutions.begin(), scoredSolutions.end());
    set<sketcherMinimizerResidue*> alreadyPlaced;
    for (auto residue : SSE) {
        if (residue->coordinatesSet) {
            continue; // placed in a previous crown
        }
        float f = getResidueDistance(bestResult->second.first,
                                     bestResult->second.second, residue, SSE);
        int index = getShapeIndex(shape, f);
        bool alreadyAResidueHere = penalties.at(index);
        sketcherMinimizerPointF position = shape.at(index);
        if (alreadyAResidueHere ||
            (placeOnlyInteracting &&
             residue->residueInteractionPartners.empty())) {
            outliers.insert(residue);
        } else {
            residue->coordinates = position;
            alreadyPlaced.insert(residue);
        }
    }
    // mark the current solution to prevent other residues from being placed on
    // top of these
    markSolution(bestResult->second, SSE, shape, penalties, outliers);
    for (auto res : alreadyPlaced) {
        res->coordinatesSet = true;
    }
    for (auto res : SSE) {
        if (res->m_isWaterMap && res->m_isClashing && res->coordinatesSet &&
            res->m_closestLigandAtom != nullptr) {
            sketcherMinimizerPointF directionToLigand =
                res->m_closestLigandAtom->coordinates - res->coordinates;
            directionToLigand.normalize();
            float displacement = BONDLENGTH * 0.3;
            res->coordinates = res->m_closestLigandAtom->coordinates -
                               directionToLigand * displacement;
        }
    }
}

void sketcherMinimizer::markSolution(
    const pair<float, float>& solution,
    const vector<sketcherMinimizerResidue*>& SSE,
    const vector<sketcherMinimizerPointF>& shape, vector<bool>& penalties,
    set<sketcherMinimizerResidue*>& outliers)
{
    float padding = abs(solution.second) * 0.5f;
    sketcherMinimizerResidue* lastRes = nullptr;
    float lastF = 0.f;
    for (auto res : SSE) {
        if (res->coordinatesSet || (res->m_isWaterMap && res->m_isClashing) ||
            outliers.find(res) != outliers.end()) {
            lastRes = nullptr;
            lastF = 0.f;
            continue;
        }
        float f = getResidueDistance(solution.first, solution.second, res, SSE);
        int startIndex = getShapeIndex(shape, f - padding);
        int endIndex = getShapeIndex(shape, f + padding);
        for (int index = startIndex; index != endIndex;
             index = (index + 1) % shape.size()) {
            penalties.at(index) = true;
        }
        if (lastRes) {
            if (solution.second < 0) {
                std::swap(lastF, f);
            }
            int startIndex = getShapeIndex(shape, lastF);
            int endIndex = getShapeIndex(shape, f);
            for (int index = startIndex; index != endIndex;
                 index = (index + 1) % shape.size()) {
                penalties.at(index) = true;
            }
        }
        lastRes = res;
        lastF = f;
    }
}

int sketcherMinimizer::getShapeIndex(
    const vector<sketcherMinimizerPointF>& shape, float floatPosition) const
{
    float normalizedF = floatPosition;
    while (normalizedF < 0) {
        normalizedF += 1.f;
    }
    while (normalizedF >= 1.f) {
        normalizedF -= 1.f;
    }
    int counter = static_cast<int>(shape.size() * normalizedF);
    return counter;
}

vector<vector<sketcherMinimizerResidue*>>
sketcherMinimizer::groupResiduesInSSEs(
    const vector<sketcherMinimizerResidue*>& residues)
{
    // divide residues by chain
    map<string, vector<sketcherMinimizerResidue*>> chainsMap;
    for (auto res : residues) {
        chainsMap[res->chain].push_back(res);
    }
    // order each chain by residue number
    for (auto& pair : chainsMap) {
        sort(pair.second.begin(), pair.second.end(),
             [](const sketcherMinimizerResidue* firstRes,
                const sketcherMinimizerResidue* secondRes) {
                 return firstRes->resnum < secondRes->resnum;
             });
    }
    int gap = 3;
    // split chains in smaller chunks whenever more than gap consecutive
    // residues are missing
    vector<vector<sketcherMinimizerResidue*>> returnValue;
    for (auto& pair : chainsMap) {
        vector<sketcherMinimizerResidue*> growingChain;
        for (auto res : pair.second) {
            if (!growingChain.empty() &&
                (res->resnum - growingChain.back()->resnum > gap ||
                 res->chain == " " || res->chain.empty())) {
                returnValue.push_back(growingChain);
                growingChain.clear();
            }
            growingChain.push_back(res);
        }
        if (!growingChain.empty()) {
            returnValue.push_back(growingChain);
        }
    }
    return returnValue;
}

vector<sketcherMinimizerPointF> sketcherMinimizer::shapeAroundLigand(int crownN)
{
    // return crownN-th crown around ligand
    float distanceOfFirstCrown = 60;
    float distanceBetweenCrowns = 60;
    // find limits
    const auto& atoms = m_atoms;
    const auto& bonds = m_bonds;
    float border = distanceBetweenCrowns * crownN + distanceOfFirstCrown;
    float minX = atoms[0]->coordinates.x();
    float maxX = atoms[0]->coordinates.x();
    float minY = atoms[0]->coordinates.y();
    float maxY = atoms[0]->coordinates.y();
    float maxPocketD = 0;
    for (auto atom : atoms) {
        float distance = atom->m_pocketDistance;
        if (distance > maxPocketD) {
            maxPocketD = distance;
        }
        float newX = atom->coordinates.x();
        float newY = atom->coordinates.y();

        if (minX > newX) {
            minX = newX;
        }
        if (maxX < newX) {
            maxX = newX;
        }
        if (minY > newY) {
            minY = newY;
        }
        if (maxY < newY) {
            maxY = newY;
        }
    }
    maxPocketD += 10; // to account for cutoffs at borders
    minX -= border + maxPocketD;
    maxX += border + maxPocketD;
    minY -= border + maxPocketD;
    maxY += border + maxPocketD;

    /* run a marching square algorithm on a grid of x_interval spacing */
    sketcherMinimizerMarchingSquares ms;
    float x_interval = 20.f;
    ms.initialize(minX, maxX, minY, maxY, x_interval);
    for (unsigned int j = 0; j < ms.getYN(); j++) {
        for (unsigned int i = 0; i < ms.getXN(); i++) {
            float pointX = ms.toRealx(static_cast<float>(i));
            float pointY = ms.toRealy(static_cast<float>(j));
            sketcherMinimizerPointF p(pointX, pointY);

            float shortestD = -1;
            for (auto a : atoms) {
                float dist = a->m_pocketDistance;

                auto vect = a->coordinates - p;
                float D = vect.length();
                D -= dist + border;
                if (D < shortestD || shortestD < 0) {
                    shortestD = D;
                }
            }

            for (auto b : bonds) {

                sketcherMinimizerPointF sp1(b->startAtom->coordinates);
                sketcherMinimizerPointF sp2(b->endAtom->coordinates);

                float distancePercentage = 1.f;
                float D2 = sketcherMinimizerMaths::squaredDistancePointSegment(
                    p, sp1, sp2, &distancePercentage);
                float D = sqrt(D2);
                float distance2 = b->startAtom->m_pocketDistance;
                float distance1 = b->endAtom->m_pocketDistance;
                float dist = distance1 * distancePercentage +
                             distance2 * (1 - distancePercentage);
                D -= dist + border;

                if (D < shortestD) {
                    shortestD = D;
                }
            }

            ms.setValue(shortestD, i, j);
        }
    }
    ms.setThreshold(0);
    ms.run();
    auto result = ms.getOrderedCoordinatesPoints();
    sort(result.begin(), result.end(),
         [](const vector<float>& firstContour,
            const vector<float>& secondContour) {
             return firstContour.size() > secondContour.size();
         });
    vector<sketcherMinimizerPointF> returnValue;
    if (!result.empty()) {
        for (unsigned int i = 0; i < result.at(0).size(); i += 2) {
            returnValue.emplace_back(result.at(0).at(i),
                                     result.at(0).at(i + 1));
        }
    }
    return returnValue;
}

float sketcherMinimizer::scoreResiduePosition(
    int index, const vector<sketcherMinimizerPointF>& shape, int shapeN,
    vector<bool>&, sketcherMinimizerResidue* residue)
{
    auto position = shape.at(index);
    float distancePenalty = 0.01f;
    float clashingLigandAtomsPenalty = 100.f;
    vector<sketcherMinimizerAtom*> targets;
    for (auto interactionPartner : residue->residueInteractionPartners) {
        if (interactionPartner->coordinatesSet) {
            targets.push_back(interactionPartner);
        }
    }
    float interactionsF = 1.f;
    if (targets.empty() && residue->m_closestLigandAtom != nullptr) {
        interactionsF = 0.2f;
        targets.push_back(residue->m_closestLigandAtom);
    }
    float score = 0.f;
    for (auto target : targets) {
        int clashingLigandAtoms = 0;
        for (auto ligandAtom : m_atoms) {
            if (ligandAtom == target) {
                continue;
            }
            auto ligandAtomPos = ligandAtom->coordinates;
            float squareDist =
                sketcherMinimizerMaths::squaredDistancePointSegment(
                    ligandAtomPos, position, target->coordinates);
            if (squareDist < 40 * 40) {
                clashingLigandAtoms++;
            }
        }
        auto distance = sketcherMinimizerMaths::squaredDistance(
                            target->coordinates, position) -
                        (shapeN * 50) * (shapeN * 50);
        score +=
            interactionsF * (distancePenalty * distance +
                             clashingLigandAtoms * clashingLigandAtomsPenalty);
    }
    return score;
}

void sketcherMinimizer::placeResidues(
    const vector<sketcherMinimizerAtom*>& atoms)
{
    if (m_residues.empty()) {
        return;
    }
    if (atoms.empty()) {
        placeResiduesProteinOnlyMode();
        return;
    }
    findClosestAtomToResidues(atoms);

    placeResiduesInCrowns();
    m_minimizer.minimizeResidues();
}

/*
 move mol around to avoid clashes with other already placed molecules. Explore a
 grid of @levels concentric levels, with #gridD resolution. @distanceFromAtoms
 is the minimum clash distance to reject a position.
 */

sketcherMinimizerPointF
sketcherMinimizer::exploreMolPosition(sketcherMinimizerMolecule* mol,
                                      unsigned int levels, float gridD,
                                      float distanceFromAtoms)
{

    sketcherMinimizerPointF v(0, 0), centerOfGrid(0, 0);
    for (unsigned int i = 0; i < levels; i++) {

        vector<sketcherMinimizerPointF> pointstoTest;
        sketcherMinimizerPointF top =
            centerOfGrid + sketcherMinimizerPointF(0.f, (1 + i) * gridD);
        sketcherMinimizerPointF bottom =
            centerOfGrid + sketcherMinimizerPointF(0.f, -((1 + i) * gridD));
        sketcherMinimizerPointF right =
            centerOfGrid + sketcherMinimizerPointF((1 + i) * gridD, 0.f);
        sketcherMinimizerPointF left =
            centerOfGrid + sketcherMinimizerPointF(-((1 + i) * gridD), 0.f);

        pointstoTest.push_back(centerOfGrid);
        pointstoTest.push_back(right);
        pointstoTest.push_back(left);
        pointstoTest.push_back(bottom);
        pointstoTest.push_back(top);

        for (unsigned int j = 0; j < i; j++) {
            pointstoTest.push_back(
                right + sketcherMinimizerPointF(0.f, gridD * (j + 1)));
            pointstoTest.push_back(
                right - sketcherMinimizerPointF(0.f, gridD * (j + 1)));
            pointstoTest.push_back(
                left + sketcherMinimizerPointF(0.f, gridD * (j + 1)));
            pointstoTest.push_back(
                left - sketcherMinimizerPointF(0.f, gridD * (j + 1)));
            pointstoTest.push_back(
                bottom + sketcherMinimizerPointF(gridD * (j + 1), 0.f));
            pointstoTest.push_back(
                bottom - sketcherMinimizerPointF(gridD * (j + 1), 0.f));
            pointstoTest.push_back(
                top + sketcherMinimizerPointF(gridD * (j + 1), 0.f));
            pointstoTest.push_back(
                top - sketcherMinimizerPointF(gridD * (j + 1), 0.f));
        }

        pointstoTest.push_back(
            centerOfGrid +
            sketcherMinimizerPointF((1 + i) * gridD, (1 + i) * gridD));
        pointstoTest.push_back(
            centerOfGrid +
            sketcherMinimizerPointF((1 + i) * gridD, -((1 + i) * gridD)));
        pointstoTest.push_back(
            centerOfGrid +
            sketcherMinimizerPointF(-((1 + i) * gridD), (1 + i) * gridD));
        pointstoTest.push_back(
            centerOfGrid +
            sketcherMinimizerPointF(-((1 + i) * gridD), -((1 + i) * gridD)));

        bool noClash = true;
        if (distanceFromAtoms < 0) {
            distanceFromAtoms = bondLength * 1.8;
        }
        float dist = distanceFromAtoms;

        for (const auto& pc : pointstoTest) {
            noClash = true;
            v = pc;
            for (sketcherMinimizerAtom* at : mol->_atoms) {
                sketcherMinimizerPointF placeNextTo = at->coordinates + v;
                for (sketcherMinimizerMolecule* m : m_molecules) {
                    if (!m->isPlaced) {
                        continue;
                    }
                    if (m == mol) {
                        continue;
                    }
                    for (sketcherMinimizerAtom* a : m->_atoms) {
                        dist = distanceFromAtoms;

                        if (((a->coordinates.x() < placeNextTo.x() + dist) &&
                             (a->coordinates.x() > placeNextTo.x() - dist)) &&
                            ((a->coordinates.y() < placeNextTo.y() + dist) &&
                             (a->coordinates.y() > placeNextTo.y() - dist))) {
                            noClash = false;
                            break;
                        }
                    }

                    if (!noClash) {
                        break;
                    }
                }

                if (!noClash) {
                    break;
                }
            }
            if (noClash) {
                break;
            }
        }
        if (noClash) {
            break;
        }
    }
    return v;
}

sketcherMinimizerPointF sketcherMinimizer::exploreGridAround(
    const sketcherMinimizerPointF& centerOfGrid, unsigned int levels,
    float gridD, float dx, float dy, float distanceFromAtoms, bool watermap,
    sketcherMinimizerResidue* residueForInteractions,
    const sketcherMinimizerPointF& direction)
{
    sketcherMinimizerPointF placeNextTo = centerOfGrid;

    for (unsigned int i = 0; i < levels; i++) {

        vector<sketcherMinimizerPointF> pointstoTest;
        sketcherMinimizerPointF top =
            centerOfGrid + sketcherMinimizerPointF(0.f, (1 + i) * gridD);
        sketcherMinimizerPointF bottom =
            centerOfGrid + sketcherMinimizerPointF(0.f, -((1 + i) * gridD));
        sketcherMinimizerPointF right =
            centerOfGrid + sketcherMinimizerPointF((1 + i) * gridD, 0.f);
        sketcherMinimizerPointF left =
            centerOfGrid + sketcherMinimizerPointF(-((1 + i) * gridD), 0.f);

        pointstoTest.push_back(centerOfGrid);
        pointstoTest.push_back(right);
        pointstoTest.push_back(left);
        pointstoTest.push_back(bottom);
        pointstoTest.push_back(top);

        for (unsigned int j = 0; j < i; j++) {
            pointstoTest.push_back(
                right + sketcherMinimizerPointF(0.f, gridD * (j + 1)));
            pointstoTest.push_back(
                right - sketcherMinimizerPointF(0.f, gridD * (j + 1)));
            pointstoTest.push_back(
                left + sketcherMinimizerPointF(0.f, gridD * (j + 1)));
            pointstoTest.push_back(
                left - sketcherMinimizerPointF(0.f, gridD * (j + 1)));
            pointstoTest.push_back(
                bottom + sketcherMinimizerPointF(gridD * (j + 1), 0.f));
            pointstoTest.push_back(
                bottom - sketcherMinimizerPointF(gridD * (j + 1), 0.f));
            pointstoTest.push_back(
                top + sketcherMinimizerPointF(gridD * (j + 1), 0.f));
            pointstoTest.push_back(
                top - sketcherMinimizerPointF(gridD * (j + 1), 0.f));
        }

        pointstoTest.push_back(
            centerOfGrid +
            sketcherMinimizerPointF((1 + i) * gridD, (1 + i) * gridD));
        pointstoTest.push_back(
            centerOfGrid +
            sketcherMinimizerPointF((1 + i) * gridD, -((1 + i) * gridD)));
        pointstoTest.push_back(
            centerOfGrid +
            sketcherMinimizerPointF(-((1 + i) * gridD), (1 + i) * gridD));
        pointstoTest.push_back(
            centerOfGrid +
            sketcherMinimizerPointF(-((1 + i) * gridD), -((1 + i) * gridD)));

        bool noClash = true;
        if (distanceFromAtoms < 0) {
            distanceFromAtoms = bondLength * 1.8;
        }
        float distanceFromResidues = bondLength * 1.3;
        float watermapDistance = 10;
        float dist = distanceFromAtoms;

        sketcherMinimizerPointF directionNormal(-direction.y(), direction.x());
        directionNormal.normalize();
        for (const auto& pc : pointstoTest) {
            noClash = true;
            sketcherMinimizerPointF point = pc - centerOfGrid;
            placeNextTo = point.y() * direction + point.x() * directionNormal +
                          centerOfGrid;
            for (sketcherMinimizerMolecule* m : m_molecules) {
                if (!m->isPlaced) {
                    continue;
                }
                for (sketcherMinimizerAtom* a : m->_atoms) {
                    if (a->isResidue()) {
                        dist = distanceFromResidues;
                    } else {
                        dist = distanceFromAtoms;
                    }
                    if (watermap) {
                        if (!a->isResidue()) {
                            continue;
                        }
                        dist = watermapDistance;
                    }
                    if (((a->coordinates.x() < placeNextTo.x() + dist + dx) &&
                         (a->coordinates.x() > placeNextTo.x() - dist - dx)) &&
                        ((a->coordinates.y() < placeNextTo.y() + dist + dy) &&
                         (a->coordinates.y() > placeNextTo.y() - dist - dy))) {
                        noClash = false;
                        break;
                    }
                    if (residueForInteractions) {
                        for (auto partnerOfA : a->residueInteractionPartners) {
                            if (a == residueForInteractions ||
                                partnerOfA == residueForInteractions ||
                                !a->coordinatesSet ||
                                !partnerOfA->coordinatesSet) {
                                continue;
                            }

                            float squareD = sketcherMinimizerMaths::
                                squaredDistancePointSegment(
                                    placeNextTo, a->coordinates,
                                    partnerOfA->coordinates);
                            if (squareD <
                                distanceFromResidues * distanceFromResidues) {
                                noClash = false;
                                break;
                            }

                            for (auto partner :
                                 residueForInteractions
                                     ->residueInteractionPartners) {
                                if (!partner->coordinatesSet) {
                                    continue;
                                }
                                if (sketcherMinimizerMaths::
                                        intersectionOfSegments(
                                            placeNextTo, partner->coordinates,
                                            a->coordinates,
                                            partnerOfA->coordinates)) {
                                    noClash = false;
                                    break;
                                }
                            }
                        }
                    }
                }

                if (!noClash) {
                    break;
                }
            }

            if (noClash) {
                break;
            }
        }
        if (noClash) {
            break;
        }
    }
    return placeNextTo;
}

vector<proximityData> sketcherMinimizer::buildProximityDataVector(
    vector<sketcherMinimizerMolecule*>& proximityMols,
    map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*>& molMap)
{

    vector<proximityData> proximityDataVector;

    for (sketcherMinimizerMolecule* mol : proximityMols) {
        proximityData data;
        sketcherMinimizerAtom* metaAtom = molMap[mol];
        vector<sketcherMinimizerPointF> additionVectors(
            metaAtom->neighbors.size(), sketcherMinimizerPointF(0.f, 0.f));
        vector<sketcherMinimizerPointF> centers(
            metaAtom->neighbors.size(), sketcherMinimizerPointF(0.f, 0.f));
        vector<int> counters(metaAtom->neighbors.size(), 0);
        for (sketcherMinimizerBond* pr : mol->m_proximityRelations) {
            sketcherMinimizerAtom* otherMetaAtom = nullptr;
            sketcherMinimizerAtom* targetAtom = nullptr;
            if (pr->startAtom->molecule == mol &&
                !(pr->endAtom->molecule == mol)) {
                otherMetaAtom = molMap[pr->endAtom->molecule];
                targetAtom = pr->startAtom;
            } else if (pr->endAtom->molecule == mol &&
                       !(pr->startAtom->molecule == mol)) {
                otherMetaAtom = molMap[pr->startAtom->molecule];
                targetAtom = pr->endAtom;
            }
            if (otherMetaAtom) {
                for (unsigned int i = 0; i < metaAtom->neighbors.size(); i++) {
                    if (metaAtom->neighbors[i] == otherMetaAtom) {
                        additionVectors[i] +=
                            targetAtom->getSingleAdditionVector();
                        centers[i] += targetAtom->coordinates;
                        counters[i]++;
                    }
                }
            }
        }
        for (unsigned int i = 0; i < centers.size(); i++) {
            if (counters[i] > 0) {
                centers[i] /= counters[i];
            }
            additionVectors[i].normalize();
        }
        data.additionVectors = additionVectors;
        data.centers = centers;
        data.counters = counters;
        proximityDataVector.push_back(data);
    }
    return proximityDataVector;
}

void sketcherMinimizer::rotateMoleculesWithProximityRelations(
    vector<sketcherMinimizerMolecule*>& proximityMols,
    map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*>& molMap,
    vector<proximityData>& proximityDataVector)
{

    for (unsigned int m = 0; m < proximityMols.size(); m++) {
        sketcherMinimizerMolecule* mol = proximityMols[m];
        sketcherMinimizerAtom* metaAtom = molMap[mol];
        vector<sketcherMinimizerPointF> additionVectors =
            proximityDataVector[m].additionVectors;
        vector<sketcherMinimizerPointF> centers =
            proximityDataVector[m].centers;
        if (mol->_atoms.size() < 2) {
            continue;
        }

        sketcherMinimizerPointF direction(1, 0);
        if (metaAtom->bonds.size() == 1) {
            direction =
                metaAtom->coordinates - metaAtom->neighbors[0]->coordinates;
            sketcherMinimizerPointF p1 = additionVectors[0];
            p1 *= -1;
            sketcherMinimizerPointF p3 = direction;
            float rotationAngle = sketcherMinimizerMaths::signedAngle(
                p1, sketcherMinimizerPointF(0, 0), p3);
            rotationAngle *= static_cast<float>(-M_PI / 180.f);
            float s = sin(rotationAngle);
            float c = cos(rotationAngle);

            for (sketcherMinimizerAtom* a : mol->_atoms) {
                sketcherMinimizerPointF coords = a->coordinates - centers[0];
                coords.rotate(s, c);
                a->coordinates = coords + centers[0];
            }
        } else if (metaAtom->bonds.size() > 1) {
            vector<sketcherMinimizerPointF> v1, v2;
            for (sketcherMinimizerAtom* n : metaAtom->neighbors) {
                v1.push_back(n->coordinates - metaAtom->coordinates);
            }
            v2 = additionVectors;
            float rotMat[4];
            alignmentMatrix(v1, v2, rotMat);
            sketcherMinimizerPointF center = mol->center();
            for (sketcherMinimizerAtom* a : mol->_atoms) {
                sketcherMinimizerPointF coords = a->coordinates - center;
                float x = coords.x();
                float y = coords.y();
                sketcherMinimizerPointF newCoords(x * rotMat[0] + y * rotMat[1],
                                                  x * rotMat[2] +
                                                      y * rotMat[3]);
                a->coordinates = center + newCoords;
            }
        }
    }
}

void sketcherMinimizer::translateMoleculesWithProximityRelations(
    vector<sketcherMinimizerMolecule*>& proximityMols,
    map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*>& molMap,
    map<sketcherMinimizerMolecule*, sketcherMinimizerPointF>& templateCenters,
    vector<proximityData>&)
{

    // placing
    int counterN = 1;
    bool cleverPlacing = false;
    do {
        cleverPlacing = !cleverPlacing; // alternatively try to be smart aboout
                                        // single atom mols
        if (!cleverPlacing) {
            counterN++;
        }

        for (auto mol : proximityMols) {
            bool residue = false;
            if (mol->_atoms.size() == 1) {
                if (mol->_atoms[0]->isResidue()) {
                    residue = true;
                }
            }
            if (!residue) {
                if (mol->hasConstrainedFragments) {
                    mol->isPlaced = true;
                    continue;
                }
            }
            if (mol->hasFixedFragments) {
                mol->isPlaced = true;
                continue;
            }
            if (!mol->m_proximityRelations.empty()) {
                sketcherMinimizerPointF atomsCenter =
                    sketcherMinimizerPointF(0, 0);
                int atomsN = 0;
                for (sketcherMinimizerBond* pr : mol->m_proximityRelations) {
                    if (pr->startAtom->molecule == mol &&
                        pr->endAtom->molecule != mol) {
                        atomsCenter += pr->startAtom->coordinates;
                        atomsN++;
                    } else if (pr->endAtom->molecule == mol &&
                               pr->startAtom->molecule != mol) {
                        atomsCenter += pr->endAtom->coordinates;
                        atomsN++;
                    }
                }
                if (atomsN > 0) {
                    atomsCenter /= atomsN;
                }

                /* positioning */
                sketcherMinimizerPointF placeNextTo = templateCenters[mol];
                placeNextTo *= counterN;

                for (sketcherMinimizerAtom* a : mol->_atoms) {
                    a->coordinates += placeNextTo - atomsCenter;
                }

                mol->isPlaced = true;
            }
        }
        if (cleverPlacing) { // replace single terminal atoms
            for (auto mol : proximityMols) {
                if (mol->_atoms.size() == 1) {
                    sketcherMinimizerAtom* metaAtom = molMap[mol];
                    if (metaAtom->neighbors.size() == 1) {

                        int bondsN = 0;
                        sketcherMinimizerPointF coords(0, 0);
                        for (sketcherMinimizerBond* pr :
                             mol->m_proximityRelations) {
                            if (pr->startAtom->molecule == mol &&
                                pr->endAtom->molecule != mol) {
                                sketcherMinimizerPointF addV =
                                    pr->endAtom->getSingleAdditionVector();
                                if (addV.length() < SKETCHER_EPSILON) {
                                    continue;
                                }
                                addV.normalize();
                                addV *= bondLength * counterN;
                                coords += pr->endAtom->coordinates + addV;
                                bondsN++;
                            } else if (pr->endAtom->molecule == mol &&
                                       pr->startAtom->molecule != mol) {
                                sketcherMinimizerPointF addV =
                                    pr->startAtom->getSingleAdditionVector();
                                if (addV.length() < SKETCHER_EPSILON) {
                                    continue;
                                }

                                addV.normalize();
                                addV *= bondLength * counterN;
                                coords += pr->startAtom->coordinates + addV;
                                bondsN++;
                            }
                        }
                        if (bondsN > 0) {
                            coords /= bondsN;

                            mol->_atoms[0]->coordinates = coords;
                        } else { // a suitable addition Vector could not be
                                 // found, try positionings a bondlength to the
                                 // right of a proximity partner
                            for (sketcherMinimizerBond* pr :
                                 mol->m_proximityRelations) {
                                if (pr->startAtom->molecule == mol &&
                                    pr->endAtom->molecule != mol) {
                                    mol->_atoms[0]->coordinates =
                                        pr->endAtom->coordinates;
                                    break;
                                } else if (pr->endAtom->molecule == mol &&
                                           pr->startAtom->molecule != mol) {
                                    mol->_atoms[0]->coordinates =
                                        pr->startAtom->coordinates;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

    } while (m_minimizer.findIntermolecularClashes(proximityMols,
                                                   bondLength * 0.5) &&
             counterN < 10);
}

void sketcherMinimizer::placeMoleculesWithProximityRelations(
    vector<sketcherMinimizerMolecule*> proximityMols)
{
    map<sketcherMinimizerMolecule*, sketcherMinimizerAtom*> molMap;

    auto* metaMol = new sketcherMinimizerMolecule;
    for (sketcherMinimizerMolecule* mol : m_molecules) {
        if (!mol->m_proximityRelations.empty()) {
            auto* a = new sketcherMinimizerAtom;
            a->molecule = metaMol;
            metaMol->_atoms.push_back(a);
            molMap[mol] = a;
        }
    }

    for (sketcherMinimizerBond* b : m_proximityRelations) {
        if (b->startAtom->molecule == b->endAtom->molecule) {
            continue;
        }
        sketcherMinimizerAtom* at1 = molMap[b->startAtom->molecule];
        sketcherMinimizerAtom* at2 = molMap[b->endAtom->molecule];
        bool found = false;
        for (sketcherMinimizerBond* b : metaMol->_bonds) {
            if ((b->startAtom == at1 && b->endAtom == at2) ||
                (b->startAtom == at2 && b->endAtom == at1)) {
                found = true;
            }
        }
        if (!found) {
            auto* newBond = new sketcherMinimizerBond;
            newBond->startAtom = at1;
            newBond->endAtom = at2;

            metaMol->_bonds.push_back(newBond);
        }
    }
    sketcherMinimizer min(m_minimizer.getPrecision());

    if (!metaMol->_atoms.empty()) {
        min.setEvenAngles(true);

        min.initialize(metaMol);
        min.findFragments();
        min.buildFromFragments(true);
        min.m_minimizer.avoidClashes();
        min.bestRotation();
        min.maybeFlip();
        min.arrangeMultipleMolecules();
    }

    bool ligandResidueStyle = true; // positions of molecules are determined
                                    // more by a bigger central molecule than by
                                    // a bonding pattern

    for (auto molecule : min.m_molecules) {
        if (!molecule->_rings.empty()) {
            // if at least three molecules are connected to each other
            // (i.e. two residues are connected to each other and both to
            // the ligand) abort the ligandResidue display style)
            ligandResidueStyle = false;
        }
    }
    sketcherMinimizerMolecule* centralMol = proximityMols[0];
    for (sketcherMinimizerMolecule* mol : proximityMols) {
        if (mol->m_proximityRelations.size() >
            centralMol->m_proximityRelations.size()) {
            centralMol = mol;
        } else if (mol->m_proximityRelations.size() ==
                       centralMol->m_proximityRelations.size() &&
                   mol->_atoms.size() > centralMol->_atoms.size()) {
            centralMol = mol;
        }
    }
    if (centralMol->_atoms.size() < MINIMUM_LIGAND_ATOMS) {
        ligandResidueStyle = false;
    }
    map<sketcherMinimizerMolecule*, sketcherMinimizerPointF> templateCenters;

    for (sketcherMinimizerMolecule* mol : proximityMols) {
        sketcherMinimizerPointF point(0, 0);

        sketcherMinimizerAtom* at = molMap[mol];
        if (at) {
            point = at->coordinates;
        }
        templateCenters[mol] = point;
    }

    if (ligandResidueStyle) {
        queue<sketcherMinimizerMolecule*> q;
        map<sketcherMinimizerMolecule*, sketcherMinimizerMolecule*> getParent;

        q.push(centralMol);
        while (!q.empty()) {
            sketcherMinimizerMolecule* mol = q.front();
            q.pop();
            if (mol->isPlaced) {
                continue;
            }
            if (mol == centralMol) {
                mol->isPlaced = true;
            } else {
                sketcherMinimizerMolecule* parent = getParent[mol];
                if (parent != nullptr) {
                    placeMolResidueLigandStyle(mol, parent);
                }
            }
            for (sketcherMinimizerBond* b : mol->m_proximityRelations) {
                if (!b->startAtom->molecule
                         ->isPlaced) { // will place a molecule twice if it has
                                       // two relations with mol. This is safe
                                       // cause of the continue for mol->isPlace
                                       // when looping the second time
                    q.push(b->startAtom->molecule);
                    getParent[b->startAtom->molecule] = mol;
                }
                if (!b->endAtom->molecule->isPlaced) {
                    q.push(b->endAtom->molecule);
                    getParent[b->endAtom->molecule] = mol;
                }
            }
        }

    } else {

        vector<proximityData> proximityDataVector =
            buildProximityDataVector(proximityMols, molMap);

        rotateMoleculesWithProximityRelations(proximityMols, molMap,
                                              proximityDataVector);

        translateMoleculesWithProximityRelations(
            proximityMols, molMap, templateCenters, proximityDataVector);
    }
}

void sketcherMinimizer::placeMolResidueLigandStyle(
    sketcherMinimizerMolecule* mol, sketcherMinimizerMolecule* parent)
{

    int n = 0;
    sketcherMinimizerPointF parentV(0, 0);
    sketcherMinimizerPointF parentAdditionV(0, 0);
    sketcherMinimizerPointF v(0, 0);
    sketcherMinimizerPointF additionV(0,
                                      0); // actually using line to centroid, to
                                          // orient the molecule away from the
                                          // ligand
    sketcherMinimizerPointF cent = mol->center();

    for (sketcherMinimizerBond* b : mol->m_proximityRelations) {
        sketcherMinimizerAtom *at = nullptr, *parentAt = nullptr;
        if (b->startAtom->molecule == parent) {
            parentAt = b->startAtom;
            at = b->endAtom;
        } else if (b->endAtom->molecule == parent) {
            at = b->startAtom;
            parentAt = b->endAtom;
        }
        if (at == nullptr || parentAt == nullptr) {
            continue;
        }
        n++;
        sketcherMinimizerPointF paddV = parentAt->getSingleAdditionVector();
        if (b->isResidueInteraction()) {
            auto* ri = static_cast<sketcherMinimizerResidueInteraction*>(b);
            if (ri->startAtom->molecule == parent &&
                !ri->m_otherStartAtoms.empty()) {
                paddV = sketcherMinimizerAtom::getSingleAdditionVector(
                    ri->getAllStartAtoms());
            } else if (ri->endAtom->molecule == parent &&
                       !ri->m_otherEndAtoms.empty()) {
                paddV = sketcherMinimizerAtom::getSingleAdditionVector(
                    ri->getAllEndAtoms());
            }
        }

        paddV.normalize();
        paddV *= bondLength * 3;
        parentV += parentAt->coordinates;
        parentAdditionV += paddV;
        additionV += at->coordinates - cent;
        v += at->coordinates;
    }
    if (n > 0) {
        v /= n;
        parentV /= n;
        parentAdditionV /= n;
        additionV /= n;
        sketcherMinimizerPointF startingPos = parentV + parentAdditionV;
        startingPos = exploreGridAround(startingPos, 15, 10);

        auto signedAngle = sketcherMinimizerMaths::signedAngle(
            startingPos - parentV, sketcherMinimizerPointF(0, 0), -additionV);
        auto angle = static_cast<float>(signedAngle / 180 * M_PI);
        float s = sin(angle);
        float c = cos(angle);

        for (sketcherMinimizerAtom* a : mol->_atoms) {
            a->coordinates -= v;
            a->coordinates.rotate(s, c);
            a->coordinates += startingPos;
            a->coordinates.round();
        }

        flipIfCrossingInteractions(mol);

        sketcherMinimizerPointF avoidClashV = exploreMolPosition(
            mol, 15,
            bondLength *
                0.5); // explore positions on a grid of points to solve clashes
        for (sketcherMinimizerAtom* a : mol->_atoms) {
            a->coordinates += avoidClashV;
        }
    }
    mol->isPlaced = true;
}

void sketcherMinimizer::flipIfCrossingInteractions(
    sketcherMinimizerMolecule* mol)
{

    for (unsigned int bb = 0; bb < mol->m_proximityRelations.size() - 1; bb++) {
        bool out = false;
        sketcherMinimizerBond* pr1 = mol->m_proximityRelations[bb];
        if (pr1->startAtom->molecule == pr1->endAtom->molecule) {
            continue;
        }
        if (!(pr1->startAtom->molecule->isPlaced ||
              pr1->startAtom->molecule == mol)) {
            continue;
        }
        if (!(pr1->endAtom->molecule->isPlaced ||
              pr1->endAtom->molecule == mol)) {
            continue;
        }

        for (unsigned int bb2 = bb + 1; bb2 < mol->m_proximityRelations.size();
             bb2++) {
            sketcherMinimizerBond* pr2 = mol->m_proximityRelations[bb2];
            if (pr2->startAtom->molecule == pr2->endAtom->molecule) {
                continue;
            }
            if (!(pr2->startAtom->molecule->isPlaced ||
                  pr2->startAtom->molecule == mol)) {
                continue;
            }
            if (!(pr2->endAtom->molecule->isPlaced ||
                  pr2->endAtom->molecule == mol)) {
                continue;
            }

            if (sketcherMinimizerMaths::intersectionOfSegments(
                    pr1->startAtom->coordinates, pr1->endAtom->coordinates,
                    pr2->startAtom->coordinates, pr2->endAtom->coordinates)) {
                /* mirror the coordinates */
                sketcherMinimizerAtom* p1 = nullptr;
                sketcherMinimizerAtom* p2 = nullptr;
                if (pr1->startAtom->molecule == mol) {
                    p1 = pr1->startAtom;
                } else if (pr1->endAtom->molecule == mol) {
                    p1 = pr1->endAtom;
                }

                if (pr2->startAtom->molecule == mol) {
                    p2 = pr2->startAtom;
                } else if (pr2->endAtom->molecule == mol) {
                    p2 = pr2->endAtom;
                }
                if (p1 && p2) {
                    sketcherMinimizerPointF middleP =
                        p1->coordinates + p2->coordinates;
                    middleP *= 0.5;
                    sketcherMinimizerPointF p1p2V =
                        p1->coordinates - p2->coordinates;
                    p1p2V.normalize();

                    for (sketcherMinimizerAtom* a : mol->_atoms) {

                        sketcherMinimizerPointF v2 = a->coordinates - middleP;
                        float dot =
                            sketcherMinimizerMaths::dotProduct(p1p2V, v2);
                        sketcherMinimizerPointF parallel = p1p2V;
                        parallel *= dot; // parallel component of v2

                        a->coordinates -= 2 * parallel;
                        a->coordinates.round();
                    }

                    out = true;
                    break;
                }
            }
        }
        if (out) {
            break;
        }
    }
}

void sketcherMinimizer::arrangeMultipleMolecules()
{
    for (auto& _residue : m_residues) { // replace residues
        _residue->coordinatesSet = false;
    }
    if (m_molecules.size() > 1) {
        // find centers for molecules bound by proximity relations
        vector<sketcherMinimizerMolecule*> proximityMols;
        for (sketcherMinimizerMolecule* mol : m_molecules) {
            if (!mol->m_proximityRelations.empty()) {
                proximityMols.push_back(mol);
            }
        }
        sketcherMinimizerPointF center(0.f, 0.f);
        if (!proximityMols.empty()) {
            placeMoleculesWithProximityRelations(proximityMols);
        } else {
            int maxI = 0;
            size_t maxSize = m_molecules[0]->_atoms.size();
            for (unsigned int i = 0; i < m_molecules.size(); i++) {
                sketcherMinimizerMolecule* m = m_molecules[i];
                size_t size = m->_atoms.size();
                if (size > maxSize) {
                    maxI = i;
                    maxSize = size;
                }
            }

            sketcherMinimizerMolecule* centralMol = m_molecules[maxI];
            centralMol->isPlaced = true;
            for (sketcherMinimizerAtom* a : m_atoms)
                a->_generalUseVisited =
                    false; // using _generalUseVisited to keep track of charged
                           // atoms that have already been used for counterions

            center = centralMol->center();
        }

        // placing non counterions
        bool foundCounterion = true;
        while (foundCounterion) {
            foundCounterion = false;
            for (sketcherMinimizerMolecule* mol : m_molecules) {
                bool residue = false;
                if (mol->_atoms.size() == 1) {
                    if (mol->_atoms[0]->isResidue()) {
                        residue = true;
                    }
                }
                if (!residue) {
                    if (mol->hasConstrainedFragments) {
                        mol->isPlaced = true;
                    }
                }
                if (mol->hasFixedFragments) {
                    mol->isPlaced = true;
                }

                if (mol->isPlaced) {
                    continue;
                }
                if (residue) {
                    continue;
                }
                int charge = mol->totalCharge();
                if (charge == 0) {
                    sketcherMinimizerPointF counterionMin, counterionMax;
                    sketcherMinimizerPointF placeNextTo = center;
                    mol->boundingBox(counterionMin, counterionMax);
                    float counteriondx =
                        (counterionMax.x() - counterionMin.x()) * .5f;
                    float counteriondy =
                        (counterionMax.y() - counterionMin.y()) * .5f;
                    sketcherMinimizerPointF counterionCenter =
                        (counterionMax + counterionMin) * .5f;

                    foundCounterion = true;
                    // explore a grid around placeNextTo to find a suitable
                    // place

                    sketcherMinimizerPointF centerOfGrid = placeNextTo;

                    float gridD = bondLength;
                    placeNextTo = exploreGridAround(centerOfGrid, 10, gridD,
                                                    counteriondx, counteriondy);

                    for (sketcherMinimizerAtom* a : mol->_atoms) {
                        a->coordinates += placeNextTo - counterionCenter;
                    }

                    mol->isPlaced = true;
                }
            }
        }

        // placing counterions
        foundCounterion = true;
        while (foundCounterion) {
            foundCounterion = false;
            for (sketcherMinimizerMolecule* mol : m_molecules) {
                if (mol->isPlaced) {
                    continue;
                }
                int charge = mol->totalCharge();
                if (charge != 0) {
                    sketcherMinimizerPointF counterionMin, counterionMax;
                    sketcherMinimizerPointF placeNextTo = center;
                    mol->boundingBox(counterionMin, counterionMax);
                    float counteriondx =
                        (counterionMax.x() - counterionMin.x()) * .5f;
                    float counteriondy =
                        (counterionMax.y() - counterionMin.y()) * .5f;
                    sketcherMinimizerPointF counterionCenter =
                        (counterionMax + counterionMin) * .5f;

                    foundCounterion = true;

                    // find an already placed charged atom to place the
                    // counterion next to

                    for (sketcherMinimizerMolecule* m : m_molecules) {
                        bool found = false;
                        if (!m->isPlaced) {
                            continue;
                        }
                        for (sketcherMinimizerAtom* a : m->_atoms) {
                            if (a->charge == 0) {
                                continue;
                            }
                            if (a->_generalUseVisited) {
                                continue;
                            }
                            if (a->charge * charge < 0) {
                                a->_generalUseVisited = true;
                                placeNextTo = a->coordinates;
                                found = true;
                                break;
                            }
                        }
                        if (found) {
                            break;
                        }
                    }

                    // explore a grid around placeNextTo to find a suitable
                    // place

                    sketcherMinimizerPointF centerOfGrid = placeNextTo;
                    float gridD = bondLength;

                    placeNextTo =
                        exploreGridAround(centerOfGrid, 10, gridD, counteriondx,
                                          counteriondy, bondLength * 0.8);

                    for (sketcherMinimizerAtom* a : mol->_atoms) {
                        a->coordinates += placeNextTo - counterionCenter;
                    }

                    mol->isPlaced = true;
                }
            }
        }
        vector<sketcherMinimizerAtom*> ligandAtoms;
        for (sketcherMinimizerAtom* a : m_atoms)
            if (a->m_isLigand) {
                ligandAtoms.push_back(a);
            }
        placeResidues(ligandAtoms);
    }
}

void sketcherMinimizer::initializeFragments()
{

    if (_fragments.empty()) {
        cerr << "Sketcherlibs warning: no fragments to initialize" << endl;
        return;
    }

    for (sketcherMinimizerFragment* indf : m_independentFragments) {
        // recursively assign it to children
        assignNumberOfChildrenAtomsFromHere(indf);
    }

    for (sketcherMinimizerFragment* f : _fragments) {
        m_fragmentBuilder.initializeCoordinates(f);
    }

    for (sketcherMinimizerFragment* indf : m_independentFragments) {
        // recursively assign it to children
        assignLongestChainFromHere(indf);
    }
}

bool sketcherMinimizer::alignWithParentDirectionConstrained(
    sketcherMinimizerFragment* fragment,
    const sketcherMinimizerPointF& position, float angle)
{
    vector<sketcherMinimizerPointF> templates, plainCoordinates,
        flippedCoordinates;
    float sine = sin(angle);
    float cosine = cos(angle);
    for (const auto& atom : fragment->_coordinates) {
        if (atom.first->constrained) {
            sketcherMinimizerPointF plainCoordinatesAtom = atom.second;
            sketcherMinimizerPointF flippedCoordinatesAtom(
                plainCoordinatesAtom.x(), -plainCoordinatesAtom.y());
            plainCoordinatesAtom.rotate(sine, cosine);
            flippedCoordinatesAtom.rotate(sine, cosine);
            templates.push_back(atom.first->templateCoordinates);
            plainCoordinates.push_back(plainCoordinatesAtom + position);
            flippedCoordinates.push_back(flippedCoordinatesAtom + position);
        }
    }
    float scorePlain =
        roundToTwoDecimalDigits(RMSD(templates, plainCoordinates));
    float scoreFlipped =
        roundToTwoDecimalDigits(RMSD(templates, flippedCoordinates));
    return (scoreFlipped < scorePlain);
}

vector<sketcherMinimizerBond*>
sketcherMinimizer::getAllTerminalBonds(sketcherMinimizerFragment* fragment)
{
    vector<sketcherMinimizerBond*> bonds;
    for (auto bond : fragment->getBonds()) {
        if (bond->isResidueInteraction()) {
            continue;
        }
        if (bond->startAtom->neighbors.size() == 1 ||
            bond->endAtom->neighbors.size() == 1) {
            bonds.push_back(bond);
        }
    }
    for (auto child : fragment->_children) {
        bonds.push_back(child->_bondToParent);
    }
    if (fragment->getParent()) {
        bonds.push_back(fragment->_bondToParent);
    }
    return bonds;
}

vector<std::pair<sketcherMinimizerPointF, float>>
sketcherMinimizer::findDirectionsToAlignWith(
    sketcherMinimizerFragment* fragment)
{
    vector<std::pair<sketcherMinimizerPointF, float>> chainDirs;

    sketcherMinimizerPointF origin =
        (fragment->_bondToParent->startAtom->coordinates +
         fragment->_bondToParent->endAtom->coordinates) *
        0.5;
    vector<sketcherMinimizerBond*> parentEndBonds =
        getAllTerminalBonds(fragment->getParent());
    for (auto bond : parentEndBonds) {
        if (bond->endAtom->fragment == fragment) {
            continue;
        }
        sketcherMinimizerPointF direction =
            origin -
            (bond->startAtom->coordinates + bond->endAtom->coordinates) * 0.5;
        direction.normalize();
        float score = 1.f;
        if (bond->bondOrder == 2) {
            score *= SCORE_MULTIPLIER_FOR_DOUBLE_BONDS;
        }
        if ((bond->startAtom->neighbors.size() == 1 &&
             bond->startAtom->atomicNumber != 6) ||
            (bond->endAtom->neighbors.size() == 1 &&
             bond->endAtom->atomicNumber != 6)) {
            score *= SCORE_MULTIPLIER_FOR_SINGLE_BONDED_HETEROATOMS;
        }
        if (bond->endAtom->fragment != fragment->getParent() ||
            bond->startAtom->fragment != fragment->getParent()) {
            score = bond->endAtom->fragment->longestChainFromHere *
                    SCORE_MULTIPLIER_FOR_FRAGMENTS;
            if (fragment->getParent()->getParent() &&
                bond->startAtom->fragment ==
                    fragment->getParent()->getParent()) {
                score *= 100;
            }
        }
        chainDirs.emplace_back(direction, score);
    }
    return chainDirs;
}

float sketcherMinimizer::testAlignment(
    const sketcherMinimizerPointF& direction,
    const std::pair<sketcherMinimizerPointF, float>& templat)
{
    float dot = sketcherMinimizerMaths::dotProduct(direction, templat.first);
    if (dot < 0) {
        dot = 0;
    }
    float score = dot * dot;
    if (dot > 1 - SKETCHER_EPSILON) {
        score += 1000;
    }
    score *= templat.second;
    return score;
}

sketcherMinimizerPointF sketcherMinimizer::scoreDirections(
    sketcherMinimizerFragment* fragment, float angle,
    const vector<std::pair<sketcherMinimizerPointF, float>>& directions,
    bool& invert)
{
    float sine = sin(angle);
    float cosine = cos(angle);
    float bestScore = 0.f;
    sketcherMinimizerPointF bestDirection(1.f, 0.f);
    vector<sketcherMinimizerBond*> terminalBonds =
        getAllTerminalBonds(fragment);
    for (auto bond : terminalBonds) {
        if (bond->startAtom->fragment != fragment) {
            continue;
        }
        sketcherMinimizerPointF bondDirectionPlain =
            (fragment->_coordinates[bond->startAtom] +
             fragment->_coordinates[bond->endAtom]) *
                0.5 -
            sketcherMinimizerPointF(-bondLength * 0.5, 0);
        bondDirectionPlain.normalize();
        sketcherMinimizerPointF bondDirectionInverted(bondDirectionPlain.x(),
                                                      -bondDirectionPlain.y());
        bondDirectionPlain.rotate(sine, cosine);
        bondDirectionInverted.rotate(sine, cosine);
        float scoreModifier = 1.f;
        if (bond->bondOrder == 2) {
            scoreModifier *= SCORE_MULTIPLIER_FOR_DOUBLE_BONDS;
        }
        if ((bond->startAtom->neighbors.size() == 1 &&
             bond->startAtom->atomicNumber != 6) ||
            (bond->endAtom->neighbors.size() == 1 &&
             bond->endAtom->atomicNumber != 6)) {
            scoreModifier *= SCORE_MULTIPLIER_FOR_SINGLE_BONDED_HETEROATOMS;
        }
        if (bond->endAtom->fragment != fragment) {
            scoreModifier = bond->endAtom->fragment->longestChainFromHere *
                            SCORE_MULTIPLIER_FOR_FRAGMENTS;
        }
        for (const auto& direction : directions) {
            float scorePlain =
                testAlignment(bondDirectionPlain, direction) * scoreModifier;
            if (scorePlain > bestScore) {
                bestScore = scorePlain;
                bestDirection = direction.first;
                invert = false;
            }
            float scoreInverted =
                testAlignment(bondDirectionInverted, direction) * scoreModifier;
            if (scoreInverted > bestScore) {
                bestScore = scoreInverted;
                bestDirection = direction.first;
                invert = true;
            }
        }
    }
    return bestDirection;
}

bool sketcherMinimizer::alignWithParentDirectionUnconstrained(
    sketcherMinimizerFragment* fragment, float angle)
{
    vector<std::pair<sketcherMinimizerPointF, float>> directions =
        findDirectionsToAlignWith(fragment);
    bool invert = false;
    scoreDirections(fragment, angle, directions, invert);
    return invert;
}

void sketcherMinimizer::alignWithParentDirection(
    sketcherMinimizerFragment* f, const sketcherMinimizerPointF& position,
    float angle)
{
    // deciding which "side" the fragment will be drawn, rotating 180° around
    // the axis of its bond to parent
    if (f->fixed) {
        return;
    }
    bool invert = (f->constrained
                       ? alignWithParentDirectionConstrained(f, position, angle)
                       : alignWithParentDirectionUnconstrained(f, angle));
    if (invert) {
        for (auto& atom : f->_coordinates) {
            atom.second.setY(-atom.second.y());
        }
        for (auto atom : f->getAtoms()) {
            if (atom->hasStereochemistrySet) {
                for (auto bond : atom->bonds) {
                    bond->isWedge = !bond->isWedge;
                }
            }
        }
    }
}

void sketcherMinimizer::assignNumberOfChildrenAtomsFromHere(
    sketcherMinimizerFragment* f)
{
    size_t cumulatedNumberOfAtoms = 0;
    float cumulatedNumberOfAtomsRanks = 0;
    size_t childrenAtoms = 0;
    for (sketcherMinimizerFragment* child : f->_children) {
        assignNumberOfChildrenAtomsFromHere(child);
        cumulatedNumberOfAtoms += child->numberOfChildrenAtoms;
        cumulatedNumberOfAtomsRanks += child->numberOfChildrenAtomsRank;
        childrenAtoms += child->getAtoms().size();
    }
    f->numberOfChildrenAtoms = cumulatedNumberOfAtoms + childrenAtoms;
    f->numberOfChildrenAtomsRank =
        static_cast<float>(0.01f * cumulatedNumberOfAtomsRanks + childrenAtoms);
}

void sketcherMinimizer::assignLongestChainFromHere(sketcherMinimizerFragment* f)
{
    float longestDist = 0;
    for (sketcherMinimizerFragment* child : f->_children) {
        assignLongestChainFromHere(child);
        if (child->longestChainFromHere > longestDist) {
            longestDist = child->longestChainFromHere;
        }
    }
    sketcherMinimizerPointF positionFromParent(0.f, 0.f);
    if (f->getParent()) {
        positionFromParent =
            f->getParent()->_coordinates[f->_bondToParent->endAtom];
    }
    f->longestChainFromHere = longestDist + positionFromParent.length();
}

sketcherMinimizerBond*
sketcherMinimizer::getBond(const sketcherMinimizerAtom* a1,
                           const sketcherMinimizerAtom* a2)
{
    for (unsigned int i = 0; i < a1->neighbors.size(); i++) {
        if (a1->neighbors[i] == a2) {
            return a1->bonds[i];
        }
    }
    return nullptr;
}

sketcherMinimizerAtom*
sketcherMinimizer::pickBestAtom(vector<sketcherMinimizerAtom*>& atoms)
{

    vector<sketcherMinimizerAtom*> candidates, oldCandidates;

    {
        size_t biggestSize = atoms[0]->fragment->numberOfChildrenAtoms;
        for (sketcherMinimizerAtom* a : atoms) {
            size_t size = a->fragment->numberOfChildrenAtoms;
            if (size == biggestSize) {
                candidates.push_back(a);
            } else if (size > biggestSize) {
                biggestSize = size;
                candidates.clear();
                candidates.push_back(a);
            }
        }
        if (candidates.size() == 1) {
            return candidates[0];
        }
        oldCandidates = candidates;
        candidates.clear();
    }

    {
        float biggestSize =
            oldCandidates[0]->fragment->numberOfChildrenAtomsRank;
        for (sketcherMinimizerAtom* a : oldCandidates) {
            float size = a->fragment->numberOfChildrenAtomsRank;
            if (size == biggestSize) {
                candidates.push_back(a);
            } else if (size > biggestSize) {
                biggestSize = size;
                candidates.clear();
                candidates.push_back(a);
            }
        }
        if (candidates.size() == 1) {
            return candidates[0];
        }
        oldCandidates = candidates;
        candidates.clear();
    }

    {
        int biggestSize = oldCandidates[0]->atomicNumber;
        for (sketcherMinimizerAtom* a : oldCandidates) {
            int size = a->atomicNumber;
            if (size == biggestSize) {
                candidates.push_back(a);
            } else if (size > biggestSize) {
                biggestSize = size;
                candidates.clear();
                candidates.push_back(a);
            }
        }
        if (candidates.size() == 1) {
            return candidates[0];
        }
        oldCandidates = candidates;
        candidates.clear();
    }
    // give up
    return oldCandidates[0];
}

void sketcherMinimizer::constrainAllAtoms()
{
    //   cerr << "sketcherMinimizer::constrainAllAtoms ()"<<endl;
    for (sketcherMinimizerAtom* a : m_atoms)
        a->constrained = true;
}

void sketcherMinimizer::constrainAtoms(const vector<bool>& constrained)
{
    if (constrained.size() == m_referenceAtoms.size()) {
        for (unsigned int i = 0; i < constrained.size(); i++) {
            if (constrained[i]) {
                m_referenceAtoms[i]->constrained = true;
            }
        }
    } else {
        cerr << "warning, wrong size of vector for constrained atoms. Ignoring"
             << endl;
    }
}

void sketcherMinimizer::fixAtoms(const vector<bool>& fixed)
{
    if (fixed.size() == m_referenceAtoms.size()) {
        for (unsigned int i = 0; i < fixed.size(); i++) {
            if (fixed[i]) {
                m_referenceAtoms[i]->fixed = true;
            }
        }
    } else {
        cerr << "warning, wrong size of vector for fixed atoms. Ignoring"
             << endl;
    }
}

void sketcherMinimizer::findClosestAtomToResidues(
    const vector<sketcherMinimizerAtom*>& catoms)
{
    const vector<sketcherMinimizerAtom*>& atoms =
        catoms.empty() ? m_atoms : catoms;

    for (sketcherMinimizerAtom* r : m_residues) {
        float squareD = 9999999.f;
        sketcherMinimizerAtom* closestA = nullptr;
        for (sketcherMinimizerAtom* a : atoms) {
            if (!a->isResidue()) {
                float diffx = a->m_x3D - r->m_x3D;
                float diffy = a->m_y3D - r->m_y3D;
                float diffz = a->m_z3D - r->m_z3D;

                float newSquareD =
                    diffx * diffx + diffy * diffy + diffz * diffz;
                if (newSquareD < squareD) {
                    squareD = newSquareD;
                    closestA = a;
                }
            }
        }
        static_cast<sketcherMinimizerResidue*>(r)->m_closestLigandAtom =
            closestA;
        if (!r->m_isClashing) {
            r->m_isClashing = (squareD < RESIDUE_CLASH_DISTANCE_SQUARED);
        }
    }
    for (sketcherMinimizerBond* b : m_bonds) {
        if (b->startAtom->isResidue()) {
            static_cast<sketcherMinimizerResidue*>(b->startAtom)
                ->m_closestLigandAtom = b->endAtom;
        }
        if (b->endAtom->isResidue()) {
            static_cast<sketcherMinimizerResidue*>(b->endAtom)
                ->m_closestLigandAtom = b->startAtom;
        }
    }
}

float sketcherMinimizer::RMSD(const vector<sketcherMinimizerPointF>& templates,
                              const vector<sketcherMinimizerPointF>& points)
{
    assert(templates.size() == points.size());
    size_t counter = templates.size();
    float total = 0.f;
    for (unsigned int i = 0; i < templates.size(); i++) {
        //        cerr << templates[i].x () << ", "<< templates[i].y () << "
        //        " << points[i].x()<<", "<<points[i].y ()<<endl;
        sketcherMinimizerPointF diff = templates[i] - points[i];
        total += diff.x() * diff.x() + diff.y() * diff.y();
    }
    if (counter > 0) {
        total /= counter;
    }
    return sqrt(total);
}

void sketcherMinimizer::alignmentMatrix(
    const vector<sketcherMinimizerPointF>& ref,
    const vector<sketcherMinimizerPointF>& points, float* m)
{
    float U[4];
    float Sig[4];

    float V[4];
    float a[4];

    a[0] = 0.f;
    a[1] = 0.f;
    a[2] = 0.f;
    a[3] = 0.f;
    assert(ref.size() == points.size());
    for (unsigned int i = 0; i < ref.size(); i++) {
        a[0] += ref[i].x() * points[i].x();
        a[1] += ref[i].y() * points[i].x();
        a[2] += ref[i].x() * points[i].y();
        a[3] += ref[i].y() * points[i].y();
    }
    svd(a, U, Sig, V);
    m[0] = V[0] * U[0] + V[1] * U[1];
    m[1] = V[0] * U[2] + V[1] * U[3];
    m[2] = V[2] * U[0] + V[3] * U[1];
    m[3] = V[2] * U[2] + V[3] * U[3];
}

void sketcherMinimizer::svd(float* a, float* U, float* Sig, float* V)
{
    float a1[4];
    a1[0] = a[0];
    a1[1] = a[2];
    a1[2] = a[1];
    a1[3] = a[3];

    float Su[4];
    Su[0] = a[0] * a1[0] + a[1] * a1[2];
    Su[1] = a[0] * a1[1] + a[1] * a1[3];
    Su[2] = a[2] * a1[0] + a[3] * a1[2];
    Su[3] = a[2] * a1[1] + a[3] * a1[3];

    auto phi = static_cast<float>(0.5 * atan2(Su[1] + Su[2], Su[0] - Su[3]));
    float cphi = cos(phi);
    cphi = roundToTwoDecimalDigits(cphi);
    float sphi = sin(phi);
    sphi = roundToTwoDecimalDigits(sphi);

    U[0] = cphi * (-1);
    U[1] = -sphi;
    U[2] = sphi * (-1);
    U[3] = cphi;

    float Sw[4];
    Sw[0] = a1[0] * a[0] + a1[1] * a[2];
    Sw[1] = a1[0] * a[1] + a1[1] * a[3];
    Sw[2] = a1[2] * a[0] + a1[3] * a[2];
    Sw[3] = a1[2] * a[1] + a1[3] * a[3];

    auto theta = static_cast<float>(0.5 * atan2(Sw[1] + Sw[2], Sw[0] - Sw[3]));
    float ctheta = cos(theta);
    float stheta = sin(theta);

    float W[4];
    W[0] = ctheta;
    W[1] = -stheta;
    W[2] = stheta;
    W[3] = ctheta;

    float SUsum = Su[0] + Su[3];
    float SUdif = sqrt((Su[0] - Su[3]) * (Su[0] - Su[3]) + 4.f * Su[1] * Su[2]);

    Sig[0] = sqrt((SUsum + SUdif) * 0.5f);
    Sig[1] = 0.f;
    Sig[2] = 0.f;
    Sig[3] = sqrt((SUsum - SUdif) * 0.5f);

    float U1[4];
    U1[0] = U[0];
    U1[1] = U[2];
    U1[2] = U[1];
    U1[3] = U[3];

    float U1A[4];

    U1A[0] = U1[0] * a[0] + U1[1] * a[2];
    U1A[1] = U1[0] * a[1] + U1[1] * a[3];
    U1A[2] = U1[2] * a[0] + U1[3] * a[2];
    U1A[3] = U1[2] * a[1] + U1[3] * a[3];

    float S[4];

    S[0] = U1A[0] * W[0] + U1A[1] * W[2];
    S[1] = U1A[0] * W[1] + U1A[1] * W[3];
    S[2] = U1A[2] * W[0] + U1A[3] * W[2];
    S[3] = U1A[2] * W[1] + U1A[3] * W[3];
    S[0] = roundToTwoDecimalDigits(S[0]);
    S[1] = roundToTwoDecimalDigits(S[1]);
    S[2] = roundToTwoDecimalDigits(S[2]);
    S[3] = roundToTwoDecimalDigits(S[3]);

    float sign1 = 1.f, sign2 = 1.f;
    if (S[0] < 0.f) {
        sign1 = -1.f;
    }
    if (S[3] < 0.f) {
        sign2 = -1.f;
    }
    float C[4];

    C[0] = sign1;
    C[1] = 0.f;
    C[2] = 0.f;
    C[3] = sign2;

    V[0] = W[0] * C[0] + W[1] * C[2];
    V[1] = W[0] * C[1] + W[1] * C[3];
    V[2] = W[2] * C[0] + W[3] * C[2];
    V[3] = W[2] * C[1] + W[3] * C[3];
    V[0] = roundToTwoDecimalDigits(V[0]);
    V[1] = roundToTwoDecimalDigits(V[1]);
    V[2] = roundToTwoDecimalDigits(V[2]);
    V[3] = roundToTwoDecimalDigits(V[3]);
}

bool sketcherMinimizer::compare(const vector<sketcherMinimizerAtom*>& atoms,
                                const vector<sketcherMinimizerBond*>& bonds,
                                sketcherMinimizerMolecule* templ,
                                vector<unsigned int>& mapping)
{
    if (atoms.size() != templ->_atoms.size()) {
        return false;
    }
    vector<int> molScores, tempScores;
    int molIter = morganScores(atoms, bonds, molScores);
    int templateIter = morganScores(templ->_atoms, templ->_bonds, tempScores);

    if (molIter != templateIter) {
        return false;
    }

    size_t size = atoms.size();
    vector<bool> matrix(size * size, false);

    vector<sketcherMinimizerPointF> templateCoordinates;
    vector<vector<size_t>> molBonds;
    vector<vector<size_t>> templateBonds;

    // for double bond chirality
    //   vector < vector < int > > templateAllowedZChains;

    vector<vector<size_t>> molCisTransChains;
    vector<bool> molIsCis;

    for (unsigned int ma = 0; ma < size; ma++) {
        vector<size_t> vec;
        molBonds.push_back(vec);
    }

    for (unsigned int ta = 0; ta < size; ta++) {
        vector<size_t> vec;
        templateBonds.push_back(vec);
    }

    for (unsigned int ta = 0; ta < size; ta++) {
        templateCoordinates.push_back(templ->_atoms[ta]->coordinates);
    }
    for (auto bond : bonds) {
        if (bond->bondOrder == 2) {

            sketcherMinimizerBond* b = bond;
            if (b->m_ignoreZE) {
                continue;
            }

            bool ok = false;
            if (b->rings.empty()) {
                ok = true;
            } else {
                ok = true;
                for (auto& ring : b->rings) {
                    if (ring->_atoms.size() < MACROCYCLE) {
                        ok = false;
                    }
                }
            }
            if (!ok) {
                continue;
            }
            sketcherMinimizerAtom* startA = b->startAtomCIPFirstNeighbor();
            sketcherMinimizerAtom* endA = b->endAtomCIPFirstNeighbor();

            sketcherMinimizerAtom* startN = nullptr;
            sketcherMinimizerAtom* endN = nullptr;
            for (unsigned int i = 0; i < b->startAtom->neighbors.size(); i++) {
                if (b->startAtom->neighbors[i] == b->endAtom) {
                    continue;
                }
                bool found = false;
                for (auto atom : atoms) {
                    if (atom == b->startAtom->neighbors[i]) {
                        startN = b->startAtom->neighbors[i];
                        found = true;
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }
            for (unsigned int i = 0; i < b->endAtom->neighbors.size(); i++) {
                if (b->endAtom->neighbors[i] == b->startAtom) {
                    continue;
                }
                bool found = false;
                for (auto atom : atoms) {
                    if (atom == b->endAtom->neighbors[i]) {
                        endN = b->endAtom->neighbors[i];
                        found = true;
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }

            if (startN && endN && startA && endA) {
                bool isCis = b->isZ;
                if (startN != startA) {
                    isCis = !isCis;
                }
                if (endN != endA) {
                    isCis = !isCis;
                }
                vector<size_t> chain;
                chain.push_back(startN->_generalUseN);
                chain.push_back(b->startAtom->_generalUseN);
                chain.push_back(b->endAtom->_generalUseN);
                chain.push_back(endN->_generalUseN);
                molCisTransChains.push_back(chain);
                molIsCis.push_back(isCis);
                // cerr << "chiral double bond" <<isCis<<"    "<<b->isZ<<endl;
            }
        }
    }
    // assuming that _generalUseN is set as the index of each atom
    for (auto bond : bonds) {
        size_t in1 = bond->startAtom->_generalUseN;
        size_t in2 = bond->endAtom->_generalUseN;

        if (in1 < in2) {
            molBonds[in2].push_back(in1);
        } else {
            molBonds[in1].push_back(in2);
        }
    }
    for (auto& _bond : templ->_bonds) {
        size_t in1 = _bond->startAtom->_generalUseN;
        size_t in2 = _bond->endAtom->_generalUseN;
        if (in1 < in2) {
            templateBonds[in2].push_back(in1);
        } else {
            templateBonds[in1].push_back(in2);
        }
    }
    for (unsigned int ma = 0; ma < atoms.size(); ma++) {
        for (unsigned int ta = 0; ta < templ->_atoms.size(); ta++) {
            if (molScores[ma] == tempScores[ta]) {
                matrix[ma * size + ta] = true;
            }
        }
    }
    bool found = false;
    vector<unsigned int> solution;
    for (unsigned int i = 0; i < size; i++) {
        if (!matrix[i]) {
            continue;
        }
        checkIdentity(solution, i, matrix, templateCoordinates, molBonds,
                      templateBonds, molCisTransChains, molIsCis, size, found,
                      mapping);
        if (found) {
            break;
        }
    }
    return found;
}

void sketcherMinimizer::checkIdentity(
    vector<unsigned int> solution, int newSol, vector<bool>& matrix,
    vector<sketcherMinimizerPointF>& templateCoordinates,
    vector<vector<size_t>>& molBonds, vector<vector<size_t>>& templateBonds,
    vector<vector<size_t>>& molCisTransChains, vector<bool>& molIsCis,
    size_t size, bool& found, vector<unsigned int>& mapping)
{
    solution.push_back(newSol);
    if (solution.size() == size) {
        // check double bonds stereo
        bool doubleBondsStereoOk = true;

        for (size_t i = 0; i < molCisTransChains.size(); i++) {
            sketcherMinimizerPointF p1 =
                templateCoordinates[solution[molCisTransChains[i][0]]];
            sketcherMinimizerPointF p2 =
                templateCoordinates[solution[molCisTransChains[i][1]]];
            sketcherMinimizerPointF p3 =
                templateCoordinates[solution[molCisTransChains[i][2]]];
            sketcherMinimizerPointF p4 =
                templateCoordinates[solution[molCisTransChains[i][3]]];

            if (molIsCis[i] !=
                sketcherMinimizerMaths::sameSide(p1, p4, p2, p3)) {
                doubleBondsStereoOk = false;
                break;
            }
        }

        if (doubleBondsStereoOk) {
            found = true;
            mapping = solution;
        }
    } else {
        for (unsigned int i = 0; i < size; i++) {
            if (found) {
                break;
            }
            if (!(matrix[solution.size() * size + i])) {
                continue;
            }
            bool check = true;
            for (unsigned int ss : solution) {
                if (ss == i) {
                    check = false;
                    break;
                }
            }
            if (check) {
                for (const auto& molBond : molBonds[solution.size()]) {
                    check = false;
                    unsigned int high = i;
                    unsigned int low = solution[molBond];
                    if (low > high) {
                        std::swap(low, high);
                    }
                    for (auto ch : templateBonds[high]) {
                        if (ch == low) {
                            check = true;
                            break;
                        }
                    }
                    if (!check) {
                        break;
                    }
                }
            }

            if (check) {
                checkIdentity(solution, i, matrix, templateCoordinates,
                              molBonds, templateBonds, molCisTransChains,
                              molIsCis, size, found, mapping);
            }
        }
    }
}

void sketcherMinimizer::setTemplateFileDir(string dir)
{
    sketcherMinimizer::m_templates.setTemplateDir(std::move(dir));
}

static void normalizeTemplate(sketcherMinimizerMolecule* mol)
{
    // normalize bond length and set _generalUseN on atoms.
    vector<float> dds;
    vector<int> ns;
    for (auto& _bond : mol->_bonds) {
        sketcherMinimizerPointF v =
            _bond->startAtom->coordinates - _bond->endAtom->coordinates;
        float dd = v.x() * v.x() + v.y() * v.y();
        bool found = false;
        for (unsigned int j = 0; j < dds.size(); ++j) {
            if (dd * 0.9 < dds[j] && dd * 1.1 > dds[j]) {
                ++ns[j];
                found = true;
                break;
            }
        }
        if (!found) {
            dds.push_back(dd);
            ns.push_back(1);
        }
    }

    if (!dds.empty()) {
        int maxI = 0;
        for (unsigned int i = 0; i < ns.size(); ++i) {
            if (ns[i] > ns[maxI]) {
                maxI = i;
            }
        }

        float f = sqrt(dds[maxI]);
        for (unsigned int i = 0; i < mol->_atoms.size(); ++i) {
            mol->_atoms[i]->coordinates /= f;
            mol->_atoms[i]->_generalUseN = i;
        }
    }
}

#ifdef USE_MAEPARSER
static string getTempFileProjDir()
{
    return sketcherMinimizer::m_templates.getTemplateDir();
}

static string getUserTemplateFileName()
{
    const string suffix = "user_templates.mae"; // whatever you wish
    return getTempFileProjDir() + suffix;
}

static void loadTemplate(const string& filename,
                         vector<sketcherMinimizerMolecule*>& templates)
{
    auto pFile = fopen(filename.c_str(), "r");
    if (pFile == nullptr) {
        return;
    }
    schrodinger::mae::Reader r(pFile);
    std::shared_ptr<schrodinger::mae::Block> b;
    while ((b = r.next(mae::CT_BLOCK)) != nullptr) {
        auto mol = mol_from_mae_block(*b);
        normalizeTemplate(mol);
        templates.push_back(mol);
    }
    fclose(pFile);
}
#endif

void sketcherMinimizer::loadTemplates()
{
    static bool loaded = false;
    if (loaded || !m_templates.getTemplates().empty()) {
        return;
    }
    m_templates.getTemplates() = coordgen_templates();
    for (auto* m : m_templates.getTemplates()) {
        normalizeTemplate(m);
    }
#ifdef USE_MAEPARSER

    string filename = getUserTemplateFileName();
    loadTemplate(filename, m_templates.getTemplates());

#endif
    loaded = true;
}

int sketcherMinimizer::morganScores(const vector<sketcherMinimizerAtom*>& atoms,
                                    const vector<sketcherMinimizerBond*>& bonds,
                                    vector<int>& oldScores)
{
    // assumes that atom[i]->_generalUseN = i
    if (atoms.size() < 2) {
        return 0;
    }
    oldScores = vector<int>(atoms.size(), 1);
    vector<int> newScores(atoms.size(), 0);
    vector<int> orderedScores;
    bool goOn = false;
    int n = 0;
    size_t idx1, idx2;
    size_t oldTies = atoms.size();
    unsigned int i = 0, j = 0;
    do {
        ++n;
        for (i = 0; i < bonds.size(); ++i) {
            idx1 = bonds[i]->startAtom->_generalUseN;
            idx2 = bonds[i]->endAtom->_generalUseN;
            newScores[idx1] += oldScores[idx2];
            newScores[idx2] += oldScores[idx1];
        }
        orderedScores = newScores;
        stable_sort(orderedScores.begin(), orderedScores.end());
        size_t newTies = 0;
        for (j = 1; j < orderedScores.size(); j++) {
            if (orderedScores[j] == orderedScores[j - 1]) {
                newTies++;
            }
        }
        if (newTies < oldTies) {
            goOn = true;
            oldTies = newTies;
            oldScores = newScores;
        } else {
            goOn = false;
        }
    } while (goOn);
    return n;
}

// interactions with m_minimizer and m_fragmentBuilder
const std::vector<sketcherMinimizerBendInteraction*>&
sketcherMinimizer::getBendInteractions() const
{
    return m_minimizer.getBendInteractions();
}

const std::vector<sketcherMinimizerStretchInteraction*>&
sketcherMinimizer::getStretchInteractions() const
{
    return m_minimizer.getStretchInteractions();
}

const std::vector<sketcherMinimizerInteraction*>&
sketcherMinimizer::getInteractions() const
{
    return m_minimizer.getInteractions();
}

std::set<sketcherMinimizerAtom*> sketcherMinimizer::getChetoCs(
    const std::vector<sketcherMinimizerAtom*>& allAtoms) const
{
    return m_minimizer.getChetoCs(allAtoms);
}

std::set<sketcherMinimizerAtom*> sketcherMinimizer::getAminoNs(
    const std::vector<sketcherMinimizerAtom*>& allAtoms) const
{
    return m_minimizer.getAminoNs(allAtoms);
}

std::set<sketcherMinimizerAtom*> sketcherMinimizer::getAlphaCs(
    const std::vector<sketcherMinimizerAtom*>& allAtoms,
    const std::set<sketcherMinimizerAtom*>& chetoCs,
    const std::set<sketcherMinimizerAtom*>& aminoNs) const
{
    return m_minimizer.getAlphaCs(allAtoms, chetoCs, aminoNs);
}

void sketcherMinimizer::clearInteractions()
{
    m_minimizer.clearInteractions();
}

void sketcherMinimizer::addInteractionsOfMolecule(
    sketcherMinimizerMolecule* molecule, bool intrafragmentClashes)
{
    m_minimizer.addInteractionsOfMolecule(molecule, intrafragmentClashes);
}

void sketcherMinimizer::minimizeMolecule(sketcherMinimizerMolecule* molecule)
{
    CoordgenFragmenter::splitIntoFragments(molecule);
    m_minimizer.minimizeMolecule(molecule);
};

void sketcherMinimizer::forceFieldMinimize()
{
    m_minimizer.run();
}

void sketcherMinimizer::addExtraInteraction(
    sketcherMinimizerMolecule* molecule,
    sketcherMinimizerInteraction* interaction)
{
    m_minimizer.addExtraInteraction(molecule, interaction);
}

void sketcherMinimizer::setEvenAngles(bool b)
{
    m_fragmentBuilder.m_evenAngles = b;
    m_minimizer.m_evenAngles = b;
}
void sketcherMinimizer::setSkipMinimization(bool b)
{
    m_minimizer.skipMinimization = b;
}
void sketcherMinimizer::setForceOpenMacrocycles(bool b)
{
    m_fragmentBuilder.setForceOpenMacrocycles(b);
}

void sketcherMinimizer::buildFromFragments(bool b)
{
    m_minimizer.buildFromFragments(b);
}
bool sketcherMinimizer::avoidClashesOfMolecule(
    sketcherMinimizerMolecule* molecule,
    const std::vector<sketcherMinimizerInteraction*>& extraInteractions)
{
    return m_minimizer.avoidClashesOfMolecule(molecule, extraInteractions);
}

void sketcherMinimizer::addExtraBond(sketcherMinimizerBond* bond)
{
    m_extraBonds.push_back(bond);
}

CoordgenTemplates sketcherMinimizer::m_templates;
