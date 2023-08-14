//
// Copyright (C) David Cosgrove 2023
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <regex>
#include <set>

#include <GraphMol/MolOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "RascalResult.h"

namespace RDKit {

namespace RascalMCES {

RascalResult::RascalResult(const RDKit::ROMol &mol1, const RDKit::ROMol &mol2,
                           const std::vector<std::vector<int>> &adjMatrix1,
                           const std::vector<std::vector<int>> &adjMatrix2,
                           const std::vector<unsigned int> &clique,
                           const std::vector<std::pair<int, int>> &vtx_pairs,
                           bool timedOut, bool swapped, double tier1Sim,
                           double tier2Sim, bool ringMatchesRingOnly,
                           bool singleLargestFrag, int maxFragSep)
    : d_timedOut(timedOut),
      d_tier1Sim(tier1Sim),
      d_tier2Sim(tier2Sim),
      d_ringMatchesRingOnly(ringMatchesRingOnly),
      d_maxFragSep(maxFragSep) {
  const std::vector<std::vector<int>> *mol1AdjMatrix;
  if (swapped) {
    d_mol1.reset(new RDKit::ROMol(mol2));
    d_mol2.reset(new RDKit::ROMol(mol1));
    mol1AdjMatrix = &adjMatrix2;
  } else {
    d_mol1.reset(new RDKit::ROMol(mol1));
    d_mol2.reset(new RDKit::ROMol(mol2));
    mol1AdjMatrix = &adjMatrix1;
  }

  extractClique(clique, vtx_pairs, swapped, d_bondMatches);
  matchCliqueAtoms(*mol1AdjMatrix);
  if (d_maxFragSep != -1) {
    applyMaxFragSep();
  }
  if (singleLargestFrag) {
    largestFragOnly();
  }
}

RascalResult::RascalResult(double tier1Sim, double tier2Sim)
    : d_tier1Sim(tier1Sim), d_tier2Sim(tier2Sim) {}

RascalResult::RascalResult(const RascalResult &other)
    : d_bondMatches(other.d_bondMatches),
      d_atomMatches(other.d_atomMatches),
      d_smarts(other.d_smarts),
      d_timedOut(other.d_timedOut),
      d_tier1Sim(other.d_tier1Sim),
      d_tier2Sim(other.d_tier2Sim),
      d_numFrags(other.d_numFrags),
      d_ringNonRingBondScore(other.d_ringNonRingBondScore),
      d_atomMatchScore(other.d_atomMatchScore),
      d_maxDeltaAtomAtomDist(other.d_maxDeltaAtomAtomDist),
      d_largestFragSize(other.d_largestFragSize) {
  if (other.d_mol1) {
    d_mol1.reset(new ROMol(*other.d_mol1));
  }
  if (other.d_mol2) {
    d_mol2.reset(new ROMol(*other.d_mol2));
  }
  if (other.d_mcesMol) {
    d_mcesMol.reset(new ROMol(*other.d_mcesMol));
  }
}

RascalResult &RascalResult::operator=(const RascalResult &other) {
  if (this == &other) {
    return *this;
  }
  d_bondMatches = other.d_bondMatches;
  d_atomMatches = other.d_atomMatches;
  d_smarts = other.d_smarts;
  d_timedOut = other.d_timedOut;
  d_numFrags = other.d_numFrags;
  d_ringNonRingBondScore = other.d_ringNonRingBondScore;
  d_atomMatchScore = other.d_atomMatchScore;
  d_maxDeltaAtomAtomDist = other.d_maxDeltaAtomAtomDist;
  d_largestFragSize = other.d_largestFragSize;
  if (other.d_mol1) {
    d_mol1.reset(new ROMol(*other.d_mol1));
  }
  if (other.d_mol2) {
    d_mol2.reset(new ROMol(*other.d_mol2));
  }
  if (other.d_mcesMol) {
    d_mcesMol.reset(new ROMol(*other.d_mcesMol));
  }
  return *this;
}

void RascalResult::largestFragOnly() { largestFragsOnly(1); }

void RascalResult::largestFragsOnly(int numFrags) {
  std::unique_ptr<RDKit::ROMol> mol1_frags(makeMolFrags(1));
  // getMolFrags() returns boost::shared_ptr.  Ho-hum.
  auto frags = RDKit::MolOps::getMolFrags(*mol1_frags, false);
  if (numFrags < 1 || frags.size() < numFrags) {
    return;
  }
  std::sort(frags.begin(), frags.end(),
            [](const boost::shared_ptr<ROMol> &f1,
               const boost::shared_ptr<ROMol> &f2) -> bool {
              return f1->getNumAtoms() > f2->getNumAtoms();
            });
  frags.erase(frags.begin() + numFrags, frags.end());
  rebuildFromFrags(frags);
}

void RascalResult::trimSmallFrags(int minFragSize) {
  std::unique_ptr<RDKit::ROMol> mol1_frags(makeMolFrags(1));
  // getMolFrags() returns boost::shared_ptr.  Ho-hum.
  auto frags = RDKit::MolOps::getMolFrags(*mol1_frags, false);
  for (auto &frag : frags) {
    if (frag->getNumAtoms() < minFragSize) {
      frag.reset();
    }
  }
  frags.erase(std::remove_if(
                  frags.begin(), frags.end(),
                  [](const boost::shared_ptr<ROMol> &f) -> bool { return !f; }),
              frags.end());
  rebuildFromFrags(frags);
}

double RascalResult::similarity() const {
  if (!d_mol1 || !d_mol2) {
    return 0.0;
  }
  return johnsonSimilarity(d_bondMatches, d_atomMatches, *d_mol1, *d_mol2);
}

void RascalResult::rebuildFromFrags(
    const std::vector<boost::shared_ptr<ROMol>> &frags) {
  // Force the re-creation of the SMARTS and other properties next time
  // they-re needed.
  d_smarts = "";
  d_maxFragSep = -1;
  d_ringNonRingBondScore = -1;
  d_maxDeltaAtomAtomDist = -1;
  d_largestFragSize = -1;

  std::set<int> fragAtoms, fragBonds;
  for (const auto &f : frags) {
    for (auto atom : f->atoms()) {
      if (atom->hasProp("ORIG_INDEX")) {
        fragAtoms.insert(atom->getProp<int>("ORIG_INDEX"));
      }
    }
    for (auto bond : f->bonds()) {
      if (bond->hasProp("ORIG_INDEX")) {
        fragBonds.insert(bond->getProp<int>("ORIG_INDEX"));
      }
    }
  }
  std::vector<std::pair<int, int>> newAtomMatches;
  for (const auto &am : d_atomMatches) {
    if (fragAtoms.find(am.first) != fragAtoms.end()) {
      newAtomMatches.push_back(am);
    }
  }
  d_atomMatches = newAtomMatches;
  std::vector<std::pair<int, int>> new_bond_matches;
  for (const auto &bm : d_bondMatches) {
    if (fragBonds.find(bm.first) != fragBonds.end()) {
      new_bond_matches.push_back(bm);
    }
  }
  d_bondMatches = new_bond_matches;
  d_numFrags = frags.size();
  d_largestFragSize = frags.empty() ? 0 : frags.front()->getNumAtoms();
}

std::string RascalResult::createSmartsString() const {
  if (!d_mol1 || !d_mol2) {
    return "";
  }
  std::unique_ptr<RDKit::RWMol> smartsMol(new RDKit::RWMol);
  std::map<int, unsigned int> atomMap;
  auto mol1Rings = d_mol1->getRingInfo();
  auto mol2Rings = d_mol2->getRingInfo();
  for (const auto &am : d_atomMatches) {
    RDKit::QueryAtom a;
    auto mol1Atom = d_mol1->getAtomWithIdx(am.first);
    a.setQuery(RDKit::makeAtomNumQuery(mol1Atom->getAtomicNum()));
    auto mol2Atom = d_mol2->getAtomWithIdx(am.second);
    if (mol1Atom->getAtomicNum() != mol2Atom->getAtomicNum()) {
      a.expandQuery(RDKit::makeAtomNumQuery(mol2Atom->getAtomicNum()),
                    Queries::COMPOSITE_OR);
    }
    if (mol1Atom->getIsAromatic() && mol2Atom->getIsAromatic()) {
      a.expandQuery(RDKit::makeAtomAromaticQuery(), Queries::COMPOSITE_AND,
                    true);
    } else if (!mol1Atom->getIsAromatic() && !mol2Atom->getIsAromatic()) {
      a.expandQuery(RDKit::makeAtomAliphaticQuery(), Queries::COMPOSITE_AND,
                    true);
    }
    if (d_ringMatchesRingOnly && !mol1Atom->getIsAromatic() &&
        !mol2Atom->getIsAromatic() &&
        mol1Rings->numAtomRings(mol1Atom->getIdx()) &&
        mol2Rings->numAtomRings(mol2Atom->getIdx())) {
      a.expandQuery(RDKit::makeAtomInRingQuery(), Queries::COMPOSITE_AND, true);
    }
    auto ai = smartsMol->addAtom(&a);
    atomMap.insert(std::make_pair(am.first, ai));
  }

  for (const auto &bm : d_bondMatches) {
    RDKit::QueryBond b;
    auto mol1Bond = d_mol1->getBondWithIdx(bm.first);
    b.setBeginAtomIdx(atomMap[mol1Bond->getBeginAtomIdx()]);
    b.setEndAtomIdx(atomMap[mol1Bond->getEndAtomIdx()]);
    b.setQuery(makeBondOrderEqualsQuery(mol1Bond->getBondType()));
    auto mol2Bond = d_mol2->getBondWithIdx(bm.second);
    if (mol1Bond->getBondType() != mol2Bond->getBondType()) {
      b.expandQuery(makeBondOrderEqualsQuery(mol2Bond->getBondType()),
                    Queries::COMPOSITE_OR);
    }
    if (d_ringMatchesRingOnly && !mol1Bond->getIsAromatic() &&
        !mol2Bond->getIsAromatic() &&
        mol1Rings->numBondRings(mol1Bond->getIdx()) &&
        mol2Rings->numBondRings(mol2Bond->getIdx())) {
      b.expandQuery(RDKit::makeBondIsInRingQuery(), Queries::COMPOSITE_AND,
                    true);
    }
    smartsMol->addBond(&b, false);
  }
  std::string smt = RDKit::MolToSmarts(*smartsMol, true);
  cleanSmarts(smt);
  return smt;
}

namespace {
// Return the atom common to the two bonds, -1 if there isn't one.
int common_atom_in_bonds(const RDKit::Bond *bond1, const RDKit::Bond *bond2) {
  int commonAtom = -1;
  if (bond1->getBeginAtomIdx() == bond2->getBeginAtomIdx()) {
    commonAtom = bond1->getBeginAtomIdx();
  } else if (bond1->getEndAtomIdx() == bond2->getBeginAtomIdx()) {
    commonAtom = bond1->getEndAtomIdx();
  } else if (bond1->getBeginAtomIdx() == bond2->getEndAtomIdx()) {
    commonAtom = bond1->getBeginAtomIdx();
  } else if (bond1->getEndAtomIdx() == bond2->getEndAtomIdx()) {
    commonAtom = bond1->getEndAtomIdx();
  }
  return commonAtom;
}
}  // namespace

void RascalResult::matchCliqueAtoms(
    const std::vector<std::vector<int>> &mol1_adj_matrix) {
  if (d_bondMatches.empty()) {
    return;
  }
  std::vector<int> mol1Matches(d_mol1->getNumAtoms(), -1);
  // set the clique atoms to -2 in mol1Matches, to mark them as yet undecided.
  for (const auto &bm : d_bondMatches) {
    auto bond1 = d_mol1->getBondWithIdx(bm.first);
    mol1Matches[bond1->getBeginAtomIdx()] = -2;
    mol1Matches[bond1->getEndAtomIdx()] = -2;
  }

  // First, use the line graphs to match atoms that have 2 matching bonds
  // incident on them.
  for (size_t i = 0; i < d_bondMatches.size() - 1; ++i) {
    const auto &pair1 = d_bondMatches[i];
    auto bond1_1 = d_mol1->getBondWithIdx(pair1.first);
    auto bond2_1 = d_mol2->getBondWithIdx(pair1.second);
    for (size_t j = i + 1; j < d_bondMatches.size(); ++j) {
      const auto &pair2 = d_bondMatches[j];
      if (mol1_adj_matrix[pair1.first][pair2.first]) {
        // the 2 bonds are incident on the same atom, so the 2 atoms must match
        auto bond1_2 = d_mol1->getBondWithIdx(pair2.first);
        auto bond2_2 = d_mol2->getBondWithIdx(pair2.second);
        auto mol1Atom = common_atom_in_bonds(bond1_1, bond1_2);
        auto mol2Atom = common_atom_in_bonds(bond2_1, bond2_2);
        if (mol1Atom != -1) {
          mol1Matches[mol1Atom] = mol2Atom;
          auto omol1Atom = bond1_1->getOtherAtomIdx(mol1Atom);
          auto omol2Atom = bond2_1->getOtherAtomIdx(mol2Atom);
          mol1Matches[omol1Atom] = omol2Atom;
          omol1Atom = bond1_2->getOtherAtomIdx(mol1Atom);
          omol2Atom = bond2_2->getOtherAtomIdx(mol2Atom);
          mol1Matches[omol1Atom] = omol2Atom;
        }
      }
    }
  }
  // if there are -2 entries in mol1Matches there's more to do.
  if (std::count(mol1Matches.begin(), mol1Matches.end(), -2)) {
    // Any -2 entries in mol1Matches are down to isolated bonds, which are a bit
    // tricky.
    for (size_t i = 0; i < d_bondMatches.size(); ++i) {
      const auto &pair1 = d_bondMatches[i];
      auto bond1_1 = d_mol1->getBondWithIdx(pair1.first);
      if (mol1Matches[bond1_1->getBeginAtomIdx()] == -2 &&
          mol1Matches[bond1_1->getEndAtomIdx()] == -2) {
        auto bond2_1 = d_mol2->getBondWithIdx(pair1.second);
        if (bond1_1->getBeginAtom()->getAtomicNum() !=
            bond1_1->getEndAtom()->getAtomicNum()) {
          // it's fairly straightforward:
          if (bond1_1->getBeginAtom()->getAtomicNum() ==
              bond2_1->getBeginAtom()->getAtomicNum()) {
            mol1Matches[bond1_1->getBeginAtomIdx()] =
                bond2_1->getBeginAtomIdx();
            mol1Matches[bond1_1->getEndAtomIdx()] = bond2_1->getEndAtomIdx();
          } else {
            mol1Matches[bond1_1->getBeginAtomIdx()] = bond2_1->getEndAtomIdx();
            mol1Matches[bond1_1->getEndAtomIdx()] = bond2_1->getBeginAtomIdx();
          }
        } else if (bond1_1->getBeginAtom()->getTotalNumHs() !=
                   bond1_1->getEndAtom()->getTotalNumHs()) {
          // try it on number of hydrogens
          if (bond1_1->getBeginAtom()->getTotalNumHs() >
              bond1_1->getEndAtom()->getTotalNumHs()) {
            mol1Matches[bond1_1->getBeginAtomIdx()] =
                bond2_1->getBeginAtomIdx();
            mol1Matches[bond1_1->getEndAtomIdx()] = bond2_1->getEndAtomIdx();
          } else {
            mol1Matches[bond1_1->getBeginAtomIdx()] = bond2_1->getEndAtomIdx();
            mol1Matches[bond1_1->getEndAtomIdx()] = bond2_1->getBeginAtomIdx();
          }
        } else {
          // it probably doesn't matter
          mol1Matches[bond1_1->getBeginAtomIdx()] = bond2_1->getBeginAtomIdx();
          mol1Matches[bond1_1->getEndAtomIdx()] = bond2_1->getEndAtomIdx();
        }
      }
    }
  }
  for (size_t i = 0u; i < d_mol1->getNumAtoms(); ++i) {
    if (mol1Matches[i] >= 0) {
      d_atomMatches.push_back(std::make_pair(i, mol1Matches[i]));
    }
  }
}

void RascalResult::applyMaxFragSep() {
  std::unique_ptr<RDKit::ROMol> mol1_frags(makeMolFrags(1));
  auto frags1 = RDKit::MolOps::getMolFrags(*mol1_frags, false);
  if (frags1.size() < 2) {
    return;
  }
  auto fragFragDist = [](const boost::shared_ptr<RDKit::ROMol> &frag1,
                         const boost::shared_ptr<RDKit::ROMol> &frag2,
                         const double *pathMatrix, int num_atoms) -> double {
    int minDist = std::numeric_limits<int>::max();
    for (auto at1 : frag1->atoms()) {
      int at1Idx = at1->getProp<int>("ORIG_INDEX");
      for (auto at2 : frag2->atoms()) {
        int at2Idx = at2->getProp<int>("ORIG_INDEX");
        int dist = std::nearbyint(pathMatrix[at1Idx * num_atoms + at2Idx]);
        if (dist < minDist) {
          minDist = dist;
        }
      }
    }
    return minDist;
  };

  std::unique_ptr<RDKit::ROMol> mol2Frags(makeMolFrags(2));
  auto frags2 = RDKit::MolOps::getMolFrags(*mol2Frags, false);
  // These arrays must not be deleted - they are cached in the molecule and
  // deleted when it is. The distance matrix will be re-calculated in case
  // something's been copied over somewhere.
  auto mol1Dists = RDKit::MolOps::getDistanceMat(*d_mol1, false, false, true);
  auto mol2Dists = RDKit::MolOps::getDistanceMat(*d_mol2, false, false, true);

  bool deletedFrag = false;
  for (size_t i = 0; i < frags1.size() - 1; ++i) {
    if (!frags1[i]) {
      continue;
    }
    for (size_t j = i + 1; j < frags1.size(); ++j) {
      if (!frags1[j]) {
        continue;
      }
      int mol1Dist =
          fragFragDist(frags1[i], frags1[j], mol1Dists, d_mol1->getNumAtoms());
      int mol2Dist =
          fragFragDist(frags2[i], frags2[j], mol2Dists, d_mol2->getNumAtoms());
      if (mol1Dist > d_maxFragSep || mol2Dist > d_maxFragSep) {
        deletedFrag = true;
        if (frags1[i]->getNumAtoms() < frags1[j]->getNumAtoms()) {
          frags1[i].reset();
          frags2[i].reset();
        } else {
          frags1[j].reset();
          frags2[j].reset();
        }
      }
    }
  }

  if (deletedFrag) {
    // rebuild the d_bondMatches
    std::vector<std::pair<int, int>> new_bond_matches;
    for (size_t i = 0; i < frags1.size(); ++i) {
      if (!frags1[i]) {
        continue;
      }
      for (auto b : frags1[i]->bonds()) {
        int b_idx = b->getProp<int>("ORIG_INDEX");
        for (auto &bm : d_bondMatches) {
          if (b_idx == bm.first) {
            new_bond_matches.push_back(bm);
            break;
          }
        }
      }
    }
    d_bondMatches = new_bond_matches;
    // and the d_atomMatches
    std::vector<std::pair<int, int>> new_atom_matches;
    for (size_t i = 0; i < frags1.size(); ++i) {
      if (!frags1[i]) {
        continue;
      }
      for (auto a : frags1[i]->atoms()) {
        int a_idx = a->getProp<int>("ORIG_INDEX");
        for (auto &am : d_atomMatches) {
          if (a_idx == am.first) {
            new_atom_matches.push_back(am);
            break;
          }
        }
      }
    }
    d_atomMatches = new_atom_matches;
  }
}

// Return a molecule with the clique in it.  Each atom will have the property
// ORIG_INDEX giving its index in the original molecule.
RDKit::ROMol *RascalResult::makeMolFrags(int molNum) const {
  std::shared_ptr<RDKit::ROMol> theMol;
  if (molNum == 1) {
    theMol = d_mol1;
  } else if (molNum == 2) {
    theMol = d_mol2;
  } else {
    return nullptr;
  }
  if (!theMol) {
    return nullptr;
  }
  auto *molFrags = new RDKit::RWMol(*theMol);
  std::vector<char> ainClique(theMol->getNumAtoms(), 0);
  for (const auto &am : d_atomMatches) {
    if (molNum == 1) {
      ainClique[am.first] = 1;
    } else {
      ainClique[am.second] = 1;
    }
  }
  std::vector<char> binClique(theMol->getNumBonds(), 0);
  for (const auto &bm : d_bondMatches) {
    if (molNum == 1) {
      binClique[bm.first] = 1;
    } else {
      binClique[bm.second] = 1;
    }
  }
  molFrags->beginBatchEdit();
  for (auto &a : molFrags->atoms()) {
    if (!ainClique[a->getIdx()]) {
      molFrags->removeAtom(a);
    } else {
      a->setProp<int>("ORIG_INDEX", a->getIdx());
    }
  }
  for (auto &b : molFrags->bonds()) {
    if (!binClique[b->getIdx()]) {
      molFrags->removeBond(b->getBeginAtomIdx(), b->getEndAtomIdx());
    } else {
      b->setProp<int>("ORIG_INDEX", b->getIdx());
    }
  }
  molFrags->commitBatchEdit();
  return molFrags;
}

// Calculate a score for how many bonds in the clique don't match
// cyclic/non-cyclic
int RascalResult::calcRingNonRingScore() const {
  if (!d_mol1 || !d_mol2) {
    return 0;
  }

  int score = 0;
  for (const auto &bm : d_bondMatches) {
    auto nbr1 = d_mol1->getRingInfo()->numBondRings(bm.first);
    auto nbr2 = d_mol2->getRingInfo()->numBondRings(bm.second);

    if ((nbr1 && !nbr2) || (!nbr1 && nbr2)) {
      ++score;
    }
  }
  return score;
}

// Calculate a score for how well the atoms in the clique from mol1 match the
// atoms for the clique in mol2.  The atom scores are made up of H count and
// summed for the molecule. Its so that, for example, an OH in mol1 that could
// match an OH or OMe matches the OH for preference.
int RascalResult::calcAtomMatchScore() const {
  if (!d_mol1 || !d_mol2) {
    return 0;
  }
  int score = 0;
  for (const auto &am : d_atomMatches) {
    int num_h_1 = d_mol1->getAtomWithIdx(am.first)->getTotalNumHs();
    int num_h_2 = d_mol2->getAtomWithIdx(am.second)->getTotalNumHs();
    score += std::abs(num_h_1 - num_h_2);
  }
  return score;
}

int RascalResult::calcMaxDeltaAtomAtomDistScore() const {
  // Possibly this could be improved, to be the total of the minimum distances
  // between each fragment.
  if (d_atomMatches.empty()) {
    return 0;
  }
  // These arrays are cached so shouldn't be deleted.  The final 'true' in the
  // call is to force recalculation, just in case there's some other type copied
  // over from the input molecule.
  const auto *mol1Dists =
      RDKit::MolOps::getDistanceMat(*d_mol1, false, false, true);
  const auto *mol2Dists =
      RDKit::MolOps::getDistanceMat(*d_mol2, false, false, true);

  int score = 0;
  auto dist = [](int idx1, int idx2, const double *dists,
                 int num_atoms) -> int {
    return int(std::nearbyint(dists[idx1 * num_atoms + idx2]));
  };
  for (size_t i = 0; i < d_atomMatches.size() - 1; ++i) {
    for (size_t j = i + 1; j < d_atomMatches.size(); ++j) {
      auto d1 = dist(d_atomMatches[i].first, d_atomMatches[j].first, mol1Dists,
                     d_mol1->getNumAtoms());
      auto d2 = dist(d_atomMatches[i].second, d_atomMatches[j].second,
                     mol2Dists, d_mol2->getNumAtoms());
      auto deltaDist = abs(d1 - d2);
      if (deltaDist > score) {
        score = deltaDist;
      }
    }
  }
  return score;
}

int RascalResult::calcLargestFragSize() const {
  if (!d_mol1 || !d_mol2) {
    return 0;
  }
  std::unique_ptr<RDKit::ROMol> mol1_frags(makeMolFrags(1));
  std::vector<int> mapping;
  auto numFrags = RDKit::MolOps::getMolFrags(*mol1_frags, mapping);
  int lfs = -1;
  for (unsigned int i = 0; i < numFrags; ++i) {
    auto fragSize = std::count(mapping.begin(), mapping.end(), i);
    if (fragSize > lfs) {
      lfs = fragSize;
    }
  }
  return lfs;
}

int RascalResult::numFrags() const {
  if (!d_mol1 || !d_mol2) {
    return 0;
  }
  if (d_numFrags == -1) {
    std::unique_ptr<RDKit::ROMol> mol1_frags(makeMolFrags(1));
    std::vector<int> mol1_frag_mapping;
    d_numFrags = RDKit::MolOps::getMolFrags(*mol1_frags, mol1_frag_mapping);
  }
  return d_numFrags;
}

int RascalResult::ringNonRingBondScore() const {
  if (!d_mol1 || !d_mol2) {
    return 0;
  }
  if (d_ringNonRingBondScore == -1) {
    d_ringNonRingBondScore = calcRingNonRingScore();
  }
  return d_ringNonRingBondScore;
}

int RascalResult::atomMatchScore() const {
  if (!d_mol1 || !d_mol2) {
    return 0;
  }
  if (d_atomMatchScore == -1) {
    d_atomMatchScore = calcAtomMatchScore();
  }
  return d_atomMatchScore;
}

int RascalResult::maxDeltaAtomAtomDist() const {
  if (!d_mol1 || !d_mol2) {
    return 0;
  }
  if (d_maxDeltaAtomAtomDist == -1) {
    d_maxDeltaAtomAtomDist = calcMaxDeltaAtomAtomDistScore();
  }
  return d_maxDeltaAtomAtomDist;
}

int RascalResult::largestFragSize() const {
  if (!d_mol1 || !d_mol2) {
    return 0;
  }
  if (d_largestFragSize == -1) {
    d_largestFragSize = calcLargestFragSize();
  }
  return d_largestFragSize;
}

std::string RascalResult::smarts() const {
  if (!d_mol1 || !d_mol2) {
    return "";
  }
  if (d_smarts.empty()) {
    d_smarts = createSmartsString();
  }
  return d_smarts;
}

const std::shared_ptr<ROMol> RascalResult::mcesMol() const {
  if (d_mcesMol || !d_mol1) {
    return d_mcesMol;
  }

  std::set<int> mol1Bonds;
  for (const auto &bm : d_bondMatches) {
    mol1Bonds.insert(bm.first);
  }
  std::set<int> mol1Atoms;
  for (const auto &am : d_atomMatches) {
    mol1Atoms.insert(am.first);
  }
  RWMol tmpMol(*d_mol1);
  MolOps::KekulizeIfPossible(tmpMol);
  tmpMol.beginBatchEdit();
  for (auto &bond : tmpMol.bonds()) {
    if (mol1Bonds.find(bond->getIdx()) == mol1Bonds.end()) {
      auto bo = bond->getBondType();
      bond->getBeginAtom()->setNumExplicitHs(
          bond->getBeginAtom()->getNumExplicitHs() + bo);
      bond->getEndAtom()->setNumExplicitHs(
          bond->getEndAtom()->getNumExplicitHs() + bo);
      tmpMol.removeBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
    }
  }
  for (auto atom : tmpMol.atoms()) {
    if (mol1Atoms.find(atom->getIdx()) == mol1Atoms.end()) {
      tmpMol.removeAtom(atom);
    }
  }
  tmpMol.commitBatchEdit();
  MolOps::removeHs(tmpMol);
  MolOps::sanitizeMol(tmpMol);
  d_mcesMol.reset(new ROMol(tmpMol));
  return d_mcesMol;
}

bool resultSort(const RascalResult &res1, const RascalResult &res2) {
  if (res1.bondMatches().size() == res2.bondMatches().size()) {
    if (res1.numFrags() == res2.numFrags()) {
      if (res1.largestFragSize() == res2.largestFragSize()) {
        if (res1.ringNonRingBondScore() == res2.ringNonRingBondScore()) {
          if (res1.atomMatchScore() == res2.atomMatchScore()) {
            if (res1.maxDeltaAtomAtomDist() == res2.maxDeltaAtomAtomDist()) {
              return res1.smarts() < res2.smarts();
            } else {
              return res1.maxDeltaAtomAtomDist() < res2.maxDeltaAtomAtomDist();
            }
          } else {
            return res1.atomMatchScore() < res2.atomMatchScore();
          }
        } else {
          return res1.ringNonRingBondScore() < res2.ringNonRingBondScore();
        }
      } else {
        return res1.largestFragSize() > res2.largestFragSize();
      }
    } else {
      return res1.numFrags() < res2.numFrags();
    }
  } else {
    return res1.bondMatches().size() > res2.bondMatches().size();
  }
}

void extractClique(const std::vector<unsigned int> &clique,
                   const std::vector<std::pair<int, int>> &vtxPairs,
                   bool swapped,
                   std::vector<std::pair<int, int>> &bondMatches) {
  bondMatches.clear();
  for (auto mem : clique) {
    if (swapped) {
      bondMatches.push_back(
          std::make_pair(vtxPairs[mem].second, vtxPairs[mem].first));
    } else {
      bondMatches.push_back(
          std::make_pair(vtxPairs[mem].first, vtxPairs[mem].second));
    }
  }
  std::sort(bondMatches.begin(), bondMatches.end());
}

void cleanSmarts(std::string &smarts) {
  const static std::vector<std::pair<std::regex, std::string>> repls{
      {std::regex(R"(\[#6&A\])"), "C"},
      {std::regex(R"(\[#6&A&R\])"), "[C&R]"},
      {std::regex(R"(\[#6&a\])"), "c"},
      {std::regex(R"(\[#7&A\])"), "N"},
      {std::regex(R"(\[#7&A&R\])"), "[N&R]"},
      {std::regex(R"(\[#7&a\])"), "n"},
      {std::regex(R"(\[#8&A\])"), "O"},
      {std::regex(R"(\[#8&A&R\])"), "[O&R]"},
      {std::regex(R"(\[#8&a\])"), "o"},
      {std::regex(R"(\[#9&A\])"), "F"},
      {std::regex(R"(\[#16&A\])"), "S"},
      {std::regex(R"(\[#16&a\])"), "s"},
      {std::regex(R"(\[#17&A\])"), "Cl"},
      {std::regex(R"(\[#35&A\])"), "Br"},
      {std::regex(R"(\[#53&A\])"), "I"},
      {std::regex(R"(([cnops][1-9]*):([cnops]))"), "$1$2"},
      {std::regex(R"(([cnops]):([1-9]))"), "$1$2"},
      {std::regex(R"(([A-Z])-([cnops]))"), "$1$2"},
      {std::regex(R"(([cnops][1-9]*)-([A-Z]))"), "$1$2"},
      {std::regex(R"(([A-Z][1-9]*)-([A-Z]))"), "$1$2"},
      {std::regex(R"(([A-Z])-([1-9]))"), "$1$2"}};
  // Sometimes it needs more than 1 pass through
  for (int i = 0; i < 10; ++i) {
    std::string start_smt = smarts;
    for (const auto &repl : repls) {
      smarts = std::regex_replace(smarts, std::get<0>(repl), std::get<1>(repl));
    }
    if (start_smt == smarts) {
      break;
    }
  }
}

void printBondMatches(const RascalResult &res, std::ostream &os) {
  os << "Bond 1 matches : " << res.bondMatches().size() << " : [";
  for (const auto &bm : res.bondMatches()) {
    os << bm.first << ",";
  }
  os << "]" << std::endl;
  os << "Bond 2 matches : " << res.bondMatches().size() << " : [";
  for (const auto &bm : res.bondMatches()) {
    os << bm.second << ",";
  }
  os << "]" << std::endl;
}

void printAtomMatches(const RascalResult &res, std::ostream &os) {
  os << "Atom 1 matches : " << res.atomMatches().size() << " : [";
  for (const auto &am : res.atomMatches()) {
    os << am.first << ",";
  }
  os << "]" << std::endl;
  os << "Atom 2 matches : " << res.atomMatches().size() << " : [";
  for (const auto &am : res.atomMatches()) {
    os << am.second << ",";
  }
  os << "]" << std::endl;
}

void printScores(const RascalResult &res, std::ostream &os) {
  os << res.bondMatches().size() << " : " << res.numFrags() << " : "
     << res.largestFragSize() << " : " << res.ringNonRingBondScore() << " : "
     << res.atomMatchScore() << " : " << res.maxDeltaAtomAtomDist() << " : "
     << res.smarts() << std::endl;
}

double johnsonSimilarity(const std::vector<std::pair<int, int>> &bondMatches,
                         const std::vector<std::pair<int, int>> &atomMatches,
                         const RDKit::ROMol &mol1, const RDKit::ROMol &mol2) {
  double num = (bondMatches.size() + atomMatches.size()) *
               (bondMatches.size() + atomMatches.size());
  double denom = (mol1.getNumAtoms() + mol1.getNumBonds()) *
                 (mol2.getNumAtoms() + mol2.getNumBonds());
  return num / denom;
}

}  // namespace RascalMCES
}  // namespace RDKit
