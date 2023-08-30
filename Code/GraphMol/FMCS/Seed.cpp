//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MaximumCommonSubgraph.h"
#include "Composition2N.h"
#include "Seed.h"

#include "DebugTrace.h"
#include "../SmilesParse/SmilesWrite.h"

#include <set>

namespace RDKit {
namespace FMCS {

unsigned int Seed::addAtom(const Atom *atom) {
  unsigned int i = MoleculeFragment.Atoms.size();
  unsigned int aqi = atom->getIdx();
  MoleculeFragment.Atoms.push_back(atom);
  MoleculeFragment.SeedAtomIdxMap[aqi] = i;
  Topology.addAtom(aqi);
#ifdef DUP_SUBSTRUCT_CACHE
  DupCacheKey.addAtom(aqi);
#endif
  return i;
}

unsigned int Seed::addBond(const Bond *bond) {
  unsigned int bi = bond->getIdx();
  CHECK_INVARIANT(!ExcludedBonds.test(bi), "");
  ExcludedBonds.set(bi);
  MoleculeFragment.Bonds.push_back(bond);
  // remap idx to seed's indices:
  unsigned int i = MoleculeFragment.SeedAtomIdxMap.at(bond->getBeginAtomIdx());
  unsigned int j = MoleculeFragment.SeedAtomIdxMap.at(bond->getEndAtomIdx());
  Topology.addBond(bi, i, j);
#ifdef DUP_SUBSTRUCT_CACHE
  DupCacheKey.addBond(bi);
#endif
  return getNumBonds();
}

boost::dynamic_bitset<> Seed::addNewBondsToSeed(const ROMol &qmol,
                                                Seed &seed) const {
  boost::dynamic_bitset<> newAtomsSet(qmol.getNumAtoms());
  for (const auto &newBond : NewBonds) {
    unsigned int aIdx = newBond.EndAtomIdx;
    if (NotSet == aIdx) {  // new atom
      // check if new bonds simultaneously close a ring
      if (!newAtomsSet.test(newBond.NewAtomIdx)) {
        const auto end_atom = newBond.NewAtom;
        aIdx = seed.addAtom(end_atom);
        newAtomsSet.set(newBond.NewAtomIdx);
      }
    }
    const auto src_bond = qmol.getBondWithIdx(newBond.BondIdx);
    seed.addBond(src_bond);
  }
  seed.RemainingBonds = RemainingBonds - NewBonds.size();  // Added ALL !!!
  seed.RemainingAtoms =
      RemainingAtoms - newAtomsSet.count();  // new atoms added to seed
  return newAtomsSet;
}

bool Seed::canAddAllNonFusedRingBondsConnectedToBond(
    const Atom &srcAtom, const Bond &bond, MaximumCommonSubgraph &mcs) const {
  const auto &mol = bond.getOwningMol();
  const auto ri = mol.getRingInfo();
  int bondIdx = bond.getIdx();
  const auto &bondRings = ri->bondRings().at(ri->bondMembers(bondIdx).front());
  std::set<unsigned int> nonFusedRingBondIndices;
  boost::dynamic_bitset<> connectedAtomIndices(mol.getNumAtoms());
  Seed seed;
  seed.createFromParent(this);
  for (const auto &bi : bondRings) {
    if (bi != bondIdx && ri->numBondRings(bi) == 1) {
      nonFusedRingBondIndices.insert(bi);
    }
  }
  auto currAtom = &srcAtom;
  auto currBond = &bond;
  auto excludedBonds = seed.ExcludedBonds;
  while (currBond) {
    if (!seed.ExcludedBonds.test(currBond->getIdx())) {
      connectedAtomIndices.set(currBond->getBeginAtomIdx());
      connectedAtomIndices.set(currBond->getEndAtomIdx());
      seed.addNewBondFromAtom(*currAtom, *currBond);
    }
    currBond = nullptr;
    for (const auto &candBondIdx : nonFusedRingBondIndices) {
      const auto candBond = mol.getBondWithIdx(candBondIdx);
      if (connectedAtomIndices.test(candBond->getBeginAtomIdx())) {
        currAtom = candBond->getBeginAtom();
      } else if (connectedAtomIndices.test(candBond->getEndAtomIdx())) {
        currAtom = candBond->getEndAtom();
      } else {
        continue;
      }
      nonFusedRingBondIndices.erase(candBondIdx);
      currBond = candBond;
      break;
    }
  }
  if (seed.NewBonds.empty()) {
    return false;
  }
  seed.addNewBondsToSeed(mol, seed);
  seed.MatchResult = MatchResult;
  seed.ExcludedBonds = excludedBonds;
  MCSParameters &p = mcs.parameters();
  bool origMatchFusedRings = p.BondCompareParameters.MatchFusedRings;
  bool origMatchFusedRingsStrict =
      p.BondCompareParameters.MatchFusedRingsStrict;
  p.BondCompareParameters.MatchFusedRings = false;
  p.BondCompareParameters.MatchFusedRingsStrict = false;
  bool res = mcs.match(seed);
  p.BondCompareParameters.MatchFusedRings = origMatchFusedRings;
  p.BondCompareParameters.MatchFusedRingsStrict = origMatchFusedRingsStrict;
  return res;
}

void Seed::addNewBondFromAtom(const Atom &srcAtom, const Bond &bond) const {
  const auto end_atom = bond.getOtherAtom(&srcAtom);
  unsigned int end_atom_idx = NotSet;
  for (unsigned int i = 0; i < getNumAtoms(); ++i) {
    // already exists in this seed
    if (end_atom == MoleculeFragment.Atoms.at(i)) {
      end_atom_idx = i;
      break;
    }
  }
  NewBonds.emplace_back(bond.getIdx(), end_atom->getIdx(), end_atom_idx,
                        NotSet == end_atom_idx ? end_atom : nullptr);
}

void Seed::fillNewBonds(const ROMol &qmol, MaximumCommonSubgraph *mcs) const {
  auto excludedBonds = ExcludedBonds;
  const auto ri = qmol.getRingInfo();
  // all atoms added on previous growing only
  for (unsigned int srcAtomIdx = LastAddedAtomsBeginIdx;
       srcAtomIdx < getNumAtoms(); ++srcAtomIdx) {
    const auto atom = MoleculeFragment.Atoms.at(srcAtomIdx);
    for (const auto &nbri :
         boost::make_iterator_range(qmol.getAtomBonds(atom))) {
      const auto bond = qmol[nbri];
      const auto bi = bond->getIdx();
      // already in the seed or NewBonds list from another atom in a RING
      if (excludedBonds.test(bi)) {
        continue;
      }
      excludedBonds.set(bi);
      if (mcs && mcs->parameters().BondCompareParameters.CompleteRingsOnly &&
          ri->numBondRings(bi) == 1 &&
          !canAddAllNonFusedRingBondsConnectedToBond(*atom, *bond, *mcs)) {
        continue;
      }
      addNewBondFromAtom(*atom, *bond);
    }
  }
}

void Seed::grow(MaximumCommonSubgraph &mcs) const {
  if (!canGrowBiggerThan(mcs.getMaxNumberBonds(), mcs.getMaxNumberAtoms())) {
    GrowingStage = NotSet;  // finished
#ifdef VERBOSE_STATISTICS_ON
    ++mcs.VerboseStatistics.RemainingSizeRejected;
#endif
    return;
  }

  const auto &qmol = mcs.getQueryMolecule();
  boost::dynamic_bitset<> newAtomsSet(
      qmol.getNumAtoms());  // keep track of newly added atoms

  if (0 == GrowingStage) {
    // 0. Fill out list of all directly connected outgoing bonds
    // non const method, multistage growing optimisation
    fillNewBonds(qmol, &mcs);
    if (NewBonds.empty()) {
      GrowingStage = NotSet;  // finished
      return;
    }
    // 1. Check and add the biggest child seed with all outgoing bonds added:
    // Add all bonds at first (build the biggest child seed). All new atoms are
    // already in the seed
    Seed seed;
    seed.createFromParent(this);
    newAtomsSet = addNewBondsToSeed(qmol, seed);
#ifdef VERBOSE_STATISTICS_ON
    ++mcs.VerboseStatistics.Seed;
#endif
    if (!seed.canGrowBiggerThan(mcs.getMaxNumberBonds(),
                                mcs.getMaxNumberAtoms())) {
      GrowingStage = NotSet;
#ifdef VERBOSE_STATISTICS_ON
      ++mcs.VerboseStatistics.RemainingSizeRejected;
#endif
      return;  // the biggest possible subgraph from this seed is too small for
               // future growing. So, skip ALL children !
    }
    seed.MatchResult = MatchResult;
    // this seed + all extern bonds is a part of MCS
    bool allMatched = mcs.checkIfMatchAndAppend(seed);
    GrowingStage = 1;
    if (allMatched && NewBonds.size() > 1) {
      return;  // grow deep first. postpone next growing steps
    }
  }
  // 2. Check and add all 2^N-1-1 other possible seeds:
  if (1 == NewBonds.size()) {
    GrowingStage = NotSet;
    return;  // everything has been done
  }
  // OPTIMISATION:
  // check each individual bond first: if (this seed + individual bond) does not
  // exist in MCS, exclude this new bond from growing this seed.
  unsigned int numErasedNewBonds = 0;
  for (auto &newBond : NewBonds) {
#ifdef VERBOSE_STATISTICS_ON
    { ++mcs.VerboseStatistics.Seed; }
#endif
    Seed seed;
    seed.createFromParent(this);

    // existing in this parent seed (ring) or -1
    unsigned int aIdx = newBond.EndAtomIdx;
    if (NotSet == aIdx) {  // new atom
      const auto end_atom = newBond.NewAtom;
      aIdx = seed.addAtom(end_atom);
    }
    const auto src_bond = qmol.getBondWithIdx(newBond.BondIdx);
    seed.addBond(src_bond);
    seed.computeRemainingSize(qmol);

    if (seed.canGrowBiggerThan(mcs.getMaxNumberBonds(),
                               mcs.getMaxNumberAtoms())) {
      if (!MatchResult.empty()) {
        seed.MatchResult = MatchResult;
      }
      if (!mcs.checkIfMatchAndAppend(seed)) {
        // exclude this new bond from growing this seed
        // - decrease 2^^N-1 to 2^^k-1, k<N.
        newBond.BondIdx = NotSet;
        ++numErasedNewBonds;
#ifdef VERBOSE_STATISTICS_ON
        ++mcs.VerboseStatistics.IndividualBondExcluded;
#endif
      }
    } else {  // seed too small
#ifdef VERBOSE_STATISTICS_ON
      ++mcs.VerboseStatistics.RemainingSizeRejected;
#endif
    }
  }

  if (numErasedNewBonds > 0) {
    std::vector<NewBond> dirtyNewBonds;
    dirtyNewBonds.reserve(NewBonds.size());
    dirtyNewBonds.swap(NewBonds);
    for (const auto &dirtyNewBond : dirtyNewBonds) {
      if (NotSet != dirtyNewBond.BondIdx) {
        NewBonds.push_back(dirtyNewBond);
      }
    }
  }
  // add all other from 2^k-1 possible seeds, where k=newBonds.size()
  // if just one new bond, then seed has already been created
  if (NewBonds.size() > 1) {
    if (sizeof(unsigned long long) * 8 < NewBonds.size()) {
      throw std::runtime_error(
          "Max number of new external bonds of a seed >64");
    }
    BitSet maxCompositionValue;
    Composition2N::compute2N(NewBonds.size(), maxCompositionValue);
    --maxCompositionValue;  // 2^N-1
    Composition2N composition(maxCompositionValue, maxCompositionValue);

#ifdef EXCLUDE_WRONG_COMPOSITION
    std::vector<BitSet> failedCombinations;
    BitSet failedCombinationsMask = 0uLL;
#endif
    while (composition.generateNext()) {
      // exclude already processed individual external bond combinations
      if (composition.is2Power()) {
        continue;
      }
      if (0 == numErasedNewBonds &&
          composition.getBitSet() == maxCompositionValue) {
        continue;  // exclude already processed all external bonds combination
      }
      // 2N-1
#ifdef EXCLUDE_WRONG_COMPOSITION
      // OPTIMISATION. reduce amount of generated seeds and match calls
      // 2120 instead of 2208 match calls on small test. 43 wrongComp-s, 83
      // rejected
      if (failedCombinationsMask & composition.getBitSet()) {
        // possibly exists in the list
        bool compositionWrong = false;
        for (std::vector<BitSet>::const_iterator failed =
                 failedCombinations.begin();
             failed != failedCombinations.end(); failed++)
          if (*failed == (*failed & composition.getBitSet())) {
            // combination includes failed combination
            compositionWrong = true;
            break;
          }
        if (compositionWrong) {
#ifdef VERBOSE_STATISTICS_ON
          ++mcs.VerboseStatistics.WrongCompositionRejected;
#endif
          continue;
        }
      }
#endif
#ifdef VERBOSE_STATISTICS_ON
      { ++mcs.VerboseStatistics.Seed; }
#endif
      Seed seed;
      seed.createFromParent(this);
      newAtomsSet.reset();

      for (const auto &newBond : NewBonds) {
        const auto i = &newBond - &NewBonds.front();
        if (composition.isSet(i)) {
          // existing in this parent seed (ring) or -1
          unsigned int aIdx = newBond.EndAtomIdx;
          if (NotSet == aIdx) {  // new atom
            if (!newAtomsSet.test(newBond.NewAtomIdx)) {
              const auto end_atom = newBond.NewAtom;
              aIdx = seed.addAtom(end_atom);
              newAtomsSet.set(newBond.NewAtomIdx);
            }
          }
          const auto src_bond = qmol.getBondWithIdx(newBond.BondIdx);
          seed.addBond(src_bond);
        }
      }
      seed.computeRemainingSize(qmol);
      if (!seed.canGrowBiggerThan(mcs.getMaxNumberBonds(),
                                  mcs.getMaxNumberAtoms())) {  // seed too small
#ifdef VERBOSE_STATISTICS_ON
        ++mcs.VerboseStatistics.RemainingSizeRejected;
#endif
      } else {
        seed.MatchResult = MatchResult;
        bool found = mcs.checkIfMatchAndAppend(seed);

        if (!found) {
#ifdef EXCLUDE_WRONG_COMPOSITION  // if seed does not match it is possible to
                                  // exclude this mismatched combination for
                                  // performance improvement
          failedCombinations.push_back(composition.getBitSet());
          failedCombinationsMask &= composition.getBitSet();
#ifdef VERBOSE_STATISTICS_ON
          ++mcs.VerboseStatistics.WrongCompositionDetected;
#endif
#endif
        }
      }
    }
  }
  GrowingStage = NotSet;  // finished
}

void Seed::computeRemainingSize(const ROMol &qmol) {
  RemainingBonds = RemainingAtoms = 0;

  std::vector<unsigned int> end_atom_stack;
  auto visitedBonds = ExcludedBonds;
  boost::dynamic_bitset<> visitedAtoms(qmol.getNumAtoms());

  std::for_each(
      MoleculeFragment.Atoms.begin(), MoleculeFragment.Atoms.end(),
      [&visitedAtoms](const auto &atom) { visitedAtoms.set(atom->getIdx()); });

  // SDF all paths
  // 1. direct neighbours
  for (unsigned int seedAtomIdx = LastAddedAtomsBeginIdx;
       seedAtomIdx < getNumAtoms(); ++seedAtomIdx) {
    // just now added new border vertices (candidates for
    // future growing)
    const auto atom = MoleculeFragment.Atoms.at(seedAtomIdx);
    for (const auto &nbri :
         boost::make_iterator_range(qmol.getAtomBonds(atom))) {
      const auto bond = qmol[nbri];
      if (!visitedBonds.test(bond->getIdx())) {
        ++RemainingBonds;
        visitedBonds.set(bond->getIdx());
        unsigned int end_atom_idx = (atom == bond->getBeginAtom())
                                        ? bond->getEndAtomIdx()
                                        : bond->getBeginAtomIdx();
        if (!visitedAtoms.test(end_atom_idx)) {  // check RING/CYCLE
          ++RemainingAtoms;
          visitedAtoms.set(end_atom_idx);
          end_atom_stack.push_back(end_atom_idx);
        }
      }
    }
  }
  // 2. go deep
  while (!end_atom_stack.empty()) {
    unsigned int ai = end_atom_stack.back();
    end_atom_stack.pop_back();
    const auto atom = qmol.getAtomWithIdx(ai);
    for (const auto &nbri :
         boost::make_iterator_range(qmol.getAtomBonds(atom))) {
      const auto bond = qmol[nbri];
      if (!visitedBonds.test(bond->getIdx())) {
        ++RemainingBonds;
        visitedBonds.set(bond->getIdx());
        unsigned int end_atom_idx = (ai == bond->getBeginAtomIdx())
                                        ? bond->getEndAtomIdx()
                                        : bond->getBeginAtomIdx();
        if (!visitedAtoms.test(end_atom_idx)) {  // check RING/CYCLE
          ++RemainingAtoms;
          visitedAtoms.set(end_atom_idx);
          end_atom_stack.push_back(end_atom_idx);
        }
      }
    }
  }
}
}  // namespace FMCS
}  // namespace RDKit
