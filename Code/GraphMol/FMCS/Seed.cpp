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

unsigned int Seed::addAtom(const Atom* atom) {
  unsigned int i = MoleculeFragment.AtomsIdx.size();
  unsigned int aqi = atom->getIdx();
  MoleculeFragment.Atoms.push_back(atom);
  MoleculeFragment.AtomsIdx.push_back(aqi);
  MoleculeFragment.SeedAtomIdxMap[aqi] = i;
  Topology.addAtom(aqi);
#ifdef DUP_SUBSTRUCT_CACHE
  DupCacheKey.addAtom(aqi);
#endif
  return i;
}

unsigned int Seed::addBond(const Bond* bond) {
  unsigned int b = bond->getIdx();
  CHECK_INVARIANT(!ExcludedBonds[b], "");
  ExcludedBonds[b] = true;
  MoleculeFragment.BondsIdx.push_back(b);
  MoleculeFragment.Bonds.push_back(bond);
  // remap idx to seed's indices:
  unsigned int i = MoleculeFragment.SeedAtomIdxMap[bond->getBeginAtomIdx()];
  unsigned int j = MoleculeFragment.SeedAtomIdxMap[bond->getEndAtomIdx()];
  Topology.addBond(b, i, j);
#ifdef DUP_SUBSTRUCT_CACHE
  DupCacheKey.addBond(b);
#endif
  return getNumBonds();
}

void Seed::fillNewBonds(const ROMol& qmol) {
  auto excludedBonds = ExcludedBonds;
  // all atoms added on previous growing only
  for (unsigned int srcAtomIdx = LastAddedAtomsBeginIdx;
       srcAtomIdx < getNumAtoms(); srcAtomIdx++) {
    const auto atom = MoleculeFragment.Atoms[srcAtomIdx];
    for (const auto& nbri :
         boost::make_iterator_range(qmol.getAtomBonds(atom))) {
      const auto bond = qmol[nbri];
      if (!excludedBonds.test(bond->getIdx())) {
        // already in the seed or NewBonds list from another atom in a RING
        excludedBonds.set(bond->getIdx());
        unsigned int ai = (atom == bond->getBeginAtom())
                              ? bond->getEndAtomIdx()
                              : bond->getBeginAtomIdx();
        const auto end_atom = qmol.getAtomWithIdx(ai);
        unsigned int end_atom_idx = NotSet;
        for (unsigned int i = 0; i < getNumAtoms(); i++) {
          // already exists in this seed
          if (end_atom == MoleculeFragment.Atoms[i]) {
            end_atom_idx = i;
            break;
          }
        }
        NewBonds.emplace_back(srcAtomIdx, bond->getIdx(), ai, end_atom_idx,
                              NotSet == end_atom_idx ? end_atom : nullptr);
      }
    }
  }
}

void Seed::grow(MaximumCommonSubgraph& mcs) const {
  const auto& qmol = mcs.getQueryMolecule();
  std::set<unsigned int> newAtomsSet;  // keep track of newly added atoms

  if (!canGrowBiggerThan(mcs.getMaxNumberBonds(),
                         mcs.getMaxNumberAtoms())) {  // prune() parent
    GrowingStage = NotSet;                            // finished
#ifdef VERBOSE_STATISTICS_ON
    ++mcs.VerboseStatistics.RemainingSizeRejected;
#endif
    return;
  }

  if (0 == GrowingStage) {
    // 0. Fill out list of all directly connected outgoing bonds
    // non const method, multistage growing optimisation
    const_cast<Seed*>(this)->fillNewBonds(qmol);

    if (NewBonds.empty()) {
      GrowingStage = NotSet;  // finished
      return;
    }

    // 1. Check and add the biggest child seed with all outgoing bonds added:
    // Add all bonds at first (build the biggest child seed). All new atoms are
    // already in the seed
    Seed seed;
    seed.createFromParent(this);

    for (const auto& newBond : NewBonds) {
      unsigned int aIdx = newBond.EndAtomIdx;
      if (NotSet == aIdx) {  // new atom
        // check if new bonds simultaneously close a ring
        if (newAtomsSet.find(newBond.NewAtomIdx) == newAtomsSet.end()) {
          const auto end_atom = newBond.NewAtom;
          aIdx = seed.addAtom(end_atom);
          newAtomsSet.insert(newBond.NewAtomIdx);
        }
      }
      const auto src_bond = qmol.getBondWithIdx(newBond.BondIdx);
      seed.addBond(src_bond);
    }
#ifdef VERBOSE_STATISTICS_ON
    { ++mcs.VerboseStatistics.Seed; }
#endif
    seed.RemainingBonds = RemainingBonds - NewBonds.size();  // Added ALL !!!
    seed.RemainingAtoms =
        RemainingAtoms - newAtomsSet.size();  // new atoms added to seed

    // prune() Best Sizes
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
  // check each single bond first: if (this seed + single bond) does not exist
  // in MCS, exclude this new bond from growing this seed.
  unsigned int numErasedNewBonds = 0;
  for (auto& nbi : NewBonds) {
#ifdef VERBOSE_STATISTICS_ON
    { ++mcs.VerboseStatistics.Seed; }
#endif
    Seed seed;
    seed.createFromParent(this);

    // existed in this parent seed (ring) or -1
    unsigned int aIdx = nbi.EndAtomIdx;
    if (NotSet == aIdx) {  // new atom
      const auto end_atom = nbi.NewAtom;
      aIdx = seed.addAtom(end_atom);
    }
    const auto src_bond = qmol.getBondWithIdx(nbi.BondIdx);
    seed.addBond(src_bond);
    seed.computeRemainingSize(qmol);

    if (seed.canGrowBiggerThan(mcs.getMaxNumberBonds(),
                               mcs.getMaxNumberAtoms())) {  // prune()
      if (!MatchResult.empty()) {
        seed.MatchResult = MatchResult;
      }
      if (!mcs.checkIfMatchAndAppend(seed)) {
        nbi.BondIdx = NotSet;  // exclude this new bond from growing this seed
                               // - decrease 2^^N-1 to 2^^k-1, k<N.
        ++numErasedNewBonds;
#ifdef VERBOSE_STATISTICS_ON
        ++mcs.VerboseStatistics.SingleBondExcluded;
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
    for (const auto& dirtyNewBond : dirtyNewBonds) {
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
      // exclude already processed single external bond combinations
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
      newAtomsSet.clear();

      for (unsigned int i = 0; i < NewBonds.size(); i++) {
        if (composition.isSet(i)) {
          const auto nbi = &NewBonds[i];
          // existed in this parent seed (ring) or -1
          unsigned int aIdx = nbi->EndAtomIdx;
          if (NotSet == aIdx) {  // new atom
            if (newAtomsSet.find(nbi->NewAtomIdx) == newAtomsSet.end()) {
              const auto end_atom = nbi->NewAtom;
              aIdx = seed.addAtom(end_atom);
              newAtomsSet.insert(nbi->NewAtomIdx);
            }
          }
          const auto src_bond = qmol.getBondWithIdx(nbi->BondIdx);
          seed.addBond(src_bond);
        }
      }
      seed.computeRemainingSize(qmol);
      if (!seed.canGrowBiggerThan(
              mcs.getMaxNumberBonds(),
              mcs.getMaxNumberAtoms())) {  // prune(). // seed too small
#ifdef VERBOSE_STATISTICS_ON
        ++mcs.VerboseStatistics.RemainingSizeRejected;
#endif
      } else {
        seed.MatchResult = MatchResult;
        bool found = mcs.checkIfMatchAndAppend(seed);

        if (!found) {
#ifdef EXCLUDE_WRONG_COMPOSITION  // if seed does not matched it is possible to
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

void Seed::computeRemainingSize(const ROMol& qmol) {
  RemainingBonds = RemainingAtoms = 0;

  std::vector<unsigned int> end_atom_stack;
  auto visitedBonds = ExcludedBonds;
  boost::dynamic_bitset<> visitedAtoms(qmol.getNumAtoms());

  std::for_each(MoleculeFragment.AtomsIdx.begin(),
                MoleculeFragment.AtomsIdx.end(),
                [&visitedAtoms](unsigned int i) { visitedAtoms.set(i); });

  // SDF all paths
  // 1. direct neighbours
  for (unsigned int seedAtomIdx = LastAddedAtomsBeginIdx;
       seedAtomIdx < getNumAtoms();
       seedAtomIdx++) {  // just now added new border vertices (candidates for
                         // future growing)
    const auto atom = MoleculeFragment.Atoms[seedAtomIdx];
    for (const auto& nbri :
         boost::make_iterator_range(qmol.getAtomBonds(atom))) {
      const auto bond = qmol[nbri];
      if (!visitedBonds.test(bond->getIdx())) {
        ++RemainingBonds;
        visitedBonds.set(bond->getIdx());
        unsigned int end_atom_idx =
            (MoleculeFragment.AtomsIdx[seedAtomIdx] == bond->getBeginAtomIdx())
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
    for (const auto& nbri :
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
