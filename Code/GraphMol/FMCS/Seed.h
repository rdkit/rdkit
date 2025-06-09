//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#pragma once
#include <map>
#include <boost/dynamic_bitset.hpp>
#include "../RDKitBase.h"
// algorithm optimisation definitions
#include "DebugTrace.h"
#include "Graph.h"
#include "DuplicatedSeedCache.h"
#include "SubstructMatchCustom.h"
#include "TargetMatch.h"

namespace RDKit {
namespace FMCS {
class MaximumCommonSubgraph;
struct TargetMatch;

// Reference to a fragment of source molecule
struct RDKIT_FMCS_EXPORT MolFragment {
  std::vector<const Atom *> Atoms;
  std::vector<const Bond *> Bonds;
  // Full Query Molecule to Seed indices backward conversionmap
  std::map<unsigned int, unsigned int> SeedAtomIdxMap;
};

struct RDKIT_FMCS_EXPORT NewBond {
  // index in qmol of new bond scheduled to be added into
  // seed. This is outgoing bond from SourceAtomIdx
  unsigned int BondIdx{0};
  // index in qmol of new atom scheduled to be
  // added into seed. Another end of new bond
  unsigned int NewAtomIdx{0};
  // index in the seed. RING. "New" Atom on the another
  // end of new bond if it already exists in the seed.
  unsigned int EndAtomIdx{0};
  // pointer to qmol's new atom scheduled to be
  // added into seed. Another end of new bond
  const Atom *NewAtom{nullptr};

  NewBond()

  {}

  NewBond(unsigned int bond_idx, unsigned int new_atom, unsigned int to_atom,
          const Atom *a)
      : BondIdx(bond_idx),
        NewAtomIdx(new_atom),
        EndAtomIdx(to_atom),
        NewAtom(a) {}
};

class RDKIT_FMCS_EXPORT Seed {
 private:
  boost::dynamic_bitset<> addNewBondsToSeed(const ROMol &qmol,
                                            Seed &seed) const;
  bool canAddAllNonFusedRingBondsConnectedToBond(
      const Atom &srcAtom, const Bond &bond, MaximumCommonSubgraph &mcs) const;
  void addNewBondFromAtom(const Atom &srcAtom, const Bond &bond) const;
  // for multistage growing. all directly connected outgoing bonds
  mutable std::vector<NewBond> NewBonds;
  bool StoreAllDegenerateMCS = false;

 public:
  // this seed has been completely copied into list.
  // postponed non-locked copy for MULTI_THREAD
  bool CopyComplete{false};
  // 0 new seed; -1 finished; n>0 in
  // progress, exact stage of growing for SDF
  mutable unsigned int GrowingStage{0};
  // Reference to a fragment of source molecule
  MolFragment MoleculeFragment;
  // seed topology with references to source molecule
  Graph Topology;

  boost::dynamic_bitset<> ExcludedBonds;
  // in this subgraph for improving performance of future growing
  unsigned int LastAddedAtomsBeginIdx{0};
  // in this subgraph for DEBUG ONLY
  unsigned int LastAddedBondsBeginIdx{0};
  unsigned int RemainingBonds{0};
  unsigned int RemainingAtoms{0};
#ifdef DUP_SUBSTRUCT_CACHE
  DuplicatedSeedCache::TKey DupCacheKey;
#endif
  // for each target
  std::vector<TargetMatch> MatchResult;

 public:
  Seed()

  {}

  void setMoleculeFragment(const Seed &src) {
    MoleculeFragment = src.MoleculeFragment;
  }
  Seed &operator=(const Seed &src) {
    NewBonds = src.NewBonds;
    GrowingStage = src.GrowingStage;
    MoleculeFragment = src.MoleculeFragment;
    Topology = src.Topology;
    ExcludedBonds = src.ExcludedBonds;
    LastAddedAtomsBeginIdx = src.LastAddedAtomsBeginIdx;
    LastAddedBondsBeginIdx = src.LastAddedBondsBeginIdx;
    RemainingBonds = src.RemainingBonds;
    RemainingAtoms = src.RemainingAtoms;
    StoreAllDegenerateMCS = src.StoreAllDegenerateMCS;
#ifdef DUP_SUBSTRUCT_CACHE
    DupCacheKey = src.DupCacheKey;
#endif
    MatchResult = src.MatchResult;
    CopyComplete = true;  // LAST
    return *this;
  }
  void createFromParent(const Seed *parent) {
    MoleculeFragment = parent->MoleculeFragment;
    Topology = parent->Topology;
    ExcludedBonds = parent->ExcludedBonds;
    RemainingBonds = parent->RemainingBonds;
    RemainingAtoms = parent->RemainingAtoms;
    StoreAllDegenerateMCS = parent->StoreAllDegenerateMCS;
#ifdef DUP_SUBSTRUCT_CACHE
    DupCacheKey = parent->DupCacheKey;
#endif
    LastAddedAtomsBeginIdx = getNumAtoms();  // previous size
    LastAddedBondsBeginIdx = getNumBonds();  // previous size
    GrowingStage = 0;
  }

  unsigned int getNumAtoms() const { return MoleculeFragment.Atoms.size(); }
  unsigned int getNumBonds() const { return MoleculeFragment.Bonds.size(); }

  void grow(MaximumCommonSubgraph &mcs) const;
  bool canGrowBiggerThan(unsigned int maxBonds, unsigned int maxAtoms) const {
    return RemainingBonds + getNumBonds() > maxBonds ||
           (RemainingBonds + getNumBonds() == maxBonds &&
            (RemainingAtoms + getNumAtoms() > maxAtoms ||
             (StoreAllDegenerateMCS &&
              RemainingAtoms + getNumAtoms() == maxAtoms)));
  }
  void computeRemainingSize(const ROMol &qmol);

  unsigned int addAtom(const Atom *atom);
  unsigned int addBond(const Bond *bond);
  void fillNewBonds(const ROMol &qmol,
                    MaximumCommonSubgraph *mcs = nullptr) const;
  void setStoreAllDegenerateMCS(bool value) { StoreAllDegenerateMCS = value; }
};
}  // namespace FMCS
}  // namespace RDKit
