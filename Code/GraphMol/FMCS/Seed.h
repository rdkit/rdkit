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
#include "../RDKitBase.h"
#include "DebugTrace.h"  // algorithm optimisation definitions
#include "Graph.h"
#include "DuplicatedSeedCache.h"
#include "SubstructMatchCustom.h"

namespace RDKit {
namespace FMCS {
class MaximumCommonSubgraph;
struct TargetMatch;

struct RDKIT_FMCS_EXPORT
    MolFragment {  // Reference to a fragment of source molecule
  std::vector<const Atom*> Atoms;
  std::vector<const Bond*> Bonds;
  std::vector<unsigned> AtomsIdx;
  std::vector<unsigned> BondsIdx;  // need for results and size() only !
  std::map<unsigned, unsigned> SeedAtomIdxMap;  // Full Query Molecule to Seed
                                                // indices backward conversion
                                                // map
};

struct RDKIT_FMCS_EXPORT NewBond {
  unsigned SourceAtomIdx{0};  // index in the seed. Atom is already in the seed
  unsigned BondIdx{0};  // index in qmol of new bond scheduled to be added into
                        // seed. This is outgoing bond from SourceAtomIdx
  unsigned NewAtomIdx{0};  // index in qmol of new atom scheduled to be added
                           // into seed. Another end of new bond
  const Atom* NewAtom{nullptr};  // pointer to qmol's new atom scheduled to be
                                 // added into seed. Another end of new bond
  unsigned EndAtomIdx{0};  // index in the seed. RING. "New" Atom on the another
                           // end of new bond is already exists in the seed.

  NewBond()

  {}

  NewBond(unsigned from_atom, unsigned bond_idx, unsigned new_atom,
          unsigned to_atom, const Atom* a)
      : SourceAtomIdx(from_atom),
        BondIdx(bond_idx),
        NewAtomIdx(new_atom),
        NewAtom(a),
        EndAtomIdx(to_atom) {}
};

class RDKIT_FMCS_EXPORT Seed {
 private:
  mutable std::vector<NewBond> NewBonds;  // for multistage growing. all
                                          // directly connected outgoing bonds
 public:
  bool CopyComplete{false};  // this seed has been completely copied into list.
                             // postponed non-locked copy for MULTI_THREAD
  mutable unsigned GrowingStage{0};  // 0 new seed; -1 finished; n>0 in
                                     // progress, exact stage of growing for SDF
  MolFragment MoleculeFragment;  // Reference to a fragment of source molecule
  Graph Topology;  // seed topology with references to source molecule

  std::vector<bool> ExcludedBonds;
  unsigned LastAddedAtomsBeginIdx{0};  // in this subgraph for improving
                                       // performance of future growing
  unsigned LastAddedBondsBeginIdx{0};  // in this subgraph for DEBUG ONLY
  unsigned RemainingBonds{0};
  unsigned RemainingAtoms{0};
#ifdef DUP_SUBSTRUCT_CACHE
  DuplicatedSeedCache::TKey DupCacheKey;
#endif
  std::vector<TargetMatch> MatchResult;  // for each target
 public:
  Seed()

  {}

  void setMoleculeFragment(const Seed& src) {
    MoleculeFragment = src.MoleculeFragment;
  }
  Seed& operator=(const Seed& src) {
    NewBonds = src.NewBonds;
    GrowingStage = src.GrowingStage;
    MoleculeFragment = src.MoleculeFragment;
    Topology = src.Topology;
    ExcludedBonds = src.ExcludedBonds;
    LastAddedAtomsBeginIdx = src.LastAddedAtomsBeginIdx;
    LastAddedBondsBeginIdx = src.LastAddedBondsBeginIdx;
    RemainingBonds = src.RemainingBonds;
    RemainingAtoms = src.RemainingAtoms;
#ifdef DUP_SUBSTRUCT_CACHE
    DupCacheKey = src.DupCacheKey;
#endif
    MatchResult = src.MatchResult;
    CopyComplete = true;  // LAST
    return *this;
  }
  void createFromParent(const Seed* parent) {
    MoleculeFragment = parent->MoleculeFragment;
    Topology = parent->Topology;
    ExcludedBonds = parent->ExcludedBonds;
    RemainingBonds = parent->RemainingBonds;
    RemainingAtoms = parent->RemainingAtoms;
#ifdef DUP_SUBSTRUCT_CACHE
    DupCacheKey = parent->DupCacheKey;
#endif
    LastAddedAtomsBeginIdx = getNumAtoms();  // previous size
    LastAddedBondsBeginIdx = getNumBonds();  // previous size
    GrowingStage = 0;
  }

  unsigned getNumAtoms() const { return MoleculeFragment.AtomsIdx.size(); }
  unsigned getNumBonds() const { return MoleculeFragment.BondsIdx.size(); }

  void grow(MaximumCommonSubgraph& mcs) const;
  bool canGrowBiggerThan(unsigned maxBonds,
                         unsigned maxAtoms) const {  // prune()
    return RemainingBonds + getNumBonds() > maxBonds ||
           (RemainingBonds + getNumBonds() == maxBonds &&
            RemainingAtoms + getNumAtoms() > maxAtoms);
  }
  void computeRemainingSize(const ROMol& qmol);

  unsigned addAtom(const Atom* atom);
  unsigned addBond(const Bond* bond);
  void fillNewBonds(const ROMol& qmol);
};
}  // namespace FMCS
}  // namespace RDKit
