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
#include <vector>
#include <string>
#include <stdexcept>
#include "../RDKitBase.h"
#include "FMCS.h"
#include "DebugTrace.h"  // algorithm filter definitions
#include "SeedSet.h"
#include "Target.h"
#include "SubstructureCache.h"
#include "DuplicatedSeedCache.h"
#include "MatchTable.h"
#include "TargetMatch.h"

namespace RDKit {

inline bool FinalChiralityCheckFunction(
    const std::uint32_t c1[], const std::uint32_t c2[], const ROMol& mol1,
    const FMCS::Graph& query, const ROMol& mol2, const FMCS::Graph& target,
    const MCSParameters* p);

bool FinalMatchCheckFunction(const std::uint32_t c1[], const std::uint32_t c2[],
                             const ROMol& mol1, const FMCS::Graph& query,
                             const ROMol& mol2, const FMCS::Graph& target,
                             const MCSParameters* p);

namespace FMCS {
class RDKIT_FMCS_EXPORT MaximumCommonSubgraph {
  // current result. Reference to a fragment of source molecule
  struct MCS {
    std::vector<const Atom*> Atoms;
    std::vector<const Bond*> Bonds;
    const ROMol* QueryMolecule;
    std::vector<Target> Targets;
  };
  unsigned long long To;
  MCSProgressData Stat;
  detail::MCSParametersInternal Parameters;
  // min number of matches
  unsigned int ThresholdCount;
  std::vector<const ROMol*> Molecules;
#ifdef FAST_SUBSTRUCT_CACHE
  // for Morgan code. Value based on current functor and parameters
  std::vector<unsigned int> QueryAtomLabels;
  // for Morgan code. Value based on current functor and parameters
  std::vector<unsigned int> QueryBondLabels;
  SubstructureCache HashCache;
  MatchTable QueryAtomMatchTable;
  MatchTable QueryBondMatchTable;
#endif
#ifdef DUP_SUBSTRUCT_CACHE
  DuplicatedSeedCache DuplicateCache;
#endif
  const ROMol* QueryMolecule;
  unsigned int QueryMoleculeMatchedBonds;
  unsigned int QueryMoleculeMatchedAtoms;
  const Atom* QueryMoleculeSingleMatchedAtom;
  std::vector<Target> Targets;
  SeedSet Seeds;
  MCS McsIdx;
  std::map<std::vector<unsigned int>, MCS> DegenerateMcsMap;

 public:
#ifdef VERBOSE_STATISTICS_ON
  ExecStatistics VerboseStatistics;
#endif

  MaximumCommonSubgraph(const MCSParameters* params);
  ~MaximumCommonSubgraph() { clear(); }
  MCSResult find(const std::vector<ROMOL_SPTR>& mols);
  const ROMol& getQueryMolecule() const { return *QueryMolecule; }
  unsigned int getMaxNumberBonds() const { return McsIdx.Bonds.size(); }

  unsigned int getMaxNumberAtoms() const { return McsIdx.Atoms.size(); }
  bool checkIfMatchAndAppend(Seed& seed);
  bool match(Seed& seed);
  const MCSParameters& parameters() const { return Parameters; }
  MCSParameters& parameters() { return Parameters; }

 private:
  void clear() {
    Targets.clear();
    Molecules.clear();
    To = nanoClock();
  }
  void init(size_t startIdx);
  void makeInitialSeeds();
  bool createSeedFromMCS(size_t newQueryTarget, Seed& seed);
  bool growSeeds();  // returns false if canceled
  std::pair<std::string, ROMOL_SPTR> generateResultSMARTSAndQueryMol(
      const MCS& mcsIdx) const;

  bool matchIncrementalFast(Seed& seed, unsigned int itarget);
};
}  // namespace FMCS
}  // namespace RDKit
