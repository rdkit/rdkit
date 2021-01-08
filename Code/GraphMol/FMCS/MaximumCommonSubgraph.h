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
#include "RingMatchTableSet.h"

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
  struct MCS {  // current result. Reference to a fragment of source molecule
    std::vector<const Atom*> Atoms;
    std::vector<const Bond*> Bonds;
    std::vector<unsigned> AtomsIdx;
    std::vector<unsigned> BondsIdx;  // need for results and size() only !
    const ROMol* QueryMolecule;
    std::vector<Target> Targets;
  };

  unsigned long long To;
  MCSProgressData Stat;
  MCSParameters Parameters;
  unsigned ThresholdCount;  // min number of matching
  std::vector<const ROMol*> Molecules;
  MCSCompareFunctionsData CompareFunctionsData;
#ifdef FAST_SUBSTRUCT_CACHE
  std::vector<unsigned> QueryAtomLabels;  // for code Morgan. Value based on
                                          // current functor and parameters
  std::vector<unsigned> QueryBondLabels;  // for code Morgan. Value based on
                                          // current functor and parameters
  SubstructureCache HashCache;
  MatchTable QueryAtomMatchTable;
  MatchTable QueryBondMatchTable;
  RingMatchTableSet RingMatchTables;
#endif
#ifdef DUP_SUBSTRUCT_CACHE
  DuplicatedSeedCache DuplicateCache;
#endif
  const ROMol* QueryMolecule;
  unsigned QueryMoleculeMatchedBonds;
  unsigned QueryMoleculeMatchedAtoms;
  const Atom* QueryMoleculeSingleMatchedAtom;
  std::vector<Target> Targets;
  SeedSet Seeds;
  MCS McsIdx;

 public:
#ifdef VERBOSE_STATISTICS_ON
  ExecStatistics VerboseStatistics;
#endif

  MaximumCommonSubgraph(const MCSParameters* params, void* userData);
  ~MaximumCommonSubgraph() { clear(); }
  MCSResult find(const std::vector<ROMOL_SPTR>& mols);
  const ROMol& getQueryMolecule() const { return *QueryMolecule; }
  unsigned getMaxNumberBonds() const { return McsIdx.BondsIdx.size(); }

  unsigned getMaxNumberAtoms() const { return McsIdx.AtomsIdx.size(); }
  // internal:
  bool checkIfMatchAndAppend(Seed& seed);

 private:
  void clear() {
    Targets.clear();
    Molecules.clear();
    CompareFunctionsData.clear();
    To = nanoClock();
  }
  std::pair<bool, bool> evaluateQueryMolecule(size_t i);
  void init();
  void loadMatchTables(const ROMol* targetMolecule, MatchTable& atomMatchTable, MatchTable& bondMatchTable);
  void loadBondLabels();
  void buildTargetTopology(Target& target);
  void printVerboseStatistics(const MCSResult& res) const;
  void makeInitialSeeds();
  void loadInitialSeedParameter(std::vector<bool>& excludedBonds);
  void loadSeedsFromQueryBonds(std::vector<bool>& excludedBonds);
  bool createSeedFromMCS(size_t newQueryTarget, Seed& seed);
  bool growSeeds();  // returns false if canceled
  std::pair<std::string, RWMol*> generateResultSMARTSAndQueryMol(
      const MCS& mcsIdx) const;
  bool addFusedBondQueries(const MCS& McsIdx, RWMol* rwMol) const;
  void recordChosenConformerIdx(const ROMol *mol, std::vector<unsigned>::const_iterator *conformerIterPtr);
  bool match(Seed& seed);
  bool matchIncrementalFast(Seed& seed, unsigned itarget);
};
}  // namespace FMCS
}  // namespace RDKit
