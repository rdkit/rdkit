//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once
#include <vector>
#include <string>
#include <stdexcept>
#include "../RDKitBase.h"
#include "Graph.h"
#ifndef RDKIT_WRAP_DECL
#define RDKIT_WRAP_DECL
#endif

namespace RDKit {
struct MCSParameters;

struct MCSAtomCompareParameters {
  bool MatchValences;
  bool MatchChiralTag;
  bool MatchFormalCharge;

 public:
  MCSAtomCompareParameters()
      : MatchValences(false), MatchChiralTag(false), MatchFormalCharge(false) {}
};

struct MCSBondCompareParameters {
  bool RingMatchesRingOnly;
  bool CompleteRingsOnly;
  bool MatchStereo;

 public:
  MCSBondCompareParameters()
      : RingMatchesRingOnly(false),
        CompleteRingsOnly(false),
        MatchStereo(false) {}
};

typedef bool (*MCSFinalMatchCheckFunction)(
    const short unsigned c1[], const short unsigned c2[], const ROMol& mol1,
    const FMCS::Graph& query, const ROMol& mol2, const FMCS::Graph& target,
    const MCSParameters* p);
typedef bool (*MCSAtomCompareFunction)(const MCSAtomCompareParameters& p,
                                       const ROMol& mol1, unsigned int atom1,
                                       const ROMol& mol2, unsigned int atom2,
                                       void* userData);
typedef bool (*MCSBondCompareFunction)(const MCSBondCompareParameters& p,
                                       const ROMol& mol1, unsigned int bond1,
                                       const ROMol& mol2, unsigned int bond2,
                                       void* userData);

// Some predefined functors:
RDKIT_WRAP_DECL bool MCSAtomCompareAny(const MCSAtomCompareParameters& p,
                                       const ROMol& mol1, unsigned int atom1,
                                       const ROMol& mol2, unsigned int atom2,
                                       void* userData);

RDKIT_WRAP_DECL bool MCSAtomCompareElements(const MCSAtomCompareParameters& p,
                                            const ROMol& mol1,
                                            unsigned int atom1,
                                            const ROMol& mol2,
                                            unsigned int atom2, void* userData);
RDKIT_WRAP_DECL bool MCSAtomCompareIsotopes(const MCSAtomCompareParameters& p,
                                            const ROMol& mol1,
                                            unsigned int atom1,
                                            const ROMol& mol2,
                                            unsigned int atom2, void* userData);

RDKIT_WRAP_DECL bool MCSBondCompareAny(const MCSBondCompareParameters& p,
                                       const ROMol& mol1, unsigned int bond1,
                                       const ROMol& mol2, unsigned int bond2,
                                       void* userData);
RDKIT_WRAP_DECL bool MCSBondCompareOrder(
    const MCSBondCompareParameters& p, const ROMol& mol1, unsigned int bond1,
    const ROMol& mol2, unsigned int bond2,
    void* userData);  // ignore Aromatization
RDKIT_WRAP_DECL bool MCSBondCompareOrderExact(
    const MCSBondCompareParameters& p, const ROMol& mol1, unsigned int bond1,
    const ROMol& mol2, unsigned int bond2, void* userData);

struct MCSProgressData {
  unsigned NumAtoms;
  unsigned NumBonds;
  unsigned SeedProcessed;

 public:
  MCSProgressData() : NumAtoms(0), NumBonds(0), SeedProcessed(0) {}
};

typedef bool (*MCSProgressCallback)(const MCSProgressData& stat,
                                    const MCSParameters& params,
                                    void* userData);
bool MCSProgressCallbackTimeout(const MCSProgressData& stat,
                                const MCSParameters& params, void* userData);

struct MCSParameters {
  bool MaximizeBonds;
  double Threshold;
  unsigned Timeout;  // in seconds
  bool Verbose;
  MCSAtomCompareParameters AtomCompareParameters;
  MCSBondCompareParameters BondCompareParameters;
  MCSAtomCompareFunction AtomTyper;
  MCSBondCompareFunction BondTyper;
  void* CompareFunctionsUserData;
  MCSProgressCallback ProgressCallback;  // return false to interrupt execution
  void* ProgressCallbackUserData;
  MCSFinalMatchCheckFunction
      FinalMatchChecker;    // FinalChiralityCheckFunction() to check chirality
  std::string InitialSeed;  // user defined or empty string (default)
 public:
  MCSParameters()
      : MaximizeBonds(true),
        Threshold(1.0)  // match to all
        ,
        Timeout(-1),
        Verbose(false),
        AtomTyper(MCSAtomCompareElements),
        BondTyper(MCSBondCompareOrder),
        CompareFunctionsUserData(0),
        ProgressCallback(0),
        ProgressCallbackUserData(0),
        FinalMatchChecker(0),
        InitialSeed("") {}
};

struct MCSResult {
  unsigned NumAtoms;
  unsigned NumBonds;
  std::string SmartsString;
  bool Canceled;  // interrupted by timeout or user defined progress callback.
                  // Contains valid current MCS !
 public:
  MCSResult() : NumAtoms(0), NumBonds(0), Canceled(false) {}
  bool isCompleted() const { return !Canceled; }
};

void parseMCSParametersJSON(const char* json, MCSParameters* params);

MCSResult findMCS(const std::vector<ROMOL_SPTR>& mols,
                  const MCSParameters* params = 0);
MCSResult findMCS_P(const std::vector<ROMOL_SPTR>& mols,
                    const char* params_json);

typedef enum {
  AtomCompareAny,
  AtomCompareElements,
  AtomCompareIsotopes
} AtomComparator;
typedef enum {
  BondCompareAny,
  BondCompareOrder,
  BondCompareOrderExact
} BondComparator;
MCSResult findMCS(const std::vector<ROMOL_SPTR>& mols, bool maximizeBonds,
                  double threshold = 1.0, unsigned timeout = 3600,
                  bool verbose = false, bool matchValences = false,
                  bool ringMatchesRingOnly = false,
                  bool completeRingsOnly = false, bool matchChiralTag = false,
                  AtomComparator atomComp = AtomCompareElements,
                  BondComparator bondComp = BondCompareOrder);

}  // namespace RDKit
