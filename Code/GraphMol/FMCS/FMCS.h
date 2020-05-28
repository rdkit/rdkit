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
#include "Graph.h"

namespace RDKit {
struct MCSParameters;

typedef enum {
  AtomCompareAny,
  AtomCompareElements,
  AtomCompareIsotopes,
  AtomCompareAnyHeavyAtom
} AtomComparator;

typedef enum {
  BondCompareAny,
  BondCompareOrder,
  BondCompareOrderExact
} BondComparator;

typedef enum {
  IgnoreRingFusion,
  PermissiveRingFusion,
  StrictRingFusion
} RingComparator;

struct RDKIT_FMCS_EXPORT MCSAtomCompareParameters {
  bool MatchValences = false;
  bool MatchChiralTag = false;
  bool MatchFormalCharge = false;
  bool RingMatchesRingOnly = false;
  bool MatchIsotope = false;
};

struct RDKIT_FMCS_EXPORT MCSBondCompareParameters {
  bool RingMatchesRingOnly = false;
  bool CompleteRingsOnly = false;
  bool MatchFusedRings = false;
  bool MatchFusedRingsStrict = false;
  bool MatchStereo = false;
};

typedef bool (*MCSFinalMatchCheckFunction)(
    const std::uint32_t c1[], const std::uint32_t c2[], const ROMol& mol1,
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
RDKIT_FMCS_EXPORT bool checkAtomRingMatch(const MCSAtomCompareParameters& p,
                                          const ROMol& mol1, unsigned int atom1,
                                          const ROMol& mol2,
                                          unsigned int atom2);
RDKIT_FMCS_EXPORT bool checkAtomCharge(const MCSAtomCompareParameters& p,
                                       const ROMol& mol1, unsigned int atom1,
                                       const ROMol& mol2, unsigned int atom2);
RDKIT_FMCS_EXPORT bool checkAtomChirality(const MCSAtomCompareParameters& p,
                                          const ROMol& mol1, unsigned int atom1,
                                          const ROMol& mol2,
                                          unsigned int atom2);

RDKIT_FMCS_EXPORT bool MCSAtomCompareAny(const MCSAtomCompareParameters& p,
                                         const ROMol& mol1, unsigned int atom1,
                                         const ROMol& mol2, unsigned int atom2,
                                         void* userData);
RDKIT_FMCS_EXPORT bool MCSAtomCompareAnyHeavyAtom(
    const MCSAtomCompareParameters& p, const ROMol& mol1, unsigned int atom1,
    const ROMol& mol2, unsigned int atom2, void* userData);

RDKIT_FMCS_EXPORT bool MCSAtomCompareElements(
    const MCSAtomCompareParameters& p, const ROMol& mol1, unsigned int atom1,
    const ROMol& mol2, unsigned int atom2, void* userData);
RDKIT_FMCS_EXPORT bool MCSAtomCompareIsotopes(
    const MCSAtomCompareParameters& p, const ROMol& mol1, unsigned int atom1,
    const ROMol& mol2, unsigned int atom2, void* userData);

RDKIT_FMCS_EXPORT bool checkBondStereo(const MCSBondCompareParameters& p,
                                       const ROMol& mol1, unsigned int bond1,
                                       const ROMol& mol2, unsigned int bond2);
RDKIT_FMCS_EXPORT bool checkBondRingMatch(const MCSBondCompareParameters& p,
                                          const ROMol& mol1, unsigned int bond1,
                                          const ROMol& mol2, unsigned int bond2,
                                          void* v_ringMatchMatrixSet);

RDKIT_FMCS_EXPORT bool MCSBondCompareAny(const MCSBondCompareParameters& p,
                                         const ROMol& mol1, unsigned int bond1,
                                         const ROMol& mol2, unsigned int bond2,
                                         void* userData);
RDKIT_FMCS_EXPORT bool MCSBondCompareOrder(
    const MCSBondCompareParameters& p, const ROMol& mol1, unsigned int bond1,
    const ROMol& mol2, unsigned int bond2,
    void* userData);  // ignore Aromatization
RDKIT_FMCS_EXPORT bool MCSBondCompareOrderExact(
    const MCSBondCompareParameters& p, const ROMol& mol1, unsigned int bond1,
    const ROMol& mol2, unsigned int bond2, void* userData);

struct RDKIT_FMCS_EXPORT MCSProgressData {
  unsigned NumAtoms{0};
  unsigned NumBonds{0};
  unsigned SeedProcessed{0};

 public:
  MCSProgressData() {}
};

typedef bool (*MCSProgressCallback)(const MCSProgressData& stat,
                                    const MCSParameters& params,
                                    void* userData);
RDKIT_FMCS_EXPORT bool MCSProgressCallbackTimeout(const MCSProgressData& stat,
                                                  const MCSParameters& params,
                                                  void* userData);

struct RDKIT_FMCS_EXPORT MCSParameters {
  bool MaximizeBonds = true;
  double Threshold = 1.0;  // match all molecules
  unsigned Timeout = -1;   // in seconds
  bool Verbose = false;
  MCSAtomCompareParameters AtomCompareParameters;
  MCSBondCompareParameters BondCompareParameters;
  MCSAtomCompareFunction AtomTyper = MCSAtomCompareElements;
  MCSBondCompareFunction BondTyper = MCSBondCompareOrder;
  void* CompareFunctionsUserData = nullptr;
  MCSProgressCallback ProgressCallback =
      nullptr;  // return false to interrupt execution
  void* ProgressCallbackUserData = nullptr;
  MCSFinalMatchCheckFunction FinalMatchChecker =
      nullptr;  // FinalMatchCheckFunction() to check chirality and ring fusion
  std::string InitialSeed = "";  // user defined or empty string (default)
  void setMCSAtomTyperFromEnum(AtomComparator atomComp);
  void setMCSAtomTyperFromConstChar(const char* atomComp);
  void setMCSBondTyperFromEnum(BondComparator bondComp);
  void setMCSBondTyperFromConstChar(const char* bondComp);
};

struct RDKIT_FMCS_EXPORT MCSResult {
  unsigned NumAtoms{0};
  unsigned NumBonds{0};
  std::string SmartsString;
  bool Canceled{false};  // interrupted by timeout or user defined progress
                         // callback. Contains valid current MCS !
  ROMOL_SPTR QueryMol;

 public:
  MCSResult() {}
  bool isCompleted() const { return !Canceled; }
};

RDKIT_FMCS_EXPORT void parseMCSParametersJSON(const char* json,
                                              MCSParameters* params);

RDKIT_FMCS_EXPORT MCSResult findMCS(const std::vector<ROMOL_SPTR>& mols,
                                    const MCSParameters* params = nullptr);
RDKIT_FMCS_EXPORT MCSResult findMCS_P(const std::vector<ROMOL_SPTR>& mols,
                                      const char* params_json);

RDKIT_FMCS_EXPORT MCSResult findMCS(
    const std::vector<ROMOL_SPTR>& mols, bool maximizeBonds, double threshold,
    unsigned timeout, bool verbose, bool matchValences,
    bool ringMatchesRingOnly, bool completeRingsOnly, bool matchChiralTag,
    AtomComparator atomComp, BondComparator bondComp, RingComparator ringComp);
RDKIT_FMCS_EXPORT MCSResult
findMCS(const std::vector<ROMOL_SPTR>& mols, bool maximizeBonds,
        double threshold = 1.0, unsigned timeout = 3600, bool verbose = false,
        bool matchValences = false, bool ringMatchesRingOnly = false,
        bool completeRingsOnly = false, bool matchChiralTag = false,
        AtomComparator atomComp = AtomCompareElements,
        BondComparator bondComp = BondCompareOrder);

}  // namespace RDKit
