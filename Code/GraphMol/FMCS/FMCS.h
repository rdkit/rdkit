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
  bool CompleteRingsOnly = false;
  bool MatchIsotope = false;
  double MaxDistance = -1.0;
};

struct RDKIT_FMCS_EXPORT MCSBondCompareParameters {
  bool RingMatchesRingOnly = false;
  bool CompleteRingsOnly = false;
  bool MatchFusedRings = false;
  bool MatchFusedRingsStrict = false;
  bool MatchStereo = false;
};

typedef bool (*MCSAtomCompareFunction)(const MCSAtomCompareParameters &,
                                       const ROMol &, unsigned int,
                                       const ROMol &, unsigned int, void *);
typedef bool (*MCSBondCompareFunction)(const MCSBondCompareParameters &,
                                       const ROMol &, unsigned int,
                                       const ROMol &, unsigned int, void *);
typedef bool (*MCSAcceptanceFunction)(const ROMol &, const ROMol &,
                                      const std::vector<std::pair<int, int>> &,
                                      const std::vector<std::pair<int, int>> &,
                                      const MCSParameters *);
typedef bool (*MCSFinalMatchCheckFunction)(const std::uint32_t[],
                                           const std::uint32_t[], const ROMol &,
                                           const FMCS::Graph &, const ROMol &,
                                           const FMCS::Graph &,
                                           const MCSParameters *);

// Some predefined functors:
RDKIT_FMCS_EXPORT bool checkAtomRingMatch(const MCSAtomCompareParameters &p,
                                          const ROMol &mol1, unsigned int atom1,
                                          const ROMol &mol2,
                                          unsigned int atom2);
RDKIT_FMCS_EXPORT bool checkAtomCharge(const MCSAtomCompareParameters &p,
                                       const ROMol &mol1, unsigned int atom1,
                                       const ROMol &mol2, unsigned int atom2);
RDKIT_FMCS_EXPORT bool checkAtomChirality(const MCSAtomCompareParameters &p,
                                          const ROMol &mol1, unsigned int atom1,
                                          const ROMol &mol2,
                                          unsigned int atom2);
RDKIT_FMCS_EXPORT bool checkAtomDistance(const MCSAtomCompareParameters &p,
                                         const ROMol &mol1, unsigned int atom1,
                                         const ROMol &mol2, unsigned int atom2);

RDKIT_FMCS_EXPORT bool MCSAtomCompareAny(const MCSAtomCompareParameters &p,
                                         const ROMol &mol1, unsigned int atom1,
                                         const ROMol &mol2, unsigned int atom2,
                                         void *userData);
RDKIT_FMCS_EXPORT bool MCSAtomCompareAnyHeavyAtom(
    const MCSAtomCompareParameters &p, const ROMol &mol1, unsigned int atom1,
    const ROMol &mol2, unsigned int atom2, void *userData);

RDKIT_FMCS_EXPORT bool MCSAtomCompareElements(
    const MCSAtomCompareParameters &p, const ROMol &mol1, unsigned int atom1,
    const ROMol &mol2, unsigned int atom2, void *userData);
RDKIT_FMCS_EXPORT bool MCSAtomCompareIsotopes(
    const MCSAtomCompareParameters &p, const ROMol &mol1, unsigned int atom1,
    const ROMol &mol2, unsigned int atom2, void *userData);

RDKIT_FMCS_EXPORT bool checkBondStereo(const MCSBondCompareParameters &p,
                                       const ROMol &mol1, unsigned int bond1,
                                       const ROMol &mol2, unsigned int bond2);

RDKIT_FMCS_EXPORT bool havePairOfCompatibleRings(
    const MCSBondCompareParameters &p, const ROMol &mol1, unsigned int bond1,
    const ROMol &mol2, unsigned int bond2);

RDKIT_FMCS_EXPORT bool checkBondRingMatch(const MCSBondCompareParameters &p,
                                          const ROMol &mol1, unsigned int bond1,
                                          const ROMol &mol2,
                                          unsigned int bond2);

RDKIT_FMCS_EXPORT bool MCSBondCompareAny(const MCSBondCompareParameters &p,
                                         const ROMol &mol1, unsigned int bond1,
                                         const ROMol &mol2, unsigned int bond2,
                                         void *userData);
RDKIT_FMCS_EXPORT bool MCSBondCompareOrder(
    const MCSBondCompareParameters &p, const ROMol &mol1, unsigned int bond1,
    const ROMol &mol2, unsigned int bond2,
    void *userData);  // ignore Aromatization
RDKIT_FMCS_EXPORT bool MCSBondCompareOrderExact(
    const MCSBondCompareParameters &p, const ROMol &mol1, unsigned int bond1,
    const ROMol &mol2, unsigned int bond2, void *userData);

struct RDKIT_FMCS_EXPORT MCSProgressData {
  unsigned int NumAtoms{0};
  unsigned int NumBonds{0};
  unsigned int SeedProcessed{0};

 public:
  MCSProgressData() {}
};

typedef bool (*MCSProgressCallback)(const MCSProgressData &stat,
                                    const MCSParameters &params,
                                    void *userData);
RDKIT_FMCS_EXPORT bool MCSProgressCallbackTimeout(const MCSProgressData &stat,
                                                  const MCSParameters &params,
                                                  void *userData);

struct RDKIT_FMCS_EXPORT MCSParameters {
  MCSParameters() {}
  MCSParameters(const MCSParameters *other) : MCSParameters() {
    if (other) {
      *this = *other;
    }
  }
  MCSParameters(const MCSParameters &other) = default;
  MCSParameters &operator=(const MCSParameters &other) = default;
  virtual ~MCSParameters() {}
  bool StoreAll = false;
  bool MaximizeBonds = true;
  double Threshold = 1.0;    // match all molecules
  unsigned int Timeout = 0;  // in seconds
  bool Verbose = false;
  MCSAtomCompareParameters AtomCompareParameters;
  MCSBondCompareParameters BondCompareParameters;
  MCSAtomCompareFunction AtomTyper = MCSAtomCompareElements;
  MCSBondCompareFunction BondTyper = MCSBondCompareOrder;
  void *CompareFunctionsUserData = nullptr;
  MCSProgressCallback ProgressCallback =
      nullptr;  // return false to interrupt execution
  void *ProgressCallbackUserData = nullptr;
  // FinalMatchCheckFunction() to accept/reject a growing MCS candidate based on
  // user-defined criteria
  MCSFinalMatchCheckFunction FinalMatchChecker = nullptr;
  void *FinalMatchCheckerUserData = nullptr;
  // ShouldAcceptMCS() to accept/reject a fully-grown MCS candidate based on
  // user-defined criteria
  MCSAcceptanceFunction ShouldAcceptMCS = nullptr;
  void *ShouldAcceptMCSUserData = nullptr;
  std::string InitialSeed = "";  // user defined or empty string (default)
  void setMCSAtomTyperFromEnum(AtomComparator atomComp);
  void setMCSAtomTyperFromConstChar(const char *atomComp);
  void setMCSBondTyperFromEnum(BondComparator bondComp);
  void setMCSBondTyperFromConstChar(const char *bondComp);
};

namespace detail {
struct MCSParametersInternal : public MCSParameters {
  MCSParametersInternal() {}
  ~MCSParametersInternal() {}
  MCSParametersInternal(const MCSParameters *params);
  MCSFinalMatchCheckFunction UserFinalMatchChecker = nullptr;
};
}  // end namespace detail

struct RDKIT_FMCS_EXPORT MCSResult {
  unsigned int NumAtoms{0};
  unsigned int NumBonds{0};
  std::string SmartsString;
  bool Canceled{false};  // interrupted by timeout or user defined progress
                         // callback. Contains valid current MCS !
  ROMOL_SPTR QueryMol;
  std::map<std::string, ROMOL_SPTR> DegenerateSmartsQueryMolDict;

 public:
  MCSResult() {}
  bool isCompleted() const { return !Canceled; }
};

RDKIT_FMCS_EXPORT void parseMCSParametersJSON(const char *json,
                                              MCSParameters *params);

RDKIT_FMCS_EXPORT MCSResult findMCS(const std::vector<ROMOL_SPTR> &mols,
                                    const MCSParameters *params = nullptr);
RDKIT_FMCS_EXPORT MCSResult findMCS_P(const std::vector<ROMOL_SPTR> &mols,
                                      const char *params_json);

RDKIT_FMCS_EXPORT MCSResult findMCS(
    const std::vector<ROMOL_SPTR> &mols, bool maximizeBonds, double threshold,
    unsigned int timeout, bool verbose, bool matchValences,
    bool ringMatchesRingOnly, bool completeRingsOnly, bool matchChiralTag,
    AtomComparator atomComp, BondComparator bondComp, RingComparator ringComp);
RDKIT_FMCS_EXPORT MCSResult findMCS(
    const std::vector<ROMOL_SPTR> &mols, bool maximizeBonds,
    double threshold = 1.0, unsigned int timeout = 3600, bool verbose = false,
    bool matchValences = false, bool ringMatchesRingOnly = false,
    bool completeRingsOnly = false, bool matchChiralTag = false,
    AtomComparator atomComp = AtomCompareElements,
    BondComparator bondComp = BondCompareOrder);

}  // namespace RDKit
