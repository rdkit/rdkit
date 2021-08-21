//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file StructChecker.h

\brief Contains the public API of the StructChecker

\b Note that this should be considered beta and that the API may change in
future
releases.

*/
#include <RDGeneral/export.h>
#pragma once
#ifndef RD_STRUCTCHECKER_H_Oct2016
#define RD_STRUCTCHECKER_H_Oct2016

#include <string>
#include <vector>
#include "../RDKitBase.h"

namespace RDKit {
namespace StructureCheck {

// Flags for the return values of the StructureChecker

// TypeDefs for translating augmented atom pairs
static const int ANY_CHARGE = 8;
enum RadicalType {
  RT_NONE = 0,
  SINGLET = 1,
  DOUBLET = 2,
  TRIPLET = 3,
  ANY_RADICAL = 0xFF
};

enum AABondType {  // MDL CTFile bond types plus extensions
  BT_NONE = 0,     // means REMOVE Bond
  SINGLE = 1,
  DOUBLE = 2,
  TRIPLE = 3,
  AROMATIC = 4,
  SINGLE_DOUBLE = 5,
  SINGLE_AROMATIC = 6,
  DOUBLE_AROMATIC = 7,
  ANY_BOND = 8,
  ALL_BOND_TYPES = 0xF
};

enum AATopology {
  TP_NONE = 0,  // Don't care
  RING = 1,     // Ring
  CHAIN = 2     // Chain
};

struct RDKIT_STRUCTCHECKER_EXPORT Ligand {
  std::string AtomSymbol;
  int Charge;
  RadicalType Radical;
  unsigned SubstitutionCount;  // substitution count 0 = don't care
  AABondType BondType;
  Ligand()
      : Charge(ANY_CHARGE),
        Radical(ANY_RADICAL),
        SubstitutionCount(0),
        BondType(ANY_BOND) {}
};

struct RDKIT_STRUCTCHECKER_EXPORT AugmentedAtom {
  std::string AtomSymbol;
  std::string ShortName;
  int Charge;
  RadicalType Radical;
  AATopology Topology;
  std::vector<Ligand> Ligands;

  AugmentedAtom()
      : Charge(ANY_CHARGE), Radical(ANY_RADICAL), Topology(TP_NONE) {}

  AugmentedAtom(const std::string &symbol, const std::string &name, int charge,
                RadicalType radical, AATopology topology)
      : AtomSymbol(symbol),
        ShortName(name),
        Charge(charge),
        Radical(radical),
        Topology(topology) {}
};

struct RDKIT_STRUCTCHECKER_EXPORT IncEntry {
  std::string AtomSymbol;
  double LocalInc;
  double AlphaInc;
  double BetaInc;
  double MultInc;

  // Used for logging
  int local_inc_used;
  int alpha_inc_used;
  int beta_inc_used;
  int mult_inc_used;
};

struct RDKIT_STRUCTCHECKER_EXPORT PathEntry {
  AugmentedAtom Path;
  double Cond;
  // Used for logging
  int cond_used;
};
//-------------

//! Structure Check Options
///   Holds all the user options for the StructureChecking.
///   Can be initialized from factory functions, perhaps serialized
struct RDKIT_STRUCTCHECKER_EXPORT StructCheckerOptions {
  double AcidityLimit;
  bool RemoveMinorFragments;
  int DesiredCharge;
  bool CheckCollisions;
  int CollisionLimitPercent;
  unsigned MaxMolSize;
  bool ConvertSText;
  bool SqueezeIdentifiers;
  bool StripZeros;
  bool CheckStereo;
  bool ConvertAtomTexts;
  bool GroupsToSGroups;
  bool Verbose;

  // Internal data for struchk
  std::vector<std::pair<AugmentedAtom, AugmentedAtom>> AugmentedAtomPairs;
  std::vector<AugmentedAtom> AcidicAtoms;
  std::vector<AugmentedAtom> GoodAtoms;
  std::vector<ROMOL_SPTR> Patterns;
  std::vector<ROMOL_SPTR> RotatePatterns;
  std::vector<ROMOL_SPTR> StereoPatterns;
  std::vector<ROMOL_SPTR> FromTautomer;
  std::vector<ROMOL_SPTR> ToTautomer;

  double Elneg0;                          // elneg_table[0].value;
  std::map<unsigned, double> ElnegTable;  // AtomicNumber -> eleng
  std::vector<IncEntry> AtomAcidity;      // atom_acidity_table[]
  std::vector<IncEntry> ChargeIncTable;
  // std::map AtomSymbol(or AtomicNumber) -> IncEntry
  /* [ReadTransformation() ]
   * The alpha, beta coefficients of the transfomation function used
   * to stretch the preliminary pKa values to the actual predictions.
   * The function is pKa = 7 + (pKa'-7)*beta + ((pKa'-7)*alpha)^3.
   */

  double Alpha, Beta;
  std::vector<PathEntry> AlphaPathTable, BetaPathTable;

 public:
  StructCheckerOptions();

  void clear() { *this = StructCheckerOptions(); }

  bool loadAugmentedAtomTranslations(const std::string &path);
  void setAugmentedAtomTranslations(
      const std::vector<std::pair<AugmentedAtom, AugmentedAtom>> &aaPairs);

  bool loadAcidicAugmentedAtoms(const std::string &path);
  void setAcidicAugmentedAtoms(const std::vector<AugmentedAtom> &acidicAtoms);

  bool loadGoodAugmentedAtoms(const std::string &path);
  void setGoodAugmentedAtoms(const std::vector<AugmentedAtom> &acidicAtoms);

  bool loadPatterns(const std::string &path);  // file with clean patterns
  void parsePatterns(
      const std::vector<std::string> &smarts);  // can throw RDKit exceptions
  void setPatterns(const std::vector<ROMOL_SPTR> &p);

  bool loadRotatePatterns(
      const std::string &path);  // file with rotate patterns
  void parseRotatePatterns(
      const std::vector<std::string> &smarts);  // can throw RDKit exceptions
  void setRotatePatterns(const std::vector<ROMOL_SPTR> &p);

  bool loadStereoPatterns(
      const std::string &path);  // file with stereo patterns
  void parseStereoPatterns(
      const std::vector<std::string> &smarts);  // can throw RDKit exceptions
  void setStereoPatterns(const std::vector<ROMOL_SPTR> &p);

  bool loadTautomerData(const std::string &path);  // file path
  void parseTautomerData(const std::vector<std::string> &smartsFrom,
                         const std::vector<std::string> &smartsTo);
  void setTautomerData(const std::vector<ROMOL_SPTR> &from,
                       const std::vector<ROMOL_SPTR> &to);
  bool loadChargeDataTables(const std::string &path);  // file path
};

RDKIT_STRUCTCHECKER_EXPORT bool parseOptionsJSON(const std::string &json,
                                                 StructCheckerOptions &op);

RDKIT_STRUCTCHECKER_EXPORT bool loadOptionsFromFiles(
    StructCheckerOptions &op,
    const std::string &augmentedAtomTranslationsFile = "",
    // ?? AcidicAtoms;
    // ?? GoodAtoms;
    const std::string &patternFile = "",        // file with clean patterns
    const std::string &rotatePatternFile = "",  // file with rotate patterns
    const std::string &stereoPatternFile = "",  // file with stereo patterns
    const std::string &tautomerFile = "");

//! \brief Class for performing structure validation and cleanup
/*! \b NOTE: This class should be considered beta. The API may change in future
releases.

Examples of Usage

\code
  StructChecker chk;
  int flags = StructureCheck::checkMolStructure( mol ); // use defaults
\endcode

or

\code
    StructureCheck::StructCheckerOptions options;   // use defaults
    // To use external data
    StructureCheck::loadOptionsFromFiles(options, file1, file2);
    StructChecker chk(options);

    for( mol in mols ) {
        int flags = StructureCheck::checkMolStructure( mol, &options);
        if (0!=(flags & StructureCheck::StructureFlags::BAD_SET)) {
        // write to error file
        } else if (0!=(flags & StructureCheck::StructureFlags::TRANSFORMED_SET))
{
        // input molecule was transformed
        } else { // flag == NO_CHANGE
        // no change
        }
    }
\endcode
*/
class RDKIT_STRUCTCHECKER_EXPORT StructChecker {
 public:
  typedef enum StructureFlags {
    NO_CHANGE = 0,
    BAD_MOLECULE = 0x0001,
    ALIAS_CONVERSION_FAILED = 0x0002,
    STEREO_ERROR = 0x0004,
    STEREO_FORCED_BAD = 0x0008,
    ATOM_CLASH = 0x0010,
    ATOM_CHECK_FAILED = 0x0020,
    SIZE_CHECK_FAILED = 0x0040,
    // reserved error = 0x0080,
    TRANSFORMED = 0x0100,
    FRAGMENTS_FOUND = 0x0200,
    EITHER_WARNING = 0x0400,
    DUBIOUS_STEREO_REMOVED = 0x0800,
    RECHARGED = 0x1000,
    STEREO_TRANSFORMED = 0x2000,
    TEMPLATE_TRANSFORMED = 0x4000,
    TAUTOMER_TRANSFORMED = 0x8000,
    // mask:
    BAD_SET = (BAD_MOLECULE | ALIAS_CONVERSION_FAILED | STEREO_ERROR |
               STEREO_FORCED_BAD | ATOM_CLASH | ATOM_CHECK_FAILED |
               SIZE_CHECK_FAILED),

    TRANSFORMED_SET = (TRANSFORMED | FRAGMENTS_FOUND | EITHER_WARNING |
                       DUBIOUS_STEREO_REMOVED | STEREO_TRANSFORMED |
                       TEMPLATE_TRANSFORMED | TAUTOMER_TRANSFORMED | RECHARGED),
  } StructureFlags;
  // attributes:
 private:
  StructCheckerOptions Options;

 public:
  inline StructChecker() {}
  inline StructChecker(const StructCheckerOptions &options)
      : Options(options) {}

  const StructCheckerOptions &GetOptions() const { return Options; }
  void SetOptions(const StructCheckerOptions &options) { Options = options; }

  // Check and fix (if need) molecule structure and return a set of
  // StructureFlags
  // that describes what have been done
  unsigned checkMolStructure(RWMol &mol) const;

  // an instance independent helper methods:
  // Converts structure property flags to a comma separated string
  static std::string StructureFlagsToString(unsigned flags);
  // Converts a comma separated string to a StructureFlag unsigned integer
  static unsigned StringToStructureFlags(const std::string &str);
  // internal implementation:
 private:
};
}  // namespace StructureCheck
}  // namespace RDKit
#endif
