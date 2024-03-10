//
//  Copyright (C) 2001-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_SMILESPARSE_H
#define RD_SMILESPARSE_H

#include <GraphMol/SanitException.h>
#include <string>
#include <exception>
#include <map>
#include <memory>

namespace RDKit {
class RWMol;
class Atom;
class Bond;

namespace SmilesParse {
class RDKIT_SMILESPARSE_EXPORT SmilesParseException : public std::exception {
 public:
  SmilesParseException(const char *msg) : _msg(msg) {}
  SmilesParseException(const std::string msg) : _msg(msg) {}
  const char *what() const noexcept override { return _msg.c_str(); }
  ~SmilesParseException() noexcept override = default;

 private:
  std::string _msg;
};

}  // namespace SmilesParse

namespace v2 {
namespace SmilesParse {
using RDKit::SmilesParse::SmilesParseException;

struct RDKIT_SMILESPARSE_EXPORT SmilesParserParams {
  bool sanitize = true;      /**< sanitize the molecule after building it */
  bool allowCXSMILES = true; /**< recognize and parse CXSMILES*/
  bool strictCXSMILES =
      true; /**< throw an exception if the CXSMILES parsing fails */
  bool parseName = true;    /**< parse (and set) the molecule name as well */
  bool removeHs = true;     /**< remove Hs after constructing the molecule */
  bool skipCleanup = false; /**<  skip the final cleanup stage */
  bool debugParse = false;  /**< enable debugging in the SMILES parser*/
  std::map<std::string, std::string>
      replacements; /**< allows SMILES "macros" */
};

struct RDKIT_SMILESPARSE_EXPORT SmartsParserParams {
  bool allowCXSMILES = true; /**< recognize and parse CXSMILES extensions */
  bool strictCXSMILES =
      true; /**< throw an exception if the CXSMILES parsing fails */
  bool parseName = true; /**< parse (and set) the molecule name as well */
  bool mergeHs =
      false; /**< toggles merging H atoms in the SMARTS into neighboring atoms*/
  bool skipCleanup = false; /**<  skip the final cleanup stage */
  bool debugParse = false;  /**< enable debugging in the SMARTS parser*/
  std::map<std::string, std::string>
      replacements; /**< allows SMARTS "macros" */
};

RDKIT_SMILESPARSE_EXPORT std::unique_ptr<RDKit::RWMol> MolFromSmiles(
    const std::string &smi,
    const SmilesParserParams &params = SmilesParserParams());
RDKIT_SMILESPARSE_EXPORT std::unique_ptr<RDKit::RWMol> MolFromSmarts(
    const std::string &sma,
    const SmartsParserParams &params = SmartsParserParams());

RDKIT_SMILESPARSE_EXPORT std::unique_ptr<RDKit::Atom> AtomFromSmiles(
    const std::string &smi);
RDKIT_SMILESPARSE_EXPORT std::unique_ptr<RDKit::Bond> BondFromSmiles(
    const std::string &smi);

RDKIT_SMILESPARSE_EXPORT std::unique_ptr<RDKit::Atom> AtomFromSmarts(
    const std::string &sma);
RDKIT_SMILESPARSE_EXPORT std::unique_ptr<RDKit::Bond> BondFromSmarts(
    const std::string &sma);

}  // namespace SmilesParse
}  // namespace v2

inline namespace v1 {
using RDKit::SmilesParse::SmilesParseException;

struct RDKIT_SMILESPARSE_EXPORT SmilesParserParams {
  int debugParse = 0;   /**< enable debugging in the SMILES parser*/
  bool sanitize = true; /**< sanitize the molecule after building it */
  std::map<std::string, std::string> *replacements =
      nullptr;               /**< allows SMILES "macros" */
  bool allowCXSMILES = true; /**< recognize and parse CXSMILES*/
  bool strictCXSMILES =
      true; /**< throw an exception if the CXSMILES parsing fails */
  bool parseName = true;    /**< parse (and set) the molecule name as well */
  bool removeHs = true;     /**< remove Hs after constructing the molecule */
  bool skipCleanup = false; /**<  skip the final cleanup stage */
};

struct RDKIT_SMILESPARSE_EXPORT SmartsParserParams {
  int debugParse = 0; /**< enable debugging in the SMARTS parser*/
  std::map<std::string, std::string> *replacements =
      nullptr;               /**< allows SMARTS "macros" */
  bool allowCXSMILES = true; /**< recognize and parse CXSMILES extensions */
  bool strictCXSMILES =
      true; /**< throw an exception if the CXSMILES parsing fails */
  bool parseName = true; /**< parse (and set) the molecule name as well */
  bool mergeHs =
      false; /**< toggles merging H atoms in the SMARTS into neighboring atoms*/
  bool skipCleanup = false; /**<  skip the final cleanup stage */
};

inline RDKit::RWMol *SmilesToMol(const std::string &smi,
                                 const SmilesParserParams &ps) {
  RDKit::v2::SmilesParse::SmilesParserParams v2ps;
  v2ps.debugParse = ps.debugParse;
  v2ps.sanitize = ps.sanitize;

  if (ps.replacements) {
    v2ps.replacements = *ps.replacements;
  }
  v2ps.allowCXSMILES = ps.allowCXSMILES;
  v2ps.strictCXSMILES = ps.strictCXSMILES;
  v2ps.parseName = ps.parseName;
  v2ps.removeHs = ps.removeHs;
  v2ps.skipCleanup = ps.skipCleanup;
  return RDKit::v2::SmilesParse::MolFromSmiles(smi, v2ps).release();
}

inline Atom *SmilesToAtom(const std::string &smi) {
  auto res = RDKit::v2::SmilesParse::AtomFromSmiles(smi).release();
  return res;
}

inline Bond *SmilesToBond(const std::string &smi) {
  return RDKit::v2::SmilesParse::BondFromSmiles(smi).release();
}

//! Construct a molecule from a SMILES string
/*!
 \param smi           the SMILES to convert
 \param debugParse    toggles verbose debugging information from the parser
 \param sanitize      toggles H removal and sanitization of the molecule
 \param replacements  a string->string map of replacement strings. See below
                      for more information about replacements.

 \return a pointer to the new molecule; the caller is responsible for free'ing
 this.

 The optional replacements map can be used to do string substitution of
 abbreviations
 in the input SMILES. The set of substitutions is repeatedly looped through
 until
 the string no longer changes. It is the responsibility of the caller to make
 sure
 that substitutions results in legal and sensible SMILES.

 Examples of substitutions:
 \code
   CC{Q}C with {"{Q}":"OCCO"} -> CCOCCOC
   C{A}C{Q}C with {"{Q}":"OCCO", "{A}":"C1(CC1)"} -> CC1(CC1)COCCOC
   C{A}C{Q}C with {"{Q}":"{X}CC{X}", "{A}":"C1CC1", "{X}":"N"} -> CC1CC1CNCCNC
 \endcode

 */
inline RWMol *SmilesToMol(
    const std::string &smi, int debugParse = 0, bool sanitize = true,
    std::map<std::string, std::string> *replacements = nullptr) {
  RDKit::v2::SmilesParse::SmilesParserParams params;
  params.debugParse = debugParse;
  if (replacements) {
    params.replacements = *replacements;
  }
  if (sanitize) {
    params.sanitize = true;
    params.removeHs = true;
  } else {
    params.sanitize = false;
    params.removeHs = false;
  }
  return RDKit::v2::SmilesParse::MolFromSmiles(smi, params).release();
};

inline RWMol *SmartsToMol(const std::string &sma,
                          const SmartsParserParams &ps) {
  RDKit::v2::SmilesParse::SmartsParserParams v2ps;
  v2ps.debugParse = ps.debugParse;
  if (ps.replacements) {
    v2ps.replacements = *ps.replacements;
  }
  v2ps.allowCXSMILES = ps.allowCXSMILES;
  v2ps.strictCXSMILES = ps.strictCXSMILES;
  v2ps.parseName = ps.parseName;
  v2ps.mergeHs = ps.mergeHs;
  v2ps.skipCleanup = ps.skipCleanup;

  return RDKit::v2::SmilesParse::MolFromSmarts(sma, v2ps).release();
}

//! Construct a molecule from a SMARTS string
/*!
 \param sma           the SMARTS to convert
 \param debugParse    toggles verbose debugging information from the parser
 \param mergeHs       toggles merging H atoms in the SMARTS into neighboring
 atoms
 \param replacements  a string->string map of replacement strings.
                      \see SmilesToMol for more information about replacements

 \return a pointer to the new molecule; the caller is responsible for free'ing
 this.
 */
inline RWMol *SmartsToMol(
    const std::string &sma, int debugParse = 0, bool mergeHs = false,
    std::map<std::string, std::string> *replacements = nullptr) {
  RDKit::v2::SmilesParse::SmartsParserParams ps;
  ps.debugParse = debugParse;
  ps.mergeHs = mergeHs;
  if (replacements) {
    ps.replacements = *replacements;
  }
  return RDKit::v2::SmilesParse::MolFromSmarts(sma, ps).release();
};

inline Atom *SmartsToAtom(const std::string &sma) {
  return RDKit::v2::SmilesParse::AtomFromSmarts(sma).release();
}
inline Bond *SmartsToBond(const std::string &sma) {
  return RDKit::v2::SmilesParse::BondFromSmarts(sma).release();
}
}  // namespace v1

inline std::unique_ptr<RDKit::RWMol> operator"" _smiles(const char *text,
                                                        size_t len) {
  std::string smi(text, len);
  try {
    return v2::SmilesParse::MolFromSmiles(smi);
  } catch (const RDKit::MolSanitizeException &) {
    return nullptr;
  }
}
inline std::unique_ptr<RDKit::RWMol> operator"" _smarts(const char *text,
                                                        size_t len) {
  std::string smi(text, len);
  return v2::SmilesParse::MolFromSmarts(smi);
}

}  // namespace RDKit

#endif
