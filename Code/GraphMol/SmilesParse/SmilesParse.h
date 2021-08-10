//
//  Copyright (C) 2001-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_SMILESPARSE_H
#define RD_SMILESPARSE_H

#include "SmilesParsev2.h"

namespace RDKit {

inline RWMol *SmilesToMol(const std::string &smi,
                          const SmilesParserParams &params) {
  return SmilesParser::SmilesToMol(smi, params).release();
};

inline Atom *SmilesToAtom(const std::string &smi) {
  return SmilesParser::SmilesToAtom(smi).release();
};
inline Bond *SmilesToBond(const std::string &smi) {
  return SmilesParser::SmilesToBond(smi).release();
};

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
  SmilesParserParams params;
  params.debugParse = debugParse;
  params.replacements = replacements;
  if (sanitize) {
    params.sanitize = true;
    params.removeHs = true;
  } else {
    params.sanitize = false;
    params.removeHs = false;
  }
  return SmilesParser::SmilesToMol(smi, params).release();
};

inline RWMol *SmartsToMol(std::string sma, SmartsParserParams ps) {
  return SmilesParser::SmartsToMol(sma, ps).release();
};
;

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
  SmartsParserParams ps;
  ps.debugParse = debugParse;
  ps.mergeHs = mergeHs;
  ps.replacements = replacements;
  return SmilesParser::SmartsToMol(sma, ps).release();
}

inline Atom *SmartsToAtom(const std::string &sma) {
  return SmilesParser::SmartsToAtom(sma).release();
}
inline Bond *SmartsToBond(const std::string &sma) {
  return SmilesParser::SmartsToBond(sma).release();
}

}  // namespace RDKit

#endif
