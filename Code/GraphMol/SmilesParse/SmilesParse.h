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

//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *SmilesToMol(const std::string &smi,
                          const SmilesParserParams &params) {
  return SmilesParser::SmilesToMol(smi, params).release();
}
//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline Atom *SmilesToAtom(const std::string &smi) {
  return SmilesParser::SmilesToAtom(smi).release();
}
//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline Bond *SmilesToBond(const std::string &smi) {
  return SmilesParser::SmilesToBond(smi).release();
}
//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
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
}
//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *SmartsToMol(std::string sma, SmartsParserParams ps) {
  return SmilesParser::SmartsToMol(sma, ps).release();
}
//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *SmartsToMol(
    const std::string &sma, int debugParse = 0, bool mergeHs = false,
    std::map<std::string, std::string> *replacements = nullptr) {
  SmartsParserParams ps;
  ps.debugParse = debugParse;
  ps.mergeHs = mergeHs;
  ps.replacements = replacements;
  return SmilesParser::SmartsToMol(sma, ps).release();
}
//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline Atom *SmartsToAtom(const std::string &sma) {
  return SmilesParser::SmartsToAtom(sma).release();
}
//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline Bond *SmartsToBond(const std::string &sma) {
  return SmilesParser::SmartsToBond(sma).release();
}

}  // namespace RDKit

#endif
