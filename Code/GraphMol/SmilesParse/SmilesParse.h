//
//  Copyright (C) 2001-2016 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_SMILESPARSE_H_
#define _RD_SMILESPARSE_H_

#include <GraphMol/RWMol.h>
#include <GraphMol/SanitException.h>
#include <string>
#include <exception>
#include <map>

namespace RDKit {

struct RDKIT_SMILESPARSE_EXPORT SmilesParserParams {
  int debugParse;
  bool sanitize;
  std::map<std::string, std::string> *replacements;
  bool allowCXSMILES;
  bool parseName;
  bool removeHs;
  SmilesParserParams()
      : debugParse(0),
        sanitize(true),
        replacements(NULL),
        allowCXSMILES(true),
        parseName(false),
        removeHs(true){};
};
RDKIT_SMILESPARSE_EXPORT RWMol *SmilesToMol(const std::string &smi,
                                            const SmilesParserParams &params);

RDKIT_SMILESPARSE_EXPORT Atom *SmilesToAtom(const std::string &smi);
RDKIT_SMILESPARSE_EXPORT Bond *SmilesToBond(const std::string &smi);

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
 the string no longer changes. It is the responsiblity of the caller to make
 sure
 that substitutions results in legal and sensible SMILES.

 Examples of substitutions:
 \code
   CC{Q}C with {"{Q}":"OCCO"} -> CCOCCOC
   C{A}C{Q}C with {"{Q}":"OCCO", "{A}":"C1(CC1)"} -> CC1(CC1)COCCOC
   C{A}C{Q}C with {"{Q}":"{X}CC{X}", "{A}":"C1CC1", "{X}":"N"} -> CC1CC1CCNCCNC
 \endcode

 */
inline RWMol *SmilesToMol(
    const std::string &smi, int debugParse = 0, bool sanitize = true,
    std::map<std::string, std::string> *replacements = 0) {
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
  return SmilesToMol(smi, params);
};

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
RDKIT_SMILESPARSE_EXPORT RWMol *SmartsToMol(
    const std::string &sma, int debugParse = 0, bool mergeHs = false,
    std::map<std::string, std::string> *replacements = 0);

RDKIT_SMILESPARSE_EXPORT Atom *SmartsToAtom(const std::string &sma);
RDKIT_SMILESPARSE_EXPORT Bond *SmartsToBond(const std::string &sma);

class RDKIT_SMILESPARSE_EXPORT SmilesParseException : public std::exception {
 public:
  SmilesParseException(const char *msg) : _msg(msg){};
  SmilesParseException(const std::string msg) : _msg(msg){};
  const char *message() const { return _msg.c_str(); };
  ~SmilesParseException() throw(){};

 private:
  std::string _msg;
};

inline std::unique_ptr<RDKit::RWMol> operator"" _smiles(const char *text,
                                                        size_t len) {
  std::string smi(text, len);
  RWMol *ptr = nullptr;
  try {
    ptr = SmilesToMol(smi);
  } catch (const RDKit::MolSanitizeException &e) {
    ptr = nullptr;
  }
  return std::unique_ptr<RWMol>(ptr);
}
inline std::unique_ptr<RDKit::RWMol> operator"" _smarts(const char *text,
                                                        size_t len) {
  std::string smi(text, len);
  RWMol *ptr = nullptr;
  try {
    ptr = SmartsToMol(smi);
  } catch (const RDKit::MolSanitizeException &e) {
    ptr = nullptr;
  }
  return std::unique_ptr<RWMol>(ptr);
}

}  // namespace RDKit

#endif
