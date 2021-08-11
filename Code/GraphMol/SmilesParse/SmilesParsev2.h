//
//  Copyright (C) 2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_SMILESPARSEv2_H
#define RD_SMILESPARSEv2_H

#include <GraphMol/RWMol.h>
#include <GraphMol/SanitException.h>
#include <string>
#include <exception>
#include <map>

namespace RDKit {

struct RDKIT_SMILESPARSE_EXPORT SmilesParserParams {
  int debugParse = 0;   /**< enable debugging in the SMILES parser*/
  bool sanitize = true; /**< sanitize the molecule after building it */
  std::map<std::string, std::string> *replacements =
      nullptr;               /**< allows SMILES "macros" */
  bool allowCXSMILES = true; /**< recognize and parse CXSMILES*/
  bool strictCXSMILES =
      true; /**< throw an exception if the CXSMILES parsing fails */
  bool parseName = false; /**< parse (and set) the molecule name as well */
  bool removeHs = true;   /**< remove Hs after constructing the molecule */
  bool useLegacyStereo =
      true; /**< use the legacy stereochemistry perception code */
};

struct RDKIT_SMILESPARSE_EXPORT SmartsParserParams {
  int debugParse = 0; /**< enable debugging in the SMARTS parser*/
  std::map<std::string, std::string> *replacements =
      nullptr;               /**< allows SMARTS "macros" */
  bool allowCXSMILES = true; /**< recognize and parse CXSMILES extensions */
  bool strictCXSMILES =
      true; /**< throw an exception if the CXSMILES parsing fails */
  bool parseName = false; /**< parse (and set) the molecule name as well */
  bool mergeHs =
      false; /**< toggles merging H atoms in the SMARTS into neighboring atoms*/
};

namespace SmilesParser {
//! Construct a molecule from a SMILES string
RDKIT_SMILESPARSE_EXPORT std::unique_ptr<RWMol> SmilesToMol(
    const std::string &smi,
    const SmilesParserParams &params = SmilesParserParams());

//! Construct an atom from a SMILES string
RDKIT_SMILESPARSE_EXPORT std::unique_ptr<Atom> SmilesToAtom(
    const std::string &smi);
//! Construct a bond from a SMILES string
RDKIT_SMILESPARSE_EXPORT std::unique_ptr<Bond> SmilesToBond(
    const std::string &smi);

//! Construct a molecule from a SMARTS string
RDKIT_SMILESPARSE_EXPORT std::unique_ptr<RWMol> SmartsToMol(
    const std::string &sma,
    const SmartsParserParams &params = SmartsParserParams());
//! Construct an atom from a SMARTS string
RDKIT_SMILESPARSE_EXPORT std::unique_ptr<Atom> SmartsToAtom(
    const std::string &sma);
//! Construct a bond from a SMARTS string
RDKIT_SMILESPARSE_EXPORT std::unique_ptr<Bond> SmartsToBond(
    const std::string &sma);

}  // namespace SmilesParser

class RDKIT_SMILESPARSE_EXPORT SmilesParseException : public std::exception {
 public:
  SmilesParseException(const char *msg) : _msg(msg) {}
  SmilesParseException(const std::string &msg) : _msg(msg) {}
  const char *what() const noexcept override { return _msg.c_str(); }
  ~SmilesParseException() noexcept override = default;

 private:
  std::string _msg;
};
RDKIT_SMILESPARSE_EXPORT std::unique_ptr<RDKit::RWMol> operator"" _smiles(
    const char *text, size_t len);
RDKIT_SMILESPARSE_EXPORT std::unique_ptr<RDKit::RWMol> operator"" _smarts(
    const char *text, size_t len);

}  // namespace RDKit

#endif
