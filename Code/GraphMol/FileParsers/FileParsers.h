//
//  Copyright (C) 2002-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_FILEPARSERS_H
#define RD_FILEPARSERS_H

#include "FileParsersv2.h"

namespace RDKit {

//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *MolDataStreamToMol(std::istream *inStream, unsigned int &line,
                                 bool sanitize = true, bool removeHs = true,
                                 bool strictParsing = true) {
  return FileParsers::MolDataStreamToMol(inStream, line, sanitize, removeHs,
                                         strictParsing)
      .release();
}
//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *MolDataStreamToMol(std::istream &inStream, unsigned int &line,
                                 bool sanitize = true, bool removeHs = true,
                                 bool strictParsing = true) {
  return FileParsers::MolDataStreamToMol(&inStream, line, sanitize, removeHs,
                                         strictParsing)
      .release();
}
//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *MolBlockToMol(const std::string &molBlock, bool sanitize = true,
                            bool removeHs = true, bool strictParsing = true) {
  return FileParsers::MolBlockToMol(molBlock, sanitize, removeHs, strictParsing)
      .release();
}

//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *MolFileToMol(const std::string &fName, bool sanitize = true,
                           bool removeHs = true, bool strictParsing = true) {
  return FileParsers::MolFileToMol(fName, sanitize, removeHs, strictParsing)
      .release();
}
//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *TPLDataStreamToMol(std::istream *inStream, unsigned int &line,
                                 bool sanitize = true,
                                 bool skipFirstConf = false) {
  return FileParsers::TPLDataStreamToMol(inStream, line, sanitize,
                                         skipFirstConf)
      .release();
}

//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *TPLFileToMol(const std::string &fName, bool sanitize = true,
                           bool skipFirstConf = false) {
  return FileParsers::TPLFileToMol(fName, sanitize, skipFirstConf).release();
}

//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *Mol2FileToMol(const std::string &fName, bool sanitize = true,
                            bool removeHs = true,
                            Mol2Type variant = Mol2Type::CORINA,
                            bool cleanupSubstructures = true) {
  return FileParsers::Mol2FileToMol(fName, sanitize, removeHs, variant,
                                    cleanupSubstructures)
      .release();
}

//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *Mol2DataStreamToMol(std::istream *inStream, bool sanitize = true,
                                  bool removeHs = true,
                                  Mol2Type variant = Mol2Type::CORINA,
                                  bool cleanupSubstructures = true) {
  return FileParsers::Mol2DataStreamToMol(inStream, sanitize, removeHs, variant,
                                          cleanupSubstructures)
      .release();
}
//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *Mol2DataStreamToMol(std::istream &inStream, bool sanitize = true,
                                  bool removeHs = true,
                                  Mol2Type variant = Mol2Type::CORINA,
                                  bool cleanupSubstructures = true) {
  return FileParsers::Mol2DataStreamToMol(&inStream, sanitize, removeHs,
                                          variant, cleanupSubstructures)
      .release();
}

//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *Mol2BlockToMol(const std::string &molBlock, bool sanitize = true,
                             bool removeHs = true,
                             Mol2Type variant = Mol2Type::CORINA,
                             bool cleanupSubstructures = true) {
  return FileParsers::Mol2BlockToMol(molBlock, sanitize, removeHs, variant,
                                     cleanupSubstructures)
      .release();
}

//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *PDBBlockToMol(const char *str, bool sanitize = true,
                            bool removeHs = true, unsigned int flavor = 0,
                            bool proximityBonding = true) {
  return FileParsers::PDBBlockToMol(str, sanitize, removeHs, flavor,
                                    proximityBonding)
      .release();
}

//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *PDBBlockToMol(const std::string &str, bool sanitize = true,
                            bool removeHs = true, unsigned int flavor = 0,
                            bool proximityBonding = true) {
  return FileParsers::PDBBlockToMol(str, sanitize, removeHs, flavor,
                                    proximityBonding)
      .release();
}
//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *PDBDataStreamToMol(std::istream *inStream, bool sanitize = true,
                                 bool removeHs = true, unsigned int flavor = 0,
                                 bool proximityBonding = true) {
  return FileParsers::PDBDataStreamToMol(inStream, sanitize, removeHs, flavor,
                                         proximityBonding)
      .release();
}
//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *PDBDataStreamToMol(std::istream &inStream, bool sanitize = true,
                                 bool removeHs = true, unsigned int flavor = 0,
                                 bool proximityBonding = true) {
  return FileParsers::PDBDataStreamToMol(&inStream, sanitize, removeHs, flavor,
                                         proximityBonding)
      .release();
}
//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *PDBFileToMol(const std::string &fname, bool sanitize = true,
                           bool removeHs = true, unsigned int flavor = 0,
                           bool proximityBonding = true) {
  return FileParsers::PDBFileToMol(fname, sanitize, removeHs, flavor,
                                   proximityBonding)
      .release();
}

//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *RDKitSVGToMol(const std::string &svg, bool sanitize = true,
                            bool removeHs = true) {
  return FileParsers::RDKitSVGToMol(svg, sanitize, removeHs).release();
}
//! Backwards-compatibility function. Please use the version in the SmilesParser
/// namespace instead
inline RWMol *RDKitSVGToMol(std::istream *instream, bool sanitize = true,
                            bool removeHs = true) {
  return FileParsers::RDKitSVGToMol(instream, sanitize, removeHs).release();
}

}  // namespace RDKit

#endif
