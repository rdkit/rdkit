//
//  Copyright (C) 2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_PNGPARSER_H
#define RD_PNGPARSER_H

#include <RDGeneral/types.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <boost/format.hpp>

#include <string>
#include <fstream>
#include <sstream>

namespace RDKit {

namespace PNGData {
RDKIT_FILEPARSERS_EXPORT extern const std::string smilesTag;
RDKIT_FILEPARSERS_EXPORT extern const std::string molTag;
RDKIT_FILEPARSERS_EXPORT extern const std::string pklTag;
}  // namespace PNGData

RDKIT_FILEPARSERS_EXPORT std::vector<std::pair<std::string, std::string>>
PNGStreamToMetadata(std::istream &inStream);
inline std::vector<std::pair<std::string, std::string>> PNGFileToMetadata(
    const std::string &fname) {
  std::ifstream inStream(fname.c_str(), std::ios::binary);
  if (!inStream || (inStream.bad())) {
    throw BadFileException((boost::format("Bad input file %s") % fname).str());
  }
  return PNGStreamToMetadata(inStream);
}
inline std::vector<std::pair<std::string, std::string>> PNGStringToMetadata(
    const std::string &data) {
  std::stringstream inStream(data);
  return PNGStreamToMetadata(inStream);
}

RDKIT_FILEPARSERS_EXPORT ROMol *PNGStreamToMol(
    std::istream &inStream,
    const SmilesParserParams &params = SmilesParserParams());
inline ROMol *PNGFileToMol(
    const std::string &fname,
    const SmilesParserParams &params = SmilesParserParams()) {
  std::ifstream inStream(fname.c_str(), std::ios::binary);
  if (!inStream || (inStream.bad())) {
    throw BadFileException((boost::format("Bad input file %s") % fname).str());
  }
  return PNGStreamToMol(inStream, params);
}
inline ROMol *PNGStringToMol(
    const std::string &data,
    const SmilesParserParams &params = SmilesParserParams()) {
  std::stringstream inStream(data);
  return PNGStreamToMol(inStream, params);
}

RDKIT_FILEPARSERS_EXPORT std::vector<std::unique_ptr<ROMol>> PNGStreamToMols(
    std::istream &inStream, const std::string &tagToUse = PNGData::pklTag,
    const SmilesParserParams &params = SmilesParserParams());
inline std::vector<std::unique_ptr<ROMol>> PNGFileToMols(
    const std::string &fname, const std::string &tagToUse = PNGData::pklTag,
    const SmilesParserParams &params = SmilesParserParams()) {
  std::ifstream inStream(fname.c_str(), std::ios::binary);
  if (!inStream || (inStream.bad())) {
    throw BadFileException((boost::format("Bad input file %s") % fname).str());
  }
  return PNGStreamToMols(inStream, tagToUse, params);
}
inline std::vector<std::unique_ptr<ROMol>> PNGStringToMols(
    const std::string &data, const std::string &tagToUse = PNGData::pklTag,
    const SmilesParserParams &params = SmilesParserParams()) {
  std::stringstream inStream(data);
  return PNGStreamToMols(inStream, tagToUse, params);
}

//! The compressed flag is ignored if the RDKit is not built with
//! boost::iostreams support
RDKIT_FILEPARSERS_EXPORT std::string addMetadataToPNGStream(
    std::istream &iStream,
    const std::vector<std::pair<std::string, std::string>> &metadata,
    bool compressed = true);
inline std::string addMetadataToPNGString(
    const std::string &pngString,
    const std::vector<std::pair<std::string, std::string>> &metadata,
    bool compressed = true) {
  std::stringstream inStream(pngString);
  return addMetadataToPNGStream(inStream, metadata, compressed);
}
inline std::string addMetadataToPNGFile(
    const std::string &fname,
    const std::vector<std::pair<std::string, std::string>> &metadata,
    bool compressed = true) {
  std::ifstream inStream(fname.c_str(), std::ios::binary);
  return addMetadataToPNGStream(inStream, metadata, compressed);
}

RDKIT_FILEPARSERS_EXPORT std::string addMolToPNGStream(
    const ROMol &mol, std::istream &iStream, bool includePkl = true,
    bool includeSmiles = true, bool includeMol = false);
inline std::string addMolToPNGString(const ROMol &mol,
                                     const std::string &pngString,
                                     bool includePkl = true,
                                     bool includeSmiles = true,
                                     bool includeMol = false) {
  std::stringstream inStream(pngString);
  return addMolToPNGStream(mol, inStream, includePkl, includeSmiles,
                           includeMol);
}
inline std::string addMolToPNGFile(const ROMol &mol, const std::string &fname,
                                   bool includePkl = true,
                                   bool includeSmiles = true,
                                   bool includeMol = false) {
  std::ifstream inStream(fname.c_str(), std::ios::binary);
  return addMolToPNGStream(mol, inStream, includePkl, includeSmiles,
                           includeMol);
}

}  // namespace RDKit

#endif
