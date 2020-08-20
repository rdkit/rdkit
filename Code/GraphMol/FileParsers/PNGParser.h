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

//! Tags used for PNG metadata
namespace PNGData {
RDKIT_FILEPARSERS_EXPORT extern const std::string smilesTag;
RDKIT_FILEPARSERS_EXPORT extern const std::string molTag;
RDKIT_FILEPARSERS_EXPORT extern const std::string pklTag;
}  // namespace PNGData

//! \name metadata to/from PNG
//@{

//! \brief returns the metadata (tEXt and zTXt chunks) from PNG data
RDKIT_FILEPARSERS_EXPORT std::vector<std::pair<std::string, std::string>>
PNGStreamToMetadata(std::istream &inStream);

//! \brief returns the metadata (tEXt and zTXt chunks) from PNG data
inline std::vector<std::pair<std::string, std::string>> PNGFileToMetadata(
    const std::string &fname) {
  std::ifstream inStream(fname.c_str(), std::ios::binary);
  if (!inStream || (inStream.bad())) {
    throw BadFileException((boost::format("Bad input file %s") % fname).str());
  }
  return PNGStreamToMetadata(inStream);
}

//! \brief returns the metadata (tEXt and zTXt chunks) from PNG data
inline std::vector<std::pair<std::string, std::string>> PNGStringToMetadata(
    const std::string &data) {
  std::stringstream inStream(data);
  return PNGStreamToMetadata(inStream);
}

//! \brief adds metadata to a PNG stream.
//! The modified PNG data is returned.
/*!

The compressed flag is ignored if the RDKit is not built with
boost::iostreams support

*/
RDKIT_FILEPARSERS_EXPORT std::string addMetadataToPNGStream(
    std::istream &iStream,
    const std::vector<std::pair<std::string, std::string>> &metadata,
    bool compressed = true);

//! \brief adds metadata to a PNG string.
//! The modified PNG data is returned.
inline std::string addMetadataToPNGString(
    const std::string &pngString,
    const std::vector<std::pair<std::string, std::string>> &metadata,
    bool compressed = true) {
  std::stringstream inStream(pngString);
  return addMetadataToPNGStream(inStream, metadata, compressed);
}

//! \brief adds metadata to a PNG file.
//! The modified PNG data is returned.
inline std::string addMetadataToPNGFile(
    const std::string &fname,
    const std::vector<std::pair<std::string, std::string>> &metadata,
    bool compressed = true) {
  std::ifstream inStream(fname.c_str(), std::ios::binary);
  return addMetadataToPNGStream(inStream, metadata, compressed);
}
//@}

//! \name molecules to/from PNG
//@{

//! \brief constructs an ROMol from the metadata in a PNG stream
/*!

Looks through the metadata in the PNG to find the first tag that matches one of
the tags in \c RDKit::PNGData. A molecule is constructed from this chunk.

Throws a \c FileParseException if no suitable tag is found.

The caller is responsible for the returned pointer.

 */
RDKIT_FILEPARSERS_EXPORT ROMol *PNGStreamToMol(
    std::istream &inStream,
    const SmilesParserParams &params = SmilesParserParams());
//! \brief constructs an ROMol from the metadata in a PNG file.
//! See \c PNGStreamToMol() for more details.
inline ROMol *PNGFileToMol(
    const std::string &fname,
    const SmilesParserParams &params = SmilesParserParams()) {
  std::ifstream inStream(fname.c_str(), std::ios::binary);
  if (!inStream || (inStream.bad())) {
    throw BadFileException((boost::format("Bad input file %s") % fname).str());
  }
  return PNGStreamToMol(inStream, params);
}
//! \brief constructs an ROMol from the metadata in a PNG string.
//! See \c PNGStreamToMol() for more details.
inline ROMol *PNGStringToMol(
    const std::string &data,
    const SmilesParserParams &params = SmilesParserParams()) {
  std::stringstream inStream(data);
  return PNGStreamToMol(inStream, params);
}

//! \brief constructs a vector of ROMol from the metadata in a PNG stream
/*!

Looks through the metadata in the PNG to find tags that start with tagToUse
(must be one of the tags in \c RDKit::PNGData). The molecules constructed from
these data are returned.

 */
RDKIT_FILEPARSERS_EXPORT std::vector<std::unique_ptr<ROMol>> PNGStreamToMols(
    std::istream &inStream, const std::string &tagToUse = PNGData::pklTag,
    const SmilesParserParams &params = SmilesParserParams());
//! \brief constructs a vector of ROMol from the metadata in a PNG file.
//! See \c PNGStreamToMols() for more details.
inline std::vector<std::unique_ptr<ROMol>> PNGFileToMols(
    const std::string &fname, const std::string &tagToUse = PNGData::pklTag,
    const SmilesParserParams &params = SmilesParserParams()) {
  std::ifstream inStream(fname.c_str(), std::ios::binary);
  if (!inStream || (inStream.bad())) {
    throw BadFileException((boost::format("Bad input file %s") % fname).str());
  }
  return PNGStreamToMols(inStream, tagToUse, params);
}
//! \brief constructs a vector of ROMol from the metadata in a PNG string.
//! See \c PNGStreamToMols() for more details.
inline std::vector<std::unique_ptr<ROMol>> PNGStringToMols(
    const std::string &data, const std::string &tagToUse = PNGData::pklTag,
    const SmilesParserParams &params = SmilesParserParams()) {
  std::stringstream inStream(data);
  return PNGStreamToMols(inStream, tagToUse, params);
}

//! \brief adds metadata for an ROMol to the data from a PNG stream.
//! The modified PNG data is returned.
/*!

  \param mol            the molecule to add
  \param iStream        the stream to read from
  \param includePkl     include a molecule pickle
  \param includeSmiles  include CXSMILES for the molecule
  \param includeMol     include a mol block for the molecule

*/
RDKIT_FILEPARSERS_EXPORT std::string addMolToPNGStream(
    const ROMol &mol, std::istream &iStream, bool includePkl = true,
    bool includeSmiles = true, bool includeMol = false);

//! \brief adds metadata for an ROMol to a PNG string.
//! The modified PNG data is returned.
//! See \c addMolToPNGStream() for more details.
inline std::string addMolToPNGString(const ROMol &mol,
                                     const std::string &pngString,
                                     bool includePkl = true,
                                     bool includeSmiles = true,
                                     bool includeMol = false) {
  std::stringstream inStream(pngString);
  return addMolToPNGStream(mol, inStream, includePkl, includeSmiles,
                           includeMol);
}
//! \brief adds metadata for an ROMol to the data from a PNG file.
//! The modified PNG data is returned.
//! See \c addMolToPNGStream() for more details.
inline std::string addMolToPNGFile(const ROMol &mol, const std::string &fname,
                                   bool includePkl = true,
                                   bool includeSmiles = true,
                                   bool includeMol = false) {
  std::ifstream inStream(fname.c_str(), std::ios::binary);
  return addMolToPNGStream(mol, inStream, includePkl, includeSmiles,
                           includeMol);
}
//@}

}  // namespace RDKit

#endif
