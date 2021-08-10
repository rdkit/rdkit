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

// \brief construct a molecule from MDL mol data in a stream
/*!
 *   \param inStream - stream containing the data
 *   \param line     - current line number (used for error reporting)
 *   \param sanitize - toggles sanitization and stereochemistry
 *                     perception of the molecule
 *   \param removeHs - toggles removal of Hs from the molecule. H removal
 *                     is only done if the molecule is sanitized
 *   \param line     - current line number (used for error reporting)
 *   \param strictParsing - if set to false, the parser is more lax about
 * correctness of the contents.
 *
 */
inline RWMol *MolDataStreamToMol(std::istream *inStream, unsigned int &line,
                                 bool sanitize = true, bool removeHs = true,
                                 bool strictParsing = true) {
  return FileParsers::MolDataStreamToMol(inStream, line, sanitize, removeHs,
                                         strictParsing)
      .release();
}
// \overload
inline RWMol *MolDataStreamToMol(std::istream &inStream, unsigned int &line,
                                 bool sanitize = true, bool removeHs = true,
                                 bool strictParsing = true) {
  return FileParsers::MolDataStreamToMol(inStream, line, sanitize, removeHs,
                                         strictParsing)
      .release();
}
// \brief construct a molecule from an MDL mol block
/*!
 *   \param molBlock - string containing the mol block
 *   \param sanitize - toggles sanitization and stereochemistry
 *                     perception of the molecule
 *   \param removeHs - toggles removal of Hs from the molecule. H removal
 *                     is only done if the molecule is sanitized
 *   \param strictParsing - if set to false, the parser is more lax about
 * correctness of the contents.
 */
inline RWMol *MolBlockToMol(const std::string &molBlock, bool sanitize = true,
                            bool removeHs = true, bool strictParsing = true) {
  return FileParsers::MolBlockToMol(molBlock, sanitize, removeHs, strictParsing)
      .release();
}

// \brief construct a molecule from an MDL mol file
/*!
 *   \param fName    - string containing the file name
 *   \param sanitize - toggles sanitization and stereochemistry
 *                     perception of the molecule
 *   \param removeHs - toggles removal of Hs from the molecule. H removal
 *                     is only done if the molecule is sanitized
 *   \param strictParsing - if set to false, the parser is more lax about
 * correctness of the contents.
 */
inline RWMol *MolFileToMol(const std::string &fName, bool sanitize = true,
                           bool removeHs = true, bool strictParsing = true) {
  return FileParsers::MolFileToMol(fName, sanitize, removeHs, strictParsing)
      .release();
}
//-----
//  TPL handling:
//-----

//! \brief translate TPL data (BioCad format) into a multi-conf molecule
/*!
  \param inStream:      the stream from which to read
  \param line:          used to track the line number of errors
  \param sanitize:      toggles sanitization and stereochemistry
                        perception of the molecule
  \param skipFirstConf: according to the TPL format description, the atomic
                        coords in the atom-information block describe the first
                        conformation and the first conf block describes second
                        conformation. The CombiCode, on the other hand, writes
                        the first conformation data both to the atom-information
                        block and to the first conf block. We want to be able to
                        read CombiCode-style tpls, so we'll allow this
  mis-feature
                        to be parsed when this flag is set.
*/
inline RWMol *TPLDataStreamToMol(std::istream *inStream, unsigned int &line,
                                 bool sanitize = true,
                                 bool skipFirstConf = false) {
  return FileParsers::TPLDataStreamToMol(inStream, line, sanitize,
                                         skipFirstConf)
      .release();
}

//! \brief construct a multi-conf molecule from a TPL (BioCad format) file
/*!
  \param fName:         the name of the file from which to read
  \param sanitize:      toggles sanitization and stereochemistry
                        perception of the molecule
  \param skipFirstConf: according to the TPL format description, the atomic
                        coords in the atom-information block describe the first
                        conformation and the first conf block describes second
                        conformation. The CombiCode, on the other hand, writes
                        the first conformation data both to the atom-information
                        block and to the first conf block. We want to be able to
                        read CombiCode-style tpls, so we'll allow this
  mis-feature
                        to be parsed when this flag is set.
*/
inline RWMol *TPLFileToMol(const std::string &fName, bool sanitize = true,
                           bool skipFirstConf = false) {
  return FileParsers::TPLFileToMol(fName, sanitize, skipFirstConf).release();
}

//-----
//  MOL2 handling
//-----

// \brief construct a molecule from a Tripos mol2 file
/*!
 *
 *   \param fName    - string containing the file name
 *   \param sanitize - toggles sanitization of the molecule
 *   \param removeHs - toggles removal of Hs from the molecule. H removal
 *                     is only done if the molecule is sanitized
 *   \param variant  - the atom type definitions to use
 *   \param cleanupSubstructures - toggles recognition and cleanup of common
 *                                 substructures
 */
inline RWMol *Mol2FileToMol(const std::string &fName, bool sanitize = true,
                            bool removeHs = true,
                            Mol2Type variant = Mol2Type::CORINA,
                            bool cleanupSubstructures = true) {
  return FileParsers::Mol2FileToMol(fName, sanitize, removeHs, variant,
                                    cleanupSubstructures)
      .release();
}

// \brief construct a molecule from Tripos mol2 data in a stream
/*!
 *   \param inStream - stream containing the data
 *   \param sanitize - toggles sanitization of the molecule
 *   \param removeHs - toggles removal of Hs from the molecule. H removal
 *                     is only done if the molecule is sanitized
 *   \param variant  - the atom type definitions to use
 *   \param cleanupSubstructures - toggles recognition and cleanup of common
 *                                 substructures
 */
inline RWMol *Mol2DataStreamToMol(std::istream *inStream, bool sanitize = true,
                                  bool removeHs = true,
                                  Mol2Type variant = Mol2Type::CORINA,
                                  bool cleanupSubstructures = true) {
  return FileParsers::Mol2DataStreamToMol(inStream, sanitize, removeHs, variant,
                                          cleanupSubstructures)
      .release();
}
// \overload
inline RWMol *Mol2DataStreamToMol(std::istream &inStream, bool sanitize = true,
                                  bool removeHs = true,
                                  Mol2Type variant = Mol2Type::CORINA,
                                  bool cleanupSubstructures = true) {
  return FileParsers::Mol2DataStreamToMol(inStream, sanitize, removeHs, variant,
                                          cleanupSubstructures)
      .release();
}

// \brief construct a molecule from a Tripos mol2 block
/*!
 *   \param molBlock - string containing the mol block
 *   \param sanitize - toggles sanitization of the molecule
 *   \param removeHs - toggles removal of Hs from the molecule. H removal
 *                     is only done if the molecule is sanitized
 *   \param variant  - the atom type definitions to use
 *   \param cleanupSubstructures - toggles recognition and cleanup of common
 *                                 substructures
 */
inline RWMol *Mol2BlockToMol(const std::string &molBlock, bool sanitize = true,
                             bool removeHs = true,
                             Mol2Type variant = Mol2Type::CORINA,
                             bool cleanupSubstructures = true) {
  return FileParsers::Mol2BlockToMol(molBlock, sanitize, removeHs, variant,
                                     cleanupSubstructures)
      .release();
}

inline RWMol *PDBBlockToMol(const char *str, bool sanitize = true,
                            bool removeHs = true, unsigned int flavor = 0,
                            bool proximityBonding = true) {
  return FileParsers::PDBBlockToMol(str, sanitize, removeHs, flavor,
                                    proximityBonding)
      .release();
}

inline RWMol *PDBBlockToMol(const std::string &str, bool sanitize = true,
                            bool removeHs = true, unsigned int flavor = 0,
                            bool proximityBonding = true) {
  return FileParsers::PDBBlockToMol(str, sanitize, removeHs, flavor,
                                    proximityBonding)
      .release();
}
inline RWMol *PDBDataStreamToMol(std::istream *inStream, bool sanitize = true,
                                 bool removeHs = true, unsigned int flavor = 0,
                                 bool proximityBonding = true) {
  return FileParsers::PDBDataStreamToMol(inStream, sanitize, removeHs, flavor,
                                         proximityBonding)
      .release();
}
inline RWMol *PDBDataStreamToMol(std::istream &inStream, bool sanitize = true,
                                 bool removeHs = true, unsigned int flavor = 0,
                                 bool proximityBonding = true) {
  return FileParsers::PDBDataStreamToMol(inStream, sanitize, removeHs, flavor,
                                         proximityBonding)
      .release();
}
inline RWMol *PDBFileToMol(const std::string &fname, bool sanitize = true,
                           bool removeHs = true, unsigned int flavor = 0,
                           bool proximityBonding = true) {
  return FileParsers::PDBFileToMol(fname, sanitize, removeHs, flavor,
                                   proximityBonding)
      .release();
}

// \brief reads a molecule from the metadata in an RDKit-generated SVG file
/*!
 *   \param svg      - string containing the SVG
 *   \param sanitize - toggles sanitization of the molecule
 *   \param removeHs - toggles removal of Hs from the molecule. H removal
 *                     is only done if the molecule is sanitized
 *
 *   **NOTE** This functionality should be considered beta.
 */
inline RWMol *RDKitSVGToMol(const std::string &svg, bool sanitize = true,
                            bool removeHs = true) {
  return FileParsers::RDKitSVGToMol(svg, sanitize, removeHs).release();
}
/*! \overload
 */
inline RWMol *RDKitSVGToMol(std::istream *instream, bool sanitize = true,
                            bool removeHs = true) {
  return FileParsers::RDKitSVGToMol(instream, sanitize, removeHs).release();
}

}  // namespace RDKit

#endif
