//
//  Copyright (C) 2002-2024 Greg Landrum and other RDKit contributors
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

#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileWriters.h>
#include "CDXMLParser.h"
#include <string>
#include <string_view>
#include <iostream>
#include <vector>
#include <exception>

#include <boost/shared_ptr.hpp>

namespace RDKit {

RDKIT_FILEPARSERS_EXPORT std::string strip(const std::string &orig);

namespace v2 {
namespace FileParsers {
class RDKIT_FILEPARSERS_EXPORT MolFileUnhandledFeatureException
    : public std::exception {
 public:
  //! construct with an error message
  explicit MolFileUnhandledFeatureException(const char *msg) : _msg(msg) {}
  //! construct with an error message
  explicit MolFileUnhandledFeatureException(const std::string msg)
      : _msg(msg) {}
  //! get the error message
  const char *what() const noexcept override { return _msg.c_str(); }
  ~MolFileUnhandledFeatureException() noexcept override = default;

 private:
  std::string _msg;
};

struct RDKIT_FILEPARSERS_EXPORT MolFileParserParams {
  bool sanitize = true;      /**< sanitize the molecule after building it */
  bool removeHs = true;      /**< remove Hs after constructing the molecule */
  bool strictParsing = true; /**< if set to false, the parser is more lax about
                                correctness of the contents. */
  bool expandAttachmentPoints =
      false; /**< toggle conversion of attachment points into dummy atoms */
};
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromMolDataStream(
    std::istream &inStream, unsigned int &line,
    const MolFileParserParams &params = MolFileParserParams());
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromMolBlock(
    const std::string &molBlock,
    const MolFileParserParams &params = MolFileParserParams());
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromMolFile(
    const std::string &fName,
    const MolFileParserParams &params = MolFileParserParams());

}  // namespace FileParsers
}  // namespace v2

inline namespace v1 {
using RDKit::v2::FileParsers::MolFileUnhandledFeatureException;
//-----
// mol files
//-----
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
  v2::FileParsers::MolFileParserParams ps;
  ps.sanitize = sanitize;
  ps.removeHs = removeHs;
  ps.strictParsing = strictParsing;
  return v2::FileParsers::MolFromMolDataStream(*inStream, line, ps).release();
};
// \overload
inline RWMol *MolDataStreamToMol(std::istream &inStream, unsigned int &line,
                                 bool sanitize = true, bool removeHs = true,
                                 bool strictParsing = true) {
  return MolDataStreamToMol(&inStream, line, sanitize, removeHs, strictParsing);
};
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
  v2::FileParsers::MolFileParserParams ps;
  ps.sanitize = sanitize;
  ps.removeHs = removeHs;
  ps.strictParsing = strictParsing;
  return v2::FileParsers::MolFromMolBlock(molBlock, ps).release();
};

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
  v2::FileParsers::MolFileParserParams ps;
  ps.sanitize = sanitize;
  ps.removeHs = removeHs;
  ps.strictParsing = strictParsing;
  return v2::FileParsers::MolFromMolFile(fName, ps).release();
};
}  // namespace v1

//-----
//  TPL handling:
//-----

namespace v2 {
namespace FileParsers {
struct RDKIT_FILEPARSERS_EXPORT TPLParserParams {
  bool sanitize = true; /**< sanitize the molecule after building it */
  bool skipFirstConf =
      false; /**< if set to true, the first conformer will be skipped */
};
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromTPLDataStream(
    std::istream &inStream, unsigned int &line,
    const TPLParserParams &params = TPLParserParams());
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromTPLFile(
    const std::string &fName,
    const TPLParserParams &params = TPLParserParams());

}  // namespace FileParsers
}  // namespace v2

inline namespace v1 {
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
  v2::FileParsers::TPLParserParams ps;
  ps.sanitize = sanitize;
  ps.skipFirstConf = skipFirstConf;
  return v2::FileParsers::MolFromTPLDataStream(*inStream, line, ps).release();
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
  v2::FileParsers::TPLParserParams ps;
  ps.sanitize = sanitize;
  ps.skipFirstConf = skipFirstConf;
  return v2::FileParsers::MolFromTPLFile(fName, ps).release();
}
}  // namespace v1

namespace v2 {
namespace FileParsers {

//-----
//  MOL2 handling
//-----

typedef enum {
  CORINA = 0  //!< supports output from Corina and some dbtranslate output
} Mol2Type;

struct Mol2ParserParams {
  bool sanitize = true; /**< sanitize the molecule after building it */
  bool removeHs = true; /**< remove Hs after constructing the molecule */
  Mol2Type variant = Mol2Type::CORINA; /**< the atom type definitions to use */
  bool cleanupSubstructures =
      true; /**< toggles recognition and cleanup of common substructures */
};

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromMol2DataStream(
    std::istream &inStream,
    const Mol2ParserParams &params = Mol2ParserParams());
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromMol2Block(
    const std::string &molBlock,
    const Mol2ParserParams &params = Mol2ParserParams());
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromMol2File(
    const std::string &fName,
    const Mol2ParserParams &params = Mol2ParserParams());

}  // namespace FileParsers
}  // namespace v2

inline namespace v1 {
using RDKit::v2::FileParsers::Mol2Type;

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
  v2::FileParsers::Mol2ParserParams ps;
  ps.sanitize = sanitize;
  ps.removeHs = removeHs;
  ps.variant = variant;
  ps.cleanupSubstructures = cleanupSubstructures;
  return v2::FileParsers::MolFromMol2File(fName, ps).release();
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
inline RWMol *Mol2DataStreamToMol(std::istream &inStream, bool sanitize = true,
                                  bool removeHs = true,
                                  Mol2Type variant = Mol2Type::CORINA,
                                  bool cleanupSubstructures = true) {
  v2::FileParsers::Mol2ParserParams ps;
  ps.sanitize = sanitize;
  ps.removeHs = removeHs;
  ps.variant = variant;
  ps.cleanupSubstructures = cleanupSubstructures;
  return v2::FileParsers::MolFromMol2DataStream(inStream, ps).release();
}
// \overload
inline RWMol *Mol2DataStreamToMol(std::istream *inStream, bool sanitize = true,
                                  bool removeHs = true,
                                  Mol2Type variant = Mol2Type::CORINA,
                                  bool cleanupSubstructures = true) {
  return Mol2DataStreamToMol(*inStream, sanitize, removeHs, variant,
                             cleanupSubstructures);
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
  v2::FileParsers::Mol2ParserParams ps;
  ps.sanitize = sanitize;
  ps.removeHs = removeHs;
  ps.variant = variant;
  ps.cleanupSubstructures = cleanupSubstructures;
  return v2::FileParsers::MolFromMol2Block(molBlock, ps).release();
}
}  // namespace v1

namespace v2 {
namespace FileParsers {

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromXYZDataStream(
    std::istream &inStream);
// \brief construct a molecule from an xyz block
/*!
 *   \param xyzBlock    - string containing the xyz block
 */
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromXYZBlock(
    const std::string &xyzBlock);
// \brief construct a molecule from an xyz file
/*!
 *   \param fName    - string containing the file name
 */
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromXYZFile(
    const std::string &fName);
}  // namespace FileParsers
}  // namespace v2
inline namespace v1 {
inline RWMol *XYZDataStreamToMol(std::istream &inStream) {
  return v2::FileParsers::MolFromXYZDataStream(inStream).release();
}
// \brief construct a molecule from an xyz block
/*!
 *   \param xyzBlock    - string containing the xyz block
 */
inline RWMol *XYZBlockToMol(const std::string &xyzBlock) {
  return v2::FileParsers::MolFromXYZBlock(xyzBlock).release();
}
// \brief construct a molecule from an xyz file
/*!
 *   \param fName    - string containing the file name
 */
inline RWMol *XYZFileToMol(const std::string &fName) {
  return v2::FileParsers::MolFromXYZFile(fName).release();
}

}  // namespace v1

namespace v2 {
namespace FileParsers {
struct RDKIT_FILEPARSERS_EXPORT PDBParserParams {
  bool sanitize = true; /**< sanitize the molecule after building it */
  bool removeHs = true; /**< remove Hs after constructing the molecule */
  bool proximityBonding = true; /**< if set to true, proximity bonding will be
                                   performed */
  unsigned int flavor = 0;      /**< flavor to use */
};

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromPDBDataStream(
    std::istream &inStream, const PDBParserParams &params = PDBParserParams());
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromPDBFile(
    const std::string &fname,
    const PDBParserParams &params = PDBParserParams());
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromPDBBlock(
    const std::string &str, const PDBParserParams &params = PDBParserParams());
}  // namespace FileParsers
}  // namespace v2

inline namespace v1 {
using RDKit::v2::FileParsers::PDBParserParams;
inline RWMol *PDBBlockToMol(const std::string &str, bool sanitize = true,
                            bool removeHs = true, unsigned int flavor = 0,
                            bool proximityBonding = true) {
  v2::FileParsers::PDBParserParams ps;
  ps.sanitize = sanitize;
  ps.removeHs = removeHs;
  ps.flavor = flavor;
  ps.proximityBonding = proximityBonding;
  return v2::FileParsers::MolFromPDBBlock(str, ps).release();
}
inline RWMol *PDBBlockToMol(const char *str, bool sanitize = true,
                            bool removeHs = true, unsigned int flavor = 0,
                            bool proximityBonding = true) {
  return PDBBlockToMol(std::string(str), sanitize, removeHs, flavor,
                       proximityBonding);
}
inline RWMol *PDBFileToMol(const std::string &fname, bool sanitize = true,
                           bool removeHs = true, unsigned int flavor = 0,
                           bool proximityBonding = true) {
  v2::FileParsers::PDBParserParams ps;
  ps.sanitize = sanitize;
  ps.removeHs = removeHs;
  ps.flavor = flavor;
  ps.proximityBonding = proximityBonding;
  return v2::FileParsers::MolFromPDBFile(fname, ps).release();
}
inline RWMol *PDBDataStreamToMol(std::istream &inStream, bool sanitize = true,
                                 bool removeHs = true, unsigned int flavor = 0,
                                 bool proximityBonding = true) {
  v2::FileParsers::PDBParserParams ps;
  ps.sanitize = sanitize;
  ps.removeHs = removeHs;
  ps.flavor = flavor;
  ps.proximityBonding = proximityBonding;
  return v2::FileParsers::MolFromPDBDataStream(inStream, ps).release();
}
inline RWMol *PDBDataStreamToMol(std::istream *inStream, bool sanitize = true,
                                 bool removeHs = true, unsigned int flavor = 0,
                                 bool proximityBonding = true) {
  return PDBDataStreamToMol(*inStream, sanitize, removeHs, flavor,
                            proximityBonding);
}
}  // namespace v1

// \brief reads a molecule from the metadata in an RDKit-generated SVG file
/*!
 *   \param svg      - string containing the SVG
 *   \param sanitize - toggles sanitization of the molecule
 *   \param removeHs - toggles removal of Hs from the molecule. H removal
 *                     is only done if the molecule is sanitized
 *
 *   **NOTE** This functionality should be considered beta.
 */
RDKIT_FILEPARSERS_EXPORT RWMol *RDKitSVGToMol(const std::string &svg,
                                              bool sanitize = true,
                                              bool removeHs = true);
/*! \overload
 */
RDKIT_FILEPARSERS_EXPORT RWMol *RDKitSVGToMol(std::istream *instream,
                                              bool sanitize = true,
                                              bool removeHs = true);

inline std::unique_ptr<RDKit::RWMol> operator"" _ctab(const char *text,
                                                      size_t len) {
  std::string data(text, len);
  try {
    return v2::FileParsers::MolFromMolBlock(data);
  } catch (const RDKit::MolSanitizeException &) {
    return nullptr;
  }
}
inline std::unique_ptr<RDKit::RWMol> operator"" _mol2(const char *text,
                                                      size_t len) {
  std::string data(text, len);
  try {
    return v2::FileParsers::MolFromMol2Block(data);
  } catch (const RDKit::MolSanitizeException &) {
    return nullptr;
  }
}

inline std::unique_ptr<RDKit::RWMol> operator"" _pdb(const char *text,
                                                     size_t len) {
  std::string data(text, len);
  try {
    return v2::FileParsers::MolFromPDBBlock(data);
  } catch (const RDKit::MolSanitizeException &) {
    return nullptr;
  }
}

}  // namespace RDKit

#endif
