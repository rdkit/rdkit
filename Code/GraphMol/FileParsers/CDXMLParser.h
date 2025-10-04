//
//  Copyright (c) 2022 Brian P Kelley
//  All rights reserved.
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_CDXML_FILEPARSERS_H
#define RD_CDXML_FILEPARSERS_H

#include <RDGeneral/types.h>
#include <string>
#include <vector>

namespace RDKit {
class RWMol;

namespace v2 {
namespace CDXMLParser {

enum class CDXMLFormat {
  CDXML = 0,
  CDX = 1,
  Auto = 2
};

//! \brief Returns true if the RDKit was build with ChemDraw CDX support
RDKIT_FILEPARSERS_EXPORT bool hasChemDrawCDXSupport();

struct RDKIT_FILEPARSERS_EXPORT CDXMLParserParams {
  bool sanitize = true;
  bool removeHs = true;
  CDXMLFormat format = CDXMLFormat::Auto;

  CDXMLParserParams() = default;
  CDXMLParserParams(bool sanitize, bool removeHs, CDXMLFormat format)
      : sanitize(sanitize), removeHs(removeHs), format(format) {}
};

//! \brief construct molecules from a CDXML file
//! The RDKit is optionally built with the Revvity ChemDraw parser
//! If this is available, CDX and CDXML can be read, see CDXMLParserParams
//!   Note that the CDXML format is large and complex, the RDKit doesn't
//!   support full functionality, just the base ones required for molecule and
//!   reaction parsing.
//! Note: If the ChemDraw extensions are available, this auto detects between
//!  CDXML and CDX
/*!
 *   \param inStream - string containing the mol block
 *   \param params - parameters controlling the parsing and post-processing
 */
RDKIT_FILEPARSERS_EXPORT std::vector<std::unique_ptr<RWMol>>
MolsFromCDXMLDataStream(std::istream &inStream,
                        const CDXMLParserParams &params = CDXMLParserParams());
//! \brief construct molecules from a CDXML file
//! The RDKit is optionally built with the Revvity ChemDraw parser
//! If this is available, CDX and CDXML can be read, see CDXMLParserParams
//!   Note that the CDXML format is large and complex, the RDKit doesn't
//!   support full functionality, just the base ones required for molecule and
//!   reaction parsing.
/*!
 *   \param fileName - cdxml fileName
 *   \param params - parameters controlling the parsing and post-processing
 */
RDKIT_FILEPARSERS_EXPORT std::vector<std::unique_ptr<RWMol>> MolsFromCDXMLFile(
    const std::string &filename,
    const CDXMLParserParams &params = CDXMLParserParams(true, true,
                                                        CDXMLFormat::Auto));

//! \brief construct molecules from a CDXML block
//! The RDKit is optionally built with the Revvity ChemDraw parser
//! If this is available, CDX and CDXML can be read, see CDXMLParserParams
//!   Note that the CDXML format is large and complex, the RDKit doesn't
//!   support full functionality, just the base ones required for molecule and
//!   reaction parsing.
//! Note: If the ChemDraw extensions are available,
//!   CDXMLFormat::Auto attempts to see if the input string is CDXML or CDX
/*!
 *   \param cdxml - string containing the mol block
 *   \param params - parameters controlling the parsing and post-processing
 */
RDKIT_FILEPARSERS_EXPORT std::vector<std::unique_ptr<RWMol>> MolsFromCDXML(
    const std::string &cdxml,
    const CDXMLParserParams &params =
        CDXMLParserParams(true, true, v2::CDXMLParser::CDXMLFormat::Auto));
}  // namespace CDXMLParser
}  // namespace v2

inline namespace v1 {

//! \brief construct molecules from a CDXML file
//! Note that the CDXML format is large and complex, the RDKit doesn't support
//!  full functionality, just the base ones required for molecule and
//!  reaction parsing.
//! Note: If the ChemDraw extensions are available, this auto detects between
//!  CDXML and CDX
/*!
 *   \param inStream - string containing the mol block
 *   \param sanitize - toggles sanitization and stereochemistry
 *                     perception of the molecule
 *   \param removeHs - toggles removal of Hs from the molecule. H removal
 *                     is only done if the molecule is sanitized
 * correctness of the contents.
 */
inline std::vector<std::unique_ptr<RWMol>> CDXMLDataStreamToMols(
    std::istream &inStream, bool sanitize = true, bool removeHs = true) {
  v2::CDXMLParser::CDXMLParserParams params(sanitize, removeHs,
                                            v2::CDXMLParser::CDXMLFormat::Auto);
  return v2::CDXMLParser::MolsFromCDXMLDataStream(inStream, params);
}

//! \brief construct molecules from a CDXML file
//! Note that the CDXML format is large and complex, the RDKit doesn't support
//!  full functionality, just the base ones required for molecule and
//!  reaction parsing.
//! Note: If the ChemDraw extensions are available,
//!   This function uses the file extension to determine the file type, .cdx or
//!   .cdxml If not, it defaults to CDXML
/*!
 *   \param fileName - cdxml fileName
 *   \param sanitize - toggles sanitization and stereochemistry
 *                     perception of the molecule
 *   \param removeHs - toggles removal of Hs from the molecule. H removal
 *                     is only done if the molecule is sanitized
 * correctness of the contents.
 */
inline std::vector<std::unique_ptr<RWMol>> CDXMLFileToMols(
    const std::string &filename, bool sanitize = true, bool removeHs = true) {
  v2::CDXMLParser::CDXMLParserParams params;
  params.sanitize = sanitize;
  params.removeHs = removeHs;
  params.format = v2::CDXMLParser::CDXMLFormat::Auto;
  return v2::CDXMLParser::MolsFromCDXMLFile(filename, params);
}

//! \brief construct molecules from a CDXML block
//! Note that the CDXML format is large and complex, the RDKit doesn't support
//!  full functionality, just the base ones required for molecule and
//!  reaction parsing.
//! Note: to parse CDX files see the CDXParserParams variant of this function
/*!
 *   \param cdxml - string containing the mol block
 *   \param sanitize - toggles sanitization and stereochemistry
 *                     perception of the molecule
 *   \param removeHs - toggles removal of Hs from the molecule. H removal
 *                     is only done if the molecule is sanitized
 * correctness of the contents.
 */
inline std::vector<std::unique_ptr<RWMol>> CDXMLToMols(const std::string &cdxml,
                                                       bool sanitize = true,
                                                       bool removeHs = true) {
  v2::CDXMLParser::CDXMLParserParams params;
  params.sanitize = sanitize;
  params.removeHs = removeHs;
  params.format = v2::CDXMLParser::CDXMLFormat::Auto;
  return v2::CDXMLParser::MolsFromCDXML(cdxml, params);
}
}  // namespace v1

}  // namespace RDKit
#endif  // RD_CDXML_FILEPARSERS_H
