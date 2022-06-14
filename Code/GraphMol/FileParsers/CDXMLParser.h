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
#ifndef _RD_CDXML_FILEPARSERS_H
#define _RD_CDXML_FILEPARSERS_H

#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>

#include <string>
#include <iostream>
#include <vector>
#include <exception>

#include <boost/shared_ptr.hpp>

namespace RDKit {
// \brief construct molecules from a CDXML file
/*!
 *   \param inStream - string containing the mol block
 *   \param sanitize - toggles sanitization and stereochemistry
 *                     perception of the molecule
 *   \param removeHs - toggles removal of Hs from the molecule. H removal
 *                     is only done if the molecule is sanitized
 * correctness of the contents.
 */
RDKIT_FILEPARSERS_EXPORT std::vector<std::unique_ptr<RWMol>> CDXMLToMols(
			       std::istream &inStream,
			       bool sanitize = true,
			       bool removeHs = true);
// \brief construct molecules from a CDXML file
/*!
 *   \param fileName - cdxml fileName
 *   \param sanitize - toggles sanitization and stereochemistry
 *                     perception of the molecule
 *   \param removeHs - toggles removal of Hs from the molecule. H removal
 *                     is only done if the molecule is sanitized
 * correctness of the contents.
 */
RDKIT_FILEPARSERS_EXPORT std::vector<std::unique_ptr<RWMol>> CDXMLFileToMols(
			       const std::string &filename,
			       bool sanitize = true,
			       bool removeHs = true);

// \brief construct molecules from a CDXML file
/*!
 *   \param cdxml - string containing the mol block
 *   \param sanitize - toggles sanitization and stereochemistry
 *                     perception of the molecule
 *   \param removeHs - toggles removal of Hs from the molecule. H removal
 *                     is only done if the molecule is sanitized
 * correctness of the contents.
 */
RDKIT_FILEPARSERS_EXPORT std::vector<std::unique_ptr<RWMol>> CDXMLStringToMols(
			       const std::string &cdxml,
			       bool sanitize = true,
			       bool removeHs = true);

}
#endif // _RD_CDXML_FILEPARSERS_H
