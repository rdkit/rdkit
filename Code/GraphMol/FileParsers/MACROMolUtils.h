//
//  Copyright (C) 2010-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_MACROMOLRUTILS_H
#define RD_MACROMOLRUTILS_H

#include <string>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MACROMol.h>

#include <string_view>

namespace RDKit {
class RWMol;
class Conformer;

enum class MACROTemplateNames {
  AsEntered,      //<! use the name of the temlate as entered in the MACRO Mol
  UseFirstName,   //<!Use the first name in the template
                  // def (For AA, the 3 letter code
  UseSecondName,  //<!use the second name in the tempate def (
                  // For AA, the 1 letter code)
  All             //<! use all names in the template def
};

struct RDKIT_FILEPARSERS_EXPORT MolFromMACROMolParams {
  bool includeLeavingGroups =
      true; /**< when true, leaving groups on atoms that are not exo-bonded
               are retained.  When false, no leaving groups are retained */
  MACROTemplateNames macroTemplateNames = MACROTemplateNames::All;
};

enum class MACROUseTemplateName {
  UseFirstName,   //<!Use the first name in the template
                  // def (For AA, the 3 letter code
  UseSecondName,  //<!use the second name in the tempate def (
                  // For AA, the 1 letter code)
};

struct RDKIT_FILEPARSERS_EXPORT MolToMACROParams {
  MACROUseTemplateName macroUseTemplateName =
      MACROUseTemplateName::UseFirstName;
};

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::RWMol> MolFromMACROMol(
    MACROMol *macroMol,
    const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams,
    const RDKit::MolFromMACROMolParams &molFromMACROMolParams);

RDKIT_FILEPARSERS_EXPORT void MACROMolToSCSRMolFile(
    RDKit::MACROMol &macroMol, const std::string &fName,
    const RDKit::MolWriterParams &params, int confId);

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::MACROMol> MolToMACROMol(
    const ROMol &mol, RDKit::MACROMolTemplateLib &templates,
    MolToMACROParams molToMACROMolParams = MolToMACROParams());

}  // namespace RDKit

#endif
