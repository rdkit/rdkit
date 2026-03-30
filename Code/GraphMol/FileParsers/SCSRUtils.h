//
//  Copyright (C) 2010-2025 Tad Hurst and other RDKit contributor
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_SCSRUTILS_H
#define RD_SCSRUTILS_H

#include <string>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include "FileParsers.h"
#include <GraphMol/MACROMol.h>
#include <GraphMol/FileParsers/MACROMolUtils.h>
#include "SCSRUtils.h"
#include <string_view>

namespace RDKit {
class RWMol;
class Conformer;

enum class SCSRBaseHbondOptions {
  Ignore,     //<! Do not include base Hbonds in expanded output
  UseSapAll,  //<!use all hbonds defined in SAPs
              // can be more than one per base
  UseSapOne,  //<!use only one SAP hbond per base
              // If multiple SAPs are defined, use the first
              // even if it is not the best
              //(this just maintains the relationship between
              // the to base pairs)
  Auto        //<!For bases that are C,G,A,T,U,In (and
              // derivatives) use the standard Watson-Crick
              // Hbonding.  No SAPs need to be defined, and if
              // defined, they are ignored.
};

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::RWMol> MolFromSCSRDataStream(
    std::istream &inStream, unsigned int &line,
    const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams =
        RDKit::v2::FileParsers::MolFileParserParams(),
    const RDKit::MolFromMACROMolParams &molFromMACROMolParams =
        RDKit::MolFromMACROMolParams(),
    const SCSRBaseHbondOptions scsrBaseHbondOptions =
        SCSRBaseHbondOptions::Auto);

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::RWMol> MolFromSCSRBlock(
    const std::string &molBlock,
    const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams =
        RDKit::v2::FileParsers::MolFileParserParams(),
    const RDKit::MolFromMACROMolParams &molFromMACROMolParams =
        RDKit::MolFromMACROMolParams(),
    const SCSRBaseHbondOptions scsrBaseHbondOptions =
        SCSRBaseHbondOptions::Auto);

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::RWMol> MolFromSCSRFile(
    const std::string &fName,
    const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams =
        RDKit::v2::FileParsers::MolFileParserParams(),
    const RDKit::MolFromMACROMolParams &molFromMACROMolParams =
        RDKit::MolFromMACROMolParams(),
    const SCSRBaseHbondOptions scsrBaseHbondOptions =
        SCSRBaseHbondOptions::Auto);

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::MACROMol>
MACROMolFromSCSRDataStream(
    std::istream &inStream, unsigned int &line,
    const RDKit::v2::FileParsers::MolFileParserParams &params =
        RDKit::v2::FileParsers::MolFileParserParams(),
    const SCSRBaseHbondOptions scsrBaseHbondOptions =
        SCSRBaseHbondOptions::Auto);

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::MACROMol> MACROMolFromSCSRBlock(
    const std::string &molBlock,
    const RDKit::v2::FileParsers::MolFileParserParams &params =
        RDKit::v2::FileParsers::MolFileParserParams(),
    const SCSRBaseHbondOptions scsrBaseHbondOptions =
        SCSRBaseHbondOptions::Auto);

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::MACROMol> MACROMolFromSCSRFile(
    const std::string &fName,
    const RDKit::v2::FileParsers::MolFileParserParams &params =
        RDKit::v2::FileParsers::MolFileParserParams(),
    const SCSRBaseHbondOptions scsrBaseHbondOptions =
        SCSRBaseHbondOptions::Auto);

RDKIT_FILEPARSERS_EXPORT std::string MACROMolToSCSRMolBlock(
    MACROMol &macroMol,
    const RDKit::MolWriterParams &params = RDKit::MolWriterParams(),
    int confId = -1, bool keepBaseHbondInfo = true);

RDKIT_FILEPARSERS_EXPORT void MACROMolToSCSRMolFile(
    RDKit::MACROMol &macroMol, const std::string &fName,
    const RDKit::MolWriterParams &params = RDKit::MolWriterParams(),
    int confId = -1);
}  // namespace RDKit

#endif
