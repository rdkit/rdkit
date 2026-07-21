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
#include <GraphMol/MacroMol.h>
#include <GraphMol/FileParsers/MacroMolUtils.h>
#include "SCSRUtils.h"
#include <string_view>

namespace RDKit {
class RWMol;
class Conformer;
class SCSRUtils { // this class is a utility container - only static methods
public:
    // Delete the constructor to prevent instantiation
    SCSRUtils() = delete;  //


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


private: 
    class SCSRHbondData {
    public:
    std::string oldAttachLabel;
    std::vector<bool> donorFlags;
    SCSRHbondData() : oldAttachLabel("") {}
    SCSRHbondData(const std::string &oldAttachLabel)
        : oldAttachLabel(oldAttachLabel) {}

    SCSRHbondData(const SCSRHbondData &oldOne)
        : oldAttachLabel(oldOne.oldAttachLabel), donorFlags(oldOne.donorFlags) {}

    SCSRHbondData &operator=(const SCSRHbondData &other) {
        if (this == &other) {
        return *this;
        }
        oldAttachLabel = other.oldAttachLabel;
        donorFlags = other.donorFlags;
        return *this;
    }
    };

    struct HydrogenBondConnection {
    unsigned int d_templateAtomIdx;
    bool d_isDonor;

    HydrogenBondConnection(unsigned int templateAtomIdx, bool isDonor)
        : d_templateAtomIdx(templateAtomIdx), d_isDonor(isDonor) {}
    };

    struct BondToAdd {
    unsigned int d_beginAtomIdx;
    unsigned int d_endAtomIdx;
    std::string d_attachPt1;
    std::string d_attachPt2;
    BondToAdd(unsigned int beginAtomIdx, unsigned int endAtomIdx,
                const std::string &attachPt1, const std::string &attachPt2)
        : d_beginAtomIdx(beginAtomIdx),
            d_endAtomIdx(endAtomIdx),
            d_attachPt1(attachPt1),
            d_attachPt2(attachPt2) {}
    auto getBeginAtomIdx() const { return d_beginAtomIdx; }
    auto getEndAtomIdx() const { return d_endAtomIdx; }
    auto getAttachPt1() const { return d_attachPt1; }
    auto getAttachPt2() const { return d_attachPt2; }
    };


private:
    static void skipSpaces(const char *&linePtr) ;

    static std::string getToken(const char *&linePtr, char delim);


    static void makeHydrogenBonds(unsigned int atom1Idx, unsigned int atom2Idx,
                        std::vector<bool> donorFlags1,
                        std::vector<bool> donorFlags2,
                        std::vector<BondToAdd> &hBondsToAdd);

    static std::string getQuotedToken(const char *&linePtr);

    static void parseTemplateLine(std::string lineStr,
                       unsigned int &line,  std::string &templateClass, std::vector<std::string> &templateNames, std::vector<std::pair<std::string,std::string>> &otherTokens);


public:

    static RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::RWMol> MolFromSCSRDataStream(
        std::istream &inStream, unsigned int &line,
        const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams =
            RDKit::v2::FileParsers::MolFileParserParams(),
        const RDKit::MolFromMacroMolParams &molFromMacroMolParams =
            RDKit::MolFromMacroMolParams(),
        const SCSRBaseHbondOptions scsrBaseHbondOptions =
            SCSRBaseHbondOptions::Auto);

    static RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::RWMol> MolFromSCSRBlock(
        const std::string &molBlock,
        const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams =
            RDKit::v2::FileParsers::MolFileParserParams(),
        const RDKit::MolFromMacroMolParams &molFromMacroMolParams =
            RDKit::MolFromMacroMolParams(),
        const SCSRBaseHbondOptions scsrBaseHbondOptions =
            SCSRBaseHbondOptions::Auto);

    static RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::RWMol> MolFromSCSRFile(
        const std::string &fName,
        const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams =
            RDKit::v2::FileParsers::MolFileParserParams(),
        const RDKit::MolFromMacroMolParams &molFromMacroMolParams =
            RDKit::MolFromMacroMolParams(),
        const SCSRBaseHbondOptions scsrBaseHbondOptions =
            SCSRBaseHbondOptions::Auto);

    static RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::MacroMol>
    MacroMolFromSCSRDataStream(
        std::istream &inStream, unsigned int &line,
        const RDKit::v2::FileParsers::MolFileParserParams &params =
            RDKit::v2::FileParsers::MolFileParserParams(),
        const SCSRBaseHbondOptions scsrBaseHbondOptions =
            SCSRBaseHbondOptions::Auto);

    static RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::MacroMol> MacroMolFromSCSRBlock(
        const std::string &molBlock,
        const RDKit::v2::FileParsers::MolFileParserParams &params =
            RDKit::v2::FileParsers::MolFileParserParams(),
        const SCSRBaseHbondOptions scsrBaseHbondOptions =
            SCSRBaseHbondOptions::Auto);

    static RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::MacroMol> MacroMolFromSCSRFile(
        const std::string &fName,
        const RDKit::v2::FileParsers::MolFileParserParams &params =
            RDKit::v2::FileParsers::MolFileParserParams(),
        const SCSRBaseHbondOptions scsrBaseHbondOptions =
            SCSRBaseHbondOptions::Auto);

    static RDKIT_FILEPARSERS_EXPORT std::string MacroMolToSCSRMolBlock(
        MacroMol &macroMol,
        const RDKit::MolWriterParams &params = RDKit::MolWriterParams(),
        int confId = -1, bool keepBaseHbondInfo = true);

    static RDKIT_FILEPARSERS_EXPORT void MacroMolToSCSRMolFile(
        RDKit::MacroMol &macroMol, const std::string &fName,
        const RDKit::MolWriterParams &params = RDKit::MolWriterParams(),
        int confId = -1);

}; // end of class SCSRUtils

}  // namespace RDKit

#endif
