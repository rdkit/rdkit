//
//  Copyright (C) 20126 Tad Hurst and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "FileParsers.h"
#include "FileParserUtils.h"
#include "MolSGroupParsing.h"
#include <RDGeneral/StreamOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SubstanceGroup.h>

#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/FileParsers/SCSRUtils.h>
#include <GraphMol/MACROMol.h>
#include <GraphMol/FileParsers/MACROMolUtils.h>

namespace RDKit {

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

struct HbondQueryData {
  std::string d_smarts;
  std::string d_name;
  std::vector<HydrogenBondConnection> d_hBondConnections;
  std::shared_ptr<ROMol> dp_mol;
  HbondQueryData(const std::string &smarts, const std::string &name,
                 std::vector<HydrogenBondConnection> hBondConnections)
      : d_smarts(smarts),
        d_name(name),
        d_hBondConnections(hBondConnections),
        dp_mol{SmartsToMol(smarts)} {}
};

void makeHydrogenBonds(unsigned int atom1Idx, unsigned int atom2Idx,
                       std::vector<bool> donorFlags1,
                       std::vector<bool> donorFlags2,
                       std::vector<BondToAdd> &hBondsToAdd) {
  if (donorFlags1.size() == 3 && donorFlags2.size() == 3) {
    // see if the two donor flag lists are
    // complementary

    bool complimentary = true;
    for (auto index = 0; index < 3; ++index) {
      // .second is the bool for donors
      if (donorFlags1[index] == donorFlags2[index]) {
        complimentary = false;  // both are donors or both are
                                // acceptors
        break;
      }
    }

    if (complimentary) {
      for (auto index = 0; index < 3; ++index) {
        hBondsToAdd.emplace_back(atom1Idx, atom2Idx,
                                 "Hb" + std::to_string(index + 1),
                                 "Hb" + std::to_string(index + 1));
      }

      return;
    }

    // not complimentary, but could be a wobble pair - check for G-U`}

    // Check for G-U type pairs in a wobble bond
    // configuration first see if one has the
    // configuration DDA  (like G)

    bool foundDDA = false, foundADA = false;
    if (donorFlags1[0] && donorFlags1[1] && !donorFlags1[2]) {
      foundDDA = true;
    } else if (donorFlags2[0] && donorFlags2[1] && !donorFlags2[2]) {
      foundDDA = true;
      std::swap(donorFlags1, donorFlags2);
      std::swap(atom1Idx, atom2Idx);
    }

    if (foundDDA) {
      // see if the other side has the
      // configuration  ADA (like U)
      if (!donorFlags2[0] && donorFlags2[1] && !donorFlags2[2]) {
        foundADA = true;
      }

      if (foundDDA && foundADA) {
        // use connection 1 and 2 from DDA and 0 and 1 from ADA
        // add the two bonds

        for (auto index1 = 1, index2 = 0; index1 < 3; ++index1, ++index2) {
          hBondsToAdd.emplace_back(atom1Idx, atom2Idx,
                                   "Hb" + std::to_string(index1 + 1),
                                   "Hb" + std::to_string(index2 + 1));
        }

        return;
      }
    }
  }

  // check for one of the ends having only two connections
  if (donorFlags1.size() * donorFlags2.size() == 6) /* one is 2 and one is 3*/ {
    if (donorFlags2.size() == 2) {
      // make sure the first set has two atoms, and the second has 3
      std::swap(donorFlags1, donorFlags2);
      std::swap(atom1Idx, atom2Idx);
    }

    if (donorFlags1[0] && !donorFlags1[1]) {
      if (!donorFlags2[1] && donorFlags2[2]) {
        // first end is DA
        // Second end 2nd and 3rd conns are AD like A (DAD) or C (AAD)

        for (auto index1 = 0, index2 = 1; index1 < 2; ++index1, ++index2) {
          hBondsToAdd.emplace_back(atom1Idx, atom2Idx,
                                   "Hb" + std::to_string(index1 + 1),
                                   "Hb" + std::to_string(index2 + 1));
        }

        return;

      } else if (!donorFlags2[0] && donorFlags2[1]) {
        // like U (ADA)
        for (auto index = 0; index < 2; ++index) {
          hBondsToAdd.emplace_back(atom1Idx, atom2Idx,
                                   "Hb" + std::to_string(index + 1),
                                   "Hb" + std::to_string(index + 1));
        }
        return;
      }
    }
  }

  // if here, we have no wobble bond, so just add
  // the bond between the first atoms on both sides
  hBondsToAdd.emplace_back(atom1Idx, atom2Idx, "Hb1", "Hb1");
  return;
}

void skipSpaces(const char *&linePtr) {
  while (*linePtr == ' ') {
    ++linePtr;  // skip spaces
  }
}

std::string getToken(const char *&linePtr, char delim) {
  skipSpaces(linePtr);
  unsigned int charCount = 0;
  while (linePtr[charCount] && linePtr[charCount] != delim) {
    ++charCount;
  }
  std::string res(linePtr, charCount);
  linePtr += charCount;
  if (*linePtr == delim) {
    ++linePtr;  // skip delim
  }
  return res;
}

std::string getQuotedToken(const char *&linePtr) {
  skipSpaces(linePtr);
  std::string res;
  if (*linePtr != '"') {
    return res;
  }
  ++linePtr;  // skip opening quote
  while (*linePtr && *linePtr != '"') {
    if (*linePtr == '\\') {
      // escaped char
      ++linePtr;
      if (*linePtr == '\0') {
        break;
      }
    }

    res += *linePtr++;
  }
  if (*linePtr != '"') {
    res = "";  // error: no closing quote
  } else {
    ++linePtr;  // skip closing quote
  }

  return res;
}

void parseTemplateLine(std::string lineStr,
                       unsigned int &line,  std::string &templateClass, std::vector<std::string> &templateNames, std::vector<std::pair<std::string,std::string>> &otherTokens) {

  // TEMPLATE 1 AA/Cya/Cya/ NATREPLACE=AA/A COMMENT=comment
  // FULLNAME=fullname CATEGORY=cat UNIQUEID=uniqueid CASNUMBER=xxxx,
  // COLLABORATOR=col, PROTECTION=prot

  // other attributes are allowed.  We capture them and ignore them, except
  // for writing them back out

  if (lineStr.substr(0, 9) != "TEMPLATE ") {
    std::ostringstream errout;
    errout << "Expected \"TEMPLATE\" at line  " << line;
    throw FileParseException(errout.str());
  }

  const char *linePtr = lineStr.c_str() + 9;
  std::string token = getToken(linePtr, ' ');  // Template ID
  if (token.empty()) {
    std::ostringstream errout;
    errout << "Expected a Template ID at line  " << line;
    throw FileParseException(errout.str());
  }
  // get the class and template names

  token = getToken(linePtr, ' ');  // Template ID
  if (token.empty()) {
    std::ostringstream errout;
    errout
        << "Type/Name(s) string of the form \"AA/Gly/G/\" was not found at line  "
        << line;
    throw FileParseException(errout.str());
  }

  // get the class and template names from the token
  std::vector<std::string> subTokens;
  boost::algorithm::split(subTokens, token, boost::algorithm::is_any_of("/"));
  if (subTokens.size() < 3) {
    std::ostringstream errout;
    errout << "Type/Name(s) string is not of the form \"AA/Gly/G/\" at line  "
           << line;
    throw FileParseException(errout.str());
  }

  templateClass = subTokens[0];
  //templateMol->setProp(common_properties::molAtomClass, subTokens[0]);

  templateNames.clear();
  std::copy_if(subTokens.begin() + 1, subTokens.end(),
               std::back_inserter(templateNames),
               [](const std::string &s) { return !s.empty(); });

  //templateMol->setProp(common_properties::templateNames, templateNames);

  // now parse attrs of the form ATTRNAME=VALUE or ATTRNAME="VALUE WITH
  // SPACES"

  otherTokens.clear();
  while (true) {
    std::string attrName = getToken(linePtr, '=');
    if (attrName.empty()) {
      break;
    }
    std::string attrValue;
    if (*linePtr == '"') {
      attrValue = getQuotedToken(linePtr);
    } else {
      attrValue = getToken(linePtr, ' ');
    }

    otherTokens.push_back(std::pair(attrName, attrValue));
    //templateMol->setProp(attrName, attrValue);
  }

  if (*linePtr != '\0') {
    std::ostringstream errout;
    errout
        << "extra characters at the end of a TEMPLATE definition line at line  "
        << line;
    throw FileParseException(errout.str());
  }

  return;
}



std::unique_ptr<RDKit::MACROMol> MACROMolFromSCSRDataStream(
    std::istream &inStream, unsigned int &line,
    const RDKit::v2::FileParsers::MolFileParserParams &params,
    const SCSRBaseHbondOptions scsrBaseHbondOptions) {
  bool chiralityPossible = false;
  if (inStream.eof()) {
    return nullptr;
  }

  auto localParams = params;
  localParams.parsingSCSRMol = true;

  auto tempMol =
      RDKit::v2::FileParsers::MolFromMolDataStream(inStream, line, localParams);
  auto res = std::unique_ptr<RDKit::MACROMol>(new RDKit::MACROMol(tempMol));
  res->updatePropertyCache(false);

  // now get all of the templates

  auto tempStr = RDKit::FileParserUtils::getV3000Line(&inStream, line);

  if (tempStr != "BEGIN TEMPLATE") {
    std::ostringstream errout;
    errout << "BEGIN TEMPLATE not found at line  " << line;
    throw RDKit::FileParseException(errout.str());
  }
  tempStr = RDKit::FileParserUtils::getV3000Line(&inStream, line);

  // TEMPLATE 1 AA/Cya/Cya/ NATREPLACE=AA/A COMMENT=comment
  // FULLNAME=fullname CATEGORY=cat UNIQUEID=uniqueid CASNUMBER=xxxx,
  // COLLABORATOR=col, PROTECTION=prot

  // other attributes are allowed.  We capture them and ignore them, except
  // for writing them back out

  while (tempStr.substr(0, 8) == "TEMPLATE") {
    auto templateMol = std::unique_ptr<RDKit::RWMol>(new RDKit::RWMol());
    std::string templateClass;
    std::vector<std::string> templateNames;
    std::vector<std::pair<std::string,std::string>> otherTokens;

    parseTemplateLine( tempStr.c_str(), line, templateClass, templateNames, otherTokens);

    auto molComplete = false;
    RDKit::Conformer *conf = nullptr;
    try {
      unsigned int nAtoms = 0, nBonds = 0;
      bool expectMEND = false;
      molComplete = RDKit::FileParserUtils::ParseV3000CTAB(
          &inStream, line, templateMol.get(), conf, chiralityPossible, nAtoms,
          nBonds, params.strictParsing, expectMEND);
    } catch (RDKit::MolFileUnhandledFeatureException &e) {
      // unhandled mol file feature, show an error
      res.reset();
      delete conf;
      conf = nullptr;
      BOOST_LOG(rdErrorLog) << " Unhandled CTAB feature: '" << e.what()
                            << "'. Molecule skipped." << std::endl;

      if (!inStream.eof()) {
        tempStr = RDKit::getLine(inStream);
      }
      ++line;
      while (!inStream.eof() && !inStream.fail() &&
             tempStr.substr(0, 6) != "M  END" &&
             tempStr.substr(0, 4) != "$$$$") {
        tempStr = RDKit::getLine(inStream);
        ++line;
      }
      molComplete = !inStream.eof() || tempStr.substr(0, 6) == "M  END" ||
                    tempStr.substr(0, 4) == "$$$$";
    } catch (RDKit::FileParseException &e) {
      // catch our exceptions and throw them back after cleanup
      delete conf;
      conf = nullptr;
      throw e;
    }

    if (!molComplete) {
      delete conf;
      conf = nullptr;
      std::ostringstream errout;
      errout
          << "Problems encountered parsing Mol data, M  END or BEGIN TEMPLATE missing around line "
          << line;
      throw RDKit::FileParseException(errout.str());
    }

    if (templateMol) {
      auto tempParams = params;
      tempParams.sanitize = false;
      tempParams.removeHs = false;
      RDKit::FileParserUtils::finishMolProcessing(
          templateMol.get(), chiralityPossible, tempParams);

      std::unique_ptr<MACROMolTemplate> newTemplate(new MACROMolTemplate(
          templateMol, templateClass, templateNames, otherTokens));
      newTemplate->updatePropertyCache(false);

      res->addTemplate(newTemplate);
      templateMol = nullptr;

    }

    tempStr = RDKit::FileParserUtils::getV3000Line(&inStream, line);
  }

  if (tempStr != "END TEMPLATE") {
    std::ostringstream errout;
    errout << "END TEMPLATE not found at line  " << line;
    throw RDKit::FileParseException(errout.str());
  }

  tempStr = RDKit::getLine(inStream);
  ++line;
  if (tempStr.substr(0, 6) != "M  END") {
    std::ostringstream errout;
    errout << "M  END not found at line  " << line;
    throw RDKit::FileParseException(errout.str());
  }

  // now validate the templates to the main mol
  // each atom in the mol that requires a template should match one,
  // and the ATTCHORD entries of each main atom must be consistent
  // with the template mol it matches.

  unsigned int atomCount = res.get()->getNumAtoms();
  for (unsigned int atomIdx = 0; atomIdx < atomCount; ++atomIdx) {
    auto atom = res->getAtomWithIdx(atomIdx);

    auto templateMol = res->atomIdxToTemplatePtr(atomIdx);
    if (templateMol == nullptr) {
      continue;  // no template for this atom - regular atom
    }

    auto mainSup = templateMol->getMainSgroup();

    std::vector<std::pair<unsigned int, std::string>> attchOrds;
    if (atom->hasProp(RDKit::common_properties::molAttachOrderTemplate)) {
      atom->getProp(RDKit::common_properties::molAttachOrderTemplate,
                    attchOrds);

      // check the ATTCHORD entries
      auto supAttachPoints = mainSup->getAttachPoints();
      for (const auto &attachOrd : attchOrds) {
        if (find_if(supAttachPoints.begin(), supAttachPoints.end(),
                    [&attachOrd](auto a) {
                      return a.id == attachOrd.second;
                    }) == supAttachPoints.end()) {
          std::ostringstream errout;
          errout << "No attachment point found in the template for atom "
                 << atom->getIdx() << " and attachment label "
                 << attachOrd.second;
          throw RDKit::FileParseException(errout.str());
        }
      }
    }
  }

  // now fix the bonds.   SCSR format has the connection points for macro
  // atoms defined in the main mol atoms, while MACROMol has that
  // information in the bonds.
  // H-bonds are even trickier.  SCSR can have a hydrogen
  // bond to a BASE tempate instance, but not have the connection
  // point in the atom def nor in the template definition.
  // It also might
  // have connection info in the macro atom in the main mol, and have one
  // or more SAP entries in the template with the same name.
  //
  // MACROMols require that every connection from or to
  // a template atom have a single uniquely nameed connection point in the
  // template definition (e.g. Hb1, Hb2, Hb3).
  //

  std::vector<std::pair<unsigned int, unsigned int>> bondsToRemove;
  std::vector<BondToAdd> hBondsToAdd;

  // first find the h-bond names for each template that is
  // involved in an h-bond. fix the template to have the standard hbond
  // attachment names (Hb1, Hb2, Hb3)

  std::map<MACROMolTemplate *, SCSRHbondData> baseTemplateHbondData;

  for (auto bond : res->bonds()) {
    if (bond->hasProp(common_properties::_MolFileBondAttachPt1) ||
        bond->hasProp(common_properties::_MolFileBondAttachPt2)) {
      std::ostringstream errout;
      errout
          << "MACROMol style bond connections found in an SCSR file/block for bond "
          << bond->getIdx();
      throw RDKit::FileParseException(errout.str());
    }

    if (bond->getBondType() != Bond::BondType::HYDROGEN) {
      continue;
    }

    auto atom1 = bond->getBeginAtom();
    auto atom2 = bond->getEndAtom();

    for (auto atomToCheck : {atom1, atom2}) {
      auto otherAtomId =
          (atomToCheck == atom1) ? atom2->getIdx() : atom1->getIdx();

      auto templatePtr = res->atomIdxToTemplatePtr(atomToCheck->getIdx());
      if (templatePtr == nullptr) {
        continue;  // not a template atom
      }

      auto atomClass =
          atomToCheck->getProp<std::string>(common_properties::molAtomClass);
      if (atomClass != "BASE") {
        continue;  // only process H-bonds to BASE templates
      }

      std::vector<std::pair<unsigned int, std::string>> attchOrds;
      if (atomToCheck->hasProp(
              RDKit::common_properties::molAttachOrderTemplate)) {
        atomToCheck->getProp(RDKit::common_properties::molAttachOrderTemplate,
                             attchOrds);
      }

      // look for an attachment point on the atom that matches the other
      // atom of the bond.
      std::string hbondAttachName = "";
      std::vector<std::pair<unsigned int, std::string>> newAttachOrds;
      for (auto &attchOrd : attchOrds) {
        if (attchOrd.first == otherAtomId) {
          if (hbondAttachName != "") {
            std::ostringstream errout;
            errout << "Multiple attachment points found for an h-bond to atom "
                   << atomToCheck->getIdx();
            throw RDKit::FileParseException(errout.str());
          }

          hbondAttachName = attchOrd.second;
        } else {
          newAttachOrds.push_back(attchOrd);
        }
      }

      if (hbondAttachName != "") {
        baseTemplateHbondData[templatePtr] =
            SCSRHbondData(hbondAttachName);  // no donor flags yet
        atomToCheck->setProp(RDKit::common_properties::molAttachOrderTemplate,
                             newAttachOrds);  // remove the references to the
                                              // h-bond attach points
      } else {  // it is an h-bond to a template ref, but no atom attach point
                // references it
        baseTemplateHbondData[templatePtr] = SCSRHbondData(
            std::string(""));  // no attach point name (nor donor flags yet)
      }
    }
  }

  for (auto &hbondData : baseTemplateHbondData) {
    auto templateMol = hbondData.first;
    auto &saps = templateMol->getMainSgroup()->getAttachPoints();
    std::vector<RDKit::SubstanceGroup::AttachPoint> newSaps;
    bool sapsChanged = false;
    for (auto &sap : saps) {
      if (sap.id == hbondData.second.oldAttachLabel) {
        if (scsrBaseHbondOptions == SCSRBaseHbondOptions::UseSapOne ||
            scsrBaseHbondOptions == SCSRBaseHbondOptions::UseSapAll) {
          auto isDonor =
              templateMol->getAtomWithIdx(sap.aIdx)->getTotalNumHs() > 0;

          hbondData.second.donorFlags.push_back(isDonor);
          RDKit::SubstanceGroup::AttachPoint newSap(sap);

          newSaps.emplace_back(
              sap.aIdx, -1,
              "Hb" + std::to_string(hbondData.second.donorFlags.size()));
          sapsChanged = true;
        }
      } else {
        newSaps.emplace_back(sap);
      }
    }

    // now for automatic h-bond site detection

    if (scsrBaseHbondOptions == SCSRBaseHbondOptions::Auto) {
      static const std::vector<HbondQueryData> hbondQueries = {
          {"[NH]1C=NC2=C1N=C([NH2])[NH]C2=O",
           "Guanine",
           {HydrogenBondConnection(7, true), HydrogenBondConnection(8, true),
            HydrogenBondConnection(10, false)}},
          {"[NH]1C=NC2=C1N=[CH]N=C2[NH2]",
           "Adenine",
           {HydrogenBondConnection(6, true), HydrogenBondConnection(7, false),
            HydrogenBondConnection(9, true)}},
          {"[NH]1C=CC([NH2])=NC1=O",
           "Cytosine",
           {HydrogenBondConnection(7, false), HydrogenBondConnection(5, false),
            HydrogenBondConnection(4, true)}},
          {"[NH]1C=CC(=O)[NH]C1=O",
           "Thyamine-Uracil",
           {HydrogenBondConnection(7, false), HydrogenBondConnection(5, true),
            HydrogenBondConnection(4, false)}},
          {"[NH]1C=NC2=C1N=C[NH]C2=O",
           "Inosine",
           {HydrogenBondConnection(7, true),
            HydrogenBondConnection(9, false)}}};

      templateMol->updatePropertyCache(false);

      for (const auto &querySmi : hbondQueries) {
        RDKit::SubstructMatchParameters params;
        auto match = SubstructMatch(*templateMol, *querySmi.dp_mol, params);

        if (match.empty()) {
          continue;
        }

        unsigned int hBondCount = 0;
        hbondData.second.donorFlags.clear();
        for (const auto &hBondConnection : querySmi.d_hBondConnections) {
          newSaps.emplace_back(
              match[0][hBondConnection.d_templateAtomIdx].second, -1,
              "Hb" + std::to_string(++hBondCount));
          hbondData.second.donorFlags.push_back(hBondConnection.d_isDonor);
          sapsChanged = true;
        }

        break;  // only take the first match
      }
    }

    if (sapsChanged) {
      templateMol->getMainSgroup()->clearAttachPoints();
      for (auto &newSap : newSaps) {
        templateMol->getMainSgroup()->addAttachPoint(newSap.aIdx, newSap.lvIdx,
                                                     newSap.id);
      }
    }
  }

  //  go back over the bonds and fix them - if either atom is a template atom,
  //  we have to figure out which attachment point to bond

  for (auto bond : res->bonds()) {
    auto atom1 = bond->getBeginAtom();
    auto atom2 = bond->getEndAtom();
    auto atom1Idx = atom1->getIdx();
    auto atom2Idx = atom2->getIdx();

    auto templatePtr1 = res->atomIdxToTemplatePtr(atom1Idx);
    auto templatePtr2 = res->atomIdxToTemplatePtr(atom2Idx);

    if (bond->getBondType() == Bond::BondType::HYDROGEN &&
        templatePtr1 != nullptr && templatePtr2 != nullptr &&
        baseTemplateHbondData.contains(templatePtr1) &&
        baseTemplateHbondData.contains(templatePtr2)) {
      // hydrogen bond between base template atoms

      if (SCSRBaseHbondOptions::Ignore != scsrBaseHbondOptions) {
        auto donorFlags1 = baseTemplateHbondData[templatePtr1].donorFlags;
        auto donorFlags2 = baseTemplateHbondData[templatePtr2].donorFlags;

        makeHydrogenBonds(atom1Idx, atom2Idx, donorFlags1, donorFlags2,
                          hBondsToAdd);
      }

      bondsToRemove.emplace_back(bond->getBeginAtom()->getIdx(),
                                 bond->getEndAtom()->getIdx());
    }

    else {  // regular bonds  or an H-bond where one or both sides don't have
            // template attach point info

      for (auto atomToCheck : {atom1, atom2}) {
        auto otherAtomId =
            (atomToCheck == atom1) ? atom2->getIdx() : atom1->getIdx();

        if (res->atomIdxToTemplatePtr(atomToCheck->getIdx()) == nullptr) {
          continue;  // not a template atom
        }

        std::vector<std::pair<unsigned int, std::string>> attchOrds;
        if (atomToCheck->hasProp(
                RDKit::common_properties::molAttachOrderTemplate)) {
          atomToCheck->getProp(RDKit::common_properties::molAttachOrderTemplate,
                               attchOrds);
        }

        // look for an attachment point on the atom that matches the
        // other atom of the bond.
        std::pair<unsigned int, std::string> *foundAttchOrd = nullptr;
        for (auto &attchOrd : attchOrds) {
          if (attchOrd.first == otherAtomId) {
            if (foundAttchOrd != nullptr) {
              std::ostringstream errout;
              errout << "Multiple attachment points found for a bond to atom "
                     << atomToCheck->getIdx();
              throw RDKit::FileParseException(errout.str());
            }

            foundAttchOrd = &attchOrd;
            break;
          }
        }
        if (foundAttchOrd == nullptr) {
          std::ostringstream errout;
          errout << "Attachment point not found for bond " << bond->getIdx()
                 << " and atom " << atom1->getIdx();
          throw RDKit::FileParseException(errout.str());
        }
        // add the attachment point info to the bond

        auto propToAdd = (atomToCheck == atom1)
                             ? common_properties::_MolFileBondAttachPt1
                             : common_properties::_MolFileBondAttachPt2;
        bond->setProp(propToAdd, (*foundAttchOrd).second);

        // remove the attachment point from the atom
        std::erase(attchOrds, *foundAttchOrd);
        atomToCheck->setProp(RDKit::common_properties::molAttachOrderTemplate,
                             attchOrds);
      }
    }
  }

  res->beginBatchEdit();
  for (auto &hbondToRemove : bondsToRemove) {
    // the following call will delelte all bonds that have the same begin
    // and end atoms.
    res->removeBond(hbondToRemove.first, hbondToRemove.second);
  }
  res->commitBatchEdit();

  res->beginBatchEdit();
  for (auto &hBondToAdd : hBondsToAdd) {
    auto newBondIdx =
        res->addBond(hBondToAdd.d_beginAtomIdx, hBondToAdd.d_endAtomIdx,
                     Bond::BondType::HYDROGEN);
    auto newBond = res->getBondWithIdx(newBondIdx - 1);
    if (!hBondToAdd.d_attachPt1.empty()) {
      newBond->setProp(common_properties::_MolFileBondAttachPt1,
                       hBondToAdd.d_attachPt1);
    }
    if (!hBondToAdd.d_attachPt2.empty()) {
      newBond->setProp(common_properties::_MolFileBondAttachPt2,
                       hBondToAdd.d_attachPt2);
    }
  }
  res->commitBatchEdit();

  return res;
}

//------------------------------------------------
//
//  Read an SCSR molecule from a string
//
//------------------------------------------------
std::unique_ptr<RDKit::MACROMol> MACROMolFromSCSRBlock(
    const std::string &molBlock,
    const RDKit::v2::FileParsers::MolFileParserParams &params,
    const SCSRBaseHbondOptions scsrBaseHbondOptions) {
  std::istringstream inStream(molBlock);
  unsigned int line = 0;
  return MACROMolFromSCSRDataStream(inStream, line, params,
                                    scsrBaseHbondOptions);
}

//------------------------------------------------
//
//  Read an SCSR molecule from a file
//
//------------------------------------------------
std::unique_ptr<RDKit::MACROMol> MACROMolFromSCSRFile(
    const std::string &fName,
    const RDKit::v2::FileParsers::MolFileParserParams &params,
    const SCSRBaseHbondOptions scsrBaseHbondOptions) {
  std::ifstream inStream(fName.c_str());
  if (!inStream || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fName;
    throw RDKit::BadFileException(errout.str());
  }
  if (!inStream.eof()) {
    unsigned int line = 0;
    return MACROMolFromSCSRDataStream(inStream, line, params,
                                      scsrBaseHbondOptions);
  } else {
    return std::unique_ptr<RDKit::MACROMol>();
  }
}

std::unique_ptr<RDKit::RWMol> MolFromSCSRDataStream(
    std::istream &inStream, unsigned int &line,
    const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams,
    const RDKit::MolFromMACROMolParams &molFromMACROMolParams,
    const SCSRBaseHbondOptions scsrBaseHbondOptions) {
  auto macroMol = RDKit::MACROMolFromSCSRDataStream(
      inStream, line, molFileParserParams, scsrBaseHbondOptions);
  return MolFromMACROMol(macroMol.get(), molFileParserParams,
                         molFromMACROMolParams);
}

std::unique_ptr<RDKit::RWMol> MolFromSCSRBlock(
    const std::string &molBlock,
    const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams,
    const RDKit::MolFromMACROMolParams &molFromMACROMolParams,
    const SCSRBaseHbondOptions scsrBaseHbondOptions) {
  auto macroMol = RDKit::MACROMolFromSCSRBlock(molBlock, molFileParserParams,
                                               scsrBaseHbondOptions);
  return MolFromMACROMol(macroMol.get(), molFileParserParams,
                         molFromMACROMolParams);
}

std::unique_ptr<RDKit::RWMol> MolFromSCSRFile(
    const std::string &fName,
    const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams,
    const RDKit::MolFromMACROMolParams &molFromMACROMolParams,
    const SCSRBaseHbondOptions scsrBaseHbondOptions) {
  auto macroMol = RDKit::MACROMolFromSCSRFile(fName, molFileParserParams,
                                              scsrBaseHbondOptions);
  auto res = MolFromMACROMol(macroMol.get(), molFileParserParams,
                             molFromMACROMolParams);
  return res;
}

std::string MACROMolToSCSRMolBlock(MACROMol &macroMol,
                                   const RDKit::MolWriterParams &params,
                                   int confId, bool keepBaseHbondInfo) {
  RDKit::Utils::LocaleSwitcher switcher;
  std::string res = "";
  auto localParams = params;
  localParams.forceV3000 = true;
  localParams.writeEndMolLine = false;

  // make a copy of the main macromol so we can modify it for output without
  // affecting the original

  auto localMol = std::unique_ptr<RWMol>(new RWMol(macroMol));

  // move the connection point info from the bonds to the atoms for the main
  // mol for H-bond connections in templates
  //  1) retain only one connection point in each atom named "Hb"
  //  2) retain only one bond connecting the same two template atoms
  //  3) change the template so that all hbond saps are called Hb

  // for each BASE tempate, the list of
  // h-bond connection point names to fix
  std::map<MACROMolTemplate *, std::vector<std::string>> templateHbondConnections;

  // for each BASE tempate, the list of
  // h-bond connection point name BASES to fix
  // for sets of hbond connections like Hb1, Hb2, Hb3, this would be "Hb"
  // THis allows us to change all such connections even when only one or two
  // are used in the main mol
  std::map<MACROMolTemplate *, std::vector<std::string>> templateHbondConnectionBases;

  std::vector<std::pair<unsigned int, unsigned int>>
      hBonds;  // the bonds to keep (re-add).

  // RDKit will remove all bonds with the same begin and end atoms when called
  // to remove one, so we must keep track of the FIRST hbond between any two
  // atoms We will remove all of the hbonds with those two atoms, then re-add
  // one hbond back.

  for (auto bond : localMol->bonds()) {
    // see if this is an h-bond and either of the atoms are macro refs
    if (bond->getBondType() == Bond::BondType::HYDROGEN &&
        (bond->getBeginAtom()->hasProp(common_properties::molAtomClass) ||
         bond->getEndAtom()->hasProp(common_properties::molAtomClass))) {
      if (std::find_if(hBonds.begin(), hBonds.end(),
                       [bond](std::pair<unsigned int, unsigned int> &hBond) {
                         return (hBond.first == bond->getBeginAtomIdx() &&
                                 hBond.second == bond->getEndAtomIdx()) ||
                                (hBond.second == bond->getBeginAtomIdx() &&
                                 hBond.first == bond->getEndAtomIdx());
                       }) != hBonds.end()) {
        continue;  // next bond
      }

      hBonds.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());

      for (auto atomToCheck : {bond->getBeginAtom(), bond->getEndAtom()}) {
        std::string attrToCheck;
        unsigned int otherAtomId;
        if (atomToCheck == bond->getBeginAtom()) {
          attrToCheck = common_properties::_MolFileBondAttachPt1;
          otherAtomId = bond->getEndAtomIdx();
        } else {
          attrToCheck = common_properties::_MolFileBondAttachPt2;
          otherAtomId = bond->getBeginAtomIdx();
        }

        if (bond->hasProp(attrToCheck)) {
          std::vector<std::pair<unsigned int, std::string>> attachPoints;
          atomToCheck->getProp(common_properties::molAttachOrderTemplate,
                               attachPoints);

          // make sure the atom have one and only one  "Hb" attach point.
          if (std::find_if(attachPoints.begin(), attachPoints.end(),
                           [](const auto &attachPoint) {
                             return attachPoint.second == "Hb";
                           }) == attachPoints.end()) {
            attachPoints.emplace_back(otherAtomId, "Hb");
            atomToCheck->setProp(common_properties::molAttachOrderTemplate,
                                 attachPoints);
          }
        }

        // mark which templates need to be fixed so the attach points are "Hb"

        auto templatePtr = macroMol.atomIdxToTemplatePtr(atomToCheck->getIdx());
        if (!templateHbondConnections.contains(templatePtr)) {
          templateHbondConnections[templatePtr] = std::vector<std::string>();
        }
        if (!templateHbondConnectionBases.contains(templatePtr)) {
          templateHbondConnectionBases[templatePtr] =
              std::vector<std::string>();
        }

        auto attrName = bond->getProp<std::string>(attrToCheck);
        templateHbondConnections[templatePtr].push_back(attrName);
        auto attrNameBase = attrName;
        if (isdigit(attrName[attrName.size() - 1])) {
          // change Hb1, Hb2, etc to Hb
          templateHbondConnectionBases[templatePtr].push_back(
              attrName.substr(0, attrName.size() - 1));
        }

        bond->clearProp(attrToCheck);
      }
    } else {  // here for regular bonds - not  hbonds between templates

      if (bond->hasProp(common_properties::_MolFileBondAttachPt1)) {
        auto attachPt = bond->getProp<std::string>(
            common_properties::_MolFileBondAttachPt1);

        std::vector<std::pair<unsigned int, std::string>> attachPts;
        bond->getBeginAtom()->getPropIfPresent(
            RDKit::common_properties::molAttachOrderTemplate, attachPts);
        attachPts.push_back({bond->getEndAtomIdx(), attachPt});
        bond->getBeginAtom()->setProp(
            RDKit::common_properties::molAttachOrderTemplate, attachPts);
        bond->clearProp(common_properties::_MolFileBondAttachPt1);
      }

      if (bond->hasProp(common_properties::_MolFileBondAttachPt2)) {
        auto attachPt = bond->getProp<std::string>(
            common_properties::_MolFileBondAttachPt2);

        std::vector<std::pair<unsigned int, std::string>> attachPts;
        bond->getEndAtom()->getPropIfPresent(
            RDKit::common_properties::molAttachOrderTemplate, attachPts);
        attachPts.push_back({bond->getBeginAtomIdx(), attachPt});
        bond->getEndAtom()->setProp(
            RDKit::common_properties::molAttachOrderTemplate, attachPts);
        bond->clearProp(common_properties::_MolFileBondAttachPt2);
      }
    }
  }

  // remove the duplicate h-bonds

  if (!hBonds.empty()) {
    localMol->beginBatchEdit();
    for (auto &hbondToRemove : hBonds) {
      // the following call will delete all bonds that have the same begin
      // and end atoms.
      localMol->removeBond(hbondToRemove.first, hbondToRemove.second);
    }
    localMol->commitBatchEdit();

    localMol->beginBatchEdit();
    for (auto &hBond : hBonds) {
      // the following call adds one hbond back for each set that was removed
      localMol->addBond(hBond.first, hBond.second, Bond::BondType::HYDROGEN);
    }
    localMol->commitBatchEdit();
    hBonds.clear();
  }
  res += RDKit::MolToMolBlock(*(localMol), localParams, confId);

  // now write out all of the templates

  res += "M  V30 BEGIN TEMPLATE\n";
  for (unsigned int templateId = 0; templateId < macroMol.getNumTemplates();
       ++templateId) {
    auto macroMolTemplate = macroMol.getTemplate(templateId);
    std::string templateClass;
    std::vector<std::string> templateNames;
    std::string natReplace = "";
    std::string comment = "";
    std::string fullname = "";
    std::string category = "";
    std::string uniqueId = "";
    std::string casNumber = "";
    std::string collaborator = "";
    std::string protection = "";

    macroMolTemplate->getProp(common_properties::molAtomClass, templateClass);

    macroMolTemplate->getProp(common_properties::templateNames, templateNames);

    res += "M  V30 TEMPLATE " + std::to_string(templateId + 1) + " " +
           templateClass;
    for (auto templateName : templateNames) {
      res += "/" + templateName;
    }
    res += "/";

    for (auto &propName : macroMolTemplate->getPropList()) {
      if (propName != common_properties::molAtomClass &&
          propName != common_properties::templateNames && propName[0] != '_') {
        std::string propVal;
        macroMolTemplate->getProp<std::string>(propName, propVal);
        res += " " + propName + "=" + propVal;
      }
    }

    res += "\n";

    // copy template mol so we can modify it for output without affecting the
    // original, if needed

    RWMol *molToUse = macroMolTemplate;
    std::unique_ptr<RWMol> tMol;

    // check to see if the template need to be modified for hbonds
    // connections. if we are keeping them, they are changed to "Hb" if not
    // keeping them, they are removed from the template definition

    if (templateHbondConnections.contains(macroMolTemplate) ||
        templateHbondConnectionBases.contains(macroMolTemplate)) {
      tMol.reset(new RWMol(*macroMolTemplate));
      molToUse = tMol.get();

      auto &hbondConnections = templateHbondConnections[macroMolTemplate];
      auto &hbondConnectionBases = templateHbondConnectionBases[macroMolTemplate];
      auto &sgroups = RDKit::getSubstanceGroups(*(tMol.get()));
      auto &localMainSgroup =
          sgroups[macroMolTemplate->getMainSgroup()->getProp<unsigned int>(
                      "index") -
                  1];
      auto &attachPoints = localMainSgroup.getAttachPoints();

      if (keepBaseHbondInfo) {
        for (auto &attachPoint : attachPoints) {
          if (std::find(hbondConnections.begin(), hbondConnections.end(),
                        attachPoint.id) != hbondConnections.end()) {
            attachPoint.id = "Hb";
          } else if (isdigit(attachPoint.id[attachPoint.id.size() - 1])) {
            auto baseId = attachPoint.id.substr(0, attachPoint.id.size() - 1);
            if (std::find(hbondConnectionBases.begin(),
                          hbondConnectionBases.end(),
                          baseId) != hbondConnectionBases.end()) {
              attachPoint.id = "Hb";
            }
          }
        }
      } else  // remove hbond connection points for this template
      {
        std::erase_if(attachPoints, [hbondConnections, hbondConnectionBases](
                                        const SubstanceGroup::AttachPoint &ap) {
          if (std::find(hbondConnections.begin(), hbondConnections.end(),
                        ap.id) != hbondConnections.end()) {
            return true;
          }
          if (isdigit(ap.id[ap.id.size() - 1])) {
            auto baseId = ap.id.substr(0, ap.id.size() - 1);
            if (std::find(hbondConnectionBases.begin(),
                          hbondConnectionBases.end(),
                          baseId) != hbondConnectionBases.end()) {
              return true;
            }
          }
          return false;
        });
      }
    }

    boost::dynamic_bitset<> aromaticBonds(macroMolTemplate->getNumBonds());
    prepareMol(*molToUse, localParams, aromaticBonds);
    res += FileParserUtils::getV3000CTAB(*molToUse, aromaticBonds, confId,
                                         localParams.precision);
  }

  res += "M  V30 END TEMPLATE\n";
  res += "M  END\n";

  return res;
}

void MACROMolToSCSRMolFile(RDKit::MACROMol &macroMol, const std::string &fName,
                           const RDKit::MolWriterParams &params, int confId) {
  auto *outStream = new std::ofstream(fName.c_str());
  if (!(*outStream) || outStream->bad()) {
    delete outStream;
    std::ostringstream errout;
    errout << "Bad output file " << fName;
    throw BadFileException(errout.str());
  }
  std::string outString = MACROMolToSCSRMolBlock(macroMol, params, confId);
  *outStream << outString;
  delete outStream;
}

}  // namespace RDKit
