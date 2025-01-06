//
//  Copyright (C) 2002-2024 Tad Hurst and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "FileParserUtils.h"
#include "MolSGroupParsing.h"
#include <RDGeneral/StreamOps.h>

#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>

using namespace RDKit::SGroupParsing;

namespace RDKit {

namespace v2 {
namespace FileParsers {

//------------------------------------------------
//
//  Read a SCVSR molecule from a stream
//
//------------------------------------------------
std::unique_ptr<RDKit::SCSRMol> ScsrFromScsrDataStream(
    std::istream &inStream, unsigned int &line,
    const RDKit::v2::FileParsers::MolFileParserParams &params) {
  std::string tempStr;
  bool molComplete = false;
  bool chiralityPossible = false;
  Utils::LocaleSwitcher ls;
  // mol name
  line++;
  tempStr = getLine(inStream);
  if (inStream.eof()) {
    return nullptr;
  }
  auto res = std::unique_ptr<RDKit::SCSRMol>(new RDKit::SCSRMol());

  res->getMol()->setProp(common_properties::_Name, tempStr);

  // info
  line++;
  tempStr = getLine(inStream);
  res->getMol()->setProp("_MolFileInfo", tempStr);
  if (tempStr.length() >= 22) {
    std::string dimLabel = tempStr.substr(20, 2);
    // Unless labelled as 3D we assume 2D
    if (dimLabel == "3d" || dimLabel == "3D") {
      res->getMol()->setProp(common_properties::_3DConf, 1);
    }
  }
  // comments
  line++;
  tempStr = getLine(inStream);
  res->getMol()->setProp("_MolFileComments", tempStr);

  unsigned int nAtoms = 0, nBonds = 0, nLists = 0, chiralFlag = 0, nsText = 0,
               nRxnComponents = 0;
  int nReactants = 0, nProducts = 0, nIntermediates = 0;
  (void)nLists;  // read from the file but unused
  (void)nsText;
  (void)nRxnComponents;
  (void)nReactants;
  (void)nProducts;
  (void)nIntermediates;
  // counts line, this is where we really get started
  line++;
  tempStr = getLine(inStream);

  if (tempStr.size() < 6) {
    if (res) {
      res = nullptr;
    }
    std::ostringstream errout;
    errout << "Counts line too short: '" << tempStr << "' on line" << line;
    throw FileParseException(errout.str());
  }

  unsigned int spos = 0;
  // this needs to go into a try block because if the lexical_cast throws an
  // exception we want to catch throw a different exception
  try {
    nAtoms = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    spos = 3;
    nBonds = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    spos = 6;
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << tempStr.substr(spos, 3)
           << "' to unsigned int on line " << line;
    throw FileParseException(errout.str());
  }
  try {
    spos = 6;
    if (tempStr.size() >= 9) {
      nLists = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    }

    spos = 12;
    if (tempStr.size() >= spos + 3) {
      chiralFlag = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    }

    spos = 15;
    if (tempStr.size() >= spos + 3) {
      nsText = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    }

    spos = 18;
    if (tempStr.size() >= spos + 3) {
      nRxnComponents =
          FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    }

    spos = 21;
    if (tempStr.size() >= spos + 3) {
      nReactants = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    }

    spos = 24;
    if (tempStr.size() >= spos + 3) {
      nProducts = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    }

    spos = 27;
    if (tempStr.size() >= spos + 3) {
      nIntermediates =
          FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    }

  } catch (boost::bad_lexical_cast &) {
    // some SD files (such as some from NCI) lack all the extra information
    // on the header line, so ignore problems parsing there.
  }

  if (tempStr.size() > 35) {
    if (tempStr.size() < 39 || tempStr[34] != 'V') {
      std::ostringstream errout;
      errout << "CTAB version string invalid at line " << line;
      if (params.strictParsing) {
        throw FileParseException(errout.str());
      } else {
        BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
      }
    } else if (tempStr.substr(34, 5) != "V3000") {
      std::ostringstream errout;
      errout << "SCSR Mol files in not V3000 at line" << line;
      throw FileParseException(errout.str());
    }
  }

  res->getMol()->setProp(common_properties::_MolFileChiralFlag, chiralFlag);

  Conformer *conf = nullptr;
  try {
    if (nAtoms != 0 || nBonds != 0) {
      std::ostringstream errout;
      errout << "V3000 mol blocks should have 0s in the initial counts line. "
                "(line: "
             << line << ")";
      if (params.strictParsing) {
        throw FileParseException(errout.str());
      } else {
        BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
      }
    }
    bool expectMEND = false;
    bool expectMacroAtoms = true;
    molComplete = FileParserUtils::ParseV3000CTAB(
        &inStream, line, (RWMol *)res->getMol(), conf, chiralityPossible,
        nAtoms, nBonds, params.strictParsing, expectMEND, expectMacroAtoms);
  } catch (MolFileUnhandledFeatureException &e) {
    // unhandled mol file feature, show an error
    res.reset();
    delete conf;
    conf = nullptr;
    BOOST_LOG(rdErrorLog) << " Unhandled CTAB feature: '" << e.what()
                          << "'. Molecule skipped." << std::endl;

    if (!inStream.eof()) {
      tempStr = getLine(inStream);
    }
    ++line;
    while (!inStream.eof() && !inStream.fail() &&
           tempStr.substr(0, 6) != "M  END" && tempStr.substr(0, 4) != "$$$$") {
      tempStr = getLine(inStream);
      ++line;
    }
    molComplete = !inStream.eof() || tempStr.substr(0, 6) == "M  END" ||
                  tempStr.substr(0, 4) == "$$$$";
  } catch (FileParseException &e) {
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
    throw FileParseException(errout.str());
  }

  if (res) {
    FileParserUtils::finishMolProcessing((RWMol *)res->getMol(),
                                         chiralityPossible, params);
  }

  // now get all of the templates

  tempStr = FileParserUtils::getV3000Line(&inStream, line);

  if (tempStr != "BEGIN TEMPLATE") {
    delete conf;
    conf = nullptr;
    std::ostringstream errout;
    errout << "BEGIN TEMPLATE not found at line  " << line;
    throw FileParseException(errout.str());
  }
  tempStr = FileParserUtils::getV3000Line(&inStream, line);

  // TEMPLATE 1 AA/Cya/Cya/ NATREPLACE=AA/A
  while (tempStr.substr(0, 8) == "TEMPLATE") {
    std::vector<std::string> tokens;
    boost::algorithm::split(tokens, tempStr, boost::algorithm::is_space());

    std::string natReplace = "";
    if (tokens.size() == 4) {
      if (tokens[3].size() < 12 || tokens[3].substr(0, 11) != "NATREPLACE=") {
        delete conf;
        conf = nullptr;
        std::ostringstream errout;
        errout << "Bad NATREPLACE entry at line  " << line;
        throw FileParseException(errout.str());
      }
      natReplace = tokens[3].substr(12);
    } else if (tokens.size() != 3) {
      delete conf;
      conf = nullptr;
      std::ostringstream errout;
      errout << "Bad TEMPLATE at line  " << line;
      throw FileParseException(errout.str());
    }
    boost::algorithm::split(tokens, tokens[2],
                            boost::algorithm::is_any_of("/"));
    if (tokens.size() < 3) {
      delete conf;
      conf = nullptr;
      std::ostringstream errout;
      errout << "Type/Name(s) string is not of the form \"AA/Gly/G/\" at line  "
             << line;
      throw FileParseException(errout.str());
    }

    res->addTemplate(new ROMol());
    auto templateMol = (RWMol *)res->getTemplate(res->getTemplateCount() - 1);

    templateMol->setProp(common_properties::natReplace, natReplace);
    templateMol->setProp(common_properties::molAtomClass, tokens[0]);

    std::vector<std::string> templateNames;
    for (unsigned int i = 1; i < tokens.size(); ++i) {
      if (tokens[i] != "") {
        templateNames.push_back(tokens[i]);
      }
    }

    templateMol->setProp(common_properties::templateNames, templateNames);

    Conformer *conf = nullptr;
    try {
      unsigned int nAtoms = 0, nBonds = 0;
      bool expectMEND = false;
      molComplete = FileParserUtils::ParseV3000CTAB(
          &inStream, line, templateMol, conf, chiralityPossible, nAtoms, nBonds,
          params.strictParsing, expectMEND);
    } catch (MolFileUnhandledFeatureException &e) {
      // unhandled mol file feature, show an error
      res.reset();
      delete conf;
      conf = nullptr;
      BOOST_LOG(rdErrorLog) << " Unhandled CTAB feature: '" << e.what()
                            << "'. Molecule skipped." << std::endl;

      if (!inStream.eof()) {
        tempStr = getLine(inStream);
      }
      ++line;
      while (!inStream.eof() && !inStream.fail() &&
             tempStr.substr(0, 6) != "M  END" &&
             tempStr.substr(0, 4) != "$$$$") {
        tempStr = getLine(inStream);
        ++line;
      }
      molComplete = !inStream.eof() || tempStr.substr(0, 6) == "M  END" ||
                    tempStr.substr(0, 4) == "$$$$";
    } catch (FileParseException &e) {
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
      throw FileParseException(errout.str());
    }

    if (templateMol) {
      FileParserUtils::finishMolProcessing(templateMol, chiralityPossible,
                                           params);
    }

    tempStr = FileParserUtils::getV3000Line(&inStream, line);
  }

  if (tempStr != "END TEMPLATE") {
    delete conf;
    conf = nullptr;
    std::ostringstream errout;
    errout << "END TEMPLATE not found at line  " << line;
    throw FileParseException(errout.str());
  }

  tempStr = getLine(inStream);
  ++line;
  if (tempStr.substr(0, 6) != "M  END") {
    delete conf;
    conf = nullptr;
    std::ostringstream errout;
    errout << "M  END not found at line  " << line;
    throw FileParseException(errout.str());
  }

  // now validate the templates to the main mol
  // each atom in the mol that requires a template should match one,
  // and the ATTCHORD entries of each main atom must be consistent
  // with the template mol it matches.

  unsigned int atomCount = res.get()->getMol()->getNumAtoms();
  for (unsigned int atomIdx = 0; atomIdx < atomCount; ++atomIdx) {
    auto atom = res->getMol()->getAtomWithIdx(atomIdx);

    std::string dummyLabel = "";
    std::string atomClass = "";
    std::vector<RDKit::AtomAttchOrd> attchOrds;
    if (atom->hasProp(common_properties::dummyLabel)) {
      dummyLabel = atom->getProp<std::string>(common_properties::dummyLabel);
    }
    if (atom->hasProp(common_properties::molAtomClass)) {
      atomClass = atom->getProp<std::string>(common_properties::molAtomClass);
    }
    if (atom->hasProp(common_properties::molAttachOrderTemplate)) {
      attchOrds = atom->getProp<std::vector<RDKit::AtomAttchOrd>>(
          common_properties::molAttachOrderTemplate);
    }
    if (dummyLabel != "" && atomClass != "" && attchOrds.size() > 0) {
      // find the template that matches the class and label

      bool templateFound = false;
      for (unsigned int templateIdx = 0;
           !templateFound && templateIdx < res.get()->getTemplateCount();
           ++templateIdx) {
        auto templateMol = res->getTemplate(templateIdx);
        if (templateMol->getProp<std::string>(
                common_properties::molAtomClass) == atomClass) {
          for (auto templateName :
               templateMol->getProp<std::vector<std::string>>(
                   common_properties::templateNames)) {
            if (templateName == dummyLabel) {
              // now check the ATTCHORD entries

              for (auto attchOrd : attchOrds) {
                for (auto sgroup : RDKit::getSubstanceGroups(*templateMol)) {
                  if (sgroup.getProp<std::string>("TYPE") == "SUP" &&
                      sgroup.getProp<std::string>("CLASS") == atomClass) {
                    for (auto sap : sgroup.getAttachPoints()) {
                      if (sap.id == attchOrd.getLabel()) {
                        templateFound = true;
                        break;
                      }
                    }
                  }
                  if (templateFound) {
                    break;
                  }
                }
                if (templateFound) {
                  break;
                }
              }
            }
            if (templateFound) {
              break;
            }
          }
        }
      }

      if (!templateFound) {
        std::ostringstream errout;
        errout << "No template found for atom " << atom->getIdx() << " class "
               << atomClass << " label " << dummyLabel;
        throw FileParseException(errout.str());
      }
    }
  }

  return res;
}

//------------------------------------------------
//
//  Read a molecule from a string
//
//------------------------------------------------
std::unique_ptr<RDKit::SCSRMol> ScsrFromScsrBlock(
    const std::string &molBlock,
    const RDKit::v2::FileParsers::MolFileParserParams &params) {
  std::istringstream inStream(molBlock);
  unsigned int line = 0;
  return ScsrFromScsrDataStream(inStream, line, params);
}

//------------------------------------------------
//
//  Read a molecule from a file
//
//------------------------------------------------
std::unique_ptr<RDKit::SCSRMol> ScsrFromScsrFile(
    const std::string &fName, const MolFileParserParams &params) {
  std::ifstream inStream(fName.c_str());
  if (!inStream || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fName;
    throw BadFileException(errout.str());
  }
  if (!inStream.eof()) {
    unsigned int line = 0;
    return ScsrFromScsrDataStream(inStream, line, params);
  } else {
    return std::unique_ptr<RDKit::SCSRMol>();
  }
}

class OriginAtomDef {
 public:
  OriginAtomDef(unsigned int mainAtomIdx, unsigned int templateAtomIdx)
      : mainAtomIdx(mainAtomIdx), templateAtomIdx(templateAtomIdx) {}

 private:
  unsigned int mainAtomIdx;
  unsigned int templateAtomIdx;

 public:
  bool operator<(const OriginAtomDef &other) const {
    if (mainAtomIdx != other.mainAtomIdx) {
      return mainAtomIdx < other.mainAtomIdx;
    }
    return templateAtomIdx < other.templateAtomIdx;
  }
};

class OriginAtomConnection {
 public:
  OriginAtomConnection(unsigned int mainAtomIdx, std::string attachLabel)
      : mainAtomIdx(mainAtomIdx), attachLabel(attachLabel) {}

 private:
  unsigned int mainAtomIdx;
  std::string attachLabel;

 public:
  unsigned int getMainAtomIdx() const { return mainAtomIdx; }
  std::string getAttachLabel() const { return attachLabel; }

  bool operator<(const OriginAtomConnection &other) const {
    if (mainAtomIdx != other.getMainAtomIdx()) {
      return mainAtomIdx < other.mainAtomIdx;
    }
    return attachLabel < other.attachLabel;
  }
};

unsigned int getNewAtomForBond(
    const Atom *atom, unsigned int otherAtomIdx,
    const std::map<OriginAtomDef, unsigned int> &originAtomMap,
    const std::map<OriginAtomConnection, unsigned int> &attachMap) {
  unsigned int beginSeq = 0;
  unsigned int atomIdx = atom->getIdx();
  if (!atom->getPropIfPresent<unsigned int>(common_properties::molAtomSeqId,
                                            beginSeq)) {
    return originAtomMap.at(OriginAtomDef(atomIdx, UINT_MAX));
  }

  // if here , it is a template atom

  auto attchOrds = atom->getProp<std::vector<RDKit::AtomAttchOrd>>(
      common_properties::molAttachOrderTemplate);
  for (auto attchOrd : attchOrds) {
    if (attchOrd.getAtomIdx() == otherAtomIdx) {
      auto attachMapIt =
          attachMap.find(OriginAtomConnection(atomIdx, attchOrd.getLabel()));
      if (attachMapIt == attachMap.end()) {
        return UINT_MAX;
      }
      return originAtomMap.at(OriginAtomDef(atomIdx, attachMapIt->second));
      // return originAtomMap.at(OriginAtomDef(
      //     atomIdx,
      //     attachMap.at(OriginAtomConnection(atomIdx, attchOrd.getLabel()))));
      break;
    }
  }

  // error attachment ord not found
  return UINT_MAX;
}

std::unique_ptr<RDKit::RWMol> MolFromScsr(
    const RDKit::SCSRMol &scsrMol, const MolFromScsrParams &molFromScsrParams) {
  auto resMol = std::unique_ptr<RWMol>(new RWMol());
  auto mol = scsrMol.getMol();
  // first get some information from the templates to be used when creating the
  // coords for the new atoms. this is a dirty approach that simply expands the
  // orginal macro atom coords to be big enough to hold any expanded macro atom.
  // No attempt is made to make this look nice, or to avoid overlaps.
  std::vector<RDGeom::Point3D> templateCentroids;
  double maxSize = 0.0;

  bool skipCoords = true;
  const Conformer *conf = nullptr;
  std::unique_ptr<Conformer> newConf(nullptr);
  if (mol->getNumConformers() != 0) {
    skipCoords = false;
    conf = &mol->getConformer(0);
    newConf.reset(new Conformer(scsrMol.getMol()->getNumAtoms()));
    newConf->set3D(conf->is3D());
  }

  for (unsigned int templateIdx = 0;
       !skipCoords && templateIdx < scsrMol.getTemplateCount(); ++templateIdx) {
    auto templateMol = scsrMol.getTemplate(templateIdx);
    RDGeom::Point3D sumOfCoords;
    auto confCount = templateMol->getNumConformers();
    if (confCount == 0) {
      skipCoords = true;
      break;
    }
    auto templateConf = templateMol->getConformer(0);
    RDGeom::Point3D maxCoord = templateConf.getAtomPos(0);
    RDGeom::Point3D minCoord = maxCoord;
    for (unsigned int atomIdx = 0; atomIdx < templateMol->getNumAtoms();
         ++atomIdx) {
      auto atomCoord = templateConf.getAtomPos(atomIdx);
      sumOfCoords += atomCoord;

      if (atomCoord.x > maxCoord.x) {
        maxCoord.x = atomCoord.x;
      }
      if (atomCoord.y > maxCoord.y) {
        maxCoord.y = atomCoord.y;
      }
      if (atomCoord.z > maxCoord.z) {
        maxCoord.z = atomCoord.z;
      }
      if (atomCoord.x < minCoord.x) {
        minCoord.x = atomCoord.x;
      }
      if (atomCoord.y < minCoord.y) {
        minCoord.y = atomCoord.y;
      }
      if (atomCoord.z < minCoord.z) {
        minCoord.z = atomCoord.z;
      }
    }
    templateCentroids.push_back(sumOfCoords / templateMol->getNumAtoms());
    if (maxCoord.x - minCoord.x > maxSize) {
      maxSize = maxCoord.x - minCoord.x;
    }
    if (maxCoord.y - minCoord.y > maxSize) {
      maxSize = maxCoord.y - minCoord.y;
    }
    if (maxCoord.z - minCoord.z > maxSize) {
      maxSize = maxCoord.z - minCoord.z;
    }
  }

  // for each atom in the main mol, expand it to full atom form

  std::vector<StereoGroup> newStereoGroups;
  std::vector<Atom *> absoluteAtoms;
  std::vector<Bond *> absoluteBonds;

  std::map<OriginAtomDef, unsigned int> originAtomMap;
  // maps main atom# and attach label to template atom#
  std::map<OriginAtomConnection, unsigned int> attachMap;
  std::vector<std::vector<unsigned int>> atomMap;

  unsigned int atomCount = mol->getNumAtoms();
  for (unsigned int atomIdx = 0; atomIdx < atomCount; ++atomIdx) {
    auto atom = mol->getAtomWithIdx(atomIdx);
    atomMap.push_back(std::vector<unsigned int>());

    std::string dummyLabel = "";
    std::string atomClass = "";

    if (!atom->getPropIfPresent(common_properties::dummyLabel, dummyLabel) ||
        !atom->getPropIfPresent(common_properties::molAtomClass, atomClass) ||
        dummyLabel == "" || atomClass == "") {
      // NOT a template atom - just copy it
      auto newAtom = new Atom(*atom);
      resMol->addAtom(newAtom, true, true);
      atomMap[atomIdx].push_back(newAtom->getIdx());

      // templatesFound.push_back(nullptr);
      originAtomMap[OriginAtomDef(atomIdx, UINT_MAX)] = newAtom->getIdx();
      if (!skipCoords) {
        newConf->setAtomPos(newAtom->getIdx(),
                            conf->getAtomPos(atomIdx) * maxSize);
      }
    } else {  // it is a macro atom - expand it

      unsigned int seqId = 0;
      atom->getPropIfPresent(common_properties::molAtomSeqId, seqId);

      //  find the template that matches the class and label

      const ROMol *templateMol = nullptr;
      unsigned int templateIdx;
      bool templateFound = false;
      for (templateIdx = 0; templateIdx < scsrMol.getTemplateCount();
           ++templateIdx) {
        templateMol = scsrMol.getTemplate(templateIdx);
        std::vector<std::string> templateNames;
        std::string templateName;
        std::string templateAtomClass;
        if (templateMol->getPropIfPresent<std::string>(
                common_properties::molAtomClass, templateAtomClass) &&
            templateAtomClass == atomClass &&
            templateMol->getPropIfPresent<std::vector<std::string>>(
                common_properties::templateNames, templateNames) &&
            std::find(templateNames.begin(), templateNames.end(), dummyLabel) !=
                templateNames.end()) {
          templateFound = true;
          break;
        }
      }

      if (!templateFound) {
        std::ostringstream errout;
        errout << "No template found for atom " << atom->getIdx() << " class "
               << atomClass << " label " << dummyLabel;
        throw FileParseException(errout.str());
      }

      // find all the atoms in the leaving groups for this atom's
      // expansion
      // also record the attachment atoms in the template

      auto attchOrds = atom->getProp<std::vector<RDKit::AtomAttchOrd>>(
          common_properties::molAttachOrderTemplate);

      std::vector<unsigned int> leavingGroupAtoms;

      // first find the sgroup that is the base for this atom's template

      const SubstanceGroup *sgroup = nullptr;
      for (auto &sgroupToTest : RDKit::getSubstanceGroups(*templateMol)) {
        std::string sup;
        std::string sgroupAtomClass;
        if (sgroupToTest.getPropIfPresent<std::string>("TYPE", sup) &&
            sup == "SUP" &&
            sgroupToTest.getPropIfPresent<std::string>("CLASS",
                                                       sgroupAtomClass) &&
            sgroupAtomClass == atomClass) {
          sgroup = &sgroupToTest;
          break;
        }
      }

      if (sgroup == nullptr) {
        throw FileParseException("No SUP sgroup found for atom");
      }

      // find  the attachment points used by this atom and record them in
      // attachMap.
      // also, if we are including leaving groups that are not
      // substitued out, find the leaving group atoms that will not be included
      // in the final molecule

      for (auto attchOrd : attchOrds) {
        for (auto attachPoint : sgroup->getAttachPoints()) {
          if (attachPoint.id == attchOrd.getLabel()) {
            attachMap[OriginAtomConnection(atomIdx, attachPoint.id)] =
                attachPoint.aIdx;
            if (molFromScsrParams.includeLeavingGroups) {
              // only those leaving groups that are substitued are NOT included
              for (auto lGrp : RDKit::getSubstanceGroups(*templateMol)) {
                auto atomsInSAP = lGrp.getAtoms();
                std::string lgrpClass;
                if (lGrp.getProp<std::string>("TYPE") == "SUP" &&
                    lGrp.getPropIfPresent<std::string>("CLASS", lgrpClass) &&
                    lgrpClass == "LGRP" &&
                    std::find(atomsInSAP.begin(), atomsInSAP.end(),
                              attachPoint.lvIdx) != atomsInSAP.end()) {
                  leavingGroupAtoms.insert(leavingGroupAtoms.end(),
                                           atomsInSAP.begin(),
                                           atomsInSAP.end());
                  break;
                }
              }
            }
            break;
          }
        }
      }

      // if we are not including atoms from any leaving groups even if they are
      // not substitued, find all the leaving group atoms

      if (!molFromScsrParams.includeLeavingGroups) {
        // all leaving groups are not to be included
        for (auto lGrp : RDKit::getSubstanceGroups(*templateMol)) {
          auto atomsInSAP = lGrp.getAtoms();
          std::string lgrpClass;
          if (lGrp.getProp<std::string>("TYPE") == "SUP" &&
              lGrp.getPropIfPresent<std::string>("CLASS", lgrpClass) &&
              lgrpClass == "LGRP") {
            leavingGroupAtoms.insert(leavingGroupAtoms.end(),
                                     atomsInSAP.begin(), atomsInSAP.end());
          }
        }
      }

      // add an superatom sgroup to mark the atoms from this macro atom
      // expansion
      std::string typ = "SUP";
      std::unique_ptr<SubstanceGroup> newSgroup(
          new SubstanceGroup(resMol.get(), typ));
      newSgroup->setProp<unsigned int>("index", seqId);
      newSgroup->setProp("LABEL", dummyLabel + "_" + std::to_string(seqId));

      // copy the atoms of the template into the new molecule
      // leaving group atoms are not copied if they are part of the main
      // mol's atom's attachord
      if (!skipCoords) {
        newConf->resize(newConf->getNumAtoms() + templateMol->getNumAtoms());
      }

      RDGeom::Point3D newCentroid = conf->getAtomPos(atomIdx) * maxSize;

      for (auto templateAtom : templateMol->atoms()) {
        if (std::find(leavingGroupAtoms.begin(), leavingGroupAtoms.end(),
                      templateAtom->getIdx()) == leavingGroupAtoms.end()) {
          auto newAtom = new Atom(*templateAtom);

          resMol->addAtom(newAtom, true, true);
          newSgroup->addAtomWithIdx(newAtom->getIdx());

          originAtomMap[OriginAtomDef(atomIdx, templateAtom->getIdx())] =
              newAtom->getIdx();

          atomMap[atomIdx].push_back(newAtom->getIdx());
          if (!skipCoords) {
            newConf->setAtomPos(newAtom->getIdx(),
                                newCentroid +
                                    templateMol->getConformer().getAtomPos(
                                        templateAtom->getIdx()) -
                                    templateCentroids[templateIdx]);
          }
        }
      }

      // copy the bonds of the template into the new molecule
      // if the bonds are between atoms NOT in the leaving group

      for (auto bond : templateMol->bonds()) {
        if (std::find(leavingGroupAtoms.begin(), leavingGroupAtoms.end(),
                      bond->getBeginAtomIdx()) == leavingGroupAtoms.end() &&
            std::find(leavingGroupAtoms.begin(), leavingGroupAtoms.end(),
                      bond->getEndAtomIdx()) == leavingGroupAtoms.end()) {
          auto newBeginAtomIdx =
              originAtomMap[OriginAtomDef(atomIdx, bond->getBeginAtomIdx())];
          auto newEndAtomIdx =
              originAtomMap[OriginAtomDef(atomIdx, bond->getEndAtomIdx())];

          auto newBond = new Bond(bond->getBondType());
          newBond->setBeginAtomIdx(newBeginAtomIdx);
          newBond->setEndAtomIdx(newEndAtomIdx);
          newBond->updateProps(*bond, false);
          resMol->addBond(newBond, true);
          newSgroup->addBondWithIdx(newBond->getIdx());
        }
      }

      if (newSgroup->getIsValid()) {
        addSubstanceGroup(*resMol, *newSgroup.get());
      }

      // take care of stereo groups in the template

      for (auto &sg : templateMol->getStereoGroups()) {
        std::vector<Atom *> newGroupAtoms;
        std::vector<Bond *> newGroupBonds;

        for (auto sgAtom : sg.getAtoms()) {
          auto originAtom = OriginAtomDef(atomIdx, sgAtom->getIdx());

          auto newAtomPtr = originAtomMap.find(originAtom);
          if (newAtomPtr != originAtomMap.end()) {
            newGroupAtoms.push_back(resMol->getAtomWithIdx(newAtomPtr->second));
          }
        }

        for (auto sgBond : sg.getBonds()) {
          auto originBeginAtom =
              OriginAtomDef(atomIdx, sgBond->getBeginAtomIdx());
          auto originEndAtom = OriginAtomDef(atomIdx, sgBond->getEndAtomIdx());

          auto newBeginAtomPtr = originAtomMap.find(originBeginAtom);
          auto newEndAtomPtr = originAtomMap.find(originEndAtom);
          if (newBeginAtomPtr != originAtomMap.end() &&
              newEndAtomPtr != originAtomMap.end()) {
            auto newBond = resMol->getBondBetweenAtoms(newBeginAtomPtr->second,
                                                       newEndAtomPtr->second);
            if (newBond != nullptr) {
              newGroupBonds.push_back(newBond);
            }
          }
        }

        if (sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
          absoluteAtoms.insert(absoluteAtoms.end(), newGroupAtoms.begin(),
                               newGroupAtoms.end());
          absoluteBonds.insert(absoluteBonds.end(), newGroupBonds.begin(),
                               newGroupBonds.end());
        } else {
          // make a new group

          newStereoGroups.emplace_back(sg.getGroupType(), newGroupAtoms,
                                       newGroupBonds);
        }
      }
    }
  }

  if (resMol->getNumAtoms() == 0) {
    return resMol;
  }

  newConf->resize(resMol->getNumAtoms());
  resMol->addConformer(newConf.release());

  // now deal with the bonds from the original mol.

  for (auto bond : scsrMol.getMol()->bonds()) {
    auto newBeginAtomIdx = getNewAtomForBond(
        bond->getBeginAtom(), bond->getEndAtomIdx(), originAtomMap, attachMap);
    if (newBeginAtomIdx == UINT_MAX) {
      continue;
    }
    auto newEndAtomIdx = getNewAtomForBond(
        bond->getEndAtom(), bond->getBeginAtomIdx(), originAtomMap, attachMap);
    if (newEndAtomIdx == UINT_MAX) {
      continue;
    }

    auto newBond = new Bond(bond->getBondType());
    newBond->setBeginAtomIdx(newBeginAtomIdx);
    newBond->setEndAtomIdx(newEndAtomIdx);
    newBond->updateProps(*bond, false);

    resMol->addBond(newBond, true);
  }

  // copy any attrs from the main mol

  for (auto &prop : mol->getPropList(false, false)) {
    std::string propVal;
    if (mol->getPropIfPresent(prop, propVal)) {
      resMol->setProp(prop, propVal);
    }
  }

  // copy the sgroups from the main mol

  for (auto &sg : getSubstanceGroups(*mol)) {
    if (sg.getIsValid()) {
      RDKit::SubstanceGroup newSg(sg);
      auto &atoms = sg.getAtoms();
      std::vector<unsigned int> newAtoms;
      for (auto atom : atoms) {
        auto originAtom = OriginAtomDef(atom, UINT_MAX);
        auto newAtomPtr = originAtomMap.find(originAtom);
        if (newAtomPtr != originAtomMap.end()) {
          newAtoms.push_back(newAtomPtr->second);
        }
      }

      newSg.setAtoms(newAtoms);
      addSubstanceGroup(*resMol, newSg);
    }
  }

  // take care of stereo groups in the main mol - for atoms that are NOT
  // template refs

  for (auto &sg : mol->getStereoGroups()) {
    std::vector<Atom *> newGroupAtoms;
    std::vector<Bond *> newGroupBonds;

    for (auto sgAtom : sg.getAtoms()) {
      auto originAtom = OriginAtomDef(sgAtom->getIdx(), UINT_MAX);
      auto newAtomPtr = originAtomMap.find(originAtom);
      if (newAtomPtr != originAtomMap.end()) {
        newGroupAtoms.push_back(resMol->getAtomWithIdx(newAtomPtr->second));
      }
    }

    for (auto sgBond : sg.getBonds()) {
      auto originBeginAtom = OriginAtomDef(sgBond->getBeginAtomIdx(), UINT_MAX);
      auto originEndAtom = OriginAtomDef(sgBond->getEndAtomIdx(), UINT_MAX);
      auto newBeginAtomPtr = originAtomMap.find(originBeginAtom);
      auto newEndAtomPtr = originAtomMap.find(originEndAtom);

      if (newBeginAtomPtr != originAtomMap.end() &&
          newEndAtomPtr != originAtomMap.end()) {
        auto newBond = resMol->getBondBetweenAtoms(newBeginAtomPtr->second,
                                                   newEndAtomPtr->second);
        if (newBond != nullptr) {
          newGroupBonds.push_back(newBond);
        }
      }

      if (sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
        absoluteAtoms.insert(absoluteAtoms.end(), newGroupAtoms.begin(),
                             newGroupAtoms.end());
        absoluteBonds.insert(absoluteBonds.end(), newGroupBonds.begin(),
                             newGroupBonds.end());
      } else {
        // make a new group

        newStereoGroups.emplace_back(sg.getGroupType(), newGroupAtoms,
                                     newGroupBonds);
      }
    }
  }

  // make an absolute group that contains any absolute atoms or bonds from
  // either the main mol or the template instantiations

  if (!absoluteAtoms.empty() || !absoluteBonds.empty()) {
    newStereoGroups.emplace_back(StereoGroupType::STEREO_ABSOLUTE,
                                 absoluteAtoms, absoluteBonds);
  }

  if (newStereoGroups.size() > 0) {
    resMol->setStereoGroups(newStereoGroups);
  }

  return resMol;
}

}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
