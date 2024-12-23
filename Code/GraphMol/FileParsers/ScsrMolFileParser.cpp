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
std::unique_ptr<SCSRMol> ScsrMolFromScsrDataStream(
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
  auto res = std::unique_ptr<SCSRMol>(new SCSRMol());

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
    molComplete = FileParserUtils::ParseV3000CTAB(
        &inStream, line, (RWMol *)res->getMol(), conf, chiralityPossible,
        nAtoms, nBonds, params.strictParsing, expectMEND);
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
std::unique_ptr<SCSRMol> ScsrMolFromScsrBlock(
    const std::string &molBlock,
    const RDKit::v2::FileParsers::MolFileParserParams &params) {
  std::istringstream inStream(molBlock);
  unsigned int line = 0;
  return ScsrMolFromScsrDataStream(inStream, line, params);
}

//------------------------------------------------
//
//  Read a molecule from a file
//
//------------------------------------------------
std::unique_ptr<SCSRMol> ScsrMolFromScsrFile(
    const std::string &fName,
    const RDKit::v2::FileParsers::MolFileParserParams &params) {
  std::ifstream inStream(fName.c_str());
  if (!inStream || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fName;
    throw BadFileException(errout.str());
  }
  if (!inStream.eof()) {
    unsigned int line = 0;
    return ScsrMolFromScsrDataStream(inStream, line, params);
  } else {
    return std::unique_ptr<SCSRMol>();
  }
}
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
