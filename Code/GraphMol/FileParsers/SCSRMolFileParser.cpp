//
//  Copyright (C) 2024 Tad Hurst and other RDKit contributors
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
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>

using namespace RDKit::SGroupParsing;

namespace RDKit {

namespace v2 {
namespace FileParsers {

class SCSRMol {
 private:
  std::unique_ptr<ROMol> p_mol;
  std::vector<std::unique_ptr<ROMol>> p_templates;

 public:
  SCSRMol() {};
  SCSRMol(const SCSRMol &other) = delete;
  SCSRMol(SCSRMol &&other) noexcept = delete;
  SCSRMol &operator=(SCSRMol &&other) noexcept = delete;

  SCSRMol &operator=(const SCSRMol &) = delete;  // disable assignment
  ~SCSRMol() {}

  void addTemplate(std::unique_ptr<ROMol> templateMol) {
    PRECONDITION(templateMol, "bad template molecule");
    p_templates.push_back(std::move(templateMol));
  }

  unsigned int getTemplateCount() const { return p_templates.size(); }

  ROMol *getTemplate(unsigned int index) { return p_templates[index].get(); };

  const ROMol *getMol() const { return p_mol.get(); }

  ROMol *getMol() { return p_mol.get(); }

  void setMol(std::unique_ptr<ROMol> mol) {
    PRECONDITION(mol, "bad molecule");
    p_mol = std::move(mol);
  }
};

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
  *linePtr++;  // skip opening quote
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
    linePtr++;  // skip closing quote
  }

  return res;
}

void parseTemplateLine(RWMol *templateMol, std::string lineStr,
                       unsigned int &line) {
  PRECONDITION(templateMol, "bad template molecule");

  // TEMPLATE 1 AA/Cya/Cya/ NATREPLACE=AA/A COMMENT=comment FULLNAME=fullname
  // CATEGORY=cat UNIQUEID=uniqueid CASNUMBER=xxxx, COLLABORATOR=col,
  // PROTECTION=prot

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

  std::vector<std::string> subTokens;

  // get the class and template names from the token

  boost::algorithm::split(subTokens, token, boost::algorithm::is_any_of("/"));
  if (subTokens.size() < 3) {
    std::ostringstream errout;
    errout << "Type/Name(s) string is not of the form \"AA/Gly/G/\" at line  "
           << line;
    throw FileParseException(errout.str());
  }

  templateMol->setProp(common_properties::molAtomClass, subTokens[0]);

  std::vector<std::string> templateNames;

  for (unsigned int i = 1; i < subTokens.size(); ++i) {
    if (subTokens[i] != "") {
      templateNames.push_back(subTokens[i]);
    }
  }
  templateMol->setProp(common_properties::templateNames, templateNames);

  // now parse attrs of the form ATTRNAME=VALUE or ATTRNAME="VALUE WITH
  // SPACES"

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

    templateMol->setProp(attrName, attrValue);
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

//------------------------------------------------
//
//  Read a SCVSR molecule from a stream
//
//------------------------------------------------
static std::unique_ptr<SCSRMol> SCSRMolFromSCSRDataStream(
    std::istream &inStream, unsigned int &line,
    const RDKit::v2::FileParsers::MolFileParserParams &params) {
  bool chiralityPossible = false;
  if (inStream.eof()) {
    return nullptr;
  }
  auto res = std::unique_ptr<SCSRMol>(new SCSRMol());
  auto localParams = params;
  localParams.parsingSCSRMol = true;
  res->setMol(RDKit::v2::FileParsers::MolFromMolDataStream(inStream, line,
                                                           localParams));

  // now get all of the templates

  auto tempStr = FileParserUtils::getV3000Line(&inStream, line);

  if (tempStr != "BEGIN TEMPLATE") {
    std::ostringstream errout;
    errout << "BEGIN TEMPLATE not found at line  " << line;
    throw FileParseException(errout.str());
  }
  tempStr = FileParserUtils::getV3000Line(&inStream, line);

  // TEMPLATE 1 AA/Cya/Cya/ NATREPLACE=AA/A COMMENT=comment FULLNAME=fullname
  // CATEGORY=cat UNIQUEID=uniqueid CASNUMBER=xxxx, COLLABORATOR=col,
  // PROTECTION=prot

  // other attributes are allowed.  We capture them and ignore them, except
  // for writing them back out

  while (tempStr.substr(0, 8) == "TEMPLATE") {
    res->addTemplate(std::unique_ptr<ROMol>(new ROMol()));
    auto templateMol = (RWMol *)res->getTemplate(res->getTemplateCount() - 1);

    parseTemplateLine(templateMol, tempStr.c_str(), line);

    // res->addTemplate(std::move(templateMol));

    // while (tempStr.substr(0, 8) == "TEMPLATE") {
    //   std::vector<std::string> tokens;
    //   std::vector<std::string> subTokens;
    //   boost::algorithm::split(tokens, tempStr, boost::algorithm::is_space());

    //   if (tokens.size() < 3) {
    //     std::ostringstream errout;
    //     errout << "Bad Template entry at line  " << line;
    //     throw FileParseException(errout.str());
    //   }

    //   // get the class and template names

    //   boost::algorithm::split(subTokens, tokens[2],
    //                           boost::algorithm::is_any_of("/"));
    //   if (subTokens.size() < 3) {
    //     std::ostringstream errout;
    //     errout << "Type/Name(s) string is not of the form \"AA/Gly/G/\" at
    //     line
    //     "
    //            << line;
    //     throw FileParseException(errout.str());
    //   }

    //   res->addTemplate(std::unique_ptr<ROMol>(new ROMol()));
    //   auto templateMol = (RWMol *)res->getTemplate(res->getTemplateCount() -
    //   1);

    //   templateMol->setProp(common_properties::molAtomClass, subTokens[0]);

    //   std::vector<std::string> templateNames;

    //   for (unsigned int i = 1; i < subTokens.size(); ++i) {
    //     if (subTokens[i] != "") {
    //       templateNames.push_back(subTokens[i]);
    //     }
    //   }
    //   templateMol->setProp(common_properties::templateNames, templateNames);

    //   for (unsigned int i = 3; i < tokens.size(); ++i) {
    //     boost::algorithm::split(subTokens, tokens[i],
    //                             boost::algorithm::is_any_of("="));
    //     if (subTokens.size() != 2) {
    //       std::ostringstream errout;
    //       errout
    //           << "Attribute  string is not of the form \"AttrName=value\" at
    //           line  "
    //           << line;
    //       throw FileParseException(errout.str());
    //     }
    //     templateMol->setProp(subTokens[0], subTokens[1]);
    //   }

    auto molComplete = false;
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

    // if (templateMol) {
    auto tempParams = params;
    tempParams.sanitize = false;
    tempParams.removeHs = false;
    FileParserUtils::finishMolProcessing(templateMol, chiralityPossible,
                                         tempParams);
    //}

    tempStr = FileParserUtils::getV3000Line(&inStream, line);
  }

  if (tempStr != "END TEMPLATE") {
    std::ostringstream errout;
    errout << "END TEMPLATE not found at line  " << line;
    throw FileParseException(errout.str());
  }

  tempStr = getLine(inStream);
  ++line;
  if (tempStr.substr(0, 6) != "M  END") {
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
    std::vector<std::pair<unsigned int, std::string>> attchOrds;
    if (atom->hasProp(common_properties::dummyLabel)) {
      dummyLabel = atom->getProp<std::string>(common_properties::dummyLabel);
    }
    if (atom->hasProp(common_properties::molAtomClass)) {
      atomClass = atom->getProp<std::string>(common_properties::molAtomClass);
    }
    if (atom->hasProp(common_properties::molAttachOrderTemplate)) {
      atom->getProp(common_properties::molAttachOrderTemplate, attchOrds);
    }
    if (dummyLabel != "" && atomClass != "") {
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
            if (templateFound) {
              break;
            }
            if (templateName != dummyLabel) {
              continue;  // check next template name
            }
            // find the SUP group for the main connections (not leaving
            // groups)

            const RDKit::SubstanceGroup *mainSUP = nullptr;
            for (auto &sgroup : RDKit::getSubstanceGroups(*templateMol)) {
              if (sgroup.getProp<std::string>("TYPE") == "SUP" &&
                  sgroup.getProp<std::string>("CLASS") == atomClass) {
                mainSUP = &sgroup;
                break;
              }
            }

            if (mainSUP == nullptr) {
              break;  // breaks out of the loop over template names -
                      // continues to next template
            }

            // now check the ATTCHORD entries
            templateFound = true;  // tentative - until some attachment point
                                   // in the atom is not found in the template
            for (const auto &attachOrd : attchOrds) {
              auto supAttachPoints = mainSUP->getAttachPoints();
              if (std::find_if(supAttachPoints.begin(), supAttachPoints.end(),
                               [&attachOrd](auto a) {
                                 return a.id == attachOrd.second;
                               }) == supAttachPoints.end()) {
                templateFound = false;
                break;
              }
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
static std::unique_ptr<SCSRMol> SCSRMolFromSCSRBlock(
    const std::string &molBlock,
    const RDKit::v2::FileParsers::MolFileParserParams &params) {
  std::istringstream inStream(molBlock);
  unsigned int line = 0;
  return SCSRMolFromSCSRDataStream(inStream, line, params);
}

//------------------------------------------------
//
//  Read a molecule from a file
//
//------------------------------------------------
static std::unique_ptr<SCSRMol> SCSRMolFromSCSRFile(
    const std::string &fName, const MolFileParserParams &params) {
  std::ifstream inStream(fName.c_str());
  if (!inStream || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fName;
    throw BadFileException(errout.str());
  }
  if (!inStream.eof()) {
    unsigned int line = 0;
    return SCSRMolFromSCSRDataStream(inStream, line, params);
  } else {
    return std::unique_ptr<SCSRMol>();
  }
}

class MolFromSCSRMolConverter {
 private:
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

  struct HydrogenBondConnection {
    unsigned int d_templateAtomIdx;
    bool df_isDonor;

    HydrogenBondConnection(unsigned int templateAtomIdx, bool isDonor)
        : d_templateAtomIdx(templateAtomIdx), df_isDonor(isDonor) {}
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

  SCSRMol *scsrMol;
  std::unique_ptr<RWMol> resMol;
  const ROMol *mol;
  const MolFileParserParams molFileParserParams;
  const MolFromSCSRParams molFromSCSRParams;

  std::map<unsigned int, std::vector<HydrogenBondConnection>>
      hBondSitesForTemplate;
  std::map<unsigned int, unsigned int> atomIdxToTemplateIdx;

  ROMol *atomIdxToTemplateMol(unsigned int atomIdx) {
    return scsrMol->getTemplate(atomIdxToTemplateIdx[atomIdx]);
  }

  // maps main atom# and template atom# to new atom#
  std::map<OriginAtomDef, unsigned int> originAtomMap;

  // maps main atom# and attach label to template atom#
  std::map<OriginAtomConnection, std::vector<unsigned int>> attachMap;

  unsigned int getNewAtomForBond(const Atom *atom, unsigned int otherAtomIdx) {
    std::string atomClass = "";
    unsigned int atomIdx = atom->getIdx();
    if (!atom->getPropIfPresent<std::string>(common_properties::molAtomClass,
                                             atomClass)) {
      return originAtomMap.at(OriginAtomDef(atomIdx, UINT_MAX));
    }

    // if here , it is a template atom
    // this routine is NOT called for H-bonds, so there is only one atom to
    // attach
    std::vector<std::pair<unsigned int, std::string>> attchOrds;
    atom->getProp(common_properties::molAttachOrderTemplate, attchOrds);
    for (const auto &[idx, lbl] : attchOrds) {
      if (idx == otherAtomIdx) {
        auto attachMapIt = attachMap.find(OriginAtomConnection(atomIdx, lbl));
        if (attachMapIt == attachMap.end() || attachMapIt->second.empty()) {
          throw FileParseException("Attachment ord not found");
        }
        return originAtomMap.at(OriginAtomDef(atomIdx, attachMapIt->second[0]));
      }
    }

    // error attachment ord not found
    throw FileParseException("Attachment ord not found");
  }

  void getNewAtomsForHydrogenBond(
      Atom *atom, unsigned int otherAtomIdx,
      std::vector<HydrogenBondConnection> &hydrogenBondConnections) {
    std::string atomClass = "";
    unsigned int atomIdx = atom->getIdx();
    auto templateMol = atomIdxToTemplateMol(atomIdx);
    hydrogenBondConnections.clear();
    if (!atom->getPropIfPresent<std::string>(common_properties::molAtomClass,
                                             atomClass)) {
      return;  // hatoms only allowed in templates
    }

    // if here , it is a template atom
    // there could be more than one atom in the template that matches the atom
    // for instance, the hydrogen bonds to the template can result in multiple
    // hydrogen bonds in the expanded molecule

    std::vector<std::pair<unsigned int, std::string>> attchOrds;
    atom->getProp(common_properties::molAttachOrderTemplate, attchOrds);
    for (const auto &[idx, lbl] : attchOrds) {
      if (idx == otherAtomIdx) {
        auto attachMapIt = attachMap.find(OriginAtomConnection(atomIdx, lbl));
        if (attachMapIt == attachMap.end()) {
          return;  // error, attachment ord not found
        }

        for (auto templateAtomIdx : attachMapIt->second) {
          auto templateAtom = templateMol->getAtomWithIdx(templateAtomIdx);
          templateAtom->updatePropertyCache();
          auto isDonor = templateAtom->getTotalNumHs() > 0;
          hydrogenBondConnections.emplace_back(templateAtomIdx, isDonor);
          if (molFromSCSRParams.scsrBaseHbondOptions ==
              SCSRBaseHbondOptions::UseSapOne) {
            return;
          }
        }
        return;
      }
    }

    // error attachment ord not found
    return;
  }

  void copySgroupIntoResult(
      const unsigned int atomIdx, const RDKit::SubstanceGroup &sgroup,
      std::string sgroupName,
      std::vector<std::unique_ptr<SubstanceGroup>> &newSgroups,
      RDKit::Conformer *newConf, const RDGeom::Point3D &coordOffset) {
    // add a superatom sgroup to mark the atoms from this macro atom
    // expansion. These new superatom sgroups are not put into the output mol
    // yet, because the bonds have not be added to the mol nor to the sgroup.
    // They are saved in an array to be added later

    const std::string typ = "SUP";
    newSgroups.emplace_back(new SubstanceGroup((ROMol *)resMol.get(), typ));
    auto newSgroup = newSgroups.back().get();
    newSgroup->setProp("LABEL", sgroupName);

    // copy the atoms of the sgroup into the new molecule

    if (newConf) {
      newConf->resize(newConf->getNumAtoms() +
                      atomIdxToTemplateMol(atomIdx)->getNumAtoms());
    }

    for (auto templateAtomIdx : sgroup.getAtoms()) {
      auto templateAtom =
          atomIdxToTemplateMol(atomIdx)->getAtomWithIdx(templateAtomIdx);
      auto newAtom = new Atom(*templateAtom);

      resMol->addAtom(newAtom, true, true);
      newSgroup->addAtomWithIdx(newAtom->getIdx());

      originAtomMap[OriginAtomDef(atomIdx, templateAtomIdx)] =
          newAtom->getIdx();

      // atomMap[atomIdx].push_back(newAtom->getIdx());
      if (newConf) {
        newConf->setAtomPos(
            newAtom->getIdx(),
            coordOffset +
                atomIdxToTemplateMol(atomIdx)->getConformer().getAtomPos(
                    templateAtomIdx));
      }
    }
  }

  void addBond(const Bond::BondType bondType, unsigned int beginAtom,
               unsigned int endAtom) {
    auto newBond = new Bond(bondType);
    newBond->setBeginAtomIdx(beginAtom);
    newBond->setEndAtomIdx(endAtom);
    // newBond->updateProps(*bond, false);
    resMol->addBond(newBond, true);
  }

  void addBonds(const Bond::BondType bondType,
                const unsigned int mainBeginAtomIdx,
                const unsigned int mainEndAtomIdx,
                const std::vector<HydrogenBondConnection> &beginAtoms,
                const std::vector<HydrogenBondConnection> &endAtoms) {
    for (unsigned int i = 0; i < beginAtoms.size(); i++) {
      addBond(bondType,
              originAtomMap[OriginAtomDef(mainBeginAtomIdx,
                                          beginAtoms[i].d_templateAtomIdx)],
              originAtomMap[OriginAtomDef(mainEndAtomIdx,
                                          endAtoms[i].d_templateAtomIdx)]);
    }
  }

  void processBondInMainMol(const Bond *bond) {
    unsigned int newBeginAtom;
    unsigned int newEndAtom;

    if (bond->getBondType() != Bond::HYDROGEN) {
      newBeginAtom =
          getNewAtomForBond(bond->getBeginAtom(), bond->getEndAtomIdx());
      if (newBeginAtom == UINT_MAX) {
        throw FileParseException("Error getting new atom for bond");
      }

      newEndAtom =
          getNewAtomForBond(bond->getEndAtom(), bond->getBeginAtomIdx());
      if (newEndAtom == UINT_MAX) {
        throw FileParseException("Error getting new atom for bond");
      }

      addBond(bond->getBondType(), newBeginAtom, newEndAtom);
      return;
    }

    // if here it is a hydrogen bond

    // the processing of H-bonds is contorlled by the SCSRBaseHbondOptions
    // member of the SCSRMolFileParserParams parameter.  The options are:

    //  Ignore - if this is selected, all H-bonds are
    //  ignored and not processed.

    // UseSapAll = 1,  // if this is selected, all SAPs
    // for the hbond are used.  They must be defined in the template in the
    // correct order, which starts with the first atom nearest the Al
    // connection, and continues sequentially

    // If there are 3 sites on eash side and they are complimentary (the
    // donors match acceptors and vice versa), we add the bonds. and
    // we are done.

    // If not, we check to check to see if they comply to a wobble bond
    // configuration. There are four generally accepted wobble bonds, and we
    // deal with these four types only.  The four known wobble bonds are:
    // 1.  I-C
    // 2.  I-U
    // 3.  I-A
    // 4.  G-U
    // "I" stands for inosine - it has only two available hbond sites . The
    // pair G-U has three sites on each end, but they are not complimentary.
    // (https://en.wikipedia.org/wiki/Wobble_base_pair#:~:text=A%20wobble%20base%20pair%20is,hypoxanthine%2Dcytosine%20(I%2DC).

    // For pairs that have I or something like it, the configuration must be
    // DA, so the two atoms form the two H bonds. The other side (C,U or A)
    // has confiuration of AAD (C), ADA (U), or DAD (A). For C-type bases and
    // A type bases, the second and third atoms are used (AD), and for U
    // types, the first two atoms are used (AD).

    // for the GU pair, both sides have 3 atoms, but they are not
    // complimentary.  The second and third sites on the G side are used (DA),
    // and the first two sites on the U side are used (AD).

    // in any other case, we punt and add just one bond,  between the first
    // atom on both sides even if they are NOT complemenary.  This just
    // indicates that the sides are h-bonded somehow and keeps the overall
    // pairing straight.

    // UseSapOne
    // if this is selected, use only one SAP hbond per base
    // If multiple SAPs are defined, use the first
    // even if it is not the best
    //(this just maintains the relationship between
    // the to base pairs)

    // Auto
    // For bases that are C,G,A,T,U,In (and
    // derivatives) use the standard Watson-Crick
    // Hbonding.  No SAPs need to be defined, and if
    // defined, they are ignored.
    // the definitions of the binding sites are determined by substructure
    // matching

    std::vector<HydrogenBondConnection> beginHatomConnections;
    std::vector<HydrogenBondConnection> endHatomConnections;

    if (molFromSCSRParams.scsrBaseHbondOptions ==
        SCSRBaseHbondOptions::Ignore) {
      return;
    }
    if (molFromSCSRParams.scsrBaseHbondOptions ==
            SCSRBaseHbondOptions::UseSapAll ||
        molFromSCSRParams.scsrBaseHbondOptions ==
            SCSRBaseHbondOptions::UseSapOne) {
      // get the hydrogen bond sites from the template SAPs
      getNewAtomsForHydrogenBond(bond->getBeginAtom(), bond->getEndAtomIdx(),
                                 beginHatomConnections);

      if (beginHatomConnections.empty()) {
        throw FileParseException("No Attach Point for bond");
      }

      getNewAtomsForHydrogenBond(bond->getEndAtom(), bond->getBeginAtomIdx(),
                                 endHatomConnections);
      if (endHatomConnections.empty()) {
        throw FileParseException("No Attach Point for bond");
      }
    } else if (molFromSCSRParams.scsrBaseHbondOptions ==
               SCSRBaseHbondOptions::Auto) {
      // get the hbond sites from the map already created for each
      // base template
      auto templateIdx = atomIdxToTemplateIdx[bond->getBeginAtomIdx()];
      beginHatomConnections = hBondSitesForTemplate[templateIdx];

      templateIdx = atomIdxToTemplateIdx[bond->getEndAtomIdx()];
      endHatomConnections = hBondSitesForTemplate[templateIdx];
    }
    unsigned int mainBeginAtomIdx = bond->getBeginAtomIdx();
    unsigned int mainEndAtomIdx = bond->getEndAtomIdx();

    if (beginHatomConnections.empty() || endHatomConnections.empty()) {
      // no hydrogen bond sites found, so just add the bond between the
      // first atoms on both sides
      return;  // no hydrogen bond sites found, so skip the bond
    }
    if (beginHatomConnections.size() == 3 && endHatomConnections.size() == 3) {
      // see if the two donor flag lists are
      // complementary

      bool complimentary = true;
      for (auto index = 0; index < 3; ++index) {
        // .second is the bool for donors
        if (beginHatomConnections[index].df_isDonor ==
            endHatomConnections[index].df_isDonor) {
          complimentary = false;  // both are donors or both are
                                  // acceptors
          break;
        }
      }

      if (complimentary) {
        addBonds(bond->getBondType(), mainBeginAtomIdx, mainEndAtomIdx,
                 beginHatomConnections, endHatomConnections);
        return;
      }

      // Check for G-U type pairs in a wobble bond
      // configuratio first see if one has the
      // configuration DDA  (like G)

      bool foundDDA = false;
      if (beginHatomConnections[0].df_isDonor &&
          beginHatomConnections[1].df_isDonor &&
          !beginHatomConnections[2].df_isDonor) {
        foundDDA = true;
      } else if (endHatomConnections[0].df_isDonor &&
                 endHatomConnections[1].df_isDonor &&
                 !endHatomConnections[2].df_isDonor) {
        foundDDA = true;
        std::swap(beginHatomConnections, endHatomConnections);
        std::swap(mainBeginAtomIdx, mainEndAtomIdx);
      }

      if (foundDDA) {
        std::vector<HydrogenBondConnection> wobbleBeginAtoms;
        std::vector<HydrogenBondConnection> wobbleEndAtoms;

        wobbleBeginAtoms.push_back(beginHatomConnections[1]);
        wobbleBeginAtoms.push_back(beginHatomConnections[2]);

        // see if the other side has the
        // configuration  ADA (like U)
        if (!endHatomConnections[0].df_isDonor &&
            endHatomConnections[1].df_isDonor &&
            !endHatomConnections[2].df_isDonor) {
          // ADA use the first two atoms
          wobbleEndAtoms.push_back(endHatomConnections[0]);
          wobbleEndAtoms.push_back(endHatomConnections[1]);
          addBonds(bond->getBondType(), mainBeginAtomIdx, mainEndAtomIdx,
                   wobbleBeginAtoms, wobbleEndAtoms);
          return;
        }
      }
    } else if (beginHatomConnections.size() * endHatomConnections.size() ==
               6) /* one is 2 and one is 3*/ {
      if (endHatomConnections.size() == 2) {
        // make sure the first set has two atoms, and the second has 3
        std::swap(beginHatomConnections, endHatomConnections);
        std::swap(mainBeginAtomIdx, mainEndAtomIdx);
      }

      std::vector<HydrogenBondConnection> wobbleEndAtoms;
      if (beginHatomConnections[0].df_isDonor &&
          !beginHatomConnections[1].df_isDonor) {  // like I (DA)
        if (!endHatomConnections[1].df_isDonor &&
            endHatomConnections[2].df_isDonor) {  // 2nd and 3rd are AD
          // like A (DAD) or C (AAD)
          wobbleEndAtoms.push_back(endHatomConnections[1]);
          wobbleEndAtoms.push_back(endHatomConnections[2]);
          addBonds(bond->getBondType(), mainBeginAtomIdx, mainEndAtomIdx,
                   beginHatomConnections, wobbleEndAtoms);
          return;

        } else if (!endHatomConnections[0].df_isDonor &&
                   endHatomConnections[1].df_isDonor) {
          // like U (ADA)
          wobbleEndAtoms.push_back(endHatomConnections[0]);
          wobbleEndAtoms.push_back(endHatomConnections[1]);
          addBonds(bond->getBondType(), mainBeginAtomIdx, mainEndAtomIdx,
                   beginHatomConnections, wobbleEndAtoms);
          return;
        }
      }
    }

    // if here, we have no wobble bond, so just add
    // the bond between the first atoms on both sides

    addBond(bond->getBondType(),
            originAtomMap[OriginAtomDef(
                mainBeginAtomIdx, beginHatomConnections[0].d_templateAtomIdx)],
            originAtomMap[OriginAtomDef(
                mainEndAtomIdx, endHatomConnections[0].d_templateAtomIdx)]);
    return;
  }

 public:
  MolFromSCSRMolConverter(SCSRMol *scsrMolInit,
                          const MolFileParserParams &molFileParserParamsInit,
                          const MolFromSCSRParams &molFromSCSRParamsInit)
      : scsrMol(scsrMolInit),
        molFileParserParams(molFileParserParamsInit),
        molFromSCSRParams(molFromSCSRParamsInit) {}

  std::unique_ptr<RDKit::RWMol> convert() {
    resMol.reset(new RWMol());
    mol = scsrMol->getMol();

    // first get some information from the templates to be used when
    // creating the coords for the new atoms. this is a dirty approach
    // that simply expands the orginal macro atom coords to be big
    // enough to hold any expanded macro atom. No attempt is made to
    // make this look nice, or to avoid overlaps.
    std::vector<RDGeom::Point3D> templateCentroids;
    double maxSize = 0.0;

    const Conformer *conf = nullptr;
    std::unique_ptr<Conformer> newConf(nullptr);
    if (mol->getNumConformers() != 0) {
      conf = &mol->getConformer(0);
      newConf.reset(new Conformer(scsrMol->getMol()->getNumAtoms()));
      newConf->set3D(conf->is3D());

      for (unsigned int templateIdx = 0;
           templateIdx < scsrMol->getTemplateCount(); ++templateIdx) {
        auto templateMol = scsrMol->getTemplate(templateIdx);
        RDGeom::Point3D sumOfCoords;
        const RDKit::Conformer *templateConf = nullptr;
        auto confCount = templateMol->getNumConformers();
        if (confCount == 0) {
          conf = nullptr;
          break;
        }
        templateConf = &templateMol->getConformer(0);
        RDGeom::Point3D maxCoord = templateConf->getAtomPos(0);
        RDGeom::Point3D minCoord = maxCoord;
        for (unsigned int atomIdx = 0; atomIdx < templateMol->getNumAtoms();
             ++atomIdx) {
          auto atomCoord = templateConf->getAtomPos(atomIdx);
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
    }

    // if we are perceiving the H-bonding locations for base templates, do
    // that here

    if (molFromSCSRParams.scsrBaseHbondOptions == SCSRBaseHbondOptions::Auto) {
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
      for (unsigned int templateIdx = 0;
           templateIdx < scsrMol->getTemplateCount(); ++templateIdx) {
        auto templateMol = scsrMol->getTemplate(templateIdx);
        templateMol->updatePropertyCache(false);

        if (templateMol->getProp<std::string>(
                common_properties::molAtomClass) != "BASE") {
          continue;
        }

        hBondSitesForTemplate[templateIdx] =
            std::vector<HydrogenBondConnection>();
        for (const auto &querySmi : hbondQueries) {
          RDKit::SubstructMatchParameters params;
          auto match = SubstructMatch(*templateMol, *querySmi.dp_mol, params);

          if (match.empty()) {
            continue;
          }

          for (const auto &hbondData : querySmi.d_hBondConnections) {
            hBondSitesForTemplate[templateIdx].emplace_back(
                match[0][hbondData.d_templateAtomIdx].second,
                hbondData.df_isDonor);
          }

          break;  // only take the first match
        }
      }
    }

    // for each atom in the main mol, expand it to full atom form

    std::vector<StereoGroup> newStereoGroups;
    std::vector<Atom *> absoluteAtoms;
    std::vector<Bond *> absoluteBonds;

    std::vector<std::unique_ptr<SubstanceGroup>> newSgroups;

    unsigned int atomCount = mol->getNumAtoms();
    for (unsigned int atomIdx = 0; atomIdx < atomCount; ++atomIdx) {
      auto atom = mol->getAtomWithIdx(atomIdx);
      std::string dummyLabel = "";
      std::string atomClass = "";

      if (!atom->getPropIfPresent(common_properties::dummyLabel, dummyLabel) ||
          !atom->getPropIfPresent(common_properties::molAtomClass, atomClass) ||
          dummyLabel == "" || atomClass == "") {
        // NOT a template atom - just copy it
        auto newAtom = new Atom(*atom);
        resMol->addAtom(newAtom, true, true);
        // atomMap[atomIdx].push_back(newAtom->getIdx());

        // templatesFound.push_back(nullptr);
        originAtomMap[OriginAtomDef(atomIdx, UINT_MAX)] = newAtom->getIdx();
        if (conf != nullptr) {
          newConf->setAtomPos(newAtom->getIdx(),
                              conf->getAtomPos(atomIdx) * maxSize);
        }
      } else {  // it is a macro atom - expand it

        unsigned int seqId = 0;
        std::string seqName = "";
        atom->getPropIfPresent(common_properties::molAtomSeqId, seqId);
        atom->getPropIfPresent(common_properties::molAtomSeqName, seqName);

        //  find the template that matches the class and label

        ROMol *templateMol = nullptr;
        unsigned int templateIdx;
        bool templateFound = false;
        std::string templateNameToUse = "";
        for (templateIdx = 0; templateIdx < scsrMol->getTemplateCount();
             ++templateIdx) {
          templateMol = scsrMol->getTemplate(templateIdx);
          std::vector<std::string> templateNames;
          std::string templateName;
          std::string templateAtomClass;
          if (templateMol->getPropIfPresent<std::string>(
                  common_properties::molAtomClass, templateAtomClass) &&
              templateAtomClass == atomClass &&
              templateMol->getPropIfPresent<std::vector<std::string>>(
                  common_properties::templateNames, templateNames) &&
              std::find(templateNames.begin(), templateNames.end(),
                        dummyLabel) != templateNames.end()) {
            templateFound = true;
            switch (molFromSCSRParams.scsrTemplateNames) {
              case SCSRTemplateNames::UseFirstName:
                templateNameToUse = templateNames[0];
                break;
              case SCSRTemplateNames::UseSecondName:
                templateNameToUse = templateNames.back();
                break;
              case SCSRTemplateNames::AsEntered:
                templateNameToUse = dummyLabel;
                break;
            }
            break;
          }
        }

        if (!templateFound) {
          std::ostringstream errout;
          errout << "No template found for atom " << atomIdx << " class "
                 << atomClass << " label " << dummyLabel;
          throw FileParseException(errout.str());
        }

        atomIdxToTemplateIdx[atomIdx] = templateIdx;
        std::vector<std::pair<unsigned int, std::string>> attchOrds;

        atom->getPropIfPresent(common_properties::molAttachOrderTemplate,
                               attchOrds);

        // first find the sgroup that is the base for this atom's
        // template

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

        // add the atoms from the main template to the new molecule
        std::string sgroupName = templateNameToUse;
        if (seqId != 0) {
          sgroupName += "_" + std::to_string(seqId);
        }
        if (seqName != "") {
          sgroupName += "_" + seqName;
        }

        auto coordOffset = (conf->getAtomPos(atomIdx) * maxSize) -
                           templateCentroids[templateIdx];

        copySgroupIntoResult(atomIdx, *sgroup, sgroupName, newSgroups,
                             newConf.get(), coordOffset);

        // find  the attachment points used by this atom and record them
        // in attachMap.

        for (const auto &[idx, lbl] : attchOrds) {
          attachMap[OriginAtomConnection(atomIdx, lbl)] =
              std::vector<unsigned int>();
          auto &attachAtoms = attachMap[OriginAtomConnection(atomIdx, lbl)];
          for (auto attachPoint : sgroup->getAttachPoints()) {
            if (attachPoint.id == lbl) {
              attachAtoms.push_back(attachPoint.aIdx);
            }
          }
          if (attachAtoms.size() == 0) {
            throw FileParseException("Attachment point not found");
          }
        }

        // if we are including atoms from leaving groups, go through the
        // attachment points of the main sgroup. If the attach point is
        // not found in the attachords, then find the sgroup for that
        // attach point and add its atoms the molecule

        if (molFromSCSRParams.includeLeavingGroups) {
          for (auto attachPoint : sgroup->getAttachPoints()) {
            if (attachMap.find(OriginAtomConnection(atomIdx, attachPoint.id)) ==
                attachMap.end()) {
              bool foundLgSgroup = false;

              for (auto lgSgroup : getSubstanceGroups(*templateMol)) {
                std::string lgSup;
                std::string lgSgroupAtomClass;
                if (lgSgroup.getPropIfPresent<std::string>("TYPE", lgSup) &&
                    lgSup == "SUP" &&
                    lgSgroup.getPropIfPresent<std::string>("CLASS",
                                                           lgSgroupAtomClass) &&
                    lgSgroupAtomClass == "LGRP") {
                  auto lgSgroupAtoms = lgSgroup.getAtoms();
                  if (std::find(lgSgroupAtoms.begin(), lgSgroupAtoms.end(),
                                attachPoint.lvIdx) != lgSgroupAtoms.end()) {
                    std::string sgroupName = dummyLabel;
                    if (seqId != 0) {
                      sgroupName += "_" + std::to_string(seqId);
                    }
                    if (seqName != "") {
                      sgroupName += "_" + seqName;
                    }
                    sgroupName += "_" + attachPoint.id;
                    foundLgSgroup = true;
                    copySgroupIntoResult(atomIdx, lgSgroup, sgroupName,
                                         newSgroups, newConf.get(),
                                         coordOffset);

                    break;
                  }
                }
              }

              if (!foundLgSgroup) {
                throw FileParseException("Leaving group SGroup " +
                                         attachPoint.id + " not found");
              }
            }
          }
        }

        // copy the bonds of the template into the new molecule
        // if the bonds are between atoms in the new molecule
        // Bonds to atoms in leaving groups that "left" are NOT copied

        for (auto bond : templateMol->bonds()) {
          if (originAtomMap.find(OriginAtomDef(
                  atomIdx, bond->getBeginAtomIdx())) == originAtomMap.end() ||
              originAtomMap.find(OriginAtomDef(
                  atomIdx, bond->getEndAtomIdx())) == originAtomMap.end()) {
            continue;  // bond not in the new molecule
          }
          auto newBeginAtomIdx =
              originAtomMap[OriginAtomDef(atomIdx, bond->getBeginAtomIdx())];
          auto newEndAtomIdx =
              originAtomMap[OriginAtomDef(atomIdx, bond->getEndAtomIdx())];

          auto newBond = new Bond(bond->getBondType());
          newBond->setBeginAtomIdx(newBeginAtomIdx);
          newBond->setEndAtomIdx(newEndAtomIdx);
          newBond->updateProps(*bond, false);
          resMol->addBond(newBond, true);
        }

        // take care of stereo groups in the template
        // abs groups are added to the list of abs atoms and bonds, so
        // that we can add ONE abs group later

        for (auto &sg : templateMol->getStereoGroups()) {
          std::vector<Atom *> newGroupAtoms;
          std::vector<Bond *> newGroupBonds;

          for (auto sgAtom : sg.getAtoms()) {
            auto originAtom = OriginAtomDef(atomIdx, sgAtom->getIdx());

            auto newAtomPtr = originAtomMap.find(originAtom);
            if (newAtomPtr != originAtomMap.end()) {
              newGroupAtoms.push_back(
                  resMol->getAtomWithIdx(newAtomPtr->second));
            }
          }

          for (auto sgBond : sg.getBonds()) {
            auto originBeginAtom =
                OriginAtomDef(atomIdx, sgBond->getBeginAtomIdx());
            auto originEndAtom =
                OriginAtomDef(atomIdx, sgBond->getEndAtomIdx());

            auto newBeginAtomPtr = originAtomMap.find(originBeginAtom);
            auto newEndAtomPtr = originAtomMap.find(originEndAtom);
            if (newBeginAtomPtr != originAtomMap.end() &&
                newEndAtomPtr != originAtomMap.end()) {
              auto newBond = resMol->getBondBetweenAtoms(
                  newBeginAtomPtr->second, newEndAtomPtr->second);
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
      return std::move(resMol);
    }

    if (conf != nullptr) {
      newConf->resize(resMol->getNumAtoms());
      resMol->addConformer(newConf.release());
    }

    // now deal with the bonds from the original mol.

    for (auto bond : scsrMol->getMol()->bonds()) {
      processBondInMainMol(bond);
    }

    // copy any attrs from the main mol

    for (auto &prop : mol->getPropList(false, false)) {
      std::string propVal;
      if (mol->getPropIfPresent(prop, propVal)) {
        resMol->setProp(prop, propVal);
      }
    }

    // copy the sgroups from the main mol for atoms not in a template

    for (auto &sg : getSubstanceGroups(*mol)) {
      if (sg.getIsValid()) {
        auto &atoms = sg.getAtoms();
        std::vector<unsigned int> newAtoms;
        for (auto atom : atoms) {
          auto originAtom = OriginAtomDef(atom, UINT_MAX);
          auto newAtomPtr = originAtomMap.find(originAtom);
          if (newAtomPtr != originAtomMap.end()) {
            newAtoms.push_back(newAtomPtr->second);
          } else {
            // some atoms were in templates and others were not - cannot
            // add this sgroup
            newAtoms.clear();
            break;
          }
        }
        if (newAtoms.size() > 0) {
          const std::string type = "SUP";
          newSgroups.emplace_back(new SubstanceGroup(resMol.get(), type));
          auto newSg = newSgroups.back().get();
          // RDKit::SubstanceGroup newSg(sg);

          newSg->setAtoms(newAtoms);
        }
      }
    }

    // now that we have all substance groups from the template and from
    // the non-template atoms, and we have all the bonds, find the
    // Xbonds for each substance group and add them

    for (auto bond : resMol->bonds()) {
      for (auto &sg : newSgroups) {
        // if one atom of the bond is found and the other is not in the
        // sgroup, this is a Xbond
        auto sgAtoms = sg->getAtoms();
        if ((std::find(sgAtoms.begin(), sgAtoms.end(),
                       bond->getBeginAtomIdx()) == sgAtoms.end()) !=
            (std::find(sgAtoms.begin(), sgAtoms.end(), bond->getEndAtomIdx()) ==
             sgAtoms.end())) {
          sg->addBondWithIdx(bond->getIdx());
        }
      }
    }

    if (newSgroups.size() > 0) {
      for (auto &sg : newSgroups) {
        addSubstanceGroup(*resMol, *sg.get());
      }
    }
    newSgroups.clear();  // just tidy cleanup

    // take care of stereo groups in the main mol - for atoms that are
    // NOT template refs

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
        auto originBeginAtom =
            OriginAtomDef(sgBond->getBeginAtomIdx(), UINT_MAX);
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

    // make an absolute group that contains any absolute atoms or bonds
    // from either the main mol or the template instantiations

    if (!absoluteAtoms.empty() || !absoluteBonds.empty()) {
      newStereoGroups.emplace_back(StereoGroupType::STEREO_ABSOLUTE,
                                   absoluteAtoms, absoluteBonds);
    }

    if (newStereoGroups.size() > 0) {
      resMol->setStereoGroups(newStereoGroups);
    }

    bool chiralityPossible = false;
    FileParserUtils::finishMolProcessing(resMol.get(), chiralityPossible,
                                         molFileParserParams);

    return std::move(resMol);
  }
};

static std::unique_ptr<RDKit::RWMol> MolFromSCSRMol(
    SCSRMol *scsrMol, const MolFileParserParams &molFileParserParams,
    const MolFromSCSRParams &molFromSCSRParams) {
  MolFromSCSRMolConverter converter(scsrMol, molFileParserParams,
                                    molFromSCSRParams);
  return converter.convert();
}

std::unique_ptr<RDKit::RWMol> MolFromSCSRDataStream(
    std::istream &inStream, unsigned int &line,
    const MolFileParserParams &molFileParserParams,
    const MolFromSCSRParams &molFromSCSRParams) {
  auto scsr = SCSRMolFromSCSRDataStream(inStream, line, molFileParserParams);
  return MolFromSCSRMol(scsr.get(), molFileParserParams, molFromSCSRParams);
}

std::unique_ptr<RDKit::RWMol> MolFromSCSRBlock(
    const std::string &molBlock, const MolFileParserParams &molFileParserParams,
    const MolFromSCSRParams &molFromSCSRParams) {
  auto scsr = SCSRMolFromSCSRBlock(molBlock, molFileParserParams);
  return MolFromSCSRMol(scsr.get(), molFileParserParams, molFromSCSRParams);
}

std::unique_ptr<RDKit::RWMol> MolFromSCSRFile(
    const std::string &fName, const MolFileParserParams &molFileParserParams,
    const MolFromSCSRParams &molFromSCSRParams) {
  auto scsr = SCSRMolFromSCSRFile(fName, molFileParserParams);
  return MolFromSCSRMol(scsr.get(), molFileParserParams, molFromSCSRParams);
}

}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
