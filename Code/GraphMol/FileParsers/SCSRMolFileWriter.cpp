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
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/SCSRMol.h>

using namespace RDKit::SGroupParsing;

namespace RDKit {
//------------------------------------------------
//
//  Read a SCVSR molecule from a stream
//
//------------------------------------------------

std::string SCSRMolToSCSRMolBlock(SCSRMol &scsrMol,
                                  const RDKit::MolWriterParams &params,
                                  int confId) {
  RDKit::Utils::LocaleSwitcher switcher;
  std::string res = "";
  auto localParams = params;
  localParams.forceV3000 = true;
  localParams.writeEndMolLine = false;

  res += RDKit::MolToMolBlock(*(scsrMol.getMol()), localParams, confId);

  // now write out all of the templates

  res += "M  V30 BEGIN TEMPLATE\n";
  for (unsigned int templateId = 0; templateId < scsrMol.getTemplateCount();
       ++templateId) {
    auto scsrTemplate = scsrMol.getTemplate(templateId);

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

    // scsrTemplate->getPropIfPresent<std::string>(common_properties::natReplace,
    //                                             natReplace);
    // scsrTemplate->getPropIfPresent<std::string>(
    //     common_properties::molTemplateComment, comment);
    // scsrTemplate->getPropIfPresent<std::string>(
    //     common_properties::molTemplateFullName, fullname);
    // scsrTemplate->getPropIfPresent<std::string>(
    //     common_properties::molTemplateCategory, category);
    // scsrTemplate->getPropIfPresent<std::string>(
    //     common_properties::molTemplateUniqueId, uniqueId);
    // scsrTemplate->getPropIfPresent<std::string>(
    //     common_properties::molTemplateCasNumber, casNumber);
    // scsrTemplate->getPropIfPresent<std::string>(
    //     common_properties::molTemplateCollaborator, collaborator);
    // scsrTemplate->getPropIfPresent<std::string>(
    //     common_properties::molTemplateProtection, protection);

    scsrTemplate->getProp(common_properties::molAtomClass, templateClass);

    scsrTemplate->getProp(common_properties::templateNames, templateNames);

    res += "M  V30 TEMPLATE " + std::to_string(templateId + 1) + " " +
           templateClass;
    for (auto templateName : templateNames) {
      res += "/" + templateName;
    }
    res += "/";

    for (auto &propName : scsrTemplate->getPropList()) {
      if (propName != common_properties::molAtomClass &&
          propName != common_properties::templateNames && propName[0] != '_') {
        std::string propVal;
        scsrTemplate->getProp<std::string>(propName, propVal);
        res += " " + propName + "=" + propVal;
      }
    }

    res += "\n";

    boost::dynamic_bitset<> aromaticBonds(scsrTemplate->getNumBonds());
    RWMol tMol(*scsrTemplate);
    prepareMol(tMol, localParams, aromaticBonds);
    res += FileParserUtils::getV3000CTAB(tMol, aromaticBonds, confId,
                                         localParams.precision);
  }

  res += "M  V30 END TEMPLATE\n";
  res += "M  END\n";

  return res;
}

RDKIT_FILEPARSERS_EXPORT void SCSRMolToSCSRMolFile(
    RDKit::SCSRMol &scsrMol, const std::string &fName,
    const RDKit::MolWriterParams &params, int confId) {
  auto *outStream = new std::ofstream(fName.c_str());
  if (!(*outStream) || outStream->bad()) {
    delete outStream;
    std::ostringstream errout;
    errout << "Bad output file " << fName;
    throw BadFileException(errout.str());
  }
  std::string outString = SCSRMolToSCSRMolBlock(scsrMol, params, confId);
  *outStream << outString;
  delete outStream;
}

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::SCSRMol>
SCSRAtributedMolToTemplate(const ROMol &mol) {
  auto res = std::unique_ptr<RDKit::SCSRMol>(new RDKit::SCSRMol());

  boost::dynamic_bitset<> atomsDone(mol.getNumAtoms());
  boost::dynamic_bitset<> bondsDone(mol.getNumBonds());

  // got through the SGROUPS and find all of the SUP groups that define a
  // template

  for (auto &sgroup : RDKit::getSubstanceGroups(mol)) {
    if (sgroup.getProp<std::string>("TYPE") == "SUP") {
      std::string atomClass = sgroup.getProp<std::string>("CLASS");

      // // make a new template molecule
      // auto templateMol = std::unique_ptr<ROMol>(new ROMol());

      // std::map<unsigned int, unsigned int> oldToNewAtomIdx;
      // std::map<unsigned int, unsigned int> newToOldAtomIdx;

      // for (auto match : matchVect) {
      //   // add the atoms
      //   for (auto atomIdx : sgroup.getAtoms()) {
      //     if (!atomsDone[atomIdx]) {
      //       auto atom = mol.getAtomWithIdx(atomIdx);
      //       auto newAtom = new Atom(*atom);
      //       unsigned int newAtomIdx =
      //           templateMol->addAtom(newAtom, false, true);
      //       oldToNewAtomIdx[atomIdx] = newAtomIdx;
      //       newToOldAtomIdx[newAtomIdx] = atomIdx;
      //       atomsDone[atomIdx] = 1;
      //     }
      //   }

      //   // add the bonds
      //   for (auto bondIdx : sgroup.getBonds()) {
      //     if (!bondsDone[bondIdx]) {
      //       auto bond = mol.getBondWithIdx(bondIdx);
      //       unsigned int begAtomIdx = bond->getBeginAtomIdx();
      //       unsigned int endAtomIdx = bond->getEndAtomIdx();
      //       PRECONDITION(
      //           oldToNewAtomIdx.find(begAtomIdx) != oldToNewAtomIdx.end(),
      //           "no mapping for bond begin atom");
      //       PRECONDITION(
      //           oldToNewAtomIdx.find(endAtomIdx) != oldToNewAtomIdx.end(),
      //           "no mapping for bond end atom");
      //       unsigned int newBegAtomIdx = oldToNewAtomIdx[begAtomIdx];
      //       unsigned int newEndAtomIdx = oldToNew

      //           return res;
    }
  }
  return res;
}

bool isTemplateMatchAHit(
    const MatchVectType &match, const ROMol &mol, const ROMol *templateMol,
    const RWMol *queryMol, std::map<unsigned int, unsigned int> &atomMap,
    const std::vector<SubstanceGroup::AttachPoint> &supAttachPoints,
    std::vector<std::pair<unsigned int, std::string>> &newAttachOrds,
    std::vector<unsigned int> &atomsInMatch) {
  // get a new set of match points to be used for the new tempalte ref atoms
  atomsInMatch.clear();
  newAttachOrds.clear();
  std::map<unsigned int, unsigned int> molToQueryMap;

  // check to see if any of the atoms or bonds are already used
  for (auto pair : match) {
    if (atomMap.contains(pair.second)) {
      // atom already used in another match - skip this match
      return false;
    }

    molToQueryMap[pair.second] = pair.first;
    atomsInMatch.push_back(pair.second);
  }
  // for each atom in the match, check its bonds.  If it has a bond to
  // another atom in the hit,
  //  make sure that bond is the query.  If not, skip this match.

  // If the atom is an attachment point for the template, see if there is
  // a bond to an atom NOT in the match. If so, record the attachment
  // order for later.

  // if the atom is not an attachment point,  it has an external bond:
  //  if the external bond is a hydrogen bond, record that attachment.
  //  if the external bond is NOT a hydrogen bond, skip this match.

  for (auto molToQueryMapItem : molToQueryMap) {
    auto atom = mol.getAtomWithIdx(molToQueryMapItem.first);
    auto queryAtom = queryMol->getAtomWithIdx(molToQueryMapItem.second);
    unsigned int templateAtomIdx =
        queryAtom->getProp<unsigned int>("origAtomId");
    for (auto nbri : boost::make_iterator_range(mol.getAtomBonds(atom))) {
      auto bond = mol[nbri];
      unsigned int nbrAtomIdx = bond->getOtherAtomIdx(atom->getIdx());
      if (molToQueryMap.contains(nbrAtomIdx)) {
        // make sure the bond is also in the query mol - if not skip this
        // match

        auto queryBond = queryMol->getBondBetweenAtoms(
            molToQueryMap[bond->getBeginAtomIdx()],
            molToQueryMap[bond->getEndAtomIdx()]);

        if (!queryBond) {
          return false;
        }

        // bond to another atom is in the query - ok
        continue;
      } else {
        // bond to an atom NOT in the main SUP - could be a leaving group

        bool sapFound = false;
        for (auto &attachPoint : supAttachPoints) {
          if (attachPoint.aIdx == templateAtomIdx) {
            // see if the attachment point is an H or OH, and the nbr atom is
            // also an H or OH if so, this leaving group can be considered part
            // of the template

            sapFound = true;
            auto nbrAtom = mol.getAtomWithIdx(nbrAtomIdx);
            const Atom *templateNbrAtom = nullptr;
            if (attachPoint.lvIdx >= 0) {
              templateNbrAtom = templateMol->getAtomWithIdx(attachPoint.lvIdx);
            }

            if (attachPoint.lvIdx >= 0 &&
                ((nbrAtom->getAtomicNum() == 1 &&
                  nbrAtom->getTotalDegree() == 1) ||
                 (nbrAtom->getAtomicNum() == 8 &&
                  nbrAtom->getTotalDegree() == 2)) &&
                templateNbrAtom->getAtomicNum() == nbrAtom->getAtomicNum() &&
                templateNbrAtom->getTotalDegree() ==
                    nbrAtom->getTotalDegree()) {
              // ok - treat this as part of the template
              atomsInMatch.push_back(nbrAtomIdx);
            } else {
              newAttachOrds.push_back(
                  std::make_pair(nbrAtomIdx, attachPoint.id));
            }
            break;
          }
        }

        if (!sapFound) {
          // not an attachment point - see if it is a hydrogen bond - a H-bond
          // is OK, but any other type is NOT OK
          if (bond->getBondType() != Bond::HYDROGEN) {
            {
              // external bond that is not allowed - skip this match
              return false;
            }
          }
        }
      }
    }
  }

  return true;
}

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::SCSRMol> MolToScsrMol(
    const ROMol &mol, RDKit::SCSRMol &templates,
    MolToSCSRParams molToSCSRParams) {
  auto res = std::unique_ptr<SCSRMol>(new SCSRMol());
  auto resMol = std::unique_ptr<RWMol>(new RWMol());

  std::map<unsigned int, unsigned int> atomMap;

  for (unsigned int templateIndex = 0;
       templateIndex < templates.getTemplateCount(); ++templateIndex) {
    auto templateMol = templates.getTemplate(templateIndex);
    templateMol->updatePropertyCache();
    std::vector<std::string> templateNames;

    std::string templateAtomClass;
    std::string templateNameToUse;
    templateMol->getPropIfPresent<std::string>(common_properties::molAtomClass,
                                               templateAtomClass);
    templateMol->getPropIfPresent<std::vector<std::string>>(
        common_properties::templateNames, templateNames);
    switch (molToSCSRParams.scsrUseTemplateName) {
      case SCSRUseTemplateName::UseFirstName:
        templateNameToUse = templateNames[0];
        break;
      case SCSRUseTemplateName::UseSecondName:
        templateNameToUse = templateNames.back();
        break;
    }
    std::vector<std::pair<unsigned int, std::string>> attachOrds;
    std::vector<std::pair<unsigned int, std::string>> newAttachOrds;

    // first find the sgroup that is the base for this template

    const SubstanceGroup *sgroup = nullptr;
    for (auto &sgroupToTest : RDKit::getSubstanceGroups(*templateMol)) {
      std::string sup;
      std::string sgroupAtomClass;
      if ((!sgroupToTest.getPropIfPresent<std::string>("TYPE", sup) &&
           sup != "SUP") ||
          (!sgroupToTest.getPropIfPresent<std::string>("CLASS",
                                                       sgroupAtomClass) ||
           std::find(ScsrClasses.begin(), ScsrClasses.end(), sgroupAtomClass) ==
               ScsrClasses.end())) {
        continue;  // sgroup is not a sup for a template
      }
      if (sgroupAtomClass != "LGRP") {
        sgroup = &sgroupToTest;
        break;
      }
    }

    if (sgroup == nullptr) {
      continue;  // skip this template
    }
    auto supAttachPoints = sgroup->getAttachPoints();

    std::unique_ptr<RWMol> queryMol(new RWMol(*templateMol));

    queryMol->beginBatchEdit();
    for (auto atom : queryMol->atoms()) {
      atom->setProp("origAtomId", atom->getIdx());

      if (std::find(sgroup->getAtoms().begin(), sgroup->getAtoms().end(),
                    atom->getIdx()) == sgroup->getAtoms().end()) {
        queryMol->removeAtom(atom->getIdx());
      }
    }
    queryMol->commitBatchEdit();
    auto queryOut = RDKit::MolToMolBlock(*queryMol);

    // find all occurances of the template in the molecule

    MolOps::sanitizeMol(*queryMol);

    queryMol->updatePropertyCache();

    std::vector<MatchVectType> matchVect;

    SubstructMatchParameters params;
    params.recursionPossible = false;
    params.useChirality = true;
    params.useQueryQueryMatches = false;
    params.maxMatches = 0;  // find all matches

    matchVect = SubstructMatch(mol, *queryMol, params);

    if (!matchVect.size()) {
      continue;
    }

    bool templateCopied = false;
    for (const auto &match : matchVect) {
      std::vector<unsigned int> atomsInMatch;

      // add this match to the SCSRMol if it is a valid hit

      if (!isTemplateMatchAHit(match, mol, templateMol, queryMol.get(), atomMap,
                               supAttachPoints, newAttachOrds, atomsInMatch)) {
        continue;
      }

      if (!templateCopied) {
        // add the template to the SCSRMol
        res->addTemplate(std::unique_ptr<ROMol>(new ROMol(*templateMol)));
        templateCopied = true;
      }

      auto newAtomIdx = resMol->getNumAtoms();
      for (const auto atomMatch : atomsInMatch) {
        atomMap[atomMatch] = newAtomIdx;
      }

      // create the macro atom reference

      resMol->addAtom(new Atom(0), true, true);
      auto resAtom = resMol->getAtomWithIdx(newAtomIdx);

      resAtom->setAtomicNum(0);
      resAtom->setProp(common_properties::dummyLabel, templateNameToUse);
      resAtom->setProp(common_properties::molAtomClass, templateAtomClass);

      // add the attach ords

      if (newAttachOrds.size() > 0) {
        std::vector<std::pair<unsigned int, std::string>> attchOrds;
        for (const auto &newSap : newAttachOrds) {
          attchOrds.emplace_back(newSap.first, newSap.second);
        }
        resAtom->setProp(common_properties::molAttachOrderTemplate, attchOrds);
      }
    }
  }

  // add all atoms that are not in a template

  for (const auto atom : mol.atoms()) {
    if (!atomMap.contains(atom->getIdx())) {
      auto newAtomIdx = resMol->getNumAtoms();
      resMol->addAtom(new Atom(*atom), true, true);
      atomMap[atom->getIdx()] = newAtomIdx;
    }
  }

  // fix the saps - they were previouisly added to include the atomIds from
  // the original mol, but should be the new atom Ids in the SCSRMol

  for (const auto atom : resMol->atoms()) {
    if (!atom->hasProp(common_properties::molAttachOrderTemplate)) {
      continue;  // not a template reference atom
    }
    std::vector<std::pair<unsigned int, std::string>> attachOrds;
    std::vector<std::pair<unsigned int, std::string>> newAttachOrds;

    atom->getProp(common_properties::molAttachOrderTemplate, attachOrds);

    for (const auto &attachOrd : attachOrds) {
      PRECONDITION(atomMap.contains(attachOrd.first),
                   "no mapping for attachment point atom");
      auto attachOrdFound = false;
      for (auto newAttachOrd : newAttachOrds) {
        if (newAttachOrd.first == atomMap[attachOrd.first]) {
          attachOrdFound = true;
          break;
        }
      }
      if (!attachOrdFound) {
        newAttachOrds.emplace_back(atomMap[attachOrd.first], attachOrd.second);
      }
    }

    atom->setProp(common_properties::molAttachOrderTemplate, newAttachOrds);
  }

  // bonds:  if the atoms of the original bonds are mapped to the same
  // atom, they are in the template def

  for (auto bond : mol.bonds()) {
    unsigned int begAtomIdx = bond->getBeginAtomIdx();
    unsigned int endAtomIdx = bond->getEndAtomIdx();

    if (atomMap.contains(begAtomIdx) && atomMap.contains(endAtomIdx) &&
        atomMap[begAtomIdx] == atomMap[endAtomIdx]) {
      continue;
    }

    // bonds are in different new atoms - add the bond
    // check for Hydrogen bonds - they can be multiple hbonds in the atomistic
    // mol, but only one hbond in the SCSRMol
    if (bond->getBondType() == Bond::HYDROGEN) {
      if (resMol->getBondBetweenAtoms(atomMap[begAtomIdx],
                                      atomMap[endAtomIdx])) {
        continue;
      }
    }

    resMol->addBond(atomMap[begAtomIdx], atomMap[endAtomIdx],
                    bond->getBondType());
  }

  res->setMol(std::move(resMol));

  return res;
}

}  // namespace RDKit
