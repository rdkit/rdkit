
//
//  Copyright (C) 2025 Tad Hurst and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/FileParsers/MACROMolUtils.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <GraphMol/MACROMol.h>
#include <GraphMol/FileParsers/MACROMolUtils.h>

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <RDGeneral/BoostEndInclude.h>
#endif

namespace RDKit {

class MolFromMACROMolConverter {
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

  RDKit::MACROMol *macroMol;
  std::unique_ptr<RWMol> resMol;
  const RDKit::v2::FileParsers::MolFileParserParams molFileParserParams;
  const MolFromMACROMolParams molFromMACROMolParams;

  // maps main atom# and template atom# to new atom#
  std::map<OriginAtomDef, unsigned int> originAtomMap;

  // maps main atom# and attach label to template atom#
  std::map<OriginAtomConnection, unsigned int> attachMap;

  unsigned int getNewAtomForBond(const Atom *atom, unsigned int otherAtomIdx) {
    std::string atomClass = "";
    unsigned int atomIdx = atom->getIdx();
    if (!atom->getPropIfPresent<std::string>(common_properties::molAtomClass,
                                             atomClass)) {
      return originAtomMap.at(OriginAtomDef(atomIdx, UINT_MAX));
    }

    // if here , it is a template atom

    std::vector<std::pair<unsigned int, std::string>> attchOrds;
    atom->getProp(common_properties::molAttachOrderTemplate, attchOrds);
    for (const auto &[idx, lbl] : attchOrds) {
      if (idx == otherAtomIdx) {
        auto attachMapIt = attachMap.find(OriginAtomConnection(atomIdx, lbl));
        if (attachMapIt == attachMap.end()) {
          throw FileParseException("Attachment ord not found");
        }
        return originAtomMap.at(OriginAtomDef(atomIdx, attachMapIt->second));
      }
    }

    // error attachment ord not found
    throw FileParseException("Attachment ord not found");
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
                      macroMol->atomIdxToTemplateMol(atomIdx)->getNumAtoms());
    }

    for (auto templateAtomIdx : sgroup.getAtoms()) {
      auto templateAtom =
          macroMol->atomIdxToTemplateMol(atomIdx)->getAtomWithIdx(
              templateAtomIdx);
      auto newAtom = new Atom(*templateAtom);

      resMol->addAtom(newAtom, true, true);
      newSgroup->addAtomWithIdx(newAtom->getIdx());

      originAtomMap[OriginAtomDef(atomIdx, templateAtomIdx)] =
          newAtom->getIdx();

      if (newConf) {
        newConf->setAtomPos(
            newAtom->getIdx(),
            coordOffset + macroMol->atomIdxToTemplateMol(atomIdx)
                              ->getConformer()
                              .getAtomPos(templateAtomIdx));
      }
    }
  }

  void addBond(const Bond::BondType bondType, unsigned int beginAtom,
               unsigned int endAtom) {
    resMol->addBond(beginAtom, endAtom, bondType);
  }

  void processBondInMainMol(const Bond *bond) {
    unsigned int newBeginAtom;
    unsigned int newEndAtom;

    newBeginAtom =
        getNewAtomForBond(bond->getBeginAtom(), bond->getEndAtomIdx());
    if (newBeginAtom == UINT_MAX) {
      throw FileParseException("Error getting new atom for bond");
    }

    newEndAtom = getNewAtomForBond(bond->getEndAtom(), bond->getBeginAtomIdx());
    if (newEndAtom == UINT_MAX) {
      throw FileParseException("Error getting new atom for bond");
    }

    addBond(bond->getBondType(), newBeginAtom, newEndAtom);
    return;
  }

 public:
  MolFromMACROMolConverter(
      MACROMol *macroMolInit,
      const RDKit::v2::FileParsers::MolFileParserParams
          &molFileParserParamsInit,
      const MolFromMACROMolParams &molFromMACROMolParamsInit)
      : macroMol(macroMolInit),
        molFileParserParams(molFileParserParamsInit),
        molFromMACROMolParams(molFromMACROMolParamsInit) {}

  std::unique_ptr<RDKit::RWMol> convert() {
    resMol.reset(new RWMol());

    // first get some information from the templates to be used when
    // creating the coords for the new atoms. this is a dirty approach
    // that simply expands the orginal macro atom coords to be big
    // enough to hold any expanded macro atom. No attempt is made to
    // make this look nice, or to avoid overlaps.
    std::vector<RDGeom::Point3D> templateCentroids;
    double maxSize = 0.0;

    const Conformer *conf = nullptr;
    std::unique_ptr<Conformer> newConf(nullptr);
    if (macroMol->getNumConformers() != 0) {
      conf = &macroMol->getConformer(0);
      newConf.reset(new Conformer(macroMol->getNumAtoms()));
      newConf->set3D(conf->is3D());

      for (unsigned int templateIdx = 0;
           templateIdx < macroMol->getTemplateCount(); ++templateIdx) {
        auto templateMol = macroMol->getTemplate(templateIdx);
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

    // loop over the bonds in the main mol.  For each bond, reocrd the ATTACHMAP
    // entries.  We do not know the mapped template atom yet, so add a
    // placeholder for that.  We will fill in the real template atom after we
    // have processed the main atoms and know the mapping of main atoms to
    // template atoms.

    for (auto &bond : macroMol->bonds()) {
      std::string lbl;
      if (bond->getPropIfPresent(common_properties::_MolFileBondAttachPt1,
                                 lbl)) {
        auto key = OriginAtomConnection(bond->getBeginAtomIdx(), lbl);
        if (attachMap.find(key) != attachMap.end()) {
          throw FileParseException("Duplicate attachment point label");
        }
        attachMap[key] = UINT_MAX;
      }
      if (bond->getPropIfPresent(common_properties::_MolFileBondAttachPt2,
                                 lbl)) {
        auto key = OriginAtomConnection(bond->getEndAtomIdx(), lbl);
        if (attachMap.find(key) != attachMap.end()) {
          throw FileParseException("Duplicate attachment point label");
        }
        attachMap[key] = UINT_MAX;  // placeholder for now, will fill in real
                                    // template atom idx later
      }
    }

    // for each atom in the main mol, expand it to full atom form

    std::vector<StereoGroup> newStereoGroups;
    std::vector<Atom *> absoluteAtoms;
    std::vector<Bond *> absoluteBonds;

    std::vector<std::unique_ptr<SubstanceGroup>> newSgroups;

    unsigned int atomCount = macroMol->getNumAtoms();
    for (unsigned int atomIdx = 0; atomIdx < atomCount; ++atomIdx) {
      auto atom = macroMol->getAtomWithIdx(atomIdx);
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
        for (templateIdx = 0; templateIdx < macroMol->getTemplateCount();
             ++templateIdx) {
          templateMol = macroMol->getTemplate(templateIdx);
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
            switch (molFromMACROMolParams.macroTemplateNames) {
              case MACROTemplateNames::UseFirstName:
                templateNameToUse = templateNames[0];
                break;
              case MACROTemplateNames::UseSecondName:
                templateNameToUse = templateNames.back();
                break;
              case MACROTemplateNames::AsEntered:
                templateNameToUse = dummyLabel;
                break;
              case MACROTemplateNames::All:
                templateNameToUse = "";
                for (const auto &nm : templateNames) {
                  if (templateNameToUse != "") {
                    templateNameToUse += "+";
                  }
                  templateNameToUse += nm;
                }
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

        std::string sgroupName = atomClass + "_";
        if (seqId != 0) {
          sgroupName += std::to_string(seqId) + "_";
        } else {
          sgroupName += "na_";
        }
        if (seqName != "") {
          sgroupName += "_" + seqName;
        }

        sgroupName += templateNameToUse;

        auto coordOffset = (conf->getAtomPos(atomIdx) * maxSize) -
                           templateCentroids[templateIdx];

        copySgroupIntoResult(atomIdx, *sgroup, sgroupName, newSgroups,
                             newConf.get(), coordOffset);

        // if we are including atoms from leaving groups, go through the
        // attachment points of the main sgroup. If the attach point is
        // not found in the attachMap, then find the sgroup for that
        // attach point and add its atoms the molecule

        if (molFromMACROMolParams.includeLeavingGroups) {
          for (auto attachPoint : sgroup->getAttachPoints()) {
            auto key = OriginAtomConnection(atomIdx, attachPoint.id);
            if (attachMap.find(key) != attachMap.end()) {
              // fill in the tempate atom id for this attachPoint

              attachMap[key] = attachPoint.aIdx;
            } else {
              // this attach point was not found, so the leaving group is
              // included in the output molecule (if there is one).

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
                    copySgroupIntoResult(atomIdx, lgSgroup, sgroupName,
                                         newSgroups, newConf.get(),
                                         coordOffset);

                    break;
                  }
                }
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

    for (auto bond : macroMol->bonds()) {
      processBondInMainMol(bond);
    }

    // copy any attrs from the main mol

    for (auto &prop : macroMol->getPropList(false, false)) {
      std::string propVal;
      if (macroMol->getPropIfPresent(prop, propVal)) {
        resMol->setProp(prop, propVal);
      }
    }

    // copy the sgroups from the main mol for atoms not in a template

    for (auto &sg : getSubstanceGroups(*macroMol)) {
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

          newSg->updateProps(sg, false);
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

    for (auto &sg : macroMol->getStereoGroups()) {
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
    RDKit::FileParserUtils::finishMolProcessing(resMol.get(), chiralityPossible,
                                                molFileParserParams);

    return std::move(resMol);
  }
};

static std::unique_ptr<RDKit::RWMol> MolFromMACROMol(
    MACROMol *macroMol,
    const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams,
    const RDKit::MolFromMACROMolParams &molFromMACROMolParams) {
  MolFromMACROMolConverter converter(macroMol, molFileParserParams,
                                     molFromMACROMolParams);
  return converter.convert();
}

}  // end of namespace RDKit