//
//  Copyright (C) 2002-2024 Tad Hurst and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "MolSGroupParsing.h"
#include <GraphMol/QueryOps.h>

using namespace RDKit::SGroupParsing;

namespace RDKit {

namespace v2 {
namespace FileParsers {

//------------------------------------------------
//
//  Read a SCVSR molecule from a stream
//
//------------------------------------------------

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
      return originAtomMap.at(OriginAtomDef(
          atomIdx,
          attachMap.at(OriginAtomConnection(atomIdx, attchOrd.getLabel()))));
      break;
    }
  }

  // error attachment ord not found
  throw FileParseException("Attachment ord not found");
}

std::unique_ptr<RWMol> MolFromScsr(const RDKit::SCSRMol &scsrMol) {
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
      for (auto attchOrd : attchOrds) {
        for (auto sgroup : RDKit::getSubstanceGroups(*templateMol)) {
          std::string sup;
          std::string sgroupAtomClass;
          if (sgroup.getPropIfPresent<std::string>("TYPE", sup) &&
              sup == "SUP" &&
              sgroup.getPropIfPresent<std::string>("CLASS", sgroupAtomClass) &&
              sgroupAtomClass == atomClass) {
            for (auto attachPoint : sgroup.getAttachPoints()) {
              if (attachPoint.id == attchOrd.getLabel()) {
                attachMap[OriginAtomConnection(atomIdx, attachPoint.id)] =
                    attachPoint.aIdx;
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
            }
          }
        }
      }

      // add an superatom sgroup to mark the atoms from this macro atom
      // expansion
      std::string typ = "SUP";
      std::unique_ptr<SubstanceGroup> sgroup(
          new SubstanceGroup(resMol.get(), typ));
      sgroup->setProp<unsigned int>("index", seqId);
      sgroup->setProp("LABEL", dummyLabel + "_" + std::to_string(seqId));

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
          sgroup->addAtomWithIdx(newAtom->getIdx());

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
          sgroup->addBondWithIdx(newBond->getIdx());
        }
      }

      if (sgroup->getIsValid()) {
        addSubstanceGroup(*resMol, *sgroup.get());
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
    auto newEndAtomIdx = getNewAtomForBond(
        bond->getEndAtom(), bond->getBeginAtomIdx(), originAtomMap, attachMap);

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
