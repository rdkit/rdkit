//
//  Copyright (C) 2022-2023 Tad Hurst, Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// this file parses MRV file for molecules and reactions

#include <string>
#include <exception>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cfloat>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSGroupParsing.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include "MarvinParser.h"
#include "MarvinDefs.h"
#include <GraphMol/Conformer.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Atropisomers.h>
#include <GraphMol/Chirality.h>

#include <GraphMol/RDKitQueries.h>
#include <GraphMol/StereoGroup.h>
#include <GraphMol/SubstanceGroup.h>

#include "MarvinParser.h"

#include <RDGeneral/StreamOps.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>

namespace RDKit {

namespace v2 {
namespace MarvinParser {

/*
      Imports the Marvin-specific dialect of CML (Chemical Markup Language) and
   converts it to datastructures that are compatible with Molfile, RXNfile, and
   Molfile complemented with canvas objects.
*/
class MarvinCMLReader {
 public:
  MarvinCMLReader(){};

  ~MarvinCMLReader(){};

  std::unique_ptr<RWMol> parseMolecule(boost::property_tree::ptree molTree,
                       bool sanitize = false, bool removeHs = false) {
    boost::property_tree::ptree molSection;

    try {
      molSection = molTree.get_child("cml.MDocument.MChemicalStruct.molecule");
    } catch (const std::exception &) {
      try {
        molSection = molTree.get_child("cml.MDocument");
        return std::make_unique<RWMol>();
      } catch (const std::exception &) {
        try {
          molSection = molTree.get_child("cml");
          return std::make_unique<RWMol>();
        } catch (const std::exception &) {
          throw FileParseException("Expected \"molecule\" in MRV file");
        }
      }
    }

    std::unique_ptr<MarvinMol> marvinMol{
        (MarvinMol *)parseMarvinMolecule(molSection)};

    marvinMol->prepSgroupsForRDKit();

    return parseMolecule(marvinMol.get(), sanitize, removeHs);
  }

  std::unique_ptr<ChemicalReaction> parseReaction(
      boost::property_tree::ptree rxnTree,
      boost::property_tree::ptree documentTree, bool sanitize = false,
      bool removeHs = false) {
    std::unique_ptr<ChemicalReaction> rxn{new ChemicalReaction()};
    rxnTree = rxnTree.get_child("cml.MDocument.MChemicalStruct.reaction");
    std::unique_ptr<MarvinReaction> marvinReaction{
        parseMarvinReaction(rxnTree, documentTree)};
    marvinReaction->prepSgroupsForRDKit();

    // get each reactant
    for (auto &mol : marvinReaction->reactants) {
      rxn->addReactantTemplate(
          ROMOL_SPTR(parseMolecule(mol.get(), sanitize, removeHs)));
    }

    // get each agent
    for (auto &mol : marvinReaction->agents) {
      rxn->addAgentTemplate(
          ROMOL_SPTR(parseMolecule(mol.get(), sanitize, removeHs)));
    }

    // get each product
    for (auto &mol : marvinReaction->products) {
      rxn->addProductTemplate(
          ROMOL_SPTR(parseMolecule(mol.get(), sanitize, removeHs)));
    }

    // convert atoms to queries:
    for (auto mol : rxn->getReactants()) {
      for (auto atom : mol->atoms()) {
        QueryOps::replaceAtomWithQueryAtom(static_cast<RWMol *>(mol.get()),
                                           atom);
      }
    }
    for (auto mol : rxn->getProducts()) {
      for (auto atom : mol->atoms()) {
        QueryOps::replaceAtomWithQueryAtom(static_cast<RWMol *>(mol.get()),
                                           atom);
      }
    }

    marvinReaction->prepSgroupsForRDKit();

    return rxn;
  }

  Atom *molAtomFromMarvinAtom(const MarvinAtom *marvinAtom,
                              const MarvinMolBase *marvinMolBase) {
    PRECONDITION(marvinAtom, "bad marvin atom");
    PRECONDITION(marvinMolBase, "bad marvin mol");
    Atom *res = nullptr;

    try {
      std::string symb = marvinAtom->elementType;
      boost::trim(symb);

      if (symb == "*") {
        auto *query = new QueryAtom(0);  // wHAT IS AUTO *
        res = query;
        query->setQuery(makeAtomNullQuery());
        query->setProp(common_properties::dummyLabel, "*");
        // queries have no implicit Hs:
        res->setNoImplicit(true);
      } else if (symb == "R") {
        // It must have an rgroupRef/>

        if (marvinAtom->rgroupRef < 0) {
          throw FileParseException(
              "Expected an R group atom to have an rgroupRef designation");
        }

        auto *query = new QueryAtom(0);  // wHAT IS AUTO *
        res = query;

        query->setProp(common_properties::_MolFileRLabel,
                       marvinAtom->rgroupRef);
        std::string dLabel = "R" + std::to_string(marvinAtom->rgroupRef);
        query->setProp(common_properties::dummyLabel, dLabel);
        if (marvinAtom->mrvAlias != "") {
          query->setProp(common_properties::molFileAlias, marvinAtom->mrvAlias);
        }
        query->setQuery(makeAtomNullQuery());
      }

      else if (symb.size() <= 2) {
        // must be a regular atom

        if (symb.size() == 2 && symb[1] >= 'A' && symb[1] <= 'Z') {
          symb[1] = static_cast<char>(tolower(symb[1]));
        }
        res = new Atom();

        try {
          res->setAtomicNum(PeriodicTable::getTable()->getAtomicNumber(symb));
        } catch (const Invar::Invariant &e) {
          delete res;
          res = nullptr;
          throw FileParseException(e.what());
        }
      } else {
        throw FileParseException("Unrecognized atom type in MRV file");
      }

      // res->setPos(marvinAtom.x2,marvinAtom->y2,0.0);

      if (marvinAtom->formalCharge != 0) {
        res->setFormalCharge(marvinAtom->formalCharge);
      }

      if (marvinAtom->isotope != 0) {
        res->setIsotope(marvinAtom->isotope);
      }

      // mrvValence and hydrogenCount both set the number of explicit hydrogens
      // if both are present, they must agreee

      if (marvinAtom->hydrogenCount >= 0) {
        res->setNumExplicitHs(marvinAtom->hydrogenCount);
      }

      if (marvinAtom->mrvValence >= 0) {
        int explicitValence = marvinMolBase->getExplicitValence(*marvinAtom);
        // if (explicitValence > marvinAtom->mrvValence)
        //   throw FileParseException("atom has specified valence (" +
        //   std::to_string(marvinAtom->mrvValence) +  ") smaller than the drawn
        //   valence");

        int hCount = marvinAtom->mrvValence - explicitValence;
        if (hCount >= 0) {
          if (marvinAtom->hydrogenCount >= 0) {
            if (marvinAtom->hydrogenCount != hCount) {
              throw FileParseException(
                  "mrvValence and hydrogenCount both specified for an atom, and they do not agree");
            }
          } else {
            res->setNoImplicit(true);
            res->setNumExplicitHs(hCount);
          }
        }
      }

      if (marvinAtom->radical != "") {
        res->setNumRadicalElectrons(
            marvinRadicalToRadicalElectrons.at(marvinAtom->radical));
      }

      if (marvinAtom->mrvMap != 0) {
        res->setProp(common_properties::molAtomMapNumber, marvinAtom->mrvMap);
      }

      return res;
    } catch (const std::exception &) {
      delete res;
      res = nullptr;
      throw;
    }
  }

  void molBondFromMarvinBond(const MarvinBond *marvinBond,
                             const MarvinMol *marvinMol, RWMol *mol,
                             bool &chiralityPossible) {
    PRECONDITION(marvinBond, "bad marvin bond");
    PRECONDITION(marvinMol, "bad marvin mol");

    unsigned int bType = UINT_MAX;
    Bond *bond = nullptr;

    try {
      int idx1 = marvinMol->getAtomIndex(marvinBond->atomRefs2[0]);
      int idx2 = marvinMol->getAtomIndex(marvinBond->atomRefs2[1]);

      Bond::BondType type = Bond::UNSPECIFIED;
      ;

      // if there is a query type, that takes precidence over the bond type

      std::string marvinBondType = marvinBond->getBondType();

      if (marvinBondType == "SD") {
        bType = 5;
        type = Bond::UNSPECIFIED;
        bond = new QueryBond;
        bond->setQuery(makeSingleOrDoubleBondQuery());
      } else if (marvinBondType == "SA") {
        bType = 6;
        type = Bond::UNSPECIFIED;
        bond = new QueryBond;
        bond->setQuery(makeSingleOrAromaticBondQuery());
      } else if (marvinBondType == "DA") {
        bType = 7;
        type = Bond::UNSPECIFIED;
        bond = new QueryBond;
        bond->setQuery(makeDoubleOrAromaticBondQuery());
      } else if (marvinBondType == "ANY") {
        bType = 8;
        type = Bond::UNSPECIFIED;
        bond = new QueryBond;
        bond->setQuery(makeBondNullQuery());
      } else if (marvinBondType == "DATIVE") {
        bond = new Bond();
        type = Bond::DATIVE;
        bType = 1;
      } else if (marvinBondType == "1") {
        type = Bond::SINGLE;
        bond = new Bond;
        bType = 1;
      } else if (marvinBondType == "2") {
        type = Bond::DOUBLE;
        bond = new Bond;
        bType = 2;
      } else if (marvinBondType == "3") {
        type = Bond::TRIPLE;
        bond = new Bond;
        bType = 3;
      } else if (marvinBondType == "A") {
        bType = 4;
        type = Bond::AROMATIC;
        bond = new Bond;
      }

      bond->setBeginAtomIdx(idx1);
      bond->setEndAtomIdx(idx2);
      Bond::BondDir bondDir = Bond::NONE;
      Bond::BondStereo bondStereo = Bond::STEREONONE;
      unsigned int cfg = 0;

      std::string temp =
          boost::algorithm::to_lower_copy(marvinBond->bondStereo.value);
      if (temp != "") {
        if (temp == "w") {
          bondDir = Bond::BEGINWEDGE;
        } else if (temp == "h") {
          bondDir = Bond::BEGINDASH;
        } else {
          std::ostringstream err;
          err << "unrecognized bond stereo value"
              << marvinBond->bondStereo.value;
          throw FileParseException(err.str());
        }
      } else if (marvinBond->bondStereo.convention != "") {
        if (marvinBond->bondStereo.convention != "MDL") {
          std::ostringstream err;
          err << "unrecognized bond stereo conventrion"
              << marvinBond->bondStereo.convention;
          throw FileParseException(err.str());
        }
        int mdlStereoVal = 0;
        if (!getCleanNumber(marvinBond->bondStereo.conventionValue,
                            mdlStereoVal)) {
          throw FileParseException(
              "MDL Convention Value must be one of: 1, 3, 4, 6");
        }
        switch (mdlStereoVal) {
          case 1:
            bondDir = Bond::BEGINWEDGE;
            break;
          case 6:
            bondDir = Bond::BEGINDASH;
            break;
          case 3:  // "either" double bond
            bondDir = Bond::EITHERDOUBLE;
            break;
          case 4:  // "either" single bond
            bondDir = Bond::UNKNOWN;
            break;
          default:
            throw FileParseException(
                "MDL Convention Value must be one of: 1, 3, 4, 6");
        }
      } else if (marvinBond->bondStereo.dictRef != "") {
        if (marvinBond->bondStereo.dictRef == "cml:W") {
          bondDir = Bond::BEGINWEDGE;
        } else if (marvinBond->bondStereo.dictRef == "cml:H") {
          bondDir = Bond::BEGINDASH;
        } else {
          throw FileParseException("dictRef must be one of: cml:W or cml:H");
        }
      } else {
        // nothing to do, no stereo information
      }

      switch (bondDir) {
        case Bond::BEGINWEDGE:
          bType = 1;  // overrides the type specification
          type = Bond::SINGLE;
          chiralityPossible = true;
          cfg = 1;
          break;
        case Bond::BEGINDASH:
          bType = 1;  // overrides the type specification
          type = Bond::SINGLE;
          cfg = 3;
          chiralityPossible = true;
          break;
        case Bond::EITHERDOUBLE:
          bondStereo = Bond::STEREOANY;
          bType = 2;  // overrides the type specification
          type = Bond::DOUBLE;
          cfg = 2;
          break;
        case Bond::UNKNOWN:
          bType = 1;  // overrides the type specification
          type = Bond::SINGLE;
          cfg = 2;
          break;
        default:
          // nothing to do
          break;
      }

      bond->setBondType(type);
      bond->setProp(common_properties::_MolFileBondType, bType);
      bond->setBondDir(bondDir);

      if (cfg != 0) {
        bond->setProp(common_properties::_MolFileBondCfg, cfg);
      }
      if (bondStereo != Bond::STEREONONE) {
        bond->setStereo(bondStereo);
      }

      // if we got an aromatic bond set the flag on the bond and the connected
      // atoms)
      if (bond->getBondType() == Bond::AROMATIC) {
        bond->setIsAromatic(true);
      }

      mol->addBond(bond, true);
    } catch (const std::exception &) {
      delete bond;
      throw;
    }
  }

  std::unique_ptr<RWMol> parseMolecule(MarvinMol *marvinMol,
                                       bool sanitize = false,
                                       bool removeHs = false) {
    PRECONDITION(marvinMol, "no molecule");
    std::vector<MarvinStereoGroup *> stereoGroups;
    std::unique_ptr<Conformer> confPtr;
    Conformer *conf = nullptr;
    std::unique_ptr<Conformer> conf3dPtr;
    Conformer *conf3d = nullptr;

    std::unique_ptr<RWMol> mol{new RWMol()};

    try {
      mol->setProp("_MolFileComments", "Generated by RDKit");

      // set the atoms

      if (marvinMol->hasAny2dCoords()) {
        conf = new Conformer(marvinMol->atoms.size());
        confPtr = std::unique_ptr<Conformer>(conf);
        confPtr->set3D(false);
      }
      if (marvinMol->hasAny3dCoords()) {
        conf3d = new Conformer(marvinMol->atoms.size());
        conf3dPtr = std::unique_ptr<Conformer>(conf3d);
        conf3d->set3D(true);
      }

      for (auto atomPtr : marvinMol->atoms) {
        Atom *atom = molAtomFromMarvinAtom(atomPtr, marvinMol);
        unsigned int aid = mol->addAtom(atom, false, true);

        if (conf != nullptr) {
          RDGeom::Point3D pos;
          if (atomPtr->x2 != DBL_MAX && atomPtr->y2 != DBL_MAX) {
            pos.x = atomPtr->x2;
            pos.y = atomPtr->y2;
          } else {
            pos.x = 0.0;
            pos.y = 0.0;
          }
          pos.z = 0.0;

          conf->setAtomPos(aid, pos);
        }

        if (conf3d != nullptr) {
          RDGeom::Point3D pos;
          if (atomPtr->x3 != DBL_MAX && atomPtr->y3 != DBL_MAX &&
              atomPtr->z3 != DBL_MAX) {
            pos.x = atomPtr->x3;
            pos.y = atomPtr->y3;
            pos.z = atomPtr->z3;
          } else {
            pos.x = 0.0;
            pos.y = 0.0;
            pos.z = 0.0;
          }

          conf3d->setAtomPos(aid, pos);
        }

        // also collect the stereo groups here

        if (atomPtr->mrvStereoGroup != "") {
          RDKit::StereoGroupType groupType;
          int groupNumber;

          // get the group parts
          std::string temp =
              boost::algorithm::to_lower_copy(atomPtr->mrvStereoGroup);
          if (temp == "abs") {
            groupType = RDKit::StereoGroupType::STEREO_ABSOLUTE;
            groupNumber = (-1);
          } else if (boost::starts_with(temp, "and")) {
            groupType = RDKit::StereoGroupType::STEREO_AND;
            if (!getCleanNumber(temp.substr(3), groupNumber)) {
              throw FileParseException(
                  "Group Number must be an integer in a stereo group AND# in a MRV file");
            }
          } else if (boost::starts_with(temp, "or")) {
            groupType = RDKit::StereoGroupType::STEREO_OR;
            if (!getCleanNumber(temp.substr(2), groupNumber)) {
              throw FileParseException(
                  "Group Number must be an integer in a stereo group OR# in a MRV file");
            }
          } else {
            throw FileParseException("Unrecognized group definition");
          }

          // see if the group already exists

          MarvinStereoGroup *marvinStereoGroup;
          auto groupIter =
              find_if(stereoGroups.begin(), stereoGroups.end(),
                      [groupType, groupNumber](const MarvinStereoGroup *arg) {
                        return arg->groupNumber == groupNumber &&
                               arg->groupType == groupType;
                      });
          if (groupIter != stereoGroups.end()) {
            marvinStereoGroup = *groupIter;
          } else {
            marvinStereoGroup = new MarvinStereoGroup(groupType, groupNumber);
            stereoGroups.push_back(marvinStereoGroup);
          }

          // add this atom to the group

          marvinStereoGroup->atoms.push_back(
              (unsigned int)marvinMol->getAtomIndex(atomPtr->id));
        }
      }

      if (conf != nullptr) {
        mol->addConformer(conf, true);
        confPtr.release();
      }

      if (conf3d != nullptr) {
        mol->addConformer(conf3d, true);
        conf3dPtr.release();
      }

      // set the bonds

      bool chiralityPossible = false;

      for (auto bondPtr : marvinMol->bonds) {
        molBondFromMarvinBond(bondPtr, marvinMol, mol.get(), chiralityPossible);
      }

      //  add the stereo groups

      std::vector<StereoGroup> groups;

      for (auto groupPtr : stereoGroups) {
        std::vector<Atom *> atoms;
        std::vector<Bond *> bonds;
        for (auto atomPtr : groupPtr->atoms) {
          atoms.push_back(mol->getAtomWithIdx(atomPtr));
        }

        groups.emplace_back(groupPtr->groupType, std::move(atoms),
                            std::move(bonds), groupPtr->groupNumber);
      }
      if (!groups.empty()) {
        mol->setStereoGroups(std::move(groups));
      }

      // // add the sgroup records

      int sequenceId = 0;
      // now the SuperatomSgroupsExpanded
      for (auto &marvinSgroup : marvinMol->sgroups) {
        auto sgroup = std::unique_ptr<SubstanceGroup>();

        marvinSgroup->parseMoleculeSpecific(mol.get(), sgroup, sequenceId);

        if (sgroup->getIsValid()) {
          addSubstanceGroup(*mol, *sgroup.get());
        }
        // increment the sequenceId

        sequenceId++;
      }

      mol->clearAllAtomBookmarks();
      mol->clearAllBondBookmarks();

      // calculate explicit valence on each atom:
      for (RWMol::AtomIterator atomIt = mol->beginAtoms();
           atomIt != mol->endAtoms(); ++atomIt) {
        (*atomIt)->calcExplicitValence(false);
        (*atomIt)->calcImplicitValence(false);
      }

      // update the chirality and stereo-chemistry
      //
      // NOTE: we detect the stereochemistry before sanitizing/removing
      // hydrogens because the removal of H atoms may actually remove
      // the wedged bond from the molecule.  This wipes out the only
      // sign that chirality ever existed and makes us sad... so first
      // perceive chirality, then remove the Hs and sanitize.
      //

      if (chiralityPossible && conf != nullptr) {
        DetectAtomStereoChemistry(*mol, conf);
      } else if (conf3d != nullptr) {
        mol->updatePropertyCache(false);
        MolOps::assignChiralTypesFrom3D(*mol, conf3d->getId(), true);
      }

      if (conf) {
        Atropisomers::detectAtropisomerChirality(*mol, conf);
      } else if (conf3d) {
        Atropisomers::detectAtropisomerChirality(*mol, conf3d);
      } else {
        Atropisomers::detectAtropisomerChirality(*mol, nullptr);
      }

      ClearSingleBondDirFlags(*mol);

      if (sanitize) {
        if (removeHs) {
          // Bond stereo detection must happen before H removal, or
          // else we might be removing stereogenic H atoms in double
          // bonds (e.g. imines). But before we run stereo detection,
          // we need to run mol cleanup so don't have trouble with
          // e.g. nitro groups. Sadly, this a;; means we will find
          // run both cleanup and ring finding twice (a fast find
          // rings in bond stereo detection, and another in
          // sanitization's SSSR symmetrization).
          unsigned int failedOp = 0;
          MolOps::sanitizeMol(*mol, failedOp, MolOps::SANITIZE_CLEANUP);
          MolOps::detectBondStereochemistry(*mol);
          MolOps::removeHs(*mol, false, false);
        } else {
          MolOps::sanitizeMol(*mol);
          MolOps::detectBondStereochemistry(*mol);
        }

        MolOps::assignStereochemistry(*mol, true, true, true);
      } else {
        MolOps::detectBondStereochemistry(*mol);
      }

      if (mol->hasProp(common_properties::_NeedsQueryScan)) {
        mol->clearProp(common_properties::_NeedsQueryScan);
        QueryOps::completeMolQueries(mol.get());
      }

      // clean up

      for (auto &stereoGroup : stereoGroups) {
        delete stereoGroup;
      }

      return mol;
    }

    catch (const std::exception &) {
      for (auto &stereoGroup : stereoGroups) {
        delete stereoGroup;
      }

      throw;
    }
  }

  MarvinMolBase *parseMarvinMolecule(
      boost::property_tree::ptree molTree,
      MarvinMol *parentMol = nullptr)  // parent is for sub-mols
  {
    MarvinMolBase *res = nullptr;

    try {
      std::string role = "";

      if (parentMol == nullptr) {
        res = new MarvinMol(molTree);
      } else  // is is a sub-mol - used for sGroups
      {
        role = molTree.get<std::string>("<xmlattr>.role", "");
        if (role == "") {
          throw FileParseException(
              "Expected a role for a sub-molecule in MRV file");
        }

        if (role == "SuperatomSgroup") {
          // there are two types of superatomSgroups - regular and expanded.
          // Expanded has a molecule attr of atomRefs

          std::string atomRefs =
              molTree.get<std::string>("<xmlattr>.atomRefs", "");
          auto atomArray = molTree.get_child_optional("atomArray");

          if (atomRefs == "" &&
              atomArray)  // no atomRefs means regular superatom
                          // - the atoms are in this superatom
          {
            res = new MarvinSuperatomSgroup(parentMol, molTree);

          } else  // if atomRefs and atomArray does not exist, this is an
                  // expanded superatom - the atoms are already in the parent
                  // mol
          {
            res = new MarvinSuperatomSgroupExpanded(parentMol, molTree);
          }
        } else if (role == "SruSgroup" || role == "CopolymerSgroup" ||
                   role == "ModificationSgroup") {
          res = new MarvinSruCoModSgroup(parentMol, role, molTree);
        } else if (role == "MultipleSgroup") {
          res = new MarvinMultipleSgroup(parentMol, molTree);
        } else if (role == "MulticenterSgroup") {
          res = new MarvinMulticenterSgroup(parentMol, molTree);
        } else if (role == "GenericSgroup") {
          res = new MarvinGenericSgroup(parentMol, molTree);
        } else if (role == "MonomerSgroup") {
          res = new MarvinMonomerSgroup(parentMol, molTree);
        } else if (role == "DataSgroup") {
          res = new MarvinDataSgroup(parentMol, molTree);
        } else {
          throw FileParseException("Unexpected role " + role + " in MRV file");
        }
      }

      for (auto &v : molTree) {
        if (v.first == "molecule") {
          MarvinMolBase *subMol =
              parseMarvinMolecule(v.second, (MarvinMol *)res);
          res->sgroups.push_back(std::unique_ptr<MarvinMolBase>(subMol));
        }
      }

      return res;
    } catch (const std::exception &) {
      delete res;

      throw;
    }
  };

  MarvinReaction *parseMarvinReaction(
      boost::property_tree::ptree rxnTree,
      boost::property_tree::ptree documentTree,
      bool parseArrowPlusesAndConditions =
          false)  // Arrows etc are not used in RDKIT) {
  {
    auto *res = new MarvinReaction();

    try {
      boost::property_tree::ptree childTree;
      bool foundChild = false;
      try {
        childTree = rxnTree.get_child("reactantList");
        foundChild = true;
      } catch (const std::exception &) {
        foundChild = false;
      }

      if (foundChild) {
        for (auto &v : childTree) {
          res->reactants.push_back(std::unique_ptr<MarvinMol>(
              (MarvinMol *)parseMarvinMolecule(v.second)));
        }
      }

      try {
        childTree = rxnTree.get_child("agentList");
        foundChild = true;
      } catch (const std::exception &) {
        foundChild = false;
      }
      if (foundChild) {
        for (auto &v : childTree) {
          res->agents.push_back(std::unique_ptr<MarvinMol>(
              (MarvinMol *)parseMarvinMolecule(v.second)));
        }
      }

      try {
        childTree = rxnTree.get_child("productList");
        foundChild = true;
      } catch (const std::exception &) {
        foundChild = false;
      }
      if (foundChild) {
        for (auto &v : childTree) {
          res->products.push_back(std::unique_ptr<MarvinMol>(
              (MarvinMol *)parseMarvinMolecule(v.second)));
        }
      }

      if (parseArrowPlusesAndConditions) {
        // <arrow type="DEFAULT" x1="-11.816189911577812"
        // y1="-10.001443743444021" x2="-8.401759471454618"
        // y2="-10.001443743444021"/>
        boost::property_tree::ptree arrow = rxnTree.get_child("arrow");
        res->arrow.type = arrow.get<std::string>("<xmlattr>.type", "");
        if (!getCleanNumber(arrow.get<std::string>("<xmlattr>.x1", ""),
                            res->arrow.x1) ||
            !getCleanNumber(arrow.get<std::string>("<xmlattr>.y1", ""),
                            res->arrow.y1) ||
            !getCleanNumber(arrow.get<std::string>("<xmlattr>.x2", ""),
                            res->arrow.x2) ||
            !getCleanNumber(arrow.get<std::string>("<xmlattr>.y2", ""),
                            res->arrow.y1)) {
          throw FileParseException(
              "Arrow coordinates must all be large floating point numbers in MRV file");
        }

        for (auto &v : documentTree) {
          if (v.first != "MReactionSign") {
            continue;
          }
          auto *marvinPlus = new MarvinPlus();
          res->pluses.push_back(std::unique_ptr<MarvinPlus>(marvinPlus));
          marvinPlus->id = v.second.get<std::string>("<xmlattr>.id", "");
          int pointCount = 0;
          for (auto &v2 : v.second) {
            if (v2.first == "MPoint") {
              double x;
              double y;
              if (!getCleanNumber(v2.second.get<std::string>("<xmlattr>.x", ""),
                                  x) ||
                  !getCleanNumber(v2.second.get<std::string>("<xmlattr>.y", ""),
                                  y)) {
                throw FileParseException(
                    "Plus sign  coordinates must all be large floating point numbers in MRV file");
              }

              switch (pointCount) {
                case 0:  //  first point - x1 and y1 are set
                  marvinPlus->x1 = x;
                  marvinPlus->y1 = y;
                  break;
                case 1:  // x2 is set, y1 is checked
                  marvinPlus->x2 = x;
                  if (marvinPlus->y1 != y) {
                    throw FileParseException(
                        "Plus sign coordinate Y in 2nd MPoint must be the same as that from the 1st MPoint in MRV file");
                  }
                  break;
                case 2:  // y2 is set, x2 is checked
                  marvinPlus->y2 = y;
                  if (marvinPlus->x2 != x) {
                    throw FileParseException(
                        "Plus sign coordinate X in 3rd MPoint must be the same as that from the 2nd MPoint in MRV file");
                  }
                  break;
                case 3:  // x2 and y2 are checked
                  if (marvinPlus->x1 != x) {
                    throw FileParseException(
                        "Plus sign coordinate X in 4th MPoint must be the same as that from the 1st MPoint in MRV file");
                  }

                  if (marvinPlus->y2 != y) {
                    throw FileParseException(
                        "Plus sign coordinate Y in 4th MPoint must be the same as that from the 3rd MPoint in MRV file");
                  }
                  break;

                default:
                  throw FileParseException(
                      "Plus sign coordinate must have 4 MPoints in MRV file");
              }
              ++pointCount;
            }
          }
        }

        for (auto &v : documentTree) {
          if (v.first != "MTextBox") {
            continue;
          }

          auto *marvinCondition = new MarvinCondition();
          res->conditions.push_back(
              std::unique_ptr<MarvinCondition>(marvinCondition));
          marvinCondition->id = v.second.get<std::string>("<xmlattr>.id", "");
          marvinCondition->halign =
              v.second.get<std::string>("<xmlattr>.halign", "");
          marvinCondition->valign =
              v.second.get<std::string>("<xmlattr>.valign", "");
          double fontScale;
          std::string fontScaleStr =
              v.second.get<std::string>("<xmlattr>.fontScale", "");
          if (fontScaleStr != "") {
            if (!getCleanNumber(fontScaleStr, fontScale)) {
              throw FileParseException(
                  "Condition font scale must be a positive integer in MRV file");
            }
          } else {
            fontScale = 0.0;
          }
          marvinCondition->fontScale = fontScale;

          marvinCondition->text = v.second.get<std::string>("Field", "");

          int pointCount = 0;
          for (auto &v2 : v.second) {
            if (v2.first == "MPoint") {
              double x, y;
              std::string xStr = v2.second.get<std::string>("<xmlattr>.x", "");
              std::string yStr = v2.second.get<std::string>("<xmlattr>.y", "");
              if (!getCleanNumber(xStr, x) || !getCleanNumber(yStr, y)) {
                throw FileParseException(
                    "Condition coordinate must valid integers in MRV file");
              }

              switch (pointCount) {
                case 0:  //  first point - x1 and y1 are set
                  marvinCondition->x = x;
                  marvinCondition->y = y;
                  break;
                case 1:
                case 2:
                case 3:  // x and Y are checked - must be the same as point 1
                  break;

                default:
                  throw FileParseException(
                      "Condition defs must have 4 MPoints in MRV file");
              }
              ++pointCount;
            }
          }
        }
      }  // end of if (parseArrowPlusesAndConditions)

      return res;
    } catch (const std::exception &) {
      delete res;

      throw;
    }
  }
};

bool MrvDataStreamIsReaction(std::istream &inStream) {
  PRECONDITION(inStream, "no stream");
  Utils::LocaleSwitcher ls;
  ptree tree;

  // Parse the XML into the property tree.

  read_xml(inStream, tree);

  // see if the reaction header is present
  try {
    auto rxn = tree.get_child("cml.MDocument.MChemicalStruct.reaction");
  } catch (const std::exception &) {
    return false;
  }

  return true;
}

bool MrvDataStreamIsReaction(std::istream *inStream) {
  PRECONDITION(inStream, "no stream");
  return MrvDataStreamIsReaction(*inStream);
}

bool MrvBlockIsReaction(const std::string &molmrvText) {
  std::istringstream inStream(molmrvText);
  return MrvDataStreamIsReaction(inStream);
}

bool MrvFileIsReaction(const std::string &fName) {
  std::ifstream inStream(fName.c_str());
  if (!inStream || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fName;
    throw BadFileException(errout.str());
  }
  if (!inStream.eof()) {
    return MrvDataStreamIsReaction(inStream);
  }
  return false;
}

//------------------------------------------------
//
//  Read a RWMol from a stream
//
//------------------------------------------------
std::unique_ptr<RWMol> MolFromMrvDataStream(std::istream &inStream,
                                            const MrvParserParams &params) {
  ptree tree;

  // Parse the XML into the property tree.

  read_xml(inStream, tree);

  MarvinCMLReader reader;
  return reader.parseMolecule(tree, params.sanitize, params.removeHs);
}
//------------------------------------------------
//
//  Read a RWMol from a string
//
//------------------------------------------------
std::unique_ptr<RWMol> MolFromMrvBlock(const std::string &molmrvText,
                                       const MrvParserParams &params) {
  std::istringstream inStream(molmrvText);
  return MolFromMrvDataStream(inStream, params);
}

//------------------------------------------------
//
//  Read an RWMol from a file
//
//------------------------------------------------
std::unique_ptr<RWMol> MolFromMrvFile(const std::string &fName,
                                      const MrvParserParams &params) {
  std::ifstream inStream(fName.c_str());
  if (!inStream || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fName;
    throw BadFileException(errout.str());
  }
  if (inStream.eof()) {
    return nullptr;
  }
  return MolFromMrvDataStream(inStream, params);
}

//------------------------------------------------
//
//  Read a ChemicalReaction from a stream
//
//------------------------------------------------
std::unique_ptr<ChemicalReaction> ReactionFromMrvDataStream(
    std::istream &inStream, const MrvParserParams &params) {
  PRECONDITION(inStream, "no stream");

  Utils::LocaleSwitcher ls;

  ptree tree;

  // Parse the XML into the property tree.

  read_xml(inStream, tree);

  MarvinCMLReader reader;
  return reader.parseReaction(tree, tree.get_child("cml.MDocument"),
                              params.sanitize, params.removeHs);
}

//------------------------------------------------
//
//  Read a ChemicalReaction from a string
//
//------------------------------------------------
std::unique_ptr<ChemicalReaction> ReactionFromMrvString(
    const std::string &molmrvText, const MrvParserParams &params) {
  std::istringstream inStream(molmrvText);
  return ReactionFromMrvDataStream(inStream, params);
}

//------------------------------------------------
//
//  Read a ChemicalReaction from a file
//
//------------------------------------------------
std::unique_ptr<ChemicalReaction> ReactionFromMrvFile(
    const std::string &fName, const MrvParserParams &params) {
  std::ifstream inStream(fName.c_str());
  if (!inStream || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fName;
    throw BadFileException(errout.str());
  }
  return ReactionFromMrvDataStream(inStream, params);
}
}  // namespace MarvinParser
}  // namespace v2
}  // namespace RDKit
