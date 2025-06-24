//
//  Copyright (C) 2022-2023 Tad Hurst, Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// this file WRITES MRV file for molecules and reactions

#include <string>
#include <exception>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cstdio>

#include <RDGeneral/LocaleSwitcher.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/GenericGroups/GenericGroups.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/Atropisomers.h>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSGroupWriting.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include "MarvinParser.h"
#include "MarvinDefs.h"

#include <GraphMol/RDKitQueries.h>
#include <GraphMol/StereoGroup.h>

#include <RDGeneral/StreamOps.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BoostEndInclude.h>

using namespace RDKit::SGroupWriting;

#define ARROW_MIN_LENGTH 1.0
#define ARROW_SPACE 0.5
#define PLUS_SPACE 1.0

namespace RDKit {
class MarvinCMLWriter {
  bool hasComplexQuery(const Atom *atom) {
    PRECONDITION(atom, "bad atom");
    bool res = false;
    if (atom->hasQuery()) {
      res = true;
      // counter examples:
      //  1) atomic number
      //  2) the smarts parser inserts AtomAnd queries
      //     for "C" or "c":
      //
      std::string descr = atom->getQuery()->getDescription();
      if (descr == "AtomAtomicNum") {
        res = false;
      } else if (descr == "AtomAnd") {
        if ((*atom->getQuery()->beginChildren())->getDescription() ==
            "AtomAtomicNum") {
          res = false;
        }
      }
    }
    return res;
  }

  void GetMarvinAtomInfo(const Atom *atom, MarvinAtom *marvinAtom) {
    PRECONDITION(atom, "bad atom");
    PRECONDITION(marvinAtom, "bad marvinAtom");

    if (atom->hasProp(common_properties::_MolFileRLabel)) {
      marvinAtom->elementType = "R";

      unsigned int rgroupRef;
      atom->getProp(common_properties::_MolFileRLabel, rgroupRef);
      marvinAtom->rgroupRef = (int)rgroupRef;

      std::string alias;
      if (atom->hasProp(common_properties::molFileAlias)) {
        atom->getProp(common_properties::molFileAlias, marvinAtom->mrvAlias);
      }
    } else if (atom->getAtomicNum()) {
      marvinAtom->elementType = atom->getSymbol();
    } else {
      if (!atom->hasProp(common_properties::dummyLabel)) {
        if (atom->hasQuery() &&
            (atom->getQuery()->getTypeLabel() == "A" ||
             atom->getQuery()->getTypeLabel() == "" ||
             (atom->getQuery()->getNegation() &&
              atom->getQuery()->getDescription() == "AtomAtomicNum" &&
              static_cast<ATOM_EQUALS_QUERY *>(atom->getQuery())->getVal() ==
                  1))) {
          marvinAtom->elementType = "*";
        } else if (atom->hasQuery() &&
                   (atom->getQuery()->getTypeLabel() == "Q" ||
                    (atom->getQuery()->getNegation() &&
                     atom->getQuery()->getDescription() == "AtomOr" &&
                     atom->getQuery()->endChildren() -
                             atom->getQuery()->beginChildren() ==
                         2 &&
                     (*atom->getQuery()->beginChildren())->getDescription() ==
                         "AtomAtomicNum" &&
                     static_cast<ATOM_EQUALS_QUERY *>(
                         (*atom->getQuery()->beginChildren()).get())
                             ->getVal() == 6 &&
                     (*++(atom->getQuery()->beginChildren()))
                             ->getDescription() == "AtomAtomicNum" &&
                     static_cast<ATOM_EQUALS_QUERY *>(
                         (*++(atom->getQuery()->beginChildren())).get())
                             ->getVal() == 1))) {
          throw MarvinWriterException(
              "Query atoms are not supported for MarvinWriter");
        } else if (atom->hasQuery() &&
                   atom->getQuery()->getTypeLabel() == "X") {
          throw MarvinWriterException(
              "Query atoms are not supported for MarvinWriter");
        } else if (atom->hasQuery() &&
                   atom->getQuery()->getTypeLabel() == "M") {
          throw MarvinWriterException(
              "Query atoms are not supported for MarvinWriter");
        } else if (atom->hasQuery() &&
                   atom->getQuery()->getTypeLabel() == "AH") {
          throw MarvinWriterException(
              "Query atoms are not supported for MarvinWriter");
        } else if (atom->hasQuery() &&
                   atom->getQuery()->getTypeLabel() == "QH") {
          throw MarvinWriterException(
              "Query atoms are not supported for MarvinWriter");
        } else if (atom->hasQuery() &&
                   atom->getQuery()->getTypeLabel() == "XH") {
          throw MarvinWriterException(
              "Query atoms are not supported for MarvinWriter");
        } else if (atom->hasQuery() &&
                   atom->getQuery()->getTypeLabel() == "MH") {
          throw MarvinWriterException(
              "Query atoms are not supported for MarvinWriter");
        } else if (hasComplexQuery(atom)) {
          throw MarvinWriterException(
              "Query atoms are not supported for MarvinWriter");
        } else {
          marvinAtom->elementType = "*";
        }
      } else {
        std::string symb;
        atom->getProp(common_properties::dummyLabel, symb);
        if (symb == "*") {
          marvinAtom->elementType = "*";
        } else if (symb == "X") {
          marvinAtom->elementType = "R";
          marvinAtom->rgroupRef = 0;
          marvinAtom->mrvAlias = "R";
        } else if (symb == "Xa") {
          marvinAtom->elementType = "R";
          marvinAtom->rgroupRef = 1;
          marvinAtom->mrvAlias = "R1";
        } else if (symb == "Xb") {
          marvinAtom->elementType = "R";
          marvinAtom->rgroupRef = 2;
          marvinAtom->mrvAlias = "R2";
        } else if (symb == "Xc") {
          marvinAtom->elementType = "R";
          marvinAtom->rgroupRef = 3;
          marvinAtom->mrvAlias = "R3";
        } else if (symb == "Xd") {
          marvinAtom->elementType = "R";
          marvinAtom->rgroupRef = 4;
          marvinAtom->mrvAlias = "R4";
        } else if (symb == "Xf") {
          marvinAtom->elementType = "R";
          marvinAtom->rgroupRef = 5;
          marvinAtom->mrvAlias = "R5";
        } else if (symb == "Xg") {
          marvinAtom->elementType = "R";
          marvinAtom->rgroupRef = 6;
          marvinAtom->mrvAlias = "R6";
        } else if (symb == "Xh") {
          marvinAtom->elementType = "R";
          marvinAtom->rgroupRef = 7;
          marvinAtom->mrvAlias = "R7";
        } else if (symb == "Xi") {
          marvinAtom->elementType = "R";
          marvinAtom->rgroupRef = 8;
          marvinAtom->mrvAlias = "R8";
        } else if (symb == "Xj") {
          marvinAtom->elementType = "R";
          marvinAtom->rgroupRef = 9;
          marvinAtom->mrvAlias = "R9";
        } else {
          throw MarvinWriterException(
              "Query atoms are not supported for MarvinWriter");
        }
      }
    }

    return;
  }

  bool isQueryBondInRing(const Bond *bond) {
    PRECONDITION(bond, "no bond");
    PRECONDITION(bond->hasQuery(), "no query");
    Bond::QUERYBOND_QUERY *qry = bond->getQuery();
    // start by catching combined bond order + bond topology queries

    if (qry->getDescription() == "BondAnd" && !qry->getNegation() &&
        qry->endChildren() - qry->beginChildren() == 2) {
      auto child1 = qry->beginChildren();
      auto child2 = child1 + 1;
      if (((*child1)->getDescription() == "BondInRing") !=
          ((*child2)->getDescription() == "BondInRing")) {
        if ((*child1)->getDescription() != "BondInRing") {
          qry = child2->get();
        } else {
          qry = child1->get();
        }
      }
    }
    if (qry->getDescription() == "BondInRing") {
      return true;
    }

    return false;
  }

  std::string getMarvinQueryBondSymbol(const Bond *bond) {
    PRECONDITION(bond, "no bond");
    PRECONDITION(bond->hasQuery(), "no query");

    Bond::QUERYBOND_QUERY *qry = bond->getQuery();
    if (qry->getDescription() == "BondOrder" || isQueryBondInRing(bond)) {
      return "";
    } else {
      // start by catching combined bond order + bond topology queries
      if (qry->getDescription() == "BondAnd" && !qry->getNegation() &&
          qry->endChildren() - qry->beginChildren() == 2) {
        auto child1 = qry->beginChildren();
        auto child2 = child1 + 1;
        if ((*child2)->getDescription() == "BondInRing") {
          qry = child1->get();
        } else if ((*child1)->getDescription() == "BondInRing") {
          qry = child2->get();
        }
      }
      if (qry->getDescription() == "BondOr" && !qry->getNegation()) {
        if (qry->endChildren() - qry->beginChildren() == 2) {
          auto child1 = qry->beginChildren();
          auto child2 = child1 + 1;
          if ((*child1)->getDescription() == "BondOrder" &&
              !(*child1)->getNegation() &&
              (*child2)->getDescription() == "BondOrder" &&
              !(*child2)->getNegation()) {
            // ok, it's a bond query we have a chance of dealing with
            int t1 = static_cast<BOND_EQUALS_QUERY *>(child1->get())->getVal();
            int t2 = static_cast<BOND_EQUALS_QUERY *>(child2->get())->getVal();
            if (t1 > t2) {
              std::swap(t1, t2);
            }
            if (t1 == Bond::SINGLE && t2 == Bond::DOUBLE) {
              return "SD";
            } else if (t1 == Bond::SINGLE && t2 == Bond::AROMATIC) {
              return "SA";
            } else if (t1 == Bond::DOUBLE && t2 == Bond::AROMATIC) {
              return "DA";
            }
          }
        }
      } else if (qry->getDescription() == "SingleOrAromaticBond" &&
                 !qry->getNegation()) {
        return "SA";
      } else if (qry->getDescription() == "SingleOrDoubleBond" &&
                 !qry->getNegation()) {
        return "SD";
      } else if (qry->getDescription() == "DoubleOrAromaticBond" &&
                 !qry->getNegation()) {
        return "DA";
      } else if (qry->getDescription() == "BondNull" && !qry->getNegation()) {
        return "Any";
      }
    }

    throw MarvinWriterException(
        "Only SA, DA, SD, and Any query bond are supported for MarvinWriter");
  }

  void GetMarvinBondSymbol(const Bond *bond, std::string &order,
                           std::string &queryType, std::string &convention) {
    PRECONDITION(bond, "");

    convention = "";
    order = "";
    queryType = "";

    if (bond->hasQuery()) {
      order = "1";
      queryType = getMarvinQueryBondSymbol(bond);
      if (queryType == "") {
        throw MarvinWriterException(
            "Only 1,2,3,Aromatic, and query bonds SA, DA, and SD are supported for MarvinWriter");
      }
      return;
    }

    queryType = "";  // not s query

    switch (bond->getBondType()) {
      case Bond::SINGLE:
        if (bond->getIsAromatic()) {
          order = "A";
        } else {
          order = "1";
        }
        break;

      case Bond::DOUBLE:
        if (bond->getIsAromatic()) {
          order = "A";
        } else {
          order = "2";
        }
        break;
      case Bond::TRIPLE:
        order = "3";
        break;

      case Bond::AROMATIC:
        order = "A";
        break;

      case Bond::DATIVE:
        convention = "cxn:coord";
        break;

      default:
        throw MarvinWriterException(
            "Only 1,2,3,Aromatic, and query bonds SA, DA, and SD are supported for MarvinWriter");
    }
  }

 private:
  bool hasNonDefaultValence(const Atom *atom) {
    PRECONDITION(atom, "no atom");
    if (atom->getNumRadicalElectrons() != 0) {
      return true;
    }

    if (atom->hasQuery()) {
      return false;
    }

    if (atom->getAtomicNum() == 1 ||
        SmilesWrite::inOrganicSubset(atom->getAtomicNum())) {
      return false;
    }

    return true;
  }

  MarvinMol *MolToMarvinMol(RWMol *mol, int &molCount, int &atomCount,
                            int &bondCount, int &sgCount, int confId = -1) {
    PRECONDITION(mol, "no molecule");
    // molCount is the starting and ending molCount - used when called from a
    // rxn

    MarvinMol *marvinMol = nullptr;
    const Conformer *conf = nullptr;
    const Conformer *conf3d = nullptr;
    int tempMolCount = 0, tempAtomCount = 0, tempBondCount = 0, tempSgCount = 0;
    try {
      marvinMol = new MarvinMol();

      marvinMol->molID = 'm' + std::to_string(++tempMolCount);

      // get a conformer

      int confCount = mol->getNumConformers();
      if (confCount > 0) {
        if (confId >= 0 && confId < confCount) {
          Conformer *testConf = &mol->getConformer(confId);
          if (!testConf->is3D()) {
            conf = testConf;
          } else {
            conf3d = testConf;
          }
        }

        if (conf == nullptr && conf3d == nullptr) {
          for (unsigned int confId = 0; confId < mol->getNumConformers();
               ++confId) {
            Conformer *testConf = &mol->getConformer(confId);
            if (!testConf->is3D()) {
              if (conf == nullptr) {  // only take the first 2d conf
                conf = testConf;
              }
            } else {
              if (conf3d == nullptr) {  // only take the first 3d conf
                conf3d = testConf;
              }
            }
            if (conf != nullptr && conf3d != nullptr) {
              break;
            }
          }
        }
      }

      if (mol->needsUpdatePropertyCache()) {
        mol->updatePropertyCache(false);
      }
      for (auto atom : mol->atoms()) {
        auto marvinAtom = new MarvinAtom();
        marvinMol->pushOwnedAtom(marvinAtom);

        marvinMol->atoms.push_back(marvinAtom);

        marvinAtom->id = 'a' + std::to_string(++tempAtomCount);

        GetMarvinAtomInfo(atom, marvinAtom);

        marvinAtom->formalCharge = atom->getFormalCharge();

        unsigned int nRadEs = atom->getNumRadicalElectrons();
        // value of radical electrons has to be one of the expected values or it
        // is ignored
        if (nRadEs != 0) {
          if (const auto iter = radicalElectronsToMarvinRadical.find(nRadEs);
              iter != radicalElectronsToMarvinRadical.end()) {
            marvinAtom->radical = iter->second;
          }
        }

        if (marvinAtom->isElement()) {
          marvinAtom->isotope = atom->getIsotope();

          if (marvinAtom->radical == "" && hasNonDefaultValence(atom)) {
            if (atom->getTotalDegree() == 0) {
              // Specify zero valence for elements/metals without neighbors
              // or hydrogens (degree 0) instead of writing them as radicals.
              marvinAtom->mrvValence = -1;
            } else {
              // if there are explicit Hs , mark them in the atom

              if (atom->getNoImplicit() && atom->getNumExplicitHs() > 0) {
                marvinAtom->hydrogenCount = atom->getNumExplicitHs();

              } else {
                unsigned int totalValence = atom->getTotalValence();
                if (totalValence != 15) {
                  marvinAtom->mrvValence = totalValence % 15;
                }
              }
            }
          }
        }

        if (!atom->getPropIfPresent(common_properties::molAtomMapNumber,
                                    marvinAtom->mrvMap)) {
          marvinAtom->mrvMap = 0;
        }

        if (conf != nullptr) {
          const RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
          marvinAtom->x2 = pos.x;
          marvinAtom->y2 = pos.y;
        } else {
          marvinAtom->x2 = DBL_MAX;
          marvinAtom->y2 = DBL_MAX;
        }

        if (conf3d != nullptr) {
          const RDGeom::Point3D pos = conf3d->getAtomPos(atom->getIdx());
          marvinAtom->x3 = pos.x;
          marvinAtom->y3 = pos.y;
          marvinAtom->z3 = pos.z;
        } else {
          marvinAtom->x3 = DBL_MAX;
          marvinAtom->y3 = DBL_MAX;
          marvinAtom->z3 = DBL_MAX;
        }

        //  atom maps for rxns
      }

      const Conformer *confToUse = nullptr;
      if (conf) {
        confToUse = conf;
      } else if (conf3d) {
        confToUse = conf3d;
      }
      auto wedgeBonds = Chirality::pickBondsToWedge(*mol, nullptr, confToUse);

      for (auto bond : mol->bonds()) {
        auto marvinBond = new MarvinBond();
        marvinMol->pushOwnedBond(marvinBond);
        marvinMol->bonds.push_back(marvinBond);

        marvinBond->id = 'b' + std::to_string(++tempBondCount);

        GetMarvinBondSymbol(bond, marvinBond->order, marvinBond->queryType,
                            marvinBond->convention);

        Bond::BondDir bondDirection = Bond::BondDir::NONE;
        bool reverse = false;

        if (conf) {
          Chirality::GetMolFileBondStereoInfo(bond, wedgeBonds, conf,
                                              bondDirection, reverse);
        } else if (conf3d) {
          Chirality::GetMolFileBondStereoInfo(bond, wedgeBonds, conf3d,
                                              bondDirection, reverse);
        }

        if (reverse) {
          // switch the begin and end atoms on the bond line
          marvinBond->atomRefs2[0] =
              marvinMol->atoms[bond->getEndAtomIdx()]->id;
          marvinBond->atomRefs2[1] =
              marvinMol->atoms[bond->getBeginAtomIdx()]->id;
        } else {
          marvinBond->atomRefs2[0] =
              marvinMol->atoms[bond->getBeginAtomIdx()]->id;
          marvinBond->atomRefs2[1] =
              marvinMol->atoms[bond->getEndAtomIdx()]->id;
        }

        switch (bondDirection) {
          case Bond::NONE:
            marvinBond->bondStereo.value = "";
            break;
          case Bond::BEGINWEDGE:
            marvinBond->bondStereo.value = "W";
            break;
          case Bond::BEGINDASH:
            marvinBond->bondStereo.value = "H";
            break;
          case Bond::UNKNOWN:
            marvinBond->bondStereo.value = "";
            marvinBond->bondStereo.convention = "MDL";
            marvinBond->bondStereo.conventionValue = "4";
            break;
          case Bond::EITHERDOUBLE:
            marvinBond->bondStereo.value = "";
            marvinBond->bondStereo.convention = "MDL";
            marvinBond->bondStereo.conventionValue = "3";
            break;

          default:
            marvinBond->bondStereo.value = "";  // other types are ignored
        }
      }

      // get all stereoGroups - add them to the correct atoms
      int orCount = 0;
      int andCount = 0;

      for (const StereoGroup &group : mol->getStereoGroups()) {
        std::string stereoGroupType;

        switch (group.getGroupType()) {
          case RDKit::StereoGroupType::STEREO_ABSOLUTE:
            stereoGroupType = "abs";
            break;
          case RDKit::StereoGroupType::STEREO_OR:
            stereoGroupType = "or" + std::to_string(++orCount);
            break;
          case RDKit::StereoGroupType::STEREO_AND:
            stereoGroupType = "and" + std::to_string(++andCount);
            break;
          default:
            throw MarvinWriterException("Unrecognized stereo group type");
        }

        std::vector<unsigned int> atomIds;
        Atropisomers::getAllAtomIdsForStereoGroup(*mol, group, atomIds,
                                                  wedgeBonds);

        for (auto atomId : atomIds) {
          marvinMol->atoms[atomId]->mrvStereoGroup = stereoGroupType;
        }
      }

      for (const SubstanceGroup &sgroup : getSubstanceGroups(*mol)) {
        std::string type;
        if (!sgroup.getPropIfPresent("TYPE", type)) {
          throw MarvinWriterException("TYPE not found for an Sgroup");
        }
        if (type == "SRU" || type == "MOD" || type == "COP") {
          std::string mrvType;
          if (type == "SRU") {
            mrvType = "SruSgroup";
          } else if (type == "MOD") {
            mrvType = "ModificationSgroup";
          } else if (type == "COP") {
            mrvType = "CopolymerSgroup";
          }

          auto marvinCoModSruSgroup =
              new MarvinSruCoModSgroup(mrvType, marvinMol);
          marvinMol->sgroups.push_back(
              std::unique_ptr<MarvinMolBase>(marvinCoModSruSgroup));

          if (!sgroup.getPropIfPresent("LABEL", marvinCoModSruSgroup->title)) {
            throw MarvinWriterException(
                "Expected a LABEL attribute for an SRU, MOD, or COP group");
          }

          if (!sgroup.getPropIfPresent("CONNECT",
                                       marvinCoModSruSgroup->connect)) {
            throw MarvinWriterException(
                "Expected a CONNECT attribute for an SRU, MOD, or COP group");
          }

          marvinCoModSruSgroup->id = "sg" + std::to_string(++tempSgCount);
          marvinCoModSruSgroup->molID = 'm' + std::to_string(++tempMolCount);

          for (auto atomIndex : sgroup.getAtoms()) {
            marvinCoModSruSgroup->atoms.push_back(marvinMol->atoms[atomIndex]);
            if (!marvinMol->atoms[atomIndex]->sGroupRefIsSuperatom) {
              marvinMol->atoms[atomIndex]->sgroupRef = marvinCoModSruSgroup->id;
            }
          }

          for (auto bondIndex : sgroup.getBonds()) {
            marvinCoModSruSgroup->bonds.push_back(marvinMol->bonds[bondIndex]);
          }
        }

        else if (type == "DAT") {
          auto marvinDataSgroup = new MarvinDataSgroup(marvinMol);
          marvinMol->sgroups.push_back(
              std::unique_ptr<MarvinMolBase>(marvinDataSgroup));

          marvinDataSgroup->id = "sg" + std::to_string(++tempSgCount);
          marvinDataSgroup->molID = 'm' + std::to_string(++tempMolCount);
          if (!sgroup.getPropIfPresent("FIELDNAME",
                                       marvinDataSgroup->fieldName)) {
            throw MarvinWriterException(
                "FIELDNAME not found for a SuperatomSgroup");
          }

          if (!sgroup.getPropIfPresent("QUERYTYPE",
                                       marvinDataSgroup->queryType)) {
            marvinDataSgroup->queryType = "";
          }
          if (!sgroup.getPropIfPresent("QUERYOP", marvinDataSgroup->queryOp)) {
            marvinDataSgroup->queryOp = "";
          }

          std::vector<std::string> fieldDatas;
          if (!sgroup.getPropIfPresent<std::vector<std::string>>("DATAFIELDS",
                                                                 fieldDatas)) {
            marvinDataSgroup->fieldData = "";
          } else {
            marvinDataSgroup->fieldData =
                boost::algorithm::join(fieldDatas, "\n");
          }

          if (!sgroup.getPropIfPresent("UNITS", marvinDataSgroup->units)) {
            marvinDataSgroup->units = "";
          }

          if (!sgroup.getPropIfPresent<double>("X", marvinDataSgroup->x)) {
            marvinDataSgroup->x = 0.0;
          }
          if (!sgroup.getPropIfPresent<double>("Y", marvinDataSgroup->y)) {
            marvinDataSgroup->y = 0.0;
          }
          if (!sgroup.getPropIfPresent("CONTEXT", marvinDataSgroup->context)) {
            marvinDataSgroup->context = "";
          }
          if (!sgroup.getPropIfPresent("PLACEMENT",
                                       marvinDataSgroup->placement)) {
            marvinDataSgroup->placement = "";
          }
          if (!sgroup.getPropIfPresent("UNITSDISPLAYED",
                                       marvinDataSgroup->unitsDisplayed)) {
            marvinDataSgroup->unitsDisplayed = "";
          }

          for (auto atomIndex : sgroup.getAtoms()) {
            marvinDataSgroup->atoms.push_back(marvinMol->atoms[atomIndex]);
            if (!marvinMol->atoms[atomIndex]->sGroupRefIsSuperatom) {
              marvinMol->atoms[atomIndex]->sgroupRef = marvinDataSgroup->id;
            }
          }
        }

        else if (type == "SUP") {
          auto superatomSgroupExpanded =
              new MarvinSuperatomSgroupExpanded(marvinMol);
          marvinMol->sgroups.push_back(
              std::unique_ptr<MarvinMolBase>(superatomSgroupExpanded));

          superatomSgroupExpanded->id = "sg" + std::to_string(++tempSgCount);
          superatomSgroupExpanded->molID = 'm' + std::to_string(++tempMolCount);

          if (!sgroup.getPropIfPresent("LABEL",
                                       superatomSgroupExpanded->title)) {
            throw MarvinWriterException(
                "LABEL not found for a SuperatomSgroup");
          }

          for (auto atomIndex : sgroup.getAtoms()) {
            superatomSgroupExpanded->atoms.push_back(
                marvinMol->atoms[atomIndex]);
          }
        }

        else if (type == "MUL") {
          auto marvinMultipleSgroup = new MarvinMultipleSgroup(marvinMol);
          marvinMol->sgroups.push_back(
              std::unique_ptr<MarvinMolBase>(marvinMultipleSgroup));
          marvinMultipleSgroup->id = "sg" + std::to_string(++tempSgCount);
          marvinMultipleSgroup->molID = 'm' + std::to_string(++tempMolCount);

          std::string titleValue;
          if (!sgroup.getPropIfPresent("MULT", titleValue) &&
              !sgroup.getPropIfPresent("LABEL", titleValue)) {
            throw MarvinWriterException("Title not found for a MultipleSgroup");
          }
          marvinMultipleSgroup->title = titleValue;

          for (auto atomIndex : sgroup.getAtoms()) {
            marvinMultipleSgroup->atoms.push_back(marvinMol->atoms[atomIndex]);
            if (!marvinMol->atoms[atomIndex]->sGroupRefIsSuperatom) {
              marvinMol->atoms[atomIndex]->sgroupRef = marvinMultipleSgroup->id;
            }
          }
          for (auto atomIndex : sgroup.getParentAtoms()) {
            marvinMultipleSgroup->parentAtoms.push_back(
                marvinMol->atoms[atomIndex]);
          }

          marvinMultipleSgroup->isExpanded = true;
        }

        else if (type == "GEN") {
          auto marvinGenericSgroup = new MarvinGenericSgroup(marvinMol);
          marvinMol->sgroups.push_back(
              std::unique_ptr<MarvinMolBase>(marvinGenericSgroup));
          marvinGenericSgroup->id = "sg" + std::to_string(++tempSgCount);
          marvinGenericSgroup->molID = 'm' + std::to_string(++tempMolCount);

          marvinGenericSgroup->charge =
              "onAtoms";  // RDKit has not place to put the charge location
                          // value, so we assume onAtoms here

          for (auto atomIndex : sgroup.getAtoms()) {
            MarvinAtom *marvinAtom = marvinMol->atoms[atomIndex];
            marvinGenericSgroup->atoms.push_back(marvinAtom);
            if (!marvinAtom->sGroupRefIsSuperatom) {
              marvinAtom->sgroupRef = marvinGenericSgroup->id;
            }
          }
        }

        else if (type == "MON") {
          // <molecule id="sg1" role="MonomerSgroup" title="mon"
          // charge="onAtoms" molID="m2" atomRefs="a2 a1 a3 a4">
          // </molecule>
          auto marvinMonomerSgroup = new MarvinMonomerSgroup(marvinMol);
          marvinMol->sgroups.push_back(
              std::unique_ptr<MarvinMolBase>(marvinMonomerSgroup));

          marvinMonomerSgroup->id = "sg" + std::to_string(++tempSgCount);
          marvinMonomerSgroup->molID = 'm' + std::to_string(++tempMolCount);

          std::string titleValue;
          if (!sgroup.getPropIfPresent("MULT", titleValue) &&
              !sgroup.getPropIfPresent("LABEL", titleValue)) {
            throw MarvinWriterException("Title not found for a MultipleSgroup");
          }
          marvinMonomerSgroup->title = titleValue;

          marvinMonomerSgroup->charge =
              "onAtoms";  // RDKit has not place to put the charge location
                          // value, so we assume onAtoms here
        }
      }

      // convert the superInfos to supergroups

      marvinMol->processSgroupsFromRDKit();
      std::map<std::string, std::string> sgMap;
      std::map<std::string, std::string> atomMap;
      std::map<std::string, std::string> bondMap;

      marvinMol->cleanUpNumbering(molCount, atomCount, bondCount, sgCount,
                                  sgMap, atomMap, bondMap);

      return marvinMol;
    } catch (const std::exception &e) {
      delete marvinMol;
      throw;
    }
  }

 public:
  MarvinMol *MolToMarvinMol(RWMol *mol, int confId = -1) {
    PRECONDITION(mol, "bad mol");

    int molCount = 0, atomCount = 0, bondCount = 0, sgCount = 0;

    return MolToMarvinMol(mol, molCount, atomCount, bondCount, sgCount, confId);
  }

  static bool compareRowsOfRectanglesReverse(std::vector<MarvinRectangle> &v1,
                                             std::vector<MarvinRectangle> &v2) {
    auto rect1 =
        MarvinRectangle(v1);  // composite rectangle for the previous row
    auto rect2 = MarvinRectangle(v2);  // composite rectangle for the row
    return MarvinRectangle::compareRectanglesByYReverse(
        rect1, rect2);  // just compare the first one in each row
  }

  double GetArrowPerdendicularPosition(
      std::vector<std::unique_ptr<MarvinMol>>
          &molList  // list of mols (agents) to examine
                    // for a space for the arrow
      ,
      bool verticalFlag)  // if verticalFlag, the arrow is to be placed
                          // horizonatally, so look for a vertical (y) space
  {
    // dividing the mols into rows and sorted by y value

    std::vector<MarvinRectangle> rectangleList;
    for (auto &mol : molList) {
      // see if there is horizontal overlap with any existing row

      MarvinRectangle molRect(mol->atoms);
      bool foundOverlap = false;
      for (MarvinRectangle rectangle : rectangleList) {
        if ((verticalFlag == true && molRect.overlapsVertically(rectangle)) ||
            (verticalFlag == false &&
             molRect.overlapsVHorizontally(rectangle))) {
          rectangle.extend(molRect);
          foundOverlap = true;
          break;
        }
      }

      if (!foundOverlap) {  // no overlap with a current row rectangle, so
                            // make a new one
        rectangleList.push_back(molRect);
      }
    }

    // sort the rows by X or Y, depending on vertical flag

    if (verticalFlag) {
      std::sort(rectangleList.begin(), rectangleList.end(),
                MarvinRectangle::compareRectanglesByYReverse);  // sort top down
    } else {
      std::sort(rectangleList.begin(), rectangleList.end(),
                MarvinRectangle::compareRectanglesByX);
    }

    // find a  spot for the arrow between rectangles, if possible

    for (auto rect1 = rectangleList.begin(); rect1 != rectangleList.end();
         ++rect1) {
      auto rect2 = rect1 + 1;
      if (rect2 == rectangleList.end()) {
        break;
      }

      // see if there is room between for the arrow

      if (verticalFlag) {
        if (rect2->upperLeft.y - rect1->lowerRight.y >= ARROW_SPACE) {
          return (rect2->upperLeft.y + rect1->lowerRight.y) / 2.0;
        }
      } else {
        if (rect2->upperLeft.x - rect1->lowerRight.x >= ARROW_SPACE) {
          return (rect2->upperLeft.x + rect1->lowerRight.x) / 2.0;
        }
      }
    }

    // if made it to here no spot was found, so place the arrow under the
    // bottom rectangle or left of the the leftmost one

    if (verticalFlag) {
      return rectangleList.front().lowerRight.y - ARROW_SPACE;
    } else {
      return rectangleList.front().upperLeft.x - ARROW_SPACE;
    }
  }

  void AddMarvinPluses(MarvinReaction &rxn,
                       std::vector<std::unique_ptr<MarvinMol>> &molList,
                       int &plusCount) {
    // dividing the mols into rows and sorted by y value

    std::vector<std::vector<MarvinRectangle>> rowsOfRectangles;
    for (auto &mol : molList) {
      if (!mol->hasCoords()) {
        return;  // no good way to make arrows
      }

      // see if there is horizontal overlap with any existing row

      MarvinRectangle molRect(mol->atoms);
      bool foundRow = false;
      for (std::vector<MarvinRectangle> &row : rowsOfRectangles) {
        for (MarvinRectangle rect : row) {
          if (molRect.overlapsVertically(rect)) {
            foundRow = true;
            row.push_back(molRect);
            break;
          }
        }
      }

      if (!foundRow)  // no overlap with a current row, so make a new one
      {
        std::vector<MarvinRectangle> newRow;
        newRow.push_back(molRect);
        rowsOfRectangles.push_back(newRow);
      }
    }

    // sort the members of each  row by X

    for (std::vector<MarvinRectangle> row : rowsOfRectangles) {
      std::sort(row.begin(), row.end(), MarvinRectangle::compareRectanglesByX);
    }

    // sort the rows by Y

    std::sort(rowsOfRectangles.begin(), rowsOfRectangles.end(),
              compareRowsOfRectanglesReverse);

    // make a plus between each rect on each row

    for (auto rowPtr = rowsOfRectangles.begin();
         rowPtr != rowsOfRectangles.end(); ++rowPtr) {
      for (auto rect1 = rowPtr->begin(); rect1 != rowPtr->end(); ++rect1) {
        auto rect2 = rect1 + 1;
        if (rect2 == rowPtr->end()) {
          break;
        }

        double x = (rect2->upperLeft.x + rect1->lowerRight.x) / 2.0;
        double y;

        // see if there is room between for the +
        if (rect2->lowerRight.x - rect1->upperLeft.x >= PLUS_SPACE) {
          y = (rect2->getCenter().y + rect1->getCenter().y) / 2.0;
        } else {  // put it under the two rectangles
          y = std::min<double>(rect1->lowerRight.y, rect2->lowerRight.y) -
              PLUS_SPACE / 2.0;
        }

        auto &newMarvinPlus = rxn.pluses.emplace_back(new MarvinPlus);

        newMarvinPlus->id = "o" + std::to_string(++plusCount);
        newMarvinPlus->x1 = x;
        newMarvinPlus->y1 = y;
        newMarvinPlus->x2 = x + PLUS_SPACE;
        newMarvinPlus->y2 = y + PLUS_SPACE;
      }

      // for each row after the first, add plus

      if (rowPtr != rowsOfRectangles.begin())  // 2nd and subsequent rows
      {
        // if these two rows only have enough space, put the  plus between
        // them

        auto rectPrev = MarvinRectangle(
            *(rowPtr - 1));  // composite rectangle for the previous row
        auto rectCurr =
            MarvinRectangle(*rowPtr);  // composite rectangle for the row

        auto &newMarvinPlus = rxn.pluses.emplace_back(new MarvinPlus);
        if (rectPrev.lowerRight.y - rectCurr.upperLeft.y > PLUS_SPACE) {
          double x = (rectPrev.getCenter().x + rectCurr.getCenter().x) / 2;
          double y = (rectPrev.lowerRight.y + rectCurr.upperLeft.y) / 2;
          newMarvinPlus->id = "o" + std::to_string(++plusCount);
          newMarvinPlus->x1 = x;
          newMarvinPlus->y1 = y;
          newMarvinPlus->x2 = x + (PLUS_SPACE);
          newMarvinPlus->y2 = y + (PLUS_SPACE);
        } else  // just put it in front of the current row
        {
          double x = rowPtr->front().upperLeft.x;
          double y = rowPtr->front().getCenter().y;
          newMarvinPlus->id = "o" + std::to_string(++plusCount);
          newMarvinPlus->x1 = x;
          newMarvinPlus->y1 = y;
          newMarvinPlus->x2 = x + (PLUS_SPACE);
          newMarvinPlus->y2 = y + (PLUS_SPACE);
        }
      }
    }
  }

  void SetArrow(MarvinReaction *marvinReaction) {
    PRECONDITION(marvinReaction, "bad reaction");

    // add a reaction arrow
    // get the overall rectangle for the reactants and the one for the
    // products

    // first set the bad (but parsable results) in case we cannot place the
    // arrow

    marvinReaction->arrow.x1 = 0.0;
    marvinReaction->arrow.x2 = 0.0;
    marvinReaction->arrow.y1 = 0.0;
    marvinReaction->arrow.y2 = 0.0;

    // make sure we have both reactants and products

    if (marvinReaction->reactants.size() == 0 ||
        marvinReaction->products.size() == 0) {
      return;
    }

    // make sure all atoms have coords

    for (auto &reactantPtr : marvinReaction->reactants) {
      if (!reactantPtr->hasCoords()) {
        return;
      }
    }
    for (auto &prodPtr : marvinReaction->products) {
      if (!prodPtr->hasCoords()) {
        return;
      }
    }

    marvinReaction->arrow.type = "DEFAULT";

    MarvinRectangle reactantRect(marvinReaction->reactants.front()->atoms);
    for (auto reactantPtr = marvinReaction->reactants.begin() + 1;
         reactantPtr != marvinReaction->reactants.end(); ++reactantPtr) {
      reactantRect.extend(MarvinRectangle((*reactantPtr)->atoms));
    }

    MarvinRectangle productRect(marvinReaction->products.front()->atoms);
    for (auto productPtr = marvinReaction->products.begin() + 1;
         productPtr != marvinReaction->products.end(); ++productPtr) {
      productRect.extend(MarvinRectangle((*productPtr)->atoms));
    }

    // if there is room between the reactants and products, put the arrow
    // there

    if (productRect.upperLeft.x - reactantRect.lowerRight.x >
        ARROW_MIN_LENGTH + 2.0 * ARROW_SPACE) {
      marvinReaction->arrow.x1 = reactantRect.lowerRight.x + ARROW_SPACE;
      marvinReaction->arrow.x2 = productRect.upperLeft.x - ARROW_SPACE;
      if (marvinReaction->agents.size() > 0) {
        marvinReaction->arrow.y1 =
            GetArrowPerdendicularPosition(marvinReaction->agents, true);
      } else {
        marvinReaction->arrow.y1 =
            (reactantRect.getCenter().y + productRect.getCenter().y) /
            2.0;  // no agents = just put it based on the reactant and
                  // products
      }
      marvinReaction->arrow.y2 = marvinReaction->arrow.y1;
    }
    // if not enough room between the reactants and product horizontally,
    // try vertically

    else if (reactantRect.lowerRight.y - productRect.upperLeft.y >
             ARROW_MIN_LENGTH + 2.0 * ARROW_SPACE) {
      if (marvinReaction->agents.size() > 0) {
        marvinReaction->arrow.x1 =
            GetArrowPerdendicularPosition(marvinReaction->agents, false);
      } else {
        marvinReaction->arrow.x1 =
            (reactantRect.getCenter().x + productRect.getCenter().x) /
            2.0;  // no agents = just put it based on the reactant and
                  // products
      }
      marvinReaction->arrow.x2 = marvinReaction->arrow.x1;
      marvinReaction->arrow.y1 = reactantRect.lowerRight.y - ARROW_SPACE;
      marvinReaction->arrow.y2 = productRect.upperLeft.y + ARROW_SPACE;
    }

    // if not good horizontal nor vertical place, just put it between the
    // centers (hack)

    else if ((reactantRect.getCenter() - productRect.getCenter()).length() >
             ARROW_MIN_LENGTH + 2.0 * ARROW_SPACE) {
      marvinReaction->arrow.x1 = reactantRect.getCenter().x;
      marvinReaction->arrow.x2 = productRect.getCenter().x;
      marvinReaction->arrow.y1 = reactantRect.getCenter().y;
      marvinReaction->arrow.y2 = productRect.getCenter().y;
    }

    // really no good place for the arrow

    else {
      marvinReaction->arrow.x1 = reactantRect.lowerRight.x + ARROW_SPACE;
      marvinReaction->arrow.x2 =
          reactantRect.lowerRight.x + ARROW_MIN_LENGTH + ARROW_SPACE;
      marvinReaction->arrow.y1 =
          (reactantRect.getCenter().y + productRect.getCenter().y) / 2.0;
      marvinReaction->arrow.y2 = marvinReaction->arrow.y1;
    }
  }

  MarvinReaction *ChemicalReactionToMarvinRxn(const ChemicalReaction *rxn,
                                              int confId = -1) {
    PRECONDITION(rxn, "bad reaction");

    MarvinReaction *marvinReaction = nullptr;
    try {
      auto marvinReaction = new MarvinReaction();
      int molCount = 0, atomCount = 0, bondCount = 0, sgCount = 0;
      for (const auto &mol : rxn->getReactants()) {
        RWMol rwMol(*mol);
        marvinReaction->reactants.emplace_back(MolToMarvinMol(
            &rwMol, molCount, atomCount, bondCount, sgCount, confId));
      }
      for (const auto &mol : rxn->getAgents()) {
        RWMol rwMol(*mol);
        marvinReaction->agents.emplace_back(MolToMarvinMol(
            &rwMol, molCount, atomCount, bondCount, sgCount, confId));
      }
      for (const auto &mol : rxn->getProducts()) {
        RWMol rwMol(*mol);
        marvinReaction->products.emplace_back(MolToMarvinMol(
            &rwMol, molCount, atomCount, bondCount, sgCount, confId));
      }

      // make up some pluses

      int plusCount = 0;
      AddMarvinPluses(*marvinReaction, marvinReaction->reactants, plusCount);
      AddMarvinPluses(*marvinReaction, marvinReaction->products, plusCount);

      // add a reaction arrow

      SetArrow(marvinReaction);

      return marvinReaction;
    } catch (const std::exception &e) {
      delete marvinReaction;
      throw;
    }
  }
};

std::string MolToMrvBlock(const ROMol &mol, const MrvWriterParams &params,
                          int confId) {
  Utils::LocaleSwitcher ls;

  RWMol trwmol(mol);
  // NOTE: kekulize the molecule before writing it out
  // because of the way mol files handle aromaticity
  if (trwmol.needsUpdatePropertyCache()) {
    trwmol.updatePropertyCache(false);
  }
  if (params.kekulize) {
    MolOps::Kekulize(trwmol);
  }

  if (params.includeStereo && !trwmol.getNumConformers()) {
    // generate coordinates so that the stereo we generate makes sense
    RDDepict::compute2DCoords(trwmol);
  }

  GenericGroups::convertGenericQueriesToSubstanceGroups(trwmol);

  MarvinCMLWriter marvinCMLWriter;

  auto marvinMol = marvinCMLWriter.MolToMarvinMol(&trwmol, confId);
  marvinMol->setPrecision(params.precision);
  ptree pt = marvinMol->toMolPtree();
  std::ostringstream out;
  if (params.prettyPrint) {
    write_xml(out, pt,
              boost::property_tree::xml_writer_make_settings<std::string>(
                  '\t', 1, "windows-1252"));
  } else {
    write_xml(out, pt);
  }
  std::string res = out.str();
  delete marvinMol;
  return res;
}

//------------------------------------------------
//
//  Dump a molecule to a file
//
//------------------------------------------------
void MolToMrvFile(const ROMol &mol, const std::string &fName,
                  const MrvWriterParams &params, int confId) {
  auto *outStream = new std::ofstream(fName.c_str());
  if (!(*outStream) || outStream->bad()) {
    delete outStream;
    std::ostringstream errout;
    errout << "Bad output file " << fName;
    throw BadFileException(errout.str());
  }
  std::string outString = MolToMrvBlock(mol, params, confId);
  *outStream << outString;
  delete outStream;
}

std::string ChemicalReactionToMrvBlock(const ChemicalReaction &rxn,
                                       bool prettyPrint) {
  Utils::LocaleSwitcher ls;

  MarvinCMLWriter marvinCMLWriter;

  auto marvinRxn = marvinCMLWriter.ChemicalReactionToMarvinRxn(&rxn);
  ptree pt = marvinRxn->toPtree();
  std::ostringstream out;
  if (prettyPrint) {
    write_xml(out, pt,
              boost::property_tree::xml_writer_make_settings<std::string>(
                  '\t', 1, "windows-1252"));
  } else {
    write_xml(out, pt);
  }
  delete marvinRxn;
  return out.str();
};

void ChemicalReactionToMrvFile(const ChemicalReaction &rxn,
                               const std::string &fName, bool prettyPrint) {
  auto *outStream = new std::ofstream(fName.c_str());
  if (!(*outStream) || outStream->bad()) {
    delete outStream;
    std::ostringstream errout;
    errout << "Bad output file " << fName;
    throw BadFileException(errout.str());
  }
  std::string outString = ChemicalReactionToMrvBlock(rxn, prettyPrint);
  *outStream << outString;
  delete outStream;
};
}  // namespace RDKit
