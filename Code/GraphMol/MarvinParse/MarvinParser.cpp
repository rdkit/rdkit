//
//  Copyright (C) 2002-2021 Greg Landrum and other RDKit contributors
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

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSGroupParsing.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include "MarvinParser.h"
#include "MarvinDefs.h"
#include <GraphMol/Conformer.h>

#include <GraphMol/RDKitQueries.h>
#include <GraphMol/StereoGroup.h>
#include <GraphMol/SubstanceGroup.h>

#include "MarvinParser.h"

#include <RDGeneral/StreamOps.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>

#ifdef RDKIT_USE_BOOST_REGEX
#include <boost/regex.hpp>
using boost::regex;
using boost::regex_match;
using boost::smatch;
// using boost::algorithm>;

#else
#include <regex>
using std::regex;
using std::regex_match;
using std::smatch;
#endif

using namespace RDKit::SGroupParsing;

namespace RDKit {

/*
      Imports the Marvin-specific dialect of CML (Chemical Markup Language) and
   converts it to datastructures that are compatible with Molfile, RXNfile, and
   Molfile complemented with canvas objects.
*/
class MarvinCMLReader {
 public:
  MarvinCMLReader(){};

  ~MarvinCMLReader(){};

  RWMol *parseMolecule(boost::property_tree::ptree molTree,
                       bool sanitize = false, bool removeHs = false) {
    RWMol *mol = nullptr;
    boost::property_tree::ptree molSection;

    try {
      molSection = molTree.get_child("cml.MDocument.MChemicalStruct.molecule");
    } catch (const std::exception &e) {
      try {
        molSection = molTree.get_child("cml.MDocument");
        return new RWMol();
      } catch (const std::exception &e) {
        try {
          molSection = molTree.get_child("cml");
          return new RWMol();
        } catch (const std::exception &e) {
          throw FileParseException("Expected \"molecule\" in MRV file");
        }
      }
    }

    try {
      std::unique_ptr<MarvinMol> marvinMol{
          (MarvinMol *)parseMarvinMolecule(molSection)};

      marvinMol->prepSgroupsForRDKit();

      mol = parseMolecule(marvinMol.get(), sanitize, removeHs);

       return mol;
    } catch (const std::exception &e) {
      delete mol;
      throw;
    }
  }

  ChemicalReaction *parseReaction(boost::property_tree::ptree rxnTree,
                                  boost::property_tree::ptree documentTree,
                                  bool sanitize = false,
                                  bool removeHs = false) {
    std::unique_ptr<ChemicalReaction> rxn{new ChemicalReaction()};
    rxnTree = rxnTree.get_child("cml.MDocument.MChemicalStruct.reaction");
    std::unique_ptr<MarvinReaction> marvinReaction{
        parseMarvinReaction(rxnTree, documentTree)};
    marvinReaction->prepSgroupsForRDKit();

    // get each reactant
    for (auto mol : marvinReaction->reactants) {
      rxn->addReactantTemplate(
          ROMOL_SPTR(parseMolecule(mol, sanitize, removeHs)));
    }

    // get each agent
    for (auto mol : marvinReaction->agents) {
      rxn->addAgentTemplate(ROMOL_SPTR(parseMolecule(mol, sanitize, removeHs)));
    }

    // get each product
    for (auto mol : marvinReaction->products) {
      rxn->addProductTemplate(
          ROMOL_SPTR(parseMolecule(mol, sanitize, removeHs)));
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

    return rxn.release();
  }

 private:
  bool getCleanDouble(std::string strToParse, double &outDouble) {
    if (boost::algorithm::trim_copy(strToParse) !=
        strToParse) {  // should be no white space
      return false;
    }
    try {
      outDouble = boost::lexical_cast<double>(strToParse);
    } catch (const std::exception &e) {
      return false;
    }

    return true;
  }

  bool getCleanInt(std::string strToParse, int &outInt) {
    if (boost::algorithm::trim_copy(strToParse) !=
        strToParse) {  // should be no white space
      return false;
    }
    try {
      outInt = boost::lexical_cast<int>(strToParse);
    } catch (const std::exception &e) {
      return false;
    }

    return true;
  }

  Atom *molAtomFromMarvinAtom(const MarvinAtom *marvinAtom,
                              const MarvinMolBase *marvinMolBase) {
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
        query->setIsotope(marvinAtom->rgroupRef);
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
        res->setProp(common_properties::_hasMassQuery, true);
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
    } catch (const std::exception &e) {
      delete res;
      res = nullptr;
      throw;
    }
  }

  void molBondFromMarvinBond(const MarvinBond *marvinBond,
                             const MarvinMol *marvinMol, RWMol *mol,
                             bool &chiralityPossible) {
    unsigned int bType;
    Bond *bond = nullptr;

    try {
      int idx1 = marvinMol->getAtomIndex(marvinBond->atomRefs2[0]);
      int idx2 = marvinMol->getAtomIndex(marvinBond->atomRefs2[1]);

      Bond::BondType type;

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
        if (!getCleanInt(marvinBond->bondStereo.conventionValue,
                         mdlStereoVal)) {
          throw FileParseException(
              "MDL Convention Value must be one of; 1, 3, 4, 6");
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
        }
      } else if (marvinBond->bondStereo.dictRef != "") {
        if (marvinBond->bondStereo.dictRef == "cml:W") {
          bondDir = Bond::BEGINWEDGE;
        } else if (marvinBond->bondStereo.dictRef == "cml:H") {
          bondDir = Bond::BEGINDASH;
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

      // v2k has no way to set stereoCare on bonds, so set the property if both
      // the beginning and end atoms have it set:
      int care1 = 0;
      int care2 = 0;
      if (!bond->hasProp(common_properties::molStereoCare) &&
          mol->getAtomWithIdx(bond->getBeginAtomIdx())
              ->getPropIfPresent(common_properties::molStereoCare, care1) &&
          mol->getAtomWithIdx(bond->getEndAtomIdx())
              ->getPropIfPresent(common_properties::molStereoCare, care2)) {
        if (care1 && care2) {
          bond->setProp(common_properties::molStereoCare, 1);
        }
      }
      mol->addBond(bond, true);
    } catch (const std::exception &e) {
      delete bond;
      throw;
    }
  }

  // see if the MarviMol has any coordinates

  bool hasCoords(MarvinMol *marvinMol) {
    for (auto atomPtr : marvinMol->atoms) {
      if (atomPtr->x2 != DBL_MAX && atomPtr->y2 != DBL_MAX) {
        return true;
      }
    }

    return false;
  }

  RWMol *parseMolecule(MarvinMol *marvinMol, bool sanitize = false,
                       bool removeHs = false) {
    std::vector<MarvinStereoGroup *> stereoGroups;
    Conformer *conf = nullptr;
    RWMol *mol = nullptr;

    try {
      mol = new RWMol();

      // mol->setProp(common_properties::_Name, marvinMol->molID);
      mol->setProp("_MolFileComments", "Generated by RDKit");

      // set the atoms

      if (hasCoords(marvinMol)) {
        conf = new Conformer(marvinMol->atoms.size());
        conf->set3D(false);
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
            if (!getCleanInt(temp.substr(3), groupNumber)) {
              throw FileParseException(
                  "Group Number must be an integer in a stereo group AND# in a MRV file");
            }
          } else if (boost::starts_with(temp, "or")) {
            groupType = RDKit::StereoGroupType::STEREO_OR;
            if (!getCleanInt(temp.substr(2), groupNumber)) {
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
        conf = nullptr;  // conf now owned by mol
      }

      // set the bonds

      bool chiralityPossible = false;

      for (auto bondPtr : marvinMol->bonds) {
        molBondFromMarvinBond(bondPtr, marvinMol, mol, chiralityPossible);
      }

      //  add the stereo groups

      std::vector<StereoGroup> groups;

      for (auto groupPtr : stereoGroups) {
        std::vector<Atom *> atoms;
        for (auto atomPtr : groupPtr->atoms) {
          atoms.push_back(mol->getAtomWithIdx(atomPtr));
        }

        groups.emplace_back(groupPtr->groupType, std::move(atoms));
      }
      if (!groups.empty()) {
        mol->setStereoGroups(std::move(groups));
      }

      // // add the sgroup records

      int sequenceId;

      // now the SuperatomSgroupsExpanded
      for (auto marvinSgroup : marvinMol->sgroups) {
        if (marvinSgroup->role() == "SuperatomSgroupExpanded") {
          auto subMolPtr = (MarvinSuperatomSgroupExpanded *)marvinSgroup;
          std::string typ = "SUP";
          SubstanceGroup sgroup = SubstanceGroup(mol, typ);
          sgroup.setProp<unsigned int>("index", sequenceId);

          for (auto atomPtr : subMolPtr->atoms) {
            sgroup.addAtomWithIdx(marvinMol->getAtomIndex(atomPtr->id));
          }

          for (auto bondPtr : subMolPtr->bonds) {
            sgroup.addBondWithIdx(marvinMol->getBondIndex(bondPtr->id));
          }

          sgroup.setProp("LABEL", subMolPtr->title);

          if (sgroup.getIsValid()) {
            addSubstanceGroup(*mol, sgroup);
          }
        }

        // now the MultipleSgroups
        // note: sequence continues counting from the loop above

        else if (marvinSgroup->role() == "MultipleSgroup") {
          auto subMolPtr = (MarvinMultipleSgroup *)marvinSgroup;

          std::string typ = "MUL";
          SubstanceGroup sgroup = SubstanceGroup(mol, typ);
          sgroup.setProp<unsigned int>("index", sequenceId);

          for (auto atomPtr : subMolPtr->atoms) {
            sgroup.addAtomWithIdx(marvinMol->getAtomIndex(atomPtr->id));
          }

          for (auto atomPtr : subMolPtr->parentAtoms) {
            sgroup.addParentAtomWithIdx(marvinMol->getAtomIndex(atomPtr->id));
          }

          // the connection bonds

          for (auto bondptr : subMolPtr->bondsToAtomsNotInExpandedGroup) {
            sgroup.addBondWithIdx(marvinMol->getBondIndex(bondptr->id));
          }

          sgroup.setProp("MULT", subMolPtr->title);

          if (sgroup.getIsValid()) {
            addSubstanceGroup(*mol, sgroup);
          }
        }

        else if (marvinSgroup->role() == "SruSgroup" ||
                 marvinSgroup->role() == "CopolymerSgroup" ||
                 marvinSgroup->role() == "ModificationSgroup") {
          auto subMolPtr = (MarvinSruCoModSgroup *)marvinSgroup;

          std::string typ;
          if (subMolPtr->role() == "SruSgroup") {
            typ = "SRU";
          } else if (subMolPtr->role() == "CopolymerSgroup") {
            typ = "COP";
          } else if (subMolPtr->role() == "ModificationSgroup") {
            typ = "MOD";
          } else {
            throw FileParseException(
                "Internal error: unrecognized role in a MarvinSruCoPolSgroup");
          }

          SubstanceGroup sgroup = SubstanceGroup(mol, typ);
          sgroup.setProp<unsigned int>("index", sequenceId);

          sgroup.setProp("CONNECT", subMolPtr->connect);

          for (auto atomPtr : subMolPtr->atoms) {
            sgroup.addAtomWithIdx(marvinMol->getAtomIndex(atomPtr->id));
          }

          for (auto bondPtr : subMolPtr->bonds) {
            sgroup.addBondWithIdx(marvinMol->getBondIndex(bondPtr->id));
          }

          sgroup.setProp("LABEL", subMolPtr->title);

          if (sgroup.getIsValid()) {
            addSubstanceGroup(*mol, sgroup);
          }
        }

        // Now the data groups

        else if (marvinSgroup->role() == "DataSgroup") {
          auto subMolPtr = (MarvinDataSgroup *)marvinSgroup;

          std::string typ = "DAT";
          SubstanceGroup sgroup = SubstanceGroup(mol, typ);
          sgroup.setProp<unsigned int>("index", sequenceId);

          for (auto atomPtr : subMolPtr->atoms) {
            sgroup.addAtomWithIdx(marvinMol->getAtomIndex(atomPtr->id));
          }

          sgroup.setProp("FIELDNAME", subMolPtr->fieldName);

          if (subMolPtr->queryType != "") {
            sgroup.setProp("QUERYTYPE", subMolPtr->queryType);
          }
          if (subMolPtr->queryOp != "") {
            sgroup.setProp("QUERYOP", subMolPtr->queryOp);
          }

          std::ostringstream out;
          out << std::fixed << std::setw(10) << std::setprecision(4)
              << subMolPtr->x << std::fixed << std::setw(10)
              << std::setprecision(4) << subMolPtr->y
              << "    DRU   ALL  0       0";

          sgroup.setProp("FIELDDISP", out.str());  // really not used by RDKIT

          std::vector<std::string> fieldDatas;
          fieldDatas.push_back(subMolPtr->fieldData);
          sgroup.setProp("DATAFIELDS", fieldDatas);

          // The following props are not part of the RDKit structure for MOL
          // files, but we save them so that we can round-trip the MRV

          sgroup.setProp("UNITS", subMolPtr->units);
          sgroup.setProp("UNITSDISPLAYED", subMolPtr->unitsDisplayed);
          sgroup.setProp("CONTEXT", subMolPtr->context);
          sgroup.setProp("PLACEMENT", subMolPtr->placement);
          sgroup.setProp("X", subMolPtr->x);
          sgroup.setProp("Y", subMolPtr->y);

          if (sgroup.getIsValid()) {
            addSubstanceGroup(*mol, sgroup);
          }
        }

        // the MultiCenter Sgroups

        else if (marvinSgroup->role() == "MulticenterSgroup") {
          // There is really no place to put these in RDKit.  We should have
          // removed these already

          throw FileParseException(
              "Internal error: a MarvinMulticenterSgroup had not been removed");
        }

        // the Generic groups

        else if (marvinSgroup->role() == "GenericSgroup") {
          auto subMolPtr = (MarvinGenericSgroup *)marvinSgroup;

          std::string typ = "GEN";
          SubstanceGroup sgroup = SubstanceGroup(mol, typ);
          sgroup.setProp<unsigned int>("index", sequenceId);

          for (auto atomPtr : subMolPtr->atoms) {
            sgroup.addAtomWithIdx(marvinMol->getAtomIndex(atomPtr->id));
          }

          // note: there is no place to put the change="onAtoms" flag

          if (sgroup.getIsValid()) {
            addSubstanceGroup(*mol, sgroup);
          }
        }

        // now the MonomerSgroups
        // note: sequence continues counting from the loop above

        else if (marvinSgroup->role() == "MonomerSgroup") {
          auto subMolPtr = (MarvinMonomerSgroup *)marvinSgroup;

          std::string typ = "MON";
          SubstanceGroup sgroup = SubstanceGroup(mol, typ);
          sgroup.setProp<unsigned int>("index", sequenceId);

          std::vector<MarvinAtom *>::const_iterator atomIter;
          for (auto atomPtr : subMolPtr->atoms) {
            int atomIndex = marvinMol->getAtomIndex(atomPtr->id);
            sgroup.addAtomWithIdx(atomIndex);
          }
          sgroup.setProp("LABEL", subMolPtr->title);

          // Note: RDKit does not have a place for the Bracket information nor
          // the charge="onAtoms" attr

          if (sgroup.getIsValid()) {
            addSubstanceGroup(*mol, sgroup);
          }
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

      if (mol->getNumConformers() > 0) {
        const Conformer &conf2 = mol->getConformer();
        if (chiralityPossible) {
          DetectAtomStereoChemistry(*mol, &conf2);
        }
      }

      if (sanitize) {
        if (removeHs) {
          MolOps::removeHs(*mol, false, false);
        } else {
          MolOps::sanitizeMol(*mol);
        }

        // now that atom stereochem has been perceived, the wedging
        // information is no longer needed, so we clear
        // single bond dir flags:
        ClearSingleBondDirFlags(*mol);

        // unlike DetectAtomStereoChemistry we call detectBondStereochemistry
        // here after sanitization because we need the ring information:
        MolOps::detectBondStereochemistry(*mol);

        MolOps::assignStereochemistry(*mol, true, true, true);
      } else {
        // we still need to do something about double bond stereochemistry
        // (was github issue 337)
        // now that atom stereochem has been perceived, the wedging
        // information is no longer needed, so we clear
        // single bond dir flags:
        ClearSingleBondDirFlags(*mol);
        MolOps::detectBondStereochemistry(*mol);
      }

      if (mol->hasProp(common_properties::_NeedsQueryScan)) {
        mol->clearProp(common_properties::_NeedsQueryScan);
        QueryOps::completeMolQueries(mol);
      }

      // clean up

      for (auto &stereoGroup : stereoGroups) {
        delete stereoGroup;
      }

      return mol;
    }

    catch (const std::exception &e) {
      delete mol;

      delete conf;

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
      std::string molID = molTree.get<std::string>("<xmlattr>.molID", "");

      if (molID == "") {
        throw FileParseException("Expected a molID in MRV file");
      }

      std::string role = "";

      if (parentMol == nullptr) {
        res = new MarvinMol();
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
            auto *marvinSuperatomSgroup = new MarvinSuperatomSgroup(parentMol);
            res = marvinSuperatomSgroup;

            marvinSuperatomSgroup->id =
                molTree.get<std::string>("<xmlattr>.id", "");

            marvinSuperatomSgroup->title =
                molTree.get<std::string>("<xmlattr>.title", "");
            if (marvinSuperatomSgroup->title == "") {
              throw FileParseException(
                  "Expected  title for a SuperatomSgroup definition in MRV file");
            }
          } else  // if atomRefs and atomArray does not exist, this is an
                  // expanded superatom - the atoms are already in the parent
                  // mol
          {
            auto *marvinSuperatomSgroupExpanded =
                new MarvinSuperatomSgroupExpanded(parentMol);
            res = marvinSuperatomSgroupExpanded;
            role = "SuperatomSgroupExpanded";

            marvinSuperatomSgroupExpanded->id =
                molTree.get<std::string>("<xmlattr>.id", "");

            if (atomRefs != "") {
              std::vector<std::string> atomList;
              boost::algorithm::split(atomList, atomRefs,
                                      boost::algorithm::is_space());
              for (auto &it : atomList) {
                auto atomPtr = res->parent->findAtomByRef(it);

                if (atomPtr == nullptr) {
                  throw FileParseException(
                      "AtomRef specification for an SuperatomSgroupExpanded group definition was not found in any parent");
                }
                marvinSuperatomSgroupExpanded->atoms.push_back(atomPtr);
              }
            } else  // must have no atomRefs nor AtomAray - get the atoms from
                    // the parent block - the ones that referene this
                    // superSgroup
            {
              for (MarvinAtom *atom : parentMol->atoms) {
                if (atom->sgroupRef == marvinSuperatomSgroupExpanded->id) {
                  marvinSuperatomSgroupExpanded->atoms.push_back(atom);
                }
              }
            }

            marvinSuperatomSgroupExpanded->title =
                molTree.get<std::string>("<xmlattr>.title", "");
            if (marvinSuperatomSgroupExpanded->title == "") {
              throw FileParseException(
                  "Expected  title for a SuperatomSgroupExpanded definition in MRV file");
            }
          }
        } else if (role == "SruSgroup" || role == "CopolymerSgroup" ||
                   role == "ModificationSgroup") {
          //      <molecule molID="m2" id="sg1" role="SruSgroup" atomRefs="a1 a3
          //      a4" title="n" connect="ht" correspondence="" bondList=""/>
          //      <molecule molID="m3" id="sg2" role="SruSgroup" atomRefs="a5 a6
          //      a7 a8" title="n" connect="hh" correspondence="" bondList=""/>
          //      <molecule molID="m4" id="sg3" role="SruSgroup" atomRefs="a10
          //      a11" title="n" connect="eu" correspondence=""
          //      bondList=""/></molecule>

          auto *marvinCoModSruSgroup =
              new MarvinSruCoModSgroup(role, parentMol);
          res = marvinCoModSruSgroup;

          marvinCoModSruSgroup->id =
              molTree.get<std::string>("<xmlattr>.id", "");

          std::string atomRefsStr =
              molTree.get<std::string>("<xmlattr>.atomRefs", "");
          if (atomRefsStr == "") {
            throw FileParseException(
                "Expected  atomRefs for a SruSgroup definition in MRV file");
          }

          std::vector<std::string> atomList;
          boost::algorithm::split(atomList, atomRefsStr,
                                  boost::algorithm::is_space());
          for (auto &it : atomList) {
            auto atomPtr = res->parent->findAtomByRef(it);

            if (atomPtr == nullptr) {
              throw FileParseException(
                  "AtomRef specification for an SRU, Coplomer, or Modification group definition was not found in any parent");
            }
            marvinCoModSruSgroup->atoms.push_back(atomPtr);
          }

          marvinCoModSruSgroup->title =
              molTree.get<std::string>("<xmlattr>.title", "");
          if (marvinCoModSruSgroup->title == "") {
            throw FileParseException(
                "Expected title for a SruSgroup definition in MRV file");
          }

          // the title for an SRU group must be a lower case letter, or a range
          // of two positive ints (4-6)

          marvinCoModSruSgroup->connect =
              molTree.get<std::string>("<xmlattr>.connect", "");
          if (marvinCoModSruSgroup->connect == "") {
            marvinCoModSruSgroup->connect = "ht";
          }
          if (!boost::algorithm::contains(
                  sruSgroupConnectChoices,
                  std::vector<std::string>{marvinCoModSruSgroup->connect})) {
            std::ostringstream err;
            err << "Expected a connect  string of \""
                << boost::algorithm::join(sruSgroupConnectChoices, ", ")
                << "\" for a SruSgroup definition in MRV file";
            throw FileParseException(err.str());
          }

          marvinCoModSruSgroup->correspondence =
              molTree.get<std::string>("<xmlattr>.correspondence", "");

          std::string bondListStr =
              molTree.get<std::string>("<xmlattr>.bondList", "");
          if (bondListStr != "") {
            std::vector<std::string> bondList;
            boost::algorithm::split(bondList, bondListStr,
                                    boost::algorithm::is_space());

            for (auto &it : bondList) {
              auto bondPtr = res->parent->findBondByRef(it);

              if (bondPtr == nullptr) {
                throw FileParseException(
                    "BondList specification for an SRU, Coplomer, or Modification group definition was not found in any parent");
              }
              marvinCoModSruSgroup->bonds.push_back(bondPtr);
            }
          }
        } else if (role == "MultipleSgroup") {
          // first make sure this can be made into a proper MultipleSgroup.
          // To be valid, the title must be a positive integer and there must be
          // exactly two bonds the connect to the repeating group

          auto *marvinMultipleSgroup = new MarvinMultipleSgroup(parentMol);
          res = marvinMultipleSgroup;

          marvinMultipleSgroup->id =
              molTree.get<std::string>("<xmlattr>.id", "");

          std::string atomRefsStr =
              molTree.get<std::string>("<xmlattr>.atomRefs", "");
          if (atomRefsStr == "") {
            throw FileParseException(
                "Expected  atomRefs for a MultipleSgroup definition in MRV file");
          }

          std::vector<std::string> atomList;
          boost::algorithm::split(atomList, atomRefsStr,
                                  boost::algorithm::is_space());
          for (auto &it : atomList) {
            auto atomPtr = res->parent->findAtomByRef(it);

            if (atomPtr == nullptr) {
              throw FileParseException(
                  "AtomRef specification for a MultipleSgroup definition was not found in any parent");
            }

            marvinMultipleSgroup->atoms.push_back(atomPtr);
            marvinMultipleSgroup->parentAtoms.push_back(atomPtr);
          }

          marvinMultipleSgroup->title =
              molTree.get<std::string>("<xmlattr>.title", "");
          int testInt;
          if (marvinMultipleSgroup->title == "" ||
              !getCleanInt(marvinMultipleSgroup->title, testInt) ||
              testInt <= 0) {
            throw FileParseException(
                "Expected a positive integer title for a MultipleSgroup definition in MRV file");
          }
        } else if (role == "DataSgroup") {
          auto *marvinDataSgroup = new MarvinDataSgroup(parentMol);
          res = marvinDataSgroup;

          marvinDataSgroup->id = molTree.get<std::string>("<xmlattr>.id", "");

          std::string atomRefsStr =
              molTree.get<std::string>("<xmlattr>.atomRefs", "");
          if (atomRefsStr == "") {
            throw FileParseException(
                "Expected  atomRefs for a DataSgroup definition in MRV file");
          }

          std::vector<std::string> atomList;
          boost::algorithm::split(atomList, atomRefsStr,
                                  boost::algorithm::is_space());
          for (auto &it : atomList) {
            auto atomPtr = res->parent->findAtomByRef(it);

            if (atomPtr == nullptr) {
              throw FileParseException(
                  "AtomRef specification for a DataSgroup definition was not found in any parent");
            }

            marvinDataSgroup->atoms.push_back(atomPtr);
          }

          marvinDataSgroup->context =
              molTree.get<std::string>("<xmlattr>.context", "");

          marvinDataSgroup->fieldName =
              molTree.get<std::string>("<xmlattr>.fieldName", "");
          marvinDataSgroup->placement =
              molTree.get<std::string>("<xmlattr>.placement", "");

          marvinDataSgroup->unitsDisplayed = molTree.get<std::string>(
              "<xmlattr>.unitsDisplayed", "Unit not displayed");
          if (marvinDataSgroup->unitsDisplayed == "") {
            marvinDataSgroup->unitsDisplayed = "Unit not displayed";
          }
          std::string unitsDisplayed =
              boost::algorithm::to_lower_copy(marvinDataSgroup->unitsDisplayed);

          if (unitsDisplayed != "unit displayed" &&
              unitsDisplayed != "unit not displayed") {
            throw FileParseException(
                "Expected unitsDisplayed to be either \"Unit displayed\" or \"Unit not displayed\" for a DataSgroup definition in MRV file");
          }

          marvinDataSgroup->queryType =
              molTree.get<std::string>("<xmlattr>.queryType", "");
          marvinDataSgroup->queryOp =
              molTree.get<std::string>("<xmlattr>.queryOp", "");
          marvinDataSgroup->units =
              molTree.get<std::string>("<xmlattr>.units", "");

          marvinDataSgroup->fieldData =
              molTree.get<std::string>("<xmlattr>.fieldData", "");

          std::string x = molTree.get<std::string>("<xmlattr>.x", "0.0");
          std::string y = molTree.get<std::string>("<xmlattr>.y", "0.0");

          if (x == "") {
            throw FileParseException(
                "Expected x for a DataSgroup definition in MRV file");
          }
          if (!getCleanDouble(x, marvinDataSgroup->x)) {
            throw FileParseException(
                "The value for x must be a floating point value in MRV file");
          }
          if (y == "") {
            throw FileParseException(
                "Expected y for a DataSgroup definition in MRV file");
          }
          if (!getCleanDouble(y, marvinDataSgroup->y)) {
            throw FileParseException(
                "The value for y must be a floating point value in MRV file");
          }
        } else if (role == "MulticenterSgroup") {
          auto *marvinMulticenterSgroup =
              new MarvinMulticenterSgroup(parentMol);
          res = marvinMulticenterSgroup;

          marvinMulticenterSgroup->id =
              molTree.get<std::string>("<xmlattr>.id", "");

          std::string atomRefsStr =
              molTree.get<std::string>("<xmlattr>.atomRefs", "");
          if (atomRefsStr == "") {
            throw FileParseException(
                "Expected  atomRefs for a MulticenterSgroup definition in MRV file");
          }

          std::vector<std::string> atomList;
          boost::algorithm::split(atomList, atomRefsStr,
                                  boost::algorithm::is_space());
          for (auto &it : atomList) {
            auto atomPtr = res->parent->findAtomByRef(it);

            if (atomPtr == nullptr) {
              throw FileParseException(
                  "AtomRef specification for a MulticenterSgroup definition was not found in any parent");
            }

            marvinMulticenterSgroup->atoms.push_back(atomPtr);
          }

          std::string centerId =
              molTree.get<std::string>("<xmlattr>.center", "");

          auto atomPtr = res->parent->findAtomByRef(centerId);

          if (atomPtr == nullptr) {
            throw FileParseException(
                "Center specification for a MulticenterSgroup definition was not found in any parent");
          }

          marvinMulticenterSgroup->center = atomPtr;
        } else if (role == "GenericSgroup") {
          auto marvinGenericSgroup = new MarvinGenericSgroup(parentMol);
          res = marvinGenericSgroup;

          marvinGenericSgroup->id =
              molTree.get<std::string>("<xmlattr>.id", "");

          std::string atomRefsStr =
              molTree.get<std::string>("<xmlattr>.atomRefs", "");
          if (atomRefsStr == "") {
            throw FileParseException(
                "Expected  atomRefs for a GenericSgroup definition in MRV file");
          }

          std::vector<std::string> atomList;
          boost::algorithm::split(atomList, atomRefsStr,
                                  boost::algorithm::is_space());
          for (auto &it : atomList) {
            auto atomPtr = res->parent->findAtomByRef(it);

            if (atomPtr == nullptr) {
              throw FileParseException(
                  "AtomRef specification for a GenericSgroup definition was not found in any parent");
            }

            marvinGenericSgroup->atoms.push_back(atomPtr);
          }

          marvinGenericSgroup->charge =
              molTree.get<std::string>("<xmlattr>.charge", "");
          if (marvinGenericSgroup->charge != "onAtoms" &&
              marvinGenericSgroup->charge != "onBracket") {
            throw FileParseException(
                "Expected  omAtoms or onBracket for a charge attr of a GenericSgroup definition in MRV file");
          }
        } else if (role == "MonomerSgroup") {
          auto *marvinMonomerSgroup = new MarvinMonomerSgroup(parentMol);
          res = marvinMonomerSgroup;

          marvinMonomerSgroup->id =
              molTree.get<std::string>("<xmlattr>.id", "");

          std::string atomRefsStr =
              molTree.get<std::string>("<xmlattr>.atomRefs", "");
          if (atomRefsStr == "") {
            throw FileParseException(
                "Expected  atomRefs for a MonomerSgroup definition in MRV file");
          }

          std::vector<std::string> atomList;
          boost::algorithm::split(atomList, atomRefsStr,
                                  boost::algorithm::is_space());
          for (auto &it : atomList) {
            auto atomPtr = res->parent->findAtomByRef(it);

            if (atomPtr == nullptr) {
              throw FileParseException(
                  "AtomRef specification for a MonomerSgroup definition was not found in any parent");
            }

            marvinMonomerSgroup->atoms.push_back(atomPtr);
          }

          marvinMonomerSgroup->title =
              molTree.get<std::string>("<xmlattr>.title", "");
          if (marvinMonomerSgroup->title == "") {
            throw FileParseException(
                "Expected  title for a MonomerSgroup definition in MRV file");
          }
          marvinMonomerSgroup->charge =
              molTree.get<std::string>("<xmlattr>.charge", "");
          if (marvinMonomerSgroup->charge == "") {
            throw FileParseException(
                "Expected  omAtoms or onBracket for a charge attr of a MonomerSgroup definition in MRV file");
          }
        } else {
          throw FileParseException("Unexpected role " + role + " in MRV file");
        }
      }

      res->molID = molID;

      // get atoms if this mol is supposed to have them

      if (res->hasAtomBondBlocks()) {
        boost::property_tree::ptree atomArray = molTree.get_child("atomArray");

        // there are two types of atom arrays:
        // <atomArray atomID="a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11" elementType="C
        // C C C C C Cl C N O O" formalCharge="0 0 0 0 0 0 0 0 1 0 -1"
        // lonePair="0 0 0 0 0 0 3 0 0 2 3" x2="-4.3334 -5.6670 -5.6670 -4.3334
        // -2.9997 -2.9997 -4.3335 -1.6660 -7.0007 -1.6660 -0.3323" y2="1.7693
        // 0.9993 -0.5409 -1.3109 -0.5409 0.9993 3.3093 -1.3109 -1.3109 -2.8509
        // -0.5410"></atomArray>
        //  AND
        // <atomArray>
        //       <atom id="a1" elementType="C" x2="-9.4583" y2="1.9358"
        //       mrvStereoGroup="and1"/> <atom id="a2" elementType="C"
        //       x2="-10.7921" y2="1.1658"/> <atom id="a3" elementType="C"
        //       x2="-10.7921" y2="-0.3744"/> <atom id="a8" elementType="O"
        //       x2="-12.1257" y2="-1.1444" lonePair="2"/>
        //   </atomArray>

        // See which one we have

        std::string atomID =
            atomArray.get<std::string>("<xmlattr>.atomID", "UseLongForm");
        if (atomID == "UseLongForm") {
          // long form - each atom on a line

          BOOST_FOREACH (boost::property_tree::ptree::value_type &v,
                         molTree.get_child("atomArray")) {
            if (v.first != "atom") {
              continue;
            }

            auto *mrvAtom = new MarvinAtom();
            res->atoms.push_back(mrvAtom);

            mrvAtom->id = v.second.get<std::string>("<xmlattr>.id", "");
            mrvAtom->elementType =
                v.second.get<std::string>("<xmlattr>.elementType", "");

            if (mrvAtom->id == "" || mrvAtom->elementType == "") {
              throw FileParseException(
                  "Expected id, elementType for an atom definition in MRV file");
            }

            std::string x2 = v.second.get<std::string>("<xmlattr>.x2", "");
            std::string y2 = v.second.get<std::string>("<xmlattr>.y2", "");

            // x2 and y2 are doubles

            if (x2 != "" && y2 != "" &&
                (!getCleanDouble(x2, mrvAtom->x2) ||
                 !getCleanDouble(y2, mrvAtom->y2))) {
              throw FileParseException(
                  "The values x2 and y2 must be large floats in MRV file");
            }

            std::string formalCharge =
                v.second.get<std::string>("<xmlattr>.formalCharge", "");
            if (formalCharge != "") {
              if (!getCleanInt(formalCharge, mrvAtom->formalCharge)) {
                throw FileParseException(
                    "The value for formalCharge must be an integer in MRV file");
              }
            } else {
              mrvAtom->formalCharge = 0;
            }

            mrvAtom->radical =
                v.second.get<std::string>("<xmlattr>.radical", "");
            if (mrvAtom->radical != "") {
              if (!boost::algorithm::contains(
                      marvinRadicalVals,
                      std::vector<std::string>{mrvAtom->radical})) {
                std::ostringstream err;
                err << "The value for radical must be one of "
                    << boost::algorithm::join(marvinRadicalVals, ", ")
                    << " in MRV file";
                throw FileParseException(err.str());
              }
            } else {
              mrvAtom->radical = "";
            }

            std::string isotopeStr =
                v.second.get<std::string>("<xmlattr>.isotope", "");
            if (isotopeStr != "") {
              if (!getCleanInt(isotopeStr, mrvAtom->isotope) ||
                  mrvAtom->isotope <= 0) {
                std::ostringstream err;
                err << "The value for isotope must be a positive number in MRV file";
                throw FileParseException(err.str());
              }
            } else {
              mrvAtom->isotope = 0;
            }

            std::string valenceStr =
                v.second.get<std::string>("<xmlattr>.mrvValence", "");
            if (valenceStr != "") {
              if (!getCleanInt(valenceStr, mrvAtom->mrvValence) ||
                  mrvAtom->mrvValence < 0) {
                std::ostringstream err;
                err << "The value for mrvValence must be a positive number in MRV file";
                throw FileParseException(err.str());
              }
            } else {
              mrvAtom->mrvValence = -1;
            }

            std::string hCountStr =
                v.second.get<std::string>("<xmlattr>.hydrogenCount", "");
            if (hCountStr != "") {
              if (!getCleanInt(hCountStr, mrvAtom->hydrogenCount) ||
                  mrvAtom->hydrogenCount < 0) {
                std::ostringstream err;
                err << "The value for hydrogenCount must be a non-negative number in MRV file";
                throw FileParseException(err.str());
              }
            } else {
              mrvAtom->hydrogenCount = -1;
            }

            mrvAtom->mrvAlias =
                v.second.get<std::string>("<xmlattr>.mrvAlias", "");

            mrvAtom->rgroupRef = v.second.get<int>("<xmlattr>.rgroupRef", -1);

            mrvAtom->mrvStereoGroup =
                v.second.get<std::string>("<xmlattr>.mrvStereoGroup", "");
            if (mrvAtom->mrvStereoGroup == "0") {
              mrvAtom->mrvStereoGroup = "";
            }

            std::string mrvMap =
                v.second.get<std::string>("<xmlattr>.mrvMap", "");
            if (mrvMap != "") {
              if (!getCleanInt(mrvMap, mrvAtom->mrvMap) ||
                  mrvAtom->mrvMap <= 0) {
                throw FileParseException(
                    "The value for mrvMap must be an non-=negative integer in MRV file");
              }
            } else {
              mrvAtom->mrvMap = 0;
            }

            mrvAtom->sgroupRef =
                v.second.get<std::string>("<xmlattr>.sgroupRef", "");

            if (role == "SuperatomSgroup") {
              mrvAtom->sgroupAttachmentPoint = v.second.get<std::string>(
                  "<xmlattr>.sgroupAttachmentPoint", "");

              // atom->setProp(common_properties::molAttachPoint, ival);
            }
          }
        } else  // single line form of atoms
        {
          // <atomArray atomID="a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11"
          // elementType="C C C C C C Cl C N O O" formalCharge="0 0 0 0 0 0 0 0
          // 1 0 -1" lonePair="0 0 0 0 0 0 3 0 0 2 3" x2="-4.3334 -5.6670
          // -5.6670 -4.3334 -2.9997 -2.9997 -4.3335 -1.6660 -7.0007 -1.6660
          // -0.3323" y2="1.7693 0.9993 -0.5409 -1.3109 -0.5409 0.9993 3.3093
          // -1.3109 -1.3109 -2.8509 -0.5410"></atomArray>

          std::vector<std::string> atomIds;
          size_t atomCount;
          if (atomID == "") {
            atomCount = 0;
          } else {
            boost::algorithm::split(atomIds, atomID,
                                    boost::algorithm::is_space());
            atomCount = atomIds.size();
          }

          std::vector<std::string> elementTypes;
          std::string elementType =
              atomArray.get<std::string>("<xmlattr>.elementType", "");
          boost::algorithm::split(elementTypes, elementType,
                                  boost::algorithm::is_space());

          std::vector<std::string> x2s;
          std::string x2 = atomArray.get<std::string>("<xmlattr>.x2", "");
          boost::algorithm::split(x2s, x2, boost::algorithm::is_space());

          std::vector<std::string> y2s;
          std::string y2 = atomArray.get<std::string>("<xmlattr>.y2", "");
          boost::algorithm::split(y2s, y2, boost::algorithm::is_space());

          std::vector<std::string> formalCharges;
          std::string formalCharge =
              atomArray.get<std::string>("<xmlattr>.formalCharge", "");
          boost::algorithm::split(formalCharges, formalCharge,
                                  boost::algorithm::is_space());

          std::vector<std::string> isotopes;
          std::string isotope =
              atomArray.get<std::string>("<xmlattr>.isotope", "");
          boost::algorithm::split(isotopes, isotope,
                                  boost::algorithm::is_space());

          std::vector<std::string> radicals;
          std::string radical =
              atomArray.get<std::string>("<xmlattr>.radical", "");
          boost::algorithm::split(radicals, radical,
                                  boost::algorithm::is_space());

          std::vector<std::string> hydrogenCounts;
          std::string hydrogenCount =
              atomArray.get<std::string>("<xmlattr>.hydrogenCount", "");
          boost::algorithm::split(hydrogenCounts, hydrogenCount,
                                  boost::algorithm::is_space());

          std::vector<std::string> mrvValences;
          std::string mrvValence =
              atomArray.get<std::string>("<xmlattr>.mrvValence", "");
          boost::algorithm::split(mrvValences, mrvValence,
                                  boost::algorithm::is_space());

          std::vector<std::string> mrvAliases;
          std::string mrvAlias =
              atomArray.get<std::string>("<xmlattr>.mrvAlias", "");
          boost::algorithm::split(mrvAliases, mrvAlias,
                                  boost::algorithm::is_space());

          std::vector<std::string> rgroupRefs;
          std::string rgroupRef =
              atomArray.get<std::string>("<xmlattr>.rgroupRef", "");
          boost::algorithm::split(rgroupRefs, rgroupRef,
                                  boost::algorithm::is_space());

          std::vector<std::string> mrvStereoGroups;
          std::string mrvStereoGroup =
              atomArray.get<std::string>("<xmlattr>.mrvStereoGroup", "");
          boost::algorithm::split(mrvStereoGroups, mrvStereoGroup,
                                  boost::algorithm::is_space());

          std::vector<std::string> mrvMaps;
          std::string mrvMap =
              atomArray.get<std::string>("<xmlattr>.mrvMap", "");
          boost::algorithm::split(mrvMaps, mrvMap,
                                  boost::algorithm::is_space());

          std::vector<std::string> sgroupRefs;
          std::string sgroupRef =
              atomArray.get<std::string>("<xmlattr>.sgroupRef", "");
          boost::algorithm::split(sgroupRefs, sgroupRef,
                                  boost::algorithm::is_space());

          std::vector<std::string> sgroupAttachmentPoints;
          std::string sgroupAttachmentPoint =
              atomArray.get<std::string>("<xmlattr>.sgroupAttachmentPoint", "");
          boost::algorithm::split(sgroupAttachmentPoints, sgroupAttachmentPoint,
                                  boost::algorithm::is_space());

          if (atomID != "") {
            if (elementType == "") {
              throw FileParseException(
                  "Expected an elementType array for an atomArray definition in MRV file");
            }
            if (elementTypes.size() < atomCount) {
              throw FileParseException(
                  "There must be an element type for each atom id");
            }
          }

          for (size_t i = 0; i < atomCount; ++i) {
            auto *mrvAtom = new MarvinAtom();
            res->atoms.push_back(mrvAtom);

            mrvAtom->id = atomIds[i];

            mrvAtom->elementType = elementTypes[i];

            if (x2 != "" && y2 != "" && x2s.size() > i && y2s.size() > i) {
              if (!getCleanDouble(x2s[i], mrvAtom->x2) ||
                  !getCleanDouble(y2s[i], mrvAtom->y2)) {
                throw FileParseException(
                    "The values x2 and y2 must be large floats in MRV file");
              }
            }

            if (formalCharge != "" && formalCharges.size() > i) {
              if (!getCleanInt(formalCharges[i], mrvAtom->formalCharge)) {
                throw FileParseException(
                    "The value for formalCharge must be an integer in MRV file");
              }
            } else {
              mrvAtom->formalCharge = 0;
            }

            if (isotope != "" && isotopes.size() > i) {
              if (!getCleanInt(isotopes[i], mrvAtom->isotope)) {
                throw FileParseException(
                    "The value for formalCharge must be an integer in MRV file");
              }
            } else {
              mrvAtom->isotope = 0;
            }

            if (mrvValence != "" && mrvValences.size() > i) {
              if (mrvValences[i] == "-") {
                mrvAtom->mrvValence = -1;
              } else if (!getCleanInt(mrvValences[i], mrvAtom->mrvValence)) {
                throw FileParseException(
                    "The value for mrvValences must be an integer in MRV file");
              }
            } else {
              mrvAtom->mrvValence = -1;
            }

            if (hydrogenCount != "" && hydrogenCounts.size() > i) {
              if (hydrogenCounts[i] != "-" &&
                  !getCleanInt(hydrogenCounts[i], mrvAtom->hydrogenCount)) {
                throw FileParseException(
                    "The value for hydrogenCount must be an integer in MRV file");
              }
            } else {
              mrvAtom->hydrogenCount = -1;
            }

            if (radical != "" && radicals.size() > i && radicals[i] != "0") {
              mrvAtom->radical = radicals[i];
              if (!boost::algorithm::contains(
                      marvinRadicalVals,
                      std::vector<std::string>{mrvAtom->radical})) {
                std::ostringstream err;
                err << "The value for radical must be one of "
                    << boost::algorithm::join(marvinRadicalVals, ", ")
                    << " in MRV file";
                throw FileParseException(err.str());
              }
            } else {
              mrvAtom->radical = "";
            }

            if (mrvAlias != "" && mrvAliases.size() > i) {
              mrvAtom->mrvAlias = mrvAliases[i];
            } else {
              mrvAtom->mrvAlias = "";
            }

            if (rgroupRef != "" && rgroupRefs.size() > i) {
              if (!getCleanInt(rgroupRefs[i], mrvAtom->rgroupRef)) {
                throw FileParseException(
                    "rgroupRef value must be an integer in MRV file");
              }
            } else {
              mrvAtom->rgroupRef = (-1);
            }

            if (mrvStereoGroup != "" && mrvStereoGroups.size() > i &&
                mrvStereoGroups[i] != "0")  // "0" is NOT a stereo group
            {
              mrvAtom->mrvStereoGroup = mrvStereoGroups[i];
            } else {
              mrvAtom->mrvStereoGroup = "";
            }

            if (mrvMap != "" && mrvMaps.size() > i) {
              if (!getCleanInt(mrvMaps[i], mrvAtom->mrvMap) ||
                  mrvAtom->mrvMap < 0) {
                throw FileParseException(
                    "The value for mrvMap must be an non-negative integer in MRV file");
              }
            } else {
              mrvAtom->mrvMap = 0;
            }

            if (sgroupRef != "" && sgroupRefs.size() > i &&
                sgroupRefs[i] != "0") {
              mrvAtom->sgroupRef = sgroupRefs[i];
            } else {
              mrvAtom->sgroupRef = "";
            }

            if (role == "SuperatomSgroup" && sgroupAttachmentPoint != "" &&
                sgroupAttachmentPoints.size() > i &&
                sgroupAttachmentPoints[i] != "0") {
              mrvAtom->sgroupAttachmentPoint = sgroupAttachmentPoints[i];
            } else {
              mrvAtom->sgroupAttachmentPoint = "";
            }
          }
        }

        auto bondArray = molTree.get_child_optional("bondArray");
        if (bondArray) {
          BOOST_FOREACH (boost::property_tree::ptree::value_type &v,
                         molTree.get_child("bondArray")) {
            if (v.first != "bond") {
              continue;
            }

            auto *mrvBond = new MarvinBond();
            res->bonds.push_back(mrvBond);

            mrvBond->id = v.second.get<std::string>("<xmlattr>.id", "");
            if (mrvBond->id == "") {
              throw FileParseException(
                  "Expected id for an bond definition in MRV file");
            }

            std::string atomRefs2 =
                v.second.get<std::string>("<xmlattr>.atomRefs2", "");

            std::vector<std::string> atomRefs2s;
            boost::algorithm::split(atomRefs2s, atomRefs2,
                                    boost::algorithm::is_space());
            mrvBond->atomRefs2[0] = atomRefs2s[0];
            mrvBond->atomRefs2[1] = atomRefs2s[1];
            if (atomRefs2s.size() != 2 ||
                !boost::algorithm::contains(
                    res->atoms, std::vector<std::string>{mrvBond->atomRefs2[0]},
                    MarvinMol::atomRefInAtoms) ||
                !boost::algorithm::contains(
                    res->atoms, std::vector<std::string>{mrvBond->atomRefs2[1]},
                    MarvinMol::atomRefInAtoms)) {
              throw FileParseException(
                  "atomRefs2 must contain two atom refs that must appear in the atoms array in MRV file");
            }

            mrvBond->order = v.second.get<std::string>("<xmlattr>.order", "");
            if (mrvBond->order != "") {
              if (!boost::algorithm::contains(
                      marvinBondOrders,
                      std::vector<std::string>{mrvBond->order})) {
                std::ostringstream err;
                err << "Expected one of  "
                    << boost::algorithm::join(marvinBondOrders, ", ")
                    << " for order for an bond definition in MRV file";
                throw FileParseException(err.str());
              }
            }

            mrvBond->queryType =
                v.second.get<std::string>("<xmlattr>.queryType", "");
            if (mrvBond->queryType != "") {
              if (!boost::algorithm::contains(
                      marvinQueryBondsTypes,
                      std::vector<std::string>{mrvBond->queryType})) {
                std::ostringstream err;
                err << "Expected one of  "
                    << boost::algorithm::join(marvinQueryBondsTypes, ", ")
                    << " for queryType for an bond definition in MRV file";
                throw FileParseException(err.str());
              }
            }

            mrvBond->convention =
                v.second.get<std::string>("<xmlattr>.convention", "");
            if (mrvBond->convention != "") {
              if (!boost::algorithm::contains(
                      marvinConventionTypes,
                      std::vector<std::string>{mrvBond->convention})) {
                std::ostringstream err;
                err << "Expected one of  "
                    << boost::algorithm::join(marvinConventionTypes, ", ")
                    << " for convention for an bond definition in MRV file";
                throw FileParseException(err.str());
              }
            }

            int bondStereoDeclCount = 0;
            mrvBond->bondStereo.value =
                v.second.get<std::string>("bondStereo", "");
            if (mrvBond->bondStereo.value != "") {
              bondStereoDeclCount++;
              if (boost::algorithm::to_lower_copy(mrvBond->bondStereo.value) ==
                      "w" ||
                  boost::algorithm::to_lower_copy(mrvBond->bondStereo.value) ==
                      "h") {
                // do nothing  - this is OK
              } else if (boost::algorithm::to_lower_copy(
                             mrvBond->bondStereo.value) == "c" ||
                         boost::algorithm::to_lower_copy(
                             mrvBond->bondStereo.value) == "t") {
                mrvBond->bondStereo.value = "";  // cis and trans are ignored
              } else {
                throw FileParseException(
                    "The value for bondStereo must be \"H\", \"W\", \"C\" or \"T\" in MRV file (\"C\" and \"T\" are ignored)");
              }
            }

            // see if bondstereo has a dictRef or convention

            auto bondStereoItem = v.second.get_child_optional("bondStereo");
            if (bondStereoItem) {
              BOOST_FOREACH (boost::property_tree::ptree::value_type &ww,
                             v.second) {
                mrvBond->bondStereo.convention =
                    ww.second.get<std::string>("<xmlattr>.convention", "");
                if (mrvBond->bondStereo.convention != "") {
                  bondStereoDeclCount++;
                  if (mrvBond->bondStereo.convention != "MDL") {
                    std::ostringstream err;
                    err << "Expected MDL as value for the bond convention attribute";
                    throw FileParseException(err.str());
                  }
                  mrvBond->bondStereo.conventionValue =
                      ww.second.get<std::string>("<xmlattr>.conventionValue",
                                                 "");
                  if (!boost::algorithm::contains(
                          marvinStereoConventionTypes,
                          std::vector<std::string>{
                              mrvBond->bondStereo.conventionValue})) {
                    std::ostringstream err;
                    err << "Expected one of  "
                        << boost::algorithm::join(marvinStereoConventionTypes,
                                                  ", ")
                        << " for a bond convention for an bond stereo def";
                    throw FileParseException(err.str());
                  }
                }

                mrvBond->bondStereo.dictRef =
                    ww.second.get<std::string>("<xmlattr>.dictRef", "");
                if (mrvBond->bondStereo.dictRef != "") {
                  bondStereoDeclCount++;
                  if (!boost::algorithm::contains(
                          marvinStereoDictRefTypes,
                          std::vector<std::string>{
                              mrvBond->bondStereo.dictRef})) {
                    std::ostringstream err;
                    err << "Expected one of  "
                        << boost::algorithm::join(marvinStereoDictRefTypes,
                                                  ", ")
                        << " for a discRef value for an bond stereo def";
                    throw FileParseException(err.str());
                  }
                }
              }
            }

            // check that there were not too many different declarations

            if (bondStereoDeclCount > 1) {
              throw FileParseException(
                  "bondStereo can either have only one of: a value, dictRef, or convention value");
            }
          }
        }
      }

      if (role == "SuperatomSgroup") {
        // see if there is an AttachmentPointArray
        bool found;

        try {
          boost::property_tree::ptree AttachmentPointArrayTree =
              molTree.get_child("AttachmentPointArray");
          found = true;
        } catch (const std::exception &e) {
          found = false;
        }

        if (found) {
          BOOST_FOREACH (boost::property_tree::ptree::value_type &v,
                         molTree.get_child("AttachmentPointArray")) {
            std::string bondId =
                v.second.get<std::string>("<xmlattr>.bond", "");
            if (bondId == "") {  // this can happen if the attachment point is
                                 // not actually used - as in Amino acids that
                                 // have non-used crosslink atoms
              continue;
            }

            auto *marvinAttachmentPoint = new MarvinAttachmentPoint();
            ((MarvinSuperatomSgroup *)res)
                ->attachmentPoints.push_back(marvinAttachmentPoint);

            marvinAttachmentPoint->atom =
                v.second.get<std::string>("<xmlattr>.atom", "");
            marvinAttachmentPoint->order =
                v.second.get<std::string>("<xmlattr>.order", "");
            marvinAttachmentPoint->bond = bondId;

            if (marvinAttachmentPoint->atom == "" ||
                marvinAttachmentPoint->order == "") {
              throw FileParseException(
                  "Expected atom, order and bond, for an AttachmentPoint definition in MRV file");
            }

            if (!boost::algorithm::contains(
                    res->atoms,
                    std::vector<std::string>{marvinAttachmentPoint->atom},
                    MarvinMol::atomRefInAtoms)) {
              throw FileParseException(
                  "Atom specification for an AttachmentPoint definition must be in the parent's atom array in MRV file");
            }

            // bond must be found in the  bonds vector

            if (parentMol->findBondByRef(marvinAttachmentPoint->bond) ==
                nullptr) {
              throw FileParseException(
                  "Bond specification for an AttachmentPoint definition must be in the bond array in MRV file");
            }

            // the order must be an integer

            int orderInt;
            if (!getCleanInt(marvinAttachmentPoint->order, orderInt)) {
              throw FileParseException(
                  "Order for an AttachmentPoint definition must be an integer in MRV file");
            }
          }
        }
      }

      BOOST_FOREACH (boost::property_tree::ptree::value_type &v, molTree) {
        if (v.first == "molecule") {
          MarvinMolBase *subMol =
              parseMarvinMolecule(v.second, (MarvinMol *)res);
          ((MarvinMol *)res)->sgroups.push_back(subMol);
        }
      }

      return res;
    } catch (const std::exception &e) {
      delete res;

      throw;
    }
  };

  MarvinReaction *parseMarvinReaction(
      boost::property_tree::ptree rxnTree,
      boost::property_tree::ptree documentTree) {
    auto *res = new MarvinReaction();

    try {
      boost::property_tree::ptree childTree;
      bool foundChild = false;
      try {
        childTree = rxnTree.get_child("reactantList");
        foundChild = true;
      } catch (const std::exception &e) {
        foundChild = false;
      }

      if (foundChild) {
        BOOST_FOREACH (boost::property_tree::ptree::value_type &v, childTree)
          res->reactants.push_back((MarvinMol *)parseMarvinMolecule(v.second));
      }

      try {
        childTree = rxnTree.get_child("agentList");
        foundChild = true;
      } catch (const std::exception &e) {
        foundChild = false;
      }
      if (foundChild) {
        BOOST_FOREACH (boost::property_tree::ptree::value_type &v, childTree)
          res->agents.push_back((MarvinMol *)parseMarvinMolecule(v.second));
      }

      try {
        childTree = rxnTree.get_child("productList");
        foundChild = true;
      } catch (const std::exception &e) {
        foundChild = false;
      }
      if (foundChild) {
        BOOST_FOREACH (boost::property_tree::ptree::value_type &v, childTree)
          res->products.push_back((MarvinMol *)parseMarvinMolecule(v.second));
      }

      // <arrow type="DEFAULT" x1="-11.816189911577812" y1="-10.001443743444021"
      // x2="-8.401759471454618" y2="-10.001443743444021"/>
      boost::property_tree::ptree arrow = rxnTree.get_child("arrow");
      res->arrow.type = arrow.get<std::string>("<xmlattr>.type", "");
      if (!getCleanDouble(arrow.get<std::string>("<xmlattr>.x1", ""),
                          res->arrow.x1) ||
          !getCleanDouble(arrow.get<std::string>("<xmlattr>.y1", ""),
                          res->arrow.y1) ||
          !getCleanDouble(arrow.get<std::string>("<xmlattr>.x2", ""),
                          res->arrow.x2) ||
          !getCleanDouble(arrow.get<std::string>("<xmlattr>.y2", ""),
                          res->arrow.y1)) {
        throw FileParseException(
            "Arrow coordinates must all be large floating point numbers in MRV file");
      }

      BOOST_FOREACH (boost::property_tree::ptree::value_type &v, documentTree) {
        if (v.first != "MReactionSign") {
          continue;
        }
        auto *marvinPlus = new MarvinPlus();
        res->pluses.push_back(marvinPlus);
        marvinPlus->id = v.second.get<std::string>("<xmlattr>.id", "");
        int pointCount = 0;
        BOOST_FOREACH (boost::property_tree::ptree::value_type &v2, v.second) {
          if (v2.first == "MPoint") {
            double x;
            double y;
            if (!getCleanDouble(v2.second.get<std::string>("<xmlattr>.x", ""),
                                x) ||
                !getCleanDouble(v2.second.get<std::string>("<xmlattr>.y", ""),
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

      BOOST_FOREACH (boost::property_tree::ptree::value_type &v, documentTree) {
        if (v.first != "MTextBox") {
          continue;
        }

        auto *marvinCondition = new MarvinCondition();
        res->conditions.push_back(marvinCondition);
        marvinCondition->id = v.second.get<std::string>("<xmlattr>.id", "");
        marvinCondition->halign =
            v.second.get<std::string>("<xmlattr>.halign", "");
        marvinCondition->valign =
            v.second.get<std::string>("<xmlattr>.valign", "");
        double fontScale;
        std::string fontScaleStr =
            v.second.get<std::string>("<xmlattr>.fontScale", "");
        if (fontScaleStr != "") {
          if (!getCleanDouble(fontScaleStr, fontScale)) {
            throw FileParseException(
                "Condition font scale must be a positive integer in MRV file");
          }
        } else {
          fontScale = 0.0;
        }
        marvinCondition->fontScale = fontScale;

        marvinCondition->text = v.second.get<std::string>("Field", "");

        int pointCount = 0;
        BOOST_FOREACH (boost::property_tree::ptree::value_type &v2, v.second) {
          if (v2.first == "MPoint") {
            double x, y;
            std::string xStr = v2.second.get<std::string>("<xmlattr>.x", "");
            std::string yStr = v2.second.get<std::string>("<xmlattr>.y", "");
            if (!getCleanDouble(xStr, x) || !getCleanDouble(yStr, y)) {
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

      return res;
    } catch (const std::exception &e) {
      delete res;

      throw;
    }
  }
};

//------------------------------------------------
//
//  Read a molecule from a stream
//
//------------------------------------------------

bool MrvDataStreamIsReaction(std::istream &inStream) {
  PRECONDITION(inStream, "no stream");

  using boost::property_tree::ptree;
  ptree tree;

  // Parse the XML into the property tree.

  read_xml(inStream, tree);

  // see if the reaction header is present
  try {
    auto rxn = tree.get_child("cml.MDocument.MChemicalStruct.reaction");
  } catch (const std::exception &e) {
    return false;
  }

  return true;
}

bool MrvDataStreamIsReaction(std::istream *inStream) {
  return MrvDataStreamIsReaction(*inStream);
}
//------------------------------------------------
//
//  Read a molecule from a string
//
//------------------------------------------------
bool MrvBlockIsReaction(const std::string &molmrvText) {
  std::istringstream inStream(molmrvText);
  return MrvDataStreamIsReaction(inStream);
}

//------------------------------------------------
//
//  Read a RWMOL from a file
//
//------------------------------------------------
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
RWMol *MrvMolDataStreamParser(std::istream *inStream, bool sanitize,
                              bool removeHs) {
  using boost::property_tree::ptree;
  ptree tree;

  // Parse the XML into the property tree.

  read_xml(*inStream, tree);

  MarvinCMLReader reader;
  return reader.parseMolecule(tree, sanitize, removeHs);
}
//------------------------------------------------
//
//  Read a RWMol from a stream reference
//
//------------------------------------------------
RWMol *MrvMolDataStreamParser(std::istream &inStream, bool sanitize,
                              bool removeHs) {
  return MrvMolDataStreamParser(&inStream, sanitize, removeHs);
}
//------------------------------------------------
//
//  Read a RWMol from a string
//
//------------------------------------------------
RWMol *MrvMolStringParser(const std::string &molmrvText, bool sanitize,
                          bool removeHs) {
  std::istringstream inStream(molmrvText);
  return MrvMolDataStreamParser(inStream, sanitize, removeHs);
}

//------------------------------------------------
//
//  Read an RWMol from a file
//
//------------------------------------------------
RWMol *MrvMolFileParser(const std::string &fName, bool sanitize,
                        bool removeHs) {
  std::ifstream inStream(fName.c_str());
  if (!inStream || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fName;
    throw BadFileException(errout.str());
  }
  RWMol *res = nullptr;
  if (!inStream.eof()) {
    res = MrvMolDataStreamParser(inStream, sanitize, removeHs);
  }
  return res;
}

//------------------------------------------------
//
//  Read a ChemicalReaction from a stream
//
//------------------------------------------------
ChemicalReaction *MrvRxnDataStreamParser(std::istream *inStream, bool sanitize,
                                         bool removeHs) {
  using boost::property_tree::ptree;
  ptree tree;

  // Parse the XML into the property tree.

  read_xml(*inStream, tree);

  MarvinCMLReader reader;
  return reader.parseReaction(tree, tree.get_child("cml.MDocument"), sanitize,
                              removeHs);
}

//------------------------------------------------
//
//  Read a ChemicalReaction from a stream reference
//
//------------------------------------------------
ChemicalReaction *MrvRxnDataStreamParser(std::istream &inStream, bool sanitize,
                                         bool removeHs) {
  return MrvRxnDataStreamParser(&inStream, sanitize, removeHs);
}
//------------------------------------------------
//
//  Read a ChemicalReaction from a string
//
//------------------------------------------------
ChemicalReaction *MrvRxnStringParser(const std::string &molmrvText,
                                     bool sanitize, bool removeHs) {
  std::istringstream inStream(molmrvText);
  return MrvRxnDataStreamParser(inStream, sanitize, removeHs);
}

//------------------------------------------------
//
//  Read a ChemicalReaction from a file
//
//------------------------------------------------
ChemicalReaction *MrvRxnFileParser(const std::string &fName, bool sanitize,
                                   bool removeHs) {
  std::ifstream inStream(fName.c_str());
  if (!inStream || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fName;
    throw BadFileException(errout.str());
  }
  ChemicalReaction *res = MrvRxnDataStreamParser(inStream, sanitize, removeHs);
  return res;
}
}  // namespace RDKit
