//
//  Copyright (C) 2018-2022 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifdef _MSC_VER
#pragma warning(disable : 4503)
#endif

#include <RDGeneral/Invariant.h>
#include <RDGeneral/versions.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/StereoGroup.h>

#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/MolInterchange/details.h>
#include <RDGeneral/FileParseException.h>

#include <sstream>
#include <exception>
#include <map>
#include <cmath>

#include <RDGeneral/BoostStartInclude.h>
#define BOOST_JSON_NO_LIB
#ifndef BOOST_CONTAINER_NO_LIB
#define BOOST_CONTAINER_NO_LIB
#endif
#include <boost/json.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace bj = boost::json;

namespace RDKit {

namespace MolInterchange {

namespace {

template <typename T>
void addMol(const T &imol, bj::object &bjMol, const bj::object &atomDefaults,
            const bj::object &bondDefaults, const JSONWriteParameters &params);

void initAtomDefaults(bj::object &bjDefaults) {
  // "atomDefaults": {"Z": 6, "impHs": 0, "chg": 0, "stereo": "unspecified",
  // "nrad": 0, "isotope": 0},
  bjDefaults["z"] = 6;
  bjDefaults["impHs"] = 0;
  bjDefaults["chg"] = 0;
  bjDefaults["nRad"] = 0;
  bjDefaults["isotope"] = 0;
  bjDefaults["stereo"] = "unspecified";
}
void initBondDefaults(bj::object &bjDefaults) {
  // "bondDefaults": {"bo": 1, "stereo": "unspecified", "stereoAtoms": []},
  bjDefaults["bo"] = 1;
  bjDefaults["stereo"] = "unspecified";
}
void initHeader(bj::object &bjHeader, const JSONWriteParameters &params) {
  auto vers = currentRDKitJSONVersion;
  if (!params.useRDKitExtensions) {
    vers = currentMolJSONVersion;
  }
  bjHeader["version"] = vers;
}

void addIntVal(bj::object &dest, const bj::object &defaults, const char *tag,
               int val) {
  auto it = defaults.find(tag);
  if (it != defaults.end()) {
    int dval = it->value().as_int64();
    if (dval != val) {
      dest[tag] = val;
    }
  } else {
    dest[tag] = val;
  }
}

void addIntVal(bj::object &dest, const char *tag, int val) { dest[tag] = val; }

void addStringVal(bj::object &dest, const bj::object &defaults, const char *tag,
                  const std::string &val) {
  auto it = defaults.find(tag);
  if (it != defaults.end()) {
    std::string dval = it->value().as_string().c_str();
    if (val.size() && dval != val) {
      dest[tag] = val;
    }
  } else {
    dest[tag] = val;
  }
}

void addStringVal(bj::object &dest, const char *tag, const std::string &val) {
  dest[tag] = val;
}

void addAtom(const Atom &atom, bj::object &bjAtom,
             const bj::object &bjDefaults) {
  addIntVal(bjAtom, bjDefaults, "z", atom.getAtomicNum());
  if (!atom.hasQuery()) {
    addIntVal(bjAtom, bjDefaults, "impHs", atom.getTotalNumHs());
  }
  addIntVal(bjAtom, bjDefaults, "chg", atom.getFormalCharge());
  addIntVal(bjAtom, bjDefaults, "isotope", atom.getIsotope());
  addIntVal(bjAtom, bjDefaults, "nRad", atom.getNumRadicalElectrons());
  std::string chi = "";
  if (inv_chilookup.find(atom.getChiralTag()) != inv_chilookup.end()) {
    chi = inv_chilookup.find(atom.getChiralTag())->second;
  } else {
    BOOST_LOG(rdWarningLog)
        << " unrecognized atom chirality set to default while writing"
        << std::endl;
  }
  addStringVal(bjAtom, bjDefaults, "stereo", chi);
}

template <typename Q>
void addQuery(const Q &query, bj::object &bjQuery,
              const JSONWriteParameters &params) {
  bjQuery["descr"] = query.getDescription();
  if (!query.getTypeLabel().empty()) {
    bjQuery["type"] = query.getTypeLabel();
  }
  if (query.getNegation()) {
    bjQuery["negated"] = true;
  }
  if (typeid(query) == typeid(RecursiveStructureQuery)) {
    auto rq = (const RecursiveStructureQuery *)&query;
    auto submol = rq->getQueryMol();
    PRECONDITION(submol, "bad recursive query");
    bj::object subquery;
    bj::object atomDefaults, bondDefaults;
    initAtomDefaults(atomDefaults);
    initBondDefaults(bondDefaults);
    addMol(*submol, subquery, atomDefaults, bondDefaults, params);
    bjQuery["tag"] = MolPickler::QUERY_RECURSIVE;
    bjQuery["subquery"] = subquery;
  } else {
    auto qdetails = PicklerOps::getQueryDetails(&query);
    switch (qdetails.which()) {
      case 0:
        bjQuery["tag"] = boost::get<MolPickler::Tags>(qdetails);
        break;
      case 1: {
        auto v = boost::get<std::tuple<MolPickler::Tags, int32_t>>(qdetails);
        bjQuery["tag"] = std::get<0>(v);
        bjQuery["val"] = std::get<1>(v);
      } break;
      case 2: {
        auto v = boost::get<std::tuple<MolPickler::Tags, int32_t, int32_t>>(
            qdetails);
        bjQuery["tag"] = std::get<0>(v);
        bjQuery["val"] = std::get<1>(v);
        if (std::get<2>(v)) {
          bjQuery["tol"] = std::get<2>(v);
        }
      } break;
      case 3: {
        auto v = boost::get<
            std::tuple<MolPickler::Tags, int32_t, int32_t, int32_t, char>>(
            qdetails);
        bjQuery["tag"] = std::get<0>(v);
        bjQuery["lower"] = std::get<1>(v);
        bjQuery["upper"] = std::get<2>(v);
        if (std::get<3>(v)) {
          bjQuery["tol"] = std::get<3>(v);
        }
        bjQuery["ends"] = std::get<4>(v);
      } break;
      case 4: {
        auto v = boost::get<std::tuple<MolPickler::Tags, std::set<int32_t>>>(
            qdetails);
        bjQuery["tag"] = std::get<0>(v);
        const auto &tset = std::get<1>(v);
        bj::array sval;
        for (auto val : tset) {
          sval.push_back(val);
        }
        bjQuery["set"] = sval;

      } break;
      default:
        throw MolPicklerException(
            "do not know how to pickle part of the query.");
    }
  }
  if (query.endChildren() != query.beginChildren()) {
    bj::array children;
    for (auto cit = query.beginChildren(); cit != query.endChildren(); ++cit) {
      bj::object child;
      addQuery(**cit, child, params);
      children.push_back(std::move(child));
    }
    bjQuery["children"] = std::move(children);
  }
}

void addBond(const Bond &bond, bj::object &bjBond, const bj::object &bjDefaults,
             const JSONWriteParameters &params) {
  int bo = -1;
  if (auto boIter = inv_bolookup.find(bond.getBondType());
      boIter != inv_bolookup.end()) {
    bo = boIter->second;
  }
  // commonchem only supports a few bond orders:
  if (!params.useRDKitExtensions && bo > 3) {
    bo = -1;
  }
  if (bo < 0) {
    bo = 0;
    BOOST_LOG(rdWarningLog) << " unrecognized bond type " << bond.getBondType()
                            << " set to zero while writing" << std::endl;
  }
  addIntVal(bjBond, bjDefaults, "bo", bo);
  bj::array bjAtoms;
  bjAtoms.push_back(static_cast<int>(bond.getBeginAtomIdx()));
  bjAtoms.push_back(static_cast<int>(bond.getEndAtomIdx()));
  bjBond["atoms"] = std::move(bjAtoms);

  std::string chi = "";
  if (auto sbIter = inv_stereoBondlookup.find(bond.getStereo());
      sbIter != inv_stereoBondlookup.end()) {
    chi = sbIter->second;
  } else {
    BOOST_LOG(rdWarningLog) << " unrecognized bond stereo " << bond.getStereo()
                            << " set to default while writing" << std::endl;
  }
  addStringVal(bjBond, bjDefaults, "stereo", chi);
  if (chi != "unspecified" && bond.getStereoAtoms().size() == 2) {
    bj::array bjStereoAtoms;
    bjStereoAtoms.push_back(static_cast<int>(bond.getStereoAtoms()[0]));
    bjStereoAtoms.push_back(static_cast<int>(bond.getStereoAtoms()[1]));
    bjBond["stereoAtoms"] = std::move(bjStereoAtoms);
  }
}

template <typename T>
void addProperties(const T &obj, const std::vector<std::string> &propNames,
                   bj::object &properties) {
  const auto &rd_dict = obj.getDict();

  for (const auto &rdvalue : rd_dict) {
    if (std::find(propNames.begin(), propNames.end(), rdvalue.key) ==
        propNames.end()) {
      continue;
    }
    const auto tag = rdvalue.val.getTag();
    switch (tag) {
      case RDTypeTag::IntTag:
      case RDTypeTag::UnsignedIntTag:
        properties[rdvalue.key] = from_rdvalue<int>(rdvalue.val);
        break;
      case RDTypeTag::DoubleTag:
      case RDTypeTag::FloatTag:
        properties[rdvalue.key] = from_rdvalue<double>(rdvalue.val);
        break;
      default:
        try {
          properties[rdvalue.key] = from_rdvalue<std::string>(rdvalue.val);
        } catch (const std::bad_any_cast &) {
          BOOST_LOG(rdWarningLog)
              << "Warning: Could not convert property " << rdvalue.key
              << " to a recognized type. Skipping it." << std::endl;
        }
        break;
    }
  }
}

void addSubstanceGroup(const SubstanceGroup &sg, bj::object &bjSG) {
  bool includePrivate = false, includeComputed = false;
  auto propNames = sg.getPropList(includePrivate, includeComputed);
  if (propNames.size()) {
    bj::object properties;
    addProperties(sg, propNames, properties);
    bjSG["properties"] = std::move(properties);
  }

  if (!sg.getAtoms().empty()) {
    bj::array bjArr;
    for (const auto idx : sg.getAtoms()) {
      bjArr.push_back(static_cast<int>(idx));
    }
    bjSG["atoms"] = std::move(bjArr);
  }
  if (!sg.getBonds().empty()) {
    bj::array bjArr;
    for (const auto idx : sg.getBonds()) {
      bjArr.push_back(static_cast<int>(idx));
    }
    bjSG["bonds"] = std::move(bjArr);
  }
  if (!sg.getParentAtoms().empty()) {
    bj::array bjArr;
    for (const auto idx : sg.getParentAtoms()) {
      bjArr.push_back(static_cast<int>(idx));
    }
    bjSG["parentAtoms"] = std::move(bjArr);
  }
  if (!sg.getBrackets().empty()) {
    bj::array bjArr;
    for (const auto &brk : sg.getBrackets()) {
      bj::array bjBrk;
      for (const auto &pt : brk) {
        bj::array bjPos;
        bjPos.push_back(pt.x);
        bjPos.push_back(pt.y);
        bjPos.push_back(pt.z);
        bjBrk.push_back(bjPos);
      }
      bjArr.push_back(std::move(bjBrk));
    }
    bjSG["brackets"] = std::move(bjArr);
  }

  if (!sg.getCStates().empty()) {
    bj::array bjArr;
    for (const auto &cs : sg.getCStates()) {
      bj::object bjCS;
      bjCS["bond"] = cs.bondIdx;
      if ("SUP" == sg.getProp<std::string>("TYPE")) {
        bj::array bjLoc;
        bjLoc.push_back(cs.vector.x);
        bjLoc.push_back(cs.vector.y);
        bjLoc.push_back(cs.vector.z);
        bjCS["vector"] = bjLoc;
      }
      bjArr.push_back(std::move(bjCS));
    }
    bjSG["cstates"] = std::move(bjArr);
  }

  if (!sg.getAttachPoints().empty()) {
    bj::array bjArr;
    for (const auto &ap : sg.getAttachPoints()) {
      bj::object bjAP;
      bjAP["aIdx"] = ap.aIdx;
      if (ap.lvIdx != -1) {
        bjAP["lvIdx"] = ap.lvIdx;
      }
      if (!ap.id.empty()) {
        bjAP["id"] = ap.id;
      }
      bjArr.push_back(std::move(bjAP));
    }
    bjSG["attachPoints"] = std::move(bjArr);
  }
}

void addStereoGroup(const StereoGroup &sg, bj::object &bjSG) {
  if (inv_stereoGrouplookup.find(sg.getGroupType()) ==
      inv_stereoGrouplookup.end()) {
    throw ValueErrorException("unrecognized StereoGroup type");
  }
  bjSG["type"] = inv_stereoGrouplookup.at(sg.getGroupType());

  if (sg.getGroupType() != StereoGroupType::STEREO_ABSOLUTE) {
    bjSG["id"] = sg.getWriteId();
  }

  if (!sg.getAtoms().empty()) {
    bj::array bjAtoms;
    for (const auto atm : sg.getAtoms()) {
      bjAtoms.push_back(static_cast<int>(atm->getIdx()));
    }
    bjSG["atoms"] = std::move(bjAtoms);
  }
  if (!sg.getBonds().empty()) {
    bj::array bjBonds;
    for (const auto bnd : sg.getBonds()) {
      bjBonds.push_back(static_cast<int>(bnd->getIdx()));
    }
    bjSG["bonds"] = std::move(bjBonds);
  }
}

void addConformer(const Conformer &conf, bj::object &bjConf) {
  int dim = 2;
  if (conf.is3D()) {
    dim = 3;
  }
  bjConf["dim"] = dim;
  bj::array bjCoords;
  for (const auto &pos : conf.getPositions()) {
    bj::array bjPos;
    bjPos.push_back(pos.x);
    bjPos.push_back(pos.y);
    if (dim == 3) {
      bjPos.push_back(pos.z);
    }
    bjCoords.push_back(std::move(bjPos));
  }
  bjConf["coords"] = std::move(bjCoords);
}

template <typename T>
void addMol(const T &imol, bj::object &rjMol, const bj::object &atomDefaults,
            const bj::object &bondDefaults, const JSONWriteParameters &params) {
  RWMol mol(imol);
  if (!mol.getRingInfo()->isSymmSssr()) {
    MolOps::symmetrizeSSSR(mol);
  }
  if (mol.needsUpdatePropertyCache()) {
    mol.updatePropertyCache(false);
  }
  try {
    MolOps::Kekulize(mol, false);
  } catch (const KekulizeException &) {
    mol = imol;
    if (mol.needsUpdatePropertyCache()) {
      mol.updatePropertyCache(false);
    }
  }
  if (mol.hasProp(common_properties::_Name)) {
    const std::string &nm = mol.getProp<std::string>(common_properties::_Name);
    rjMol["name"] = nm;
  }
  bool hasQueryAtoms = false;
  {
    bj::array rjAtoms;
    for (const auto &at : mol.atoms()) {
      bj::object rjAtom;
      addAtom(*at, rjAtom, atomDefaults);
      rjAtoms.push_back(std::move(rjAtom));
      if (at->hasQuery()) {
        hasQueryAtoms = true;
      }
    }
    rjMol["atoms"] = std::move(rjAtoms);
  }
  bool hasQueryBonds = false;
  {
    bj::array rjBonds;
    for (const auto &bnd : mol.bonds()) {
      bj::object rjBond;
      addBond(*bnd, rjBond, bondDefaults, params);
      rjBonds.push_back(std::move(rjBond));
      if (bnd->hasQuery()) {
        hasQueryBonds = true;
      }
    }
    rjMol["bonds"] = std::move(rjBonds);
  }
  if (params.useRDKitExtensions && !mol.getStereoGroups().empty()) {
    bj::array rjStereoGroups;
    for (const auto &sg : mol.getStereoGroups()) {
      bj::object rjSG;
      addStereoGroup(sg, rjSG);
      rjStereoGroups.push_back(std::move(rjSG));
    }
    rjMol["stereoGroups"] = std::move(rjStereoGroups);
  }

  if (params.useRDKitExtensions && !getSubstanceGroups(mol).empty()) {
    bj::array rjSubstanceGroups;
    for (const auto &sg : getSubstanceGroups(mol)) {
      bj::object rjSG;
      addSubstanceGroup(sg, rjSG);
      rjSubstanceGroups.push_back(std::move(rjSG));
    }
    rjMol["substanceGroups"] = std::move(rjSubstanceGroups);
  }

  if (mol.getNumConformers()) {
    bj::array rjConfs;
    for (auto conf = mol.beginConformers(); conf != mol.endConformers();
         ++conf) {
      bj::object rjConf;
      addConformer(**conf, rjConf);
      rjConfs.push_back(std::move(rjConf));
    }

    rjMol["conformers"] = std::move(rjConfs);
  }

  bool includePrivate = false, includeComputed = false;
  auto propNames = mol.getPropList(includePrivate, includeComputed);
  if (propNames.size()) {
    bj::object properties;
    addProperties(mol, propNames, properties);
    rjMol["properties"] = std::move(properties);
  }

  bj::object representation;
  representation["name"] = "rdkitRepresentation";
  representation["formatVersion"] = currentRDKitRepresentationVersion;
  representation["toolkitVersion"] = rdkitVersion;

  bool hasArom = false;
  for (const auto &atom : mol.atoms()) {
    if (atom->getIsAromatic()) {
      hasArom = true;
      break;
    }
  }
  if (hasArom) {
    {
      bj::array rjArr;
      for (const auto &atom : mol.atoms()) {
        if (atom->getIsAromatic()) {
          rjArr.push_back(atom->getIdx());
        }
      }
      representation["aromaticAtoms"] = std::move(rjArr);
    }
    {
      bj::array rjArr;
      for (const auto &bond : mol.bonds()) {
        if (bond->getIsAromatic()) {
          rjArr.push_back(bond->getIdx());
        }
      }
      representation["aromaticBonds"] = std::move(rjArr);
    }
  }
  {
    bj::array rjArr;
    if (mol.getAtomWithIdx(0)->hasProp(common_properties::_CIPRank)) {
      for (const auto &atom : mol.atoms()) {
        rjArr.push_back(
            atom->getProp<unsigned int>(common_properties::_CIPRank));
      }
    }
    if (!rjArr.empty()) {
      representation["cipRanks"] = std::move(rjArr);
    }
  }
  {
    bj::array rjArr;
    for (const auto &atom : mol.atoms()) {
      std::string cip;
      if (atom->getPropIfPresent(common_properties::_CIPCode, cip)) {
        bj::array rjElement;
        rjElement.push_back(atom->getIdx());
        rjElement.push_back(cip.c_str());
        rjArr.push_back(rjElement);
      }
    }
    if (!rjArr.empty()) {
      representation["cipCodes"] = std::move(rjArr);
    }
  }
  if (mol.getRingInfo()->numRings()) {
    {
      bj::array rjArr;
      for (const auto &ring : mol.getRingInfo()->atomRings()) {
        bj::array rjRing;
        for (const auto &ai : ring) {
          rjRing.push_back(ai);
        }
        rjArr.push_back(rjRing);
      }
      representation["atomRings"] = std::move(rjArr);
    }
  }

  bj::array rjReprs;
  rjReprs.push_back(representation);

  if (mol.getAtomWithIdx(0)->hasProp(common_properties::_GasteigerCharge)) {
    bj::object chargeRep;
    chargeRep["name"] = "partialCharges";
    chargeRep["generator"] = "RDKit";
    chargeRep["formatVersion"] = currentChargeRepresentationVersion;
    chargeRep["generatorVersion"] = rdkitVersion;

    bj::array rjArr;
    for (const auto &at : mol.atoms()) {
      double rjval;
      if (at->hasProp(common_properties::_GasteigerCharge)) {
        rjval = at->getProp<double>(common_properties::_GasteigerCharge);
      } else {
        rjval = 0.0;
      }
      rjArr.push_back(rjval);
    }
    chargeRep["values"] = std::move(rjArr);
    rjReprs.push_back(std::move(chargeRep));
  }

  if (hasQueryAtoms || hasQueryBonds) {
    bj::object queryRep;
    queryRep["name"] = "rdkitQueries";
    queryRep["formatVersion"] = currentQueryRepresentationVersion;
    queryRep["toolkitVersion"] = rdkitVersion;

    if (hasQueryAtoms) {
      bj::array rjArr;
      for (const auto &atom : mol.atoms()) {
        bj::object rjQ;
        if (atom->hasQuery()) {
          addQuery(*atom->getQuery(), rjQ, params);
        }
        rjArr.push_back(std::move(rjQ));
      }
      queryRep["atomQueries"] = std::move(rjArr);
    }
    if (hasQueryBonds) {
      bj::array rjArr;
      for (const auto &bond : mol.bonds()) {
        bj::object rjQ;
        if (bond->hasQuery()) {
          addQuery(*bond->getQuery(), rjQ, params);
        }
        rjArr.push_back(std::move(rjQ));
      }
      queryRep["bondQueries"] = std::move(rjArr);
    }
    rjReprs.push_back(std::move(queryRep));
  }
  rjMol["extensions"] = std::move(rjReprs);
}
}  // end of anonymous namespace

template <typename T>
std::string MolsToJSONData(const std::vector<T> &mols,
                           const JSONWriteParameters &params) {
  bj::object doc;

  bj::object header;
  initHeader(header, params);
  if (!params.useRDKitExtensions) {
    doc["commonchem"] = header;
  } else {
    doc["rdkitjson"] = header;
  }

  bj::object defaults;

  bj::object atomDefaults;
  initAtomDefaults(atomDefaults);
  defaults["atom"] = atomDefaults;

  bj::object bondDefaults;
  initBondDefaults(bondDefaults);
  defaults["bond"] = bondDefaults;
  doc["defaults"] = defaults;

  bj::array bjMols;
  for (const auto &mol : mols) {
    if (!mol) {
      throw ValueErrorException("null molecule passed to MolsToJSONData");
    }
    bj::object bjMol;
    // write mol;
    addMol(*mol, bjMol, defaults["atom"].as_object(),
           defaults["bond"].as_object(), params);
    bjMols.push_back(std::move(bjMol));
  }
  doc["molecules"] = std::move(bjMols);

  return bj::serialize(doc);
};

template RDKIT_MOLINTERCHANGE_EXPORT std::string MolsToJSONData<ROMol *>(
    const std::vector<ROMol *> &, const JSONWriteParameters &);
template RDKIT_MOLINTERCHANGE_EXPORT std::string MolsToJSONData<RWMol *>(
    const std::vector<RWMol *> &, const JSONWriteParameters &);
template RDKIT_MOLINTERCHANGE_EXPORT std::string MolsToJSONData<const ROMol *>(
    const std::vector<const ROMol *> &, const JSONWriteParameters &);
template RDKIT_MOLINTERCHANGE_EXPORT std::string MolsToJSONData<const RWMol *>(
    const std::vector<const RWMol *> &, const JSONWriteParameters &);

template RDKIT_MOLINTERCHANGE_EXPORT std::string
MolsToJSONData<std::unique_ptr<ROMol>>(
    const std::vector<std::unique_ptr<ROMol>> &, const JSONWriteParameters &);
template RDKIT_MOLINTERCHANGE_EXPORT std::string
MolsToJSONData<std::unique_ptr<RWMol>>(
    const std::vector<std::unique_ptr<RWMol>> &, const JSONWriteParameters &);

}  // end of namespace MolInterchange
}  // end of namespace RDKit
