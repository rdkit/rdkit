//
//  Copyright (C) 2018-2021 Greg Landrum
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

#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/MolInterchange/details.h>
#include <RDGeneral/FileParseException.h>

#include <sstream>
#include <exception>
#include <map>

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>
#include "rapidjson/pointer.h"

namespace rj = rapidjson;

namespace RDKit {

namespace MolInterchange {

namespace {
template <typename T>
void addMol(const T &imol, rj::Value &rjMol, rj::Document &doc,
            const rj::Value &atomDefaults, const rj::Value &bondDefaults);

void initAtomDefaults(rj::Value &rjDefaults, rj::Document &document) {
  // "atomDefaults": {"Z": 6, "impHs": 0, "chg": 0, "stereo": "unspecified",
  // "nrad": 0, "isotope": 0},
  rjDefaults.AddMember("z", 6, document.GetAllocator());
  rjDefaults.AddMember("impHs", 0, document.GetAllocator());
  rjDefaults.AddMember("chg", 0, document.GetAllocator());
  rjDefaults.AddMember("nRad", 0, document.GetAllocator());
  rjDefaults.AddMember("isotope", 0, document.GetAllocator());
  rjDefaults.AddMember("stereo", "unspecified", document.GetAllocator());
}
void initBondDefaults(rj::Value &rjDefaults, rj::Document &document) {
  // "bondDefaults": {"bo": 1, "stereo": "unspecified", "stereoAtoms": []},
  rjDefaults.AddMember("bo", 1, document.GetAllocator());
  rjDefaults.AddMember("stereo", "unspecified", document.GetAllocator());
}
void initHeader(rj::Value &rjHeader, rj::Document &document) {
  rjHeader.AddMember("version", currentMolJSONVersion, document.GetAllocator());
}

void addIntVal(rj::Value &dest, const rj::Value &defaults, const char *tag,
               int val, rj::Document &doc) {
  const auto srt = rj::StringRef(tag);
  const auto &miter = defaults.FindMember(srt);
  if (miter != defaults.MemberEnd()) {
    int dval = miter->value.GetInt();
    if (dval != val) {
      dest.AddMember(srt, val, doc.GetAllocator());
    }
  } else {
    dest.AddMember(srt, val, doc.GetAllocator());
  }
}

void addStringVal(rj::Value &dest, const rj::Value &defaults, const char *tag,
                  const std::string &val, rj::Document &doc) {
  rj::Value nmv;
  nmv.SetString(val.c_str(), val.size(), doc.GetAllocator());
  const auto srt = rj::StringRef(tag);
  const auto &miter = defaults.FindMember(srt);
  if (miter != defaults.MemberEnd()) {
    std::string dval = miter->value.GetString();
    if (val.size() && dval != val) {
      dest.AddMember(srt, nmv, doc.GetAllocator());
    }
  } else {
    dest.AddMember(srt, nmv, doc.GetAllocator());
  }
}

void addAtom(const Atom &atom, rj::Value &rjAtom, rj::Document &doc,
             const rj::Value &rjDefaults) {
  addIntVal(rjAtom, rjDefaults, "z", atom.getAtomicNum(), doc);
  if (!atom.hasQuery()) {
    addIntVal(rjAtom, rjDefaults, "impHs", atom.getTotalNumHs(), doc);
  }
  addIntVal(rjAtom, rjDefaults, "chg", atom.getFormalCharge(), doc);
  addIntVal(rjAtom, rjDefaults, "isotope", atom.getIsotope(), doc);
  addIntVal(rjAtom, rjDefaults, "nRad", atom.getNumRadicalElectrons(), doc);
  std::string chi = "";
  if (inv_chilookup.find(atom.getChiralTag()) != inv_chilookup.end()) {
    chi = inv_chilookup.find(atom.getChiralTag())->second;
  } else {
    BOOST_LOG(rdWarningLog)
        << " unrecognized atom chirality set to default while writing"
        << std::endl;
  }
  addStringVal(rjAtom, rjDefaults, "stereo", chi, doc);
}

template <typename Q>
void addQuery(const Q &query, rj::Value &rjQuery, rj::Document &doc) {
  rj::Value descr;
  descr.SetString(query.getDescription().c_str(), query.getDescription().size(),
                  doc.GetAllocator());
  rjQuery.AddMember("descr", descr, doc.GetAllocator());
  if (!query.getTypeLabel().empty()) {
    rj::Value typ;
    typ.SetString(query.getTypeLabel().c_str(), query.getTypeLabel().size(),
                  doc.GetAllocator());
    rjQuery.AddMember("type", typ, doc.GetAllocator());
  }
  if (query.getNegation()) {
    rjQuery.AddMember("negated", true, doc.GetAllocator());
  }
  if (typeid(query) == typeid(RecursiveStructureQuery)) {
    auto rq = (const RecursiveStructureQuery *)&query;
    auto submol = rq->getQueryMol();
    PRECONDITION(submol, "bad recursive query");
    rj::Value subquery(rj::kObjectType);
    addMol(*submol, subquery, doc,
           *rj::GetValueByPointer(doc, "/defaults/atom"),
           *rj::GetValueByPointer(doc, "/defaults/bond"));
    rjQuery.AddMember("tag", MolPickler::QUERY_RECURSIVE, doc.GetAllocator());
    rjQuery.AddMember("subquery", subquery, doc.GetAllocator());
  } else {
    auto qdetails = PicklerOps::getQueryDetails(&query);
    switch (qdetails.which()) {
      case 0:
        rjQuery.AddMember("tag", boost::get<MolPickler::Tags>(qdetails),
                          doc.GetAllocator());
        break;
      case 1: {
        auto v = boost::get<std::tuple<MolPickler::Tags, int32_t>>(qdetails);
        rjQuery.AddMember("tag", std::get<0>(v), doc.GetAllocator());
        rjQuery.AddMember("val", std::get<1>(v), doc.GetAllocator());
      } break;
      case 2: {
        auto v = boost::get<std::tuple<MolPickler::Tags, int32_t, int32_t>>(
            qdetails);
        rjQuery.AddMember("tag", std::get<0>(v), doc.GetAllocator());
        rjQuery.AddMember("val", std::get<1>(v), doc.GetAllocator());
        if (std::get<2>(v)) {
          rjQuery.AddMember("tol", std::get<2>(v), doc.GetAllocator());
        }
      } break;
      case 3: {
        auto v = boost::get<
            std::tuple<MolPickler::Tags, int32_t, int32_t, int32_t, char>>(
            qdetails);
        rjQuery.AddMember("tag", std::get<0>(v), doc.GetAllocator());
        rjQuery.AddMember("lower", std::get<1>(v), doc.GetAllocator());
        rjQuery.AddMember("upper", std::get<2>(v), doc.GetAllocator());
        if (std::get<3>(v)) {
          rjQuery.AddMember("tol", std::get<3>(v), doc.GetAllocator());
        }
        rjQuery.AddMember("ends", std::get<4>(v), doc.GetAllocator());
      } break;
      case 4: {
        auto v = boost::get<std::tuple<MolPickler::Tags, std::set<int32_t>>>(
            qdetails);
        rjQuery.AddMember("tag", std::get<0>(v), doc.GetAllocator());
        const auto &tset = std::get<1>(v);
        rj::Value sval(rj::kArrayType);
        for (auto val : tset) {
          sval.PushBack(val, doc.GetAllocator());
        }
        rjQuery.AddMember("set", sval, doc.GetAllocator());

      } break;
      default:
        throw MolPicklerException(
            "do not know how to pickle part of the query.");
    }
  }
  if (query.endChildren() != query.beginChildren()) {
    rj::Value children(rj::kArrayType);
    for (auto cit = query.beginChildren(); cit != query.endChildren(); ++cit) {
      rj::Value child(rj::kObjectType);
      addQuery(**cit, child, doc);
      children.PushBack(child, doc.GetAllocator());
    }
    rjQuery.AddMember("children", children, doc.GetAllocator());
  }
}

void addBond(const Bond &bond, rj::Value &rjBond, rj::Document &doc,
             const rj::Value &rjDefaults) {
  int bo = 0;
  if (inv_bolookup.find(bond.getBondType()) != inv_bolookup.end()) {
    bo = inv_bolookup.find(bond.getBondType())->second;
  } else {
    BOOST_LOG(rdWarningLog)
        << " unrecognized bond type set to zero while writing" << std::endl;
  }
  addIntVal(rjBond, rjDefaults, "bo", bo, doc);
  rj::Value rjAtoms(rj::kArrayType);
  rj::Value v1(static_cast<int>(bond.getBeginAtomIdx()));
  rj::Value v2(static_cast<int>(bond.getEndAtomIdx()));
  rjAtoms.PushBack(v1, doc.GetAllocator());
  rjAtoms.PushBack(v2, doc.GetAllocator());
  rjBond.AddMember("atoms", rjAtoms, doc.GetAllocator());

  std::string chi = "";
  if (inv_stereolookup.find(bond.getStereo()) != inv_stereolookup.end()) {
    chi = inv_stereolookup.find(bond.getStereo())->second;
  } else {
    BOOST_LOG(rdWarningLog) << " unrecognized bond stereo " << bond.getStereo()
                            << " set to default while writing" << std::endl;
  }
  addStringVal(rjBond, rjDefaults, "stereo", chi, doc);
  if (chi != "unspecified" && bond.getStereoAtoms().size() == 2) {
    rj::Value rjStereoAtoms(rj::kArrayType);
    rj::Value v1(static_cast<int>(bond.getStereoAtoms()[0]));
    rj::Value v2(static_cast<int>(bond.getStereoAtoms()[1]));
    rjStereoAtoms.PushBack(v1, doc.GetAllocator());
    rjStereoAtoms.PushBack(v2, doc.GetAllocator());
    rjBond.AddMember("stereoAtoms", rjStereoAtoms, doc.GetAllocator());
  }
}

void addConformer(const Conformer &conf, rj::Value &rjConf, rj::Document &doc) {
  int dim = 2;
  if (conf.is3D()) {
    dim = 3;
  }
  rjConf.AddMember("dim", dim, doc.GetAllocator());
  rj::Value rjCoords(rj::kArrayType);
  for (const auto &pos : conf.getPositions()) {
    rj::Value rjPos(rj::kArrayType);
    rjPos.PushBack(pos.x, doc.GetAllocator());
    rjPos.PushBack(pos.y, doc.GetAllocator());
    if (dim == 3) {
      rjPos.PushBack(pos.z, doc.GetAllocator());
    }
    rjCoords.PushBack(rjPos, doc.GetAllocator());
  }
  rjConf.AddMember("coords", rjCoords, doc.GetAllocator());
}

template <typename T>
void addMol(const T &imol, rj::Value &rjMol, rj::Document &doc,
            const rj::Value &atomDefaults, const rj::Value &bondDefaults) {
  RWMol mol(imol);
  if (!mol.getRingInfo()->isInitialized()) {
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
    rj::Value nmv;
    const std::string &nm =
        mol.getProp<std::string>(common_properties::_Name).c_str();
    nmv.SetString(nm.c_str(), nm.size(), doc.GetAllocator());
    rjMol.AddMember("name", nmv, doc.GetAllocator());
  }
  rj::Value rjAtoms(rj::kArrayType);
  bool hasQueryAtoms = false;
  for (const auto &at : mol.atoms()) {
    rj::Value rjAtom(rj::kObjectType);
    addAtom(*at, rjAtom, doc, atomDefaults);
    rjAtoms.PushBack(rjAtom, doc.GetAllocator());
    if (at->hasQuery()) {
      hasQueryAtoms = true;
    }
  }
  rjMol.AddMember("atoms", rjAtoms, doc.GetAllocator());

  rj::Value rjBonds(rj::kArrayType);
  bool hasQueryBonds = false;
  for (const auto &bnd : mol.bonds()) {
    rj::Value rjBond(rj::kObjectType);
    addBond(*bnd, rjBond, doc, bondDefaults);
    rjBonds.PushBack(rjBond, doc.GetAllocator());
    if (bnd->hasQuery()) {
      hasQueryBonds = true;
    }
  }
  rjMol.AddMember("bonds", rjBonds, doc.GetAllocator());

  if (mol.getNumConformers()) {
    rj::Value rjConfs(rj::kArrayType);
    for (auto conf = mol.beginConformers(); conf != mol.endConformers();
         ++conf) {
      rj::Value rjConf(rj::kObjectType);
      addConformer(*(conf->get()), rjConf, doc);
      rjConfs.PushBack(rjConf, doc.GetAllocator());
    }

    rjMol.AddMember("conformers", rjConfs, doc.GetAllocator());
  }

  bool includePrivate = false, includeComputed = false;
  auto propNames = mol.getPropList(includePrivate, includeComputed);
  if (propNames.size()) {
    rj::Value properties(rj::kObjectType);
    for (const auto &pN : propNames) {
      rj::Value rjv;
      try {
        auto val = mol.getProp<int>(pN);
        rjv = val;
      } catch (const boost::bad_any_cast &) {
        try {
          auto val = mol.getProp<double>(pN);
          rjv = val;
        } catch (const boost::bad_any_cast &) {
          try {
            auto val = mol.getProp<std::string>(pN);
            rjv.SetString(val.c_str(), val.size(), doc.GetAllocator());
          } catch (const boost::bad_any_cast &) {
            continue;
          }
        }
      }
      rj::Value rjpN;
      rjpN.SetString(pN.c_str(), pN.size(), doc.GetAllocator());
      properties.AddMember(rjpN, rjv, doc.GetAllocator());
    }
    rjMol.AddMember("properties", properties, doc.GetAllocator());
  }

  rj::Value representation(rj::kObjectType);
  representation.AddMember("name", "rdkitRepresentation", doc.GetAllocator());
  representation.AddMember("formatVersion", currentRDKitRepresentationVersion,
                           doc.GetAllocator());
  rj::Value toolkitVersion;
  toolkitVersion.SetString(rj::StringRef(rdkitVersion));
  representation.AddMember("toolkitVersion", toolkitVersion,
                           doc.GetAllocator());

  bool hasArom = false;
  for (const auto &atom : mol.atoms()) {
    if (atom->getIsAromatic()) {
      hasArom = true;
      break;
    }
  }
  if (hasArom) {
    {
      rj::Value rjArr(rj::kArrayType);
      for (const auto &atom : mol.atoms()) {
        if (atom->getIsAromatic()) {
          rjArr.PushBack(atom->getIdx(), doc.GetAllocator());
        }
      }
      representation.AddMember("aromaticAtoms", rjArr, doc.GetAllocator());
    }
    {
      rj::Value rjArr(rj::kArrayType);
      for (const auto &bond : mol.bonds()) {
        if (bond->getIsAromatic()) {
          rjArr.PushBack(bond->getIdx(), doc.GetAllocator());
        }
      }
      representation.AddMember("aromaticBonds", rjArr, doc.GetAllocator());
    }
  }
  {
    rj::Value rjArr(rj::kArrayType);
    if (mol.getAtomWithIdx(0)->hasProp(common_properties::_CIPRank)) {
      for (const auto &atom : mol.atoms()) {
        rjArr.PushBack(atom->getProp<unsigned int>(common_properties::_CIPRank),
                       doc.GetAllocator());
      }
    }
    if (rjArr.Size()) {
      representation.AddMember("cipRanks", rjArr, doc.GetAllocator());
    }
  }
  {
    rj::Value rjArr(rj::kArrayType);
    for (const auto &atom : mol.atoms()) {
      std::string cip;
      if (atom->getPropIfPresent(common_properties::_CIPCode, cip)) {
        rj::Value cipv;
        cipv.SetString(cip.c_str(), cip.size(), doc.GetAllocator());
        rj::Value rjElement(rj::kArrayType);
        rjElement.PushBack(rj::Value(atom->getIdx()), doc.GetAllocator());
        rjElement.PushBack(cipv, doc.GetAllocator());
        rjArr.PushBack(rjElement, doc.GetAllocator());
      }
    }
    if (rjArr.Size()) {
      representation.AddMember("cipCodes", rjArr, doc.GetAllocator());
    }
  }
  if (mol.getRingInfo()->numRings()) {
    {
      rj::Value rjArr(rj::kArrayType);
      for (const auto &ring : mol.getRingInfo()->atomRings()) {
        rj::Value rjRing(rj::kArrayType);
        for (const auto &ai : ring) {
          rjRing.PushBack(ai, doc.GetAllocator());
        }
        rjArr.PushBack(rjRing, doc.GetAllocator());
      }
      representation.AddMember("atomRings", rjArr, doc.GetAllocator());
    }
  }

  rj::Value rjReprs(rj::kArrayType);
  rjReprs.PushBack(representation, doc.GetAllocator());

  if (mol.getAtomWithIdx(0)->hasProp(common_properties::_GasteigerCharge)) {
    rj::Value representation(rj::kObjectType);
    representation.AddMember("name", "partialCharges", doc.GetAllocator());
    representation.AddMember("generator", "RDKit", doc.GetAllocator());
    representation.AddMember("formatVersion",
                             currentChargeRepresentationVersion,
                             doc.GetAllocator());
    rj::Value toolkitVersion;
    toolkitVersion.SetString(rj::StringRef(rdkitVersion));
    representation.AddMember("generatorVersion", toolkitVersion,
                             doc.GetAllocator());

    rj::Value rjArr(rj::kArrayType);
    for (const auto &at : mol.atoms()) {
      rj::Value rjval;
      if (at->hasProp(common_properties::_GasteigerCharge)) {
        rjval = at->getProp<double>(common_properties::_GasteigerCharge);
      } else {
        rjval = 0.0;
      }
      rjArr.PushBack(rjval, doc.GetAllocator());
    }
    representation.AddMember("values", rjArr, doc.GetAllocator());
    rjReprs.PushBack(representation, doc.GetAllocator());
  }

  if (hasQueryAtoms || hasQueryBonds) {
    rj::Value representation(rj::kObjectType);
    representation.AddMember("name", "rdkitQueries", doc.GetAllocator());
    representation.AddMember("formatVersion", currentQueryRepresentationVersion,
                             doc.GetAllocator());
    rj::Value toolkitVersion;
    toolkitVersion.SetString(rj::StringRef(rdkitVersion));
    representation.AddMember("toolkitVersion", toolkitVersion,
                             doc.GetAllocator());

    if (hasQueryAtoms) {
      rj::Value rjArr(rj::kArrayType);
      for (const auto &atom : mol.atoms()) {
        rj::Value rjQ(rj::kObjectType);
        if (atom->hasQuery()) {
          addQuery(*atom->getQuery(), rjQ, doc);
        }
        rjArr.PushBack(rjQ, doc.GetAllocator());
      }
      representation.AddMember("atomQueries", rjArr, doc.GetAllocator());
    }
    if (hasQueryBonds) {
      rj::Value rjArr(rj::kArrayType);
      for (const auto &bond : mol.bonds()) {
        rj::Value rjQ(rj::kObjectType);
        if (bond->hasQuery()) {
          addQuery(*bond->getQuery(), rjQ, doc);
        }
        rjArr.PushBack(rjQ, doc.GetAllocator());
      }
      representation.AddMember("bondQueries", rjArr, doc.GetAllocator());
    }
    rjReprs.PushBack(representation, doc.GetAllocator());
  }
  rjMol.AddMember("extensions", rjReprs, doc.GetAllocator());
}
}  // end of anonymous namespace

template <typename T>
std::string MolsToJSONData(const std::vector<T> &mols) {
  std::string res = "";
  rj::Document doc;
  doc.SetObject();

  rj::Value header(rj::kObjectType);
  initHeader(header, doc);
  doc.AddMember("commonchem", header, doc.GetAllocator());

  rj::Value defaults(rj::kObjectType);

  rj::Value atomDefaults(rj::kObjectType);
  initAtomDefaults(atomDefaults, doc);
  defaults.AddMember("atom", atomDefaults, doc.GetAllocator());

  rj::Value bondDefaults(rj::kObjectType);
  initBondDefaults(bondDefaults, doc);
  defaults.AddMember("bond", bondDefaults, doc.GetAllocator());
  doc.AddMember("defaults", defaults, doc.GetAllocator());

  rj::Value rjMols(rj::kArrayType);
  for (const auto &mol : mols) {
    rj::Value rjMol(rj::kObjectType);
    // write mol;
    addMol(*mol, rjMol, doc, *rj::GetValueByPointer(doc, "/defaults/atom"),
           *rj::GetValueByPointer(doc, "/defaults/bond"));
    rjMols.PushBack(rjMol, doc.GetAllocator());
  }
  doc.AddMember("molecules", rjMols, doc.GetAllocator());

  rj::StringBuffer buffer;
  rj::Writer<rj::StringBuffer> writer(buffer);
  writer.SetMaxDecimalPlaces(4);
  doc.Accept(writer);
  return buffer.GetString();
};

template RDKIT_MOLINTERCHANGE_EXPORT std::string MolsToJSONData<ROMol *>(
    const std::vector<ROMol *> &);
template RDKIT_MOLINTERCHANGE_EXPORT std::string MolsToJSONData<RWMol *>(
    const std::vector<RWMol *> &);
template RDKIT_MOLINTERCHANGE_EXPORT std::string MolsToJSONData<const ROMol *>(
    const std::vector<const ROMol *> &);
template RDKIT_MOLINTERCHANGE_EXPORT std::string MolsToJSONData<const RWMol *>(
    const std::vector<const RWMol *> &);

}  // end of namespace MolInterchange
}  // end of namespace RDKit
