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
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/StereoGroup.h>

#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/MolInterchange/details.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/MolPickler.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/format.hpp>
#include <RDGeneral/BoostEndInclude.h>
using namespace Queries;

#include <sstream>
#include <exception>
#include <map>

#include <RDGeneral/BoostStartInclude.h>
// make sure we're using boost::json header-only
#define BOOST_JSON_NO_LIB
#define BOOST_CONTAINER_NO_LIB
#include <boost/json.hpp>
#include <boost/json/src.hpp>  // only include this once in the project!
#include <RDGeneral/BoostEndInclude.h>
namespace bj = boost::json;

namespace RDKit {

namespace MolInterchange {

namespace {
struct DefaultValueCache {
  DefaultValueCache(const bj::value &defs) : bjDefaults(defs) {};
  const bj::value &bjDefaults;
  mutable std::map<const char *, int> intMap;
  mutable std::map<const char *, bool> boolMap;
  mutable std::map<const char *, std::string> stringMap;
  int getInt(const char *key) const {
    PRECONDITION(key, "no key");
    const auto &lookup = intMap.find(key);
    if (lookup != intMap.end()) {
      return lookup->second;
    }
    if (const auto fobj = bjDefaults.if_object()) {
      if (const auto kit = fobj->find(key); kit != fobj->end()) {
        const auto val = kit->value().if_int64();
        if (!val) {
          throw FileParseException(std::string("Bad format: value of ") +
                                   std::string(key) +
                                   std::string(" is not an int"));
        }
        auto res = static_cast<int>(*val);
        intMap[key] = res;
        return res;
      }
    }
    return 0;
  }
  bool getBool(const char *key) const {
    PRECONDITION(key, "no key");
    const auto &lookup = boolMap.find(key);
    if (lookup != boolMap.end()) {
      return lookup->second;
    }
    if (const auto fobj = bjDefaults.if_object()) {
      if (const auto kit = fobj->find(key); kit != fobj->end()) {
        const auto val = kit->value().if_bool();
        if (!val) {
          throw FileParseException(std::string("Bad format: value of ") +
                                   std::string(key) +
                                   std::string(" is not a bool"));
        }
        bool res = *val;
        boolMap[key] = res;
        return res;
      }
    }
    return false;
  }

  std::string getString(const char *key) const {
    PRECONDITION(key, "no key");
    const auto &lookup = stringMap.find(key);
    if (lookup != stringMap.end()) {
      return lookup->second;
    }
    if (const auto fobj = bjDefaults.if_object()) {
      if (const auto kit = fobj->find(key); kit != fobj->end()) {
        const auto val = kit->value().if_string();
        if (!val) {
          throw FileParseException(std::string("Bad format: value of ") +
                                   std::string(key) +
                                   std::string(" is not a string"));
        }
        auto res = *val;
        stringMap[key] = res;
        return res;
      }
    }
    return "";
  }
};

int getIntDefaultValue(const char *key, const bj::value &from,
                       const DefaultValueCache &defaults) {
  PRECONDITION(key, "no key");
  if (const auto fobj = from.if_object()) {
    if (const auto kit = fobj->find(key); kit != fobj->end()) {
      const auto val = kit->value().if_int64();
      if (!val) {
        throw FileParseException(std::string("Bad format: value of ") +
                                 std::string(key) +
                                 std::string(" is not an int"));
      }
      return static_cast<int>(*val);
    }
  }
  return defaults.getInt(key);
}
bool getBoolDefaultValue(const char *key, const bj::value &from,
                         const DefaultValueCache &defaults) {
  PRECONDITION(key, "no key");
  if (const auto fobj = from.if_object()) {
    if (const auto kit = fobj->find(key); kit != fobj->end()) {
      const auto val = kit->value().if_bool();
      if (!val) {
        throw FileParseException(std::string("Bad format: value of ") +
                                 std::string(key) +
                                 std::string(" is not a bool"));
      }
      return *val;
    }
  }
  return defaults.getBool(key);
}
std::string getStringDefaultValue(const char *key, const bj::value &from,
                                  const DefaultValueCache &defaults) {
  PRECONDITION(key, "no key");
  if (const auto fobj = from.if_object()) {
    if (const auto kit = fobj->find(key); kit != fobj->end()) {
      const auto val = kit->value().if_string();
      if (!val) {
        throw FileParseException(std::string("Bad format: value of ") +
                                 std::string(key) +
                                 std::string(" is not a string"));
      }
      return *val;
    }
  }
  return defaults.getString(key);
}

void readAtom(RWMol *mol, const bj::value &atomVal,
              const DefaultValueCache &atomDefaults,
              const JSONParseParameters &params) {
  PRECONDITION(mol, "no mol");
  std::string stereo = getStringDefaultValue("stereo", atomVal, atomDefaults);
  auto stereoVal = chilookup.find(stereo);
  if (stereoVal == chilookup.end()) {
    throw FileParseException("Bad Format: bad stereo value for atom");
  }
  std::unique_ptr<Atom> at(new Atom(getIntDefaultValue("z", atomVal, atomDefaults)));
  if (params.useHCounts) {
    at->setNoImplicit(true);
    at->setNumExplicitHs(getIntDefaultValue("impHs", atomVal, atomDefaults));
  }
  at->setFormalCharge(getIntDefaultValue("chg", atomVal, atomDefaults));
  at->setNumRadicalElectrons(getIntDefaultValue("nRad", atomVal, atomDefaults));
  at->setIsotope(getIntDefaultValue("isotope", atomVal, atomDefaults));
  at->setChiralTag(stereoVal->second);
  bool updateLabel = false, takeOwnership = true;
  mol->addAtom(at.release(), updateLabel, takeOwnership);
}

void readBond(RWMol *mol, const bj::value &bondVal,
              const DefaultValueCache &bondDefaults, bool &needStereoLoop) {
  PRECONDITION(mol, "no mol");
  unsigned int bo = getIntDefaultValue("bo", bondVal, bondDefaults);
  auto bondOrder = bolookup.find(bo);
  if (bondOrder == bolookup.end()) {
    throw FileParseException("Bad Format: bad bond order for bond");
  }
  const auto &aids = bondVal.at("atoms").as_array();
  std:::unique_ptr<Bond> bnd(new Bond());
  bnd->setBeginAtomIdx(static_cast<int>(aids.at(0).as_int64()));
  bnd->setEndAtomIdx(static_cast<int>(aids.at(1).as_int64()));
  bnd->setBondType(bondOrder->second);
  mol->addBond(bnd.release());
  std::string stereo = getStringDefaultValue("stereo", bondVal, bondDefaults);
  if (stereo != "unspecified") {
    needStereoLoop = true;
  }
}

template <typename T>
void parseProperties(T &obj, const bj::value &propsVal) {
  for (const auto &propVal : propsVal.as_object()) {
    if (propVal.value().is_int64()) {
      obj.setProp(propVal.key(), static_cast<int>(propVal.value().as_int64()));
    } else if (propVal.value().is_double()) {
      obj.setProp(propVal.key(), propVal.value().as_double());
    } else if (propVal.value().is_string()) {
      obj.setProp(propVal.key(), propVal.value().as_string().data());
    }
  }
}

void readStereoGroups(RWMol *mol, const bj::value &sgVals) {
  PRECONDITION(mol, "no mol");

  std::vector<StereoGroup> molSGs(mol->getStereoGroups());
  for (const auto &sgVal : sgVals.as_array()) {
    if (!sgVal.as_object().contains("type")) {
      throw FileParseException("Bad Format: stereogroup does not have a type");
    }
    if (!sgVal.as_object().contains("atoms") &&
        !sgVal.as_object().contains("bonds")) {
      throw FileParseException(
          "Bad Format: stereogroup does not have either atoms or bonds");
    }
    if (MolInterchange::stereoGrouplookup.find(
            sgVal.at("type").as_string()) ==
        MolInterchange::stereoGrouplookup.end()) {
      throw FileParseException("Bad Format: bad stereogroup type");
    }
    const auto typ = MolInterchange::stereoGrouplookup.at(
        sgVal.at("type").as_string());

    unsigned gId = 0;
    if (typ != StereoGroupType::STEREO_ABSOLUTE &&
        sgVal.as_object().contains("id")) {
      gId = static_cast<unsigned>(sgVal.at("id").as_int64());
    }

    std::vector<Atom *> atoms;
    std::vector<Bond *> bonds;
    if (sgVal.as_object().contains("atoms")) {
      const auto &aids = sgVal.at("atoms").as_array();
      for (const auto &aid : aids) {
        atoms.push_back(
            mol->getAtomWithIdx(static_cast<unsigned>(aid.as_int64())));
      }
    }
    if (sgVal.as_object().contains("bonds")) {
      const auto &bids = sgVal.at("bonds").as_array();
      for (const auto &bid : bids) {
        bonds.push_back(
            mol->getBondWithIdx(static_cast<unsigned>(bid.as_int64())));
      }
    }

    if (!atoms.empty() || !bonds.empty()) {
      molSGs.emplace_back(typ, std::move(atoms), std::move(bonds), gId);
    }
  }
  mol->setStereoGroups(std::move(molSGs));
}

void readSubstanceGroups(RWMol *mol, const bj::value &sgVals) {
  PRECONDITION(mol, "no mol");

  for (const auto &sgVal : sgVals.as_array()) {
    if (!sgVal.as_object().contains("properties") ||
        !sgVal.at("properties").as_object().contains("TYPE")) {
      throw FileParseException(
          "Bad Format: substance group does not have TYPE property");
    }

    auto sgType = sgVal.at("properties").at("TYPE").as_string();
    if (!SubstanceGroupChecks::isValidType(sgType)) {
      throw FileParseException(
          (boost::format(
               "Bad Format: substance group TYPE '%s' not recognized") %
           sgType)
              .str());
    }
    SubstanceGroup sg(mol, sgType);

    parseProperties(sg, sgVal.at("properties"));
    std::string pval;
    if (sg.getPropIfPresent("SUBTYPE", pval) &&
        !SubstanceGroupChecks::isValidSubType(pval)) {
      throw FileParseException(
          (boost::format(
               "Bad Format: substance group SUBTYPE '%s' not recognized") %
           pval)
              .str());
    }
    if (sg.getPropIfPresent("CONNECT", pval) &&
        !SubstanceGroupChecks::isValidConnectType(pval)) {
      throw FileParseException(
          (boost::format(
               "Bad Format: substance group CONNECT type '%s' not recognized") %
           pval)
              .str());
    }

    if (sgVal.as_object().contains("atoms")) {
      const auto &aids = sgVal.at("atoms").as_array();
      std::vector<unsigned int> atoms;
      for (const auto &aid : aids) {
        atoms.push_back(static_cast<unsigned>(aid.as_int64()));
      }
      sg.setAtoms(atoms);
    }

    if (sgVal.as_object().contains("bonds")) {
      const auto &aids = sgVal.at("bonds").as_array();
      std::vector<unsigned int> bonds;
      for (const auto &aid : aids) {
        bonds.push_back(static_cast<unsigned>(aid.as_int64()));
      }
      sg.setBonds(bonds);
    }

    if (sgVal.as_object().contains("parentAtoms")) {
      const auto &aids = sgVal.at("parentAtoms").as_array();
      std::vector<unsigned int> atoms;
      for (const auto &aid : aids) {
        atoms.push_back(static_cast<unsigned>(aid.as_int64()));
      }
      sg.setParentAtoms(atoms);
    }

    if (sgVal.as_object().contains("brackets")) {
      const auto &brks = sgVal.at("brackets").as_array();
      for (const auto &brk : brks) {
        SubstanceGroup::Bracket bracket;
        unsigned int idx = 0;
        for (const auto &pt : brk.as_array()) {
          const auto &pta = pt.as_array();
          if (pta.size() != 3) {
            throw FileParseException(
                "Bad Format: bracket point doesn't have three coordinates");
          }
          RDGeom::Point3D loc(pta[0].as_double(), pta[1].as_double(),
                              pta[2].as_double());
          bracket[idx++] = std::move(loc);
        }
        sg.getBrackets().push_back(std::move(bracket));
      }
    }

    if (sgVal.as_object().contains("cstates")) {
      const auto &cstats = sgVal.at("cstates").as_array();
      for (const auto &cstat : cstats) {
        SubstanceGroup::CState cstate;
        cstate.bondIdx = static_cast<unsigned>(cstat.at("bond").as_int64());
        if (cstat.as_object().contains("vector")) {
          const auto &pta = cstat.at("vector").as_array();
          if (pta.size() != 3) {
            throw FileParseException(
                "Bad Format: cstate vector doesn't have three coordinates");
          }
          RDGeom::Point3D loc(pta[0].as_double(), pta[1].as_double(),
                              pta[2].as_double());
          cstate.vector = std::move(loc);
        }
        sg.getCStates().push_back(std::move(cstate));
      }
    }

    if (sgVal.as_object().contains("attachPoints")) {
      const auto &aps = sgVal.at("attachPoints").as_array();
      for (const auto &ap : aps) {
        SubstanceGroup::AttachPoint attach;
        attach.aIdx = static_cast<unsigned>(ap.at("aIdx").as_int64());
        if (ap.as_object().contains("lvIdx")) {
          attach.lvIdx = static_cast<unsigned>(ap.at("lvIdx").as_int64());
        }
        if (ap.as_object().contains("id")) {
          attach.id = ap.at("id").as_string();
        }
        sg.getAttachPoints().push_back(std::move(attach));
      }
    }
    addSubstanceGroup(*mol, sg);
  }
}

void readBondStereo(Bond *bnd, const bj::value &bondVal,
                    const DefaultValueCache &bondDefaults) {
  PRECONDITION(bnd, "no bond");

  std::string stereo = getStringDefaultValue("stereo", bondVal, bondDefaults);
  if (stereo == "unspecified") {
    return;
  }
  if (stereoBondlookup.find(stereo) == stereoBondlookup.end()) {
    throw FileParseException("Bad Format: bond stereo value for bond");
  }
  if (bondVal.as_object().contains("stereoAtoms")) {
    const auto &aids = bondVal.at("stereoAtoms").as_array();
    bnd->setStereoAtoms(static_cast<int>(aids[0].as_int64()),
                        static_cast<int>(aids[1].as_int64()));
  } else if (stereo == "cis" || stereo == "trans") {
    throw FileParseException(
        "Bad Format: bond stereo provided without stereoAtoms");
  }
  bnd->setStereo(stereoBondlookup.find(stereo)->second);
}  // namespace

void readConformer(Conformer *conf, const bj::value &confVal) {
  PRECONDITION(conf, "no conformer");

  if (!confVal.as_object().contains("dim")) {
    throw FileParseException("Bad Format: no conformer dimension");
  }
  size_t dim = static_cast<size_t>(confVal.at("dim").as_int64());
  if (dim == 2) {
    conf->set3D(false);
  } else if (dim == 3) {
    conf->set3D(true);
  } else {
    throw FileParseException("Bad Format: conformer dimension != 2 or 3");
  }
  if (!confVal.as_object().contains("coords")) {
    throw FileParseException("Bad Format: no conformer coords");
  }
  size_t idx = 0;
  for (const auto &ptVal : confVal.at("coords").as_array()) {
    const auto &arr = ptVal.as_array();
    if (arr.size() != dim) {
      throw FileParseException("coordinate contains wrong number of values");
    }
    RDGeom::Point3D pt(arr[0].as_double(), arr[1].as_double(),
                       (dim == 3 ? arr[2].as_double() : 0.0));
    conf->setAtomPos(idx++, pt);
  }
  if (idx != conf->getNumAtoms()) {
    throw FileParseException(
        "Bad Format: conformer doesn't contain coordinates for all atoms");
  }
}

void readPartialCharges(RWMol *mol, const bj::value &repVal,
                        const JSONParseParameters &) {
  PRECONDITION(mol, "no molecule");
  PRECONDITION(repVal.at("name").as_string() == std::string("partialCharges"),
               "bad charges");
  if (!repVal.as_object().contains("formatVersion")) {
    throw FileParseException("Bad Format: missing version");
  }
  if (static_cast<int>(repVal.at("formatVersion").as_int64()) >
      currentChargeRepresentationVersion) {
    BOOST_LOG(rdWarningLog)
        << "partialCharges version "
        << static_cast<int>(repVal.at("formatVersion").as_int64())
        << " too recent. Ignoring it." << std::endl;
    return;
  }
  {
    if (repVal.as_object().contains("values")) {
      const auto &values = repVal.at("values").as_array();
      if (values.size() != mol->getNumAtoms()) {
        throw FileParseException(
            "Bad Format: size of values array != num atoms");
      }
      for (unsigned int idx = 0; idx != mol->getNumAtoms(); ++idx) {
        const auto &val = values[idx];
        if (!val.is_double()) {
          throw FileParseException("Bad Format: partial charge not double");
        }
        mol->getAtomWithIdx(idx)->setProp(common_properties::_GasteigerCharge,
                                          val.as_double());
      }
    }
  }
}
void processMol(RWMol *mol, const bj::value &molval,
                const DefaultValueCache &atomDefaults,
                const DefaultValueCache &bondDefaults,
                const JSONParseParameters &params);
Query<int, Atom const *, true> *readQuery(Atom const *owner,
                                          const bj::value &repVal,
                                          const DefaultValueCache &atomDefaults,
                                          const DefaultValueCache &bondDefaults,
                                          const JSONParseParameters &params);
Query<int, Bond const *, true> *readQuery(Bond const *owner,
                                          const bj::value &repVal,
                                          const DefaultValueCache &atomDefaults,
                                          const DefaultValueCache &bondDefaults,
                                          const JSONParseParameters &params);

template <class T>
Query<int, T const *, true> *readBaseQuery(T const *owner,
                                           const bj::value &repVal,
                                           const JSONParseParameters &) {
  PRECONDITION(owner, "no query");
  PRECONDITION(repVal.as_object().contains("tag"), "no tag");
  int tag = repVal.at("tag").as_int64();
  if (!repVal.as_object().contains("descr")) {
    throw FileParseException("Bad Format: missing query description");
  }
  Query<int, T const *, true> *res = nullptr;
  switch (tag) {
    case MolPickler::QUERY_AND:
      res = new AndQuery<int, T const *, true>();
      break;
    case MolPickler::QUERY_OR:
      res = new OrQuery<int, T const *, true>();
      break;
    case MolPickler::QUERY_XOR:
      res = new XOrQuery<int, T const *, true>();
      break;
    case MolPickler::QUERY_EQUALS:
      res = new EqualityQuery<int, T const *, true>();
      static_cast<EqualityQuery<int, T const *, true> *>(res)->setVal(
          repVal.at("val").as_int64());
      if (repVal.as_object().contains("tol")) {
        static_cast<EqualityQuery<int, T const *, true> *>(res)->setTol(
            repVal.at("tol").as_int64());
      }
      break;
    case MolPickler::QUERY_GREATER:
      res = new GreaterQuery<int, T const *, true>();
      static_cast<GreaterQuery<int, T const *, true> *>(res)->setVal(
          repVal.at("val").as_int64());
      if (repVal.as_object().contains("tol")) {
        static_cast<GreaterQuery<int, T const *, true> *>(res)->setTol(
            repVal.at("tol").as_int64());
      }
      break;
    case MolPickler::QUERY_GREATEREQUAL:
      res = new GreaterEqualQuery<int, T const *, true>();
      static_cast<GreaterEqualQuery<int, T const *, true> *>(res)->setVal(
          repVal.at("val").as_int64());
      if (repVal.as_object().contains("tol")) {
        static_cast<GreaterEqualQuery<int, T const *, true> *>(res)->setTol(
            repVal.at("tol").as_int64());
      }
      break;
    case MolPickler::QUERY_LESS:
      res = new LessQuery<int, T const *, true>();
      static_cast<LessQuery<int, T const *, true> *>(res)->setVal(
          repVal.at("val").as_int64());
      if (repVal.as_object().contains("tol")) {
        static_cast<LessQuery<int, T const *, true> *>(res)->setTol(
            repVal.at("tol").as_int64());
      }
      break;
    case MolPickler::QUERY_LESSEQUAL:
      res = new LessEqualQuery<int, T const *, true>();
      static_cast<LessEqualQuery<int, T const *, true> *>(res)->setVal(
          repVal.at("val").as_int64());
      if (repVal.as_object().contains("tol")) {
        static_cast<LessEqualQuery<int, T const *, true> *>(res)->setTol(
            repVal.at("tol").as_int64());
      }
      break;
    case MolPickler::QUERY_NULL:
      res = new Query<int, T const *, true>();
      break;
    case MolPickler::QUERY_RANGE:
      res = new RangeQuery<int, T const *, true>();
      static_cast<RangeQuery<int, T const *, true> *>(res)->setLower(
          repVal.at("lower").as_int64());
      static_cast<RangeQuery<int, T const *, true> *>(res)->setUpper(
          repVal.at("upper").as_int64());
      if (repVal.as_object().contains("tol")) {
        static_cast<RangeQuery<int, T const *, true> *>(res)->setTol(
            repVal.at("tol").as_int64());
      }
      if (repVal.as_object().contains("ends")) {
        short ends = repVal.at("ends").as_int64();
        const unsigned int lowerOpen = 1 << 1;
        const unsigned int upperOpen = 1;
        static_cast<RangeQuery<int, T const *, true> *>(res)->setEndsOpen(
            ends & lowerOpen, ends & upperOpen);
      }
      break;
    case MolPickler::QUERY_SET:
      res = new SetQuery<int, T const *, true>();
      if (repVal.as_object().contains("set")) {
        for (const auto &member : repVal.at("set").as_array()) {
          static_cast<SetQuery<int, T const *, true> *>(res)->insert(
              member.as_int64());
        }
      }
      break;
    default:
      throw FileParseException(
          (boost::format("Bad Format: unknown query tag %s") % tag).str());
  }

  return res;
}

template <class T, class U>
void finishQuery(T const *owner, U *res, const bj::value &repVal,
                 const DefaultValueCache &atomDefaults,
                 const DefaultValueCache &bondDefaults,
                 const JSONParseParameters &params) {
  PRECONDITION(owner, "no owner");
  PRECONDITION(res, "no result");
  auto descr = repVal.at("descr").as_string();
  res->setDescription(descr);
  std::string typ;
  if (repVal.as_object().contains("type")) {
    typ = repVal.at("type").as_string();
  }
  if (!typ.empty()) {
    res->setTypeLabel(typ);
  }
  bool negated = false;
  if (repVal.as_object().contains("negated")) {
    negated = repVal.at("negated").as_bool();
  }
  res->setNegation(negated);
  QueryOps::finalizeQueryFromDescription(res, owner);

  if (repVal.as_object().contains("children")) {
    for (const auto &child : repVal.at("children").as_array()) {
      typename U::CHILD_TYPE childq{
          readQuery(owner, child, atomDefaults, bondDefaults, params)};
      res->addChild(childq);
    }
  }
}

Query<int, Atom const *, true> *readQuery(Atom const *owner,
                                          const bj::value &repVal,
                                          const DefaultValueCache &atomDefaults,
                                          const DefaultValueCache &bondDefaults,
                                          const JSONParseParameters &params) {
  PRECONDITION(owner, "no owner");
  if (!repVal.as_object().contains("tag")) {
    throw FileParseException("Bad Format: missing atom query tag");
  }
  Query<int, Atom const *, true> *res = nullptr;
  int tag = static_cast<int>(repVal.at("tag").as_int64());
  if (tag == MolPickler::QUERY_RECURSIVE) {
    if (!repVal.as_object().contains("subquery")) {
      throw FileParseException("Bad Format: missing subquery");
    }
    auto *mol = new RWMol();
    processMol(mol, repVal.at("subquery"), atomDefaults, bondDefaults, params);
    res = new RecursiveStructureQuery(mol);
  } else if (tag == MolPickler::QUERY_ATOMRING) {
    res = new AtomRingQuery();
    static_cast<EqualityQuery<int, Atom const *, true> *>(res)->setVal(
        static_cast<int>(repVal.at("val").as_int64()));
    if (repVal.as_object().contains("tol")) {
      static_cast<EqualityQuery<int, Atom const *, true> *>(res)->setTol(
          static_cast<int>(repVal.at("tol").as_int64()));
    }
  } else {
    res = readBaseQuery(owner, repVal, params);
  }
  if (res) {
    finishQuery(owner, res, repVal, atomDefaults, bondDefaults, params);
  }
  return res;
}
Query<int, Bond const *, true> *readQuery(Bond const *bond,
                                          const bj::value &repVal,
                                          const DefaultValueCache &atomDefaults,
                                          const DefaultValueCache &bondDefaults,
                                          const JSONParseParameters &params) {
  PRECONDITION(bond, "no owner");
  if (!repVal.as_object().contains("tag")) {
    throw FileParseException("Bad Format: missing bond query tag");
  }
  Query<int, Bond const *, true> *res = nullptr;
  res = readBaseQuery(bond, repVal, params);
  if (res) {
    finishQuery(bond, res, repVal, atomDefaults, bondDefaults, params);
  }

  return res;
}

void readQueries(RWMol *mol, const bj::value &repVal,
                 const DefaultValueCache &atomDefaults,
                 const DefaultValueCache &bondDefaults,
                 const JSONParseParameters &params) {
  PRECONDITION(mol, "no molecule");
  PRECONDITION(
      repVal.at("name").as_string() == std::string("rdkitQueries"),
      "bad queries");
  if (!repVal.as_object().contains("formatVersion")) {
    throw FileParseException("Bad Format: missing format_version");
  }
  if (repVal.at("formatVersion").as_int64() >
      currentQueryRepresentationVersion) {
    BOOST_LOG(rdWarningLog) << "RDKit query representation format version "
                            << repVal.at("formatVersion").as_int64()
                            << " too recent. Ignoring it." << std::endl;
    return;
  }
  {
    const auto &miter = repVal.as_object().find("atomQueries");
    if (miter != repVal.as_object().end()) {
      size_t idx = 0;
      for (const auto &val : miter->value().as_array()) {
        if (!val.is_object()) {
          throw FileParseException("Bad Format: atomQuery not object");
        }
        if (!val.as_object().contains("tag")) {
          // nothing here, continue
          continue;
        }
        if (idx >= mol->getNumAtoms()) {
          throw FileParseException("too much atom data found");
        }
        auto atom = mol->getAtomWithIdx(idx);
        CHECK_INVARIANT(atom != nullptr, "no atom");
        // we need to replace the current atom with a query atom:
        QueryAtom qatom(*atom);
        // that copy created a bunch of query info by default,
        // but we want to get the info from the JSON, so delete
        // that:
        qatom.setQuery(nullptr);
        mol->replaceAtom(idx, &qatom);
        atom = mol->getAtomWithIdx(idx);
        static_cast<QueryAtom *>(atom)->setQuery(
            readQuery(atom, val, atomDefaults, bondDefaults, params));
        ++idx;
      }
    }
  }
  {
    const auto &miter = repVal.as_object().find("bondQueries");
    if (miter != repVal.as_object().end()) {
      size_t idx = 0;
      for (const auto &val : miter->value().as_array()) {
        if (!val.is_object()) {
          throw FileParseException("Bad Format: bondQuery not object");
        }
        if (!val.as_object().contains("tag")) {
          // nothing here, continue
          continue;
        }
        if (idx >= mol->getNumBonds()) {
          throw FileParseException("too much bond data found");
        }
        auto bond = mol->getBondWithIdx(idx);
        CHECK_INVARIANT(bond != nullptr, "no bond");
        QueryBond qbond(*bond);
        qbond.setQuery(nullptr);
        mol->replaceBond(idx, &qbond);
        bond = mol->getBondWithIdx(idx);
        static_cast<QueryBond *>(bond)->setQuery(
            readQuery(bond, val, atomDefaults, bondDefaults, params));
        ++idx;
      }
    }
  }
}

void readRDKitRepresentation(RWMol *mol, const bj::value &repVal,
                             const JSONParseParameters &params) {
  PRECONDITION(mol, "no molecule");
  PRECONDITION(repVal.at("name").as_string() ==
                   std::string("rdkitRepresentation"),
               "bad representation");
  if (!repVal.as_object().contains("formatVersion")) {
    throw FileParseException("Bad Format: missing format_version");
  }
  if (repVal.at("formatVersion").as_int64() >
      currentRDKitRepresentationVersion) {
    BOOST_LOG(rdWarningLog) << "RDKit representation format version "
                            << repVal.at("formatVersion").as_int64()
                            << " too recent. Ignoring it." << std::endl;
    return;
  }
  {
    const auto &miter = repVal.as_object().find("aromaticAtoms");
    if (miter != repVal.as_object().end()) {
      for (const auto &val : miter->value().as_array()) {
        if (!val.is_int64()) {
          throw FileParseException("Bad Format: aromaticAtom not int");
        }
        mol->getAtomWithIdx(val.as_int64())->setIsAromatic(true);
      }
    }
  }
  {
    const auto &miter = repVal.as_object().find("aromaticBonds");
    if (miter != repVal.as_object().end()) {
      for (const auto &val : miter->value().as_array()) {
        if (!val.is_int64()) {
          throw FileParseException("Bad Format: aromaticBond not int");
        }
        mol->getBondWithIdx(val.as_int64())->setIsAromatic(true);
        if (params.setAromaticBonds) {
          mol->getBondWithIdx(val.as_int64())->setBondType(Bond::AROMATIC);
        }
      }
    }
  }
  {
    const auto &miter = repVal.as_object().find("cipRanks");
    if (miter != repVal.as_object().end()) {
      size_t i = 0;
      for (const auto &val : miter->value().as_array()) {
        if (!val.is_int64()) {
          throw FileParseException("Bad Format: ciprank not int");
        }
        mol->getAtomWithIdx(i++)->setProp(
            common_properties::_CIPRank,
            static_cast<unsigned int>(val.as_int64()));
      }
    }
  }
  {
    const auto &miter = repVal.as_object().find("cipCodes");
    if (miter != repVal.as_object().end()) {
      for (const auto &val : miter->value().as_array()) {
        if (!val.is_array()) {
          throw FileParseException("Bad Format: CIPCode not string");
        }
        mol->getAtomWithIdx(val.at(0).as_int64())
            ->setProp(common_properties::_CIPCode,
                      val.at(1).as_string());
      }
    }
  }
  {
    const auto &miter = repVal.as_object().find("atomRings");
    if (miter != repVal.as_object().end()) {
      CHECK_INVARIANT(!mol->getRingInfo()->isInitialized(),
                      "rings already initialized");
      auto ri = mol->getRingInfo();
      ri->initialize();
      for (const auto &val : miter->value().as_array()) {
        if (!val.is_array()) {
          throw FileParseException("Bad Format: atomRing not array");
        }
        INT_VECT atomRing;
        INT_VECT bondRing;
        size_t sz = val.as_array().size();
        atomRing.reserve(sz);
        bondRing.reserve(sz);
        for (size_t i = 0; i < sz - 1; ++i) {
          int idx1 = static_cast<int>(val.as_array()[i].as_int64());
          int idx2 = static_cast<int>(val.as_array()[i + 1].as_int64());
          atomRing.push_back(idx1);
          const auto &bnd = mol->getBondBetweenAtoms(idx1, idx2);
          CHECK_INVARIANT(bnd, "no bond found for ring");
          bondRing.push_back(bnd->getIdx());
        }
        int idx1 = static_cast<int>(val.as_array()[sz - 1].as_int64());
        int idx2 = static_cast<int>(val.as_array()[0].as_int64());
        atomRing.push_back(idx1);
        const auto &bnd = mol->getBondBetweenAtoms(idx1, idx2);
        CHECK_INVARIANT(bnd, "no bond found for ring");
        bondRing.push_back(bnd->getIdx());
        ri->addRing(atomRing, bondRing);
      }
    }
  }
}

void processMol(RWMol *mol, const bj::value &molval,
                const DefaultValueCache &atomDefaults,
                const DefaultValueCache &bondDefaults,
                const JSONParseParameters &params) {
  if (molval.as_object().contains("name")) {
    mol->setProp(common_properties::_Name,
                 molval.at("name").as_string());
  }
  if (!molval.as_object().contains("atoms")) {
    throw FileParseException("Bad Format: missing atoms in JSON");
  }
  if (!molval.as_object().contains("bonds")) {
    throw FileParseException("Bad Format: missing bonds in JSON");
  }

  for (const auto &atomVal : molval.at("atoms").as_array()) {
    readAtom(mol, atomVal, atomDefaults, params);
  }
  bool needStereoLoop = false;
  for (const auto &bondVal : molval.at("bonds").as_array()) {
    readBond(mol, bondVal, bondDefaults, needStereoLoop);
  }
  if (needStereoLoop) {
    // need to set bond stereo after the bonds are there
    unsigned int bidx = 0;
    for (const auto &bondVal : molval.at("bonds").as_array()) {
      Bond *bnd = mol->getBondWithIdx(bidx++);
      readBondStereo(bnd, bondVal, bondDefaults);
    }
  }

  if (molval.as_object().contains("stereoGroups")) {
    readStereoGroups(mol, molval.at("stereoGroups"));
  }
  if (molval.as_object().contains("substanceGroups")) {
    readSubstanceGroups(mol, molval.at("substanceGroups"));
  }
  if (params.parseConformers && molval.as_object().contains("conformers")) {
    for (const auto &confVal : molval.at("conformers").as_array()) {
      auto *conf = new Conformer(mol->getNumAtoms());
      readConformer(conf, confVal);
      mol->addConformer(conf, true);
    }
  }

  if (params.parseProperties && molval.as_object().contains("properties")) {
    parseProperties(*mol, molval.at("properties"));
  }

  if (molval.as_object().contains("extensions")) {
    for (const auto &propVal : molval.at("extensions").as_array()) {
      if (!propVal.as_object().contains("name")) {
        throw FileParseException(
            "Bad Format: representation has no name member");
      }
      if (propVal.at("name").as_string() ==
          std::string("rdkitRepresentation")) {
        readRDKitRepresentation(mol, propVal, params);
      } else if (propVal.at("name").as_string() ==
                 std::string("partialCharges")) {
        readPartialCharges(mol, propVal, params);
      } else if (propVal.at("name").as_string() ==
                 std::string("rdkitQueries")) {
        readQueries(mol, propVal, atomDefaults, bondDefaults, params);
      }
    }
  }
  mol->setProp(common_properties::_StereochemDone, 1);
}

std::vector<boost::shared_ptr<ROMol>> DocToMols(
    bj::value &doc, const JSONParseParameters &params) {
  std::vector<boost::shared_ptr<ROMol>> res;

  // some error checking
  if (!doc.is_object()) {
    throw FileParseException("Bad Format: JSON should be an object");
  }

  if (doc.as_object().contains("commonchem")) {
    auto jobj = doc.at("commonchem").if_object();
    if (!jobj || !jobj->contains("version")) {
      throw FileParseException("Bad Format: missing version in JSON");
    }
    if (jobj->at("version").as_int64() != currentMolJSONVersion) {
      throw FileParseException("Bad Format: bad version in JSON");
    }
  } else if (doc.as_object().contains("rdkitjson")) {
    if (!doc.at("rdkitjson").is_object() ||
        !doc.at("rdkitjson").as_object().contains("version")) {
      throw FileParseException("Bad Format: missing version in JSON");
    }
    // FIX: we want to be backwards compatible
    // Version 10 files can be read by 11, but not vice versa.
    if (int jsonVersion =
            static_cast<int>(doc.at("rdkitjson").at("version").as_int64());
        jsonVersion > currentRDKitJSONVersion || jsonVersion < 10) {
      throw FileParseException("Bad Format: bad version in JSON");
    }
  } else {
    throw FileParseException("Bad Format: missing header in JSON");
  }

  bj::value atomDefaults_;
  if (doc.as_object().contains("defaults") &&
      doc.at("defaults").as_object().contains("atom")) {
    atomDefaults_ = doc.at("defaults").at("atom");
    if (!atomDefaults_.is_object()) {
      throw FileParseException("Bad Format: atomDefaults is not an object");
    }
  }
  const DefaultValueCache atomDefaults(atomDefaults_);

  bj::value bondDefaults_;
  if (doc.as_object().contains("defaults") &&
      doc.at("defaults").as_object().contains("bond")) {
    bondDefaults_ = doc.at("defaults").at("bond");
    if (!bondDefaults_.is_object()) {
      throw FileParseException("Bad Format: bondDefaults is not an object");
    }
  }
  const DefaultValueCache bondDefaults(bondDefaults_);

  if (doc.as_object().contains("molecules")) {
    if (!doc.at("molecules").is_array()) {
      throw FileParseException("Bad Format: molecules is not an array");
    }
    for (const auto &molval : doc.at("molecules").as_array()) {
      std::unique_ptr<RWMol> mol(new RWMol());
      processMol(mol.get(), molval, atomDefaults, bondDefaults, params);
      mol->updatePropertyCache(params.strictValenceCheck);
      mol->setProp(common_properties::_StereochemDone, 1);
      res.emplace_back(static_cast<ROMol *>(mol.release()));
    }
  }

  return res;
}

}  // namespace

std::vector<boost::shared_ptr<ROMol>> JSONDataStreamToMols(
    std::istream *inStream, const JSONParseParameters &params) {
  PRECONDITION(inStream, "no stream");

  std::string jsonString((std::istreambuf_iterator<char>(*inStream)),
                         std::istreambuf_iterator<char>());
  bj::monotonic_resource mr;
  bj::value doc = bj::parse(jsonString, &mr);

  return DocToMols(doc, params);
}
std::vector<boost::shared_ptr<ROMol>> JSONDataToMols(
    const std::string &jsonBlock, const JSONParseParameters &params) {
  bj::monotonic_resource mr;
  bj::value doc = bj::parse(jsonBlock, &mr);

  return DocToMols(doc, params);
}

}  // namespace MolInterchange
}  // end of namespace RDKit
