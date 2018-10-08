//
//  Copyright (C) 2018 Greg Landrum
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

#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/MolInterchange/details.h>
#include <RDGeneral/FileParseException.h>

#include <sstream>
#include <exception>
#include <map>

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>
#include <rapidjson/pointer.h>

namespace rj = rapidjson;

namespace RDKit {

namespace MolInterchange {

namespace {
struct DefaultValueCache {
  DefaultValueCache(const rj::Value &defs) : rjDefaults(defs){};
  const rj::Value &rjDefaults;
  mutable std::map<const char *, int> intMap;
  mutable std::map<const char *, bool> boolMap;
  mutable std::map<const char *, std::string> stringMap;
  int getInt(const char *key) const {
    PRECONDITION(key, "no key");
    const auto &lookup = intMap.find(key);
    if (lookup != intMap.end()) return lookup->second;
    const auto &miter = rjDefaults.FindMember(key);
    if (miter != rjDefaults.MemberEnd()) {
      if (!miter->value.IsInt())
        throw FileParseException(std::string("Bad format: value of ") +
                                 std::string(key) +
                                 std::string(" is not an int"));
      int res = miter->value.GetInt();
      intMap[key] = res;
      return res;
    }
    return 0;
  }
  bool getBool(const char *key) const {
    PRECONDITION(key, "no key");
    const auto &lookup = boolMap.find(key);
    if (lookup != boolMap.end()) return lookup->second;
    const auto &miter = rjDefaults.FindMember(key);
    if (miter != rjDefaults.MemberEnd()) {
      if (!miter->value.IsBool())
        throw FileParseException(std::string("Bad format: value of ") +
                                 std::string(key) +
                                 std::string(" is not a bool"));
      bool res = miter->value.GetBool();
      boolMap[key] = res;
      return res;
    }
    return false;
  }
  std::string getString(const char *key) const {
    PRECONDITION(key, "no key");
    const auto &lookup = stringMap.find(key);
    if (lookup != stringMap.end()) return lookup->second;
    const auto &miter = rjDefaults.FindMember(key);
    if (miter != rjDefaults.MemberEnd()) {
      if (!miter->value.IsString())
        throw FileParseException(std::string("Bad format: value of ") +
                                 std::string(key) +
                                 std::string(" is not a string"));
      std::string res = miter->value.GetString();
      stringMap[key] = res;
      return res;
    }
    return "";
  }
};

int getIntDefaultValue(const char *key, const rj::Value &from,
                       const DefaultValueCache &defaults) {
  PRECONDITION(key, "no key");
  auto endp = from.MemberEnd();
  auto miter = from.FindMember(key);
  if (miter != endp) {
    if (!miter->value.IsInt())
      throw FileParseException(std::string("Bad format: value of ") +
                               std::string(key) +
                               std::string(" is not an int"));
    return miter->value.GetInt();
  }
  return defaults.getInt(key);
}
bool getBoolDefaultValue(const char *key, const rj::Value &from,
                         const DefaultValueCache &defaults) {
  PRECONDITION(key, "no key");
  auto endp = from.MemberEnd();
  auto miter = from.FindMember(key);
  if (miter != endp) {
    if (!miter->value.IsBool())
      throw FileParseException(std::string("Bad format: value of ") +
                               std::string(key) +
                               std::string(" is not a bool"));
    return miter->value.GetBool();
  }
  return defaults.getBool(key);
}
std::string getStringDefaultValue(const char *key, const rj::Value &from,
                                  const DefaultValueCache &defaults) {
  PRECONDITION(key, "no key");
  auto endp = from.MemberEnd();
  auto miter = from.FindMember(key);
  if (miter != endp) {
    if (!miter->value.IsString())
      throw FileParseException(std::string("Bad format: value of ") +
                               std::string(key) +
                               std::string(" is not a string"));
    return miter->value.GetString();
  }
  return defaults.getString(key);
}

void readAtom(RWMol *mol, const rj::Value &atomVal,
              const DefaultValueCache &atomDefaults) {
  PRECONDITION(mol, "no mol");
  Atom *at = new Atom(getIntDefaultValue("z", atomVal, atomDefaults));
  at->setNoImplicit(true);
  at->setNumExplicitHs(getIntDefaultValue("impHs", atomVal, atomDefaults));
  at->setFormalCharge(getIntDefaultValue("chg", atomVal, atomDefaults));
  at->setNumRadicalElectrons(getIntDefaultValue("nRad", atomVal, atomDefaults));
  at->setIsotope(getIntDefaultValue("isotope", atomVal, atomDefaults));
  std::string stereo = getStringDefaultValue("stereo", atomVal, atomDefaults);
  if (chilookup.find(stereo) == chilookup.end())
    throw FileParseException("Bad Format: bad stereo value for atom");
  at->setChiralTag(chilookup.find(stereo)->second);
  bool updateLabel = false, takeOwnership = true;
  mol->addAtom(at, updateLabel, takeOwnership);
}

void readBond(RWMol *mol, const rj::Value &bondVal,
              const DefaultValueCache &bondDefaults, bool &needStereoLoop) {
  PRECONDITION(mol, "no mol");
  const auto &aids = bondVal["atoms"].GetArray();
  unsigned int bid = mol->addBond(aids[0].GetInt(), aids[1].GetInt()) - 1;
  Bond *bnd = mol->getBondWithIdx(bid);
  unsigned int bo = getIntDefaultValue("bo", bondVal, bondDefaults);
  if (bolookup.find(bo) == bolookup.end())
    throw FileParseException("Bad Format: bad bond order for bond");
  bnd->setBondType(bolookup.find(bo)->second);
  if (bondVal.HasMember("stereoAtoms")) needStereoLoop = true;
}

void readBondStereo(Bond *bnd, const rj::Value &bondVal,
                    const DefaultValueCache &bondDefaults) {
  PRECONDITION(bnd, "no bond");

  const auto &miter = bondVal.FindMember("stereoAtoms");
  if (miter != bondVal.MemberEnd()) {
    const auto aids = miter->value.GetArray();
    bnd->setStereoAtoms(aids[0].GetInt(), aids[1].GetInt());
    std::string stereo = getStringDefaultValue("stereo", bondVal, bondDefaults);
    if (stereolookup.find(stereo) == stereolookup.end())
      throw FileParseException("Bad Format: bond stereo value for bond");
    bnd->setStereo(stereolookup.find(stereo)->second);

  } else {
    if (bondVal.HasMember("stereo")) {
      throw FileParseException(
          "Bad Format: bond stereo provided without stereoAtoms");
    }
  }
}

void readConformer(Conformer *conf, const rj::Value &confVal) {
  PRECONDITION(conf, "no conformer");

  if (!confVal.HasMember("dim"))
    throw FileParseException("Bad Format: no conformer dimension");
  size_t dim = confVal["dim"].GetInt();
  if (dim == 2)
    conf->set3D(false);
  else if (dim == 3)
    conf->set3D(true);
  else
    throw FileParseException("Bad Format: conformer dimension != 2 or 3");
  if (!confVal.HasMember("coords"))
    throw FileParseException("Bad Format: no conformer coords");
  size_t idx = 0;
  for (const auto &ptVal : confVal["coords"].GetArray()) {
    const auto &arr = ptVal.GetArray();
    if (arr.Size() != dim)
      throw FileParseException("coordinate contains wrong number of values");
    RDGeom::Point3D pt(arr[0].GetFloat(), arr[1].GetFloat(),
                       (dim == 3 ? arr[2].GetFloat() : 0.0));
    conf->setAtomPos(idx++, pt);
  }
  if (idx != conf->getNumAtoms())
    throw FileParseException(
        "Bad Format: conformer doesn't contain coordinates for all atoms");
}

void readPartialCharges(RWMol *mol, const rj::Value &repVal,
                        const JSONParseParameters &params) {
  RDUNUSED_PARAM(params);
  PRECONDITION(mol, "no molecule");
  PRECONDITION(repVal["name"].GetString() == std::string("partialCharges"),
               "bad charges");
  if (!repVal.HasMember("formatVersion"))
    throw FileParseException("Bad Format: missing version");
  if (repVal["formatVersion"].GetInt() > currentChargeRepresentationVersion) {
    BOOST_LOG(rdWarningLog)
        << "partialCharges version " << repVal["formatVersion"].GetInt()
        << " too recent. Ignoring it." << std::endl;
    return;
  }
  {
    const auto &miter = repVal.FindMember("values");
    if (miter != repVal.MemberEnd()) {
      if (miter->value.GetArray().Size() != mol->getNumAtoms())
        throw FileParseException(
            "Bad Format: size of values array != num atoms");
      for (unsigned int idx = 0; idx != mol->getNumAtoms(); ++idx) {
        const auto &val = miter->value.GetArray()[idx];
        if (!val.IsDouble())
          throw FileParseException("Bad Format: partial charge not double");
        mol->getAtomWithIdx(idx)->setProp(common_properties::_GasteigerCharge,
                                          val.GetDouble());
      }
    }
  }
}
void readRDKitRepresentation(RWMol *mol, const rj::Value &repVal,
                             const JSONParseParameters &params) {
  PRECONDITION(mol, "no molecule");
  PRECONDITION(repVal["name"].GetString() == std::string("rdkitRepresentation"),
               "bad representation");
  if (!repVal.HasMember("formatVersion"))
    throw FileParseException("Bad Format: missing format_version");
  if (repVal["formatVersion"].GetInt() > 1) {
    BOOST_LOG(rdWarningLog) << "RDKit representation format version "
                            << repVal["formatVersion"].GetInt()
                            << " too recent. Ignoring it." << std::endl;
    return;
  }
  {
    const auto &miter = repVal.FindMember("aromaticAtoms");
    if (miter != repVal.MemberEnd()) {
      for (const auto &val : miter->value.GetArray()) {
        if (!val.IsInt())
          throw FileParseException("Bad Format: aromaticAtom not int");
        mol->getAtomWithIdx(val.GetInt())->setIsAromatic(true);
      }
    }
  }
  {
    const auto &miter = repVal.FindMember("aromaticBonds");
    if (miter != repVal.MemberEnd()) {
      for (const auto &val : miter->value.GetArray()) {
        if (!val.IsInt())
          throw FileParseException("Bad Format: aromaticBond not int");
        mol->getBondWithIdx(val.GetInt())->setIsAromatic(true);
        if (params.setAromaticBonds) {
          mol->getBondWithIdx(val.GetInt())->setBondType(Bond::AROMATIC);
        }
      }
    }
  }
  {
    const auto &miter = repVal.FindMember("cipRanks");
    if (miter != repVal.MemberEnd()) {
      size_t i = 0;
      for (const auto &val : miter->value.GetArray()) {
        if (!val.IsInt())
          throw FileParseException("Bad Format: ciprank not int");
        mol->getAtomWithIdx(i++)->setProp(
            common_properties::_CIPRank,
            static_cast<unsigned int>(val.GetInt()));
      }
    }
  }
  {
    const auto &miter = repVal.FindMember("cipCodes");
    if (miter != repVal.MemberEnd()) {
      for (const auto &val : miter->value.GetArray()) {
        if (!val.IsArray())
          throw FileParseException("Bad Format: CIPCode not string");
        mol->getAtomWithIdx(val[0].GetInt())
            ->setProp(common_properties::_CIPCode, val[1].GetString());
      }
    }
  }
  {
    const auto &miter = repVal.FindMember("atomRings");
    if (miter != repVal.MemberEnd()) {
      CHECK_INVARIANT(!mol->getRingInfo()->isInitialized(),
                      "rings already initialized");
      auto ri = mol->getRingInfo();
      ri->initialize();
      for (const auto &val : miter->value.GetArray()) {
        if (!val.IsArray())
          throw FileParseException("Bad Format: atomRing not array");
        INT_VECT atomRing;
        INT_VECT bondRing;
        size_t sz = val.Size();
        atomRing.reserve(sz);
        bondRing.reserve(sz);
        for (size_t i = 0; i < sz - 1; ++i) {
          int idx1 = val[i].GetInt();
          int idx2 = val[i + 1].GetInt();
          atomRing.push_back(idx1);
          const auto &bnd = mol->getBondBetweenAtoms(idx1, idx2);
          CHECK_INVARIANT(bnd, "no bond found for ring");
          bondRing.push_back(bnd->getIdx());
        }
        int idx1 = val[sz - 1].GetInt();
        int idx2 = val[0].GetInt();
        atomRing.push_back(idx1);
        const auto &bnd = mol->getBondBetweenAtoms(idx1, idx2);
        CHECK_INVARIANT(bnd, "no bond found for ring");
        bondRing.push_back(bnd->getIdx());
        ri->addRing(atomRing, bondRing);
      }
    }
  }
}

void processMol(RWMol *mol, const rj::Value &molval,
                const DefaultValueCache &atomDefaults,
                const DefaultValueCache &bondDefaults,
                const JSONParseParameters &params) {
  if (molval.HasMember("name")) {
    mol->setProp(common_properties::_Name, molval["name"].GetString());
  }
  if (!molval.HasMember("atoms"))
    throw FileParseException("Bad Format: missing atoms in JSON");
  if (!molval.HasMember("bonds"))
    throw FileParseException("Bad Format: missing bonds in JSON");

  for (const auto &atomVal : molval["atoms"].GetArray()) {
    readAtom(mol, atomVal, atomDefaults);
  }
  bool needStereoLoop = false;
  for (const auto &bondVal : molval["bonds"].GetArray()) {
    readBond(mol, bondVal, bondDefaults, needStereoLoop);
  }
  if (needStereoLoop) {
    // need to set bond stereo after the bonds are there
    unsigned int bidx = 0;
    for (const auto &bondVal : molval["bonds"].GetArray()) {
      Bond *bnd = mol->getBondWithIdx(bidx++);
      readBondStereo(bnd, bondVal, bondDefaults);
    }
  }
  if (params.parseConformers && molval.HasMember("conformers")) {
    for (const auto &confVal : molval["conformers"].GetArray()) {
      Conformer *conf = new Conformer(mol->getNumAtoms());
      readConformer(conf, confVal);
      mol->addConformer(conf, true);
    }
  }

  if (params.parseProperties && molval.HasMember("properties")) {
    for (const auto &propVal : molval["properties"].GetObject()) {
      if (propVal.value.IsInt())
        mol->setProp(propVal.name.GetString(), propVal.value.GetInt());
      else if (propVal.value.IsDouble())
        mol->setProp(propVal.name.GetString(), propVal.value.GetDouble());
      else if (propVal.value.IsString())
        mol->setProp(propVal.name.GetString(), propVal.value.GetString());
    }
  }

  if (molval.HasMember("extensions")) {
    for (const auto &propVal : molval["extensions"].GetArray()) {
      if (!propVal.HasMember("name"))
        throw FileParseException(
            "Bad Format: representation has no name member");
      if (propVal["name"].GetString() == std::string("rdkitRepresentation")) {
        readRDKitRepresentation(mol, propVal, params);
      }
      if (propVal["name"].GetString() == std::string("partialCharges")) {
        readPartialCharges(mol, propVal, params);
      }
    }
  }
  mol->updatePropertyCache(false);
  mol->setProp(common_properties::_StereochemDone, 1);
}

std::vector<boost::shared_ptr<ROMol>> DocToMols(
    rj::Document &doc, const JSONParseParameters &params) {
  std::vector<boost::shared_ptr<ROMol>> res;

  // some error checking
  if (!doc.IsObject())
    throw FileParseException("Bad Format: JSON should be an object");
  if (!doc.HasMember("commonchem"))
    throw FileParseException("Bad Format: missing header in JSON");
  if (!doc["commonchem"].HasMember("version"))
    throw FileParseException("Bad Format: missing version in JSON");
  if (doc["commonchem"]["version"].GetInt() != currentMolJSONVersion)
    throw FileParseException("Bad Format: bad version in JSON");

  rj::Value atomDefaults_;
  if (rj::GetValueByPointer(doc, "/defaults/atom")) {
    atomDefaults_ = *rj::GetValueByPointer(doc, "/defaults/atom");
    if (!atomDefaults_.IsObject())
      throw FileParseException("Bad Format: atomDefaults is not an object");
  }
  const DefaultValueCache atomDefaults(atomDefaults_);

  rj::Value bondDefaults_;
  if (rj::GetValueByPointer(doc, "/defaults/bond")) {
    bondDefaults_ = *rj::GetValueByPointer(doc, "/defaults/bond");
    if (!bondDefaults_.IsObject())
      throw FileParseException("Bad Format: bondDefaults is not an object");
  }
  const DefaultValueCache bondDefaults(bondDefaults_);

  if (doc.HasMember("molecules")) {
    if (!doc["molecules"].IsArray())
      throw FileParseException("Bad Format: molecules is not an array");
    for (const auto &molval : doc["molecules"].GetArray()) {
      RWMol *mol = new RWMol();
      processMol(mol, molval, atomDefaults, bondDefaults, params);
      mol->updatePropertyCache(params.strictValenceCheck);
      mol->setProp(common_properties::_StereochemDone, 1);
      res.push_back(boost::shared_ptr<ROMol>(static_cast<ROMol *>(mol)));
    }
  }

  return res;
}

}  // end of anonymous namespace

std::vector<boost::shared_ptr<ROMol>> JSONDataStreamToMols(
    std::istream *inStream, const JSONParseParameters &params) {
  PRECONDITION(inStream, "no stream");

  rj::IStreamWrapper isw(*inStream);
  rj::Document doc;
  doc.ParseStream(isw);

  return DocToMols(doc, params);
}
std::vector<boost::shared_ptr<ROMol>> JSONDataToMols(
    const std::string &jsonBlock, const JSONParseParameters &params) {
  rj::Document doc;
  doc.Parse(jsonBlock.c_str());
  return DocToMols(doc, params);
}

}  // end of namespace MolInterchange
}  // end of namespace RDKit
