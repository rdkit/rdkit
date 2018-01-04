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
#include <GraphMol/MolInterchange/MolInterchange.h>
#include <RDGeneral/FileParseException.h>

#include <sstream>
#include <exception>

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>

namespace rj = rapidjson;

namespace RDKit {

namespace MolInterchange {

namespace {
int getIntDefaultValue(const char *key, const rj::Value &from,
                       const rj::Value &defaults) {
  rj::Value::ConstMemberIterator endp = from.MemberEnd();
  rj::Value::ConstMemberIterator miter = from.FindMember(key);
  if (miter == endp) {
    miter = defaults.FindMember(key);
    endp = defaults.MemberEnd();
  }
  if (miter != endp) {
    if (!miter->value.IsInt())
      throw FileParseException(std::string("Bad format: value of ") +
                               std::string(key) +
                               std::string(" is not an int"));
    return miter->value.GetInt();
  }
  return 0;
}
bool getBoolDefaultValue(const char *key, const rj::Value &from,
                         const rj::Value &defaults) {
  rj::Value::ConstMemberIterator endp = from.MemberEnd();
  rj::Value::ConstMemberIterator miter = from.FindMember(key);
  if (miter == endp) {
    miter = defaults.FindMember(key);
    endp = defaults.MemberEnd();
  }
  if (miter != endp) {
    if (!miter->value.IsBool())
      throw FileParseException(std::string("Bad format: value of ") +
                               std::string(key) +
                               std::string(" is not a bool"));
    return miter->value.GetBool();
  }
  return false;
}
std::string getStringDefaultValue(const char *key, const rj::Value &from,
                                  const rj::Value &defaults) {
  rj::Value::ConstMemberIterator endp = from.MemberEnd();
  rj::Value::ConstMemberIterator miter = from.FindMember(key);
  if (miter == endp) {
    miter = defaults.FindMember(key);
    endp = defaults.MemberEnd();
  }
  if (miter != endp) {
    if (!miter->value.IsString())
      throw FileParseException(std::string("Bad format: value of ") +
                               std::string(key) +
                               std::string(" is not a string"));
    return miter->value.GetString();
  }
  return "";
}

void processMol(RWMol *mol, const rj::Value &molval,
                const rj::Value &atomDefaults, const rj::Value &bondDefaults) {
  static const std::map<std::string, Atom::ChiralType> chilookup = {
      {"unspecified", Atom::CHI_UNSPECIFIED},
      {"cw", Atom::CHI_TETRAHEDRAL_CW},
      {"ccw", Atom::CHI_TETRAHEDRAL_CCW},
      {"other", Atom::CHI_OTHER}};
  static const std::map<unsigned int, Bond::BondType> bolookup = {
      {1, Bond::SINGLE}, {2, Bond::DOUBLE}, {3, Bond::TRIPLE}};

  if (molval.HasMember("name")) {
    mol->setProp("_Name", molval["name"].GetString());
  }

  for (const auto &atomVal : molval["atoms"].GetArray()) {
    Atom *at = new Atom(getIntDefaultValue("Z", atomVal, atomDefaults));
    at->setNoImplicit(true);
    at->setNumExplicitHs(getIntDefaultValue("impHs", atomVal, atomDefaults));
    at->setFormalCharge(getIntDefaultValue("chg", atomVal, atomDefaults));
    at->setNumRadicalElectrons(
        getIntDefaultValue("nRad", atomVal, atomDefaults));
    std::string stereo = getStringDefaultValue("stereo", atomVal, atomDefaults);
    if (chilookup.find(stereo) == chilookup.end())
      throw FileParseException("Bad Format: bad stereo value for atom");
    at->setChiralTag(chilookup.find(stereo)->second);
    mol->addAtom(at);
  }
  for (const auto &bondVal : molval["bonds"].GetArray()) {
    const auto &aids = bondVal["atoms"].GetArray();
    unsigned int bid = mol->addBond(aids[0].GetInt(), aids[1].GetInt()) - 1;
    Bond *bnd = mol->getBondWithIdx(bid);
    unsigned int bo = getIntDefaultValue("bo", bondVal, bondDefaults);
    if (bolookup.find(bo) == bolookup.end())
      throw FileParseException("Bad Format: bad bond order for bond");
    bnd->setBondType(bolookup.find(bo)->second);
  }
  mol->updatePropertyCache(false);
  mol->setProp("_StereochemDone", 1);
}

std::vector<boost::shared_ptr<RWMol>> DocToMols(rj::Document &doc) {
  std::vector<boost::shared_ptr<RWMol>> res;

  // some error checking
  if (!doc.IsObject())
    throw FileParseException("Bad Format: JSON should be an object");
  if (!doc.HasMember("moljson-header"))
    throw FileParseException("Bad Format: missing header in JSON");
  if (!doc["moljson-header"].HasMember("version"))
    throw FileParseException("Bad Format: missing version in JSON");
  if (doc["moljson-header"]["version"].GetInt() != 10)
    throw FileParseException("Bad Format: bad version in JSON");

  rj::Value atomDefaults;
  if (doc.HasMember("atomDefaults")) {
    atomDefaults = doc["atomDefaults"];
    if (!atomDefaults.IsObject())
      throw FileParseException("Bad Format: atomDefaults is not an object");
  }

  rj::Value bondDefaults;
  if (doc.HasMember("bondDefaults")) {
    bondDefaults = doc["bondDefaults"];
    if (!bondDefaults.IsObject())
      throw FileParseException("Bad Format: bondDefaults is not an object");
  }

  if (doc.HasMember("molecules")) {
    if (!doc["molecules"].IsArray())
      throw FileParseException("Bad Format: molecules is not an array");
    for (const auto &molval : doc["molecules"].GetArray()) {
      RWMol *mol = new RWMol();
      processMol(mol, molval, atomDefaults, bondDefaults);
      res.push_back(boost::shared_ptr<RWMol>(mol));
    }
  }

  return res;
}
}  // end of anonymous namespace

std::vector<boost::shared_ptr<RWMol>> JSONDataStreamToMols(
    std::istream *inStream) {
  PRECONDITION(inStream, "no stream");

  rj::IStreamWrapper isw(*inStream);
  rj::Document doc;
  doc.ParseStream(isw);

  return (DocToMols(doc));
}
std::vector<boost::shared_ptr<RWMol>> JSONDataToMols(
    const std::string &jsonBlock) {
  rj::Document doc;
  doc.Parse(jsonBlock.c_str());
  return (DocToMols(doc));
}
}  // end of namespace MolInterchange
}  // end of namespace RDKit
