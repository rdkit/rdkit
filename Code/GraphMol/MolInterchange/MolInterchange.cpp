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
#include <RDGeneral/FileParseException.h>

#include <sstream>
#include <exception>
#include <map>

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

namespace rj = rapidjson;

namespace RDKit {

namespace MolInterchange {

static int currentMolJSONVersion = 10;
static int currentRDKitRepresentationVersion = 1;

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

static const std::map<std::string, Atom::ChiralType> chilookup = {
    {"unspecified", Atom::CHI_UNSPECIFIED},
    {"cw", Atom::CHI_TETRAHEDRAL_CW},
    {"ccw", Atom::CHI_TETRAHEDRAL_CCW},
    {"other", Atom::CHI_OTHER}};
static const std::map<Atom::ChiralType, std::string> inv_chilookup = {
    {Atom::CHI_UNSPECIFIED, "unspecified"},
    {Atom::CHI_TETRAHEDRAL_CW, "cw"},
    {Atom::CHI_TETRAHEDRAL_CCW, "ccw"},
    {Atom::CHI_OTHER, "other"}};

void readAtom(RWMol *mol, const rj::Value &atomVal,
              const DefaultValueCache &atomDefaults) {
  PRECONDITION(mol, "no mol");
  Atom *at = new Atom(getIntDefaultValue("Z", atomVal, atomDefaults));
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

static const std::map<unsigned int, Bond::BondType> bolookup = {
    {0, Bond::ZERO}, {1, Bond::SINGLE}, {2, Bond::DOUBLE}, {3, Bond::TRIPLE}};
static const std::map<Bond::BondType, unsigned int> inv_bolookup = {
    {Bond::ZERO, 0}, {Bond::SINGLE, 1}, {Bond::DOUBLE, 2}, {Bond::TRIPLE, 3}};
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

static const std::map<std::string, Bond::BondStereo> stereolookup = {
    {"unspecified", Bond::STEREONONE},
    {"cis", Bond::STEREOCIS},
    {"trans", Bond::STEREOTRANS},
    {"either", Bond::STEREOANY}};
static const std::map<Bond::BondStereo, std::string> inv_stereolookup = {
    {Bond::STEREONONE, "unspecified"}, {Bond::STEREOCIS, "cis"},
    {Bond::STEREOTRANS, "trans"},      {Bond::STEREOZ, "cis"},
    {Bond::STEREOE, "trans"},          {Bond::STEREOANY, "either"}};
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

void readRDKitRepresentation(RWMol *mol, const rj::Value &repVal,
                             const JSONParseParameters &params) {
  PRECONDITION(mol, "no molecule");
  PRECONDITION(repVal["toolkit"].GetString() == std::string("RDKit"),
               "bad representation");
  if (!repVal.HasMember("format_version"))
    throw FileParseException("Bad Format: missing format_version");
  if (repVal["format_version"].GetInt() > 1) {
    BOOST_LOG(rdWarningLog) << "RDKit representation format version "
                            << repVal["format_version"].GetInt()
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
        mol->getAtomWithIdx(i++)->setProp(common_properties::_CIPRank,
                                          val.GetInt());
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
  if (molval.HasMember("conformers")) {
    for (const auto &confVal : molval["conformers"].GetArray()) {
      Conformer *conf = new Conformer(mol->getNumAtoms());
      readConformer(conf, confVal);
      mol->addConformer(conf, true);
    }
  }

  if (molval.HasMember("molProperties")) {
    for (const auto &propVal : molval["molProperties"].GetObject()) {
      if (propVal.value.IsInt())
        mol->setProp(propVal.name.GetString(), propVal.value.GetInt());
      else if (propVal.value.IsDouble())
        mol->setProp(propVal.name.GetString(), propVal.value.GetDouble());
      else if (propVal.value.IsString())
        mol->setProp(propVal.name.GetString(), propVal.value.GetString());
    }
  }

  if (molval.HasMember("representations")) {
    for (const auto &propVal : molval["representations"].GetArray()) {
      if (!propVal.HasMember("toolkit"))
        throw FileParseException(
            "Bad Format: representation has no toolkit member");
      if (propVal["toolkit"].GetString() == std::string("RDKit")) {
        readRDKitRepresentation(mol, propVal, params);
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
  if (!doc.HasMember("moljson-header"))
    throw FileParseException("Bad Format: missing header in JSON");
  if (!doc["moljson-header"].HasMember("version"))
    throw FileParseException("Bad Format: missing version in JSON");
  if (doc["moljson-header"]["version"].GetInt() != currentMolJSONVersion)
    throw FileParseException("Bad Format: bad version in JSON");

  rj::Value atomDefaults_;
  if (doc.HasMember("atomDefaults")) {
    atomDefaults_ = doc["atomDefaults"];
    if (!atomDefaults_.IsObject())
      throw FileParseException("Bad Format: atomDefaults is not an object");
  }
  const DefaultValueCache atomDefaults(atomDefaults_);

  rj::Value bondDefaults_;
  if (doc.HasMember("bondDefaults")) {
    bondDefaults_ = doc["bondDefaults"];
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

namespace {
void initAtomDefaults(rj::Value &rjDefaults, rj::Document &document) {
  // "atomDefaults": {"Z": 6, "impHs": 0, "chg": 0, "stereo": "unspecified",
  // "nrad": 0, "isotope": 0},
  rjDefaults.AddMember("Z", 6, document.GetAllocator());
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
void initHeader(rj::Value &rjHeader, rj::Document &document, const char *nm) {
  rjHeader.AddMember("version", currentMolJSONVersion, document.GetAllocator());
  rj::Value nmv;
  nmv.SetString(rj::StringRef(nm));
  rjHeader.AddMember("name", nmv, document.GetAllocator());
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
  addIntVal(rjAtom, rjDefaults, "Z", atom.getAtomicNum(), doc);
  addIntVal(rjAtom, rjDefaults, "impHs", atom.getTotalNumHs(), doc);
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
  MolOps::Kekulize(mol, false);
  if (mol.hasProp(common_properties::_Name)) {
    rj::Value nmv;
    const std::string &nm =
        mol.getProp<std::string>(common_properties::_Name).c_str();
    nmv.SetString(nm.c_str(), nm.size(), doc.GetAllocator());
    rjMol.AddMember("name", nmv, doc.GetAllocator());
  }
  rj::Value rjAtoms(rj::kArrayType);
  for (const auto &at : mol.atoms()) {
    rj::Value rjAtom(rj::kObjectType);
    addAtom(*at, rjAtom, doc, atomDefaults);
    rjAtoms.PushBack(rjAtom, doc.GetAllocator());
  }
  rjMol.AddMember("atoms", rjAtoms, doc.GetAllocator());

  rj::Value rjBonds(rj::kArrayType);
  for (const auto &bnd : mol.bonds()) {
    rj::Value rjBond(rj::kObjectType);
    addBond(*bnd, rjBond, doc, bondDefaults);
    rjBonds.PushBack(rjBond, doc.GetAllocator());
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

  rj::Value representation(rj::kObjectType);
  representation.AddMember("toolkit", "RDKit", doc.GetAllocator());
  representation.AddMember("format_version", currentRDKitRepresentationVersion,
                           doc.GetAllocator());
  rj::Value toolkitVersion;
  toolkitVersion.SetString(rj::StringRef(rdkitVersion));
  representation.AddMember("toolkit_version", toolkitVersion,
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
  rjMol.AddMember("representations", rjReprs, doc.GetAllocator());
}
}  // end of anonymous namespace

template <typename T>
std::string MolsToJSONData(const std::vector<T> &mols, const char *name) {
  std::string res = "";
  rj::Document doc;
  doc.SetObject();

  rj::Value header(rj::kObjectType);
  initHeader(header, doc, name);
  doc.AddMember("moljson-header", header, doc.GetAllocator());

  rj::Value atomDefaults(rj::kObjectType);
  initAtomDefaults(atomDefaults, doc);
  doc.AddMember("atomDefaults", atomDefaults, doc.GetAllocator());

  rj::Value bondDefaults(rj::kObjectType);
  initBondDefaults(bondDefaults, doc);
  doc.AddMember("bondDefaults", bondDefaults, doc.GetAllocator());

  rj::Value rjMols(rj::kArrayType);
  for (const auto &mol : mols) {
    rj::Value rjMol(rj::kObjectType);
    // write mol;
    addMol(*mol, rjMol, doc, doc.FindMember("atomDefaults")->value,
           doc.FindMember("bondDefaults")->value);
    rjMols.PushBack(rjMol, doc.GetAllocator());
  }
  doc.AddMember("molecules", rjMols, doc.GetAllocator());

#if 0
  DefaultValueCache atomDefaults, bondDefaults;
  initAtomDefaults(atomDefaults);
  initBondDefaults(bondDefaults);
#endif

  rj::StringBuffer buffer;
  rj::Writer<rj::StringBuffer> writer(buffer);
  writer.SetMaxDecimalPlaces(4);
  doc.Accept(writer);
  return buffer.GetString();
};

template std::string MolsToJSONData<ROMol *>(const std::vector<ROMol *> &,
                                             const char *);
template std::string MolsToJSONData<RWMol *>(const std::vector<RWMol *> &,
                                             const char *);
template std::string MolsToJSONData<const ROMol *>(
    const std::vector<const ROMol *> &, const char *);
template std::string MolsToJSONData<const RWMol *>(
    const std::vector<const RWMol *> &, const char *);

}  // end of namespace MolInterchange
}  // end of namespace RDKit
