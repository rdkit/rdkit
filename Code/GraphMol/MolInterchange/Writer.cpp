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

namespace rj = rapidjson;

namespace RDKit {

namespace MolInterchange {

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
