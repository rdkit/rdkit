//
//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <string>
#include <iostream>
#include "minilib.h"

#include <RDGeneral/versions.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Descriptors/Property.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <DataStructs/BitOps.h>

#include <INCHI-API/inchi.h>

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

namespace rj = rapidjson;

using namespace RDKit;

namespace {
ROMol *mol_from_input(const std::string &input) {
  RWMol *res = nullptr;
  if (input.find("M  END") != std::string::npos) {
    bool sanitize = false;
    res = MolBlockToMol(input, sanitize);
  } else {
    SmilesParserParams ps;
    ps.sanitize = false;
    res = SmilesToMol(input, ps);
  }
  if (res) {
    try {
      MolOps::sanitizeMol(*res);
      MolOps::assignStereochemistry(*res, true, true, true);
    } catch (...) {
      delete res;
      res = nullptr;
    }
  }
  return res;
}

ROMol *qmol_from_input(const std::string &input) {
  RWMol *res = nullptr;
  if (input.find("M  END") != std::string::npos) {
    bool sanitize = false;
    res = MolBlockToMol(input, sanitize);
  } else {
    res = SmartsToMol(input);
  }
  return res;
}

std::string svg_(const ROMol &m, unsigned int w, unsigned int h,
                 const std::vector<unsigned int> *atomIds = nullptr,
                 const std::vector<unsigned int> *bondIds = nullptr) {
  MolDraw2DSVG drawer(w, h);
  std::vector<int> *highlight_atoms = nullptr;
  if (atomIds && atomIds->size()) {
    highlight_atoms = new std::vector<int>;
    highlight_atoms->reserve(atomIds->size());
    for (auto ai : *atomIds) {
      highlight_atoms->push_back(ai);
    }
  }

  std::vector<int> *highlight_bonds = nullptr;
  if (bondIds && bondIds->size()) {
    highlight_bonds = new std::vector<int>;
    highlight_bonds->reserve(bondIds->size());
    for (auto ai : *bondIds) {
      highlight_bonds->push_back(ai);
    }
  }

  MolDraw2DUtils::prepareAndDrawMolecule(drawer, m, "", highlight_atoms,
                                         highlight_bonds);
  drawer.finishDrawing();
  delete highlight_atoms;
  delete highlight_bonds;

  return drawer.getDrawingText();
}
}  // namespace

std::string JSMol::get_smiles() const {
  if (!d_mol) return "";
  return MolToSmiles(*d_mol);
}
std::string JSMol::get_svg(unsigned int w, unsigned int h) const {
  if (!d_mol) return "";
  return svg_(*d_mol, w, h);
}
std::string JSMol::get_svg_with_highlights(const std::string &details) const {
  if (!d_mol) return "";

  rj::Document doc;
  doc.Parse(details.c_str());
  if (!doc.IsObject()) return "Invalid JSON";

  std::vector<unsigned int> atomIds;
  if (!doc.HasMember("atoms") || !doc["atoms"].IsArray()) {
    return "JSON doesn't contain 'atoms' field, or it is not an array";
  }
  for (const auto &molval : doc["atoms"].GetArray()) {
    if (!molval.IsInt()) return ("Atom IDs should be integers");
    atomIds.push_back(static_cast<unsigned int>(molval.GetInt()));
  }
  std::vector<unsigned int> bondIds;
  if (doc.HasMember("bonds")) {
    if (!doc["bonds"].IsArray()) {
      return "JSON contain 'bonds' field, but it is not an array";
    }
    for (const auto &molval : doc["bonds"].GetArray()) {
      if (!molval.IsInt()) return ("Bond IDs should be integers");
      bondIds.push_back(static_cast<unsigned int>(molval.GetInt()));
    }
  }

  unsigned int w = d_defaultWidth;
  if (doc.HasMember("width")) {
    if (!doc["width"].IsUint()) {
      return "JSON contains 'width' field, but it is not an unsigned int";
    }
    w = doc["width"].GetUint();
  }

  unsigned int h = d_defaultHeight;
  if (doc.HasMember("height")) {
    if (!doc["height"].IsUint()) {
      return "JSON contains 'height' field, but it is not an unsigned int";
    }
    h = doc["height"].GetUint();
  }

  return svg_(*d_mol, w, h, &atomIds, &bondIds);
}
std::string JSMol::get_inchi() const {
  if (!d_mol) return "";
  ExtraInchiReturnValues rv;
  return MolToInchi(*d_mol, rv);
}
std::string JSMol::get_molblock() const {
  if (!d_mol) return "";
  return MolToMolBlock(*d_mol);
}

namespace {
void get_sss_json(const ROMol *d_mol, const ROMol *q_mol,
                  const MatchVectType &match, rj::Value &obj,
                  rj::Document &doc) {
  rj::Value rjAtoms(rj::kArrayType);
  for (const auto &pr : match) {
    rjAtoms.PushBack(pr.second, doc.GetAllocator());
  }
  obj.AddMember("atoms", rjAtoms, doc.GetAllocator());

  rj::Value rjBonds(rj::kArrayType);
  for (const auto qbond : q_mol->bonds()) {
    unsigned int idx1 = match[qbond->getBeginAtomIdx()].second;
    unsigned int idx2 = match[qbond->getEndAtomIdx()].second;
    const auto bond = d_mol->getBondBetweenAtoms(idx1, idx2);
    if (bond != nullptr) {
      rjBonds.PushBack(bond->getIdx(), doc.GetAllocator());
    }
  }
  obj.AddMember("bonds", rjBonds, doc.GetAllocator());
}
}  // namespace

std::string JSMol::get_substruct_match(const JSMol &q) const {
  std::string res = "{}";
  if (!d_mol || !q.d_mol) return res;

  MatchVectType match;
  if (SubstructMatch(*d_mol, *(q.d_mol), match)) {
    rj::Document doc;
    doc.SetObject();
    get_sss_json(d_mol.get(), q.d_mol.get(), match, doc, doc);
    rj::StringBuffer buffer;
    rj::Writer<rj::StringBuffer> writer(buffer);
    doc.Accept(writer);
    res = buffer.GetString();
  }

  return res;
}

std::string JSMol::get_substruct_matches(const JSMol &q) const {
  std::string res = "{}";
  if (!d_mol || !q.d_mol) return res;

  auto matches = SubstructMatch(*d_mol, (*q.d_mol));
  if (!matches.empty()) {
    rj::Document doc;
    doc.SetArray();

    for (const auto match : matches) {
      rj::Value rjMatch(rj::kObjectType);
      get_sss_json(d_mol.get(), q.d_mol.get(), match, rjMatch, doc);
      doc.PushBack(rjMatch, doc.GetAllocator());
    }

    rj::StringBuffer buffer;
    rj::Writer<rj::StringBuffer> writer(buffer);
    doc.Accept(writer);
    res = buffer.GetString();
  }

  return res;
}

std::string JSMol::get_descriptors() const {
  if (!d_mol) return "{}";
  rj::Document doc;
  doc.SetObject();

  Descriptors::Properties props;
  std::vector<std::string> dns = props.getPropertyNames();
  std::vector<double> dvs = props.computeProperties(*d_mol);
  for (size_t i = 0; i < dns.size(); ++i) {
    rj::Value v(dvs[i]);
    const auto srt = rj::StringRef(dns[i].c_str());
    doc.AddMember(srt, v, doc.GetAllocator());
  }

  if (std::find(dns.begin(), dns.end(), std::string("amw")) == dns.end()) {
    rj::Value v(Descriptors::calcAMW(*d_mol));
    doc.AddMember("amw", v, doc.GetAllocator());
  }

  rj::StringBuffer buffer;
  rj::Writer<rj::StringBuffer> writer(buffer);
  writer.SetMaxDecimalPlaces(5);
  doc.Accept(writer);
  return buffer.GetString();
}

std::string JSMol::get_morgan_fp(unsigned int radius,
                                 unsigned int fplen) const {
  if (!d_mol) return "";
  auto fp = MorganFingerprints::getFingerprintAsBitVect(*d_mol, radius, fplen);
  std::string res = BitVectToText(*fp);
  delete fp;
  return res;
}

std::string JSMol::get_stereo_tags() const {
  if (!d_mol) return "{}";
  rj::Document doc;
  doc.SetObject();

  bool cleanIt = true;
  bool force = true;
  bool flagPossibleStereocenters = true;
  MolOps::assignStereochemistry(*d_mol, cleanIt, force,
                                flagPossibleStereocenters);

  rj::Value rjAtoms(rj::kArrayType);
  for (const auto atom : d_mol->atoms()) {
    std::string cip;
    if (!atom->getPropIfPresent(common_properties::_CIPCode, cip)) {
      if (atom->hasProp(common_properties::_ChiralityPossible)) {
        cip = "?";
      }
    }
    if (!cip.empty()) {
      cip = "(" + cip + ")";
      rj::Value entry(rj::kArrayType);
      entry.PushBack(atom->getIdx(), doc.GetAllocator());
      rj::Value v;
      v.SetString(cip.c_str(), cip.size(), doc.GetAllocator());
      entry.PushBack(v, doc.GetAllocator());
      rjAtoms.PushBack(entry, doc.GetAllocator());
    }
  }
  doc.AddMember("CIP_atoms", rjAtoms, doc.GetAllocator());

  rj::Value rjBonds(rj::kArrayType);
  for (const auto bond : d_mol->bonds()) {
    std::string cip = "";
    if (bond->getStereo() == Bond::STEREOE)
      cip = "(E)";
    else if (bond->getStereo() == Bond::STEREOZ)
      cip = "(Z)";
    if (cip.empty()) continue;
    rj::Value entry(rj::kArrayType);
    entry.PushBack(bond->getBeginAtomIdx(), doc.GetAllocator());
    entry.PushBack(bond->getEndAtomIdx(), doc.GetAllocator());
    rj::Value v;
    v.SetString(cip.c_str(), cip.size(), doc.GetAllocator());
    entry.PushBack(v, doc.GetAllocator());
    rjBonds.PushBack(entry, doc.GetAllocator());
  }

  doc.AddMember("CIP_bonds", rjBonds, doc.GetAllocator());

  rj::StringBuffer buffer;
  rj::Writer<rj::StringBuffer> writer(buffer);
  doc.Accept(writer);
  return buffer.GetString();
}

std::string JSMol::get_aromatic_form() const {
  if (!d_mol) return "";

  RWMol molCopy(*d_mol);
  MolOps::setAromaticity(molCopy);

  bool includeStereo = true;
  int confId = -1;
  bool kekulize = false;
  return MolToMolBlock(molCopy, includeStereo, confId, kekulize);
}

std::string JSMol::get_kekule_form() const {
  if (!d_mol) return "";

  RWMol molCopy(*d_mol);
  MolOps::Kekulize(molCopy);

  bool includeStereo = true;
  int confId = -1;
  bool kekulize = true;
  return MolToMolBlock(molCopy, includeStereo, confId, kekulize);
}

std::string JSMol::get_new_coords(bool useCoordGen) const {
  if (!d_mol) return "";

  RWMol molCopy(*d_mol);

#ifdef RDK_BUILD_COORDGEN_SUPPORT
  bool oprefer = RDDepict::preferCoordGen;
  RDDepict::preferCoordGen = useCoordGen;
#endif
  RDDepict::compute2DCoords(molCopy);
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = oprefer;
#endif

  return MolToMolBlock(molCopy);
}

std::string JSMol::remove_hs() const {
  if (!d_mol) return "";

  RWMol molCopy(*d_mol);
  MolOps::removeAllHs(molCopy);

  bool includeStereo = true;
  int confId = -1;
  bool kekulize = true;
  return MolToMolBlock(molCopy, includeStereo, confId, kekulize);
}

std::string JSMol::add_hs() const {
  if (!d_mol) return "";

  RWMol molCopy(*d_mol);
  MolOps::addHs(molCopy);

  // RDDepict::generateDepictionMatching2DStructure(molCopy, *d_mol);

  bool includeStereo = true;
  int confId = -1;
  bool kekulize = true;
  return MolToMolBlock(molCopy, includeStereo, confId, kekulize);
}

std::string get_inchikey_for_inchi(const std::string &input) {
  return InchiToInchiKey(input);
}

JSMol *get_mol(const std::string &input) {
  ROMol *mol = mol_from_input(input);
  return new JSMol(mol);
}

JSMol *get_qmol(const std::string &input) {
  ROMol *mol = qmol_from_input(input);
  return new JSMol(mol);
}

std::string version() { return std::string(rdkitVersion); }
