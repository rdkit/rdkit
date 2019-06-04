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

std::string svg_(const ROMol &m,
                 const std::vector<unsigned int> *atomIds = nullptr) {
  MolDraw2DSVG drawer(250, 200);
  std::vector<int> *highlight_atoms = nullptr;
  if (atomIds) {
    highlight_atoms = new std::vector<int>;
    highlight_atoms->reserve(atomIds->size());
    for (auto ai : *atomIds) {
      highlight_atoms->push_back(ai);
    }
  }
  MolDraw2DUtils::prepareAndDrawMolecule(drawer, m, "", highlight_atoms);
  drawer.finishDrawing();
  delete highlight_atoms;

  return drawer.getDrawingText();
}
}  // namespace

std::string JSMol::get_smiles() const {
  if (!d_mol) return "";
  return MolToSmiles(*d_mol);
}
std::string JSMol::get_svg() const {
  if (!d_mol) return "";
  return svg_(*d_mol);
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
  return svg_(*d_mol, &atomIds);
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
std::string JSMol::get_substruct_match(const JSMol &q) const {
  std::string res = "{}";
  if (!d_mol) return res;

  MatchVectType match;
  if (SubstructMatch(*d_mol, *(q.d_mol), match)) {
    rj::Document doc;
    doc.SetArray();
    rj::Document::AllocatorType &allocator = doc.GetAllocator();
    for (const auto &pr : match) {
      doc.PushBack(pr.second, allocator);
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
