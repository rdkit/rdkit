//
//  Copyright (C) 2021 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include <string>
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
#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/Descriptors/Property.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Abbreviations/Abbreviations.h>
#include <DataStructs/BitOps.h>

#ifndef _MSC_VER
// shutoff some warnings from rapidjson
#if !defined(__clang__) and defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif
#endif
#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>
#ifndef _MSC_VER
#if !defined(__clang__) and defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
#endif
namespace rj = rapidjson;

namespace RDKit {
namespace MinimalLib {

static constexpr unsigned int d_defaultWidth = 250;
static constexpr unsigned int d_defaultHeight = 200;

RWMol *mol_from_input(const std::string &input) {
  RWMol *res = nullptr;
  if (input.find("M  END") != std::string::npos) {
    bool sanitize = false;
    bool removeHs = true;
    bool strictParsing = false;
    res = MolBlockToMol(input, sanitize, removeHs, strictParsing);
  } else if (input.find("commonchem") != std::string::npos) {
    auto molVect = MolInterchange::JSONDataToMols(input);
    if (!molVect.empty()) {
      res = new RWMol(*molVect[0]);
    }
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

RWMol *qmol_from_input(const std::string &input) {
  RWMol *res = nullptr;
  if (input.find("M  END") != std::string::npos) {
    bool sanitize = false;
    bool removeHs = true;
    bool strictParsing = false;
    res = MolBlockToMol(input, sanitize, removeHs, strictParsing);
  } else if (input.find("commonchem") != std::string::npos) {
    auto molVect = MolInterchange::JSONDataToMols(input);
    if (!molVect.empty()) {
      res = new RWMol(*molVect[0]);
    }
  } else {
    res = SmartsToMol(input);
  }
  return res;
}

std::string process_details(const std::string &details, unsigned int &width,
                            unsigned int &height, int &offsetx, int &offsety,
                            std::string &legend, std::vector<int> &atomIds,
                            std::vector<int> &bondIds) {
  rj::Document doc;
  doc.Parse(details.c_str());
  if (!doc.IsObject()) return "Invalid JSON";

  if (doc.HasMember("atoms")) {
    if (!doc["atoms"].IsArray()) {
      return "JSON doesn't contain 'atoms' field, or it is not an array";
    }
    for (const auto &molval : doc["atoms"].GetArray()) {
      if (!molval.IsInt()) return ("Atom IDs should be integers");
      atomIds.push_back(molval.GetInt());
    }
  }
  if (doc.HasMember("bonds")) {
    if (!doc["bonds"].IsArray()) {
      return "JSON contain 'bonds' field, but it is not an array";
    }
    for (const auto &molval : doc["bonds"].GetArray()) {
      if (!molval.IsInt()) return ("Bond IDs should be integers");
      bondIds.push_back(molval.GetInt());
    }
  }

  if (doc.HasMember("width")) {
    if (!doc["width"].IsUint()) {
      return "JSON contains 'width' field, but it is not an unsigned int";
    }
    width = doc["width"].GetUint();
  }

  if (doc.HasMember("height")) {
    if (!doc["height"].IsUint()) {
      return "JSON contains 'height' field, but it is not an unsigned int";
    }
    height = doc["height"].GetUint();
  }

  if (doc.HasMember("offsetx")) {
    if (!doc["offsetx"].IsInt()) {
      return "JSON contains 'offsetx' field, but it is not an int";
    }
    offsetx = doc["offsetx"].GetInt();
  }

  if (doc.HasMember("offsety")) {
    if (!doc["offsety"].IsInt()) {
      return "JSON contains 'offsety' field, but it is not an int";
    }
    offsety = doc["offsety"].GetInt();
  }

  if (doc.HasMember("legend")) {
    if (!doc["legend"].IsString()) {
      return "JSON contains 'legend' field, but it is not a string";
    }
    legend = doc["legend"].GetString();
  }

  return "";
}

void get_sss_json(const ROMol &d_mol, const ROMol &q_mol,
                  const MatchVectType &match, rj::Value &obj,
                  rj::Document &doc) {
  rj::Value rjAtoms(rj::kArrayType);
  for (const auto &pr : match) {
    rjAtoms.PushBack(pr.second, doc.GetAllocator());
  }
  obj.AddMember("atoms", rjAtoms, doc.GetAllocator());

  rj::Value rjBonds(rj::kArrayType);
  for (const auto qbond : q_mol.bonds()) {
    unsigned int beginIdx = qbond->getBeginAtomIdx();
    unsigned int endIdx = qbond->getEndAtomIdx();
    if (beginIdx >= match.size() || endIdx >= match.size()) {
      continue;
    }
    unsigned int idx1 = match[beginIdx].second;
    unsigned int idx2 = match[endIdx].second;
    const auto bond = d_mol.getBondBetweenAtoms(idx1, idx2);
    if (bond != nullptr) {
      rjBonds.PushBack(bond->getIdx(), doc.GetAllocator());
    }
  }
  obj.AddMember("bonds", rjBonds, doc.GetAllocator());
}

std::string mol_to_svg(const ROMol &m, unsigned int w, unsigned int h,
                       const std::string &details = "") {
  std::vector<int> atomIds;
  std::vector<int> bondIds;
  std::string legend = "";
  int offsetx = 0, offsety = 0;
  if (!details.empty()) {
    auto problems = process_details(details, w, h, offsetx, offsety, legend,
                                    atomIds, bondIds);
    if (!problems.empty()) {
      return problems;
    }
  }

  MolDraw2DSVG drawer(w, h);
  if (!details.empty()) {
    MolDraw2DUtils::updateDrawerParamsFromJSON(drawer, details);
  }
  drawer.setOffset(offsetx, offsety);

  MolDraw2DUtils::prepareAndDrawMolecule(drawer, m, legend, &atomIds, &bondIds);
  drawer.finishDrawing();

  return drawer.getDrawingText();
}

}  // namespace MinimalLib
}  // namespace RDKit