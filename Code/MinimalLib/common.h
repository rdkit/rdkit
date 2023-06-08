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
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/Descriptors/Property.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/MACCS.h>
#ifdef RDK_BUILD_AVALON_SUPPORT
#include <External/AvalonTools/AvalonTools.h>
#endif
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/Substruct/SubstructUtils.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Abbreviations/Abbreviations.h>
#include <DataStructs/BitOps.h>
#include <GraphMol/MolStandardize/MolStandardize.h>
#include <GraphMol/MolStandardize/Charge.h>
#include <GraphMol/MolStandardize/Tautomer.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/SanitizeRxn.h>

#include <sstream>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

#ifndef _MSC_VER
// shutoff some warnings from rapidjson
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif
#endif
#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>
#ifndef _MSC_VER
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
#endif

#define GET_JSON_VALUE(doc, key, type)                                     \
  const auto key##It = doc.FindMember(#key);                               \
  if (key##It != doc.MemberEnd()) {                                        \
    if (!key##It->value.Is##type()) {                                      \
      return "JSON contains '" #key "' field, but its type is not '" #type \
             "'";                                                          \
    }                                                                      \
    key = key##It->value.Get##type();                                      \
  }

#define GET_JSON_VALUE_WITH_DEFAULT(doc, key, type, defaultValue) \
  GET_JSON_VALUE(doc, key, type) else key = defaultValue;

namespace rj = rapidjson;

namespace RDKit {
namespace MinimalLib {

static constexpr int d_defaultWidth = 250;
static constexpr int d_defaultHeight = 200;

#define LPT_OPT_GET(opt) opt = pt.get(#opt, opt);
#define LPT_OPT_GET2(holder, opt) holder.opt = pt.get(#opt, holder.opt);

RWMol *mol_from_input(const std::string &input,
                      const std::string &details_json = "") {
  bool sanitize = true;
  bool kekulize = true;
  bool removeHs = true;
  bool mergeQueryHs = false;
  RWMol *res = nullptr;
  boost::property_tree::ptree pt;
  if (!details_json.empty()) {
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    LPT_OPT_GET(sanitize);
    LPT_OPT_GET(kekulize);
    LPT_OPT_GET(removeHs);
    LPT_OPT_GET(mergeQueryHs);
  }
  try {
    if (input.find("M  END") != std::string::npos) {
      bool strictParsing = false;
      LPT_OPT_GET(strictParsing);
      res = MolBlockToMol(input, false, removeHs, strictParsing);
    } else if (input.find("commonchem") != std::string::npos ||
               input.find("rdkitjson") != std::string::npos) {
      auto ps = MolInterchange::defaultJSONParseParameters;
      LPT_OPT_GET2(ps, setAromaticBonds);
      LPT_OPT_GET2(ps, strictValenceCheck);
      LPT_OPT_GET2(ps, parseProperties);
      LPT_OPT_GET2(ps, parseConformers);

      auto molVect = MolInterchange::JSONDataToMols(input, ps);
      if (!molVect.empty()) {
        res = new RWMol(*molVect[0]);
      }
    } else {
      SmilesParserParams ps;
      ps.sanitize = false;
      ps.removeHs = removeHs;
      LPT_OPT_GET2(ps, strictCXSMILES);
      res = SmilesToMol(input, ps);
    }
  } catch (...) {
    // we really don't want exceptions to be thrown in here
    res = nullptr;
  }
  if (res) {
    try {
      if (sanitize) {
        unsigned int failedOp;
        unsigned int sanitizeOps = MolOps::SANITIZE_ALL;
        if (!kekulize) {
          sanitizeOps ^= MolOps::SANITIZE_KEKULIZE;
        }
        MolOps::sanitizeMol(*res, failedOp, sanitizeOps);
      } else {
        res->updatePropertyCache(false);
      }
      MolOps::assignStereochemistry(*res, true, true, true);
      if (mergeQueryHs) {
        MolOps::mergeQueryHs(*res);
      }
    } catch (...) {
      delete res;
      res = nullptr;
    }
  }
  return res;
}

RWMol *mol_from_input(const std::string &input, const char *details_json) {
  std::string json;
  if (details_json) {
    json = details_json;
  }
  return mol_from_input(input, json);
}

RWMol *qmol_from_input(const std::string &input,
                       const std::string &details_json = "") {
  RWMol *res = nullptr;
  bool removeHs = true;
  boost::property_tree::ptree pt;
  if (!details_json.empty()) {
    // FIX: this should eventually be moved somewhere else
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    LPT_OPT_GET(removeHs);
  }
  if (input.find("M  END") != std::string::npos) {
    bool strictParsing = false;
    LPT_OPT_GET(strictParsing);
    res = MolBlockToMol(input, false, removeHs, strictParsing);
  } else if (input.find("commonchem") != std::string::npos ||
             input.find("rdkitjson") != std::string::npos) {
    auto ps = MolInterchange::defaultJSONParseParameters;
    LPT_OPT_GET2(ps, setAromaticBonds);
    LPT_OPT_GET2(ps, strictValenceCheck);
    LPT_OPT_GET2(ps, parseProperties);
    LPT_OPT_GET2(ps, parseConformers);

    auto molVect = MolInterchange::JSONDataToMols(input, ps);
    if (!molVect.empty()) {
      res = new RWMol(*molVect[0]);
    }
  } else {
    bool mergeHs = false;
    LPT_OPT_GET(mergeHs);
    res = SmartsToMol(input, 0, mergeHs);
  }
  return res;
}

RWMol *qmol_from_input(const std::string &input, const char *details_json) {
  std::string json;
  if (details_json) {
    json = details_json;
  }
  return qmol_from_input(input, json);
}

ChemicalReaction *rxn_from_input(const std::string &input,
                                 const std::string &details_json = "") {
  bool useSmiles = false;
  bool sanitize = false;
  ChemicalReaction *rxn = nullptr;
  boost::property_tree::ptree pt;
  if (!details_json.empty()) {
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    LPT_OPT_GET(sanitize);
    LPT_OPT_GET(useSmiles);
  }
  try {
    if (input.find("$RXN") != std::string::npos) {
      bool removeHs = false;
      bool strictParsing = false;
      LPT_OPT_GET(removeHs);
      LPT_OPT_GET(strictParsing);
      rxn = RxnBlockToChemicalReaction(input, false, removeHs, strictParsing);
    } else {
      rxn = RxnSmartsToChemicalReaction(input, nullptr, useSmiles);
    }
  } catch (...) {
    // we really don't want exceptions to be thrown in here
    rxn = nullptr;
  }
  if (rxn) {
    try {
      if (sanitize) {
        unsigned int failedOp;
        unsigned int sanitizeOps = RxnOps::SANITIZE_ALL;
        bool adjustReactants = true;
        bool mergeQueryHs = true;
        LPT_OPT_GET(adjustReactants);
        LPT_OPT_GET(mergeQueryHs);
        if (!adjustReactants) {
          sanitizeOps ^= RxnOps::SANITIZE_ADJUST_REACTANTS;
        }
        if (!mergeQueryHs) {
          sanitizeOps ^= RxnOps::SANITIZE_MERGEHS;
        }
        RxnOps::sanitizeRxn(*rxn, failedOp, sanitizeOps);
      }
    } catch (...) {
      delete rxn;
      rxn = nullptr;
    }
  }
  return rxn;
}

ChemicalReaction *rxn_from_input(const std::string &input,
                                 const char *details_json) {
  std::string json;
  if (details_json) {
    json = details_json;
  }
  return rxn_from_input(input, json);
}

std::string parse_int_array(const rj::Document &doc, std::vector<int> &intVec,
                            const std::string &keyName,
                            const std::string &valueName) {
  const auto it = doc.FindMember(keyName.c_str());
  if (it != doc.MemberEnd()) {
    if (!it->value.IsArray()) {
      return "JSON contains '" + keyName + "' field, but it is not an array";
    }
    for (const auto &val : it->value.GetArray()) {
      if (!val.IsInt()) {
        return valueName + " should be integers";
      }
      intVec.push_back(val.GetInt());
    }
  }
  return "";
}

std::string parse_rgba_array(const rj::Value &val, DrawColour &color,
                             const std::string &keyName) {
  if (!val.IsArray() || val.Size() < 3 || val.Size() > 4) {
    return "JSON contains '" + keyName +
           "' field, but the "
           "colors are not R,G,B[,A] arrays";
  }
  std::vector<double> rgba(4, 1.0);
  unsigned int i = 0;
  for (const auto &component : val.GetArray()) {
    if (!component.IsNumber()) {
      return "JSON contains '" + keyName +
             "' field, but the "
             "R,G,B[,A] arrays contain non-float values";
    }
    CHECK_INVARIANT(i < 4, "");
    rgba[i++] = component.GetDouble();
  }
  color.r = rgba[0];
  color.g = rgba[1];
  color.b = rgba[2];
  color.a = rgba[3];
  return "";
}

std::string parse_highlight_colors(const rj::Document &doc,
                                   std::map<int, DrawColour> &colorMap,
                                   const std::string &keyName) {
  const auto it = doc.FindMember(keyName.c_str());
  if (it != doc.MemberEnd()) {
    if (!it->value.IsObject()) {
      return "JSON contains '" + keyName + "' field, but it is not an object";
    }
    for (const auto &entry : it->value.GetObject()) {
      DrawColour color;
      auto problems = parse_rgba_array(entry.value, color, keyName);
      if (!problems.empty()) {
        return problems;
      }
      int idx = std::atoi(entry.name.GetString());
      colorMap[idx] = std::move(color);
    }
  }
  return "";
}

std::string process_details(rj::Document &doc, const std::string &details,
                            int &width, int &height, int &offsetx, int &offsety,
                            std::string &legend, std::vector<int> &atomIds,
                            std::vector<int> &bondIds, bool &kekulize,
                            bool &addChiralHs, bool &wedgeBonds,
                            bool &forceCoords, bool &wavyBonds) {
  doc.Parse(details.c_str());
  if (!doc.IsObject()) {
    return "Invalid JSON";
  }
  std::string problems;
  problems = parse_int_array(doc, atomIds, "atoms", "Atom IDs");
  if (!problems.empty()) {
    return problems;
  }

  problems = parse_int_array(doc, bondIds, "bonds", "Bond IDs");
  if (!problems.empty()) {
    return problems;
  }

  GET_JSON_VALUE(doc, width, Int)
  GET_JSON_VALUE(doc, height, Int)
  GET_JSON_VALUE(doc, offsetx, Int)
  GET_JSON_VALUE(doc, offsety, Int)
  GET_JSON_VALUE(doc, legend, String)
  GET_JSON_VALUE_WITH_DEFAULT(doc, kekulize, Bool, true)
  GET_JSON_VALUE_WITH_DEFAULT(doc, addChiralHs, Bool, true)
  GET_JSON_VALUE_WITH_DEFAULT(doc, wedgeBonds, Bool, true)
  GET_JSON_VALUE_WITH_DEFAULT(doc, forceCoords, Bool, false)
  GET_JSON_VALUE_WITH_DEFAULT(doc, wavyBonds, Bool, false)

  return "";
}

std::string process_mol_details(const std::string &details, int &width,
                                int &height, int &offsetx, int &offsety,
                                std::string &legend, std::vector<int> &atomIds,
                                std::vector<int> &bondIds,
                                std::map<int, DrawColour> &atomMap,
                                std::map<int, DrawColour> &bondMap,
                                std::map<int, double> &radiiMap, bool &kekulize,
                                bool &addChiralHs, bool &wedgeBonds,
                                bool &forceCoords, bool &wavyBonds) {
  rj::Document doc;
  auto problems = process_details(
      doc, details, width, height, offsetx, offsety, legend, atomIds, bondIds,
      kekulize, addChiralHs, wedgeBonds, forceCoords, wavyBonds);
  if (!problems.empty()) {
    return problems;
  }

  problems = parse_highlight_colors(doc, atomMap, "highlightAtomColors");
  if (!problems.empty()) {
    return problems;
  }

  problems = parse_highlight_colors(doc, bondMap, "highlightBondColors");
  if (!problems.empty()) {
    return problems;
  }

  const auto radiiMapit = doc.FindMember("highlightAtomRadii");
  if (radiiMapit != doc.MemberEnd()) {
    if (!radiiMapit->value.IsObject()) {
      return "JSON contains 'highlightAtomRadii' field, but it is not an object";
    }
    for (const auto &entry : radiiMapit->value.GetObject()) {
      if (!entry.value.IsNumber()) {
        return "JSON contains 'highlightAtomRadii' field, but the radii"
               "are not floats";
      }
      int idx = std::atoi(entry.name.GetString());
      radiiMap[idx] = entry.value.GetDouble();
    }
  }
  return "";
}

std::string process_rxn_details(
    const std::string &details, int &width, int &height, int &offsetx,
    int &offsety, std::string &legend, std::vector<int> &atomIds,
    std::vector<int> &bondIds, bool &kekulize, bool &highlightByReactant,
    std::vector<DrawColour> &highlightColorsReactants) {
  rj::Document doc;
  bool addChiralHs;
  bool wedgeBonds;
  bool forceCoords;
  bool wavyBonds;
  auto problems = process_details(
      doc, details, width, height, offsetx, offsety, legend, atomIds, bondIds,
      kekulize, addChiralHs, wedgeBonds, forceCoords, wavyBonds);
  if (!problems.empty()) {
    return problems;
  }
  GET_JSON_VALUE_WITH_DEFAULT(doc, highlightByReactant, Bool, false)
  auto highlightColorsReactantsIt = doc.FindMember("highlightColorsReactants");
  if (highlightColorsReactantsIt != doc.MemberEnd()) {
    if (!highlightColorsReactantsIt->value.IsArray()) {
      return "JSON contains 'highlightColorsReactants' field, but it is not an "
             "array";
    }
    for (const auto &rgbaArray : highlightColorsReactantsIt->value.GetArray()) {
      DrawColour color;
      problems = parse_rgba_array(rgbaArray, color, "highlightColorsReactants");
      if (!problems.empty()) {
        return problems;
      }
      highlightColorsReactants.push_back(std::move(color));
    }
  }
  return "";
}

std::string molblock_helper(RWMol &mol, const char *details_json,
                            bool forceV3000) {
  bool includeStereo = true;
  bool kekulize = true;
  bool useMolBlockWedging = false;
  bool addChiralHs = false;
  if (details_json && strlen(details_json)) {
    boost::property_tree::ptree pt;
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    LPT_OPT_GET(includeStereo);
    LPT_OPT_GET(kekulize);
    LPT_OPT_GET(useMolBlockWedging);
    LPT_OPT_GET(addChiralHs);
  }
  if (useMolBlockWedging) {
    reapplyMolBlockWedging(mol);
  }
  if (addChiralHs) {
    MolDraw2DUtils::prepareMolForDrawing(mol, false, true, false, false, false);
  }
  return MolToMolBlock(mol, includeStereo, -1, kekulize, forceV3000);
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

std::string mol_to_svg(const ROMol &m, int w, int h,
                       const std::string &details = "") {
  std::vector<int> atomIds;
  std::vector<int> bondIds;
  std::map<int, DrawColour> atomMap;
  std::map<int, DrawColour> bondMap;
  std::map<int, double> radiiMap;
  std::string legend = "";
  std::string problems;
  int offsetx = 0;
  int offsety = 0;
  bool kekulize = true;
  bool addChiralHs = true;
  bool wedgeBonds = true;
  bool forceCoords = false;
  bool wavyBonds = false;
  if (!details.empty()) {
    problems =
        process_mol_details(details, w, h, offsetx, offsety, legend, atomIds,
                            bondIds, atomMap, bondMap, radiiMap, kekulize,
                            addChiralHs, wedgeBonds, forceCoords, wavyBonds);
    if (!problems.empty()) {
      return problems;
    }
  }
  MolDraw2DSVG drawer(w, h);
  if (!details.empty()) {
    MolDraw2DUtils::updateDrawerParamsFromJSON(drawer, details);
  }
  drawer.setOffset(offsetx, offsety);

  MolDraw2DUtils::prepareAndDrawMolecule(drawer, m, legend, &atomIds, &bondIds,
                                         atomMap.empty() ? nullptr : &atomMap,
                                         bondMap.empty() ? nullptr : &bondMap,
                                         radiiMap.empty() ? nullptr : &radiiMap,
                                         -1, kekulize, addChiralHs, wedgeBonds,
                                         forceCoords, wavyBonds);
  drawer.finishDrawing();

  return drawer.getDrawingText();
}

std::string rxn_to_svg(const ChemicalReaction &rxn, int w, int h,
                       const std::string &details = "") {
  std::vector<int> atomIds;
  std::vector<int> bondIds;
  std::string legend = "";
  int offsetx = 0;
  int offsety = 0;
  bool kekulize = true;
  bool highlightByReactant = false;
  std::vector<DrawColour> highlightColorsReactants;
  if (!details.empty()) {
    auto problems = process_rxn_details(
        details, w, h, offsetx, offsety, legend, atomIds, bondIds, kekulize,
        highlightByReactant, highlightColorsReactants);
    if (!problems.empty()) {
      return problems;
    }
  }

  MolDraw2DSVG drawer(w, h);
  if (!kekulize) {
    drawer.drawOptions().prepareMolsBeforeDrawing = false;
  }
  drawer.drawReaction(rxn, highlightByReactant,
                      !highlightByReactant || highlightColorsReactants.empty()
                          ? nullptr
                          : &highlightColorsReactants);
  drawer.finishDrawing();
  return drawer.getDrawingText();
}

std::string get_descriptors(const ROMol &m) {
  rj::Document doc;
  doc.SetObject();

  Descriptors::Properties props;
  std::vector<std::string> dns = props.getPropertyNames();
  std::vector<double> dvs = props.computeProperties(m);
  for (size_t i = 0; i < dns.size(); ++i) {
    rj::Value v(dvs[i]);
    const auto srt = rj::StringRef(dns[i].c_str());
    doc.AddMember(srt, v, doc.GetAllocator());
  }

  rj::StringBuffer buffer;
  rj::Writer<rj::StringBuffer> writer(buffer);
  writer.SetMaxDecimalPlaces(5);
  doc.Accept(writer);
  return buffer.GetString();
}

namespace {
template <typename T, typename U>
std::unique_ptr<RWMol> standardize_func(T &mol, const std::string &details_json,
                                        U func) {
  MolStandardize::CleanupParameters ps =
      MolStandardize::defaultCleanupParameters;
  if (!details_json.empty()) {
    MolStandardize::updateCleanupParamsFromJSON(ps, details_json);
  }
  return std::unique_ptr<RWMol>(static_cast<RWMol *>(func(mol, ps)));
}

bool invertWedgingIfMolHasFlipped(ROMol &mol,
                                  const RDGeom::Transform3D &trans) {
  constexpr double FLIP_THRESHOLD = -0.99;
  auto zRot = trans.getVal(2, 2);
  bool shouldFlip = zRot < FLIP_THRESHOLD;
  if (shouldFlip) {
    invertMolBlockWedgingInfo(mol);
  }
  return shouldFlip;
}
}  // namespace

std::unique_ptr<RWMol> do_cleanup(RWMol &mol, const std::string &details_json) {
  auto molp = &mol;
  return standardize_func(
      molp, details_json,
      static_cast<RWMol *(*)(const RWMol *,
                             const MolStandardize::CleanupParameters &)>(
          MolStandardize::cleanup));
}

std::unique_ptr<RWMol> do_normalize(RWMol &mol,
                                    const std::string &details_json) {
  auto molp = &mol;
  return standardize_func(molp, details_json, MolStandardize::normalize);
}

std::unique_ptr<RWMol> do_reionize(RWMol &mol,
                                   const std::string &details_json) {
  auto molp = &mol;
  return standardize_func(molp, details_json, MolStandardize::reionize);
}

std::unique_ptr<RWMol> do_canonical_tautomer(RWMol &mol,
                                             const std::string &details_json) {
  MolStandardize::CleanupParameters ps =
      MolStandardize::defaultCleanupParameters;
  if (!details_json.empty()) {
    MolStandardize::updateCleanupParamsFromJSON(ps, details_json);
  }
  MolStandardize::TautomerEnumerator te(ps);
  std::unique_ptr<RWMol> res(static_cast<RWMol *>(te.canonicalize(mol)));
  return res;
}

std::unique_ptr<RWMol> do_neutralize(RWMol &mol,
                                     const std::string &details_json) {
  MolStandardize::CleanupParameters ps =
      MolStandardize::defaultCleanupParameters;
  if (!details_json.empty()) {
    MolStandardize::updateCleanupParamsFromJSON(ps, details_json);
  }
  MolStandardize::Uncharger uncharger(ps.doCanonical);
  std::unique_ptr<RWMol> res(static_cast<RWMol *>(uncharger.uncharge(mol)));
  return res;
}

std::unique_ptr<RWMol> do_charge_parent(RWMol &mol,
                                        const std::string &details_json) {
  MolStandardize::CleanupParameters ps =
      MolStandardize::defaultCleanupParameters;
  bool skipStandardize = false;
  if (!details_json.empty()) {
    MolStandardize::updateCleanupParamsFromJSON(ps, details_json);
    boost::property_tree::ptree pt;
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    LPT_OPT_GET(skipStandardize);
  }
  std::unique_ptr<RWMol> res(
      MolStandardize::chargeParent(mol, ps, skipStandardize));
  return res;
}

std::unique_ptr<RWMol> do_fragment_parent(RWMol &mol,
                                          const std::string &details_json) {
  MolStandardize::CleanupParameters ps =
      MolStandardize::defaultCleanupParameters;
  bool skipStandardize = false;
  if (!details_json.empty()) {
    MolStandardize::updateCleanupParamsFromJSON(ps, details_json);
    boost::property_tree::ptree pt;
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    LPT_OPT_GET(skipStandardize);
  }
  std::unique_ptr<RWMol> res(
      MolStandardize::fragmentParent(mol, ps, skipStandardize));
  return res;
}

std::unique_ptr<ExplicitBitVect> morgan_fp_as_bitvect(
    const RWMol &mol, const char *details_json) {
  size_t radius = 2;
  size_t nBits = 2048;
  bool useChirality = false;
  bool useBondTypes = true;
  bool includeRedundantEnvironments = false;
  bool onlyNonzeroInvariants = false;
  if (details_json && strlen(details_json)) {
    // FIX: this should eventually be moved somewhere else
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(ss, pt);
    LPT_OPT_GET(radius);
    LPT_OPT_GET(nBits);
    LPT_OPT_GET(useChirality);
    LPT_OPT_GET(useBondTypes);
    LPT_OPT_GET(includeRedundantEnvironments);
    LPT_OPT_GET(onlyNonzeroInvariants);
  }
  auto fp = MorganFingerprints::getFingerprintAsBitVect(
      mol, radius, nBits, nullptr, nullptr, useChirality, useBondTypes,
      onlyNonzeroInvariants, nullptr, includeRedundantEnvironments);
  return std::unique_ptr<ExplicitBitVect>{fp};
}

std::unique_ptr<ExplicitBitVect> rdkit_fp_as_bitvect(const RWMol &mol,
                                                     const char *details_json) {
  unsigned int minPath = 1;
  unsigned int maxPath = 7;
  unsigned int nBits = 2048;
  unsigned int nBitsPerHash = 2;
  bool useHs = true;
  bool branchedPaths = true;
  bool useBondOrder = true;
  if (details_json && strlen(details_json)) {
    // FIX: this should eventually be moved somewhere else
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(ss, pt);
    LPT_OPT_GET(minPath);
    LPT_OPT_GET(maxPath);
    LPT_OPT_GET(nBits);
    LPT_OPT_GET(nBitsPerHash);
    LPT_OPT_GET(useHs);
    LPT_OPT_GET(branchedPaths);
    LPT_OPT_GET(useBondOrder);
  }
  auto fp = RDKFingerprintMol(mol, minPath, maxPath, nBits, nBitsPerHash, useHs,
                              0, 128, branchedPaths, useBondOrder);
  return std::unique_ptr<ExplicitBitVect>{fp};
}

std::unique_ptr<ExplicitBitVect> pattern_fp_as_bitvect(
    const RWMol &mol, const char *details_json) {
  unsigned int nBits = 2048;
  bool tautomericFingerprint = false;
  if (details_json && strlen(details_json)) {
    // FIX: this should eventually be moved somewhere else
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(ss, pt);
    LPT_OPT_GET(nBits);
    LPT_OPT_GET(tautomericFingerprint);
  }
  auto fp = PatternFingerprintMol(mol, nBits, nullptr, nullptr,
                                  tautomericFingerprint);
  return std::unique_ptr<ExplicitBitVect>{fp};
}

std::unique_ptr<ExplicitBitVect> topological_torsion_fp_as_bitvect(
    const RWMol &mol, const char *details_json) {
  unsigned int nBits = 2048;
  if (details_json && strlen(details_json)) {
    // FIX: this should eventually be moved somewhere else
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(ss, pt);
    LPT_OPT_GET(nBits);
  }
  auto fp =
      AtomPairs::getHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits);
  return std::unique_ptr<ExplicitBitVect>{fp};
}

std::unique_ptr<ExplicitBitVect> atom_pair_fp_as_bitvect(
    const RWMol &mol, const char *details_json) {
  unsigned int nBits = 2048;
  unsigned int minLength = 1;
  unsigned int maxLength = 30;
  if (details_json && strlen(details_json)) {
    // FIX: this should eventually be moved somewhere else
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(ss, pt);
    LPT_OPT_GET(nBits);
    LPT_OPT_GET(minLength);
    LPT_OPT_GET(maxLength);
  }
  auto fp = AtomPairs::getHashedAtomPairFingerprintAsBitVect(
      mol, nBits, minLength, maxLength);
  return std::unique_ptr<ExplicitBitVect>{fp};
}

std::unique_ptr<ExplicitBitVect> maccs_fp_as_bitvect(const RWMol &mol) {
  auto fp = MACCSFingerprints::getFingerprintAsBitVect(mol);
  return std::unique_ptr<ExplicitBitVect>{fp};
}

#ifdef RDK_BUILD_AVALON_SUPPORT
std::unique_ptr<ExplicitBitVect> avalon_fp_as_bitvect(
    const RWMol &mol, const char *details_json) {
  unsigned int nBits = 512;
  if (details_json && strlen(details_json)) {
    // FIX: this should eventually be moved somewhere else
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(ss, pt);
    LPT_OPT_GET(nBits);
  }
  std::unique_ptr<ExplicitBitVect> fp(new ExplicitBitVect(nBits));
  AvalonTools::getAvalonFP(mol, *fp, nBits);
  return fp;
}
#endif

// If alignOnly is set to true in details_json, original molblock wedging
// information is preserved, and inverted if needed (in case the rigid-body
// alignment required a flip around the Z axis).
// If alignOnly is set to false in details_json or not specified, original
// molblock wedging information is preserved if it only involves the invariant
// core whose coordinates never change, and is cleared in case coordinates were
// changed. If acceptFailure is set to true and no substructure match is found,
// coordinates will be recomputed from scratch, hence molblock wedging
// information will be cleared.
std::string generate_aligned_coords(ROMol &mol, const ROMol &templateMol,
                                    const char *details_json) {
  std::string res;
  if (!templateMol.getNumConformers()) {
    return res;
  }
  constexpr int MAX_MATCHES = 1000;
  constexpr double RMSD_THRESHOLD = 1.e-2;
  constexpr double MSD_THRESHOLD = RMSD_THRESHOLD * RMSD_THRESHOLD;
  bool useCoordGen = false;
  bool allowRGroups = false;
  bool acceptFailure = true;
  bool alignOnly = false;
  if (details_json && strlen(details_json)) {
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(ss, pt);
    LPT_OPT_GET(useCoordGen);
    LPT_OPT_GET(allowRGroups);
    LPT_OPT_GET(acceptFailure);
    LPT_OPT_GET(alignOnly);
  }
  MatchVectType match;
  int confId = -1;
  std::unique_ptr<Conformer> origConformer;
  std::unique_ptr<ROMol> templateMolAdj;
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  bool oprefer = RDDepict::preferCoordGen;
  RDDepict::preferCoordGen = useCoordGen;
#endif
  // store the original conformer so it can be restored
  // if alignment fails and acceptFailure is false
  if (!acceptFailure && mol.getNumConformers()) {
    origConformer.reset(new Conformer(mol.getConformer()));
  }
  if (alignOnly) {
    std::vector<MatchVectType> matches;
    std::unique_ptr<ROMol> molHs;
    ROMol *prbMol = &mol;
    if (allowRGroups) {
      auto atoms = templateMol.atoms();
      allowRGroups = std::any_of(atoms.begin(), atoms.end(),
                                 isAtomTerminalRGroupOrQueryHydrogen);
    }
    if (allowRGroups) {
      molHs.reset(MolOps::addHs(mol));
      prbMol = molHs.get();
      auto queryParams = RDKit::MolOps::AdjustQueryParameters::noAdjustments();
      queryParams.adjustSingleBondsToDegreeOneNeighbors = true;
      queryParams.adjustSingleBondsBetweenAromaticAtoms = true;
      templateMolAdj.reset(
          RDKit::MolOps::adjustQueryProperties(templateMol, &queryParams));
    }
    const ROMol &templateMolRef =
        templateMolAdj ? *templateMolAdj : templateMol;
    if (SubstructMatch(*prbMol, templateMolRef, matches, false)) {
      if (allowRGroups) {
        matches = sortMatchesByDegreeOfCoreSubstitution(*prbMol, templateMolRef,
                                                        matches);
        int maxMatchedHeavies = -1;
        std::vector<MatchVectType> prunedMatches;
        prunedMatches.reserve(matches.size());
        for (const auto &match : matches) {
          int nMatchedHeavies = 0;
          MatchVectType prunedMatch;
          prunedMatch.reserve(match.size());
          for (const auto &pair : match) {
            const auto templateAtom = templateMolRef.getAtomWithIdx(pair.first);
            const auto prbAtom = prbMol->getAtomWithIdx(pair.second);
            if (isAtomTerminalRGroupOrQueryHydrogen(templateAtom)) {
              if (prbAtom->getAtomicNum() == 1) {
                continue;
              }
              ++nMatchedHeavies;
            }
            prunedMatch.push_back(std::move(pair));
          }
          if (nMatchedHeavies < maxMatchedHeavies) {
            break;
          } else {
            prunedMatches.push_back(std::move(prunedMatch));
            maxMatchedHeavies = nMatchedHeavies;
          }
        }
        matches = std::move(prunedMatches);
      }
      std::for_each(matches.begin(), matches.end(), [](auto &match) {
        std::for_each(match.begin(), match.end(),
                      [](auto &pair) { std::swap(pair.first, pair.second); });
      });
      if (!mol.getNumConformers()) {
        RDDepict::compute2DCoords(mol);
      }
      RDGeom::Transform3D trans;
      MolAlign::getBestAlignmentTransform(mol, templateMolRef, trans, match,
                                          confId, confId, matches, MAX_MATCHES);
      std::for_each(match.begin(), match.end(),
                    [](auto &pair) { std::swap(pair.first, pair.second); });
      MolTransforms::transformConformer(mol.getConformer(), trans);
      invertWedgingIfMolHasFlipped(mol, trans);
    } else if (acceptFailure) {
      RDDepict::compute2DCoords(mol);
      clearMolBlockWedgingInfo(mol);
    }
  } else {
    const ROMol *refPattern = nullptr;
    // always accept failure in the original call because
    // we detect it afterwards and, in case, restore the
    // original conformation
    const bool acceptOrigFailure = true;
    std::unique_ptr<ROMol> molOrig;
    if (mol.getNumConformers()) {
      molOrig.reset(new ROMol(mol));
    }
    match = RDDepict::generateDepictionMatching2DStructure(
        mol, templateMol, confId, refPattern, acceptOrigFailure, false,
        allowRGroups);
    // we need to clear the existing wedging information if:
    // 1. there is no match and we accept failure; this means that
    //    we rebuild coordinates from scratch, hence pre-existing
    //    wedging info is not valid anymore
    // 2. there is a match
    // 3. the original molecule has no coordinates to start with
    //    (in that case it should already have no wedging info either, anyway)
    // If there is no match and we do not accept failure, we keep
    // existing coordinates and hence also keep the wedging info
    bool shouldNeverClearWedgingInfo = match.empty() && !acceptFailure;
    bool shouldClearWedgingInfo = (match.empty() && acceptFailure) || !molOrig;
    if (!shouldNeverClearWedgingInfo) {
      if (!shouldClearWedgingInfo) {
        std::set<unsigned int> molMatchIndices;
        std::transform(match.begin(), match.end(),
                       std::inserter(molMatchIndices, molMatchIndices.begin()),
                       [](const auto &pair) { return pair.second; });
        // if any of the bonds that have wedging information from the molblock
        // has at least one atom which is not part of the scaffold, we cannot
        // preserve wedging information
        auto molBonds = mol.bonds();
        shouldClearWedgingInfo = std::any_of(
            molBonds.begin(), molBonds.end(), [&molMatchIndices](const auto b) {
              return ((b->hasProp(common_properties::_MolFileBondStereo) ||
                       b->hasProp(common_properties::_MolFileBondCfg)) &&
                      (!molMatchIndices.count(b->getBeginAtomIdx()) ||
                       !molMatchIndices.count(b->getEndAtomIdx())));
            });
      }
      if (!shouldClearWedgingInfo) {
        // check that scaffold coordinates have not changed, which may
        // happen when using CoordGen
        const auto &molPos = mol.getConformer().getPositions();
        const auto &templatePos = templateMol.getConformer().getPositions();
        shouldClearWedgingInfo = std::any_of(
            match.begin(), match.end(),
            [&molPos, &templatePos, MSD_THRESHOLD](const auto &pair) {
              return (molPos.at(pair.second) - templatePos.at(pair.first))
                         .lengthSq() > MSD_THRESHOLD;
            });
      }
      // final check: we still might need to invert wedging if the molecule
      // has flipped to match the scaffold
      if (!shouldClearWedgingInfo) {
        RDGeom::Transform3D trans;
        MatchVectType identityMatch(match.size());
        std::transform(match.begin(), match.end(), identityMatch.begin(),
                       [](const auto &pair) {
                         return std::make_pair(pair.second, pair.second);
                       });
        auto rmsd = MolAlign::getAlignmentTransform(
            *molOrig, mol, trans, confId, confId, &identityMatch);
        // this should not happen as we checked that previously, but we are
        // notoriously paranoid
        if (rmsd > RMSD_THRESHOLD) {
          shouldClearWedgingInfo = true;
        } else {
          invertWedgingIfMolHasFlipped(mol, trans);
        }
      }
    }
    if (shouldClearWedgingInfo) {
      clearMolBlockWedgingInfo(mol);
    }
  }
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = oprefer;
#endif
  if (match.empty()) {
    if (acceptFailure) {
      res = "{}";
    } else {
      if (mol.getNumConformers()) {
        mol.removeConformer(mol.getConformer().getId());
      }
      if (origConformer) {
        mol.addConformer(origConformer.release());
      }
      res = "";
    }
  } else {
    rj::Document doc;
    doc.SetObject();
    MinimalLib::get_sss_json(mol, templateMol, match, doc, doc);
    rj::StringBuffer buffer;
    rj::Writer<rj::StringBuffer> writer(buffer);
    doc.Accept(writer);
    res = buffer.GetString();
  }
  return res;
}

void get_mol_frags_details(const std::string &details_json, bool &sanitizeFrags,
                           bool &copyConformers) {
  boost::property_tree::ptree pt;
  if (!details_json.empty()) {
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    LPT_OPT_GET(sanitizeFrags);
    LPT_OPT_GET(copyConformers);
  }
}

std::string get_mol_frags_mappings(
    const std::vector<int> &frags,
    const std::vector<std::vector<int>> &fragsMolAtomMapping) {
  rj::Document doc;
  doc.SetObject();
  rj::Value rjFrags(rj::kArrayType);
  for (int fragIdx : frags) {
    rjFrags.PushBack(fragIdx, doc.GetAllocator());
  }
  doc.AddMember("frags", rjFrags, doc.GetAllocator());
  rj::Value rjFragsMolAtomMapping(rj::kArrayType);
  for (const auto &atomIdxVec : fragsMolAtomMapping) {
    rj::Value rjAtomIndices(rj::kArrayType);
    for (int atomIdx : atomIdxVec) {
      rjAtomIndices.PushBack(atomIdx, doc.GetAllocator());
    }
    rjFragsMolAtomMapping.PushBack(rjAtomIndices, doc.GetAllocator());
  }
  doc.AddMember("fragsMolAtomMapping", rjFragsMolAtomMapping,
                doc.GetAllocator());
  rj::StringBuffer buffer;
  rj::Writer<rj::StringBuffer> writer(buffer);
  doc.Accept(writer);
  return buffer.GetString();
}

}  // namespace MinimalLib
}  // namespace RDKit
#undef LPT_OPT_GET
#undef LPT_OPT_GET2
