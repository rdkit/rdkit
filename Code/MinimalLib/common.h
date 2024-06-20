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
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include <string>
#include <RDGeneral/versions.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
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
#include <GraphMol/Depictor/DepictUtils.h>
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
#include <GraphMol/RGroupDecomposition/RGroupUtils.h>
#include <RDGeneral/RDLog.h>
#include "common_defs.h"

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

#define GET_JSON_VALUE(doc, drawingDetails, key, type)                     \
  const auto key##It = doc.FindMember(#key);                               \
  if (key##It != doc.MemberEnd()) {                                        \
    if (!key##It->value.Is##type()) {                                      \
      return "JSON contains '" #key "' field, but its type is not '" #type \
             "'";                                                          \
    }                                                                      \
    drawingDetails.key = key##It->value.Get##type();                       \
  }

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
  bool setAromaticity = true;
  bool fastFindRings = true;
  bool assignStereo = true;
  bool assignCIPLabels = false;
  bool mappedDummiesAreRGroups = false;
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
    LPT_OPT_GET(setAromaticity);
    LPT_OPT_GET(fastFindRings);
    LPT_OPT_GET(assignStereo);
    LPT_OPT_GET(assignCIPLabels);
    LPT_OPT_GET(mappedDummiesAreRGroups);
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
        if (!setAromaticity) {
          sanitizeOps ^= MolOps::SANITIZE_SETAROMATICITY;
        }
        MolOps::sanitizeMol(*res, failedOp, sanitizeOps);
      } else {
        res->updatePropertyCache(false);
        if (fastFindRings) {
          MolOps::fastFindRings(*res);
        }
      }
      if (assignStereo) {
        MolOps::assignStereochemistry(*res, true, true, true);
      }
      if (assignCIPLabels) {
        CIPLabeler::assignCIPLabels(*res);
      }
      if (mergeQueryHs) {
        MolOps::mergeQueryHs(*res);
      }
      if (mappedDummiesAreRGroups) {
        relabelMappedDummies(*res);
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
                            DrawingDetails &drawingDetails) {
  doc.Parse(details.c_str());
  if (!doc.IsObject()) {
    return "Invalid JSON";
  }
  std::string problems;
  problems = parse_int_array(doc, drawingDetails.atomIds, "atoms", "Atom IDs");
  if (!problems.empty()) {
    return problems;
  }

  problems = parse_int_array(doc, drawingDetails.bondIds, "bonds", "Bond IDs");
  if (!problems.empty()) {
    return problems;
  }

  GET_JSON_VALUE(doc, drawingDetails, width, Int)
  GET_JSON_VALUE(doc, drawingDetails, height, Int)
  GET_JSON_VALUE(doc, drawingDetails, offsetx, Int)
  GET_JSON_VALUE(doc, drawingDetails, offsety, Int)
  GET_JSON_VALUE(doc, drawingDetails, panelWidth, Int)
  GET_JSON_VALUE(doc, drawingDetails, panelHeight, Int)
  GET_JSON_VALUE(doc, drawingDetails, noFreetype, Bool)
  GET_JSON_VALUE(doc, drawingDetails, legend, String)
  GET_JSON_VALUE(doc, drawingDetails, kekulize, Bool)
  GET_JSON_VALUE(doc, drawingDetails, addChiralHs, Bool)
  GET_JSON_VALUE(doc, drawingDetails, wedgeBonds, Bool)
  GET_JSON_VALUE(doc, drawingDetails, forceCoords, Bool)
  GET_JSON_VALUE(doc, drawingDetails, wavyBonds, Bool)
  GET_JSON_VALUE(doc, drawingDetails, useMolBlockWedging, Bool)

  return "";
}

[[deprecated(
    "please use the overload taking DrawingDetails& as parameter")]] std::string
process_details(rj::Document &doc, const std::string &details, int &width,
                int &height, int &offsetx, int &offsety, std::string &legend,
                std::vector<int> &atomIds, std::vector<int> &bondIds,
                bool &kekulize, bool &addChiralHs, bool &wedgeBonds,
                bool &forceCoords, bool &wavyBonds) {
  DrawingDetails drawingDetails;
  auto problems = process_details(doc, details, drawingDetails);
  width = drawingDetails.width;
  height = drawingDetails.height;
  offsetx = drawingDetails.offsetx;
  offsety = drawingDetails.offsety;
  legend = drawingDetails.legend;
  atomIds = drawingDetails.atomIds;
  bondIds = drawingDetails.bondIds;
  kekulize = drawingDetails.kekulize;
  addChiralHs = drawingDetails.addChiralHs;
  wedgeBonds = drawingDetails.wedgeBonds;
  forceCoords = drawingDetails.forceCoords;
  wavyBonds = drawingDetails.wavyBonds;
  return problems;
}

std::string process_mol_details(const std::string &details,
                                MolDrawingDetails &molDrawingDetails) {
  rj::Document doc;
  auto problems = process_details(doc, details, molDrawingDetails);
  if (!problems.empty()) {
    return problems;
  }

  problems = parse_highlight_colors(doc, molDrawingDetails.atomMap,
                                    "highlightAtomColors");
  if (!problems.empty()) {
    return problems;
  }

  problems = parse_highlight_colors(doc, molDrawingDetails.bondMap,
                                    "highlightBondColors");
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
      molDrawingDetails.radiiMap[idx] = entry.value.GetDouble();
    }
  }
  return "";
}

[[deprecated(
    "please use the overload taking MolDrawingDetails& as parameter")]] std::
    string
    process_mol_details(const std::string &details, int &width, int &height,
                        int &offsetx, int &offsety, std::string &legend,
                        std::vector<int> &atomIds, std::vector<int> &bondIds,
                        std::map<int, DrawColour> &atomMap,
                        std::map<int, DrawColour> &bondMap,
                        std::map<int, double> &radiiMap, bool &kekulize,
                        bool &addChiralHs, bool &wedgeBonds, bool &forceCoords,
                        bool &wavyBonds) {
  MolDrawingDetails molDrawingDetails;
  auto problems = process_mol_details(details, molDrawingDetails);
  width = molDrawingDetails.width;
  height = molDrawingDetails.height;
  offsetx = molDrawingDetails.offsetx;
  offsety = molDrawingDetails.offsety;
  legend = molDrawingDetails.legend;
  atomIds = molDrawingDetails.atomIds;
  bondIds = molDrawingDetails.bondIds;
  atomMap = molDrawingDetails.atomMap;
  bondMap = molDrawingDetails.bondMap;
  radiiMap = molDrawingDetails.radiiMap;
  kekulize = molDrawingDetails.kekulize;
  addChiralHs = molDrawingDetails.addChiralHs;
  wedgeBonds = molDrawingDetails.wedgeBonds;
  forceCoords = molDrawingDetails.forceCoords;
  wavyBonds = molDrawingDetails.wavyBonds;
  return problems;
}

std::string process_rxn_details(const std::string &details,
                                RxnDrawingDetails &rxnDrawingDetails) {
  rj::Document doc;
  auto problems = process_details(doc, details, rxnDrawingDetails);
  if (!problems.empty()) {
    return problems;
  }
  GET_JSON_VALUE(doc, rxnDrawingDetails, highlightByReactant, Bool)
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
      rxnDrawingDetails.highlightColorsReactants.push_back(std::move(color));
    }
  }
  return "";
}

[[deprecated(
    "please use the overload taking RxnDrawingDetails& as parameter")]] std::
    string
    process_rxn_details(const std::string &details, int &width, int &height,
                        int &offsetx, int &offsety, std::string &legend,
                        std::vector<int> &atomIds, std::vector<int> &bondIds,
                        bool &kekulize, bool &highlightByReactant,
                        std::vector<DrawColour> &highlightColorsReactants) {
  RxnDrawingDetails rxnDrawingDetails;
  auto problems = process_rxn_details(details, rxnDrawingDetails);
  width = rxnDrawingDetails.width;
  height = rxnDrawingDetails.height;
  offsetx = rxnDrawingDetails.offsetx;
  offsety = rxnDrawingDetails.offsety;
  legend = rxnDrawingDetails.legend;
  atomIds = rxnDrawingDetails.atomIds;
  bondIds = rxnDrawingDetails.bondIds;
  kekulize = rxnDrawingDetails.kekulize;
  highlightByReactant = rxnDrawingDetails.highlightByReactant;
  highlightColorsReactants = rxnDrawingDetails.highlightColorsReactants;
  return problems;
}

std::string molblock_helper(const RWMol &mol, const char *details_json,
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
  const RWMol *molPtr = &mol;
  std::unique_ptr<RWMol> molCopy;
  if (useMolBlockWedging || addChiralHs) {
    molCopy.reset(new RWMol(mol));
    molPtr = molCopy.get();
    if (useMolBlockWedging) {
      RDKit::Chirality::reapplyMolBlockWedging(*molCopy);
    }
    if (addChiralHs) {
      MolDraw2DUtils::prepareMolForDrawing(*molCopy, false, true, false, false,
                                           false);
    }
  }
  return MolToMolBlock(*molPtr, includeStereo, -1, kekulize, forceV3000);
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
  std::vector<int> invMatch(q_mol.getNumAtoms(), -1);
  for (const auto &pair : match) {
    invMatch[pair.first] = &pair - &match.front();
  }
  for (const auto qbond : q_mol.bonds()) {
    auto beginIdx = invMatch.at(qbond->getBeginAtomIdx());
    auto endIdx = invMatch.at(qbond->getEndAtomIdx());
    if (beginIdx == -1 || endIdx == -1) {
      continue;
    }
    unsigned int idx1 = match.at(beginIdx).second;
    unsigned int idx2 = match.at(endIdx).second;
    const auto bond = d_mol.getBondBetweenAtoms(idx1, idx2);
    if (bond != nullptr) {
      rjBonds.PushBack(bond->getIdx(), doc.GetAllocator());
    }
  }
  obj.AddMember("bonds", rjBonds, doc.GetAllocator());
}

namespace {
class SVGDrawerFromDetails : public DrawerFromDetails {
 public:
  SVGDrawerFromDetails(int w = -1, int h = -1,
                       const std::string &details = std::string()) {
    init(w, h, details);
  }

 private:
  MolDraw2DSVG &drawer() const {
    CHECK_INVARIANT(d_drawer, "d_drawer must not be null");
    return *d_drawer;
  };
  void initDrawer(const DrawingDetails &drawingDetails) {
    d_drawer.reset(new MolDraw2DSVG(
        drawingDetails.width, drawingDetails.height, drawingDetails.panelWidth,
        drawingDetails.panelHeight, drawingDetails.noFreetype));
    updateDrawerParamsFromJSON();
  }
  std::string finalizeDrawing() {
    CHECK_INVARIANT(d_drawer, "d_drawer must not be null");
    d_drawer->finishDrawing();
    return d_drawer->getDrawingText();
  }
  std::unique_ptr<MolDraw2DSVG> d_drawer;
};
}  // end anonymous namespace

std::string mol_to_svg(const ROMol &m, int w = -1, int h = -1,
                       const std::string &details = "") {
  SVGDrawerFromDetails svgDrawer(w, h, details);
  return svgDrawer.draw_mol(m);
}

std::string rxn_to_svg(const ChemicalReaction &rxn, int w = -1, int h = -1,
                       const std::string &details = "") {
  SVGDrawerFromDetails svgDrawer(w, h, details);
  return svgDrawer.draw_rxn(rxn);
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
    RDKit::Chirality::invertMolBlockWedgingInfo(mol);
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
  RDDepict::ConstrainedDepictionParams p;
  bool useCoordGen = false;
  std::string referenceSmarts;
  if (details_json && strlen(details_json)) {
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(ss, pt);
    LPT_OPT_GET(useCoordGen);
    LPT_OPT_GET(referenceSmarts);
    LPT_OPT_GET2(p, allowRGroups);
    LPT_OPT_GET2(p, acceptFailure);
    LPT_OPT_GET2(p, alignOnly);
  }
  int confId = -1;
  MatchVectType match;
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  bool oprefer = RDDepict::preferCoordGen;
  RDDepict::preferCoordGen = useCoordGen;
#endif
  std::unique_ptr<ROMol> refPattern;
  if (!referenceSmarts.empty()) {
    try {
      refPattern.reset(SmartsToMol(referenceSmarts));
    } catch (...) {
    }
  }
  try {
    match = RDDepict::generateDepictionMatching2DStructure(
        mol, templateMol, confId, refPattern.get(), p);
  } catch (...) {
  }
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = oprefer;
#endif
  if (match.empty()) {
    res = (p.acceptFailure ? "{}" : "");
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

std::string parse_inchi_options(const char *details_json) {
  std::string options;
  if (details_json && strlen(details_json)) {
    boost::property_tree::ptree pt;
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    LPT_OPT_GET(options);
  }
  return options;
}

struct LogHandle {
 public:
  LogHandle(const std::string &logName) : d_logName(logName) {
    d_logNameToLoggers = std::map<std::string, LoggerStateVector>{
        {"rdApp.debug", {LoggerState(rdDebugLog)}},
        {"rdApp.info", {LoggerState(rdInfoLog)}},
        {"rdApp.warning", {LoggerState(rdWarningLog)}},
        {"rdApp.error", {LoggerState(rdErrorLog)}},
        {"rdApp.*",
         {LoggerState(rdDebugLog), LoggerState(rdInfoLog),
          LoggerState(rdWarningLog), LoggerState(rdErrorLog)}}};
  }
  ~LogHandle() { close(); }
  static void enableLogging() {
    initLogsIfNeeded();
    boost::logging::enable_logs("rdApp.*");
  }
  static void disableLogging() {
    initLogsIfNeeded();
    boost::logging::disable_logs("rdApp.*");
  }
  void clearBuffer() {
    d_stream.str({});
    d_stream.clear();
  }
  std::string getBuffer() {
    d_stream.flush();
    return d_stream.str();
  }
  static LogHandle *setLogTee(const char *logNameCStr) {
    return setLogCommon(logNameCStr, true);
  }
  static LogHandle *setLogCapture(const char *logNameCStr) {
    return setLogCommon(logNameCStr, false);
  }

 private:
  struct LoggerState {
   public:
    LoggerState(RDLogger &logger) : d_logger(logger) {
      if (d_logger) {
        d_prevDest = d_logger->dp_dest;
        d_prevWasEnabled = d_logger->df_enabled;
      }
    }
    ~LoggerState() {
      if (!d_prevDest) {
        d_logger = nullptr;
      } else {
        if (d_logger->dp_dest) {
          d_logger->dp_dest->flush();
        }
        d_logger->dp_dest = d_prevDest;
        d_logger->df_enabled = d_prevWasEnabled;
      }
    }
    const RDLogger &logger() const { return d_logger; }
    std::ostream *stream() const { return d_logger->dp_dest; }
    void setStream(std::ostream &ostream) {
      if (!d_logger) {
        d_logger = std::make_shared<boost::logging::rdLogger>(&ostream);
      } else {
        if (d_logger->dp_dest) {
          d_logger->dp_dest->flush();
        }
        d_logger->dp_dest = &ostream;
      }
    }

   private:
    RDLogger &d_logger;
    std::ostream *d_prevDest = nullptr;
    bool d_prevWasEnabled = false;
  };
  typedef std::vector<LoggerState> LoggerStateVector;
#ifdef RDK_BUILD_THREADSAFE_SSS
  typedef std::atomic_bool LoggingFlag;
#else
  typedef bool LoggingFlag;
#endif
  bool open(bool setTee) {
    d_haveTee = setTee;
    auto loggerStates = getLoggerStates();
    if (!loggerStates) {
      return false;
    }
    clearBuffer();
    if (d_haveTee) {
      initLogsIfNeeded();
    }
    for (auto &loggerState : *loggerStates) {
      if (d_haveTee) {
        CHECK_INVARIANT(loggerState.logger(), "");
        loggerState.logger()->SetTee(d_stream);
      } else {
        loggerState.setStream(d_stream);
      }
      loggerState.logger()->df_enabled = true;
    }
    return true;
  }
  void close() {
    const auto loggerStates = getLoggerStates();
    if (!loggerStates) {
      return;
    }
    for (const auto &loggerState : *loggerStates) {
      if (d_haveTee) {
        CHECK_INVARIANT(loggerState.logger(), "");
        loggerState.logger()->ClearTee();
      }
    }
  }
  static LogHandle *setLogCommon(const char *logNameCStr, bool setTee) {
    const auto logName = std::string(logNameCStr);
    std::unique_ptr<MinimalLib::LogHandle> log_handle(
        new MinimalLib::LogHandle(logName));
    return (log_handle->open(setTee) ? log_handle.release() : nullptr);
  }
  // init logs if not yet initialized; returns true
  // if they were actually initialized, false if not
  static bool initLogsIfNeeded() {
    if (d_loggingNeedsInit) {
      RDLog::InitLogs();
      d_loggingNeedsInit = false;
      return true;
    }
    return false;
  }
  // returns nullptr if no loggers can be found
  LoggerStateVector *getLoggerStates() {
    const auto it = d_logNameToLoggers.find(d_logName);
    if (it == d_logNameToLoggers.end()) {
      return nullptr;
    }
    return &it->second;
  }
  bool d_haveTee;
  std::map<std::string, LoggerStateVector> d_logNameToLoggers;
  std::string d_logName;
  std::stringstream d_stream;
  static LoggingFlag d_loggingNeedsInit;
};

}  // namespace MinimalLib
}  // namespace RDKit
#undef LPT_OPT_GET
#undef LPT_OPT_GET2
