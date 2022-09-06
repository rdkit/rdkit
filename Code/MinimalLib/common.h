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
#include <RDGeneral/FileParseException.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/Descriptors/Property.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#ifdef RDK_BUILD_AVALON_SUPPORT
#include <External/AvalonTools/AvalonTools.h>
#endif
#include <GraphMol/Depictor/RDDepictor.h>
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
    } else if (input.find("commonchem") != std::string::npos) {
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
      LPT_OPT_GET2(ps, useLegacyStereo);
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
  } else if (input.find("commonchem") != std::string::npos) {
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

std::string process_details(const std::string &details, int &width, int &height,
                            int &offsetx, int &offsety, std::string &legend,
                            std::vector<int> &atomIds,
                            std::vector<int> &bondIds, bool &kekulize) {
  rj::Document doc;
  doc.Parse(details.c_str());
  if (!doc.IsObject()) {
    return "Invalid JSON";
  }

  const auto atomsIt = doc.FindMember("atoms");
  if (atomsIt != doc.MemberEnd()) {
    if (!atomsIt->value.IsArray()) {
      return "JSON contains 'atoms' field, but it is not an array";
    }
    for (const auto &molval : atomsIt->value.GetArray()) {
      if (!molval.IsInt()) {
        return ("Atom IDs should be integers");
      }
      atomIds.push_back(molval.GetInt());
    }
  }

  const auto bondsIt = doc.FindMember("bonds");
  if (bondsIt != doc.MemberEnd()) {
    if (!bondsIt->value.IsArray()) {
      return "JSON contains 'bonds' field, but it is not an array";
    }
    for (const auto &molval : bondsIt->value.GetArray()) {
      if (!molval.IsInt()) {
        return ("Bond IDs should be integers");
      }
      bondIds.push_back(molval.GetInt());
    }
  }

  const auto widthIt = doc.FindMember("width");
  if (widthIt != doc.MemberEnd()) {
    if (!widthIt->value.IsInt()) {
      return "JSON contains 'width' field, but it is not an int";
    }
    width = widthIt->value.GetInt();
  }

  const auto heightIt = doc.FindMember("height");
  if (heightIt != doc.MemberEnd()) {
    if (!heightIt->value.IsInt()) {
      return "JSON contains 'height' field, but it is not an int";
    }
    height = heightIt->value.GetInt();
  }

  const auto offsetxIt = doc.FindMember("offsetx");
  if (offsetxIt != doc.MemberEnd()) {
    if (!offsetxIt->value.IsInt()) {
      return "JSON contains 'offsetx' field, but it is not an int";
    }
    offsetx = offsetxIt->value.GetInt();
  }

  const auto offsetyIt = doc.FindMember("offsety");
  if (offsetyIt != doc.MemberEnd()) {
    if (!offsetyIt->value.IsInt()) {
      return "JSON contains 'offsety' field, but it is not an int";
    }
    offsety = offsetyIt->value.GetInt();
  }

  const auto legendIt = doc.FindMember("legend");
  if (legendIt != doc.MemberEnd()) {
    if (!legendIt->value.IsString()) {
      return "JSON contains 'legend' field, but it is not a string";
    }
    legend = legendIt->value.GetString();
  }

  const auto kekulizeIt = doc.FindMember("kekulize");
  if (kekulizeIt != doc.MemberEnd()) {
    if (!kekulizeIt->value.IsBool()) {
      return "JSON contains 'kekulize' field, but it is not a bool";
    }
    kekulize = kekulizeIt->value.GetBool();
  } else {
    kekulize = true;
  }

  return "";
}

std::string process_rxn_details(
    const std::string &details, int &width, int &height, int &offsetx,
    int &offsety, std::string &legend, std::vector<int> &atomIds,
    std::vector<int> &bondIds, bool &kekulize, bool &highlightByReactant,
    std::vector<DrawColour> &highlightColorsReactants) {
  auto problems = process_details(details, width, height, offsetx, offsety,
                                  legend, atomIds, bondIds, kekulize);
  if (!problems.empty()) {
    return problems;
  }
  rj::Document doc;
  doc.Parse(details.c_str());
  if (!doc.IsObject()) {
    return "Invalid JSON";
  }
  auto highlightByReactantIt = doc.FindMember("highlightByReactant");
  if (highlightByReactantIt != doc.MemberEnd()) {
    if (!highlightByReactantIt->value.IsBool()) {
      return "JSON contains 'highlightByReactant' field, but it is not a bool";
    }
    highlightByReactant = highlightByReactantIt->value.GetBool();
  } else {
    highlightByReactant = false;
  }
  auto highlightColorsReactantsIt = doc.FindMember("highlightColorsReactants");
  if (highlightColorsReactantsIt != doc.MemberEnd()) {
    if (!highlightColorsReactantsIt->value.IsArray()) {
      return "JSON contains 'highlightColorsReactants' field, but it is not an "
             "array";
    }
    for (const auto &color : highlightColorsReactantsIt->value.GetArray()) {
      if (!color.IsArray() || color.Size() < 3 || color.Size() > 4) {
        return "JSON contains 'highlightColorsReactants' array, but the "
               "elements are not R,G,B[,A] arrays";
      }
      std::vector<double> rgba(4, 1.0);
      unsigned int i = 0;
      for (const auto &component : color.GetArray()) {
        if (!component.IsDouble()) {
          return "JSON contains 'highlightColorsReactants' array, but the "
                 "R,G,B[,A] arrays contain non-float values";
        }
        CHECK_INVARIANT(i < 4, "");
        rgba[i++] = component.GetDouble();
      }
      highlightColorsReactants.emplace_back(rgba[0], rgba[1], rgba[2], rgba[3]);
    }
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

std::string mol_to_svg(const ROMol &m, int w, int h,
                       const std::string &details = "") {
  std::vector<int> atomIds;
  std::vector<int> bondIds;
  std::string legend = "";
  int offsetx = 0;
  int offsety = 0;
  bool kekulize = true;
  if (!details.empty()) {
    auto problems = process_details(details, w, h, offsetx, offsety, legend,
                                    atomIds, bondIds, kekulize);
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
                                         nullptr, nullptr, nullptr, -1,
                                         kekulize);
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

}  // namespace MinimalLib
}  // namespace RDKit
#undef LPT_OPT_GET
#undef LPT_OPT_GET2
