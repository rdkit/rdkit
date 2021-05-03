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
#include "common.h"

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
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Abbreviations/Abbreviations.h>
#include <DataStructs/BitOps.h>

#include <INCHI-API/inchi.h>

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

namespace rj = rapidjson;

using namespace RDKit;

std::string JSMol::get_smiles() const {
  if (!d_mol) return "";
  return MolToSmiles(*d_mol);
}
std::string JSMol::get_cxsmiles() const {
  if (!d_mol) return "";
  return MolToCXSmiles(*d_mol);
}
std::string JSMol::get_svg(unsigned int w, unsigned int h) const {
  if (!d_mol) return "";
  return MinimalLib::mol_to_svg(*d_mol, w, h);
}
std::string JSMol::get_svg_with_highlights(const std::string &details) const {
  if (!d_mol) return "";

  unsigned int w = MinimalLib::d_defaultWidth;
  unsigned int h = MinimalLib::d_defaultHeight;
  return MinimalLib::mol_to_svg(*d_mol, w, h, details);
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
std::string JSMol::get_v3Kmolblock() const {
  if (!d_mol) return "";
  return MolToV3KMolBlock(*d_mol);
}

std::string JSMol::get_substruct_match(const JSMol &q) const {
  std::string res = "{}";
  if (!d_mol || !q.d_mol) return res;

  MatchVectType match;
  if (SubstructMatch(*d_mol, *(q.d_mol), match)) {
    rj::Document doc;
    doc.SetObject();
    MinimalLib::get_sss_json(*d_mol, *(q.d_mol), match, doc, doc);
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

    for (const auto &match : matches) {
      rj::Value rjMatch(rj::kObjectType);
      MinimalLib::get_sss_json(*d_mol, *(q.d_mol), match, rjMatch, doc);
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
  return MinimalLib::get_descriptors(*d_mol);
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
  CIPLabeler::assignCIPLabels(*d_mol);

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
    std::string cip;
    if (bond->getPropIfPresent(common_properties::_CIPCode, cip)) {
      cip = "(" + cip + ")";
      rj::Value entry(rj::kArrayType);
      entry.PushBack(bond->getBeginAtomIdx(), doc.GetAllocator());
      entry.PushBack(bond->getEndAtomIdx(), doc.GetAllocator());
      rj::Value v;
      v.SetString(cip.c_str(), cip.size(), doc.GetAllocator());
      entry.PushBack(v, doc.GetAllocator());
      rjBonds.PushBack(entry, doc.GetAllocator());
    }
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

std::string JSMol::condense_abbreviations(double maxCoverage, bool useLinkers) {
  if (!d_mol) return "";
  if (!useLinkers) {
    Abbreviations::condenseMolAbbreviations(
        *d_mol, Abbreviations::Utils::getDefaultAbbreviations(), maxCoverage);
  } else {
    Abbreviations::condenseMolAbbreviations(
        *d_mol, Abbreviations::Utils::getDefaultLinkers(), maxCoverage);
  }
  return "";
}

std::string JSMol::condense_abbreviations_from_defs(
    const std::string &definitions, double maxCoverage, bool areLinkers) {
  static std::string lastDefs = "";
  static std::vector<Abbreviations::AbbreviationDefinition> abbrevs;
  if (definitions != lastDefs) {
    // yes, we are making the assumption that the "areLinkers" argument remains
    // the same if the definitions are the same
    bool removeExtraDummies = areLinkers;
    bool allowConnectionToDummies = areLinkers;
    lastDefs = definitions;
    try {
      abbrevs = Abbreviations::Utils::parseAbbreviations(
          definitions, removeExtraDummies, allowConnectionToDummies);
    } catch (...) {
      return "cannot parse abbreviations";
    }
  }
  Abbreviations::condenseMolAbbreviations(*d_mol, abbrevs, maxCoverage);
}

std::string JSMol::generate_aligned_coords(const JSMol &templateMol,
                                           bool useCoordGen,
                                           bool allowRGroups) {
  std::string res;
  if (!d_mol || !templateMol.d_mol || !templateMol.d_mol->getNumConformers())
    return res;

#ifdef RDK_BUILD_COORDGEN_SUPPORT
  bool oprefer = RDDepict::preferCoordGen;
  RDDepict::preferCoordGen = useCoordGen;
#endif
  RDKit::ROMol *refPattern = nullptr;
  bool acceptFailure = true;
  int confId = -1;
  RDKit::MatchVectType match = RDDepict::generateDepictionMatching2DStructure(
      *d_mol, *(templateMol.d_mol), confId, refPattern, acceptFailure, false,
      allowRGroups);
  if (!match.empty()) {
    rj::Document doc;
    doc.SetObject();
    MinimalLib::get_sss_json(*d_mol, *templateMol.d_mol, match, doc, doc);
    rj::StringBuffer buffer;
    rj::Writer<rj::StringBuffer> writer(buffer);
    doc.Accept(writer);
    res = buffer.GetString();
  }
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = oprefer;
#endif
  return res;
};

std::string get_inchikey_for_inchi(const std::string &input) {
  return InchiToInchiKey(input);
}

JSMol *get_mol(const std::string &input) {
  RWMol *mol = MinimalLib::mol_from_input(input);
  return new JSMol(mol);
}

JSMol *get_qmol(const std::string &input) {
  RWMol *mol = MinimalLib::qmol_from_input(input);
  return new JSMol(mol);
}

std::string version() { return std::string(rdkitVersion); }

void prefer_coordgen(bool useCoordGen) {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = useCoordGen;
#endif
}
