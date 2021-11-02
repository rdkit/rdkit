//
//
//  Copyright (C) 2019-2021 Greg Landrum and other RDKit contributors
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
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#include <GraphMol/Descriptors/Property.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/MolInterchange/MolInterchange.h>
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
std::string JSMol::get_smarts() const {
  if (!d_mol) return "";
  return MolToSmarts(*d_mol);
}
std::string JSMol::get_cxsmarts() const {
  if (!d_mol) return "";
  return MolToCXSmarts(*d_mol);
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
std::string JSMol::get_json() const {
  if (!d_mol) return "";
  return MolInterchange::MolToJSONData(*d_mol);
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

std::string JSMol::get_morgan_fp_as_binary_text(unsigned int radius,
                                                unsigned int fplen) const {
  if (!d_mol) return "";
  auto fp = MorganFingerprints::getFingerprintAsBitVect(*d_mol, radius, fplen);
  std::string res = BitVectToBinaryText(*fp);
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
                                           bool useCoordGen, bool allowRGroups,
                                           bool acceptFailure) {
  std::string res;
  if (!d_mol || !templateMol.d_mol || !templateMol.d_mol->getNumConformers())
    return res;

#ifdef RDK_BUILD_COORDGEN_SUPPORT
  bool oprefer = RDDepict::preferCoordGen;
  RDDepict::preferCoordGen = useCoordGen;
#endif
  RDKit::ROMol *refPattern = nullptr;
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
  } else {
    res = "{}";
  }
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = oprefer;
#endif
  return res;
};

JSSubstructLibrary::JSSubstructLibrary(unsigned int num_bits)
    : d_sslib(new SubstructLibrary(
          boost::shared_ptr<CachedTrustedSmilesMolHolder>(
              new CachedTrustedSmilesMolHolder()),
          boost::shared_ptr<PatternHolder>(new PatternHolder()))),
      d_num_bits(num_bits) {
  d_molHolder = dynamic_cast<CachedTrustedSmilesMolHolder *>(
      d_sslib->getMolHolder().get());
  d_fpHolder = dynamic_cast<PatternHolder *>(d_sslib->getFpHolder().get());
}

int JSSubstructLibrary::add_trusted_smiles(const std::string &smi) {
  std::unique_ptr<RWMol> mol(SmilesToMol(smi, 0, false));
  if (!mol) {
    return -1;
  }
  mol->updatePropertyCache();
  ExplicitBitVect *bv = PatternFingerprintMol(*mol, d_num_bits);
  if (!bv) {
    return -1;
  }
  d_fpHolder->addFingerprint(bv);
  auto ret = d_molHolder->addSmiles(smi);
  return ret;
}

inline int JSSubstructLibrary::add_mol_helper(const ROMol &mol) {
  std::string smi = MolToSmiles(mol);
  return add_trusted_smiles(smi);
}

int JSSubstructLibrary::add_mol(const JSMol &m) {
  return add_mol_helper(*m.d_mol);
}

int JSSubstructLibrary::add_smiles(const std::string &smi) {
  std::unique_ptr<RWMol> mol(SmilesToMol(smi));
  if (!mol) {
    return -1;
  }
  return add_mol_helper(*mol);
}

JSMol *JSSubstructLibrary::get_mol(unsigned int i) {
  return new JSMol(new RWMol(*d_sslib->getMol(i)));
}

std::string JSSubstructLibrary::get_matches(const JSMol &q, bool useChirality,
                                            int numThreads,
                                            int maxResults) const {
  if (!d_sslib->size()) {
    return "[]";
  }
  std::vector<unsigned int> indices = d_sslib->getMatches(
      *q.d_mol, true, useChirality, false, numThreads, maxResults);
  rj::Document doc;
  doc.SetArray();
  auto &alloc = doc.GetAllocator();
  for (const auto &i : indices) {
    doc.PushBack(i, alloc);
  }
  rj::StringBuffer buffer;
  rj::Writer<rj::StringBuffer> writer(buffer);
  doc.Accept(writer);
  std::string res = buffer.GetString();
  return res;
}

unsigned int JSSubstructLibrary::count_matches(const JSMol &q,
                                               bool useChirality,
                                               int numThreads) const {
  return d_sslib->countMatches(*q.d_mol, true, useChirality, false, 1);
}

std::string get_inchikey_for_inchi(const std::string &input) {
  return InchiToInchiKey(input);
}

JSMol *get_mol(const std::string &input, const std::string &details_json) {
  RWMol *mol = MinimalLib::mol_from_input(input, details_json);
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
