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
#include <GraphMol/Chirality.h>
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
#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Abbreviations/Abbreviations.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <Geometry/Transform3D.h>
#include <DataStructs/BitOps.h>

#include <INCHI-API/inchi.h>

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

namespace rj = rapidjson;

using namespace RDKit;

namespace {
std::string mappingToJsonArray(const ROMol &mol) {
  std::vector<unsigned int> atomMapping;
  std::vector<unsigned int> bondMapping;
  mol.getPropIfPresent(Abbreviations::common_properties::origAtomMapping,
                       atomMapping);
  mol.getPropIfPresent(Abbreviations::common_properties::origBondMapping,
                       bondMapping);
  rj::Document doc;
  doc.SetObject();
  auto &alloc = doc.GetAllocator();
  rj::Value rjAtoms(rj::kArrayType);
  for (auto i : atomMapping) {
    rjAtoms.PushBack(i, alloc);
  }
  doc.AddMember("atoms", rjAtoms, alloc);

  rj::Value rjBonds(rj::kArrayType);
  for (auto i : bondMapping) {
    rjBonds.PushBack(i, alloc);
  }
  doc.AddMember("bonds", rjBonds, alloc);
  rj::StringBuffer buffer;
  rj::Writer<rj::StringBuffer> writer(buffer);
  doc.Accept(writer);
  std::string res = buffer.GetString();
  return res;
}
}  // end of anonymous namespace

std::string JSMol::get_smiles() const {
  if (!d_mol) {
    return "";
  }
  return MolToSmiles(*d_mol);
}
std::string JSMol::get_cxsmiles() const {
  if (!d_mol) {
    return "";
  }
  return MolToCXSmiles(*d_mol);
}
std::string JSMol::get_smarts() const {
  if (!d_mol) {
    return "";
  }
  return MolToSmarts(*d_mol);
}
std::string JSMol::get_cxsmarts() const {
  if (!d_mol) {
    return "";
  }
  return MolToCXSmarts(*d_mol);
}
std::string JSMol::get_svg(int w, int h) const {
  if (!d_mol) {
    return "";
  }
  return MinimalLib::mol_to_svg(*d_mol, w, h);
}
std::string JSMol::get_svg_with_highlights(const std::string &details) const {
  if (!d_mol) {
    return "";
  }

  int w = d_defaultWidth;
  int h = d_defaultHeight;
  return MinimalLib::mol_to_svg(*d_mol, w, h, details);
}

std::string JSMol::get_inchi() const {
  if (!d_mol) {
    return "";
  }
  ExtraInchiReturnValues rv;
  return MolToInchi(*d_mol, rv);
}
std::string JSMol::get_molblock(const std::string &details) const {
  if (!d_mol) {
    return "";
  }
  return MinimalLib::molblock_helper(*d_mol, details.c_str(), false);
}
std::string JSMol::get_v3Kmolblock(const std::string &details) const {
  if (!d_mol) {
    return "";
  }
  return MinimalLib::molblock_helper(*d_mol, details.c_str(), true);
}
std::string JSMol::get_json() const {
  if (!d_mol) {
    return "";
  }
  return MolInterchange::MolToJSONData(*d_mol);
}

std::string JSMol::get_pickle() const {
  if (!d_mol) {
    return "";
  }
  std::string pickle;
  MolPickler::pickleMol(*d_mol, pickle,
                        PicklerOps::AllProps ^ PicklerOps::ComputedProps);
  return pickle;
}

std::string JSMol::get_substruct_match(const JSMol &q) const {
  std::string res = "{}";
  if (!d_mol || !q.d_mol) {
    return res;
  }

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
  if (!d_mol || !q.d_mol) {
    return res;
  }

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
  if (!d_mol) {
    return "{}";
  }
  return MinimalLib::get_descriptors(*d_mol);
}

std::string JSMol::get_morgan_fp(const std::string &details) const {
  if (!d_mol) {
    return "";
  }
  auto fp = MinimalLib::morgan_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToText(*fp);
  return res;
}

std::string JSMol::get_morgan_fp_as_binary_text(
    const std::string &details) const {
  if (!d_mol) {
    return "";
  }
  auto fp = MinimalLib::morgan_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToBinaryText(*fp);
  return res;
}

std::string JSMol::get_pattern_fp(const std::string &details) const {
  if (!d_mol) {
    return "";
  }
  auto fp = MinimalLib::pattern_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToText(*fp);
  return res;
}

std::string JSMol::get_pattern_fp_as_binary_text(
    const std::string &details) const {
  if (!d_mol) {
    return "";
  }
  auto fp = MinimalLib::pattern_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToBinaryText(*fp);
  return res;
}

std::string JSMol::get_topological_torsion_fp(
    const std::string &details) const {
  if (!d_mol) {
    return "";
  }
  auto fp =
      MinimalLib::topological_torsion_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToText(*fp);
  return res;
}

std::string JSMol::get_topological_torsion_fp_as_binary_text(
    const std::string &details) const {
  if (!d_mol) {
    return "";
  }
  auto fp =
      MinimalLib::topological_torsion_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToBinaryText(*fp);
  return res;
}

std::string JSMol::get_rdkit_fp(const std::string &details) const {
  if (!d_mol) {
    return "";
  }
  auto fp = MinimalLib::rdkit_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToText(*fp);
  return res;
}

std::string JSMol::get_rdkit_fp_as_binary_text(
    const std::string &details) const {
  if (!d_mol) {
    return "";
  }
  auto fp = MinimalLib::rdkit_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToBinaryText(*fp);
  return res;
}

std::string JSMol::get_atom_pair_fp(const std::string &details) const {
  if (!d_mol) {
    return "";
  }
  auto fp = MinimalLib::atom_pair_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToText(*fp);
  return res;
}

std::string JSMol::get_atom_pair_fp_as_binary_text(
    const std::string &details) const {
  if (!d_mol) {
    return "";
  }
  auto fp = MinimalLib::atom_pair_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToBinaryText(*fp);
  return res;
}

std::string JSMol::get_maccs_fp() const {
  if (!d_mol) {
    return "";
  }
  auto fp = MinimalLib::maccs_fp_as_bitvect(*d_mol);
  std::string res = BitVectToText(*fp);
  return res;
}

std::string JSMol::get_maccs_fp_as_binary_text() const {
  if (!d_mol) {
    return "";
  }
  auto fp = MinimalLib::maccs_fp_as_bitvect(*d_mol);
  std::string res = BitVectToBinaryText(*fp);
  return res;
}

#ifdef RDK_BUILD_AVALON_SUPPORT
std::string get_avalon_fp(const std::string &details) const {
  if (!d_mol) {
    return "";
  }
  auto fp = MinimalLib::avalon_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToText(*fp);
  return res;
}

std::string get_avalon_fp_as_binary_text(const std::string &details) const {
  if (!d_mol) {
    return "";
  }
  auto fp = MinimalLib::avalon_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToBinaryText(*fp);
  return res;
}
#endif

std::string JSMol::get_stereo_tags() const {
  if (!d_mol) {
    return "{}";
  }
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
  if (!d_mol) {
    return "";
  }

  RWMol molCopy(*d_mol);
  MolOps::setAromaticity(molCopy);

  bool includeStereo = true;
  int confId = -1;
  bool kekulize = false;
  return MolToMolBlock(molCopy, includeStereo, confId, kekulize);
}

std::string JSMol::get_kekule_form() const {
  if (!d_mol) {
    return "";
  }

  RWMol molCopy(*d_mol);
  MolOps::Kekulize(molCopy);

  bool includeStereo = true;
  int confId = -1;
  bool kekulize = true;
  return MolToMolBlock(molCopy, includeStereo, confId, kekulize);
}

bool JSMol::set_new_coords(bool useCoordGen) {
  if (!d_mol) {
    return false;
  }

#ifdef RDK_BUILD_COORDGEN_SUPPORT
  bool oprefer = RDDepict::preferCoordGen;
  RDDepict::preferCoordGen = useCoordGen;
#endif
  RDDepict::compute2DCoords(*d_mol);
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = oprefer;
#endif

  return true;
}

std::string JSMol::get_new_coords(bool useCoordGen) const {
  if (!d_mol) {
    return "";
  }

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

bool JSMol::has_prop(const std::string &key) const {
  if (!d_mol) return false;
  return d_mol->hasProp(key);
}

std::vector<std::string> JSMol::get_prop_list(bool includePrivate,
                                              bool includeComputed) const {
  if (!d_mol) return std::vector<std::string>();
  return d_mol->getPropList(includePrivate, includeComputed);
}

bool JSMol::set_prop(const std::string &key, const std::string &val,
                     bool computed) {
  if (!d_mol) return false;
  d_mol->setProp(key, val, computed);
  return true;
}

std::string JSMol::get_prop(const std::string &key) const {
  if (!d_mol || !d_mol->hasProp(key)) return "";
  std::string val;
  d_mol->getProp(key, val);
  return val;
}

bool JSMol::clear_prop(const std::string &key) {
  if (!d_mol) return false;
  bool res = d_mol->hasProp(key);
  if (res) {
    d_mol->clearProp(key);
  }
  return res;
}

std::string JSMol::remove_hs() const {
  if (!d_mol) {
    return "";
  }

  RWMol molCopy(*d_mol);
  MolOps::removeAllHs(molCopy);

  bool includeStereo = true;
  int confId = -1;
  bool kekulize = true;
  return MolToMolBlock(molCopy, includeStereo, confId, kekulize);
}

bool JSMol::remove_hs_in_place() {
  if (!d_mol) {
    return false;
  }

  MolOps::removeAllHs(*d_mol);
  MolOps::assignStereochemistry(*d_mol, true, true);
  return true;
}

std::string JSMol::add_hs() const {
  if (!d_mol) {
    return "";
  }

  RWMol molCopy(*d_mol);
  MolOps::addHs(molCopy);

  bool includeStereo = true;
  int confId = -1;
  bool kekulize = true;
  return MolToMolBlock(molCopy, includeStereo, confId, kekulize);
}

bool JSMol::add_hs_in_place() {
  if (!d_mol) {
    return false;
  }

  bool addCoords = (d_mol->getNumConformers() > 0);
  MolOps::addHs(*d_mol, false, addCoords);
  MolOps::assignStereochemistry(*d_mol, true, true);
  return true;
}

std::string JSMol::condense_abbreviations(double maxCoverage, bool useLinkers) {
  if (!d_mol) {
    return "";
  }
  if (!useLinkers) {
    Abbreviations::condenseMolAbbreviations(
        *d_mol, Abbreviations::Utils::getDefaultAbbreviations(), maxCoverage);
  } else {
    Abbreviations::condenseMolAbbreviations(
        *d_mol, Abbreviations::Utils::getDefaultLinkers(), maxCoverage);
  }
  return mappingToJsonArray(*d_mol);
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
  return mappingToJsonArray(*d_mol);
}

std::string JSMol::generate_aligned_coords(const JSMol &templateMol,
                                           const std::string &details) {
  if (!d_mol || !templateMol.d_mol || !templateMol.d_mol->getNumConformers()) {
    return "";
  }
  return MinimalLib::generate_aligned_coords(*d_mol, *templateMol.d_mol,
                                             details.c_str());
}

int JSMol::has_coords() const {
  if (!d_mol || !d_mol->getNumConformers()) {
    return 0;
  }
  return (d_mol->getConformer().is3D() ? 3 : 2);
}

double JSMol::normalize_depiction(int canonicalize, double scaleFactor) {
  if (!d_mol || !d_mol->getNumConformers()) {
    return -1.;
  }
  return RDDepict::normalizeDepiction(*d_mol, -1, canonicalize, scaleFactor);
}

void JSMol::straighten_depiction(bool minimizeRotation) {
  if (!d_mol || !d_mol->getNumConformers()) {
    return;
  }
  RDDepict::straightenDepiction(*d_mol, -1, minimizeRotation);
}

std::pair<JSMolIterator *, std::string> JSMol::get_frags(
    const std::string &details_json) {
  if (!d_mol) {
    return std::make_pair(nullptr, "");
  }
  std::vector<int> frags;
  std::vector<std::vector<int>> fragsMolAtomMapping;
  bool sanitizeFrags = true;
  bool copyConformers = true;
  MinimalLib::get_mol_frags_details(details_json, sanitizeFrags,
                                    copyConformers);
  auto molFrags = MolOps::getMolFrags(*d_mol, sanitizeFrags, &frags,
                                      &fragsMolAtomMapping, copyConformers);
  return std::make_pair(
      new JSMolIterator(molFrags),
      MinimalLib::get_mol_frags_mappings(frags, fragsMolAtomMapping));
}

std::string JSReaction::get_svg(int w, int h) const {
  if (!d_rxn) {
    return "";
  }
  return MinimalLib::rxn_to_svg(*d_rxn, w, h);
}
std::string JSReaction::get_svg_with_highlights(
    const std::string &details) const {
  if (!d_rxn) {
    return "";
  }

  int w = d_defaultWidth;
  int h = d_defaultHeight;
  return MinimalLib::rxn_to_svg(*d_rxn, w, h, details);
}

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
  auto fp = PatternFingerprintMol(*mol, d_num_bits);
  if (!fp) {
    return -1;
  }
  d_fpHolder->addFingerprint(fp);
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

JSMol *get_mol_copy(const JSMol &other) {
  auto mol = new RWMol(*other.d_mol);
  return new JSMol(mol);
}

JSMol *get_mol(const std::string &input, const std::string &details_json) {
  auto mol = MinimalLib::mol_from_input(input, details_json);
  return new JSMol(mol);
}

JSMol *get_mol_from_pickle(const std::string &pkl) {
  RWMol *mol = nullptr;
  if (!pkl.empty()) {
    mol = new RWMol();
    try {
      MolPickler::molFromPickle(pkl, mol);
    } catch (...) {
      delete mol;
      mol = nullptr;
    }
  }
  return new JSMol(mol);
}

JSMol *get_qmol(const std::string &input) {
  auto mol = MinimalLib::qmol_from_input(input);
  return new JSMol(mol);
}

JSReaction *get_rxn(const std::string &input, const std::string &details_json) {
  auto rxn = MinimalLib::rxn_from_input(input, details_json);
  return new JSReaction(rxn);
}

std::string version() { return std::string(rdkitVersion); }

void prefer_coordgen(bool useCoordGen) {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = useCoordGen;
#endif
}

bool use_legacy_stereo_perception(bool value) {
  bool was = Chirality::getUseLegacyStereoPerception();
  Chirality::setUseLegacyStereoPerception(value);
  return was;
}

bool allow_non_tetrahedral_chirality(bool value) {
  bool was = Chirality::getAllowNontetrahedralChirality();
  Chirality::setAllowNontetrahedralChirality(value);
  return was;
}
