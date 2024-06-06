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
#ifdef RDK_BUILD_MINIMAL_LIB_SUBSTRUCTLIBRARY
#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#endif
#ifdef RDK_BUILD_MINIMAL_LIB_MCS
#include <GraphMol/FMCS/FMCS.h>
#endif
#include <GraphMol/Descriptors/Property.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Abbreviations/Abbreviations.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <Geometry/Transform3D.h>
#include <DataStructs/BitOps.h>
#include <DataStructs/ExplicitBitVect.h>

#include <INCHI-API/inchi.h>

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

namespace rj = rapidjson;

using namespace RDKit;

namespace {
static const char *NO_SUPPORT_FOR_PATTERN_FPS =
    "This SubstructLibrary was built without support for pattern fps";

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
  assert(d_mol);
  return MolToSmiles(*d_mol);
}
std::string JSMol::get_smiles(const std::string &details) const {
  assert(d_mol);
  SmilesWriteParams params;
  updateSmilesWriteParamsFromJSON(params, details);
  return MolToSmiles(*d_mol, params);
}
std::string JSMol::get_cxsmiles() const {
  assert(d_mol);
  return MolToCXSmiles(*d_mol);
}
std::string JSMol::get_cxsmiles(const std::string &details) const {
  assert(d_mol);
  SmilesWriteParams params;
  updateSmilesWriteParamsFromJSON(params, details);
  SmilesWrite::CXSmilesFields cxSmilesFields =
      SmilesWrite::CXSmilesFields::CX_ALL;
  RestoreBondDirOption restoreBondDirs = RestoreBondDirOptionClear;
  updateCXSmilesFieldsFromJSON(cxSmilesFields, restoreBondDirs, details);
  return MolToCXSmiles(*d_mol, params, cxSmilesFields, restoreBondDirs);
}
std::string JSMol::get_smarts() const {
  assert(d_mol);
  return MolToSmarts(*d_mol);
}
std::string JSMol::get_smarts(const std::string &details) const {
  assert(d_mol);
  SmilesWriteParams params;
  updateSmilesWriteParamsFromJSON(params, details);
  return MolToSmarts(*d_mol, params.doIsomericSmiles, params.rootedAtAtom);
}
std::string JSMol::get_cxsmarts() const {
  assert(d_mol);
  return MolToCXSmarts(*d_mol);
}
std::string JSMol::get_cxsmarts(const std::string &details) const {
  assert(d_mol);
  SmilesWriteParams params;
  updateSmilesWriteParamsFromJSON(params, details);
  return MolToCXSmarts(*d_mol, params.doIsomericSmiles);
}
std::string JSMol::get_svg(int w, int h) const {
  assert(d_mol);
  return MinimalLib::mol_to_svg(*d_mol, w, h);
}
std::string JSMol::get_svg_with_highlights(const std::string &details) const {
  assert(d_mol);

  int w = d_defaultWidth;
  int h = d_defaultHeight;
  return MinimalLib::mol_to_svg(*d_mol, w, h, details);
}

std::string JSMol::get_inchi(const std::string &options) const {
  assert(d_mol);
  ExtraInchiReturnValues rv;
  return MolToInchi(*d_mol, rv, !options.empty() ? options.c_str() : nullptr);
}
std::string JSMol::get_molblock(const std::string &details) const {
  assert(d_mol);
  return MinimalLib::molblock_helper(*d_mol, details.c_str(), false);
}
std::string JSMol::get_v3Kmolblock(const std::string &details) const {
  assert(d_mol);
  return MinimalLib::molblock_helper(*d_mol, details.c_str(), true);
}
std::string JSMol::get_json() const {
  assert(d_mol);
  return MolInterchange::MolToJSONData(*d_mol);
}

std::string JSMol::get_pickle() const {
  assert(d_mol);
  std::string pickle;
  MolPickler::pickleMol(*d_mol, pickle,
                        PicklerOps::AllProps ^ PicklerOps::ComputedProps);
  return pickle;
}

std::string JSMol::get_substruct_match(const JSMol &q) const {
  assert(d_mol && q.d_mol);
  std::string res = "{}";

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
  assert(d_mol && q.d_mol);

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
  assert(d_mol);
  return MinimalLib::get_descriptors(*d_mol);
}

std::string JSMol::get_morgan_fp(const std::string &details) const {
  assert(d_mol);
  auto fp = MinimalLib::morgan_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToText(*fp);
  return res;
}

std::string JSMol::get_morgan_fp_as_binary_text(
    const std::string &details) const {
  assert(d_mol);
  auto fp = MinimalLib::morgan_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToBinaryText(*fp);
  return res;
}

std::string JSMol::get_pattern_fp(const std::string &details) const {
  assert(d_mol);
  auto fp = MinimalLib::pattern_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToText(*fp);
  return res;
}

std::string JSMol::get_pattern_fp_as_binary_text(
    const std::string &details) const {
  assert(d_mol);
  auto fp = MinimalLib::pattern_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToBinaryText(*fp);
  return res;
}

std::string JSMol::get_topological_torsion_fp(
    const std::string &details) const {
  assert(d_mol);
  auto fp =
      MinimalLib::topological_torsion_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToText(*fp);
  return res;
}

std::string JSMol::get_topological_torsion_fp_as_binary_text(
    const std::string &details) const {
  assert(d_mol);
  auto fp =
      MinimalLib::topological_torsion_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToBinaryText(*fp);
  return res;
}

std::string JSMol::get_rdkit_fp(const std::string &details) const {
  assert(d_mol);
  auto fp = MinimalLib::rdkit_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToText(*fp);
  return res;
}

std::string JSMol::get_rdkit_fp_as_binary_text(
    const std::string &details) const {
  assert(d_mol);
  auto fp = MinimalLib::rdkit_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToBinaryText(*fp);
  return res;
}

std::string JSMol::get_atom_pair_fp(const std::string &details) const {
  assert(d_mol);
  auto fp = MinimalLib::atom_pair_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToText(*fp);
  return res;
}

std::string JSMol::get_atom_pair_fp_as_binary_text(
    const std::string &details) const {
  assert(d_mol);
  auto fp = MinimalLib::atom_pair_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToBinaryText(*fp);
  return res;
}

std::string JSMol::get_maccs_fp() const {
  assert(d_mol);
  auto fp = MinimalLib::maccs_fp_as_bitvect(*d_mol);
  std::string res = BitVectToText(*fp);
  return res;
}

std::string JSMol::get_maccs_fp_as_binary_text() const {
  assert(d_mol);
  auto fp = MinimalLib::maccs_fp_as_bitvect(*d_mol);
  std::string res = BitVectToBinaryText(*fp);
  return res;
}

#ifdef RDK_BUILD_AVALON_SUPPORT
std::string get_avalon_fp(const std::string &details) const {
  assert(d_mol);
  auto fp = MinimalLib::avalon_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToText(*fp);
  return res;
}

std::string get_avalon_fp_as_binary_text(const std::string &details) const {
  assert(d_mol);
  auto fp = MinimalLib::avalon_fp_as_bitvect(*d_mol, details.c_str());
  std::string res = BitVectToBinaryText(*fp);
  return res;
}
#endif

std::string JSMol::get_stereo_tags() const {
  assert(d_mol);
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

void JSMol::convert_to_aromatic_form() {
  assert(d_mol);

  d_mol->updatePropertyCache();
  MolOps::setAromaticity(*d_mol);
}

std::string JSMol::get_aromatic_form() const {
  assert(d_mol);

  RWMol molCopy(*d_mol);
  molCopy.updatePropertyCache();
  MolOps::setAromaticity(molCopy);

  bool includeStereo = true;
  int confId = -1;
  bool kekulize = false;
  return MolToMolBlock(molCopy, includeStereo, confId, kekulize);
}

void JSMol::convert_to_kekule_form() {
  assert(d_mol);

  MolOps::Kekulize(*d_mol);
}

std::string JSMol::get_kekule_form() const {
  assert(d_mol);

  RWMol molCopy(*d_mol);
  MolOps::Kekulize(molCopy);

  bool includeStereo = true;
  int confId = -1;
  bool kekulize = true;
  return MolToMolBlock(molCopy, includeStereo, confId, kekulize);
}

bool JSMol::set_new_coords(bool useCoordGen) {
  assert(d_mol);

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
  assert(d_mol);

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
  assert(d_mol);
  return d_mol->hasProp(key);
}

std::vector<std::string> JSMol::get_prop_list(bool includePrivate,
                                              bool includeComputed) const {
  assert(d_mol);
  return d_mol->getPropList(includePrivate, includeComputed);
}

bool JSMol::set_prop(const std::string &key, const std::string &val,
                     bool computed) {
  assert(d_mol);
  d_mol->setProp(key, val, computed);
  return true;
}

std::string JSMol::get_prop(const std::string &key) const {
  assert(d_mol);
  if (!d_mol->hasProp(key)) {
    return "";
  }
  std::string val;
  d_mol->getProp(key, val);
  return val;
}

bool JSMol::clear_prop(const std::string &key) {
  assert(d_mol);
  bool res = d_mol->hasProp(key);
  if (res) {
    d_mol->clearProp(key);
  }
  return res;
}

std::string JSMol::remove_hs() const {
  assert(d_mol);

  RWMol molCopy(*d_mol);
  MolOps::removeAllHs(molCopy);

  bool includeStereo = true;
  int confId = -1;
  bool kekulize = true;
  return MolToMolBlock(molCopy, includeStereo, confId, kekulize);
}

bool JSMol::remove_hs_in_place() {
  assert(d_mol);

  MolOps::removeAllHs(*d_mol);
  MolOps::assignStereochemistry(*d_mol, true, true);
  return true;
}

std::string JSMol::add_hs() const {
  assert(d_mol);

  RWMol molCopy(*d_mol);
  MolOps::addHs(molCopy);

  bool includeStereo = true;
  int confId = -1;
  bool kekulize = true;
  return MolToMolBlock(molCopy, includeStereo, confId, kekulize);
}

bool JSMol::add_hs_in_place() {
  assert(d_mol);

  bool addCoords = (d_mol->getNumConformers() > 0);
  MolOps::addHs(*d_mol, false, addCoords);
  MolOps::assignStereochemistry(*d_mol, true, true);
  return true;
}

std::string JSMol::condense_abbreviations(double maxCoverage, bool useLinkers) {
  assert(d_mol);
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
  assert(d_mol && templateMol.d_mol);
  if (!templateMol.d_mol->getNumConformers()) {
    return "";
  }
  return MinimalLib::generate_aligned_coords(*d_mol, *templateMol.d_mol,
                                             details.c_str());
}

int JSMol::has_coords() const {
  assert(d_mol);
  if (!d_mol->getNumConformers()) {
    return 0;
  }
  return (d_mol->getConformer().is3D() ? 3 : 2);
}

double JSMol::normalize_depiction(int canonicalize, double scaleFactor) {
  assert(d_mol);
  if (!d_mol->getNumConformers()) {
    return -1.;
  }
  return RDDepict::normalizeDepiction(*d_mol, -1, canonicalize, scaleFactor);
}

void JSMol::straighten_depiction(bool minimizeRotation) {
  assert(d_mol);
  if (!d_mol->getNumConformers()) {
    return;
  }
  RDDepict::straightenDepiction(*d_mol, -1, minimizeRotation);
}

bool JSMol::is_valid() const { return true; }

std::pair<JSMolList *, std::string> JSMol::get_frags(
    const std::string &details_json) {
  assert(d_mol);
  std::vector<int> frags;
  std::vector<std::vector<int>> fragsMolAtomMapping;
  bool sanitizeFrags = true;
  bool copyConformers = true;
  MinimalLib::get_mol_frags_details(details_json, sanitizeFrags,
                                    copyConformers);
  auto molFrags = MolOps::getMolFrags(*d_mol, sanitizeFrags, &frags,
                                      &fragsMolAtomMapping, copyConformers);
  return std::make_pair(
      new JSMolList(molFrags),
      MinimalLib::get_mol_frags_mappings(frags, fragsMolAtomMapping));
}

unsigned int JSMol::get_num_atoms(bool heavyOnly) const {
  assert(d_mol);
  return heavyOnly ? d_mol->getNumHeavyAtoms() : d_mol->getNumAtoms();
}

unsigned int JSMol::get_num_bonds() const {
  assert(d_mol);
  return d_mol->getNumBonds();
}

#ifdef RDK_BUILD_MINIMAL_LIB_MMPA
namespace {
bool mmpaFragmentMol(const ROMol &mol, std::vector<RDKit::ROMOL_SPTR> &cores,
                     std::vector<RDKit::ROMOL_SPTR> &sidechains,
                     unsigned int minCuts, unsigned int maxCuts,
                     unsigned int maxCutBonds) {
  std::vector<std::pair<RDKit::ROMOL_SPTR, RDKit::ROMOL_SPTR>> mmpaFrags;
  if (!RDKit::MMPA::fragmentMol(mol, mmpaFrags, minCuts, maxCuts,
                                maxCutBonds)) {
    return false;
  }
  auto numEntries = mmpaFrags.size();
  cores.clear();
  cores.reserve(numEntries);
  sidechains.clear();
  sidechains.reserve(numEntries);
  for (const auto &mmpaFrag : mmpaFrags) {
    cores.push_back(mmpaFrag.first);
    sidechains.push_back(mmpaFrag.second);
  }
  return true;
}
}  // end of anonymous namespace

std::pair<JSMolList *, JSMolList *> JSMol::get_mmpa_frags(
    unsigned int minCuts, unsigned int maxCuts,
    unsigned int maxCutBonds) const {
  std::vector<RDKit::ROMOL_SPTR> cores;
  std::vector<RDKit::ROMOL_SPTR> sidechains;
  if (!mmpaFragmentMol(*d_mol, cores, sidechains, minCuts, maxCuts,
                       maxCutBonds)) {
    return std::make_pair(nullptr, nullptr);
  }
  return std::make_pair(new JSMolList(std::move(cores)),
                        new JSMolList(std::move(sidechains)));
}
#endif

#ifdef RDK_BUILD_MINIMAL_LIB_RXN
std::string JSReaction::get_svg(int w, int h) const {
  assert(d_rxn);
  return MinimalLib::rxn_to_svg(*d_rxn, w, h);
}
std::string JSReaction::get_svg_with_highlights(
    const std::string &details) const {
  assert(d_rxn);

  int w = d_defaultWidth;
  int h = d_defaultHeight;
  return MinimalLib::rxn_to_svg(*d_rxn, w, h, details);
}
bool JSReaction::is_valid() const { return true; }

std::vector<JSMolList *> JSReaction::run_reactants(
    const JSMolList &reactants, unsigned int maxProducts) const {
  d_rxn->initReactantMatchers();
  RDKit::MOL_SPTR_VECT reactant_vec;

  for (const auto &reactant : reactants.mols()) {
    if (!reactant) {
      throw ValueErrorException("Reactant must not be null");
    }
    reactant_vec.push_back(reactant);
  }

  std::vector<RDKit::MOL_SPTR_VECT> prods;
  prods = d_rxn->runReactants(reactant_vec, maxProducts);
  std::vector<JSMolList *> newResults;
  for (auto &mol_array : prods) {
    newResults.push_back(new JSMolList(mol_array));
  }
  return newResults;
}
#endif

JSMol *JSMolList::next() {
  JSMol *res = nullptr;
  if (d_idx < d_mols.size()) {
    res = at(d_idx++);
  }
  return res;
}

JSMol *JSMolList::at(size_t idx) const {
  JSMol *res = nullptr;
  if (idx < d_mols.size()) {
    const auto &molSptr = d_mols.at(idx);
    if (molSptr) {
      res = new JSMol(new RDKit::RWMol(*molSptr));
    }
  }
  return res;
}

JSMol *JSMolList::pop(size_t idx) {
  JSMol *res = nullptr;
  if (idx < d_mols.size()) {
    res = at(idx);
    d_mols.erase(d_mols.begin() + idx);
    if (d_idx > idx) {
      --d_idx;
    }
  }
  return res;
}

size_t JSMolList::append(const JSMol &mol) {
  d_mols.emplace_back(new ROMol(*mol.d_mol));
  return d_mols.size();
}

size_t JSMolList::insert(size_t idx, const JSMol &mol) {
  size_t res = 0;
  if (idx < d_mols.size()) {
    d_mols.emplace(d_mols.begin() + idx, new ROMol(*mol.d_mol));
    res = d_mols.size();
  }
  return res;
}

#ifdef RDK_BUILD_MINIMAL_LIB_SUBSTRUCTLIBRARY
JSSubstructLibrary::JSSubstructLibrary(unsigned int num_bits)
    : d_fpHolder(nullptr) {
  boost::shared_ptr<CachedTrustedSmilesMolHolder> molHolderSptr(
      new CachedTrustedSmilesMolHolder());
  boost::shared_ptr<PatternHolder> fpHolderSptr;
  d_molHolder = molHolderSptr.get();
  if (num_bits) {
    fpHolderSptr.reset(new PatternHolder(num_bits));
    d_fpHolder = fpHolderSptr.get();
    d_sslib.reset(new SubstructLibrary(molHolderSptr, fpHolderSptr));
  } else {
    d_sslib.reset(new SubstructLibrary(molHolderSptr));
  }
}

int JSSubstructLibrary::add_trusted_smiles(const std::string &smi) {
  SmilesParserParams ps;
  ps.sanitize = false;
  ps.removeHs = false;
  std::unique_ptr<RWMol> mol(SmilesToMol(smi, ps));
  if (!mol) {
    return -1;
  }
  mol->updatePropertyCache();
  MolOps::fastFindRings(*mol);
  int smiIdx;
  if (d_fpHolder) {
    auto fp = d_fpHolder->makeFingerprint(*mol);
    if (!fp) {
      return -1;
    }
    int fpIdx = d_fpHolder->addFingerprint(fp);
    smiIdx = d_molHolder->addSmiles(smi);
    CHECK_INVARIANT(fpIdx == smiIdx, "");
  } else {
    smiIdx = d_molHolder->addSmiles(smi);
  }
  return smiIdx;
}

int JSSubstructLibrary::add_trusted_smiles_and_pattern_fp(
    const std::string &smi, const std::string &patternFp) {
  if (!d_fpHolder) {
    throw ValueErrorException(NO_SUPPORT_FOR_PATTERN_FPS);
  }
  auto bitVect = new ExplicitBitVect(patternFp);
  if (!bitVect) {
    return -1;
  }
  int fpIdx = d_fpHolder->addFingerprint(bitVect);
  int smiIdx = d_molHolder->addSmiles(smi);
  CHECK_INVARIANT(fpIdx == smiIdx, "");
  return smiIdx;
}

std::string JSSubstructLibrary::get_trusted_smiles(unsigned int i) const {
  return d_molHolder->getMols().at(i);
}

std::string JSSubstructLibrary::get_pattern_fp(unsigned int i) const {
  if (!d_fpHolder) {
    throw ValueErrorException(NO_SUPPORT_FOR_PATTERN_FPS);
  }
  return d_fpHolder->getFingerprints().at(i)->toString();
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
  auto indices = d_sslib->getMatches(*q.d_mol, true, useChirality, false,
                                     numThreads, maxResults);
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
  return d_sslib->size() ? d_sslib->countMatches(*q.d_mol, true, useChirality,
                                                 false, numThreads)
                         : 0;
}
#endif

std::string get_inchikey_for_inchi(const std::string &input) {
  return InchiToInchiKey(input);
}

JSMol *get_mol_copy(const JSMol &other) {
  return new JSMol(new RWMol(*other.d_mol));
}

JSMol *get_mol(const std::string &input, const std::string &details_json) {
  auto mol = MinimalLib::mol_from_input(input, details_json);
  return mol ? new JSMol(mol) : nullptr;
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
  return mol ? new JSMol(mol) : nullptr;
}

JSMol *get_qmol(const std::string &input) {
  auto mol = MinimalLib::qmol_from_input(input);
  return mol ? new JSMol(mol) : nullptr;
}

#ifdef RDK_BUILD_MINIMAL_LIB_RXN
JSReaction *get_rxn(const std::string &input, const std::string &details_json) {
  auto rxn = MinimalLib::rxn_from_input(input, details_json);
  return rxn ? new JSReaction(rxn) : nullptr;
}
#endif

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

#ifdef RDK_BUILD_MINIMAL_LIB_MCS
namespace {
MCSResult getMcsResult(const JSMolList &molList,
                       const std::string &details_json) {
  MCSParameters p;
  if (!details_json.empty()) {
    parseMCSParametersJSON(details_json.c_str(), &p);
  }
  return RDKit::findMCS(molList.mols(), &p);
}
}  // namespace

std::string get_mcs_as_json(const JSMolList &molList,
                            const std::string &details_json) {
  auto mcsResult = getMcsResult(molList, details_json);
  rj::Document doc;
  doc.SetObject();
  auto &alloc = doc.GetAllocator();
  rj::Value rjSmarts;
  if (!mcsResult.DegenerateSmartsQueryMolDict.empty()) {
    rjSmarts.SetArray();
    for (const auto &pair : mcsResult.DegenerateSmartsQueryMolDict) {
      rjSmarts.PushBack(rj::Value(pair.first.c_str(), pair.first.size()),
                        alloc);
    }
  } else {
    rjSmarts.SetString(mcsResult.SmartsString.c_str(),
                       mcsResult.SmartsString.size());
  }
  doc.AddMember("smarts", rjSmarts, alloc);
  rj::Value rjCanceled;
  rjCanceled.SetBool(mcsResult.Canceled);
  doc.AddMember("canceled", rjCanceled, alloc);
  rj::Value rjNumAtoms;
  rjNumAtoms.SetInt(mcsResult.NumAtoms);
  doc.AddMember("numAtoms", rjNumAtoms, alloc);
  rj::Value rjNumBonds;
  rjNumBonds.SetInt(mcsResult.NumBonds);
  doc.AddMember("numBonds", rjNumBonds, alloc);
  rj::StringBuffer buffer;
  rj::Writer<rj::StringBuffer> writer(buffer);
  doc.Accept(writer);
  std::string res = buffer.GetString();
  return res;
}

std::string get_mcs_as_smarts(const JSMolList &molList,
                              const std::string &details_json) {
  auto res = getMcsResult(molList, details_json);
  return res.SmartsString;
}

JSMol *get_mcs_as_mol(const JSMolList &molList,
                      const std::string &details_json) {
  auto res = getMcsResult(molList, details_json);
  return new JSMol(new RWMol(*res.QueryMol));
}
#endif

RDKit::MinimalLib::LogHandle::LoggingFlag
    RDKit::MinimalLib::LogHandle::d_loggingNeedsInit = true;

JSLog::JSLog(RDKit::MinimalLib::LogHandle *logHandle) : d_logHandle(logHandle) {
  assert(d_logHandle);
}

JSLog::~JSLog() { delete d_logHandle; }

std::string JSLog::get_buffer() const { return d_logHandle->getBuffer(); }

void JSLog::clear_buffer() const { d_logHandle->clearBuffer(); }

JSLog *set_log_tee(const std::string &log_name) {
  auto logHandle = RDKit::MinimalLib::LogHandle::setLogTee(log_name.c_str());
  return logHandle ? new JSLog(logHandle) : nullptr;
}

JSLog *set_log_capture(const std::string &log_name) {
  auto logHandle =
      RDKit::MinimalLib::LogHandle::setLogCapture(log_name.c_str());
  return logHandle ? new JSLog(logHandle) : nullptr;
}

void enable_logging() { RDKit::MinimalLib::LogHandle::enableLogging(); }

void disable_logging() { RDKit::MinimalLib::LogHandle::disableLogging(); }
