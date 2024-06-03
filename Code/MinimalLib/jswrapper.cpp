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
#include <emscripten.h>
#include <emscripten/val.h>
#include <emscripten/bind.h>
#include "minilib.h"
#include "common_defs.h"
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolDraw2D/MolDraw2DJS.h>

using namespace RDKit;

namespace {
std::string draw_to_canvas_with_offset(JSMol &self, emscripten::val canvas,
                                       int offsetx, int offsety, int width,
                                       int height) {
  if (!self.d_mol) {
    return "no molecule";
  }
  auto ctx = canvas.call<emscripten::val>("getContext", std::string("2d"));
  if (width < 0) {
    width = canvas["width"].as<int>();
  }
  if (height < 0) {
    height = canvas["height"].as<int>();
  }
  std::unique_ptr<MolDraw2DJS> d2d(new MolDraw2DJS(width, height, ctx));
  d2d->setOffset(offsetx, offsety);
  MolDraw2DUtils::prepareAndDrawMolecule(*d2d, *self.d_mol);
  return "";
}

std::string draw_to_canvas(JSMol &self, emscripten::val canvas, int width,
                           int height) {
  return draw_to_canvas_with_offset(self, canvas, 0, 0, width, height);
}

std::string draw_to_canvas_with_highlights(JSMol &self, emscripten::val canvas,
                                           const std::string &details) {
  if (!self.d_mol) {
    return "no molecule";
  }

  auto ctx = canvas.call<emscripten::val>("getContext", std::string("2d"));
  MinimalLib::MolDrawingDetails molDrawingDetails;
  molDrawingDetails.width = canvas["width"].as<int>();
  molDrawingDetails.height = canvas["height"].as<int>();
  if (!details.empty()) {
    auto problems = MinimalLib::process_mol_details(details, molDrawingDetails);
    if (!problems.empty()) {
      return problems;
    }
  }

  std::unique_ptr<MolDraw2DJS> d2d(new MolDraw2DJS(
      molDrawingDetails.width, molDrawingDetails.height, ctx,
      molDrawingDetails.panelWidth, molDrawingDetails.panelHeight,
      molDrawingDetails.noFreetype));
  if (!details.empty()) {
    MolDraw2DUtils::updateDrawerParamsFromJSON(*d2d, details);
  }
  d2d->setOffset(molDrawingDetails.offsetx, molDrawingDetails.offsety);

  MolDraw2DUtils::prepareAndDrawMolecule(
      *d2d, *self.d_mol, molDrawingDetails.legend, &molDrawingDetails.atomIds,
      &molDrawingDetails.bondIds,
      molDrawingDetails.atomMap.empty() ? nullptr : &molDrawingDetails.atomMap,
      molDrawingDetails.bondMap.empty() ? nullptr : &molDrawingDetails.bondMap,
      molDrawingDetails.radiiMap.empty() ? nullptr
                                         : &molDrawingDetails.radiiMap,
      -1, molDrawingDetails.kekulize, molDrawingDetails.addChiralHs,
      molDrawingDetails.wedgeBonds, molDrawingDetails.forceCoords,
      molDrawingDetails.wavyBonds);
  return "";
}

#ifdef RDK_BUILD_MINIMAL_LIB_RXN
std::string draw_rxn_to_canvas_with_offset(JSReaction &self,
                                           emscripten::val canvas, int offsetx,
                                           int offsety, int width, int height) {
  if (!self.d_rxn) {
    return "no reaction";
  }
  auto ctx = canvas.call<emscripten::val>("getContext", std::string("2d"));
  if (width < 0) {
    width = canvas["width"].as<int>();
  }
  if (height < 0) {
    height = canvas["height"].as<int>();
  }
  std::unique_ptr<MolDraw2DJS> d2d(new MolDraw2DJS(width, height, ctx));
  d2d->setOffset(offsetx, offsety);
  d2d->drawReaction(*self.d_rxn);
  return "";
}

std::string draw_rxn_to_canvas(JSReaction &self, emscripten::val canvas,
                               int width, int height) {
  return draw_rxn_to_canvas_with_offset(self, canvas, 0, 0, width, height);
}

std::string draw_rxn_to_canvas_with_highlights(JSReaction &self,
                                               emscripten::val canvas,
                                               const std::string &details) {
  if (!self.d_rxn) {
    return "no reaction";
  }

  auto ctx = canvas.call<emscripten::val>("getContext", std::string("2d"));
  int w = canvas["width"].as<int>();
  int h = canvas["height"].as<int>();
  MinimalLib::RxnDrawingDetails rxnDrawingDetails;
  if (!details.empty()) {
    auto problems = MinimalLib::process_rxn_details(details, rxnDrawingDetails);
    if (!problems.empty()) {
      return problems;
    }
  }

  std::unique_ptr<MolDraw2DJS> d2d(new MolDraw2DJS(w, h, ctx));
  if (!details.empty()) {
    MolDraw2DUtils::updateDrawerParamsFromJSON(*d2d, details);
  }
  d2d->setOffset(rxnDrawingDetails.offsetx, rxnDrawingDetails.offsety);
  if (!rxnDrawingDetails.kekulize) {
    d2d->drawOptions().prepareMolsBeforeDrawing = false;
  }
  d2d->drawReaction(*self.d_rxn, rxnDrawingDetails.highlightByReactant,
                    !rxnDrawingDetails.highlightByReactant ||
                            rxnDrawingDetails.highlightColorsReactants.empty()
                        ? nullptr
                        : &rxnDrawingDetails.highlightColorsReactants);
  return "";
}
#endif

JSMol *get_mol_no_details(const std::string &input) {
  return get_mol(input, std::string());
}

#ifdef RDK_BUILD_MINIMAL_LIB_MCS
std::string get_mcs_as_json_no_details(const JSMolList &mols) {
  return get_mcs_as_json(mols, std::string());
}

JSMol *get_mcs_as_mol_no_details(const JSMolList &mols) {
  return get_mcs_as_mol(mols, std::string());
}

std::string get_mcs_as_smarts_no_details(const JSMolList &mols) {
  return get_mcs_as_smarts(mols, std::string());
}
#endif

emscripten::val binary_string_to_uint8array(const std::string &pkl) {
  emscripten::val view(emscripten::typed_memory_view(
      pkl.size(), reinterpret_cast<const unsigned char *>(pkl.c_str())));
  auto res = emscripten::val::global("Uint8Array").new_(pkl.size());
  res.call<void>("set", view);
  return res;
}

emscripten::val uint_vector_to_uint32array(
    const std::vector<unsigned int> &vec) {
  auto res = emscripten::val::global("Uint32Array").new_(vec.size());
  if (!vec.empty()) {
    emscripten::val view(emscripten::typed_memory_view(
        vec.size(), reinterpret_cast<const unsigned int *>(vec.data())));
    res.call<void>("set", view);
  }
  return res;
}

emscripten::val get_as_uint8array(const JSMol &self) {
  return binary_string_to_uint8array(self.get_pickle());
}

JSMol *get_mol_from_uint8array(const emscripten::val &pklAsUInt8Array) {
  return get_mol_from_pickle(pklAsUInt8Array.as<std::string>());
}

#ifdef RDK_BUILD_MINIMAL_LIB_RXN
JSReaction *get_rxn_no_details(const std::string &input) {
  return get_rxn(input, std::string());
}
#endif

std::string generate_aligned_coords_helper(JSMol &self,
                                           const JSMol &templateMol,
                                           const emscripten::val &param) {
  if (param.typeOf().as<std::string>() != "string") {
    throw std::runtime_error(
        "generate_aligned_coords expects a JSON string parameter");
  }
  return self.generate_aligned_coords(templateMol, param.as<std::string>());
}

emscripten::val get_morgan_fp_as_uint8array(const JSMol &self,
                                            const std::string &details) {
  auto fp = self.get_morgan_fp_as_binary_text(details);
  return binary_string_to_uint8array(fp);
}

emscripten::val get_morgan_fp_as_uint8array(const JSMol &self) {
  return get_morgan_fp_as_uint8array(self, "{}");
}

std::string parse_pattern_fp_param(const emscripten::val &param,
                                   const std::string &funcName) {
  std::string details;
  if (param.typeOf().as<std::string>() == "string") {
    details = param.as<std::string>();
  } else {
    throw std::runtime_error(
        (funcName + " expects a JSON string as parameter").c_str());
  }
  return details;
}

std::string get_pattern_fp_helper(const JSMol &self,
                                  const emscripten::val &param) {
  auto details = parse_pattern_fp_param(param, "get_pattern_fp");
  return self.get_pattern_fp(details);
}

emscripten::val get_pattern_fp_as_uint8array_helper(
    const JSMol &self, const emscripten::val &param) {
  auto details = parse_pattern_fp_param(param, "get_pattern_fp_as_uint8array");
  auto fp = self.get_pattern_fp_as_binary_text(details);
  return binary_string_to_uint8array(fp);
}

emscripten::val get_pattern_fp_as_uint8array(const JSMol &self) {
  auto fp = self.get_pattern_fp_as_binary_text("{}");
  return binary_string_to_uint8array(fp);
}

emscripten::val get_topological_torsion_fp_as_uint8array(
    const JSMol &self, const std::string &details) {
  auto fp = self.get_topological_torsion_fp_as_binary_text(details);
  return binary_string_to_uint8array(fp);
}

emscripten::val get_topological_torsion_fp_as_uint8array(const JSMol &self) {
  return get_topological_torsion_fp_as_uint8array(self, "{}");
}

emscripten::val get_rdkit_fp_as_uint8array(const JSMol &self,
                                           const std::string &details) {
  auto fp = self.get_rdkit_fp_as_binary_text(details);
  return binary_string_to_uint8array(fp);
}

emscripten::val get_rdkit_fp_as_uint8array(const JSMol &self) {
  return get_rdkit_fp_as_uint8array(self, "{}");
}

emscripten::val get_atom_pair_fp_as_uint8array(const JSMol &self,
                                               const std::string &details) {
  auto fp = self.get_atom_pair_fp_as_binary_text(details);
  return binary_string_to_uint8array(fp);
}

emscripten::val get_atom_pair_fp_as_uint8array(const JSMol &self) {
  return get_atom_pair_fp_as_uint8array(self, "{}");
}

emscripten::val get_maccs_fp_as_uint8array(const JSMol &self) {
  auto fp = self.get_maccs_fp_as_binary_text();
  return binary_string_to_uint8array(fp);
}

emscripten::val get_frags_helper(JSMol &self, const std::string &details) {
  auto res = self.get_frags(details);
  auto obj = emscripten::val::object();
  obj.set("molList", res.first);
  obj.set("mappings", res.second);
  return obj;
}

emscripten::val get_frags_helper(JSMol &self) {
  return get_frags_helper(self, "{}");
}

#ifdef RDK_BUILD_MINIMAL_LIB_SUBSTRUCTLIBRARY
int add_trusted_smiles_and_pattern_fp_helper(
    JSSubstructLibrary &self, const std::string &smi,
    const emscripten::val &patternFpAsUInt8Array) {
  return self.add_trusted_smiles_and_pattern_fp(
      smi, patternFpAsUInt8Array.as<std::string>());
}

emscripten::val get_pattern_fp_as_uint8array_from_sslib(
    const JSSubstructLibrary &self, unsigned int i) {
  return binary_string_to_uint8array(self.get_pattern_fp(i));
}

emscripten::val get_matches_as_uint32array(const JSSubstructLibrary &self,
                                           const JSMol &q, bool useChirality,
                                           int numThreads, int maxResults) {
  auto indices = self.d_sslib->size()
                     ? self.d_sslib->getMatches(*q.d_mol, true, useChirality,
                                                false, numThreads, maxResults)
                     : std::vector<unsigned int>();
  return uint_vector_to_uint32array(indices);
}

emscripten::val get_matches_as_uint32array(const JSSubstructLibrary &self,
                                           const JSMol &q, int maxResults) {
  return get_matches_as_uint32array(self, q, self.d_defaultUseChirality,
                                    self.d_defaultNumThreads, maxResults);
}

emscripten::val get_matches_as_uint32array(const JSSubstructLibrary &self,
                                           const JSMol &q) {
  return get_matches_as_uint32array(self, q, self.d_defaultUseChirality,
                                    self.d_defaultNumThreads,
                                    self.d_defaultMaxResults);
}
#endif

#ifdef RDK_BUILD_AVALON_SUPPORT
emscripten::val get_avalon_fp_as_uint8array(const JSMol &self,
                                            const std::string &details) {
  auto fp = self.get_avalon_fp_as_binary_text(details);
  return binary_string_to_uint8array(fp);
}

emscripten::val get_avalon_fp_as_uint8array(const JSMol &self) {
  return get_avalon_fp_as_uint8array(self, "{}");
}
#endif

#ifdef RDK_BUILD_MINIMAL_LIB_MMPA
emscripten::val get_mmpa_frags_helper(const JSMol &self, unsigned int minCuts,
                                      unsigned int maxCuts,
                                      unsigned int maxCutBonds) {
  auto obj = emscripten::val::object();
  auto pairs = self.get_mmpa_frags(minCuts, maxCuts, maxCutBonds);
  obj.set("cores", pairs.first);
  obj.set("sidechains", pairs.second);
  return obj;
}
#endif

}  // namespace

using namespace emscripten;
EMSCRIPTEN_BINDINGS(RDKit_minimal) {
  register_vector<std::string>("StringList");
  register_vector<JSMolList *>("JSMolListList");

  class_<JSMol>("Mol")
      .function("is_valid", &JSMol::is_valid)
      .function("has_coords", &JSMol::has_coords)
      .function("get_smiles",
                select_overload<std::string() const>(&JSMol::get_smiles))
      .function("get_smiles",
                select_overload<std::string(const std::string &) const>(
                    &JSMol::get_smiles))
      .function("get_cxsmiles",
                select_overload<std::string() const>(&JSMol::get_cxsmiles))
      .function("get_cxsmiles",
                select_overload<std::string(const std::string &) const>(
                    &JSMol::get_cxsmiles))
      .function("get_smarts",
                select_overload<std::string() const>(&JSMol::get_smarts))
      .function("get_smarts",
                select_overload<std::string(const std::string &) const>(
                    &JSMol::get_smarts))
      .function("get_cxsmarts",
                select_overload<std::string() const>(&JSMol::get_cxsmarts))
      .function("get_cxsmarts",
                select_overload<std::string(const std::string &) const>(
                    &JSMol::get_cxsmarts))
      .function("get_molblock",
                select_overload<std::string() const>(&JSMol::get_molblock))
      .function("get_molblock",
                select_overload<std::string(const std::string &) const>(
                    &JSMol::get_molblock))
      .function("get_v3Kmolblock",
                select_overload<std::string() const>(&JSMol::get_v3Kmolblock))
      .function("get_v3Kmolblock",
                select_overload<std::string(const std::string &) const>(
                    &JSMol::get_v3Kmolblock))
      .function("get_as_uint8array", &get_as_uint8array)
      .function("get_inchi", select_overload<std::string(const std::string&) const>(&JSMol::get_inchi))
      .function("get_inchi", select_overload<std::string() const>(&JSMol::get_inchi))
      .function("get_json", &JSMol::get_json)
      .function("get_svg",
                select_overload<std::string() const>(&JSMol::get_svg))
      .function("get_svg",
                select_overload<std::string(int, int) const>(&JSMol::get_svg))

      .function("get_svg_with_highlights", &JSMol::get_svg_with_highlights)
#ifdef __EMSCRIPTEN__
      .function("draw_to_canvas_with_offset", &draw_to_canvas_with_offset)
      .function("draw_to_canvas", &draw_to_canvas)
      .function("draw_to_canvas_with_highlights",
                &draw_to_canvas_with_highlights)
      .function("generate_aligned_coords",
                select_overload<std::string(JSMol &, const JSMol &,
                                            const emscripten::val &)>(
                    generate_aligned_coords_helper))
      .function("get_morgan_fp_as_uint8array",
                select_overload<emscripten::val(const JSMol &)>(
                    get_morgan_fp_as_uint8array))
      .function(
          "get_morgan_fp_as_uint8array",
          select_overload<emscripten::val(const JSMol &, const std::string &)>(
              get_morgan_fp_as_uint8array))
      .function(
          "get_pattern_fp",
          select_overload<std::string(const JSMol &, const emscripten::val &)>(
              get_pattern_fp_helper))
      .function("get_pattern_fp_as_uint8array",
                select_overload<emscripten::val(const JSMol &)>(
                    get_pattern_fp_as_uint8array))
      .function("get_pattern_fp_as_uint8array",
                select_overload<emscripten::val(const JSMol &,
                                                const emscripten::val &)>(
                    get_pattern_fp_as_uint8array_helper))
      .function("get_topological_torsion_fp_as_uint8array",
                select_overload<emscripten::val(const JSMol &)>(
                    get_topological_torsion_fp_as_uint8array))
      .function(
          "get_topological_torsion_fp_as_uint8array",
          select_overload<emscripten::val(const JSMol &, const std::string &)>(
              get_topological_torsion_fp_as_uint8array))
      .function("get_rdkit_fp_as_uint8array",
                select_overload<emscripten::val(const JSMol &)>(
                    get_rdkit_fp_as_uint8array))
      .function(
          "get_rdkit_fp_as_uint8array",
          select_overload<emscripten::val(const JSMol &, const std::string &)>(
              get_rdkit_fp_as_uint8array))
      .function("get_atom_pair_fp_as_uint8array",
                select_overload<emscripten::val(const JSMol &)>(
                    get_atom_pair_fp_as_uint8array))
      .function(
          "get_atom_pair_fp_as_uint8array",
          select_overload<emscripten::val(const JSMol &, const std::string &)>(
              get_atom_pair_fp_as_uint8array))
      .function("get_maccs_fp_as_uint8array", &get_maccs_fp_as_uint8array)
      .function("get_frags",
                select_overload<emscripten::val(JSMol &, const std::string &)>(
                    get_frags_helper),
                allow_raw_pointers())
      .function("get_frags",
                select_overload<emscripten::val(JSMol &)>(get_frags_helper),
                allow_raw_pointers())
#ifdef RDK_BUILD_AVALON_SUPPORT
      .function("get_avalon_fp_as_uint8array",
                select_overload<emscripten::val(const JSMol &)>(
                    get_avalon_fp_as_uint8array))
      .function(
          "get_avalon_fp_as_uint8array",
          select_overload<emscripten::val(const JSMol &, const std::string &)>(
              get_avalon_fp_as_uint8array))
#endif
#endif
      .function("get_substruct_match", &JSMol::get_substruct_match)
      .function("get_substruct_matches", &JSMol::get_substruct_matches)
      .function("get_descriptors", &JSMol::get_descriptors)
      .function("get_morgan_fp",
                select_overload<std::string() const>(&JSMol::get_morgan_fp))
      .function("get_morgan_fp",
                select_overload<std::string(const std::string &) const>(
                    &JSMol::get_morgan_fp))
      .function("get_pattern_fp",
                select_overload<std::string() const>(&JSMol::get_pattern_fp))
      .function("get_topological_torsion_fp",
                select_overload<std::string() const>(
                    &JSMol::get_topological_torsion_fp))
      .function("get_topological_torsion_fp",
                select_overload<std::string(const std::string &) const>(
                    &JSMol::get_topological_torsion_fp))
      .function("get_rdkit_fp",
                select_overload<std::string() const>(&JSMol::get_rdkit_fp))
      .function("get_rdkit_fp",
                select_overload<std::string(const std::string &) const>(
                    &JSMol::get_rdkit_fp))
      .function("get_atom_pair_fp",
                select_overload<std::string() const>(&JSMol::get_atom_pair_fp))
      .function("get_atom_pair_fp",
                select_overload<std::string(const std::string &) const>(
                    &JSMol::get_atom_pair_fp))
      .function("get_maccs_fp", &JSMol::get_maccs_fp)
#ifdef RDK_BUILD_AVALON_SUPPORT
      .function("get_avalon_fp",
                select_overload<std::string() const>(&JSMol::get_avalon_fp))
      .function("get_avalon_fp",
                select_overload<std::string(const std::string &) const>(
                    &JSMol::get_avalon_fp))
#endif

      // functionality primarily useful in ketcher
      .function("get_stereo_tags", &JSMol::get_stereo_tags)
      .function("get_aromatic_form", &JSMol::get_aromatic_form)
      .function("convert_to_aromatic_form", &JSMol::convert_to_aromatic_form)
      .function("get_kekule_form", &JSMol::get_kekule_form)
      .function("convert_to_kekule_form", &JSMol::convert_to_kekule_form)
      .function("set_new_coords",
                select_overload<bool()>(&JSMol::set_new_coords))
      .function("get_new_coords",
                select_overload<std::string() const>(&JSMol::get_new_coords))
      .function("set_new_coords",
                select_overload<bool(bool)>(&JSMol::set_new_coords))
      .function("get_new_coords", select_overload<std::string(bool) const>(
                                      &JSMol::get_new_coords))
      .function("has_prop", &JSMol::has_prop)
      .function("get_prop_list",
                select_overload<std::vector<std::string>(
                    bool includePrivate, bool includeComputed) const>(
                    &JSMol::get_prop_list))
      .function(
          "get_prop_list",
          select_overload<std::vector<std::string>(bool includePrivate) const>(
              &JSMol::get_prop_list))
      .function("get_prop_list",
                select_overload<std::vector<std::string>() const>(
                    &JSMol::get_prop_list))
      .function(
          "set_prop",
          select_overload<bool(const std::string &, const std::string &, bool)>(
              &JSMol::set_prop))
      .function("set_prop",
                select_overload<bool(const std::string &, const std::string &)>(
                    &JSMol::set_prop))
      .function("get_prop", &JSMol::get_prop)
      .function("clear_prop", &JSMol::clear_prop)
      .function("condense_abbreviations",
                select_overload<std::string()>(&JSMol::condense_abbreviations))
      .function("condense_abbreviations",
                select_overload<std::string(double, bool)>(
                    &JSMol::condense_abbreviations))
      .function("add_hs", &JSMol::add_hs)
      .function("add_hs_in_place", &JSMol::add_hs_in_place)
      .function("remove_hs", &JSMol::remove_hs)
      .function("remove_hs_in_place", &JSMol::remove_hs_in_place)
      .function("normalize_depiction",
                select_overload<double()>(&JSMol::normalize_depiction))
      .function("normalize_depiction",
                select_overload<double(int)>(&JSMol::normalize_depiction))
      .function("normalize_depiction", select_overload<double(int, double)>(
                                           &JSMol::normalize_depiction))
      .function("straighten_depiction",
                select_overload<void()>(&JSMol::straighten_depiction))
      .function("straighten_depiction",
                select_overload<void(bool)>(&JSMol::straighten_depiction))
      .function("get_num_atoms", select_overload<unsigned int(bool) const>(
                                     &JSMol::get_num_atoms))
      .function("get_num_atoms",
                select_overload<unsigned int() const>(&JSMol::get_num_atoms))
      .function("get_num_bonds", &JSMol::get_num_bonds)
#ifdef RDK_BUILD_MINIMAL_LIB_MMPA
      .function("get_mmpa_frags",
                select_overload<emscripten::val(const JSMol &, unsigned int,
                                                unsigned int, unsigned int)>(
                    get_mmpa_frags_helper))
#endif
      ;

  class_<JSMolList>("MolList")
      .constructor<>()
      .function("append", &JSMolList::append)
      .function("insert", &JSMolList::insert)
      .function("at", &JSMolList::at, allow_raw_pointers())
      .function("pop", &JSMolList::pop, allow_raw_pointers())
      .function("next", &JSMolList::next, allow_raw_pointers())
      .function("reset", &JSMolList::reset)
      .function("at_end", &JSMolList::at_end)
      .function("size", &JSMolList::size);

#ifdef RDK_BUILD_MINIMAL_LIB_RXN
  class_<JSReaction>("Reaction")
#ifdef __EMSCRIPTEN__
      .function("run_reactants", select_overload<std::vector<JSMolList *>(
                                     const JSMolList &, unsigned int) const>(
                                     &JSReaction::run_reactants))
      .function("draw_to_canvas_with_offset", &draw_rxn_to_canvas_with_offset)
      .function("draw_to_canvas", &draw_rxn_to_canvas)
      .function("draw_to_canvas_with_highlights",
                &draw_rxn_to_canvas_with_highlights)
#endif
      .function("get_svg",
                select_overload<std::string() const>(&JSReaction::get_svg))
      .function("get_svg", select_overload<std::string(int, int) const>(
                               &JSReaction::get_svg))

      .function("get_svg_with_highlights",
                &JSReaction::get_svg_with_highlights);
#endif

#ifdef RDK_BUILD_MINIMAL_LIB_SUBSTRUCTLIBRARY
  class_<JSSubstructLibrary>("SubstructLibrary")
      .constructor<>()
      .constructor<unsigned int>()
      .function("add_mol", &JSSubstructLibrary::add_mol)
      .function("add_smiles", &JSSubstructLibrary::add_smiles)
      .function("add_trusted_smiles", &JSSubstructLibrary::add_trusted_smiles)
      .function("get_trusted_smiles", &JSSubstructLibrary::get_trusted_smiles)
#ifdef __EMSCRIPTEN__
      .function("add_trusted_smiles_and_pattern_fp",
                select_overload<int(JSSubstructLibrary &, const std::string &,
                                    const emscripten::val &)>(
                    add_trusted_smiles_and_pattern_fp_helper))
      .function("get_pattern_fp_as_uint8array",
                select_overload<emscripten::val(const JSSubstructLibrary &,
                                                unsigned int)>(
                    get_pattern_fp_as_uint8array_from_sslib))
      .function("get_matches_as_uint32array",
                select_overload<emscripten::val(const JSSubstructLibrary &,
                                                const JSMol &, bool, int, int)>(
                    get_matches_as_uint32array))
      .function("get_matches_as_uint32array",
                select_overload<emscripten::val(const JSSubstructLibrary &,
                                                const JSMol &, int)>(
                    get_matches_as_uint32array))
      .function("get_matches_as_uint32array",
                select_overload<emscripten::val(const JSSubstructLibrary &,
                                                const JSMol &)>(
                    get_matches_as_uint32array))
#endif
      .function("get_mol", &JSSubstructLibrary::get_mol, allow_raw_pointers())
      .function(
          "get_matches",
          select_overload<std::string(const JSMol &, bool, int, int) const>(
              &JSSubstructLibrary::get_matches))
      .function("get_matches",
                select_overload<std::string(const JSMol &, int) const>(
                    &JSSubstructLibrary::get_matches))
      .function("get_matches",
                select_overload<std::string(const JSMol &) const>(
                    &JSSubstructLibrary::get_matches))
      .function("count_matches",
                select_overload<unsigned int(const JSMol &, bool, int) const>(
                    &JSSubstructLibrary::count_matches))
      .function("count_matches",
                select_overload<unsigned int(const JSMol &, bool) const>(
                    &JSSubstructLibrary::count_matches))
      .function("count_matches",
                select_overload<unsigned int(const JSMol &) const>(
                    &JSSubstructLibrary::count_matches))
      .function("size", &JSSubstructLibrary::size)
#endif
      ;

  class_<JSLog>("Log")
      .function("get_buffer", &JSLog::get_buffer)
      .function("clear_buffer", &JSLog::clear_buffer);

  function("version", &version);
  function("prefer_coordgen", &prefer_coordgen);
  function("use_legacy_stereo_perception", &use_legacy_stereo_perception);
  function("allow_non_tetrahedral_chirality", &allow_non_tetrahedral_chirality);
  function("get_inchikey_for_inchi", &get_inchikey_for_inchi);
  function("get_mol", &get_mol, allow_raw_pointers());
  function("get_mol", &get_mol_no_details, allow_raw_pointers());
  function("get_mol_from_uint8array", &get_mol_from_uint8array,
           allow_raw_pointers());
  function("get_mol_copy", &get_mol_copy, allow_raw_pointers());
  function("get_qmol", &get_qmol, allow_raw_pointers());
  function("enable_logging", &enable_logging);
  function("disable_logging", &disable_logging);
  function("set_log_capture", &set_log_capture, allow_raw_pointers());
  function("set_log_tee", &set_log_tee, allow_raw_pointers());
#ifdef RDK_BUILD_MINIMAL_LIB_RXN
  function("get_rxn", &get_rxn, allow_raw_pointers());
  function("get_rxn", &get_rxn_no_details, allow_raw_pointers());
#endif
#ifdef RDK_BUILD_MINIMAL_LIB_MCS
  function("get_mcs_as_json", &get_mcs_as_json);
  function("get_mcs_as_json", &get_mcs_as_json_no_details);
  function("get_mcs_as_mol", &get_mcs_as_mol, allow_raw_pointers());
  function("get_mcs_as_mol", &get_mcs_as_mol_no_details, allow_raw_pointers());
  function("get_mcs_as_smarts", &get_mcs_as_smarts);
  function("get_mcs_as_smarts", &get_mcs_as_smarts_no_details);
#endif
}
