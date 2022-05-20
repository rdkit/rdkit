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
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolDraw2D/MolDraw2DJS.h>

using namespace RDKit;

namespace RDKit {
namespace MinimalLib {
extern std::string process_details(const std::string &details, int &width,
                                   int &height, int &offsetx, int &offsety,
                                   std::string &legend,
                                   std::vector<int> &atomIds,
                                   std::vector<int> &bondIds, bool &kekulize);
}  // namespace MinimalLib
}  // namespace RDKit

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
  MolDraw2DJS *d2d = new MolDraw2DJS(width, height, ctx);
  d2d->setOffset(offsetx, offsety);
  MolDraw2DUtils::prepareAndDrawMolecule(*d2d, *self.d_mol);
  delete d2d;
  return "";
}

std::string draw_to_canvas(JSMol &self, emscripten::val canvas, int width,
                           int height) {
  return draw_to_canvas_with_offset(self, canvas, 0, 0, width, height);
}

std::string draw_to_canvas_with_highlights(JSMol &self, emscripten::val canvas,
                                           const std::string &details) {
  if (!self.d_mol) return "";

  std::vector<int> atomIds;
  std::vector<int> bondIds;

  auto ctx = canvas.call<emscripten::val>("getContext", std::string("2d"));

  int w = canvas["width"].as<int>();
  int h = canvas["height"].as<int>();
  int offsetx = 0;
  int offsety = 0;
  std::string legend = "";
  bool kekulize;
  auto problems = MinimalLib::process_details(
      details, w, h, offsetx, offsety, legend, atomIds, bondIds, kekulize);
  if (!problems.empty()) {
    return problems;
  }

  MolDraw2DJS *d2d = new MolDraw2DJS(w, h, ctx);
  if (!details.empty()) {
    MolDraw2DUtils::updateDrawerParamsFromJSON(*d2d, details);
  }
  d2d->setOffset(offsetx, offsety);

  MolDraw2DUtils::prepareAndDrawMolecule(*d2d, *self.d_mol, legend, &atomIds,
                                         &bondIds, nullptr, nullptr, nullptr,
                                         -1, kekulize);
  delete d2d;
  return "";
}

JSMol *get_mol_no_details(const std::string &input) {
  return get_mol(input, std::string());
}

emscripten::val binary_string_to_uint8array(const std::string &pkl) {
  emscripten::val view(emscripten::typed_memory_view(
      pkl.size(), reinterpret_cast<const unsigned char *>(pkl.c_str())));
  auto res = emscripten::val::global("Uint8Array").new_(pkl.size());
  res.call<void>("set", view);
  return res;
}

emscripten::val get_as_uint8array(const JSMol &self) {
  return binary_string_to_uint8array(self.get_pickle());
}

JSMol *get_mol_from_uint8array(const emscripten::val &pklAsUInt8Array) {
  return get_mol_from_pickle(pklAsUInt8Array.as<std::string>());
}

emscripten::val get_morgan_fp_as_uint8array(const JSMol &self,
                                            unsigned int radius,
                                            unsigned int fplen) {
  std::string fp = self.get_morgan_fp_as_binary_text(radius, fplen);
  return binary_string_to_uint8array(fp);
}

emscripten::val get_morgan_fp_as_uint8array(const JSMol &self) {
  return get_morgan_fp_as_uint8array(self, 2, 2048);
}

emscripten::val get_pattern_fp_as_uint8array(const JSMol &self,
                                             unsigned int fplen) {
  std::string fp = self.get_pattern_fp_as_binary_text(fplen);
  return binary_string_to_uint8array(fp);
}

emscripten::val get_pattern_fp_as_uint8array(const JSMol &self) {
  return get_pattern_fp_as_uint8array(self, 2048);
}

}  // namespace

using namespace emscripten;
EMSCRIPTEN_BINDINGS(RDKit_minimal) {
  class_<JSMol>("Mol")
      .function("is_valid", &JSMol::is_valid)
      .function("has_coords", &JSMol::has_coords)
      .function("get_smiles", &JSMol::get_smiles)
      .function("get_cxsmiles", &JSMol::get_cxsmiles)
      .function("get_smarts", &JSMol::get_smarts)
      .function("get_cxsmarts", &JSMol::get_cxsmarts)
      .function("get_molblock", &JSMol::get_molblock)
      .function("get_v3Kmolblock", &JSMol::get_v3Kmolblock)
      .function("get_as_uint8array", &get_as_uint8array)
      .function("get_inchi", &JSMol::get_inchi)
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
      .function("get_morgan_fp_as_uint8array",
                select_overload<emscripten::val(const JSMol &)>(
                    get_morgan_fp_as_uint8array))
      .function("get_morgan_fp_as_uint8array",
                select_overload<emscripten::val(const JSMol &, unsigned int,
                                                unsigned int)>(
                    get_morgan_fp_as_uint8array))
      .function("get_pattern_fp_as_uint8array",
                select_overload<emscripten::val(const JSMol &)>(
                    get_pattern_fp_as_uint8array))
      .function("get_pattern_fp_as_uint8array",
                select_overload<emscripten::val(const JSMol &, unsigned int)>(
                    get_pattern_fp_as_uint8array))
#endif
      .function("get_substruct_match", &JSMol::get_substruct_match)
      .function("get_substruct_matches", &JSMol::get_substruct_matches)
      .function("get_descriptors", &JSMol::get_descriptors)
      .function("get_morgan_fp",
                select_overload<std::string() const>(&JSMol::get_morgan_fp))
      .function("get_morgan_fp",
                select_overload<std::string(unsigned int, unsigned int) const>(
                    &JSMol::get_morgan_fp))
      .function("get_pattern_fp",
                select_overload<std::string() const>(&JSMol::get_pattern_fp))
      .function("get_pattern_fp",
                select_overload<std::string(unsigned int) const>(
                    &JSMol::get_pattern_fp))

      // functionality primarily useful in ketcher
      .function("get_stereo_tags", &JSMol::get_stereo_tags)
      .function("get_aromatic_form", &JSMol::get_aromatic_form)
      .function("get_kekule_form", &JSMol::get_kekule_form)
      .function("set_new_coords",
                select_overload<bool()>(&JSMol::set_new_coords))
      .function("get_new_coords",
                select_overload<std::string() const>(&JSMol::get_new_coords))
      .function("set_new_coords",
                select_overload<bool(bool)>(&JSMol::set_new_coords))
      .function("get_new_coords", select_overload<std::string(bool) const>(
                                      &JSMol::get_new_coords))
      .function("generate_aligned_coords",
                select_overload<std::string(const JSMol &)>(
                    &JSMol::generate_aligned_coords))
      .function("generate_aligned_coords",
                select_overload<std::string(const JSMol &, bool)>(
                    &JSMol::generate_aligned_coords))
      .function("generate_aligned_coords",
                select_overload<std::string(const JSMol &, bool, bool)>(
                    &JSMol::generate_aligned_coords))
      .function("generate_aligned_coords",
                select_overload<std::string(const JSMol &, bool, bool, bool)>(
                    &JSMol::generate_aligned_coords))
      .function("condense_abbreviations",
                select_overload<std::string()>(&JSMol::condense_abbreviations))
      .function("condense_abbreviations",
                select_overload<std::string(double, bool)>(
                    &JSMol::condense_abbreviations))
      .function("add_hs", &JSMol::add_hs)
      .function("remove_hs", &JSMol::remove_hs)
      .function("normalize_depiction",
                select_overload<double()>(&JSMol::normalize_depiction))
      .function("normalize_depiction",
                select_overload<double(int)>(&JSMol::normalize_depiction))
      .function("normalize_depiction", select_overload<double(int, double)>(
                                           &JSMol::normalize_depiction))
      .function("straighten_depiction", &JSMol::straighten_depiction);

  class_<JSSubstructLibrary>("SubstructLibrary")
      .constructor<>()
      .constructor<unsigned int>()
      .function("add_mol", &JSSubstructLibrary::add_mol)
      .function("add_smiles", &JSSubstructLibrary::add_smiles)
      .function("add_trusted_smiles", &JSSubstructLibrary::add_trusted_smiles)
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
                select_overload<unsigned int(const JSMol &) const>(
                    &JSSubstructLibrary::count_matches));

  function("version", &version);
  function("prefer_coordgen", &prefer_coordgen);
  function("use_legacy_stereo_perception", &use_legacy_stereo_perception);
  function("get_inchikey_for_inchi", &get_inchikey_for_inchi);
  function("get_mol", &get_mol, allow_raw_pointers());
  function("get_mol", &get_mol_no_details, allow_raw_pointers());
  function("get_mol_from_uint8array", &get_mol_from_uint8array,
           allow_raw_pointers());
  function("get_mol_copy", &get_mol_copy, allow_raw_pointers());
  function("get_qmol", &get_qmol, allow_raw_pointers());
}
