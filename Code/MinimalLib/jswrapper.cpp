//
//
//  Copyright (C) 2019-2020 Greg Landrum
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

extern std::string process_details(const std::string &details,
                                   unsigned int &width, unsigned int &height,
                                   int &offsetx, int &offsety,
                                   std::string &legend,
                                   std::vector<int> &atomIds,
                                   std::vector<int> &bondIds);

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

  unsigned int w = canvas["width"].as<unsigned int>();
  unsigned int h = canvas["height"].as<unsigned int>();
  int offsetx = 0;
  int offsety = 0;
  std::string legend = "";
  auto problems = process_details(details, w, h, offsetx, offsety, legend,
                                  atomIds, bondIds);
  if (!problems.empty()) {
    return problems;
  }

  MolDraw2DJS *d2d = new MolDraw2DJS(w, h, ctx);
  if (!details.empty()) {
    MolDraw2DUtils::updateDrawerParamsFromJSON(*d2d, details);
  }
  d2d->setOffset(offsetx, offsety);

  MolDraw2DUtils::prepareAndDrawMolecule(*d2d, *self.d_mol, legend, &atomIds,
                                         &bondIds);
  delete d2d;
  return "";
}

}  // namespace

using namespace emscripten;
EMSCRIPTEN_BINDINGS(RDKit_minimal) {
  class_<JSMol>("Mol")
      .function("is_valid", &JSMol::is_valid)
      .function("get_smiles", &JSMol::get_smiles)
      .function("get_cxsmiles", &JSMol::get_cxsmiles)
      .function("get_molblock", &JSMol::get_molblock)
      .function("get_v3Kmolblock", &JSMol::get_v3Kmolblock)
      .function("get_inchi", &JSMol::get_inchi)
      .function("get_svg",
                select_overload<std::string() const>(&JSMol::get_svg))
      .function("get_svg",
                select_overload<std::string(unsigned int, unsigned int) const>(
                    &JSMol::get_svg))

      .function("get_svg_with_highlights", &JSMol::get_svg_with_highlights)
#ifdef __EMSCRIPTEN__
      .function("draw_to_canvas_with_offset", &draw_to_canvas_with_offset)
      .function("draw_to_canvas", &draw_to_canvas)
      .function("draw_to_canvas_with_highlights",
                &draw_to_canvas_with_highlights)
#endif
      .function("get_substruct_match", &JSMol::get_substruct_match)
      .function("get_substruct_matches", &JSMol::get_substruct_matches)
      .function("get_descriptors", &JSMol::get_descriptors)
      .function("get_morgan_fp",
                select_overload<std::string() const>(&JSMol::get_morgan_fp))
      .function("get_morgan_fp",
                select_overload<std::string(unsigned int, unsigned int) const>(
                    &JSMol::get_morgan_fp))

      // functionality primarily useful in ketcher
      .function("get_stereo_tags", &JSMol::get_stereo_tags)
      .function("get_aromatic_form", &JSMol::get_aromatic_form)
      .function("get_kekule_form", &JSMol::get_kekule_form)
      .function("get_new_coords",
                select_overload<std::string() const>(&JSMol::get_new_coords))
      .function("get_new_coords", select_overload<std::string(bool) const>(
                                      &JSMol::get_new_coords))
      .function("generate_aligned_coords",
                select_overload<std::string(const JSMol &)>(&JSMol::generate_aligned_coords))
      .function("generate_aligned_coords", select_overload<std::string(const JSMol &,bool)>(
                                      &JSMol::generate_aligned_coords))
      .function("condense_abbreviations",
                select_overload<std::string()>(&JSMol::condense_abbreviations))
      .function("condense_abbreviations",
                select_overload<std::string(double, bool)>(
                    &JSMol::condense_abbreviations))
      .function("add_hs", &JSMol::add_hs)
      .function("remove_hs", &JSMol::remove_hs);

  function("version", &version);
  function("prefer_coordgen", &prefer_coordgen);
  function("get_inchikey_for_inchi", &get_inchikey_for_inchi);
  function("get_mol", &get_mol, allow_raw_pointers());
  function("get_qmol", &get_qmol, allow_raw_pointers());
}