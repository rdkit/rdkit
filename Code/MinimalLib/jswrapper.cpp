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
#include <emscripten/bind.h>
#include "minilib.h"

using namespace emscripten;
EMSCRIPTEN_BINDINGS(RDKit_minimal) {
  register_vector<unsigned int>("VectorUInt");

  class_<JSMol>("Mol")
      .function("is_valid", &JSMol::is_valid)
      .function("get_smiles", &JSMol::get_smiles)
      .function("get_molblock", &JSMol::get_molblock)
      .function("get_inchi", &JSMol::get_inchi)
      .function("get_svg", &JSMol::get_svg)
      .function("get_svg_with_highlights", &JSMol::get_svg_with_highlights)
      .function("get_substruct_match", &JSMol::get_substruct_match)
      .function("get_descriptors", &JSMol::get_descriptors);

  function("version", &version);
  function("get_inchikey_for_inchi", &get_inchikey_for_inchi);
  function("get_mol", &get_mol, allow_raw_pointers());
  function("get_qmol", &get_qmol, allow_raw_pointers());
}