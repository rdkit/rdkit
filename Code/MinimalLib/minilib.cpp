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
#include <string>
#include <iostream>

#if 1
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

std::string canon_smiles(const std::string &input) {
  std::unique_ptr<ROMol> mol(SmilesToMol(input));
  if (!mol) return "";
  return MolToSmiles(*mol);
}
#endif
int ping() {
  std::cerr << "alive" << std::endl;
  return 1;
}

using namespace emscripten;
EMSCRIPTEN_BINDINGS(RDKit_minimal) {
  function("ping", &ping);
  function("canon_smiles", &canon_smiles);
}