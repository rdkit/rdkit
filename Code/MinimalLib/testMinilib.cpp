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

#include <emscripten.h>
#include <iostream>
#include <string>
#include <MinimalLib/minilib.h>

EMSCRIPTEN_KEEPALIVE
int lping() {
  std::cerr << "blah blah" << std::endl;
  return 2;
}

int main() {
  std::string smi = "c1ccccc1O";
  std::string csmi = canon_smiles(smi);
  std::cerr << csmi << std::endl;
  assert(smi == "Oc1ccccc1");
}