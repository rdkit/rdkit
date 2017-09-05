//
//  Copyright (C) 2017 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include "FileParsers.h"
#include "MolFileStereochem.h"
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <clocale>
#include <cstdlib>

#include <string>
#include <fstream>

using namespace RDKit;

void test1() {
  BOOST_LOG(rdInfoLog) << "testing basic extended stereo parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/two_centers_or.mol";
  RWMol *m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 8);
  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
#if 1
  test1();
#endif

  return 0;
}
