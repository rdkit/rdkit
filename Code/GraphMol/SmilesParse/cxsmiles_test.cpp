//
//  Copyright (C) 2016 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <string>
#include <GraphMol/RDKitBase.h>
#include "SmilesParse.h"
#include "SmilesWrite.h"
#include <RDGeneral/RDLog.h>
using namespace RDKit;

void testBase() {
  BOOST_LOG(rdInfoLog) << "testing base functionality"
    << std::endl;
  { // it works when nothing is provided
    std::string smiles = "CC";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}
void testCoords2D() {
  BOOST_LOG(rdInfoLog) << "testing reading 2D coordinates"
                       << std::endl;
  {
    std::string smiles = "CC |(0,.75,;0,-.75,)|";
    SmilesParserParams params;
    params.allowCXSMILES = true; 
    ROMol *m = SmilesToMol(smiles,params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==2);
    TEST_ASSERT(m->getNumConformers()==1);

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();
  testBase();
  testCoords2D();
}
