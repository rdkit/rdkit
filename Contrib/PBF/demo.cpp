// $Id$
//
//  Copyright (C) 2012 Greg Landrum
//   @@ All Rights Reserved @@
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <fstream>

#include "PBFRDKit.h"

void calc() {
  std::string fname = "testData/egfr.sdf";
  RDKit::SDMolSupplier reader(fname, true, false);

  while (!reader.atEnd()) {
    RDKit::ROMol *m = reader.next();
    if (!m) continue;
    std::string nm;
    m->getProp("_Name", nm);
    double dpbf = PBFRD(*m);
    std::cout << nm << " " << dpbf << std::endl;
    delete m;
  }
}

void test() {
  std::string fname = "testData/egfr.sdf";
  RDKit::SDMolSupplier reader(fname, true, false);
  std::ifstream instrm("testData/egfr.out");
  int nDone = 0;
  while (!reader.atEnd()) {
    RDKit::ROMol *m = reader.next();
    if (!m) continue;
    std::string nm;
    m->getProp("_Name", nm);
    double dpbf = PBFRD(*m);

    std::string inm;
    double ref;
    instrm >> inm;
    instrm >> ref;
    if (inm != nm) {
      std::cerr << "name mismatch: " << inm << " " << nm << std::endl;
      continue;
    }
    if (fabs(ref - dpbf) > .001) {
      std::cerr << "value mismatch: " << inm << " " << ref << " " << dpbf
                << std::endl;
    }
    delete m;
    ++nDone;
  }
}

int main(int argc, char *argv[]) {
  RDLog::InitLogs();
  // calc();
  test();
  // for(unsigned int i=0;i<100;++i){
  // calc();
  // }
}
