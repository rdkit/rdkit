//
//  Copyright (C) 2012-2016 Greg Landrum
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
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>

#include <GraphMol/Descriptors/AUTOCORR3D.h>

void testautocorrelation() {
  std::cout << "=>start test rdf\n";

  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/chlorobenzene.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);

  int nDone = 0;
  while (!reader.atEnd()) {
    ++nDone;

    RDKit::ROMol *m = reader.next();
    TEST_ASSERT(m);
    std::string nm;
    m->getProp("_Name", nm);

    std::vector<double> dwhim;
    // for (int i=1;i<11;i++) {
    // std::cout << "i:" << 0.005*i << "\n";
    RDKit::Descriptors::AUTOCORR3D(*m, dwhim, -1);

    // FIX: at the moment this isn't actually testing anything, it's just
    // running the calculation. The best test would be an SDF that has some
    // molecules with 3D structures and calculated values of the individual
    // descriptors (from DRAGON for example) that you can compare against.
    // many of the tests in the test.cpp directory here (for example
    // testLipinski1()) show how to do this.
    for (int j = 0; j < 80; j++) {
      std::cout << dwhim[j] << ",";
    }
    std::cout << "\n";

    //}

    std::cout << "=>read molecule: " << nDone << std::endl;

    delete m;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

int main(int argc, char *argv[]) {
  RDLog::InitLogs();
  testautocorrelation();
}
