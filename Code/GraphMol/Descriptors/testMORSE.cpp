//  Created by Guillaume GODIN
//  Copyright (C) 2012-2016 Greg Landrum
//   @@ All Rights Reserved @@
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <fstream>

#include <GraphMol/Descriptors/MORSE.h>

void testMORSE() {
  std::cout << "=>start test MORSE\n";

  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);

  std::string fName =
      pathName + "/Code/GraphMol/Descriptors/test_data/MORSE.out";

  std::ifstream instrm(fName.c_str());

  std::string line;
  std::vector<std::vector<std::string>> data;

  while (std::getline(instrm, line)) {
    std::string phrase;
    std::vector<std::string> row;
    std::stringstream ss(line);
    while (std::getline(ss, phrase, '\t')) {
      row.push_back(phrase);
    }

    data.push_back(row);
  }

  std::cout << "=>read file\n";

  int nDone = 0;
  while (!reader.atEnd()) {
    RDKit::ROMol *m = reader.next();
    TEST_ASSERT(m);
    std::string nm;
    m->getProp("_Name", nm);

    std::vector<double> dmorse;

    RDKit::Descriptors::MORSE(*m, dmorse, -1);

    std::vector<std::string> myrow = data[nDone];
    std::string inm = myrow[0];
    TEST_ASSERT(inm == nm);

    for (size_t i = 0; i < dmorse.size(); i++) {
      double ref = atof(myrow[i + 1].c_str());

      if (fabs(ref) > 0.01) {
        if (fabs((ref - dmorse[i]) / ref) > 1) {
          std::cout << "value mismatch: pos" << i << " " << inm << " " << ref
                    << " " << dmorse[i] << std::endl;
        }
      }

      if (fabs(ref) < 0.01) {
        if (fabs(ref - dmorse[i]) > 0.02) {
          std::cout << "value mismatch: pos" << i << " " << inm << " " << ref
                    << " " << dmorse[i] << std::endl;
        }
      }
      if (ref > 1 && fabs(ref - dmorse[i]) / ref > 0.02) {
        std::cout << "value mismatch: pos" << i << " " << inm << " " << ref
                  << " " << dmorse[i] << std::endl;
      }
      // we're testing reasonably sized values and want to be sure that we're
      // within 2% of the reference.
      TEST_ASSERT(ref < 1 || fabs(ref - dmorse[i]) / ref < 0.02);
    }

    delete m;
    ++nDone;
  }

  BOOST_LOG(rdErrorLog) << "test on : " << nDone << " molecules done"
                        << std::endl;
}

int main() {
  RDLog::InitLogs();
  testMORSE();
}
