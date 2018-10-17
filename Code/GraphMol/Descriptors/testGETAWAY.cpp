//
//  Copyright (c) 2016, Guillaume GODIN
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
#include <sstream>

#include <GraphMol/Descriptors/GETAWAY.h>

void testGETAWAY() {
  std::cout << "=>start test GETAWAY\n";

  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);
  std::string fName =
      pathName + "/Code/GraphMol/Descriptors/test_data/GETAWAY.new.out";

  std::ifstream instrm(fName.c_str());

  // std::string ofName =
  //     pathName + "/Code/GraphMol/Descriptors/test_data/GETAWAY.new.out";
  // std::ofstream outstrm(ofName.c_str());

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

  int nDone = 0;
  while (!reader.atEnd()) {
    // if (nDone > 10) {
    //   break;
    // }
    RDKit::ROMol *m = reader.next();
    TEST_ASSERT(m);
    std::string nm;
    m->getProp("_Name", nm);

    std::vector<double> dgetaway;

    RDKit::Descriptors::GETAWAY(*m, dgetaway);

    std::vector<std::string> myrow = data[nDone];
    std::string inm = myrow[0];
    TEST_ASSERT(inm == nm);
    // std::cout <<  "\n";
    // int numAtoms = m->getNumAtoms();
    // std::cout << "number of Atoms : " << numAtoms << "\n";

    // outstrm << nm << "\t";
    for (int i = 0; i < 273; i++) {
      double ref = atof(myrow[i + 1].c_str());

      if (fabs(ref) > 1) {
        if (fabs((ref - dgetaway[i]) / ref) > 0.01) {
          std::cerr << "value mismatch: pos" << i << " " << inm
                    << " dragon: " << ref << " rdkit: " << dgetaway[i]
                    << std::endl;
        }
      }
      if (fabs(ref) <= 1) {
        if (fabs(ref - dgetaway[i]) > 0.02) {
          std::cerr << "value mismatch: pos" << i << " " << inm
                    << " dragon: " << ref << " rdkit: " << dgetaway[i]
                    << std::endl;
        }
      }
      // if (i != 0) outstrm << "\t";
      // outstrm << dgetaway[i];
      TEST_ASSERT(fabs(ref - dgetaway[i]) < 0.05);
    }
    // outstrm << "\n";
    delete m;
    ++nDone;
    // if (nDone > 50) break;
  }

  BOOST_LOG(rdErrorLog) << "test on : " << nDone << " molecules done"
                        << std::endl;
}

int main() {
  RDLog::InitLogs();
  testGETAWAY();
}
