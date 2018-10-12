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
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>

#include <GraphMol/Descriptors/WHIM.h>

void testWHIM2() {
  std::cout << "=>start test chlorobenzene whim from rdkit\n";

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

    RDKit::Descriptors::WHIM(*m, dwhim, -1, 0.01);
    for (int j = 0; j < 114; j++) {
      std::cout << dwhim[j] << ",";
    }
    std::cout << "\n";

    delete m;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testWHIM3() {
  std::cout << "=>start test chlorobenzene whim original\n";

  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/chlorobenzene2.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);

  int nDone = 0;
  while (!reader.atEnd()) {
    ++nDone;

    RDKit::ROMol *m = reader.next();
    TEST_ASSERT(m);
    std::string nm;
    m->getProp("_Name", nm);

    std::vector<double> dwhim;

    RDKit::Descriptors::WHIM(*m, dwhim, -1, 0.01);
    for (int j = 0; j < 114; j++) {
      std::cout << dwhim[j] << ",";
    }
    std::cout << "\n";

    delete m;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testWHIM1() {
  std::cout << "=>start test rdf\n";

  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/1mol.sdf";

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
    RDKit::Descriptors::WHIM(*m, dwhim, -1, 0.01);
    for (int j = 0; j < 114; j++) {
      std::cout << dwhim[j] << ",";
    }
    std::cout << "\n";

    //}

    std::cout << "=>read molecule: " << nDone << std::endl;

    delete m;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testWHIM() {
  std::cout << "=>start test WHIM\n";
  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";
  RDKit::SDMolSupplier reader(sdfName, true, false);
  std::string fName =
      pathName + "/Code/GraphMol/Descriptors/test_data/whim.new.out";
  std::ifstream instrm(fName.c_str());

  // std::string ofName =
  //     pathName + "/Code/GraphMol/Descriptors/test_data/whim.new.out";
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
    RDKit::ROMol *m = reader.next();
    TEST_ASSERT(m);
    std::string nm;
    m->getProp("_Name", nm);
    std::vector<double> dwhim;
    RDKit::Descriptors::WHIM(*m, dwhim, -1, 0.01);
    std::vector<std::string> myrow = data[nDone];
    std::string inm = myrow[0];
    TEST_ASSERT(inm == nm);
    // outstrm << nm;
    for (int i = 0; i < 114; i++) {
      double ref = atof(myrow[i + 1].c_str());
      if (fabs(ref - dwhim[i]) >= 0.01) {
        std::cerr << "value mismatch: pos" << i << " " << inm << " " << ref
                  << " " << dwhim[i] << std::endl;
      }
      // outstrm << "\t" << dwhim[i];
      // TEST_ASSERT((ref < 1e-3 && dwhim[i] < 1e-3) ||
      //             fabs(ref - dwhim[i]) / ref < 0.01);
      TEST_ASSERT(fabs(ref - dwhim[i]) < 0.01);
    }
    // outstrm << std::endl;
    delete m;
    ++nDone;
  }

  BOOST_LOG(rdErrorLog) << "test on : " << nDone << " molecules done"
                        << std::endl;
}

int main() {
  RDLog::InitLogs();
  testWHIM();
  // testWHIM3();

  //  testWHIM();
}
