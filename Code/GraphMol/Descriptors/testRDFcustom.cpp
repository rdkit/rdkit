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
#include "GraphMol/PartialCharges/GasteigerCharges.h"
#include "GraphMol/PartialCharges/GasteigerParams.h"
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>

#include <GraphMol/Descriptors/RDF.h>
#include <GraphMol/Descriptors/MORSE.h>
#include <GraphMol/Descriptors/AUTOCORR3D.h>
#include <GraphMol/Descriptors/WHIM.h>  // strange
#include <GraphMol/Descriptors/GETAWAY.h>
// make a 2D test of not ?
#include <GraphMol/Descriptors/AUTOCORR2D.h>

void testRDFcustom() {
  std::cout << "=>start test rdf custom\n";

  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);
  std::string fName =
      pathName + "/Code/GraphMol/Descriptors/test_data/RDFcustom.out";

  std::ifstream instrm(fName.c_str());

  std::string line;
  std::vector<std::vector<std::string>> data;

  while (std::getline(instrm, line)) {
    std::string phrase;
    std::vector<std::string> row;
    std::stringstream ss(line);
    while (std::getline(ss, phrase, ',')) {
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

    std::vector<double> drdf;
    std::vector<double> charges(m->getNumAtoms(), 0);
    RDKit::computeGasteigerCharges(*m, charges, 12, true);

    const std::string atomprop = "_GasteigerCharge";
    RDKit::Descriptors::RDF(*m, drdf, -1, atomprop);

    std::vector<std::string> myrow = data[nDone];
    std::string inm = myrow[0];
    TEST_ASSERT(inm == nm);

    // std::cerr << inm << ",";
    for (size_t i = 0; i < drdf.size(); i++) {
      double ref = atof(myrow[i + 1].c_str());
      // std::cerr << "[" << i << ":" << drdf[i] << "," << ref <<"]|";
      // if (fabs(ref - drdf[i]) > 0.01) {
      //  std::cerr << "pos" << i << " " << inm
      //    << " file: " << ref << " computed: " << drdf[i] ; // << std::endl;
      //}

      TEST_ASSERT(fabs(ref - drdf[i]) < 0.001);
    }
    // std::cerr << "\n";

    delete m;
    ++nDone;
  }

  BOOST_LOG(rdErrorLog) << "test on : " << nDone << " molecules done"
                        << std::endl;
}

void testMORSEcustom() {
  std::cout << "=>start test morse custom\n";

  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);
  std::string fName =
      pathName + "/Code/GraphMol/Descriptors/test_data/MORSEcustom.out";

  std::ifstream instrm(fName.c_str());

  std::string line;
  std::vector<std::vector<std::string>> data;

  while (std::getline(instrm, line)) {
    std::string phrase;
    std::vector<std::string> row;
    std::stringstream ss(line);
    while (std::getline(ss, phrase, ',')) {
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

    std::vector<double> dmorse;
    std::vector<double> charges(m->getNumAtoms(), 0);
    RDKit::computeGasteigerCharges(*m, charges, 12, true);

    const std::string atomprop = "_GasteigerCharge";
    RDKit::Descriptors::MORSE(*m, dmorse, -1, atomprop);

    std::vector<std::string> myrow = data[nDone];
    std::string inm = myrow[0];
    TEST_ASSERT(inm == nm);

    // std::cerr << inm << ",";
    for (size_t i = 0; i < dmorse.size(); i++) {
      double ref = atof(myrow[i + 1].c_str());
      // std::cerr << dmorse[i] << ",";
      // std::cerr << "[" << i << ":" << drdf[i] << "," << ref <<"]|";
      // if (fabs(ref - drdf[i]) > 0.01) {
      //  std::cerr << "pos" << i << " " << inm
      //    << " file: " << ref << " computed: " << drdf[i] ; // << std::endl;
      //}

      TEST_ASSERT(fabs(ref - dmorse[i]) < 0.001);
    }
    // std::cerr << "\n";

    delete m;
    ++nDone;
  }

  BOOST_LOG(rdErrorLog) << "test on : " << nDone << " molecules done"
                        << std::endl;
}

void testAUTOCORR3Dcustom() {
  std::cout << "=>start test AUTOCORR3D custom\n";

  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);
  std::string fName =
      pathName + "/Code/GraphMol/Descriptors/test_data/auto3Dcustom.out";

  std::ifstream instrm(fName.c_str());

  std::string line;
  std::vector<std::vector<std::string>> data;

  while (std::getline(instrm, line)) {
    std::string phrase;
    std::vector<std::string> row;
    std::stringstream ss(line);
    while (std::getline(ss, phrase, ',')) {
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

    std::vector<double> dauto3d;
    std::vector<double> charges(m->getNumAtoms(), 0);
    RDKit::computeGasteigerCharges(*m, charges, 12, true);

    const std::string atomprop = "_GasteigerCharge";
    RDKit::Descriptors::AUTOCORR3D(*m, dauto3d, -1, atomprop);

    std::vector<std::string> myrow = data[nDone];
    std::string inm = myrow[0];
    TEST_ASSERT(inm == nm);

    // std::cerr << inm << ",";
    for (size_t i = 0; i < dauto3d.size(); i++) {
      double ref = atof(myrow[i + 1].c_str());
      // std::cerr << dauto3d[i] << ",";
      // std::cerr << "[" << i << ":" << drdf[i] << "," << ref <<"]|";
      // if (fabs(ref - drdf[i]) > 0.01) {
      //  std::cerr << "pos" << i << " " << inm
      //    << " file: " << ref << " computed: " << drdf[i] ; // << std::endl;
      //}

      TEST_ASSERT(fabs(ref - dauto3d[i]) < 0.001);
    }
    // std::cerr << "\n";

    delete m;
    ++nDone;
  }

  BOOST_LOG(rdErrorLog) << "test on : " << nDone << " molecules done"
                        << std::endl;
}

void testWHIMcustom() {
  std::cout << "=>start test 1 WHIM custom\n";

  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);
  std::string fName =
      pathName + "/Code/GraphMol/Descriptors/test_data/WHIMcustom.out";

  std::ifstream instrm(fName.c_str());

  std::string line;
  std::vector<std::vector<std::string>> data;

  while (std::getline(instrm, line)) {
    std::string phrase;
    std::vector<std::string> row;
    std::stringstream ss(line);
    while (std::getline(ss, phrase, ',')) {
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
    std::vector<double> charges(m->getNumAtoms(), 0);
    RDKit::computeGasteigerCharges(*m, charges, 12, true);

    const std::string atomprop = "_GasteigerCharge";
    RDKit::Descriptors::WHIM(*m, dwhim, -1, 0.01, atomprop);
    /*
        for (unsigned int aix = 0; aix < m->getNumAtoms(); aix++) {
          std::cerr <<
       m->getAtomWithIdx(aix)->getProp<double>("_GasteigerCharge") << ",";
        }
        std::cerr << "\n";
    */
    std::vector<std::string> myrow = data[nDone];
    std::string inm = myrow[0];
    TEST_ASSERT(inm == nm);

    // std::cerr << inm << ",";
    for (size_t i = 0; i < dwhim.size(); i++) {
      double ref = atof(myrow[i + 1].c_str());
      // std::cerr << dwhim[i] << ",";
      // std::cerr << "[" << i << ":" << drdf[i] << "," << ref <<"]|";
      // if (fabs(ref - drdf[i]) > 0.01) {
      //  std::cerr << "pos" << i << " " << inm
      //    << " file: " << ref << " computed: " << drdf[i] ; // << std::endl;
      //}

      TEST_ASSERT(fabs(ref - dwhim[i]) < 0.01);
    }
    // std::cerr << "\n";

    delete m;
    ++nDone;
  }

  BOOST_LOG(rdErrorLog) << "test on : " << nDone << " molecules done"
                        << std::endl;
}

void testWHIMcustom1() {
  std::cout << "=>start test 2 WHIM custom\n";

  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);
  std::string fName =
      pathName + "/Code/GraphMol/Descriptors/test_data/WHIM1custom.out";

  std::ifstream instrm(fName.c_str());

  std::string line;
  std::vector<std::vector<std::string>> data;

  while (std::getline(instrm, line)) {
    std::string phrase;
    std::vector<std::string> row;
    std::stringstream ss(line);
    while (std::getline(ss, phrase, ',')) {
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

    for (unsigned int aix = 0; aix < m->getNumAtoms(); aix++) {
      m->getAtomWithIdx(aix)->setProp("_MyCharge", 0.001, true);
    }

    const std::string atomprop = "_MyCharge";
    RDKit::Descriptors::WHIM(*m, dwhim, -1, 0.01, atomprop);

    std::vector<std::string> myrow = data[nDone];
    std::string inm = myrow[0];
    TEST_ASSERT(inm == nm);

    // std::cerr << inm << ",";
    for (size_t i = 0; i < dwhim.size(); i++) {
      double ref = atof(myrow[i + 1].c_str());
      // std::cerr << dwhim[i] << ",";
      // std::cerr << "[" << i << ":" << drdf[i] << "," << ref <<"]|";
      // if (fabs(ref - drdf[i]) > 0.01) {
      //  std::cerr << "pos" << i << " " << inm
      //    << " file: " << ref << " computed: " << dwhim[i] ; // << std::endl;
      //}

      TEST_ASSERT(fabs(ref - dwhim[i]) < 0.01);
    }
    // std::cerr << "\n";

    delete m;
    ++nDone;
  }

  BOOST_LOG(rdErrorLog) << "test on : " << nDone << " molecules done"
                        << std::endl;
}

void testGETAWAYcustom() {
  std::cout << "=>start test GETAWAY custom\n";

  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);
  std::string fName =
      pathName + "/Code/GraphMol/Descriptors/test_data/GETAWAYcustom.out";

  std::ifstream instrm(fName.c_str());

  std::string line;
  std::vector<std::vector<std::string>> data;

  while (std::getline(instrm, line)) {
    std::string phrase;
    std::vector<std::string> row;
    std::stringstream ss(line);
    while (std::getline(ss, phrase, ',')) {
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

    std::vector<double> dgetaway;
    std::vector<double> charges(m->getNumAtoms(), 0);
    RDKit::computeGasteigerCharges(*m, charges, 12, true);

    const std::string atomprop = "_GasteigerCharge";

    /*
     * TO DO: fix this call: precision parameter is an *int*,
     *        but setting it to 4 digits (= 0.0001) breaks this
     *        test due to the references.
     */
    RDKit::Descriptors::GETAWAY(*m, dgetaway, -1, 0U, atomprop);

    std::vector<std::string> myrow = data[nDone];
    std::string inm = myrow[0];
    TEST_ASSERT(inm == nm);

    // std::cerr << inm << ",";
    for (size_t i = 0; i < dgetaway.size(); i++) {
      double ref = atof(myrow[i + 1].c_str());
      // std::cerr << dgetaway[i] << ",";
      // std::cerr << "[" << i << ":" << drdf[i] << "," << ref <<"]|";
      // if (fabs(ref - drdf[i]) > 0.01) {
      //  std::cerr << "pos" << i << " " << inm
      //    << " file: " << ref << " computed: " << drdf[i] ; // << std::endl;
      //}

      TEST_ASSERT(fabs(ref - dgetaway[i]) < 0.001);
    }
    // std::cerr << "\n";

    delete m;
    ++nDone;
  }

  BOOST_LOG(rdErrorLog) << "test on : " << nDone << " molecules done"
                        << std::endl;
}

int main() {
  RDLog::InitLogs();
  testRDFcustom();
  testMORSEcustom();
  testAUTOCORR3Dcustom();
  testWHIMcustom1();
  testWHIMcustom();  // case of divide by zeros (ie when the sum of weight are
                     // near zero!)
  testGETAWAYcustom();
}
