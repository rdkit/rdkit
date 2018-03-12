//  Created by Guillaume GODIN
//  "Copyright 2013-2016 Tomas Racek (tom@krab1k.net)"
//  Copyright (C) 2012-2017 Greg Landrum
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
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>

#include <GraphMol/Descriptors/EEM.h>

void testEEM() {
  std::cout << "=>start test EEM\n";
  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/set00.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);

  std::string fName =
      pathName + "/Code/GraphMol/Descriptors/test_data/eem.out";

  std::ifstream instrm(fName.c_str());

  std::string line;
  std::vector<std::vector<std::string> > data;

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
  int errorMols = 0;
  while (!reader.atEnd()) {
    RDKit::ROMol *m = reader.next();
    TEST_ASSERT(m);
    std::string nm;
    m->getProp("_Name", nm);
    int errorAtoms = 0;

    std::vector<std::string> myrow = data[nDone];
    std::string inm = myrow[0];

    TEST_ASSERT(inm == nm);

    int confId=-1;
    std::vector<double> charges;

    RDKit::Descriptors::EEM(*m, charges, confId);
    int numAtoms = m->getNumAtoms();

    for (int i = 0; i < numAtoms; i++) {
      double ref = atof(myrow[i + 2].c_str());
      if(fabs(ref - charges[i]) >= 0.01) {

         //std::cout << inm << "," << "ref: " << ref << " ,val: "<<  charges[i] << "\n";
             ++errorAtoms;
           }

      //TEST_ASSERT(fabs(ref - charges[i]) < 0.01);
    }

    if(errorAtoms>0) {
      std::cout << inm << "\n";
      ++errorMols;
    }

    delete m;
    //break;
    ++nDone;
  }
    std::cout << "Errors:" <<  errorMols << "\n";


  BOOST_LOG(rdErrorLog) << "test on : " << nDone << " molecules done"
                        << std::endl;
}

int main(int argc, char *argv[]) {
  RDLog::InitLogs();
  testEEM();
}