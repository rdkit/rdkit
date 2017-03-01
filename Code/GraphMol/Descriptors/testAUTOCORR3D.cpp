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
  std::cout << "=>start test autocorr3D\n";

  std::string pathName = getenv("RDBASE");
  
  //std::cout << "Path: " << pathName << "\n";

  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);

  std::string fName = pathName+"/Code/GraphMol/Descriptors/test_data/auto3D_dragon.out";

  std::ifstream instrm(fName.c_str());

  std::string line;
  std::vector<std::vector<std::string>> data;

  while(std::getline(instrm, line)) {
      std::string phrase;
      std::vector<std::string> row;
      std::stringstream ss(line);
      while(std::getline(ss, phrase, '\t')) {
          row.push_back(std::move(phrase));
      }

      data.push_back(std::move(row));
  }

  //std::cout << "=>read file ok\n";


  int nDone = 0;
  while (!reader.atEnd()) {

    if (nDone > 15) {
      break;
    }
    RDKit::ROMol *m = reader.next();
    TEST_ASSERT(m);
    std::string nm;
    m->getProp("_Name", nm);

    std::vector<double> da3d;

    RDKit::Descriptors::AUTOCORR3D(*m, da3d, -1);

    std::vector<std::string> myrow=data[nDone];
    std::string inm= myrow[0];

    std::cout << inm.c_str() << ":";

    TEST_ASSERT(inm==nm);


    for (int j = 0; j < 14; j++) {
      std::cout << da3d[j] << ",";
    }
    std::cout << "\n";


    /*for (int i = 0; i < 80 ; i++) {
          double ref =atof(myrow[i+1].c_str());


          if(fabs(ref-da3d[i])>0.05){
            std::cout<<"value mismatch: pos" << i <<" "<< inm <<" "<< ref <<" "<< da3d[i] << std::endl;
          }

           //TEST_ASSERT(fabs(ref-drdf[i])<0.05);
    }*/

    // FIX: at the moment this isn't actually testing anything, it's just
    // running the calculation. The best test would be an SDF that has some
    // molecules with 3D structures and calculated values of the individual
    // descriptors (from DRAGON for example) that you can compare against.
    // many of the tests in the test.cpp directory here (for example
    // testLipinski1()) show how to do this.
    //for (int j = 0; j < 80; j++) {
    //  std::cout << dwhim[j] << ",";
    //}
    //std::cout << "\n";

    //}

    //std::cout << "done \n" << nDone << std::endl;

    delete m;
     ++nDone;

  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

int main(int argc, char *argv[]) {
  RDLog::InitLogs();
  testautocorrelation();
}
