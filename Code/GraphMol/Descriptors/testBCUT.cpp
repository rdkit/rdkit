//
//  Copyright (C) 2020 Brian P. Kelley
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

#include <GraphMol/Descriptors/BCUT.h>

void test1(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Basic BCUT tests." << std::endl;

  std::string pathName = getenv("RDBASE");
  std::string sdfName = pathName+"/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";
  RDKit::SDMolSupplier reader(sdfName,true,false);
  //std::string fName = pathName+"/Code/GraphMol/Descriptors/test_data/PBF_egfr.out";
  //std::ifstream instrm(fName.c_str());
  int nDone=0;
  while(!reader.atEnd()){
    std::unique_ptr<RDKit::ROMol> m(reader.next());
    TEST_ASSERT(m);
    std::string nm;
    m->getProp("_Name",nm);
    std::vector<double> bcuts = RDKit::Descriptors::BCUT2D(*m);
    std::cerr << bcuts[0] << " " << bcuts[1] << " " << bcuts[2] << " " << bcuts[3] << std::endl;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

int main() {
  RDLog::InitLogs();
  test1();
}
