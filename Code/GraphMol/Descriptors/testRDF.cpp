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

#include <GraphMol/Descriptors/RDF.h>

void testRDF(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Basic RDF tests." << std::endl;

  std::string pathName = getenv("RDBASE");
  std::string sdfName = pathName+"/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";
  RDKit::SDMolSupplier reader(sdfName,true,false);
  std::string fName = pathName+"/Code/GraphMol/Descriptors/test_data/PBF_egfr.out";
  std::ifstream instrm(fName.c_str());
  int nDone=0;
  while(!reader.atEnd()){
    RDKit::ROMol *m=reader.next();
    TEST_ASSERT(m);
    std::string nm;
    m->getProp("_Name",nm);
    //std::vector<double> drbf=RDKit::Descriptors::RBF(*m);
    std::vector<double> dpbf=RDKit::Descriptors::RDF(*m);

    delete m;
    ++nDone;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


int main(int argc, char *argv[])
{
  RDLog::InitLogs();
  testRDF();
}
