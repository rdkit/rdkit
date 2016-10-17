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


#include <GraphMol/Descriptors/WHIM.h>

void testRDF() {
  std::cout << "=>start test rdf\n";

  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);
 
  int nDone = 0;
  while (!reader.atEnd()) {
    ++nDone;

    RDKit::ROMol *m = reader.next();
    TEST_ASSERT(m);
    std::string nm;
    m->getProp("_Name",nm);


    double drdf = RDKit::Descriptors::WHIM(*m);


  
       
    std::cout << "=>read molecule: " << nDone  << std::endl;

    delete m;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

int main(int argc, char *argv[]) {
  RDLog::InitLogs();
  testRDF();
}
