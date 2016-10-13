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

  std::cout << "=>start test rdf\n";

  std::string pathName = "/Users/mbp/Github/rdkit_mine";
  std::string sdfName = pathName+"/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";
  RDKit::SDMolSupplier reader(sdfName,true,false);
  std::cout << "=>read file\n";

  int nDone=0;
  while(!reader.atEnd()){
    ++nDone;

    RDKit::ROMol *m=reader.next();
    std::cout << "=>read molecule:\n"+nDone;

    
    std::vector<double> drdf= RDKit::Descriptors::RDF(*m);
   
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
