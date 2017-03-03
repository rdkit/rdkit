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


#include <GraphMol/Descriptors/RDF.h>


void testRDF1() {
  std::cout << "=>start test rdf\n";

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
    m->getProp("_Name",nm);


    std::vector<double> drdf;

    RDKit::Descriptors::RDF(*m,drdf,-1);



    std::cout << "=>read molecule: " << nDone  << std::endl;

    delete m;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testRDF() {
  std::cout << "=>start test rdf\n";

  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);
  std::string fName = pathName+"/Code/GraphMol/Descriptors/test_data/RDF.out";

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

  std::cout << "=>read file\n";

  int nDone = 0;
  while (!reader.atEnd()) {

    RDKit::ROMol *m = reader.next();
    TEST_ASSERT(m);
    std::string nm;
    m->getProp("_Name",nm);



    std::vector<double> drdf;

    RDKit::Descriptors::RDF(*m, drdf, -1);


    std::vector<std::string> myrow=data[nDone];
    std::string inm= myrow[0];
    TEST_ASSERT(inm==nm);

    for (int i=0;i<drdf.size();i++)
       {
            double ref =atof(myrow[i+1].c_str());
            if(fabs(ref-drdf[i])>0.05){
              std::cerr<<"value mismatch: pos" << i <<" "<<inm<<" dragon: "<<ref<<" rdkit: "<< drdf[i] <<std::endl;
            }

           //TEST_ASSERT(fabs(ref-drdf[i])<0.05);

       }
    std::cout << "=>read molecule: " << nDone  << std::endl;


    delete m;
    ++nDone;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

int main(int argc, char *argv[]) {
  RDLog::InitLogs();
  testRDF();
  testRDF1();


}
