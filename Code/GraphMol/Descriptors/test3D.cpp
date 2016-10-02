//
//  Copyright (C) 2016 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <iostream>
#include <fstream>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/StreamOps.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <GraphMol/Descriptors/MolDescriptors3D.h>
using namespace RDKit;
using namespace RDKit::Descriptors;

bool compare(const std::string &inm,double ref,double val,double tol=1e-3){
  if(fabs(ref-val)>.001){
    std::cerr<<"value mismatch: "<<inm<<" "<<ref<<" "<<val<<std::endl;
  }
  return fabs(ref-val)<tol;
}

void testPMI1(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Basic PMI tests." << std::endl;

  std::string pathName = getenv("RDBASE");
  std::string sdfName = pathName+"/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";

  RDKit::SDMolSupplier reader(sdfName,true,false);
  std::string fName = pathName+"/Code/GraphMol/Descriptors/test_data/PMI_egfr.out";
  std::ifstream instrm(fName);
  int nDone=0;
  while(!reader.atEnd()){
    RDKit::ROMol *m=reader.next();
    TEST_ASSERT(m);
    RDKit::ROMol mcpy(*m);
    std::string nm;
    m->getProp("_Name",nm);
    std::string inm;
    instrm>>inm;
    TEST_ASSERT(inm==nm);
    double val;
    double pmi1_m,pmi2_m,pmi3_m,pmi1_nom,pmi2_nom,pmi3_nom;
    instrm>>pmi1_m;
    instrm>>pmi2_m;
    instrm>>pmi3_m;
    instrm>>pmi1_nom;
    instrm>>pmi2_nom;
    instrm>>pmi3_nom;

    val = RDKit::Descriptors::PMI1(*m);
    TEST_ASSERT(compare(inm,pmi1_m,val));
    val = RDKit::Descriptors::PMI2(*m);
    TEST_ASSERT(compare(inm,pmi2_m,val));
    val = RDKit::Descriptors::PMI3(*m);
    TEST_ASSERT(compare(inm,pmi3_m,val));

    val = RDKit::Descriptors::PMI1(*m,-1,false);
    TEST_ASSERT(compare(inm,pmi1_nom,val));
    val = RDKit::Descriptors::PMI2(*m,-1,false);
    TEST_ASSERT(compare(inm,pmi2_nom,val));
    val = RDKit::Descriptors::PMI3(*m,-1,false);
    TEST_ASSERT(compare(inm,pmi3_nom,val));

    // now try doing it in the reverse order to make sure caching doesn't
    // screw up.
    val = RDKit::Descriptors::PMI1(mcpy,-1,false);
    TEST_ASSERT(compare(inm,pmi1_nom,val));
    val = RDKit::Descriptors::PMI2(mcpy,-1,false);
    TEST_ASSERT(compare(inm,pmi2_nom,val));
    val = RDKit::Descriptors::PMI3(mcpy,-1,false);
    TEST_ASSERT(compare(inm,pmi3_nom,val));
    val = RDKit::Descriptors::PMI1(mcpy);
    TEST_ASSERT(compare(inm,pmi1_m,val));
    val = RDKit::Descriptors::PMI2(mcpy);
    TEST_ASSERT(compare(inm,pmi2_m,val));
    val = RDKit::Descriptors::PMI3(mcpy);
    TEST_ASSERT(compare(inm,pmi3_m,val));


    delete m;
    ++nDone;
}
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main() {
  RDLog::InitLogs();
  testPMI1();
}
