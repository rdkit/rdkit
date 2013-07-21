// $Id$
//
//  Copyright (C) 2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ReducedGraphs/ReducedGraphs.h>
#include <DataStructs/ExplicitBitVect.h>

#include <RDGeneral/RDLog.h>
#include <string>

using namespace RDKit;

void test1(){
  BOOST_LOG(rdInfoLog) <<"testing basics" << std::endl;
  {
    std::string smi = "c1ccccc1CCO";
    ROMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    
    RDNumeric::DoubleVector *fp=ReducedGraphs::getErGFingerprint(*m1);
    std::cerr<<*fp<<std::endl;
    
    delete m1;
  }

  {
    std::string smi = "c1cnccc1CCO";
    ROMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    
    RDNumeric::DoubleVector *fp=ReducedGraphs::getErGFingerprint(*m1);
    std::cerr<<*fp<<std::endl;
    
    delete m1;
  }

  {
    std::string smi = "OCCC1=CC2=C(C=CC=C2)C=C1";
    ROMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    
    RDNumeric::DoubleVector *fp=ReducedGraphs::getErGFingerprint(*m1);
    std::cerr<<*fp<<std::endl;
    
    delete m1;
  }

  {
    std::string smi = "OCCC1=CC2=C(C=CC=C2)N=C1";
    ROMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    
    RDNumeric::DoubleVector *fp=ReducedGraphs::getErGFingerprint(*m1);
    std::cerr<<*fp<<std::endl;
    
    delete m1;
  }
  {
    std::string smi = "OCCC1=CC2=C(CCCC2)N=C1";
    ROMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    
    RDNumeric::DoubleVector *fp=ReducedGraphs::getErGFingerprint(*m1);
    std::cerr<<*fp<<std::endl;
    
    delete m1;
  }
  
  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}


int main(int argc,char *argv[]){
  RDLog::InitLogs();
  test1();
  return 0;
}
