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
#include <GraphMol/ReduceGraphs/ReducedGraphs.h>
#include <DataStructs/ExplicitBitVect.h>

#include <RDGeneral/RDLog.h>
#include <string>

using namespace RDKit;

void test1(){
  BOOST_LOG(rdInfoLog) <<"testing basics" << std::endl;
  {
    std::string smi = "C1=CC=CC=C1";
    RWMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1->getNumAtoms()==6);

    delete m1;
  }

  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}


int main(int argc,char *argv[]){
  RDLog::InitLogs();
  test1();
  return 0;
}
