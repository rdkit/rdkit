// $Id$
//
//  Copyright (C) 2007 Greg Landrum
//       All Rights Reserved
//

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include "FileParsers.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <string>

using namespace RDKit;

void test1(){
  BOOST_LOG(rdInfoLog) << "testing basic tpl parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/cmpd2.tpl";
  RWMol *m = TPLFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==12);
  TEST_ASSERT(m->getNumBonds()==12);
  TEST_ASSERT(m->getNumConformers()==2);
  
  delete m;
  m = TPLFileToMol(fName,true,true);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==12);
  TEST_ASSERT(m->getNumBonds()==12);
  TEST_ASSERT(m->getNumConformers()==1);

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


int main(int argc,char *argv[]){
  RDLog::InitLogs();
  test1();
  return 0;
}
