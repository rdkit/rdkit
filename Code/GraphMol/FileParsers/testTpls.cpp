// $Id$
//
//  Copyright (C) 2007 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
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
  RWMol *m;
  Conformer conf;
  std::string propVal;

  m = TPLFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==12);
  TEST_ASSERT(m->getNumBonds()==12);
  TEST_ASSERT(m->getNumConformers()==2);
  TEST_ASSERT(m->hasProp("_Name"));
  m->getProp("_Name",propVal);
  TEST_ASSERT(propVal=="compound 2");
  TEST_ASSERT(m->hasProp("Conf_1_Name"));
  m->getProp("Conf_1_Name",propVal);
  TEST_ASSERT(propVal=="conf 1");
  
  conf = m->getConformer(0);
  TEST_ASSERT(feq(conf.getAtomPos(0).x,-1.02))
  TEST_ASSERT(feq(conf.getAtomPos(0).y,0.96))
  TEST_ASSERT(feq(conf.getAtomPos(0).z,-0.04))
  
  conf = m->getConformer(1);
  TEST_ASSERT(feq(conf.getAtomPos(0).x,-2.02))
  TEST_ASSERT(feq(conf.getAtomPos(0).y,0.96))
  TEST_ASSERT(feq(conf.getAtomPos(0).z,-0.04))

    
  delete m;
  m = TPLFileToMol(fName,true,true);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==12);
  TEST_ASSERT(m->getNumBonds()==12);
  TEST_ASSERT(m->getNumConformers()==1);
  conf = m->getConformer(0);
  TEST_ASSERT(feq(conf.getAtomPos(0).x,-1.02))
  TEST_ASSERT(feq(conf.getAtomPos(0).y,0.96))
  TEST_ASSERT(feq(conf.getAtomPos(0).z,-0.04))

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test2(){
  BOOST_LOG(rdInfoLog) << "testing tpl writing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/cmpd1.tpl";
  RWMol *m,*m2;
  Conformer conf;
  std::string propVal;

  m = TPLFileToMol(fName,true,true);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==18);
  TEST_ASSERT(m->getNumBonds()==18);
  TEST_ASSERT(m->getNumConformers()==9);
  TEST_ASSERT(m->hasProp("_Name"));
  m->getProp("_Name",propVal);
  TEST_ASSERT(propVal=="compound 1");
  TEST_ASSERT(m->hasProp("Conf_1_Name"));
  m->getProp("Conf_1_Name",propVal);
  TEST_ASSERT(propVal=="conf9");

  std::stringstream strm;
  strm << MolToTPLText(*m,"TPLCharge");

  unsigned int line;
  m2 = TPLDataStreamToMol(&strm,line);
  TEST_ASSERT(m2);
  TEST_ASSERT(m2->getNumAtoms()==18);
  TEST_ASSERT(m2->getNumBonds()==18);
  TEST_ASSERT(m2->getNumConformers()==9);
  TEST_ASSERT(m2->hasProp("_Name"));
  m2->getProp("_Name",propVal);
  TEST_ASSERT(propVal=="compound 1");
  TEST_ASSERT(m2->hasProp("Conf_1_Name"));
  m2->getProp("Conf_1_Name",propVal);
  TEST_ASSERT(propVal=="conformer_1");

  
  std::stringstream strm2;
  strm2 << MolToTPLText(*m,"TPLCharge",true);

  delete m2;
  m2 = TPLDataStreamToMol(&strm2,line,true,true);
  TEST_ASSERT(m2);
  TEST_ASSERT(m2->getNumAtoms()==18);
  TEST_ASSERT(m2->getNumBonds()==18);
  TEST_ASSERT(m2->getNumConformers()==8);
  TEST_ASSERT(m2->hasProp("_Name"));
  m2->getProp("_Name",propVal);
  TEST_ASSERT(propVal=="compound 1");
  TEST_ASSERT(m2->hasProp("Conf_1_Name"));
  m2->getProp("Conf_1_Name",propVal);
  TEST_ASSERT(propVal=="conformer_1");
  
  delete m;
  delete m2;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


int main(int argc,char *argv[]){
  RDLog::InitLogs();
  test1();
  test2();
  return 0;
}
