// $Id: testPickler.cpp 4998 2006-02-21 00:52:50Z glandrum $
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/RDLog.h>

#include <stdlib.h>
#include <ctime>
#include <iostream>
#include <fstream>


using namespace RDKit;

void test1(bool doLong=0){
  std::string fName = getenv("RDBASE");
  fName += "/Code/GraphMol/test_data/canonSmiles.smi";
  
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing String Pickle Roundtrips." << std::endl;

  SmilesMolSupplier suppl(fName,"\t",0,1,false);
  
  int count = 0;
  while(!suppl.atEnd()){
    ROMol *m = suppl.next();
    TEST_ASSERT(m);

    std::string pickle;
    MolPickler::pickleMol(*m,pickle);
    TEST_ASSERT(pickle.size());

    ROMol m2;
    MolPickler::molFromPickle(pickle,m2);
    std::string smi1=MolToSmiles(*m);
    //BOOST_LOG(rdInfoLog) << smi << " " << smi1 << std::endl;
    std::string smi2=MolToSmiles(m2);
    if(smi1!=smi2){
      BOOST_LOG(rdInfoLog) << "Line: " << count << "\n  " << smi1 << "\n != \n  " << smi2 << std::endl;
    }
    TEST_ASSERT(smi1==smi2);
    delete m;
    count++;
    if(!doLong && count >= 100) break;
  }  
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}


void _createPickleFile(){
  std::string smiName = getenv("RDBASE");
  smiName += "/Code/GraphMol/test_data/canonSmiles.smi";
  std::string pklName = getenv("RDBASE");
#ifdef OLD_PICKLE
  pklName += "/Code/GraphMol/test_data/canonSmiles.v1.pkl";
#else
  pklName += "/Code/GraphMol/test_data/canonSmiles.v2.pkl";
#endif  
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "creating pickle file." << std::endl;

  SmilesMolSupplier suppl(smiName,"\t",0,1,false);
  std::ofstream outStream(pklName.c_str(),std::ios_base::binary);
  int count = 0;
  while(!suppl.atEnd()){
    ROMol *m = suppl.next();
    TEST_ASSERT(m);

    std::string pickle;
    MolPickler::pickleMol(*m,outStream);
    delete m;
    count++;
  }  
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;

}

void test2(bool doLong=0){
  std::string smiName = getenv("RDBASE");
  smiName += "/Code/GraphMol/test_data/canonSmiles.smi";
  std::string pklName = getenv("RDBASE");
  pklName += "/Code/GraphMol/test_data/canonSmiles.v1.pkl";
  
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing reading existing pickle file (v1)." << std::endl;

  SmilesMolSupplier suppl(smiName,"\t",0,1,false);
  std::ifstream inStream(pklName.c_str(),std::ios_base::binary);
  int count = 0;
  while(!suppl.atEnd()){
    ROMol *m1 = suppl.next();
    TEST_ASSERT(m1);
    ROMol m2;
    MolPickler::molFromPickle(inStream,m2);
    
    std::string smi1=MolToSmiles(*m1);
    std::string smi2=MolToSmiles(m2);

    if(smi1!=smi2){
      BOOST_LOG(rdInfoLog) << "Line: " << count << "\n  " << smi1 << "\n != \n  " << smi2 << std::endl;
    }
    TEST_ASSERT(smi1==smi2);
    delete m1;
    count++;
    if(!doLong && count >= 100) break;
  }  
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}
  

void test3(bool doLong=0){
  std::string smiName = getenv("RDBASE");
  smiName += "/Code/GraphMol/test_data/canonSmiles.smi";
  std::string pklName = getenv("RDBASE");
  pklName += "/Code/GraphMol/test_data/canonSmiles.v2.pkl";
  
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing reading existing pickle file (v2)." << std::endl;

  SmilesMolSupplier suppl(smiName,"\t",0,1,false);
  std::ifstream inStream(pklName.c_str(),std::ios_base::binary);
  int count = 0;
  while(!suppl.atEnd()){
    ROMol *m1 = suppl.next();
    TEST_ASSERT(m1);
    ROMol m2;
    MolPickler::molFromPickle(inStream,m2);
    
    std::string smi1=MolToSmiles(*m1,1);
    std::string smi2=MolToSmiles(m2,1);

    if(smi1!=smi2){
      BOOST_LOG(rdInfoLog) << "Line: " << count << "\n  " << smi1 << "\n != \n  " << smi2 << std::endl;
    }
    TEST_ASSERT(smi1==smi2);
    delete m1;
    count++;
    if(!doLong && count >= 100) break;
  }  
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}
  
void timeTest(bool doLong=0){
  time_t t1,t2;
  std::string smiName = getenv("RDBASE");
  smiName += "/Code/GraphMol/test_data/canonSmiles.smi";

  std::string pklName = getenv("RDBASE");
  pklName += "/Code/GraphMol/test_data/canonSmiles.v2.pkl";
  
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Timing reads." << std::endl;

  t1 = std::time(0);
  SmilesMolSupplier suppl(smiName,"\t",0,1,false);
  int count = 0;
  while(!suppl.atEnd()){
    ROMol *m1 = suppl.next();
    TEST_ASSERT(m1);
    count++;
    if(!doLong && count >= 100) break;
    delete m1;
  }
  t2 = std::time(0);
  BOOST_LOG(rdInfoLog) << " Smiles time: " << std::difftime(t2,t1) << std::endl;;

  std::ifstream inStream(pklName.c_str(),std::ios_base::binary);
  t1 = std::time(0);
  while(count>0){
    ROMol m2;
    MolPickler::molFromPickle(inStream,m2);
    count--;
    if(!doLong && count >= 100) break;
  }
  t2 = std::time(0);
  BOOST_LOG(rdInfoLog) << " Pickle time: " << std::difftime(t2,t1) << std::endl;;

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}
  


void test4(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing combining molecules using pickles." << std::endl;

  ROMol *m1 = SmilesToMol("C1CC1");
  TEST_ASSERT(m1);
  ROMol *m2 = SmilesToMol("C1COC1");
  TEST_ASSERT(m2);
  ROMol m3;
  std::string pickle;
  MolPickler::pickleMol(*m1,pickle);
  MolPickler::molFromPickle(pickle,m3);
  MolPickler::pickleMol(*m2,pickle);
  MolPickler::molFromPickle(pickle,m3);

  m3.debugMol(std::cout);
  TEST_ASSERT(m3.getNumAtoms()==7);
  TEST_ASSERT(m3.getNumBonds()==7);
  TEST_ASSERT(m3.getRingInfo()->numRings()==2);
  TEST_ASSERT(m3.getRingInfo()->isAtomInRingOfSize(0,3));
  TEST_ASSERT(!m3.getRingInfo()->isAtomInRingOfSize(0,4));
  TEST_ASSERT(!m3.getRingInfo()->isAtomInRingOfSize(4,3));
  TEST_ASSERT(m3.getRingInfo()->isAtomInRingOfSize(4,4));

  TEST_ASSERT(m3.getRingInfo()->isBondInRingOfSize(0,3));
  TEST_ASSERT(!m3.getRingInfo()->isBondInRingOfSize(0,4));
  TEST_ASSERT(!m3.getRingInfo()->isBondInRingOfSize(4,3));
  TEST_ASSERT(m3.getRingInfo()->isBondInRingOfSize(4,4));

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}


void testIssue164(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing Issue164: Roundtrip Pickle Failure." << std::endl;

  std::string smi="NCCCCC(NC(C(C(C)C)NC(=O)C1CCCN1C(C(N)CCCNC(=N)N)=O)=O)C(NC(C(C)C)C(NC(Cc2ccc(O)cc2)C(=O)N6C(C(NC(CC(N)=O)C(NCC(NC(C)C(NC(CCC(O)=O)C(NC(CC(O)=O)C(NC(C(NC(CO)C(NC(C)C(NC(CCC(O)=O)C(NC(C)C(NC(Cc3ccccc3)C(=O)N5C(C(NC(CC(C)C)C(NC(C(NC(C(O)=O)Cc4ccccc4)=O)CCC(O)=O)=O)=O)CCC5)=O)=O)=O)=O)=O)CCC(O)=O)=O)=O)=O)=O)=O)=O)CCC6)=O)=O";
  ROMol *m1 = SmilesToMol(smi);
  TEST_ASSERT(m1);
  std::string pickle;
  ROMol m2;
  MolPickler::pickleMol(*m1,pickle);
  MolPickler::molFromPickle(pickle,m2);

  TEST_ASSERT(m1->getNumAtoms()==m2.getNumAtoms());

  // the issue had to do with number of atoms, so let's make an enormous
  // molecule and try again:
  RWMol *m3 = static_cast<RWMol *>(SmilesToMol(smi));
  m3->insertMol(*m1);
  m3->insertMol(*m1);
  MolOps::sanitizeMol(*m3);

  MolPickler::pickleMol(*m3,pickle);
  ROMol m4;
  MolPickler::molFromPickle(pickle,m4);

  TEST_ASSERT(m3->getNumAtoms()==m4.getNumAtoms());
  TEST_ASSERT(m3->getNumBonds()==m4.getNumBonds());
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
  


}

void testIssue219(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing Issue219: Pickle Failure with Conformations." << std::endl;

  std::string smi="CC";
  ROMol *m1 = SmilesToMol(smi);
  TEST_ASSERT(m1);
  std::string pickle;
  ROMol *m2;
  MolPickler::pickleMol(*m1,pickle);
  m2 = new ROMol();
  MolPickler::molFromPickle(pickle,*m2);
  TEST_ASSERT(m1->getNumAtoms()==m2->getNumAtoms());

  Conformer *conf = new Conformer(2);
  conf->setId(23);
  m1->addConformer(conf);
  MolPickler::pickleMol(*m1,pickle);
  delete m2;
  m2 = new ROMol();
  MolPickler::molFromPickle(pickle,*m2);
  TEST_ASSERT(m1->getNumAtoms()==m2->getNumAtoms());
  TEST_ASSERT(m2->getConformer().getId()==23);
  
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
  
}

void testIssue220(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing Issue220: Pickle problems with bond stereochemistry." << std::endl;

  std::string smi="N/N=C/C";
  ROMol *m1 = SmilesToMol(smi);
  TEST_ASSERT(m1);
  TEST_ASSERT(m1->getBondWithIdx(1)->getStereo()!=Bond::STEREONONE);
  std::string pickle;
  ROMol *m2;
  MolPickler::pickleMol(*m1,pickle);
  m2 = new ROMol();
  MolPickler::molFromPickle(pickle,*m2);
  TEST_ASSERT(m1->getNumAtoms()==m2->getNumAtoms());
  TEST_ASSERT(m1->getNumBonds()==m2->getNumBonds());
  TEST_ASSERT(m2->getBondWithIdx(1)->getStereo()!=Bond::STEREONONE);
  for(unsigned int i=0;i<m1->getNumBonds();++i){
    TEST_ASSERT(m1->getBondWithIdx(i)->getStereo() ==
		m2->getBondWithIdx(i)->getStereo());
    TEST_ASSERT(m1->getBondWithIdx(i)->getStereoAtoms() ==
		m2->getBondWithIdx(i)->getStereoAtoms());
  }
  
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
  
}


int main(int argc, char *argv[]) {
  RDLog::InitLogs();
  bool doLong=false;
  if(argc>1) {
    if(!strncmp(argv[1],"-l",2)){
      doLong=true;
    }
  }

  //_createPickleFile();
#if 1
  test1(doLong);
  //test2(doLong);
  test3(doLong);
  test4();
  testIssue164();
#endif
  //timeTest(doLong);
  testIssue219();
  testIssue220();

  return 0;

}
      
      
