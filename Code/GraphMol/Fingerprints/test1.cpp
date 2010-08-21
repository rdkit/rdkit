// $Id$
//
//  Copyright (C) 2003-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include <RDGeneral/RDLog.h>
#include <string>
#include <boost/version.hpp>

using namespace RDKit;

void test1(){
  BOOST_LOG(rdInfoLog) <<"testing basics" << std::endl;

  std::string smi = "C1=CC=CC=C1";
  RWMol *m1 = SmilesToMol(smi);
  TEST_ASSERT(m1->getNumAtoms()==6);
  ExplicitBitVect *fp1=RDKFingerprintMol(*m1);
  ExplicitBitVect *fp2=RDKFingerprintMol(*m1);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)==1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>=0.0);
  
  smi = "C1=CC=CC=N1";
  RWMol *m2 = SmilesToMol(smi);
  delete fp2;
  fp2=RDKFingerprintMol(*m2);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);


  delete m1;
  delete m2;
  delete fp1;
  delete fp2;
  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}
void test2(){
  BOOST_LOG(rdInfoLog) <<"testing subgraph invariants" << std::endl;

  std::string smi = "CC(=O)COC";
  RWMol *m1 = SmilesToMol(smi);
  TEST_ASSERT(m1->getNumAtoms()==6);
  BOOST_LOG(rdInfoLog) <<"--------------------  fp1 " << std::endl;
  ExplicitBitVect *fp1=RDKFingerprintMol(*m1,1,4,2048,4,false);
  BOOST_LOG(rdInfoLog) <<"--------------------  fp2 " << std::endl;
  ExplicitBitVect *fp2=RDKFingerprintMol(*m1,1,4,2048,4,false);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)==1.0);
  
  RWMol *m2 = SmilesToMol("CC");
  delete fp2;
  BOOST_LOG(rdInfoLog) <<"--------------------  fp2 " << std::endl;
  fp2=RDKFingerprintMol(*m2,1,4,2048,4,false);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);
  TEST_ASSERT(OnBitProjSimilarity(*fp2,*fp1)[0]==1.0);

  delete m2;
  m2 = SmilesToMol("CC=O");
  delete fp2;
  fp2=RDKFingerprintMol(*m2,1,4,2048,4,false);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);
  TEST_ASSERT(OnBitProjSimilarity(*fp2,*fp1)[0]==1.0);

  delete m2;
  m2 = SmilesToMol("CCCOC");
  delete fp2;
  fp2=RDKFingerprintMol(*m2,1,4,2048,4,false);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);
  TEST_ASSERT(OnBitProjSimilarity(*fp2,*fp1)[0]==1.0);


  delete m1;
  delete m2;
  delete fp1;
  delete fp2;

  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}

void test3(){
  BOOST_LOG(rdInfoLog) <<"testing auto folding" << std::endl;

  RWMol *m = SmilesToMol("CCCOC");

  ExplicitBitVect *fp1,*fp2;
  fp2=RDKFingerprintMol(*m,1,4,2048,4,false,
			     0.3,256);
  TEST_ASSERT(fp2->getNumBits()==256);


  delete m;
  delete fp2;
  m=SmilesToMol("CN(C)Cc1n-2c(nn1)CN=C(c1ccccc1)c1cc(Cl)ccc12");
  fp1=RDKFingerprintMol(*m,1,4,2048,4,false);
  TEST_ASSERT(fp1->getNumBits()==2048);
  fp2=RDKFingerprintMol(*m,1,4,2048,4,false,0.3,256);  
  TEST_ASSERT(fp2->getNumBits()<fp1->getNumBits());
  
  delete m;
  delete fp1;
  delete fp2;
  BOOST_LOG(rdInfoLog) <<"done" << std::endl;

}

void test1alg2(){
  BOOST_LOG(rdInfoLog) <<"testing basics alg2" << std::endl;

  std::string smi = "C1=CC=CC=C1";
  RWMol *m1 = SmilesToMol(smi);
  TEST_ASSERT(m1->getNumAtoms()==6);
  ExplicitBitVect *fp1=RDKFingerprintMol(*m1);
  ExplicitBitVect *fp2=RDKFingerprintMol(*m1);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)==1.0);
  
  smi = "C1=CC=CC=N1";
  RWMol *m2 = SmilesToMol(smi);
  delete fp2;
  fp2=RDKFingerprintMol(*m2);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);

  delete m1;
  delete m2;
  delete fp1;
  delete fp2;
  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}

void test2alg2(){
  BOOST_LOG(rdInfoLog) <<"testing subgraph invariants alg2" << std::endl;

  std::string smi = "CC(=O)COC";
  RWMol *m1 = SmilesToMol(smi);
  TEST_ASSERT(m1->getNumAtoms()==6);
  BOOST_LOG(rdInfoLog) <<"--------------------  fp1 " << std::endl;
  ExplicitBitVect *fp1=RDKFingerprintMol(*m1,1,4,2048,4,false);
  BOOST_LOG(rdInfoLog) <<"--------------------  fp2 " << std::endl;
  ExplicitBitVect *fp2=RDKFingerprintMol(*m1,1,4,2048,4,false);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)==1.0);

  RWMol *m2 = SmilesToMol("CC");
  delete fp2;
  BOOST_LOG(rdInfoLog) <<"--------------------  fp2 " << std::endl;
  fp2=RDKFingerprintMol(*m2,1,4,2048,4,false);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);
  TEST_ASSERT(OnBitProjSimilarity(*fp2,*fp1)[0]==1.0);

  delete m2;
  m2 = SmilesToMol("CC=O");
  delete fp2;
  fp2=RDKFingerprintMol(*m2,1,4,2048,4,false);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);
  TEST_ASSERT(OnBitProjSimilarity(*fp2,*fp1)[0]==1.0);

  delete m2;
  m2 = SmilesToMol("CCCOC");
  delete fp2;
  fp2=RDKFingerprintMol(*m2,1,4,2048,4,false);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);
  TEST_ASSERT(OnBitProjSimilarity(*fp2,*fp1)[0]==1.0);


  delete m1;
  delete m2;
  delete fp1;
  delete fp2;

  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}

void test4Trends(){
  BOOST_LOG(rdInfoLog) <<"testing similarity trends" << std::endl;

  double sim1,sim2;
  RWMol *m;
  ExplicitBitVect *fp1,*fp2;
  
  m = SmilesToMol("CCC");
  fp1=RDKFingerprintMol(*m);
  delete m;
  m = SmilesToMol("C1CC1");
  fp2=RDKFingerprintMol(*m);
  sim1=TanimotoSimilarity(*fp1,*fp2);

  delete m;
  m = SmilesToMol("CCCC");
  delete fp1;
  fp1=RDKFingerprintMol(*m);
  delete m;
  m = SmilesToMol("C1CCC1");
  delete fp2;
  fp2=RDKFingerprintMol(*m);
  sim2=TanimotoSimilarity(*fp1,*fp2);
  TEST_ASSERT(sim2>sim1);
  sim2=sim1;

  delete m;
  m = SmilesToMol("CCCCC");
  delete fp1;
  fp1=RDKFingerprintMol(*m);
  delete m;
  m = SmilesToMol("C1CCCC1");
  delete fp2;
  fp2=RDKFingerprintMol(*m);
  sim2=TanimotoSimilarity(*fp1,*fp2);
  TEST_ASSERT(sim2>sim1);
  sim2=sim1;

  delete m;
  m = SmilesToMol("CCCCCC");
  delete fp1;
  fp1=RDKFingerprintMol(*m);
  delete m;
  m = SmilesToMol("C1CCCCC1");
  delete fp2;
  fp2=RDKFingerprintMol(*m);
  sim2=TanimotoSimilarity(*fp1,*fp2);
  TEST_ASSERT(sim2>sim1);
  sim2=sim1;

  delete m;
  m = SmilesToMol("CCCCCCC");
  delete fp1;
  fp1=RDKFingerprintMol(*m);
  delete m;
  m = SmilesToMol("C1CCCCCC1");
  delete fp2;
  fp2=RDKFingerprintMol(*m);
  sim2=TanimotoSimilarity(*fp1,*fp2);
  TEST_ASSERT(sim2>sim1);
  sim2=sim1;


  delete m;
  delete fp1;
  delete fp2;

  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}


void test5BackwardsCompatibility(){
  BOOST_LOG(rdInfoLog) <<"testing backwards compatibility of fingerprints" << std::endl;

  RWMol *m;
  ExplicitBitVect *fp1;
  
  m = SmilesToMol("CC");
  fp1=RDKFingerprintMol(*m);
  TEST_ASSERT(fp1->getNumOnBits()==4);
#if BOOST_VERSION/100 < 1040
  // bug fixes in the uniform_int<> distribution in version 1.40
  // necessitate this
  TEST_ASSERT((*fp1)[951]);
  TEST_ASSERT((*fp1)[961]);
  TEST_ASSERT((*fp1)[1436]);
  TEST_ASSERT((*fp1)[1590]);
#else
  TEST_ASSERT((*fp1)[419]);
  TEST_ASSERT((*fp1)[718]);
  TEST_ASSERT((*fp1)[1499]);
  TEST_ASSERT((*fp1)[1504]);
#endif
  delete fp1;

  delete m;

  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}

void test1Layers(){
  BOOST_LOG(rdInfoLog) <<"testing basics layered fps" << std::endl;
#if 1
  {
    RWMol *m1 = SmilesToMol("C1=CC=CC=C1");
    RWMol *m2 = SmilesToMol("C1=CC=CC=N1");
    RWMol *m3 = SmilesToMol("C1CCCCC1");
    RWMol *m4 = SmilesToMol("C1CCC1");

    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1);
    ExplicitBitVect *fp2=LayeredFingerprintMol(*m1);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)==1.0);
  
    delete fp2;
    fp2=LayeredFingerprintMol(*m2);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);

    ExplicitBitVect *fp3=LayeredFingerprintMol(*m3);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp3)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp3)>0.0);
    TEST_ASSERT(TanimotoSimilarity(*fp2,*fp3)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp2,*fp3)>0.0);

    ExplicitBitVect *fp4=LayeredFingerprintMol(*m4);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp4)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp4)>0.0);
    TEST_ASSERT(TanimotoSimilarity(*fp2,*fp4)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp2,*fp4)>0.0);
    TEST_ASSERT(TanimotoSimilarity(*fp3,*fp4)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp3,*fp4)>0.0);

    delete fp1;delete fp2;delete fp3;delete fp4;
    delete m1;delete m2;delete m3;delete m4;
  }
#endif
  {
    RWMol *m1 = SmilesToMol("C1=CC=CC=C1");
    RWMol *m2 = SmilesToMol("C1=CC=CC=N1");
    RWMol *m3 = SmilesToMol("C1CCCCC1");
    RWMol *m4 = SmilesToMol("CCCCCC");
    unsigned int layers=0x1;
    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1,layers,1,5);

    ExplicitBitVect *fp2=LayeredFingerprintMol(*m2,layers,1,5);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)==1.0);

    ExplicitBitVect *fp3=LayeredFingerprintMol(*m3,layers,1,5);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp3)==1.0);

    ExplicitBitVect *fp4=LayeredFingerprintMol(*m4,layers,1,5);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp4)==1.0);

    delete fp1;delete fp2;delete fp3;delete fp4;
    delete m1;delete m2;delete m3;delete m4;
  }

  {
    RWMol *m1 = SmilesToMol("C1=CC=CC=C1");
    RWMol *m2 = SmilesToMol("C1=CC=CC=N1");
    RWMol *m3 = SmilesToMol("C1CCCCC1");
    RWMol *m4 = SmilesToMol("CCCCCC");
    unsigned int layers=0x3;
    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1,layers,1,5);

    ExplicitBitVect *fp2=LayeredFingerprintMol(*m2,layers,1,5);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)==1.0);

    ExplicitBitVect *fp3=LayeredFingerprintMol(*m3,layers,1,5);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp3)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp3)>0.0);

    ExplicitBitVect *fp4=LayeredFingerprintMol(*m4,layers,1,5);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp4)<=1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp4)>0.0);
    TEST_ASSERT(TanimotoSimilarity(*fp3,*fp4)==1.0);

    delete fp1;delete fp2;delete fp3;delete fp4;
    delete m1;delete m2;delete m3;delete m4;
  }

  {
    RWMol *m1 = SmilesToMol("C1=CC=CC=C1");
    RWMol *m2 = SmilesToMol("C1=CC=CC=N1");
    RWMol *m3 = SmilesToMol("C1CCCCC1");
    RWMol *m4 = SmilesToMol("CCCCCC");
    unsigned int layers=0x7;
    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1,layers,1,5);

    ExplicitBitVect *fp2=LayeredFingerprintMol(*m2,layers,1,5);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);

    ExplicitBitVect *fp3=LayeredFingerprintMol(*m3,layers,1,5);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp3)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp3)>0.0);
    TEST_ASSERT(TanimotoSimilarity(*fp2,*fp3)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp2,*fp3)>0.0);

    ExplicitBitVect *fp4=LayeredFingerprintMol(*m4,layers,1,5);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp4)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp4)>0.0);
    TEST_ASSERT(TanimotoSimilarity(*fp2,*fp4)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp2,*fp4)>0.0);
    TEST_ASSERT(TanimotoSimilarity(*fp3,*fp4)==1.0);

    delete fp1;delete fp2;delete fp3;delete fp4;
    delete m1;delete m2;delete m3;delete m4;
  }

  {
    RWMol *m1 = SmilesToMol("c1ccccc1");
    RWMol *m2 = SmilesToMol("C1CCCCC1");
    RWMol *m3 = SmilesToMol("CCCCCC");
    RWMol *m4 = SmilesToMol("C1CCCCC1");
    unsigned int layers=0xF;  // add the "in ring" bit
    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1,layers,1,5);

    ExplicitBitVect *fp2=LayeredFingerprintMol(*m2,layers,1,5);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);
    ExplicitBitVect *fp3=LayeredFingerprintMol(*m3,layers,1,5);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp3)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp3)>0.0);
    TEST_ASSERT(TanimotoSimilarity(*fp2,*fp3)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp2,*fp3)>0.0);

    ExplicitBitVect *fp4=LayeredFingerprintMol(*m4,layers,1,5);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp4)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp4)>0.0);
    TEST_ASSERT(TanimotoSimilarity(*fp2,*fp4)==1.0);

    delete fp1;delete fp2;delete fp3;delete fp4;
    delete m1;delete m2;delete m3;delete m4;
  }

  {
    RWMol *m1 = SmilesToMol("C1=CC=CC=C1");
    RWMol *m2 = SmilesToMol("C1CCCCC1");
    RWMol *m3 = SmilesToMol("CCCCCC");
    RWMol *m4 = SmilesToMol("C1CCCCCC1");
    unsigned int layers=0x1F;  // add the ring size bit
    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1,layers,1,5);

    ExplicitBitVect *fp2=LayeredFingerprintMol(*m2,layers,1,5);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);

    ExplicitBitVect *fp3=LayeredFingerprintMol(*m3,layers,1,5);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp3)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp3)>0.0);
    TEST_ASSERT(TanimotoSimilarity(*fp2,*fp3)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp2,*fp3)>0.0);

    ExplicitBitVect *fp4=LayeredFingerprintMol(*m4,layers,1,5);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp4)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp4)>0.0);
    TEST_ASSERT(TanimotoSimilarity(*fp2,*fp4)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp2,*fp4)>0.0);
    TEST_ASSERT(TanimotoSimilarity(*fp3,*fp4)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp3,*fp4)>0.0);

    delete fp1;delete fp2;delete fp3;delete fp4;
    delete m1;delete m2;delete m3;delete m4;
  }

  {
    RWMol *m1 = SmilesToMol("Cc1ncccc1");
    RWMol *m2 = SmilesToMol("Cn1ccc2nn(C)c(=O)c-2c1C");
    unsigned int layers=0x7;
    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1,layers,1,5);

    ExplicitBitVect *fp2=LayeredFingerprintMol(*m2,layers,1,5);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);
    TEST_ASSERT(((*fp1)&(*fp2))==(*fp1));

    delete fp1;delete fp2;
    delete m1;delete m2;
  }

  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}


void test2Layers(){
  BOOST_LOG(rdInfoLog) <<"testing advanced layered fps" << std::endl;
  {
    std::vector<unsigned int> atomCounts;
    RWMol *m1 = SmilesToMol("CC(C)C");

    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1);
    ExplicitBitVect *fp2=LayeredFingerprintMol(*m1);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)==1.0);
  
    atomCounts.clear();
    atomCounts.resize(m1->getNumAtoms());
    for(unsigned int i=0;i<m1->getNumAtoms();++i) atomCounts[i]=0;

    delete fp2;
    fp2=LayeredFingerprintMol(*m1,0xFFFFFFFF,1,7,2048,0,128,&atomCounts);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)==1.0);
    TEST_ASSERT(atomCounts[0]==4);
    TEST_ASSERT(atomCounts[1]==7);
    TEST_ASSERT(atomCounts[2]==4);
    TEST_ASSERT(atomCounts[3]==4);

    delete fp2;
    fp2=LayeredFingerprintMol(*m1,0xFFFFFFFF,1,7,2048,0,128,&atomCounts);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)==1.0);
    TEST_ASSERT(atomCounts[0]==8);
    TEST_ASSERT(atomCounts[1]==14);
    TEST_ASSERT(atomCounts[2]==8);
    TEST_ASSERT(atomCounts[3]==8);


    delete fp1;delete fp2;
    delete m1;
  }


  {
    std::vector<unsigned int> atomCounts;
    RWMol *m1 = SmilesToMol("CC(C)C");
    RWMol *m2 = SmilesToMol("CCC");

    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1);
    ExplicitBitVect *fp2=LayeredFingerprintMol(*m2);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
  
    atomCounts.clear();
    atomCounts.resize(m1->getNumAtoms());
    for(unsigned int i=0;i<m1->getNumAtoms();++i) atomCounts[i]=0;

    ExplicitBitVect *fp3=LayeredFingerprintMol(*m1,0xFFFFFFFF,1,7,2048,0,128,&atomCounts,fp2);
    TEST_ASSERT(TanimotoSimilarity(*fp1,*fp3)<1.0);
    TEST_ASSERT(atomCounts[0]==3);
    TEST_ASSERT(atomCounts[1]==6);
    TEST_ASSERT(atomCounts[2]==3);
    TEST_ASSERT(atomCounts[3]==3);

    delete fp1;delete fp2;delete fp3;
    delete m1;delete m2;
  }

  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}

void test3Layers(){
  BOOST_LOG(rdInfoLog) <<"testing layered fps and queries" << std::endl;
  {
    RWMol *m1 = SmilesToMol("CC(C)C");
    RWMol *m2 = SmartsToMol("CC(C)C");

    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1,0x7,1,5,256);
    ExplicitBitVect *fp2=LayeredFingerprintMol(*m2,0x7,1,5,256);
    TEST_ASSERT(fp2->getNumOnBits());
    
    ExplicitBitVect fp3=(*fp1)&(*fp2);
#if 0
    std::cerr<<BitVectToText(*fp1)<<std::endl;
    std::cerr<<"---"<<std::endl;
    std::cerr<<BitVectToText(*fp2)<<std::endl;
    std::cerr<<BitVectToText(fp3)<<std::endl;
#endif
    TEST_ASSERT(fp3==(*fp2));
  
    delete fp1;delete fp2;
    delete m1;delete m2;
  }
 
  {
    RWMol *m1 = SmilesToMol("CC(C)C");
    RWMol *m2 = SmartsToMol("C-C(-C)-C");

    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1,0x7,1,5,256);
    ExplicitBitVect *fp2=LayeredFingerprintMol(*m2,0x7,1,5,256);
    TEST_ASSERT(fp2->getNumOnBits());
    
    ExplicitBitVect fp3=(*fp1)&(*fp2);
    TEST_ASSERT(fp3==(*fp2));
  
    delete fp1;delete fp2;
    delete m1;delete m2;
  }
 
  {
    RWMol *m1 = SmilesToMol("CC(=C)C");
    RWMol *m2 = SmartsToMol("C-C(-C)-C");

    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1,0x7,1,5,256);
    ExplicitBitVect *fp2=LayeredFingerprintMol(*m2,0x7,1,5,256);
    TEST_ASSERT(fp2->getNumOnBits());
    
    ExplicitBitVect fp3=(*fp1)&(*fp2);
    TEST_ASSERT(fp3!=(*fp2));
    TEST_ASSERT(fp3!=(*fp1));
  
    delete fp1;delete fp2;
    delete m1;delete m2;
  }
 
  {
    RWMol *m1 = SmilesToMol("CC(C)C");
    RWMol *m2 = SmartsToMol("CC([C,N])C");

    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1,0x7,1,5,256);
    ExplicitBitVect *fp2=LayeredFingerprintMol(*m2,0x7,1,5,256);
    TEST_ASSERT(fp2->getNumOnBits());
    
    ExplicitBitVect fp3=(*fp1)&(*fp2);
    TEST_ASSERT(fp3==(*fp2));
  
    delete fp1;delete fp2;
    delete m1;delete m2;
  }

  {
    RWMol *m1 = SmilesToMol("CC(C)C");
    RWMol *m2 = SmartsToMol("CC([!C])C");

    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1,0x7,1,5,256);
    ExplicitBitVect *fp2=LayeredFingerprintMol(*m2,0x7,1,5,256);
    TEST_ASSERT(fp2->getNumOnBits());
    
    // the query says "!C", but we still match in the fp...
    ExplicitBitVect fp3=(*fp1)&(*fp2);
    TEST_ASSERT(fp3==(*fp2));


    // of course we better match if the thing really isn't a C:
    delete fp1;delete m1;
    m1 = SmilesToMol("CC(O)C");
    fp1=LayeredFingerprintMol(*m1,0x7,1,5,256);
    fp3=(*fp1)&(*fp2);
    TEST_ASSERT(fp3==(*fp2));
  
    delete fp1;delete fp2;
    delete m1;delete m2;
  }

  {
    RWMol *m1 = SmilesToMol("CC(O)C");
    RWMol *m2 = SmartsToMol("CC([C,N])C");

    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1,0x7,1,5,256);
    ExplicitBitVect *fp2=LayeredFingerprintMol(*m2,0x7,1,5,256);
    TEST_ASSERT(fp2->getNumOnBits());
    
    ExplicitBitVect fp3=(*fp1)&(*fp2);
    TEST_ASSERT(fp3==(*fp2));
  
    delete fp1;delete fp2;
    delete m1;delete m2;
  }
 
  {
    RWMol *m1 = SmilesToMol("c1ccccc1");
    RWMol *m2 = SmartsToMol("cccc");
    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1,0x7,1,5,256);
    ExplicitBitVect *fp2=LayeredFingerprintMol(*m2,0x7,1,5,256);
    TEST_ASSERT(fp2->getNumOnBits());
    
    ExplicitBitVect fp3=(*fp1)&(*fp2);
    TEST_ASSERT(fp3==(*fp2));

    delete fp1;delete m1;
    // somewhat counter-intuitive result here, but this happens
    // because the aromaticity does not factor into the
    // calculation:
    m1 = SmilesToMol("C1CCCCC1");
    fp1=LayeredFingerprintMol(*m1,0x7,1,5,256);
    fp3=(*fp1)&(*fp2);
    TEST_ASSERT(fp3==(*fp2));
    
    delete fp1;delete fp2;
    delete m1;delete m2;
  }
 
  {
    RWMol *m1 = SmilesToMol("c1ccccc1");
    RWMol *m2 = SmartsToMol("c:c:c");
    ExplicitBitVect *fp1=LayeredFingerprintMol(*m1,0xFFFF,1,5,512);
    ExplicitBitVect *fp2=LayeredFingerprintMol(*m2,0x5,1,5,512);
    TEST_ASSERT(fp2->getNumOnBits());
    ExplicitBitVect fp3=(*fp1)&(*fp2);
    TEST_ASSERT(fp3==(*fp2));

    delete fp1;delete m1;
    m1 = SmilesToMol("C1CCCCC1");
    fp1=LayeredFingerprintMol(*m1,0xFFFF,1,5,512);
    fp3=(*fp1)&(*fp2);
    TEST_ASSERT(fp3==(*fp2));
    delete fp2;
    fp2=LayeredFingerprintMol(*m2,0x7,1,5,512);
    fp3=(*fp1)&(*fp2);
    TEST_ASSERT(fp3!=(*fp2));

    delete fp2;
    fp2=LayeredFingerprintMol(*m2,0x25,1,5,512);
    fp3=(*fp1)&(*fp2);
    TEST_ASSERT(fp3!=(*fp2));

    delete fp2;
    fp2=LayeredFingerprintMol(*m2,0x7,1,5,512);
    delete fp1;delete m1;
    // atom type information is factored in:
    m1 = SmilesToMol("c1ncncn1");
    fp1=LayeredFingerprintMol(*m1,0xFFFF,1,5,512);
    fp3=(*fp1)&(*fp2);
    TEST_ASSERT(fp3!=(*fp2));

    delete fp2;delete m2;
    // but queries switch us back to generic types, so we match again:
    m2 = SmartsToMol("[C,N]:[C,N]:[C,N]");
    fp2=LayeredFingerprintMol(*m2,0x7,1,5,512);
    TEST_ASSERT(fp2->getNumOnBits());
    fp3=(*fp1)&(*fp2);
    TEST_ASSERT(fp3==(*fp2));
    
    delete fp1;delete fp2;
    delete m1;delete m2;
  }
 
 
 BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}


void test1MorganFPs(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Morgan Fingerprints." << std::endl;

  {
    ROMol *mol;
    SparseIntVect<boost::uint32_t> *fp;

    mol = SmilesToMol("CCCCC");
    fp = MorganFingerprints::getFingerprint(*mol,0);
    TEST_ASSERT(fp->getNonzeroElements().size()==2);
    delete fp;
    fp = MorganFingerprints::getFingerprint(*mol,1);
    TEST_ASSERT(fp->getNonzeroElements().size()==5);
    delete fp;
    fp = MorganFingerprints::getFingerprint(*mol,2);
    TEST_ASSERT(fp->getNonzeroElements().size()==7);
    delete fp;
    fp = MorganFingerprints::getFingerprint(*mol,3);
    TEST_ASSERT(fp->getNonzeroElements().size()==7);
    delete fp;
  
  
    delete mol;
  }
  {
    ROMol *mol;
    SparseIntVect<boost::uint32_t> *fp;

    mol = SmilesToMol("O=C(O)CC1CC1");
    fp = MorganFingerprints::getFingerprint(*mol,0);
    TEST_ASSERT(fp->getNonzeroElements().size()==5);
    delete fp;
    fp = MorganFingerprints::getFingerprint(*mol,1);
    TEST_ASSERT(fp->getNonzeroElements().size()==11);
    delete fp;
    fp = MorganFingerprints::getFingerprint(*mol,2);
    TEST_ASSERT(fp->getNonzeroElements().size()==15);
    delete fp;
    fp = MorganFingerprints::getFingerprint(*mol,3);
    TEST_ASSERT(fp->getNonzeroElements().size()==16);
    delete fp;

    delete mol;
  }

  {
    // test that the results aren't order dependent, i.e. that we're
    // "canonicalizing" the fps correctly
    ROMol *mol,*mol2;
    SparseIntVect<boost::uint32_t> *fp,*fp2;

    mol = SmilesToMol("O=C(O)CC1CC1");
    mol2 = SmilesToMol("OC(=O)CC1CC1");
    fp = MorganFingerprints::getFingerprint(*mol,0);
    fp2 = MorganFingerprints::getFingerprint(*mol2,0);
    TEST_ASSERT(fp->getNonzeroElements().size()==5);
    TEST_ASSERT(fp2->getNonzeroElements().size()==5);
    TEST_ASSERT(*fp==*fp2);
    delete fp;
    delete fp2;

    fp = MorganFingerprints::getFingerprint(*mol,1);
    fp2 = MorganFingerprints::getFingerprint(*mol2,1);
    TEST_ASSERT(fp->getNonzeroElements().size()==11);
    TEST_ASSERT(fp2->getNonzeroElements().size()==11);
    TEST_ASSERT(*fp==*fp2);
    delete fp;
    delete fp2;

    fp = MorganFingerprints::getFingerprint(*mol,2);
    fp2 = MorganFingerprints::getFingerprint(*mol2,2);
    TEST_ASSERT(fp->getNonzeroElements().size()==15);
    TEST_ASSERT(fp2->getNonzeroElements().size()==15);
    TEST_ASSERT(*fp==*fp2);
    delete fp;
    delete fp2;

    fp = MorganFingerprints::getFingerprint(*mol,3);
    fp2 = MorganFingerprints::getFingerprint(*mol2,3);
    TEST_ASSERT(fp->getNonzeroElements().size()==16);
    TEST_ASSERT(fp2->getNonzeroElements().size()==16);
    TEST_ASSERT(*fp==*fp2);
    delete fp;
    delete fp2;

    delete mol;
    delete mol2;
  }

  {
    // symmetry test:
    ROMol *mol;
    SparseIntVect<boost::uint32_t> *fp;

    mol = SmilesToMol("OCCCCO");
    fp = MorganFingerprints::getFingerprint(*mol,2);
    TEST_ASSERT(fp->getNonzeroElements().size()==7);
    SparseIntVect<boost::uint32_t>::StorageType::const_iterator iter;
    for(iter=fp->getNonzeroElements().begin();
        iter!=fp->getNonzeroElements().end();++iter){
      TEST_ASSERT(iter->second==2 || iter->second==4);
    }
    
    delete fp;
    delete mol;
  }

  {
    // chirality test:
    ROMol *mol;
    SparseIntVect<boost::uint32_t> *fp;

    mol = SmilesToMol("CC(F)(Cl)C(F)(Cl)C");
    fp = MorganFingerprints::getFingerprint(*mol,0);
    TEST_ASSERT(fp->getNonzeroElements().size()==4);
    delete fp;
    fp = MorganFingerprints::getFingerprint(*mol,1);
    TEST_ASSERT(fp->getNonzeroElements().size()==8);
    delete mol;
    delete fp;

    mol = SmilesToMol("CC(F)(Cl)[C@](F)(Cl)C");
    fp = MorganFingerprints::getFingerprint(*mol,0);
    TEST_ASSERT(fp->getNonzeroElements().size()==4);
    delete fp;
    fp = MorganFingerprints::getFingerprint(*mol,1);
    TEST_ASSERT(fp->getNonzeroElements().size()==9);
    delete fp;
    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void test2MorganFPsFromAtoms(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Morgan Fingerprints using fromAtoms argument." << std::endl;

  {
    ROMol *mol;
    SparseIntVect<boost::uint32_t> *fp;
    std::vector<boost::uint32_t> atoms;
    atoms.push_back(0);
    
    mol = SmilesToMol("CCCCC");
    fp = MorganFingerprints::getFingerprint(*mol,0);
    TEST_ASSERT(fp->getNonzeroElements().size()==2);
    delete fp;

    fp = MorganFingerprints::getFingerprint(*mol,0,(std::vector<boost::uint32_t> *)NULL,&atoms);
    TEST_ASSERT(fp->getNonzeroElements().size()==1);
    delete fp;
  
    fp = MorganFingerprints::getFingerprint(*mol,1,(std::vector<boost::uint32_t> *)NULL,&atoms);
    TEST_ASSERT(fp->getNonzeroElements().size()==2);
    delete fp;
  
    fp = MorganFingerprints::getFingerprint(*mol,2,(std::vector<boost::uint32_t> *)NULL,&atoms);
    TEST_ASSERT(fp->getNonzeroElements().size()==2);
    delete fp;
  
    delete mol;
  }

  {
    ROMol *mol;
    SparseIntVect<boost::uint32_t> *fp;
    std::vector<boost::uint32_t> atoms;
    
    mol = SmilesToMol("CCCCC");
    fp = MorganFingerprints::getFingerprint(*mol,0,(std::vector<boost::uint32_t> *)0,&atoms);
    TEST_ASSERT(fp->getNonzeroElements().size()==0);
    delete fp;
  
    fp = MorganFingerprints::getFingerprint(*mol,1,(std::vector<boost::uint32_t> *)0,&atoms);
    TEST_ASSERT(fp->getNonzeroElements().size()==0);
    delete fp;
  
    delete mol;
  }

  {
    ROMol *mol;
    SparseIntVect<boost::uint32_t> *fp;
    std::vector<boost::uint32_t> atoms;
    atoms.push_back(0);
    
    mol = SmilesToMol("C(CC)CO");

    fp = MorganFingerprints::getFingerprint(*mol,0,(std::vector<boost::uint32_t> *)0,&atoms);
    TEST_ASSERT(fp->getNonzeroElements().size()==1);
    delete fp;
  
    fp = MorganFingerprints::getFingerprint(*mol,1,(std::vector<boost::uint32_t> *)0,&atoms);
    TEST_ASSERT(fp->getNonzeroElements().size()==2);
    delete fp;
  
    fp = MorganFingerprints::getFingerprint(*mol,2,(std::vector<boost::uint32_t> *)0,&atoms);
    TEST_ASSERT(fp->getNonzeroElements().size()==3);
    delete fp;
  
    fp = MorganFingerprints::getFingerprint(*mol,3,(std::vector<boost::uint32_t> *)0,&atoms);
    TEST_ASSERT(fp->getNonzeroElements().size()==3);
    delete fp;
  
    delete mol;
  }
  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void test3MorganFPs(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Morgan Fingerprints as bit vects." << std::endl;

  {
    ROMol *mol;
    ExplicitBitVect *fp;

    mol = SmilesToMol("CCCCC");
    fp = MorganFingerprints::getFingerprintAsBitVect(*mol,0,2048);
    TEST_ASSERT(fp->getNumOnBits()==2);
    delete fp;
    fp = MorganFingerprints::getFingerprintAsBitVect(*mol,1,2048);
    TEST_ASSERT(fp->getNumOnBits()==5);
    delete fp;
    fp = MorganFingerprints::getFingerprintAsBitVect(*mol,2,2048);
    TEST_ASSERT(fp->getNumOnBits()==7);
    delete fp;
    fp = MorganFingerprints::getFingerprintAsBitVect(*mol,3,2048);
    TEST_ASSERT(fp->getNumOnBits()==7);
    delete fp;
  
  
    delete mol;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testIssue2875658(){
  BOOST_LOG(rdInfoLog) <<"testing issue 2875658" << std::endl;

  RWMol *m;
  ExplicitBitVect *fp1,*fp2;
  
  m = SmilesToMol("c1ccccc1");
  fp1=RDKFingerprintMol(*m);
  delete m;
  m = SmilesToMol("c1cccnc1");
  fp2=RDKFingerprintMol(*m);
  delete m;

#if BOOST_VERSION/100 >= 1040
  // bug fixes in the uniform_int<> distribution in version 1.40
  // necessitate this
  TEST_ASSERT(fp1->getNumOnBits()==24);
  TEST_ASSERT((*fp1)[11]);
  TEST_ASSERT((*fp1)[857]);
  TEST_ASSERT((*fp1)[1786]);
  TEST_ASSERT((*fp1)[2020]);

  TEST_ASSERT(fp2->getNumOnBits()==64);
  TEST_ASSERT((*fp2)[11]);
  TEST_ASSERT((*fp2)[179]);
  TEST_ASSERT((*fp2)[1878]);
  TEST_ASSERT((*fp2)[2020]);
#else
  TEST_ASSERT(fp1->getNumOnBits()==24);
  TEST_ASSERT((*fp1)[23]);
  TEST_ASSERT((*fp1)[434]);
  TEST_ASSERT((*fp1)[1445]);
  TEST_ASSERT((*fp1)[2031]);

  TEST_ASSERT(fp2->getNumOnBits()==62);
  TEST_ASSERT((*fp2)[23]);
  TEST_ASSERT((*fp2)[173]);
  TEST_ASSERT((*fp2)[1847]);
  TEST_ASSERT((*fp2)[1975]);
#endif
  delete fp1;
  delete fp2;
}

void testAtomCodes(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Atom Codes." << std::endl;

  ROMol *mol;
  boost::uint32_t tgt;
  mol = SmilesToMol("C=C");
  TEST_ASSERT(AtomPairs::getAtomCode(mol->getAtomWithIdx(0))==AtomPairs::getAtomCode(mol->getAtomWithIdx(1)));
  tgt = 1 | (1 | 1<<AtomPairs::numPiBits)<<AtomPairs::numBranchBits;
  TEST_ASSERT(AtomPairs::getAtomCode(mol->getAtomWithIdx(0))==tgt);
  tgt = 1<<AtomPairs::numBranchBits | 1<<(AtomPairs::numBranchBits+AtomPairs::numPiBits);
  TEST_ASSERT(AtomPairs::getAtomCode(mol->getAtomWithIdx(0),1)==tgt);

  delete mol;
  mol = SmilesToMol("C#CO");
  tgt = 1 | 2<<AtomPairs::numBranchBits | 1<<(AtomPairs::numBranchBits+AtomPairs::numPiBits);
  TEST_ASSERT(AtomPairs::getAtomCode(mol->getAtomWithIdx(0))==tgt);
  tgt = 2 | 2<<AtomPairs::numBranchBits | 1<<(AtomPairs::numBranchBits+AtomPairs::numPiBits);
  TEST_ASSERT(AtomPairs::getAtomCode(mol->getAtomWithIdx(1))==tgt);
  tgt = 1 | 0<<AtomPairs::numBranchBits | 3<<(AtomPairs::numBranchBits+AtomPairs::numPiBits);
  TEST_ASSERT(AtomPairs::getAtomCode(mol->getAtomWithIdx(2))==tgt);

  delete mol;
  mol = SmilesToMol("CC(O)C(O)(O)C");
  tgt = 1 | 0<<AtomPairs::numBranchBits | 1<<(AtomPairs::numBranchBits+AtomPairs::numPiBits);
  TEST_ASSERT(AtomPairs::getAtomCode(mol->getAtomWithIdx(1),2)==tgt);
  tgt = 2 | 0<<AtomPairs::numBranchBits | 1<<(AtomPairs::numBranchBits+AtomPairs::numPiBits);
  TEST_ASSERT(AtomPairs::getAtomCode(mol->getAtomWithIdx(3),2)==tgt);
  
  delete mol;
  mol = SmilesToMol("C=CC(=O)O");
  tgt = 0 | 0<<AtomPairs::numBranchBits | 3<<(AtomPairs::numBranchBits+AtomPairs::numPiBits);
  TEST_ASSERT(AtomPairs::getAtomCode(mol->getAtomWithIdx(4),1)==tgt);
  tgt = 3 | 1<<AtomPairs::numBranchBits | 1<<(AtomPairs::numBranchBits+AtomPairs::numPiBits);
  TEST_ASSERT(AtomPairs::getAtomCode(mol->getAtomWithIdx(2))==tgt);


  delete mol;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testAtomPairs(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Atom Pairs." << std::endl;

  ROMol *mol;
  SparseIntVect<boost::int32_t> *fp;
  boost::uint32_t tgt;
  boost::uint32_t c1,c2,c3;

  mol = SmilesToMol("CCCCC");
  c1=AtomPairs::getAtomCode(mol->getAtomWithIdx(0));
  c2=AtomPairs::getAtomCode(mol->getAtomWithIdx(1));
  c3=AtomPairs::getAtomCode(mol->getAtomWithIdx(2));
  tgt = 1 | (std::min(c1,c2) | std::max(c1,c2)<<AtomPairs::codeSize)<< AtomPairs::numPathBits;
  TEST_ASSERT(AtomPairs::getAtomPairCode(c1,c2,1)==tgt);
  TEST_ASSERT(AtomPairs::getAtomPairCode(c2,c1,1)==tgt);
  tgt = 2 | (std::min(c1,c3) | std::max(c1,c3)<<AtomPairs::codeSize)<< AtomPairs::numPathBits;
  TEST_ASSERT(AtomPairs::getAtomPairCode(c1,c3,2)==tgt);
  TEST_ASSERT(AtomPairs::getAtomPairCode(c3,c1,2)==tgt);

  delete mol;
  mol = SmilesToMol("CCC");
  fp=AtomPairs::getAtomPairFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal()==3);
  TEST_ASSERT(fp->getNonzeroElements().size()==2);

  c1=AtomPairs::getAtomCode(mol->getAtomWithIdx(0));
  c2=AtomPairs::getAtomCode(mol->getAtomWithIdx(1));
  c3=AtomPairs::getAtomCode(mol->getAtomWithIdx(2));
  TEST_ASSERT(fp->getVal(AtomPairs::getAtomPairCode(c1,c2,1))==2);
  TEST_ASSERT(fp->getVal(AtomPairs::getAtomPairCode(c1,c3,2))==1);
  
  delete mol;
  delete fp;
  mol = SmilesToMol("CC=O.Cl");
  fp=AtomPairs::getAtomPairFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal()==3);
  TEST_ASSERT(fp->getNonzeroElements().size()==3);

  c1=AtomPairs::getAtomCode(mol->getAtomWithIdx(0));
  c2=AtomPairs::getAtomCode(mol->getAtomWithIdx(1));
  c3=AtomPairs::getAtomCode(mol->getAtomWithIdx(2));
  TEST_ASSERT(fp->getVal(AtomPairs::getAtomPairCode(c1,c2,1))==1);
  TEST_ASSERT(fp->getVal(AtomPairs::getAtomPairCode(c1,c2,1))==1);
  TEST_ASSERT(fp->getVal(AtomPairs::getAtomPairCode(c2,c3,1))==1);
  
  delete mol;
  delete fp;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testAtomPairs2(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Atom Pairs part 2." << std::endl;

  {
    ROMol *mol;
    SparseIntVect<boost::int32_t> *fp;

    mol = SmilesToMol("CCC");
    fp=AtomPairs::getAtomPairFingerprint(*mol,1,2);
    TEST_ASSERT(fp->getTotalVal()==3);
    TEST_ASSERT(fp->getNonzeroElements().size()==2);
  
    delete fp;
    fp=AtomPairs::getAtomPairFingerprint(*mol,2,2);
    TEST_ASSERT(fp->getTotalVal()==1);
    TEST_ASSERT(fp->getNonzeroElements().size()==1);
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testHashedAtomPairs(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Hashed Atom Pairs." << std::endl;

  {
    ROMol *mol;
    mol = SmilesToMol("c1ccccc1");
    SparseIntVect<boost::int32_t> *fp1;
    fp1=AtomPairs::getHashedAtomPairFingerprint(*mol);
    SparseIntVect<boost::int32_t> *fp2;
    fp2=AtomPairs::getHashedAtomPairFingerprint(*mol);
    TEST_ASSERT(DiceSimilarity(*fp1,*fp2)==1.0);
    TEST_ASSERT(*fp1==*fp2);

    delete mol;
    delete fp2;
    mol = SmilesToMol("c1ccccn1");
    fp2=AtomPairs::getHashedAtomPairFingerprint(*mol);
    RANGE_CHECK(0.0,DiceSimilarity(*fp1,*fp2),1.0);

    delete mol;
    delete fp1;
    delete fp2;
  }

  {
    ROMol *mol;
    mol = SmilesToMol("c1ccccc1");
    SparseIntVect<boost::int32_t> *fp1;
    fp1=AtomPairs::getHashedAtomPairFingerprint(*mol,2048);
    SparseIntVect<boost::int32_t> *fp2;
    fp2=AtomPairs::getHashedAtomPairFingerprint(*mol,2048,1,3);
    TEST_ASSERT(DiceSimilarity(*fp1,*fp2)==1.0);
    TEST_ASSERT(*fp1==*fp2);

    delete fp2;
    fp2=AtomPairs::getHashedAtomPairFingerprint(*mol,2048,1,2);
    RANGE_CHECK(0.0,DiceSimilarity(*fp1,*fp2),1.0);

    delete mol;
    delete fp1;
    delete fp2;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void testTorsions(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Topological Torsions." << std::endl;

  ROMol *mol;
  SparseIntVect<boost::int64_t> *fp;
  boost::uint64_t tgt;
  boost::uint64_t  c1,c2,c3,c4;
  std::vector<boost::uint32_t> codes;

  mol = SmilesToMol("CCCC");
  c1=AtomPairs::getAtomCode(mol->getAtomWithIdx(0))-1;
  c2=AtomPairs::getAtomCode(mol->getAtomWithIdx(1))-2;
  c3=AtomPairs::getAtomCode(mol->getAtomWithIdx(2))-2;
  c4=AtomPairs::getAtomCode(mol->getAtomWithIdx(3))-1;
  tgt = c1 | (c2 | (c3 | c4<<AtomPairs::codeSize)<<AtomPairs::codeSize)<<AtomPairs::codeSize;
  codes.clear();
  codes.push_back(static_cast<unsigned int>(c1));
  codes.push_back(static_cast<unsigned int>(c2));
  codes.push_back(static_cast<unsigned int>(c3));
  codes.push_back(static_cast<unsigned int>(c4));
  TEST_ASSERT(AtomPairs::getTopologicalTorsionCode(codes)==tgt);

  fp = AtomPairs::getTopologicalTorsionFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal()==1);
  TEST_ASSERT(fp->getNonzeroElements().size()==1);


  delete mol;
  delete fp;
  mol = SmilesToMol("CCCCO.Cl");
  fp=AtomPairs::getTopologicalTorsionFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal()==2);
  TEST_ASSERT(fp->getNonzeroElements().size()==2);

  delete fp;
  fp = AtomPairs::getTopologicalTorsionFingerprint(*mol,3);
  TEST_ASSERT(fp->getTotalVal()==3);
  TEST_ASSERT(fp->getNonzeroElements().size()==3);

  delete mol;
  delete fp;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testHashedTorsions(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Hashed torsions." << std::endl;

  {
    ROMol *mol;
    mol = SmilesToMol("c1ccccc1");
    SparseIntVect<boost::int64_t> *fp1;
    fp1=AtomPairs::getHashedTopologicalTorsionFingerprint(*mol);
    SparseIntVect<boost::int64_t> *fp2;
    fp2=AtomPairs::getHashedTopologicalTorsionFingerprint(*mol);
    TEST_ASSERT(DiceSimilarity(*fp1,*fp2)==1.0);
    TEST_ASSERT(*fp1==*fp2);

    delete mol;
    delete fp2;
    mol = SmilesToMol("c1ccccn1");
    fp2=AtomPairs::getHashedTopologicalTorsionFingerprint(*mol);
    RANGE_CHECK(0.0,DiceSimilarity(*fp1,*fp2),1.0);

    delete mol;
    delete fp1;
    delete fp2;
  }

  {
    ROMol *mol;
    mol = SmilesToMol("c1ccccc1");
    SparseIntVect<boost::int64_t> *fp1;
    fp1=AtomPairs::getHashedTopologicalTorsionFingerprint(*mol,2048,6);
    SparseIntVect<boost::int64_t> *fp2;
    fp2=AtomPairs::getHashedTopologicalTorsionFingerprint(*mol,2048,6);
    TEST_ASSERT(DiceSimilarity(*fp1,*fp2)==1.0);
    TEST_ASSERT(*fp1==*fp2);

    delete mol;
    delete fp2;
    mol = SmilesToMol("c1ccccn1");
    fp2=AtomPairs::getHashedTopologicalTorsionFingerprint(*mol,2048,6);
    RANGE_CHECK(0.0,DiceSimilarity(*fp1,*fp2),1.0);

    delete mol;
    delete fp1;
    delete fp2;
  }


  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}




void testBulkTorsions(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Bulk Topological Torsions." << std::endl;

  std::string fName= getenv("RDBASE");
  fName += "/Projects/DbCLI/testData/pubchem.200.sdf";
  SDMolSupplier suppl(fName);
  while(!suppl.atEnd()){
    ROMol *mol=suppl.next();
    SparseIntVect<boost::int64_t> *fp;
    fp = AtomPairs::getTopologicalTorsionFingerprint(*mol);
    TEST_ASSERT(fp->getTotalVal()>1);
    delete mol;
    delete fp;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testRootedAtomPairs(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Rooted Atom Pairs." << std::endl;

  ROMol *mol;
  SparseIntVect<boost::int32_t> *fp1,*fp2;
  std::vector<boost::uint32_t> roots;

  mol = SmilesToMol("OCCCCC");
  fp1=AtomPairs::getAtomPairFingerprint(*mol);
  SparseIntVect<boost::int32_t>::StorageType nz1=fp1->getNonzeroElements();
  TEST_ASSERT(nz1.size()>0);

  roots.push_back(0);
  fp2=AtomPairs::getAtomPairFingerprint(*mol,&roots);
  SparseIntVect<boost::int32_t>::StorageType nz2=fp2->getNonzeroElements();
  TEST_ASSERT(nz2.size()>0);
  TEST_ASSERT(nz2.size()<nz1.size());

  for(SparseIntVect<boost::int32_t>::StorageType::const_iterator bIt=nz2.begin();
      bIt!=nz2.end();++bIt){
    TEST_ASSERT(bIt->second<=fp2->getVal(bIt->first));
  }

  delete mol;
  delete fp1;
  delete fp2;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void testRootedTorsions(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Rooted Topological Torsions." << std::endl;

  ROMol *mol;
  SparseIntVect<boost::int64_t> *fp1,*fp2;
  std::vector<boost::uint32_t> roots;


  mol = SmilesToMol("OCCCC");
  roots.push_back(0);

  fp1 = AtomPairs::getTopologicalTorsionFingerprint(*mol);
  SparseIntVect<boost::int64_t>::StorageType nz1=fp1->getNonzeroElements();
  TEST_ASSERT(nz1.size()>0);

  fp2 = AtomPairs::getTopologicalTorsionFingerprint(*mol,4,&roots);
  SparseIntVect<boost::int64_t>::StorageType nz2=fp2->getNonzeroElements();
  TEST_ASSERT(nz2.size()>0);
  TEST_ASSERT(nz2.size()<nz1.size());

  for(SparseIntVect<boost::int64_t>::StorageType::const_iterator bIt=nz2.begin();
      bIt!=nz2.end();++bIt){
    TEST_ASSERT(bIt->second<=fp2->getVal(bIt->first));
  }

  delete mol;
  delete fp1;
  delete fp2;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}



int main(int argc,char *argv[]){
  RDLog::InitLogs();
#if 1
  test1();
  test2();
  test3();
  test1alg2();
  test2alg2();
  test4Trends();
  test1Layers();
  test2Layers();
  test1MorganFPs();
  test2MorganFPsFromAtoms();
  test3MorganFPs();
  test3Layers();
#endif
  test5BackwardsCompatibility();
  testIssue2875658();
  testAtomCodes();
  testAtomPairs();
  testAtomPairs2();
  testTorsions();
  testBulkTorsions();
  testHashedAtomPairs();
  testHashedTorsions();
  testRootedAtomPairs();
  testRootedTorsions();
  return 0;
}
