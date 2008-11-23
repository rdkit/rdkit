// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include <RDGeneral/RDLog.h>
#include <string>

using namespace RDKit;

void test1(){
  BOOST_LOG(rdInfoLog) <<"testing basics" << std::endl;

  std::string smi = "C1=CC=CC=C1";
  RWMol *m1 = SmilesToMol(smi);
  TEST_ASSERT(m1->getNumAtoms()==6);
  ExplicitBitVect *fp1=DaylightFingerprintMol(*m1);
  ExplicitBitVect *fp2=DaylightFingerprintMol(*m1);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)==1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>=0.0);
  
  smi = "C1=CC=CC=N1";
  RWMol *m2 = SmilesToMol(smi);
  delete fp2;
  fp2=DaylightFingerprintMol(*m2);
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
  ExplicitBitVect *fp1=DaylightFingerprintMol(*m1,1,4,2048,4,false);
  BOOST_LOG(rdInfoLog) <<"--------------------  fp2 " << std::endl;
  ExplicitBitVect *fp2=DaylightFingerprintMol(*m1,1,4,2048,4,false);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)==1.0);
  
  RWMol *m2 = SmilesToMol("CC");
  delete fp2;
  BOOST_LOG(rdInfoLog) <<"--------------------  fp2 " << std::endl;
  fp2=DaylightFingerprintMol(*m2,1,4,2048,4,false);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);
  TEST_ASSERT(OnBitProjSimilarity(*fp2,*fp1)[0]==1.0);

  delete m2;
  m2 = SmilesToMol("CC=O");
  delete fp2;
  fp2=DaylightFingerprintMol(*m2,1,4,2048,4,false);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);
  TEST_ASSERT(OnBitProjSimilarity(*fp2,*fp1)[0]==1.0);

  delete m2;
  m2 = SmilesToMol("CCCOC");
  delete fp2;
  fp2=DaylightFingerprintMol(*m2,1,4,2048,4,false);
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
  fp2=DaylightFingerprintMol(*m,1,4,2048,4,false,
			     0.3,256);
  TEST_ASSERT(fp2->GetNumBits()==256);


  delete m;
  delete fp2;
  m=SmilesToMol("CN(C)Cc1n-2c(nn1)CN=C(c1ccccc1)c1cc(Cl)ccc12");
  fp1=DaylightFingerprintMol(*m,1,4,2048,4,false);
  TEST_ASSERT(fp1->GetNumBits()==2048);
  fp2=DaylightFingerprintMol(*m,1,4,2048,4,false,0.3,256);  
  TEST_ASSERT(fp2->GetNumBits()<fp1->GetNumBits());
  
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
  smi = "C1=CC=CC=N1";
  RWMol *m2 = SmilesToMol(smi);

  ExplicitBitVect *fp1=RDKFingerprintMol(*m1);
  ExplicitBitVect *fp2=RDKFingerprintMol(*m1);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)==1.0);
  
  delete fp2;
  fp2=RDKFingerprintMol(*m2);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);

  delete fp1;
  delete fp2;

  delete m1;
  delete m2;
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
  TEST_ASSERT(fp1->GetNumOnBits()==4);
  TEST_ASSERT((*fp1)[951]);
  TEST_ASSERT((*fp1)[961]);
  TEST_ASSERT((*fp1)[1436]);
  TEST_ASSERT((*fp1)[1590]);
  delete fp1;

#if 1 
  // boost 1.35.0
  fp1=DaylightFingerprintMol(*m);
  CHECK_INVARIANT(fp1->GetNumOnBits()==4,"Fingerprint compatibility problem detected");
  CHECK_INVARIANT((*fp1)[28],"Fingerprint compatibility problem detected");
  CHECK_INVARIANT((*fp1)[1243],"Fingerprint compatibility problem detected");
  CHECK_INVARIANT((*fp1)[1299],"Fingerprint compatibility problem detected");
  CHECK_INVARIANT((*fp1)[1606],"Fingerprint compatibility problem detected");
  delete fp1;
#else
  // boost 1.34.1 and earlier
  fp1=DaylightFingerprintMol(*m);
  CHECK_INVARIANT(fp1->GetNumOnBits()==4,"Fingerprint compatibility problem detected");
  CHECK_INVARIANT((*fp1)[1141],"Fingerprint compatibility problem detected");
  CHECK_INVARIANT((*fp1)[1317],"Fingerprint compatibility problem detected");
  CHECK_INVARIANT((*fp1)[1606],"Fingerprint compatibility problem detected");
  CHECK_INVARIANT((*fp1)[1952],"Fingerprint compatibility problem detected");
  delete fp1;
#endif
  delete m;

  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}

void test1Layers(){
  BOOST_LOG(rdInfoLog) <<"testing basics layered fps" << std::endl;
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

  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}



int main(int argc,char *argv[]){
  RDLog::InitLogs();
  test1();
  test2();
  test3();
  test1alg2();
  test2alg2();
  test4Trends();
  test5BackwardsCompatibility();
  test1Layers();
  return 0;
}
