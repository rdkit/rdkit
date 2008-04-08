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

  INT_VECT::const_iterator i;
  INT_VECT v1,v2;
  fp1->GetOnBits(v1);
  //for(i=v1.begin();i!=v1.end();i++){
  //  BOOST_LOG(rdInfoLog) <<*i << " ";
  //}
  //BOOST_LOG(rdInfoLog) <<std::endl;
  
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
  fp1=DaylightFingerprintMol(*m,1,4,2048,4,false,
			     0.3,256);
  TEST_ASSERT(fp1->GetNumBits()==256);
  TEST_ASSERT(fp1->GetNumOnBits()==29);

  delete m;
  delete fp1;
  m=SmilesToMol("CN(C)Cc1n-2c(nn1)CN=C(c1ccccc1)c1cc(Cl)ccc12");
  fp1=DaylightFingerprintMol(*m,1,4,2048,4,false);
  TEST_ASSERT(fp1->GetNumBits()==2048);
  TEST_ASSERT(fp1->GetNumOnBits()==334);
  fp2=DaylightFingerprintMol(*m,1,4,2048,4,false,0.3,256);  
  TEST_ASSERT(fp2->GetNumBits()==1024);
  TEST_ASSERT(fp2->GetNumOnBits()==309);
  
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
  ExplicitBitVect *fp1=RDKFingerprintMol2(*m1);
  ExplicitBitVect *fp2=RDKFingerprintMol2(*m1);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)==1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>=0.0);
  
  smi = "C1=CC=CC=N1";
  RWMol *m2 = SmilesToMol(smi);
  delete fp2;
  fp2=RDKFingerprintMol2(*m2);
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
  ExplicitBitVect *fp1=RDKFingerprintMol2(*m1,1,4,2048,4,false);
  BOOST_LOG(rdInfoLog) <<"--------------------  fp2 " << std::endl;
  ExplicitBitVect *fp2=RDKFingerprintMol2(*m1,1,4,2048,4,false);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)==1.0);

  INT_VECT::const_iterator i;
  INT_VECT v1,v2;
  fp1->GetOnBits(v1);
  //for(i=v1.begin();i!=v1.end();i++){
  //  BOOST_LOG(rdInfoLog) <<*i << " ";
  //}
  //BOOST_LOG(rdInfoLog) <<std::endl;
  
  RWMol *m2 = SmilesToMol("CC");
  delete fp2;
  BOOST_LOG(rdInfoLog) <<"--------------------  fp2 " << std::endl;
  fp2=RDKFingerprintMol2(*m2,1,4,2048,4,false);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);
  TEST_ASSERT(OnBitProjSimilarity(*fp2,*fp1)[0]==1.0);

  delete m2;
  m2 = SmilesToMol("CC=O");
  delete fp2;
  fp2=RDKFingerprintMol2(*m2,1,4,2048,4,false);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);
  TEST_ASSERT(OnBitProjSimilarity(*fp2,*fp1)[0]==1.0);

  delete m2;
  m2 = SmilesToMol("CCCOC");
  delete fp2;
  fp2=RDKFingerprintMol2(*m2,1,4,2048,4,false);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)<1.0);
  TEST_ASSERT(TanimotoSimilarity(*fp1,*fp2)>0.0);
  TEST_ASSERT(OnBitProjSimilarity(*fp2,*fp1)[0]==1.0);


  delete m1;
  delete m2;
  delete fp1;
  delete fp2;

  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}


int main(int argc,char *argv[]){
  RDLog::InitLogs();
  test1();
  test2();
  test3();
  test1alg2();
  test2alg2();
  return 0;
}
