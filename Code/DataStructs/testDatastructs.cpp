// $Id$
//
//  Copyright (C) 2001-2014 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ios>
#include "BitVects.h"
#include "BitOps.h"
#include "BitVectUtils.h"
#include "base64.h"
#include <cmath>
#include "DiscreteValueVect.h"
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDBoost/Exceptions.h>
#include <DataStructs/SparseIntVect.h>

#include <stdlib.h>

using namespace std;
using namespace RDKit;
template< typename T >
inline void TXTMSG(const char *__a__,T __b__){  BOOST_LOG(rdInfoLog) << (__a__) << " " << (__b__) << std::endl; }

bool feq(double v1,double v2,double tol=1e-4){
  return fabs(v1-v2)<1e-4;
}

template<typename T> void Test(T arg){
  T t1(20);
  TXTMSG("Set 10:",t1.setBit(10));
  TXTMSG("Set 11:",t1.setBit(11));
  TXTMSG("Set 14:",t1.setBit(14));
  TXTMSG("Set 10:",t1.setBit(10));
  TXTMSG("Get 14:",t1.getBit(14));
  TXTMSG("Num:",t1.getNumBits());
  TXTMSG("NumOn:",t1.getNumOnBits());
  TXTMSG("NumOff:",t1.getNumOffBits());
  TEST_ASSERT(t1==t1);

  IntVect onBits;
  t1.getOnBits(onBits);
  std::copy(onBits.begin(),onBits.end(),std::ostream_iterator<int>(std::cout,", "));
  std::cout << std::endl;

  T t2(t1);
  //onBits = t2.getOnBits();
  TXTMSG("t2[19]:",t2[19]);
  TXTMSG("t2[14]:",t2[14]);

  TEST_ASSERT(t2==t1);
  
  t2 = t1;
  //onBits = t2.getOnBits();
  TXTMSG("t2[19]:",t2[19]);
  t2.unsetBit(14);
  TXTMSG("t2[14]:",t2[14]);
  t2.setBit(15);
  t2.setBit(17);
  TEST_ASSERT(t2!=t1);


  std::cout << "t1: ";
  t1.getOnBits(onBits);
  std::copy(onBits.begin(),onBits.end(),std::ostream_iterator<int>(std::cout,", "));
  std::cout << std::endl;

  std::cout << "t2: ";
  t2.getOnBits(onBits);
  std::copy(onBits.begin(),onBits.end(),std::ostream_iterator<int>(std::cout,", "));
  std::cout << std::endl;

  std::cout << "t1|t2: ";
  T t3=t1|t2;
  t3.getOnBits(onBits);
  std::copy(onBits.begin(),onBits.end(),std::ostream_iterator<int>(std::cout,", "));
  std::cout << std::endl;

  std::cout << "t1&t2: ";
  t3=t1 & t2;
  t3.getOnBits(onBits);
  std::copy(onBits.begin(),onBits.end(),std::ostream_iterator<int>(std::cout,", "));
  std::cout << std::endl;
  
  std::cout << "t1^t2: ";
  t3=t1 ^ t2;
  t3.getOnBits(onBits);
  std::copy(onBits.begin(),onBits.end(),std::ostream_iterator<int>(std::cout,", "));
  std::cout << std::endl;

  std::cout << "~t1: ";
  t3= ~t1;
  t3.getOnBits(onBits);
  std::copy(onBits.begin(),onBits.end(),std::ostream_iterator<int>(std::cout,", "));
  std::cout << std::endl;

  try{
    t3.getBit(4000);
  } catch (IndexErrorException) {
    std::cout << " except " << endl;
  } catch (...) {
    std::cout << " ERROR EXCEPT " << endl;
  }

  T t4(t1.toString());
  TEST_ASSERT(t4==t1);
  TEST_ASSERT(feq(TanimotoSimilarity(t1,t4),1.0));
  
  T *t5 = FoldFingerprint(t1);
  TEST_ASSERT(t5->getNumBits() == t1.getNumBits()/2);
  TEST_ASSERT(t5->getBit(0));
  TEST_ASSERT(t5->getBit(1));
  TEST_ASSERT(t5->getBit(4));
  TEST_ASSERT(!t5->getBit(2));
  TEST_ASSERT(!t5->getBit(3));
  delete t5;

  std::string pkl=t1.toString();
  const char *pkl64=Base64Encode(pkl.c_str(),pkl.size());
  T t6(t1.getNumBits());
  t6.initFromText(pkl64,strlen(pkl64),true);
  delete [] pkl64;
  TEST_ASSERT(t6==t1);
}


template<typename T> void TaniTest(T &arg){
  std::string fps[4] = {
    ".b+HHa.EgU6+ibEIr89.CpX0g8FZiXH+R0+Ps.mr6tg.2",
    ".b7HEa..ccc+gWEIr89.8lV8gOF3aXFFR.+Ps.mZ6lg.2",
    ".H+nHq2EcY09y5EIr9e.8p50h0NgiWGNx4+Hm+Gbslw.2",
    ".1.HHa..cUI6i5E2rO8.Op10d0NoiWGVx.+Hm.Gb6lo.2",
  };
  double dists[] = {
    1.0,0.788991,0.677165,0.686957,
    1.0,0.578125,0.591304,
    1.0,0.732759,
    1.0
  };
  int idx=0;
  for(int i=0;i<4;i++){
    T v1(256);
    FromDaylightString(v1,fps[i]);
    for(int j=i;j<4;j++){
      T v2(256);
      FromDaylightString(v2,fps[j]);
      double tani=TanimotoSimilarity(v1,v2);
      TEST_ASSERT(feq(tani,dists[idx]));
      tani = TverskySimilarity(v1,v2,1.,1.);
      TEST_ASSERT(feq(tani,dists[idx]));
      tani = SimilarityWrapper(v1,v2,TanimotoSimilarity<T,T>);
      TEST_ASSERT(feq(tani,dists[idx]));
      tani = SimilarityWrapper(v1,v2,1.,1.,TverskySimilarity<T,T>);
      TEST_ASSERT(feq(tani,dists[idx]));
      tani = SimilarityWrapper(v1,v2,TanimotoSimilarity<T,T>,true);
      TEST_ASSERT(feq(tani,1.-dists[idx]));
      tani = SimilarityWrapper(v1,v2,1.,1.,TverskySimilarity<T,T>,true);
      TEST_ASSERT(feq(tani,1.-dists[idx]));
      idx++;
    }
  }
}


template<typename T> void ProbeTest(T &arg){
  int sz=1000;
  T t1(sz),t2(sz);
  for(int i=0;i<sz;i+=2){
    t1.setBit(i);
    if(i<3*sz/4) t2.setBit(i);
  }
  std::string pkl=t1.toString();
  TEST_ASSERT(AllProbeBitsMatch(t1,pkl));
  TEST_ASSERT(AllProbeBitsMatch(t2,pkl));
  TEST_ASSERT(AllProbeBitsMatch(t1.toString(),pkl));
  TEST_ASSERT(AllProbeBitsMatch(t2.toString(),pkl));
  TEST_ASSERT(AllProbeBitsMatch(t1.toString().c_str(),pkl.c_str()));
  TEST_ASSERT(AllProbeBitsMatch(t2.toString().c_str(),pkl.c_str()));
  pkl = t2.toString();
  TEST_ASSERT(!AllProbeBitsMatch(t1,pkl));
  TEST_ASSERT(AllProbeBitsMatch(t2,pkl));
  TEST_ASSERT(!AllProbeBitsMatch(t1.toString(),pkl));
  TEST_ASSERT(AllProbeBitsMatch(t2.toString(),pkl));
  TEST_ASSERT(!AllProbeBitsMatch(t1.toString().c_str(),pkl.c_str()));
  TEST_ASSERT(AllProbeBitsMatch(t2.toString().c_str(),pkl.c_str()));
}


void test1DiscreteVect() {
  DiscreteValueVect vect1(DiscreteValueVect::ONEBITVALUE, 30);
  unsigned int i;
  for (i = 0; i < 15; ++i) {
    vect1.setVal(2*i, 1);
  }

  TEST_ASSERT(vect1.getLength() == 30);
  TEST_ASSERT(vect1.getTotalVal() == 15);
  for (i = 0; i < vect1.getLength(); ++i) {
    TEST_ASSERT(vect1.getVal(i) ==  (i+1)%2);
  }
  try {
    vect1.setVal(28,2);
  } catch (ValueErrorException &dexp) {
    BOOST_LOG(rdInfoLog) << "Expected failure: " << dexp.message() << "\n";
  }

  // all these tests should fail if unsigned int changes from being 
  // 32 bits
  DiscreteValueVect vect2(DiscreteValueVect::TWOBITVALUE, 30);
  for (i = 0; i < vect2.getLength(); ++i) {
    vect2.setVal(i, i%4);
  }

  for (i = 0; i < vect2.getLength(); ++i) {
    TEST_ASSERT(vect2.getVal(i) ==  i%4);
  }
  TEST_ASSERT(vect2.getTotalVal() == 43);
  try {
    vect2.setVal(28,10);
  } catch (ValueErrorException &dexp) {
    BOOST_LOG(rdInfoLog) << "Expected failure: " << dexp.message() << "\n";
  }

  DiscreteValueVect vect4(DiscreteValueVect::FOURBITVALUE, 30);
  for (i = 0; i < vect4.getLength(); ++i) {
    vect4.setVal(i, i%16);
  }

  for (i = 0; i < vect4.getLength(); ++i) {
    TEST_ASSERT(vect4.getVal(i) ==  i%16);
  }
  TEST_ASSERT(vect4.getTotalVal() == 211);
  try {
    vect4.setVal(28,16);
  } catch (ValueErrorException &dexp) {
    BOOST_LOG(rdInfoLog) << "Expected failure: " << dexp.message() << "\n";
  }

  DiscreteValueVect vect8(DiscreteValueVect::EIGHTBITVALUE, 32);
  for (i = 0; i < vect8.getLength(); ++i) {
    vect8.setVal(i, i%256);
  }

  for (i = 0; i < vect8.getLength(); ++i) {
    TEST_ASSERT(vect8.getVal(i) ==  i%256);
  }
  TEST_ASSERT(vect8.getTotalVal() == 496);
  try {
    vect8.setVal(28,257);
  } catch (ValueErrorException &dexp) {
    BOOST_LOG(rdInfoLog) << "Expected failure: " << dexp.message() << "\n";
  }

  DiscreteValueVect vect16(DiscreteValueVect::SIXTEENBITVALUE, 300);
  for (i = 0; i < vect16.getLength(); ++i) {
    vect16.setVal(i, i%300);
  }

  for (i = 0; i < vect16.getLength(); ++i) {
    TEST_ASSERT(vect16.getVal(i) ==  i%300);
  }
 
  TEST_ASSERT(vect16.getTotalVal() == 44850);
  vect16.setVal(28,65535);
  try {
    vect16.setVal(28,65536);
  } catch (ValueErrorException &dexp) {
    BOOST_LOG(rdInfoLog) << "Expected failure: " << dexp.message() << "\n";
  }
  
}

void test2DiscreteVectDists() {
  DiscreteValueVect v1(DiscreteValueVect::ONEBITVALUE, 30);
  DiscreteValueVect v2(DiscreteValueVect::ONEBITVALUE, 30);
  unsigned int i;
  for (i = 0; i < 15; ++i) {
    v1.setVal(2*i, 1);
    v2.setVal(2*i, 1);
  }
  TEST_ASSERT(computeL1Norm(v1, v2) == 0);
  for (i = 0; i < 30; ++i) {
    v2.setVal(i, i%2);
  }
  
  TEST_ASSERT(computeL1Norm(v1, v2) == 30);
  
  for (i = 0; i < 30; ++i) {
    if (i%3 == 0) {
      v2.setVal(i, 1);
    } else {
      v2.setVal(i,0);
    }
  }
  
  TEST_ASSERT(computeL1Norm(v1, v2) == 15);

  DiscreteValueVect v21(DiscreteValueVect::TWOBITVALUE, 30);
  DiscreteValueVect v22(DiscreteValueVect::TWOBITVALUE, 30);
  for (i = 0; i < 30; ++i) {
    v21.setVal(i, i%4);
    v22.setVal(i, i%4);
  }
  TEST_ASSERT(computeL1Norm(v21, v22) == 0);
  for (i = 0; i < 30; ++i) {
    v22.setVal(i, (i+1)%4);
  }
  TEST_ASSERT(computeL1Norm(v21, v22) == 44);

  DiscreteValueVect v41(DiscreteValueVect::FOURBITVALUE, 16);
  DiscreteValueVect v42(DiscreteValueVect::FOURBITVALUE, 16);
  for (i = 0; i < 16; ++i) {
    v41.setVal(i, i%16);
    v42.setVal(i, i%16);
  }
  TEST_ASSERT(computeL1Norm(v41, v42) == 0);

  for (i = 0; i < 16; ++i) {
    v42.setVal(i, i%5);
  }
  TEST_ASSERT(computeL1Norm(v41, v42) ==90);

  DiscreteValueVect v43(v42);
  TEST_ASSERT(computeL1Norm(v42, v43) == 0);

  DiscreteValueVect v81(DiscreteValueVect::EIGHTBITVALUE, 5);
  DiscreteValueVect v82(DiscreteValueVect::EIGHTBITVALUE, 5);
  v81.setVal(0, 34); v82.setVal(0, 34);
  v81.setVal(1, 167); v82.setVal(1, 167);
  v81.setVal(2, 3); v82.setVal(2, 3);
  v81.setVal(3, 56); v82.setVal(3, 56);
  v81.setVal(4, 128); v82.setVal(4, 128);
  TEST_ASSERT(computeL1Norm(v81, v82) == 0);

  v82.setVal(0, 14); v82.setVal(1, 67);
  v82.setVal(2, 103); v82.setVal(3, 6);
  v82.setVal(4, 228); 
  TEST_ASSERT(computeL1Norm(v81, v82) == 370);

  DiscreteValueVect v161(DiscreteValueVect::SIXTEENBITVALUE, 3);
  DiscreteValueVect v162(DiscreteValueVect::SIXTEENBITVALUE, 3);
  v161.setVal(0, 2345); v162.setVal(0, 2345);
  v161.setVal(1, 64578); v162.setVal(1, 64578);
  v161.setVal(2, 34); v162.setVal(2, 34);
  TEST_ASSERT(computeL1Norm(v161, v162) == 0);

  v162.setVal(0, 1345);
  v162.setVal(1, 54578);
  v162.setVal(2, 10034);
  TEST_ASSERT(computeL1Norm(v161, v162) == 21000);

}

void test3DiscreteVectPickles() {
  DiscreteValueVect v1(DiscreteValueVect::ONEBITVALUE, 30);
  unsigned int i;
  for (i = 0; i < 15; ++i) {
    v1.setVal(2*i, 1);
  }
  DiscreteValueVect v2(v1.toString());
  TEST_ASSERT(computeL1Norm(v1, v2) == 0);

  DiscreteValueVect v21(DiscreteValueVect::TWOBITVALUE, 30);
  for (i = 0; i < 30; ++i) {
    v21.setVal(i, i%4);
  }
  DiscreteValueVect v22(v21.toString());
  TEST_ASSERT(computeL1Norm(v21, v22) == 0);

  DiscreteValueVect v41(DiscreteValueVect::FOURBITVALUE, 16);
  for (i = 0; i < 16; ++i) {
    v41.setVal(i, i%16);
  }
  DiscreteValueVect v42(v41.toString());
  TEST_ASSERT(computeL1Norm(v41, v42) == 0);

  DiscreteValueVect v81(DiscreteValueVect::EIGHTBITVALUE, 5);
  v81.setVal(0, 34); 
  v81.setVal(1, 167);
  v81.setVal(2, 3); 
  v81.setVal(3, 56);
  v81.setVal(4, 128);
  DiscreteValueVect v82(v81.toString());
  TEST_ASSERT(computeL1Norm(v81, v82) == 0);

  DiscreteValueVect v161(DiscreteValueVect::SIXTEENBITVALUE, 3);
  v161.setVal(0, 2345);
  v161.setVal(1, 64578);
  v161.setVal(2, 34);
  DiscreteValueVect v162(v161.toString());
  TEST_ASSERT(computeL1Norm(v161, v162) == 0);
}


void test4DiscreteVectOps1() {
  DiscreteValueVect vect1(DiscreteValueVect::ONEBITVALUE, 8);
  for (unsigned int i = 0; i < 4; ++i) {
    vect1.setVal(2*i, 1);
  }
  TEST_ASSERT(vect1.getLength() == 8);
  TEST_ASSERT(vect1.getTotalVal() == 4);

  DiscreteValueVect vect2(DiscreteValueVect::ONEBITVALUE, 8);
  for (unsigned int i = 0; i < 4; ++i) {
    vect2.setVal(2*i+1, 1);
  }
  TEST_ASSERT(vect2.getTotalVal() == 4);

  DiscreteValueVect vect3=vect1&vect2;
  TEST_ASSERT(vect3.getLength() == 8);
  TEST_ASSERT(vect3.getTotalVal() == 0);
  
  DiscreteValueVect vect4=vect1|vect2;
  TEST_ASSERT(vect4.getLength() == 8);
  TEST_ASSERT(vect4.getTotalVal() == 8);
#if 0
  DiscreteValueVect vect5=~vect1;
  TEST_ASSERT(vect5.getLength() == 8);
  TEST_ASSERT(vect5.getTotalVal() == 4);

  TEST_ASSERT((vect5&vect1).getTotalVal()==0);
  TEST_ASSERT((vect5&vect2).getTotalVal()==4);
#endif
}

void test5DiscreteVectOps2() {
  DiscreteValueVect vect1(DiscreteValueVect::TWOBITVALUE, 8);
  for (unsigned int i = 0; i < 4; ++i) {
    vect1.setVal(2*i, 2);
  }
  TEST_ASSERT(vect1.getLength() == 8);
  TEST_ASSERT(vect1.getTotalVal() == 8);

  DiscreteValueVect vect2(DiscreteValueVect::TWOBITVALUE, 8);
  for (unsigned int i = 0; i < 4; ++i) {
    vect2.setVal(2*i+1, 2);
    vect2.setVal(2*i, 1);
  }
  TEST_ASSERT(vect2.getTotalVal() == 12);

  DiscreteValueVect vect3=vect1&vect2;
  TEST_ASSERT(vect3.getLength() == 8);
  TEST_ASSERT(vect3.getTotalVal() == 4);
  
  DiscreteValueVect vect4=vect1|vect2;
  TEST_ASSERT(vect4.getLength() == 8);
  TEST_ASSERT(vect4.getTotalVal() == 16);

  DiscreteValueVect vect5=vect1+vect2;
  TEST_ASSERT(vect5.getLength() == 8);
  TEST_ASSERT(vect5.getTotalVal() == 20);

  vect5=vect1-vect2;
  TEST_ASSERT(vect5.getTotalVal()==4);
  vect5=vect2-vect1;
  TEST_ASSERT(vect5.getTotalVal()==8);


#if 0
  DiscreteValueVect vect5=~vect1;
  TEST_ASSERT(vect5.getLength() == 8);
  TEST_ASSERT(vect5.getTotalVal() == 16);

  TEST_ASSERT((vect5&vect1).getTotalVal()==4);
  TEST_ASSERT((vect5&vect2).getTotalVal()==12);
#endif
}

void test6SparseIntVect() {
  SparseIntVect<int> iVect(255);

  TEST_ASSERT(iVect.getLength() == 255);
  TEST_ASSERT(iVect.getVal(23) ==0);
  iVect.setVal(23,14);
  TEST_ASSERT(iVect.getVal(23) ==14);
  
  SparseIntVect<int> oVect(iVect);
  TEST_ASSERT(oVect.getLength() == 255);
  TEST_ASSERT(oVect.getVal(23) ==14);

  std::vector<int> tmpV(3);
  tmpV[0]=1;
  tmpV[1]=5;
  tmpV[2]=1;
  TEST_ASSERT(iVect.getVal(1) ==0);
  TEST_ASSERT(iVect[1] ==0);
  TEST_ASSERT(iVect.getVal(5) ==0);
  TEST_ASSERT(iVect[5] ==0);
  updateFromSequence(iVect,tmpV);
  TEST_ASSERT(iVect.getVal(1) ==2);
  TEST_ASSERT(iVect[1] ==2);
  TEST_ASSERT(iVect.getVal(5) ==1);
  TEST_ASSERT(iVect[5] ==1);

  iVect.setVal(3,-4);
  TEST_ASSERT(iVect.getTotalVal()==13);
  
  try {
    iVect.setVal(-1,13);
    TEST_ASSERT(0);
  } catch (IndexErrorException &dexp) {
    ;
  }
  try {
    iVect.setVal(255,42);
    TEST_ASSERT(0);
  } catch (IndexErrorException &dexp) {
    ;
  }
  try {
    iVect.getVal(-1);
    TEST_ASSERT(0);
  } catch (IndexErrorException &dexp) {
    ;
  }
  try {
    iVect.getVal(255);
    TEST_ASSERT(0);
  } catch (IndexErrorException &dexp) {
    ;
  }
  try {
    iVect[-1];
    TEST_ASSERT(0);
  } catch (IndexErrorException &dexp) {
    ;
  }

  { 
    SparseIntVect<int> iV1(5);
    iV1.setVal(4,4);
    iV1.setVal(0,2);
    iV1.setVal(3,1);
    SparseIntVect<int>::StorageType::const_iterator iter=iV1.getNonzeroElements().begin();
    TEST_ASSERT(iter->first==0);
    TEST_ASSERT(iter->second==2);
    ++iter;
    TEST_ASSERT(iter->first==3);
    TEST_ASSERT(iter->second==1);
    ++iter;
    TEST_ASSERT(iter->first==4);
    TEST_ASSERT(iter->second==4);
    ++iter;
    TEST_ASSERT(iter==iV1.getNonzeroElements().end());
    TEST_ASSERT(feq(DiceSimilarity(iV1,iV1),1.));
  }
  
  { // iV1 &= iV2
    SparseIntVect<int> iV1(5),iV2(5);
    iV1.setVal(0,2);
    iV1.setVal(2,1);
    iV1.setVal(3,4);
    iV1.setVal(4,4);

    iV2.setVal(1,2);
    iV2.setVal(2,3);
    iV2.setVal(3,4);
    iV2.setVal(4,6);

    TEST_ASSERT(feq(DiceSimilarity(iV1,iV2),18./26.));

    iV1 &= iV2;
    TEST_ASSERT(iV1[0]==0);
    TEST_ASSERT(iV1[1]==0);
    TEST_ASSERT(iV1[2]==1);
    TEST_ASSERT(iV1[3]==4);
    TEST_ASSERT(iV1[4]==4);

    TEST_ASSERT(iV2[0]==0);
    TEST_ASSERT(iV2[1]==2);
    TEST_ASSERT(iV2[2]==3);
    TEST_ASSERT(iV2[3]==4);
    TEST_ASSERT(iV2[4]==6);

    TEST_ASSERT(feq(DiceSimilarity(iV1,iV2),18./24.));
    TEST_ASSERT(feq(TverskySimilarity(iV1,iV2,0.5,0.5,false),9./12.));
    TEST_ASSERT(feq(TverskySimilarity(iV1,iV2,1.0,1.0,false),9./15.));
    TEST_ASSERT(feq(TanimotoSimilarity(iV1,iV2),9./15.));
    TEST_ASSERT(feq(TverskySimilarity(iV1,iV2,0.333333333,0.66666666667,false),9./13.));
    TEST_ASSERT(feq(TverskySimilarity(iV1,iV2,1.0,0.0,false),9./9.));

    try {
      iV1 &= iVect;
      TEST_ASSERT(0);
    } catch (ValueErrorException &dexp) {
      ;
    }
  }
  
  { // iV3 = iv1&iV2
    SparseIntVect<int> iV1(5),iV2(5),iV3(5);
    iV1.setVal(0,2);
    iV1.setVal(2,1);
    iV1.setVal(3,4);
    iV1.setVal(4,4);

    iV2.setVal(1,2);
    iV2.setVal(2,3);
    iV2.setVal(3,4);
    iV2.setVal(4,6);

    iV3 = iV1 & iV2;
    TEST_ASSERT(iV3[0]==0);
    TEST_ASSERT(iV3[1]==0);
    TEST_ASSERT(iV3[2]==1);
    TEST_ASSERT(iV3[3]==4);
    TEST_ASSERT(iV3[4]==4);

    TEST_ASSERT(iV1[0]==2);
    TEST_ASSERT(iV1[1]==0);
    TEST_ASSERT(iV1[2]==1);
    TEST_ASSERT(iV1[3]==4);
    TEST_ASSERT(iV1[4]==4);
  
    TEST_ASSERT(iV2[0]==0);
    TEST_ASSERT(iV2[1]==2);
    TEST_ASSERT(iV2[2]==3);
    TEST_ASSERT(iV2[3]==4);
    TEST_ASSERT(iV2[4]==6);
  }
  
  { // iV2 &= iV1
    SparseIntVect<int> iV1(5),iV2(5);
    iV1.setVal(0,2);
    iV1.setVal(2,1);
    iV1.setVal(3,4);
    iV1.setVal(4,4);

    iV2.setVal(1,2);
    iV2.setVal(2,3);
    iV2.setVal(3,4);
    iV2.setVal(4,6);

    iV2 &= iV1;
    TEST_ASSERT(iV2[0]==0);
    TEST_ASSERT(iV2[1]==0);
    TEST_ASSERT(iV2[2]==1);
    TEST_ASSERT(iV2[3]==4);
    TEST_ASSERT(iV2[4]==4);

    TEST_ASSERT(iV1[0]==2);
    TEST_ASSERT(iV1[1]==0);
    TEST_ASSERT(iV1[2]==1);
    TEST_ASSERT(iV1[3]==4);
    TEST_ASSERT(iV1[4]==4);
  
    try {
      iV2 &= iVect;
      TEST_ASSERT(0);
    } catch (ValueErrorException &dexp) {
      ;
    }
  }
  
  {  // iV1 |= iV2
    SparseIntVect<int> iV1(5),iV2(5);
    iV1.setVal(0,2);
    iV1.setVal(2,1);
    iV1.setVal(3,4);
    iV1.setVal(4,4);

    iV2.setVal(1,2);
    iV2.setVal(2,3);
    iV2.setVal(3,4);
    iV2.setVal(4,6);

    iV1 |= iV2;
    TEST_ASSERT(iV1[0]==2);
    TEST_ASSERT(iV1[1]==2);
    TEST_ASSERT(iV1[2]==3);
    TEST_ASSERT(iV1[3]==4);
    TEST_ASSERT(iV1[4]==6);

    TEST_ASSERT(iV2[0]==0);
    TEST_ASSERT(iV2[1]==2);
    TEST_ASSERT(iV2[2]==3);
    TEST_ASSERT(iV2[3]==4);
    TEST_ASSERT(iV2[4]==6);
  
    try {
      iV1 |= iVect;
      TEST_ASSERT(0);
    } catch (ValueErrorException &dexp) {
      ;
    }
  }
  
  { // iV3 = iv1 |iV2
    SparseIntVect<int> iV1(5),iV2(5),iV3(5);
    iV1.setVal(0,2);
    iV1.setVal(2,1);
    iV1.setVal(3,4);
    iV1.setVal(4,4);

    iV2.setVal(1,2);
    iV2.setVal(2,3);
    iV2.setVal(3,4);
    iV2.setVal(4,6);

    iV3 = iV1 | iV2;
    TEST_ASSERT(iV3[0]==2);
    TEST_ASSERT(iV3[1]==2);
    TEST_ASSERT(iV3[2]==3);
    TEST_ASSERT(iV3[3]==4);
    TEST_ASSERT(iV3[4]==6);

    TEST_ASSERT(iV1[0]==2);
    TEST_ASSERT(iV1[1]==0);
    TEST_ASSERT(iV1[2]==1);
    TEST_ASSERT(iV1[3]==4);
    TEST_ASSERT(iV1[4]==4);
  
    TEST_ASSERT(iV2[0]==0);
    TEST_ASSERT(iV2[1]==2);
    TEST_ASSERT(iV2[2]==3);
    TEST_ASSERT(iV2[3]==4);
    TEST_ASSERT(iV2[4]==6);
  }
  
  {  // iV2 |= iV1
    SparseIntVect<int> iV1(5),iV2(5);
    iV1.setVal(0,2);
    iV1.setVal(2,1);
    iV1.setVal(3,4);
    iV1.setVal(4,4);

    iV2.setVal(1,2);
    iV2.setVal(2,3);
    iV2.setVal(3,4);
    iV2.setVal(4,6);

    iV2 |= iV1;
    TEST_ASSERT(iV2[0]==2);
    TEST_ASSERT(iV2[1]==2);
    TEST_ASSERT(iV2[2]==3);
    TEST_ASSERT(iV2[3]==4);
    TEST_ASSERT(iV2[4]==6);

    TEST_ASSERT(iV1[0]==2);
    TEST_ASSERT(iV1[1]==0);
    TEST_ASSERT(iV1[2]==1);
    TEST_ASSERT(iV1[3]==4);
    TEST_ASSERT(iV1[4]==4);
  }
  
  {  // iV1 += iV2
    SparseIntVect<int> iV1(5),iV2(5);
    iV1.setVal(0,2);
    iV1.setVal(2,1);
    iV1.setVal(3,4);
    iV1.setVal(4,4);

    iV2.setVal(1,2);
    iV2.setVal(2,3);
    iV2.setVal(3,-4);
    iV2.setVal(4,6);

    iV1 += iV2;
    TEST_ASSERT(iV1[0]==2);
    TEST_ASSERT(iV1[1]==2);
    TEST_ASSERT(iV1[2]==4);
    TEST_ASSERT(iV1[3]==0);
    TEST_ASSERT(iV1[4]==10);

    TEST_ASSERT(iV2[0]==0);
    TEST_ASSERT(iV2[1]==2);
    TEST_ASSERT(iV2[2]==3);
    TEST_ASSERT(iV2[3]==-4);
    TEST_ASSERT(iV2[4]==6);
  }
  {  // iV3 = IV1 + iV2
    SparseIntVect<int> iV1(5),iV2(5),iV3(5);
    iV1.setVal(2,1);
    iV1.setVal(0,2);
    iV1.setVal(4,4);
    iV1.setVal(3,4);

    iV2.setVal(1,2);
    iV2.setVal(2,3);
    iV2.setVal(3,-4);
    iV2.setVal(4,6);

    iV3=iV1+iV2;
    TEST_ASSERT(iV3[0]==2);
    TEST_ASSERT(iV3[1]==2);
    TEST_ASSERT(iV3[2]==4);
    TEST_ASSERT(iV3[3]==0);
    TEST_ASSERT(iV3[4]==10);

    TEST_ASSERT(iV1[0]==2);
    TEST_ASSERT(iV1[1]==0);
    TEST_ASSERT(iV1[2]==1);
    TEST_ASSERT(iV1[3]==4);
    TEST_ASSERT(iV1[4]==4);
    TEST_ASSERT(iV2[0]==0);
    TEST_ASSERT(iV2[1]==2);
    TEST_ASSERT(iV2[2]==3);
    TEST_ASSERT(iV2[3]==-4);
    TEST_ASSERT(iV2[4]==6);
  }

  {  // iV1 -= iV2
    SparseIntVect<int> iV1(5),iV2(5);
    iV1.setVal(0,2);
    iV1.setVal(2,1);
    iV1.setVal(3,4);
    iV1.setVal(4,4);

    iV2.setVal(1,2);
    iV2.setVal(2,3);
    iV2.setVal(3,4);
    iV2.setVal(4,6);

    iV1 -= iV2;
    TEST_ASSERT(iV1[0]==2);
    TEST_ASSERT(iV1[1]==-2);
    TEST_ASSERT(iV1[2]==-2);
    TEST_ASSERT(iV1[3]==0);
    TEST_ASSERT(iV1[4]==-2);

    TEST_ASSERT(iV2[0]==0);
    TEST_ASSERT(iV2[2]==3);
    TEST_ASSERT(iV2[1]==2);
    TEST_ASSERT(iV2[4]==6);
    TEST_ASSERT(iV2[3]==4);
  }
  {  // iV3 = IV1 - iV2
    SparseIntVect<int> iV1(5),iV2(5),iV3(5);
    iV1.setVal(0,2);
    iV1.setVal(2,1);
    iV1.setVal(3,4);
    iV1.setVal(4,4);

    iV2.setVal(3,4);
    iV2.setVal(4,6);
    iV2.setVal(1,2);
    iV2.setVal(2,3);

    iV3=iV1-iV2;
    TEST_ASSERT(iV3[0]==2);
    TEST_ASSERT(iV3[1]==-2);
    TEST_ASSERT(iV3[2]==-2);
    TEST_ASSERT(iV3[3]==0);
    TEST_ASSERT(iV3[4]==-2);

    TEST_ASSERT(iV1[0]==2);
    TEST_ASSERT(iV1[1]==0);
    TEST_ASSERT(iV1[2]==1);
    TEST_ASSERT(iV1[3]==4);
    TEST_ASSERT(iV1[4]==4);
    TEST_ASSERT(iV2[0]==0);
    TEST_ASSERT(iV2[1]==2);
    TEST_ASSERT(iV2[2]==3);
    TEST_ASSERT(iV2[3]==4);
    TEST_ASSERT(iV2[4]==6);
  }
  
  {  // operator== and operator!=
    SparseIntVect<int> iV1(5),iV2(5),iV3(3);
    iV1.setVal(0,2);
    iV1.setVal(2,1);
    iV1.setVal(3,4);
    iV1.setVal(4,4);

    iV2.setVal(1,2);
    iV2.setVal(2,3);
    iV2.setVal(3,4);
    iV2.setVal(4,6);

    TEST_ASSERT(iV1==iV1);
    TEST_ASSERT(iV2==iV2);
    TEST_ASSERT(iV3==iV3);
    TEST_ASSERT(iV1!=iV2);
    TEST_ASSERT(iV1!=iV3);
    TEST_ASSERT(iV2!=iV1);
    TEST_ASSERT(iV3!=iV1);
    TEST_ASSERT(iV1!=iV3);
  }

  { //test negative values (was sf.net Issue 3295215)
    SparseIntVect<int> iV1(5),iV2(5);
    iV1.setVal(0,-2);
    iV1.setVal(2,1);
    iV1.setVal(3,-4);
    iV1.setVal(4,4);

    iV2.setVal(1,-2);
    iV2.setVal(2,3);
    iV2.setVal(3,-4);
    iV2.setVal(4,6);

    TEST_ASSERT(feq(DiceSimilarity(iV1,iV2),18./26.));
  }

}

void test7SparseIntVectPickles() {
  {
    SparseIntVect<int> iV1(5),iV2(3);
    iV1.setVal(0,2);
    iV1.setVal(2,1);
    iV1.setVal(3,4);
    iV1.setVal(4,4);

    iV2.setVal(2,3);
    TEST_ASSERT(iV1!=iV2);
    std::string pkl;
    pkl = iV1.toString();
    iV2.fromString(pkl);
    TEST_ASSERT(iV1==iV2);
  }
  
  {
    SparseIntVect<char> iV1(5);
    SparseIntVect<int>iV2(3);
    iV1.setVal(0,2);
    iV1.setVal(2,1);
    iV1.setVal(3,4);

    iV2.setVal(1,1);
    std::string pkl;
    pkl = iV1.toString();
    iV2.fromString(pkl);
    TEST_ASSERT(iV2.getLength()==iV1.getLength());
    TEST_ASSERT(iV2[0]==2)
    TEST_ASSERT(iV2[1]==0)
    TEST_ASSERT(iV2[2]==1)
    TEST_ASSERT(iV2[3]==4)
  }
  
  {
    SparseIntVect<int> iV1(5);
    SparseIntVect<char>iV2(3);
    iV1.setVal(0,2);
    iV1.setVal(2,1);
    iV1.setVal(3,4);

    std::string pkl;
    pkl = iV1.toString();
    try{
      iV2.fromString(pkl);
      TEST_ASSERT(0);
    } catch (ValueErrorException &dexp) {
      ;
    }
  }
}


void test8BitVectPickles() {
#if 0
  {
    std::string dirName = getenv("RDBASE");
    dirName+="/Code/DataStructs/testData/";
    std::string pklName = dirName+"test1.bin";
    std::ofstream outS;
    outS.open(pklName.c_str(),std::ios_base::binary);

    ExplicitBitVect bv(32);
    for(int i=0;i<32;i+=2){
      bv.setBit(i);
    }
    std::string pkl=bv.toString();
    unsigned int sz=pkl.size();
    outS<<sz;
    outS<<pkl;
    outS.close();
  }
#endif


  {
    std::string dirName = getenv("RDBASE");
    dirName+="/Code/DataStructs/testData/";
    std::string pklName = dirName+"test1.bin";
    std::ifstream inS;
    inS.open(pklName.c_str(),std::ios_base::binary);
    unsigned int length;
    inS >> length;
    char *buff = new char[length];
    length=inS.readsome(buff,length);
    inS.close();
    std::string pkl(buff,length);
    delete [] buff;
    ExplicitBitVect bv(pkl);

    TEST_ASSERT(bv.getNumBits()==32);
    TEST_ASSERT(bv.getNumOnBits()==16);
    TEST_ASSERT(bv[0]);
    TEST_ASSERT(!bv[1]);
  }
}

void test9BitVectFPS() {
  {
    ExplicitBitVect bv(32);
    std::string fps;

    fps = BitVectToFPSText(bv);
    TEST_ASSERT(fps=="00000000");
    
    bv.setBit(0);
    bv.setBit(1);
    bv.setBit(17);
    bv.setBit(23);
    bv.setBit(31);

    fps = BitVectToFPSText(bv);
    TEST_ASSERT(fps=="03008280");
  }
  {
    ExplicitBitVect bv(32),bv2(32);
    std::string fps;

    fps = BitVectToFPSText(bv);
    TEST_ASSERT(fps=="00000000");
    UpdateBitVectFromFPSText(bv2,fps);
    TEST_ASSERT(bv==bv2);
    
    bv.setBit(0);
    bv.setBit(1);
    bv.setBit(4);
    bv.setBit(17);
    bv.setBit(23);
    bv.setBit(31);

    fps = BitVectToFPSText(bv);
    UpdateBitVectFromFPSText(bv2,fps);
    TEST_ASSERT(bv==bv2);
  }
  {
    ExplicitBitVect bv(33);
    std::string fps;

    fps = BitVectToFPSText(bv);
    TEST_ASSERT(fps=="0000000000");
    
    bv.setBit(0);
    bv.setBit(32);

    fps = BitVectToFPSText(bv);
    TEST_ASSERT(fps=="0100000001");
  }
  {
    ExplicitBitVect bv(33),bv2(33);
    std::string fps;

    fps = BitVectToFPSText(bv);
    TEST_ASSERT(fps=="0000000000");
    UpdateBitVectFromFPSText(bv2,fps);
    TEST_ASSERT(bv==bv2);
    
    bv.setBit(0);
    bv.setBit(1);
    bv.setBit(4);
    bv.setBit(17);
    bv.setBit(23);
    bv.setBit(32);

    fps = BitVectToFPSText(bv);
    UpdateBitVectFromFPSText(bv2,fps);
    TEST_ASSERT(bv==bv2);
  }
}

void test10BitVectBinaryText() {
  {
    ExplicitBitVect bv(32);
    std::string fps;

    fps = BitVectToBinaryText(bv);
    TEST_ASSERT(fps.size()==4);
    for(unsigned int i=0;i<fps.size();++i){
      TEST_ASSERT(fps[i]==0);
    }
    
    bv.setBit(0);
    bv.setBit(9);
    bv.setBit(17);
    bv.setBit(26);

    fps = BitVectToBinaryText(bv);
    TEST_ASSERT(fps.size()==4);
    for(unsigned int i=0;i<fps.size();++i){
      TEST_ASSERT(fps[i]!=0);
    }
  }
  {
    ExplicitBitVect bv(32),bv2(32);
    std::string fps;

    fps = BitVectToBinaryText(bv);
    TEST_ASSERT(fps.size()==4);
    for(unsigned int i=0;i<fps.size();++i){
      TEST_ASSERT(fps[i]==0);
    }
    UpdateBitVectFromBinaryText(bv2,fps);
    TEST_ASSERT(bv==bv2);
    
    bv.setBit(0);
    bv.setBit(1);
    bv.setBit(4);
    bv.setBit(17);
    bv.setBit(23);
    bv.setBit(31);

    fps = BitVectToBinaryText(bv);
    UpdateBitVectFromBinaryText(bv2,fps);
    TEST_ASSERT(bv==bv2);
  }
  {
    ExplicitBitVect bv(33);
    std::string fps;

    fps = BitVectToBinaryText(bv);
    TEST_ASSERT(fps.size()==5);
    for(unsigned int i=0;i<fps.size();++i){
      TEST_ASSERT(fps[i]==0);
    }
    
    bv.setBit(0);
    bv.setBit(32);

    fps = BitVectToBinaryText(bv);
    TEST_ASSERT(fps.size()==5);
    TEST_ASSERT(fps[0]!=0);
    for(unsigned int i=1;i<fps.size()-1;++i){
      TEST_ASSERT(fps[i]==0);
    }
    TEST_ASSERT(fps[fps.size()-1]!=0);
  }
  {
    ExplicitBitVect bv(33),bv2(33);
    std::string fps;

    fps = BitVectToBinaryText(bv);
    TEST_ASSERT(fps.size()==5);
    UpdateBitVectFromBinaryText(bv2,fps);
    TEST_ASSERT(bv==bv2);
    
    bv.setBit(0);
    bv.setBit(1);
    bv.setBit(4);
    bv.setBit(17);
    bv.setBit(23);
    bv.setBit(32);

    fps = BitVectToBinaryText(bv);
    UpdateBitVectFromBinaryText(bv2,fps);
    TEST_ASSERT(bv==bv2);
  }
}

void test11SimilaritiesBV() {
    // similarity = 1.0
	ExplicitBitVect bv(10);
	bv.setBit(0);
	bv.setBit(1);
	bv.setBit(4);
	bv.setBit(6);
	bv.setBit(9);
	ExplicitBitVect bv2 = bv;

	TEST_ASSERT(feq(TanimotoSimilarity(bv,bv2),1.));
	TEST_ASSERT(feq(DiceSimilarity(bv,bv2),1.));
	TEST_ASSERT(feq(CosineSimilarity(bv,bv2),1.));
	TEST_ASSERT(feq(KulczynskiSimilarity(bv,bv2),1.));
	TEST_ASSERT(feq(SokalSimilarity(bv,bv2),1.));
	TEST_ASSERT(feq(McConnaugheySimilarity(bv,bv2),1.));
	TEST_ASSERT(feq(BraunBlanquetSimilarity(bv,bv2),1.));
	TEST_ASSERT(feq(RusselSimilarity(bv,bv2),0.5));
	TEST_ASSERT(feq(RogotGoldbergSimilarity(bv,bv2),1.));
	TEST_ASSERT(feq(AllBitSimilarity(bv,bv2),1.));

	// similarity = 0.0
	bv2 = ExplicitBitVect(10);
	bv2.setBit(2);
	bv2.setBit(3);
	bv2.setBit(5);
	bv2.setBit(7);
	bv2.setBit(8);

	TEST_ASSERT(feq(TanimotoSimilarity(bv,bv2),0.));
	TEST_ASSERT(feq(DiceSimilarity(bv,bv2),0.));
	TEST_ASSERT(feq(CosineSimilarity(bv,bv2),0.));
	TEST_ASSERT(feq(KulczynskiSimilarity(bv,bv2),0.));
	TEST_ASSERT(feq(SokalSimilarity(bv,bv2),0.));
	TEST_ASSERT(feq(McConnaugheySimilarity(bv,bv2),-1.));
	TEST_ASSERT(feq(BraunBlanquetSimilarity(bv,bv2),0.));
	TEST_ASSERT(feq(RusselSimilarity(bv,bv2),0.));
	TEST_ASSERT(feq(RogotGoldbergSimilarity(bv,bv2),0.));
	TEST_ASSERT(feq(AllBitSimilarity(bv,bv2),0.));

	// similarity ~= 0.5
	bv.setBit(5);
	bv2 = ExplicitBitVect(10);
	bv2.setBit(0);
	bv2.setBit(2);
	bv2.setBit(4);
	bv2.setBit(5);
	bv2.setBit(8);
	bv2.setBit(9);

	TEST_ASSERT(feq(TanimotoSimilarity(bv,bv2),0.5));
	TEST_ASSERT(feq(DiceSimilarity(bv,bv2),0.6666));
	TEST_ASSERT(feq(CosineSimilarity(bv,bv2),0.6666));
	TEST_ASSERT(feq(KulczynskiSimilarity(bv,bv2),0.6666));
	TEST_ASSERT(feq(SokalSimilarity(bv,bv2),0.3333));
	TEST_ASSERT(feq(McConnaugheySimilarity(bv,bv2),0.3333));
	TEST_ASSERT(feq(BraunBlanquetSimilarity(bv,bv2),0.6666));
	TEST_ASSERT(feq(RusselSimilarity(bv,bv2),0.4));
	TEST_ASSERT(feq(RogotGoldbergSimilarity(bv,bv2),0.5833));
	TEST_ASSERT(feq(AllBitSimilarity(bv,bv2),0.6));
}

void test12SimilaritiesSparseBV() {
    // similarity = 1.0
	SparseBitVect sbv(10);
	sbv.setBit(0);
	sbv.setBit(1);
	sbv.setBit(4);
	sbv.setBit(6);
	sbv.setBit(9);
	SparseBitVect sbv2 = sbv;

	TEST_ASSERT(feq(TanimotoSimilarity(sbv,sbv2),1.));
	TEST_ASSERT(feq(DiceSimilarity(sbv,sbv2),1.));
	TEST_ASSERT(feq(CosineSimilarity(sbv,sbv2),1.));
	TEST_ASSERT(feq(KulczynskiSimilarity(sbv,sbv2),1.));
	TEST_ASSERT(feq(SokalSimilarity(sbv,sbv2),1.));
	TEST_ASSERT(feq(McConnaugheySimilarity(sbv,sbv2),1.));
	TEST_ASSERT(feq(BraunBlanquetSimilarity(sbv,sbv2),1.));
	TEST_ASSERT(feq(RusselSimilarity(sbv,sbv2),0.5));
	TEST_ASSERT(feq(RogotGoldbergSimilarity(sbv,sbv2),1.));
	TEST_ASSERT(feq(AllBitSimilarity(sbv,sbv2),1.));

	// similarity = 0.0
	sbv2 = SparseBitVect(10);
	sbv2.setBit(2);
	sbv2.setBit(3);
	sbv2.setBit(5);
	sbv2.setBit(7);
	sbv2.setBit(8);

	TEST_ASSERT(feq(TanimotoSimilarity(sbv,sbv2),0.));
	TEST_ASSERT(feq(DiceSimilarity(sbv,sbv2),0.));
	TEST_ASSERT(feq(CosineSimilarity(sbv,sbv2),0.));
	TEST_ASSERT(feq(KulczynskiSimilarity(sbv,sbv2),0.));
	TEST_ASSERT(feq(SokalSimilarity(sbv,sbv2),0.));
	TEST_ASSERT(feq(McConnaugheySimilarity(sbv,sbv2),-1.));
	TEST_ASSERT(feq(BraunBlanquetSimilarity(sbv,sbv2),0.));
	TEST_ASSERT(feq(RusselSimilarity(sbv,sbv2),0.));
	TEST_ASSERT(feq(RogotGoldbergSimilarity(sbv,sbv2),0.));
	TEST_ASSERT(feq(AllBitSimilarity(sbv,sbv2),0.));

	// similarity ~= 0.5
	sbv.setBit(5);
	sbv2 = SparseBitVect(10);
	sbv2.setBit(0);
	sbv2.setBit(2);
	sbv2.setBit(4);
	sbv2.setBit(5);
	sbv2.setBit(8);
	sbv2.setBit(9);

	TEST_ASSERT(feq(TanimotoSimilarity(sbv,sbv2),0.5));
	TEST_ASSERT(feq(DiceSimilarity(sbv,sbv2),0.6666));
	TEST_ASSERT(feq(CosineSimilarity(sbv,sbv2),0.6666));
	TEST_ASSERT(feq(KulczynskiSimilarity(sbv,sbv2),0.6666));
	TEST_ASSERT(feq(SokalSimilarity(sbv,sbv2),0.3333));
	TEST_ASSERT(feq(McConnaugheySimilarity(sbv,sbv2),0.3333));
	TEST_ASSERT(feq(BraunBlanquetSimilarity(sbv,sbv2),0.6666));
	TEST_ASSERT(feq(RusselSimilarity(sbv,sbv2),0.4));
	TEST_ASSERT(feq(RogotGoldbergSimilarity(sbv,sbv2),0.5833));
	TEST_ASSERT(feq(AllBitSimilarity(sbv,sbv2),0.6));
}

void test13BitVectAllOnes() {
  {
    ExplicitBitVect bv(32, false);
    TEST_ASSERT(bv.getNumOnBits()==0);
    TEST_ASSERT(bv[0]==0);

    ExplicitBitVect bv2(32, true);
    TEST_ASSERT(bv2.getNumOnBits()==32);
    TEST_ASSERT(bv2[0]==1);
  }
}

void test18BitVectConcatenation() {
 {
   ExplicitBitVect bv(32, false);
   ExplicitBitVect bv2(32, true);
   ExplicitBitVect bv3 = bv + bv2;
   TEST_ASSERT(bv3.getNumBits() == 64);
   TEST_ASSERT(bv3.getNumOnBits() == 32);
   TEST_ASSERT(bv3.getNumOffBits() == 32);
  }
}

int main(){
  RDLog::InitLogs();
  try{
    throw IndexErrorException(3);
  } catch (IndexErrorException) {
    BOOST_LOG(rdInfoLog) << "pass" << endl;
  }

  stringstream ss(ios_base::binary|ios_base::out|ios_base::in);
  int v1=4,v2=5,v3,v4;

  ss.write((const char *)&v1,sizeof(v1));
  ss.write((const char *)&v2,sizeof(v2));
  ss.seekp(0,ios_base::beg);
  RDKit::streamRead(ss,v3);
  RDKit::streamRead(ss,v4);
  
  TXTMSG("v3",v3);
  TXTMSG("v4",v4);
  
  BOOST_LOG(rdInfoLog) << " SPARSE -----------------------------------" << std::endl;
  SparseBitVect sparseFoo(10);
  Test(sparseFoo);
  TaniTest(sparseFoo);
  ProbeTest(sparseFoo);
  BOOST_LOG(rdInfoLog) << " Explicit ----------------------------------" << std::endl;
  ExplicitBitVect explicitFoo(10);
  Test(explicitFoo);
  TaniTest(explicitFoo);
  BOOST_LOG(rdInfoLog) << " Done" << std::endl;
  
  BOOST_LOG(rdInfoLog) << " Test DiscreteValue Vectors 1 ----------------------------" << endl;
  test1DiscreteVect();
  BOOST_LOG(rdInfoLog) << " Test DiscreteValue Vectors 2 ------------------------------" << endl;
  test2DiscreteVectDists();
  BOOST_LOG(rdInfoLog) << " Test DiscreteValue Vectors 3 ---------------------------" << endl;
  test3DiscreteVectPickles();

  BOOST_LOG(rdInfoLog) << " Test DiscreteValue Operations -----------------------------" << endl;
  test4DiscreteVectOps1();

  BOOST_LOG(rdInfoLog) << " Test DiscreteValue Operations 2 -------------------------- "<< endl;
  test5DiscreteVectOps2();

  BOOST_LOG(rdInfoLog) << " Test SparseIntVect  ------------------------------------" << std::endl;
  test6SparseIntVect();
  BOOST_LOG(rdInfoLog) << " Test SparseIntVect Serialization  --------------------------" << std::endl;
  test7SparseIntVectPickles();

  BOOST_LOG(rdInfoLog) << " Test BitVect Serialization  -------------------------------" << std::endl;
  test8BitVectPickles();
  
  BOOST_LOG(rdInfoLog) << " Test BitVect to FPS  -------------------------------" << std::endl;
  test9BitVectFPS();
  
  BOOST_LOG(rdInfoLog) << " Test BitVect to binary string  -------------------------------" << std::endl;
  test10BitVectBinaryText();
  
  BOOST_LOG(rdInfoLog) << " Test Similarity Measures BitVect -------------------------------" << std::endl;
  test11SimilaritiesBV();

  BOOST_LOG(rdInfoLog) << " Test Similarity Measures SparseBitVect -------------------------------" << std::endl;
    test12SimilaritiesSparseBV();

  BOOST_LOG(rdInfoLog) << " Test BitVect with all ones -------------------------------" << std::endl;
  test13BitVectAllOnes();

  BOOST_LOG(rdInfoLog) << " Test Explicit BitVects: Concatenation Operation  -------------------------------" << std::endl;
  test18BitVectConcatenation();

  return 0;
  
}
