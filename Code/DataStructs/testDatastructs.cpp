// $Id$
//
//  Copyright (C) 2001-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#include <iostream>
#include <sstream>
#include <fstream>
#include <ios>
#include "BitVects.h"
#include "BitOps.h"
#include "BitVectUtils.h"
#include <cmath>
#include "DiscreteValueVect.h"
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDBoost/Exceptions.h>


using namespace std;
using namespace RDKit;
template< typename T >
inline void TXTMSG(const char *__a__,T __b__){  BOOST_LOG(rdInfoLog) << (__a__) << " " << (__b__) << std::endl; }

template<typename T> void Test(T arg){
  T t1(20);
  TXTMSG("Set 10:",t1.SetBit(10));
  TXTMSG("Set 11:",t1.SetBit(11));
  TXTMSG("Set 14:",t1.SetBit(14));
  TXTMSG("Set 10:",t1.SetBit(10));
  TXTMSG("Get 14:",t1.GetBit(14));
  TXTMSG("Num:",t1.GetNumBits());
  TXTMSG("NumOn:",t1.GetNumOnBits());
  TXTMSG("NumOff:",t1.GetNumOffBits());

  IntVect onBits;
  t1.GetOnBits(onBits);
  std::copy(onBits.begin(),onBits.end(),std::ostream_iterator<int>(std::cout,", "));
  std::cout << std::endl;

  T t2(t1);
  //onBits = t2.GetOnBits();
  TXTMSG("t2[19]:",t2[19]);
  TXTMSG("t2[14]:",t2[14]);

  t2 = t1;
  //onBits = t2.GetOnBits();
  TXTMSG("t2[19]:",t2[19]);
  t2.UnSetBit(14);
  TXTMSG("t2[14]:",t2[14]);
  t2.SetBit(15);
  t2.SetBit(17);


  cout << "t1: ";
  t1.GetOnBits(onBits);
  std::copy(onBits.begin(),onBits.end(),std::ostream_iterator<int>(std::cout,", "));
  std::cout << std::endl;

  cout << "t2: ";
  t2.GetOnBits(onBits);
  std::copy(onBits.begin(),onBits.end(),std::ostream_iterator<int>(std::cout,", "));
  std::cout << std::endl;

  cout << "t1|t2: ";
  T t3=t1|t2;
  t3.GetOnBits(onBits);
  std::copy(onBits.begin(),onBits.end(),std::ostream_iterator<int>(std::cout,", "));
  std::cout << std::endl;

  cout << "t1&t2: ";
  t3=t1 & t2;
  t3.GetOnBits(onBits);
  std::copy(onBits.begin(),onBits.end(),std::ostream_iterator<int>(std::cout,", "));
  std::cout << std::endl;
  
  cout << "t1^t2: ";
  t3=t1 ^ t2;
  t3.GetOnBits(onBits);
  std::copy(onBits.begin(),onBits.end(),std::ostream_iterator<int>(std::cout,", "));
  std::cout << std::endl;

  cout << "~t1: ";
  t3= ~t1;
  t3.GetOnBits(onBits);
  std::copy(onBits.begin(),onBits.end(),std::ostream_iterator<int>(std::cout,", "));
  std::cout << std::endl;

  try{
    t3.GetBit(4000);
  } catch (IndexErrorException) {
    cout << " except " << endl;
  } catch (...) {
    cout << " ERROR EXCEPT " << endl;
  }

  T t4(t1.ToString());
  cout << "tan(t1,t4): " << TanimotoSimilarity(t1,t4) << endl;
  
  T *t5 = FoldFingerprint(t1);
  TEST_ASSERT(t5->GetNumBits() == t1.GetNumBits()/2);
  TEST_ASSERT(t5->GetBit(0));
  TEST_ASSERT(t5->GetBit(1));
  TEST_ASSERT(t5->GetBit(4));
  TEST_ASSERT(!t5->GetBit(2));
  TEST_ASSERT(!t5->GetBit(3));
  delete t5;
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
      TEST_ASSERT(abs(tani-dists[idx])<1e-4);
      idx++;
    }
  }
}


template<typename T> void ProbeTest(T &arg){
  int sz=1000;
  T t1(sz),t2(sz);
  for(int i=0;i<sz;i+=2){
    t1.SetBit(i);
    if(i<3*sz/4) t2.SetBit(i);
  }
  std::string pkl=t1.ToString();
  TEST_ASSERT(AllProbeBitsMatch(t1,pkl));
  TEST_ASSERT(AllProbeBitsMatch(t2,pkl));
  TEST_ASSERT(AllProbeBitsMatch(t1.ToString(),pkl));
  TEST_ASSERT(AllProbeBitsMatch(t2.ToString(),pkl));
  TEST_ASSERT(AllProbeBitsMatch(t1.ToString().c_str(),pkl.c_str()));
  TEST_ASSERT(AllProbeBitsMatch(t2.ToString().c_str(),pkl.c_str()));
  pkl = t2.ToString();
  TEST_ASSERT(!AllProbeBitsMatch(t1,pkl));
  TEST_ASSERT(AllProbeBitsMatch(t2,pkl));
  TEST_ASSERT(!AllProbeBitsMatch(t1.ToString(),pkl));
  TEST_ASSERT(AllProbeBitsMatch(t2.ToString(),pkl));
  TEST_ASSERT(!AllProbeBitsMatch(t1.ToString().c_str(),pkl.c_str()));
  TEST_ASSERT(AllProbeBitsMatch(t2.ToString().c_str(),pkl.c_str()));
}


void test1DiscreteVect() {
  DiscreteValueVect vect1(DiscreteValueVect::ONEBITVALUE, 30);
  unsigned int i;
  for (i = 0; i < 15; ++i) {
    vect1.setVal(2*i, 1);
  }

  CHECK_INVARIANT(vect1.getLength() == 30, "");
  CHECK_INVARIANT(vect1.getTotalVal() == 15, "");
  for (i = 0; i < vect1.getLength(); ++i) {
    CHECK_INVARIANT(vect1.getVal(i) ==  (i+1)%2, "");
  }
  try {
    vect1.setVal(28,2);
  } catch (ValueErrorException &dexp) {
    std::cout << "Expected failure: " << dexp.message() << "\n";
  }

  // all these tests should fail if unsigned int changes from being 
  // 32 bits
  DiscreteValueVect vect2(DiscreteValueVect::TWOBITVALUE, 30);
  for (i = 0; i < vect2.getLength(); ++i) {
    vect2.setVal(i, i%4);
  }

  for (i = 0; i < vect2.getLength(); ++i) {
    CHECK_INVARIANT(vect2.getVal(i) ==  i%4, "");
  }
  CHECK_INVARIANT(vect2.getTotalVal() == 43, "");
  try {
    vect2.setVal(28,10);
  } catch (ValueErrorException &dexp) {
    std::cout << "Expected failure: " << dexp.message() << "\n";
  }

  DiscreteValueVect vect4(DiscreteValueVect::FOURBITVALUE, 30);
  for (i = 0; i < vect4.getLength(); ++i) {
    vect4.setVal(i, i%16);
  }

  for (i = 0; i < vect4.getLength(); ++i) {
    CHECK_INVARIANT(vect4.getVal(i) ==  i%16, "");
  }
  CHECK_INVARIANT(vect4.getTotalVal() == 211, "");
  try {
    vect4.setVal(28,16);
  } catch (ValueErrorException &dexp) {
    std::cout << "Expected failure: " << dexp.message() << "\n";
  }

  DiscreteValueVect vect8(DiscreteValueVect::EIGHTBITVALUE, 32);
  for (i = 0; i < vect8.getLength(); ++i) {
    vect8.setVal(i, i%256);
  }

  for (i = 0; i < vect8.getLength(); ++i) {
    CHECK_INVARIANT(vect8.getVal(i) ==  i%256, "");
  }
  CHECK_INVARIANT(vect8.getTotalVal() == 496, "");
  try {
    vect8.setVal(28,257);
  } catch (ValueErrorException &dexp) {
    std::cout << "Expected failure: " << dexp.message() << "\n";
  }

  DiscreteValueVect vect16(DiscreteValueVect::SIXTEENBITVALUE, 300);
  for (i = 0; i < vect16.getLength(); ++i) {
    vect16.setVal(i, i%300);
  }

  for (i = 0; i < vect16.getLength(); ++i) {
    CHECK_INVARIANT(vect16.getVal(i) ==  i%300, "");
  }
 
  CHECK_INVARIANT(vect16.getTotalVal() == 44850, "");
  vect16.setVal(28,65535);
  try {
    vect16.setVal(28,65536);
  } catch (ValueErrorException &dexp) {
    std::cout << "Expected failure: " << dexp.message() << "\n";
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
  CHECK_INVARIANT(computeL1Norm(v1, v2) == 0, " ");
  for (i = 0; i < 30; ++i) {
    v2.setVal(i, i%2);
  }
  
  CHECK_INVARIANT(computeL1Norm(v1, v2) == 30, " ");
  
  for (i = 0; i < 30; ++i) {
    if (i%3 == 0) {
      v2.setVal(i, 1);
    } else {
      v2.setVal(i,0);
    }
  }
  
  CHECK_INVARIANT(computeL1Norm(v1, v2) == 15, " ");

  DiscreteValueVect v21(DiscreteValueVect::TWOBITVALUE, 30);
  DiscreteValueVect v22(DiscreteValueVect::TWOBITVALUE, 30);
  for (i = 0; i < 30; ++i) {
    v21.setVal(i, i%4);
    v22.setVal(i, i%4);
  }
  CHECK_INVARIANT(computeL1Norm(v21, v22) == 0, " ");
  for (i = 0; i < 30; ++i) {
    v22.setVal(i, (i+1)%4);
  }
  CHECK_INVARIANT(computeL1Norm(v21, v22) == 44, " ");

  DiscreteValueVect v41(DiscreteValueVect::FOURBITVALUE, 16);
  DiscreteValueVect v42(DiscreteValueVect::FOURBITVALUE, 16);
  for (i = 0; i < 16; ++i) {
    v41.setVal(i, i%16);
    v42.setVal(i, i%16);
  }
  CHECK_INVARIANT(computeL1Norm(v41, v42) == 0, " ");

  for (i = 0; i < 16; ++i) {
    v42.setVal(i, i%5);
  }
  CHECK_INVARIANT(computeL1Norm(v41, v42) ==90, " ");

  DiscreteValueVect v43(v42);
  CHECK_INVARIANT(computeL1Norm(v42, v43) == 0, " ");

  DiscreteValueVect v81(DiscreteValueVect::EIGHTBITVALUE, 5);
  DiscreteValueVect v82(DiscreteValueVect::EIGHTBITVALUE, 5);
  v81.setVal(0, 34); v82.setVal(0, 34);
  v81.setVal(1, 167); v82.setVal(1, 167);
  v81.setVal(2, 3); v82.setVal(2, 3);
  v81.setVal(3, 56); v82.setVal(3, 56);
  v81.setVal(4, 128); v82.setVal(4, 128);
  CHECK_INVARIANT(computeL1Norm(v81, v82) == 0, " ");

  v82.setVal(0, 14); v82.setVal(1, 67);
  v82.setVal(2, 103); v82.setVal(3, 6);
  v82.setVal(4, 228); 
  CHECK_INVARIANT(computeL1Norm(v81, v82) == 370, "");

  DiscreteValueVect v161(DiscreteValueVect::SIXTEENBITVALUE, 3);
  DiscreteValueVect v162(DiscreteValueVect::SIXTEENBITVALUE, 3);
  v161.setVal(0, 2345); v162.setVal(0, 2345);
  v161.setVal(1, 64578); v162.setVal(1, 64578);
  v161.setVal(2, 34); v162.setVal(2, 34);
  CHECK_INVARIANT(computeL1Norm(v161, v162) == 0, " ");

  v162.setVal(0, 1345);
  v162.setVal(1, 54578);
  v162.setVal(2, 10034);
  CHECK_INVARIANT(computeL1Norm(v161, v162) == 21000, " ");

}

void test3DiscreteVectPickles() {
  DiscreteValueVect v1(DiscreteValueVect::ONEBITVALUE, 30);
  unsigned int i;
  for (i = 0; i < 15; ++i) {
    v1.setVal(2*i, 1);
  }
  DiscreteValueVect v2(v1.toString());
  CHECK_INVARIANT(computeL1Norm(v1, v2) == 0, " ");

  DiscreteValueVect v21(DiscreteValueVect::TWOBITVALUE, 30);
  for (i = 0; i < 30; ++i) {
    v21.setVal(i, i%4);
  }
  DiscreteValueVect v22(v21.toString());
  CHECK_INVARIANT(computeL1Norm(v21, v22) == 0, " ");

  DiscreteValueVect v41(DiscreteValueVect::FOURBITVALUE, 16);
  for (i = 0; i < 16; ++i) {
    v41.setVal(i, i%16);
  }
  DiscreteValueVect v42(v41.toString());
  CHECK_INVARIANT(computeL1Norm(v41, v42) == 0, " ");

  DiscreteValueVect v81(DiscreteValueVect::EIGHTBITVALUE, 5);
  v81.setVal(0, 34); 
  v81.setVal(1, 167);
  v81.setVal(2, 3); 
  v81.setVal(3, 56);
  v81.setVal(4, 128);
  DiscreteValueVect v82(v81.toString());
  CHECK_INVARIANT(computeL1Norm(v81, v82) == 0, " ");

  DiscreteValueVect v161(DiscreteValueVect::SIXTEENBITVALUE, 3);
  v161.setVal(0, 2345);
  v161.setVal(1, 64578);
  v161.setVal(2, 34);
  DiscreteValueVect v162(v161.toString());
  CHECK_INVARIANT(computeL1Norm(v161, v162) == 0, " ");
}

int main(){
  try{
    throw IndexErrorException(3);
  } catch (IndexErrorException) {
    cerr << "pass" << endl;
  }

  stringstream ss(ios_base::binary|ios_base::out|ios_base::in);
  int v1=4,v2=5,v3,v4;

  ss.write((const char *)&v1,sizeof(v1));
  ss.write((const char *)&v2,sizeof(v2));
#if 0
  ss.close();
  fstream ss2("blah.bin",ios_base::binary|ios_base::in);
  ss2.read((char *)&v3,sizeof(v3));
  ss2.read((char *)&v4,sizeof(v4));
#endif
  ss.seekp(0,ios_base::beg);
  ss.read((char *)&v3,sizeof(v3));
  ss.read((char *)&v4,sizeof(v4));
  
  TXTMSG("v3",v3);
  TXTMSG("v4",v4);
  
  cerr << " SPARSE -----------------------------------------" << endl;
  SparseBitVect sparseFoo(10);
  Test(sparseFoo);
  TaniTest(sparseFoo);
  ProbeTest(sparseFoo);
  cerr << " Explicit -----------------------------------------" << endl;
  ExplicitBitVect explicitFoo(10);
  Test(explicitFoo);
  TaniTest(explicitFoo);
  cerr << " Done" << endl;
  
  std::cout << " Test DiscreteValue Vectors 1 ------------------------------------" << endl;
  test1DiscreteVect();
  std::cout << " Test DiscreteValue Vectors 2 ------------------------------------" << endl;
  test2DiscreteVectDists();
  std::cout << " Test DiscreteValue Vectors 3 ------------------------------------" << endl;
  test3DiscreteVectPickles();

  return 0;
  
}
