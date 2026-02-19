//
//  Copyright (C) 2001-2024 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>

#include <RDGeneral/test.h>
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
#include "RealValueVect.h"
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Exceptions.h>
#include <DataStructs/SparseIntVect.h>
#include <DataStructs/DatastructsStreamOps.h>
#include <cstdlib>

using namespace std;
using namespace RDKit;

template <typename T>
inline void TXTMSG(const char *__a__, T __b__) {
  BOOST_LOG(rdInfoLog) << (__a__) << " " << (__b__) << std::endl;
}

template <typename T>
void Test(T arg) {
  (void)arg;
  T t1(20);
  TXTMSG("Set 10:", t1.setBit(10));
  TXTMSG("Set 11:", t1.setBit(11));
  TXTMSG("Set 14:", t1.setBit(14));
  TXTMSG("Set 10:", t1.setBit(10));
  TXTMSG("Get 14:", t1.getBit(14));
  TXTMSG("Num:", t1.getNumBits());
  TXTMSG("NumOn:", t1.getNumOnBits());
  TXTMSG("NumOff:", t1.getNumOffBits());
  REQUIRE(t1 == t1);

  IntVect onBits;
  t1.getOnBits(onBits);
  std::copy(onBits.begin(), onBits.end(),
            std::ostream_iterator<int>(std::cout, ", "));
  std::cout << std::endl;

  T t2(t1);
  // onBits = t2.getOnBits();
  TXTMSG("t2[19]:", t2[19]);
  TXTMSG("t2[14]:", t2[14]);

  REQUIRE(t2 == t1);

  t2 = t1;
  // onBits = t2.getOnBits();
  TXTMSG("t2[19]:", t2[19]);
  t2.unsetBit(14);
  TXTMSG("t2[14]:", t2[14]);
  t2.setBit(15);
  t2.setBit(17);
  REQUIRE(t2 != t1);

  std::cout << "t1: ";
  t1.getOnBits(onBits);
  std::copy(onBits.begin(), onBits.end(),
            std::ostream_iterator<int>(std::cout, ", "));
  std::cout << std::endl;

  std::cout << "t2: ";
  t2.getOnBits(onBits);
  std::copy(onBits.begin(), onBits.end(),
            std::ostream_iterator<int>(std::cout, ", "));
  std::cout << std::endl;

  std::cout << "t1|t2: ";
  T t3 = t1 | t2;
  t3.getOnBits(onBits);
  std::copy(onBits.begin(), onBits.end(),
            std::ostream_iterator<int>(std::cout, ", "));
  std::cout << std::endl;

  std::cout << "t1&t2: ";
  t3 = t1 & t2;
  t3.getOnBits(onBits);
  std::copy(onBits.begin(), onBits.end(),
            std::ostream_iterator<int>(std::cout, ", "));
  std::cout << std::endl;

  std::cout << "t1^t2: ";
  t3 = t1 ^ t2;
  t3.getOnBits(onBits);
  std::copy(onBits.begin(), onBits.end(),
            std::ostream_iterator<int>(std::cout, ", "));
  std::cout << std::endl;

  std::cout << "~t1: ";
  t3 = ~t1;
  t3.getOnBits(onBits);
  std::copy(onBits.begin(), onBits.end(),
            std::ostream_iterator<int>(std::cout, ", "));
  std::cout << std::endl;

  try {
    t3.getBit(4000);
  } catch (IndexErrorException &) {
    std::cout << " except " << endl;
  } catch (...) {
    std::cout << " ERROR EXCEPT " << endl;
  }

  T t4(t1.toString());
  REQUIRE(t4 == t1);
  REQUIRE_THAT(TanimotoSimilarity(t1, t4),
               Catch::Matchers::WithinAbs(1.0, 1e-4));

  T *t5 = FoldFingerprint(t1);
  REQUIRE(t5->getNumBits() == t1.getNumBits() / 2);
  REQUIRE(t5->getBit(0));
  REQUIRE(t5->getBit(1));
  REQUIRE(t5->getBit(4));
  REQUIRE(!t5->getBit(2));
  REQUIRE(!t5->getBit(3));
  delete t5;

  std::string pkl = t1.toString();
  const char *pkl64 = Base64Encode(pkl.c_str(), pkl.size());
  T t6(t1.getNumBits());
  t6.initFromText(pkl64, strlen(pkl64), true);
  delete[] pkl64;
  REQUIRE(t6 == t1);
}

template <typename T>
void TaniTest(T &arg) {
  (void)arg;  // unused var;
  std::string fps[4] = {
      ".b+HHa.EgU6+ibEIr89.CpX0g8FZiXH+R0+Ps.mr6tg.2",
      ".b7HEa..ccc+gWEIr89.8lV8gOF3aXFFR.+Ps.mZ6lg.2",
      ".H+nHq2EcY09y5EIr9e.8p50h0NgiWGNx4+Hm+Gbslw.2",
      ".1.HHa..cUI6i5E2rO8.Op10d0NoiWGVx.+Hm.Gb6lo.2",
  };
  double dists[] = {1.0,      0.788991, 0.677165, 0.686957, 1.0,
                    0.578125, 0.591304, 1.0,      0.732759, 1.0};
  int idx = 0;
  for (int i = 0; i < 4; i++) {
    T v1(256);
    FromDaylightString(v1, fps[i]);
    for (int j = i; j < 4; j++) {
      T v2(256);
      FromDaylightString(v2, fps[j]);
      double tani = TanimotoSimilarity(v1, v2);
      REQUIRE_THAT(tani, Catch::Matchers::WithinAbs(dists[idx], 1e-4));
      tani = TverskySimilarity(v1, v2, 1., 1.);
      REQUIRE_THAT(tani, Catch::Matchers::WithinAbs(dists[idx], 1e-4));
      tani = SimilarityWrapper(v1, v2, TanimotoSimilarity<T, T>);
      REQUIRE_THAT(tani, Catch::Matchers::WithinAbs(dists[idx], 1e-4));
      tani = SimilarityWrapper(v1, v2, 1., 1., TverskySimilarity<T, T>);
      REQUIRE_THAT(tani, Catch::Matchers::WithinAbs(dists[idx], 1e-4));
      tani = SimilarityWrapper(v1, v2, TanimotoSimilarity<T, T>, true);
      REQUIRE_THAT(tani, Catch::Matchers::WithinAbs(1. - dists[idx], 1e-4));
      tani = SimilarityWrapper(v1, v2, 1., 1., TverskySimilarity<T, T>, true);
      REQUIRE_THAT(tani, Catch::Matchers::WithinAbs(1. - dists[idx], 1e-4));
      idx++;
    }
  }
}

template <typename T>
void ProbeTest(T &arg) {
  (void)arg;  // unused var
  int sz = 1000;
  T t1(sz), t2(sz);
  for (int i = 0; i < sz; i += 2) {
    t1.setBit(i);
    if (i < 3 * sz / 4) {
      t2.setBit(i);
    }
  }
  std::string pkl = t1.toString();
  REQUIRE(AllProbeBitsMatch(t1, pkl));
  REQUIRE(AllProbeBitsMatch(t2, pkl));
  REQUIRE(AllProbeBitsMatch(t1.toString(), pkl));
  REQUIRE(AllProbeBitsMatch(t2.toString(), pkl));
  REQUIRE(AllProbeBitsMatch(t1.toString().c_str(), pkl.c_str()));
  REQUIRE(AllProbeBitsMatch(t2.toString().c_str(), pkl.c_str()));
  pkl = t2.toString();
  REQUIRE(!AllProbeBitsMatch(t1, pkl));
  REQUIRE(AllProbeBitsMatch(t2, pkl));
  REQUIRE(!AllProbeBitsMatch(t1.toString(), pkl));
  REQUIRE(AllProbeBitsMatch(t2.toString(), pkl));
  REQUIRE(!AllProbeBitsMatch(t1.toString().c_str(), pkl.c_str()));
  REQUIRE(AllProbeBitsMatch(t2.toString().c_str(), pkl.c_str()));
}

TEST_CASE("test1DiscreteVect") {
  DiscreteValueVect vect1(DiscreteValueVect::ONEBITVALUE, 30);
  unsigned int i;
  for (i = 0; i < 15; ++i) {
    vect1.setVal(2 * i, 1);
  }

  REQUIRE(vect1.getLength() == 30);
  REQUIRE(vect1.getTotalVal() == 15);
  for (i = 0; i < vect1.getLength(); ++i) {
    REQUIRE(vect1.getVal(i) == (i + 1) % 2);
  }
  CHECK_THROWS_AS(vect1.setVal(28, 2), ValueErrorException);

  // all these tests should fail if unsigned int changes from being
  // 32 bits
  DiscreteValueVect vect2(DiscreteValueVect::TWOBITVALUE, 30);
  for (i = 0; i < vect2.getLength(); ++i) {
    vect2.setVal(i, i % 4);
  }

  for (i = 0; i < vect2.getLength(); ++i) {
    REQUIRE(vect2.getVal(i) == i % 4);
  }
  REQUIRE(vect2.getTotalVal() == 43);
  CHECK_THROWS_AS(vect2.setVal(28, 10), ValueErrorException);

  DiscreteValueVect vect4(DiscreteValueVect::FOURBITVALUE, 30);
  for (i = 0; i < vect4.getLength(); ++i) {
    vect4.setVal(i, i % 16);
  }

  for (i = 0; i < vect4.getLength(); ++i) {
    REQUIRE(vect4.getVal(i) == i % 16);
  }
  REQUIRE(vect4.getTotalVal() == 211);
  CHECK_THROWS_AS(vect4.setVal(28, 16), ValueErrorException);

  DiscreteValueVect vect8(DiscreteValueVect::EIGHTBITVALUE, 32);
  for (i = 0; i < vect8.getLength(); ++i) {
    vect8.setVal(i, i % 256);
  }

  for (i = 0; i < vect8.getLength(); ++i) {
    REQUIRE(vect8.getVal(i) == i % 256);
  }
  REQUIRE(vect8.getTotalVal() == 496);
  CHECK_THROWS_AS(vect8.setVal(28, 257), ValueErrorException);

  DiscreteValueVect vect16(DiscreteValueVect::SIXTEENBITVALUE, 300);
  for (i = 0; i < vect16.getLength(); ++i) {
    vect16.setVal(i, i % 300);
  }

  for (i = 0; i < vect16.getLength(); ++i) {
    REQUIRE(vect16.getVal(i) == i % 300);
  }

  REQUIRE(vect16.getTotalVal() == 44850);
  vect16.setVal(28, 65535);
  CHECK_THROWS_AS(vect16.setVal(28, 65536), ValueErrorException);
}

TEST_CASE("test2DiscreteVectDists") {
  DiscreteValueVect v1(DiscreteValueVect::ONEBITVALUE, 30);
  DiscreteValueVect v2(DiscreteValueVect::ONEBITVALUE, 30);
  unsigned int i;
  for (i = 0; i < 15; ++i) {
    v1.setVal(2 * i, 1);
    v2.setVal(2 * i, 1);
  }
  REQUIRE(computeL1Norm(v1, v2) == 0);
  for (i = 0; i < 30; ++i) {
    v2.setVal(i, i % 2);
  }

  REQUIRE(computeL1Norm(v1, v2) == 30);

  for (i = 0; i < 30; ++i) {
    if (i % 3 == 0) {
      v2.setVal(i, 1);
    } else {
      v2.setVal(i, 0);
    }
  }

  REQUIRE(computeL1Norm(v1, v2) == 15);

  DiscreteValueVect v21(DiscreteValueVect::TWOBITVALUE, 30);
  DiscreteValueVect v22(DiscreteValueVect::TWOBITVALUE, 30);
  for (i = 0; i < 30; ++i) {
    v21.setVal(i, i % 4);
    v22.setVal(i, i % 4);
  }
  REQUIRE(computeL1Norm(v21, v22) == 0);
  for (i = 0; i < 30; ++i) {
    v22.setVal(i, (i + 1) % 4);
  }
  REQUIRE(computeL1Norm(v21, v22) == 44);

  DiscreteValueVect v41(DiscreteValueVect::FOURBITVALUE, 16);
  DiscreteValueVect v42(DiscreteValueVect::FOURBITVALUE, 16);
  for (i = 0; i < 16; ++i) {
    v41.setVal(i, i % 16);
    v42.setVal(i, i % 16);
  }
  REQUIRE(computeL1Norm(v41, v42) == 0);

  for (i = 0; i < 16; ++i) {
    v42.setVal(i, i % 5);
  }
  REQUIRE(computeL1Norm(v41, v42) == 90);

  DiscreteValueVect v43(v42);
  REQUIRE(computeL1Norm(v42, v43) == 0);

  DiscreteValueVect v81(DiscreteValueVect::EIGHTBITVALUE, 5);
  DiscreteValueVect v82(DiscreteValueVect::EIGHTBITVALUE, 5);
  v81.setVal(0, 34);
  v82.setVal(0, 34);
  v81.setVal(1, 167);
  v82.setVal(1, 167);
  v81.setVal(2, 3);
  v82.setVal(2, 3);
  v81.setVal(3, 56);
  v82.setVal(3, 56);
  v81.setVal(4, 128);
  v82.setVal(4, 128);
  REQUIRE(computeL1Norm(v81, v82) == 0);

  v82.setVal(0, 14);
  v82.setVal(1, 67);
  v82.setVal(2, 103);
  v82.setVal(3, 6);
  v82.setVal(4, 228);
  REQUIRE(computeL1Norm(v81, v82) == 370);

  DiscreteValueVect v161(DiscreteValueVect::SIXTEENBITVALUE, 3);
  DiscreteValueVect v162(DiscreteValueVect::SIXTEENBITVALUE, 3);
  v161.setVal(0, 2345);
  v162.setVal(0, 2345);
  v161.setVal(1, 64578);
  v162.setVal(1, 64578);
  v161.setVal(2, 34);
  v162.setVal(2, 34);
  REQUIRE(computeL1Norm(v161, v162) == 0);

  v162.setVal(0, 1345);
  v162.setVal(1, 54578);
  v162.setVal(2, 10034);
  REQUIRE(computeL1Norm(v161, v162) == 21000);
}

TEST_CASE("test3DiscreteVectPickles") {
  DiscreteValueVect v1(DiscreteValueVect::ONEBITVALUE, 30);
  unsigned int i;
  for (i = 0; i < 15; ++i) {
    v1.setVal(2 * i, 1);
  }
  DiscreteValueVect v2(v1.toString());
  REQUIRE(computeL1Norm(v1, v2) == 0);

  DiscreteValueVect v21(DiscreteValueVect::TWOBITVALUE, 30);
  for (i = 0; i < 30; ++i) {
    v21.setVal(i, i % 4);
  }
  DiscreteValueVect v22(v21.toString());
  REQUIRE(computeL1Norm(v21, v22) == 0);

  DiscreteValueVect v41(DiscreteValueVect::FOURBITVALUE, 16);
  for (i = 0; i < 16; ++i) {
    v41.setVal(i, i % 16);
  }
  DiscreteValueVect v42(v41.toString());
  REQUIRE(computeL1Norm(v41, v42) == 0);

  DiscreteValueVect v81(DiscreteValueVect::EIGHTBITVALUE, 5);
  v81.setVal(0, 34);
  v81.setVal(1, 167);
  v81.setVal(2, 3);
  v81.setVal(3, 56);
  v81.setVal(4, 128);
  DiscreteValueVect v82(v81.toString());
  REQUIRE(computeL1Norm(v81, v82) == 0);

  DiscreteValueVect v161(DiscreteValueVect::SIXTEENBITVALUE, 3);
  v161.setVal(0, 2345);
  v161.setVal(1, 64578);
  v161.setVal(2, 34);
  DiscreteValueVect v162(v161.toString());
  REQUIRE(computeL1Norm(v161, v162) == 0);
}

TEST_CASE("test4DiscreteVectOps1") {
  DiscreteValueVect vect1(DiscreteValueVect::ONEBITVALUE, 8);
  for (unsigned int i = 0; i < 4; ++i) {
    vect1.setVal(2 * i, 1);
  }
  REQUIRE(vect1.getLength() == 8);
  REQUIRE(vect1.getTotalVal() == 4);

  DiscreteValueVect vect2(DiscreteValueVect::ONEBITVALUE, 8);
  for (unsigned int i = 0; i < 4; ++i) {
    vect2.setVal(2 * i + 1, 1);
  }
  REQUIRE(vect2.getTotalVal() == 4);

  DiscreteValueVect vect3 = vect1 & vect2;
  REQUIRE(vect3.getLength() == 8);
  REQUIRE(vect3.getTotalVal() == 0);

  DiscreteValueVect vect4 = vect1 | vect2;
  REQUIRE(vect4.getLength() == 8);
  REQUIRE(vect4.getTotalVal() == 8);
}

TEST_CASE("test5DiscreteVectOps2") {
  DiscreteValueVect vect1(DiscreteValueVect::TWOBITVALUE, 8);
  for (unsigned int i = 0; i < 4; ++i) {
    vect1.setVal(2 * i, 2);
  }
  REQUIRE(vect1.getLength() == 8);
  REQUIRE(vect1.getTotalVal() == 8);

  DiscreteValueVect vect2(DiscreteValueVect::TWOBITVALUE, 8);
  for (unsigned int i = 0; i < 4; ++i) {
    vect2.setVal(2 * i + 1, 2);
    vect2.setVal(2 * i, 1);
  }
  REQUIRE(vect2.getTotalVal() == 12);

  DiscreteValueVect vect3 = vect1 & vect2;
  REQUIRE(vect3.getLength() == 8);
  REQUIRE(vect3.getTotalVal() == 4);

  DiscreteValueVect vect4 = vect1 | vect2;
  REQUIRE(vect4.getLength() == 8);
  REQUIRE(vect4.getTotalVal() == 16);

  DiscreteValueVect vect5 = vect1 + vect2;
  REQUIRE(vect5.getLength() == 8);
  REQUIRE(vect5.getTotalVal() == 20);

  vect5 = vect1 - vect2;
  REQUIRE(vect5.getTotalVal() == 4);
  vect5 = vect2 - vect1;
  REQUIRE(vect5.getTotalVal() == 8);
}

TEST_CASE("test6SparseIntVect") {
  SparseIntVect<int> iVect(255);

  REQUIRE(iVect.getLength() == 255);
  REQUIRE(iVect.getVal(23) == 0);
  iVect.setVal(23, 14);
  REQUIRE(iVect.getVal(23) == 14);

  SparseIntVect<int> oVect(iVect);
  REQUIRE(oVect.getLength() == 255);
  REQUIRE(oVect.getVal(23) == 14);

  std::vector<int> tmpV(3);
  tmpV[0] = 1;
  tmpV[1] = 5;
  tmpV[2] = 1;
  REQUIRE(iVect.getVal(1) == 0);
  REQUIRE(iVect[1] == 0);
  REQUIRE(iVect.getVal(5) == 0);
  REQUIRE(iVect[5] == 0);
  updateFromSequence(iVect, tmpV);
  REQUIRE(iVect.getVal(1) == 2);
  REQUIRE(iVect[1] == 2);
  REQUIRE(iVect.getVal(5) == 1);
  REQUIRE(iVect[5] == 1);

  iVect.setVal(3, -4);
  REQUIRE(iVect.getTotalVal() == 13);

  REQUIRE_THROWS_AS(iVect.setVal(-1, 13), IndexErrorException);
  REQUIRE_THROWS_AS(iVect.setVal(255, 42), IndexErrorException);
  REQUIRE_THROWS_AS(iVect.getVal(-1), IndexErrorException);
  REQUIRE_THROWS_AS(iVect.getVal(255), IndexErrorException);
  REQUIRE_THROWS_AS(iVect[-1], IndexErrorException);

  {
    SparseIntVect<int> iV1(5);
    iV1.setVal(4, 4);
    iV1.setVal(0, 2);
    iV1.setVal(3, 1);
    auto iter = iV1.getNonzeroElements().begin();
    REQUIRE(iter->first == 0);
    REQUIRE(iter->second == 2);
    ++iter;
    REQUIRE(iter->first == 3);
    REQUIRE(iter->second == 1);
    ++iter;
    REQUIRE(iter->first == 4);
    REQUIRE(iter->second == 4);
    ++iter;
    REQUIRE(iter == iV1.getNonzeroElements().end());
    REQUIRE_THAT(DiceSimilarity(iV1, iV1),
                 Catch::Matchers::WithinAbs(1., 1e-4));
  }

  {  // iV1 &= iV2
    SparseIntVect<int> iV1(5), iV2(5);
    iV1.setVal(0, 2);
    iV1.setVal(2, 1);
    iV1.setVal(3, 4);
    iV1.setVal(4, 4);

    iV2.setVal(1, 2);
    iV2.setVal(2, 3);
    iV2.setVal(3, 4);
    iV2.setVal(4, 6);

    REQUIRE_THAT(DiceSimilarity(iV1, iV2),
                 Catch::Matchers::WithinAbs(18. / 26., 1e-4));

    iV1 &= iV2;
    REQUIRE(iV1[0] == 0);
    REQUIRE(iV1[1] == 0);
    REQUIRE(iV1[2] == 1);
    REQUIRE(iV1[3] == 4);
    REQUIRE(iV1[4] == 4);

    REQUIRE(iV2[0] == 0);
    REQUIRE(iV2[1] == 2);
    REQUIRE(iV2[2] == 3);
    REQUIRE(iV2[3] == 4);
    REQUIRE(iV2[4] == 6);

    REQUIRE_THAT(DiceSimilarity(iV1, iV2),
                 Catch::Matchers::WithinAbs(18. / 24., 1e-4));
    REQUIRE_THAT(TverskySimilarity(iV1, iV2, 0.5, 0.5, false),
                 Catch::Matchers::WithinAbs(9. / 12., 1e-4));
    REQUIRE_THAT(TverskySimilarity(iV1, iV2, 1.0, 1.0, false),
                 Catch::Matchers::WithinAbs(9. / 15., 1e-4));
    REQUIRE_THAT(TanimotoSimilarity(iV1, iV2),
                 Catch::Matchers::WithinAbs(9. / 15., 1e-4));
    REQUIRE_THAT(TverskySimilarity(iV1, iV2, 0.333333333, 0.66666666667, false),
                 Catch::Matchers::WithinAbs(9. / 13., 1e-4));
    REQUIRE_THAT(TverskySimilarity(iV1, iV2, 1.0, 0.0, false),
                 Catch::Matchers::WithinAbs(9. / 9., 1e-4));

    REQUIRE_THROWS_AS(iV1 &= iVect, ValueErrorException);
  }

  {  // iV3 = iv1&iV2
    SparseIntVect<int> iV1(5), iV2(5), iV3(5);
    iV1.setVal(0, 2);
    iV1.setVal(2, 1);
    iV1.setVal(3, 4);
    iV1.setVal(4, 4);

    iV2.setVal(1, 2);
    iV2.setVal(2, 3);
    iV2.setVal(3, 4);
    iV2.setVal(4, 6);

    iV3 = iV1 & iV2;
    REQUIRE(iV3[0] == 0);
    REQUIRE(iV3[1] == 0);
    REQUIRE(iV3[2] == 1);
    REQUIRE(iV3[3] == 4);
    REQUIRE(iV3[4] == 4);

    REQUIRE(iV1[0] == 2);
    REQUIRE(iV1[1] == 0);
    REQUIRE(iV1[2] == 1);
    REQUIRE(iV1[3] == 4);
    REQUIRE(iV1[4] == 4);

    REQUIRE(iV2[0] == 0);
    REQUIRE(iV2[1] == 2);
    REQUIRE(iV2[2] == 3);
    REQUIRE(iV2[3] == 4);
    REQUIRE(iV2[4] == 6);
  }

  {  // iV2 &= iV1
    SparseIntVect<int> iV1(5), iV2(5);
    iV1.setVal(0, 2);
    iV1.setVal(2, 1);
    iV1.setVal(3, 4);
    iV1.setVal(4, 4);

    iV2.setVal(1, 2);
    iV2.setVal(2, 3);
    iV2.setVal(3, 4);
    iV2.setVal(4, 6);

    iV2 &= iV1;
    REQUIRE(iV2[0] == 0);
    REQUIRE(iV2[1] == 0);
    REQUIRE(iV2[2] == 1);
    REQUIRE(iV2[3] == 4);
    REQUIRE(iV2[4] == 4);

    REQUIRE(iV1[0] == 2);
    REQUIRE(iV1[1] == 0);
    REQUIRE(iV1[2] == 1);
    REQUIRE(iV1[3] == 4);
    REQUIRE(iV1[4] == 4);
    REQUIRE_THROWS_AS(iV2 &= iVect, ValueErrorException);
  }

  {  // iV1 |= iV2
    SparseIntVect<int> iV1(5), iV2(5);
    iV1.setVal(0, 2);
    iV1.setVal(2, 1);
    iV1.setVal(3, 4);
    iV1.setVal(4, 4);

    iV2.setVal(1, 2);
    iV2.setVal(2, 3);
    iV2.setVal(3, 4);
    iV2.setVal(4, 6);

    iV1 |= iV2;
    REQUIRE(iV1[0] == 2);
    REQUIRE(iV1[1] == 2);
    REQUIRE(iV1[2] == 3);
    REQUIRE(iV1[3] == 4);
    REQUIRE(iV1[4] == 6);

    REQUIRE(iV2[0] == 0);
    REQUIRE(iV2[1] == 2);
    REQUIRE(iV2[2] == 3);
    REQUIRE(iV2[3] == 4);
    REQUIRE(iV2[4] == 6);

    REQUIRE_THROWS_AS(iV1 |= iVect, ValueErrorException);
  }

  {  // iV3 = iv1 |iV2
    SparseIntVect<int> iV1(5), iV2(5), iV3(5);
    iV1.setVal(0, 2);
    iV1.setVal(2, 1);
    iV1.setVal(3, 4);
    iV1.setVal(4, 4);

    iV2.setVal(1, 2);
    iV2.setVal(2, 3);
    iV2.setVal(3, 4);
    iV2.setVal(4, 6);

    iV3 = iV1 | iV2;
    REQUIRE(iV3[0] == 2);
    REQUIRE(iV3[1] == 2);
    REQUIRE(iV3[2] == 3);
    REQUIRE(iV3[3] == 4);
    REQUIRE(iV3[4] == 6);

    REQUIRE(iV1[0] == 2);
    REQUIRE(iV1[1] == 0);
    REQUIRE(iV1[2] == 1);
    REQUIRE(iV1[3] == 4);
    REQUIRE(iV1[4] == 4);

    REQUIRE(iV2[0] == 0);
    REQUIRE(iV2[1] == 2);
    REQUIRE(iV2[2] == 3);
    REQUIRE(iV2[3] == 4);
    REQUIRE(iV2[4] == 6);
  }

  {  // iV2 |= iV1
    SparseIntVect<int> iV1(5), iV2(5);
    iV1.setVal(0, 2);
    iV1.setVal(2, 1);
    iV1.setVal(3, 4);
    iV1.setVal(4, 4);

    iV2.setVal(1, 2);
    iV2.setVal(2, 3);
    iV2.setVal(3, 4);
    iV2.setVal(4, 6);

    iV2 |= iV1;
    REQUIRE(iV2[0] == 2);
    REQUIRE(iV2[1] == 2);
    REQUIRE(iV2[2] == 3);
    REQUIRE(iV2[3] == 4);
    REQUIRE(iV2[4] == 6);

    REQUIRE(iV1[0] == 2);
    REQUIRE(iV1[1] == 0);
    REQUIRE(iV1[2] == 1);
    REQUIRE(iV1[3] == 4);
    REQUIRE(iV1[4] == 4);
  }

  {  // iV1 += iV2
    SparseIntVect<int> iV1(5), iV2(5);
    iV1.setVal(0, 2);
    iV1.setVal(2, 1);
    iV1.setVal(3, 4);
    iV1.setVal(4, 4);

    iV2.setVal(1, 2);
    iV2.setVal(2, 3);
    iV2.setVal(3, -4);
    iV2.setVal(4, 6);

    iV1 += iV2;
    REQUIRE(iV1[0] == 2);
    REQUIRE(iV1[1] == 2);
    REQUIRE(iV1[2] == 4);
    REQUIRE(iV1[3] == 0);
    REQUIRE(iV1[4] == 10);

    REQUIRE(iV2[0] == 0);
    REQUIRE(iV2[1] == 2);
    REQUIRE(iV2[2] == 3);
    REQUIRE(iV2[3] == -4);
    REQUIRE(iV2[4] == 6);
  }
  {  // iV3 = IV1 + iV2
    SparseIntVect<int> iV1(5), iV2(5), iV3(5);
    iV1.setVal(2, 1);
    iV1.setVal(0, 2);
    iV1.setVal(4, 4);
    iV1.setVal(3, 4);

    iV2.setVal(1, 2);
    iV2.setVal(2, 3);
    iV2.setVal(3, -4);
    iV2.setVal(4, 6);

    iV3 = iV1 + iV2;
    REQUIRE(iV3[0] == 2);
    REQUIRE(iV3[1] == 2);
    REQUIRE(iV3[2] == 4);
    REQUIRE(iV3[3] == 0);
    REQUIRE(iV3[4] == 10);

    REQUIRE(iV1[0] == 2);
    REQUIRE(iV1[1] == 0);
    REQUIRE(iV1[2] == 1);
    REQUIRE(iV1[3] == 4);
    REQUIRE(iV1[4] == 4);
    REQUIRE(iV2[0] == 0);
    REQUIRE(iV2[1] == 2);
    REQUIRE(iV2[2] == 3);
    REQUIRE(iV2[3] == -4);
    REQUIRE(iV2[4] == 6);
  }

  {  // iV1 -= iV2
    SparseIntVect<int> iV1(5), iV2(5);
    iV1.setVal(0, 2);
    iV1.setVal(2, 1);
    iV1.setVal(3, 4);
    iV1.setVal(4, 4);

    iV2.setVal(1, 2);
    iV2.setVal(2, 3);
    iV2.setVal(3, 4);
    iV2.setVal(4, 6);

    iV1 -= iV2;
    REQUIRE(iV1[0] == 2);
    REQUIRE(iV1[1] == -2);
    REQUIRE(iV1[2] == -2);
    REQUIRE(iV1[3] == 0);
    REQUIRE(iV1[4] == -2);

    REQUIRE(iV2[0] == 0);
    REQUIRE(iV2[2] == 3);
    REQUIRE(iV2[1] == 2);
    REQUIRE(iV2[4] == 6);
    REQUIRE(iV2[3] == 4);
  }
  {  // iV3 = IV1 - iV2
    SparseIntVect<int> iV1(5), iV2(5), iV3(5);
    iV1.setVal(0, 2);
    iV1.setVal(2, 1);
    iV1.setVal(3, 4);
    iV1.setVal(4, 4);

    iV2.setVal(3, 4);
    iV2.setVal(4, 6);
    iV2.setVal(1, 2);
    iV2.setVal(2, 3);

    iV3 = iV1 - iV2;
    REQUIRE(iV3[0] == 2);
    REQUIRE(iV3[1] == -2);
    REQUIRE(iV3[2] == -2);
    REQUIRE(iV3[3] == 0);
    REQUIRE(iV3[4] == -2);

    REQUIRE(iV1[0] == 2);
    REQUIRE(iV1[1] == 0);
    REQUIRE(iV1[2] == 1);
    REQUIRE(iV1[3] == 4);
    REQUIRE(iV1[4] == 4);
    REQUIRE(iV2[0] == 0);
    REQUIRE(iV2[1] == 2);
    REQUIRE(iV2[2] == 3);
    REQUIRE(iV2[3] == 4);
    REQUIRE(iV2[4] == 6);
  }

  {  // operator== and operator!=
    SparseIntVect<int> iV1(5), iV2(5), iV3(3), iV4(5);
    iV1.setVal(0, 2);
    iV1.setVal(2, 1);
    iV1.setVal(3, 4);
    iV1.setVal(4, 4);

    iV2.setVal(1, 2);
    iV2.setVal(2, 3);
    iV2.setVal(3, 4);
    iV2.setVal(4, 6);

    iV4.setVal(1, 2);
    iV4.setVal(2, 3);
    iV4.setVal(3, 4);
    iV4.setVal(4, 6);

    REQUIRE(iV1 == iV1);
    REQUIRE(iV2 == iV2);
    REQUIRE(iV3 == iV3);
    REQUIRE(iV1 != iV2);
    REQUIRE(iV1 != iV3);
    REQUIRE(iV2 != iV1);
    REQUIRE(iV3 != iV1);
    REQUIRE(iV1 != iV3);
    REQUIRE(iV1 != iV4);
    REQUIRE(iV2 == iV4);
  }

  {  // test negative values (was sf.net Issue 3295215)
    SparseIntVect<int> iV1(5), iV2(5);
    iV1.setVal(0, -2);
    iV1.setVal(2, 1);
    iV1.setVal(3, -4);
    iV1.setVal(4, 4);

    iV2.setVal(1, -2);
    iV2.setVal(2, 3);
    iV2.setVal(3, -4);
    iV2.setVal(4, 6);

    REQUIRE_THAT(DiceSimilarity(iV1, iV2),
                 Catch::Matchers::WithinAbs(18. / 26., 1e-4));
  }
}

TEST_CASE("test7SparseIntVectPickles") {
  {
    SparseIntVect<int> iV1(5), iV2(3);
    iV1.setVal(0, 2);
    iV1.setVal(2, 1);
    iV1.setVal(3, 4);
    iV1.setVal(4, 4);

    iV2.setVal(2, 3);
    REQUIRE(iV1 != iV2);
    std::string pkl;
    pkl = iV1.toString();
    iV2.fromString(pkl);
    REQUIRE(iV1 == iV2);
  }

  {
    SparseIntVect<char> iV1(5);
    SparseIntVect<int> iV2(3);
    iV1.setVal(0, 2);
    iV1.setVal(2, 1);
    iV1.setVal(3, 4);

    iV2.setVal(1, 1);
    std::string pkl;
    pkl = iV1.toString();
    iV2.fromString(pkl);
    REQUIRE(iV2.getLength() == iV1.getLength());
    REQUIRE(iV2[0] == 2);
    REQUIRE(iV2[1] == 0);
    REQUIRE(iV2[2] == 1);
    REQUIRE(iV2[3] == 4);
  }

  {
    SparseIntVect<int> iV1(5);
    SparseIntVect<char> iV2(3);
    iV1.setVal(0, 2);
    iV1.setVal(2, 1);
    iV1.setVal(3, 4);

    std::string pkl;
    pkl = iV1.toString();
    REQUIRE_THROWS_AS(iV2.fromString(pkl), ValueErrorException);
  }
}

TEST_CASE("test8BitVectPickles") {
  {
    std::string dirName = getenv("RDBASE");
    dirName += "/Code/DataStructs/testData/";
    std::string pklName = dirName + "test1.bin";
    std::ifstream inS;
    inS.open(pklName.c_str(), std::ios_base::binary);
    unsigned int length;
    inS >> length;
    auto *buff = new char[length];
    length = inS.readsome(buff, length);
    inS.close();
    std::string pkl(buff, length);
    delete[] buff;
    ExplicitBitVect bv(pkl);

    REQUIRE(bv.getNumBits() == 32);
    REQUIRE(bv.getNumOnBits() == 16);
    REQUIRE(bv[0]);
    REQUIRE(!bv[1]);
  }
}

TEST_CASE("test9BitVectFPS") {
  {
    ExplicitBitVect bv(32);
    std::string fps;

    fps = BitVectToFPSText(bv);
    REQUIRE(fps == "00000000");

    bv.setBit(0);
    bv.setBit(1);
    bv.setBit(17);
    bv.setBit(23);
    bv.setBit(31);

    fps = BitVectToFPSText(bv);
    REQUIRE(fps == "03008280");
  }
  {
    ExplicitBitVect bv(32), bv2(32);
    std::string fps;

    fps = BitVectToFPSText(bv);
    REQUIRE(fps == "00000000");
    UpdateBitVectFromFPSText(bv2, fps);
    REQUIRE(bv == bv2);

    bv.setBit(0);
    bv.setBit(1);
    bv.setBit(4);
    bv.setBit(17);
    bv.setBit(23);
    bv.setBit(31);

    fps = BitVectToFPSText(bv);
    UpdateBitVectFromFPSText(bv2, fps);
    REQUIRE(bv == bv2);
  }
  {
    ExplicitBitVect bv(33);
    std::string fps;

    fps = BitVectToFPSText(bv);
    REQUIRE(fps == "0000000000");

    bv.setBit(0);
    bv.setBit(32);

    fps = BitVectToFPSText(bv);
    REQUIRE(fps == "0100000001");
  }
  {
    ExplicitBitVect bv(33), bv2(33);
    std::string fps;

    fps = BitVectToFPSText(bv);
    REQUIRE(fps == "0000000000");
    UpdateBitVectFromFPSText(bv2, fps);
    REQUIRE(bv == bv2);

    bv.setBit(0);
    bv.setBit(1);
    bv.setBit(4);
    bv.setBit(17);
    bv.setBit(23);
    bv.setBit(32);

    fps = BitVectToFPSText(bv);
    UpdateBitVectFromFPSText(bv2, fps);
    REQUIRE(bv == bv2);
  }
}

TEST_CASE("test10BitVectBinaryText") {
  {
    ExplicitBitVect bv(32);
    std::string fps;

    fps = BitVectToBinaryText(bv);
    REQUIRE(fps.size() == 4);
    for (char fp : fps) {
      REQUIRE(fp == 0);
    }

    bv.setBit(0);
    bv.setBit(9);
    bv.setBit(17);
    bv.setBit(26);

    fps = BitVectToBinaryText(bv);
    REQUIRE(fps.size() == 4);
    for (char fp : fps) {
      REQUIRE(fp != 0);
    }
  }
  {
    ExplicitBitVect bv(32), bv2(32);
    std::string fps;

    fps = BitVectToBinaryText(bv);
    REQUIRE(fps.size() == 4);
    for (char fp : fps) {
      REQUIRE(fp == 0);
    }
    UpdateBitVectFromBinaryText(bv2, fps);
    REQUIRE(bv == bv2);

    bv.setBit(0);
    bv.setBit(1);
    bv.setBit(4);
    bv.setBit(17);
    bv.setBit(23);
    bv.setBit(31);

    fps = BitVectToBinaryText(bv);
    UpdateBitVectFromBinaryText(bv2, fps);
    REQUIRE(bv == bv2);
  }
  {
    ExplicitBitVect bv(33);
    std::string fps;

    fps = BitVectToBinaryText(bv);
    REQUIRE(fps.size() == 5);
    for (char fp : fps) {
      REQUIRE(fp == 0);
    }

    bv.setBit(0);
    bv.setBit(32);

    fps = BitVectToBinaryText(bv);
    REQUIRE(fps.size() == 5);
    REQUIRE(fps[0] != 0);
    for (unsigned int i = 1; i < fps.size() - 1; ++i) {
      REQUIRE(fps[i] == 0);
    }
    REQUIRE(fps[fps.size() - 1] != 0);
  }
  {
    ExplicitBitVect bv(33), bv2(33);
    std::string fps;

    fps = BitVectToBinaryText(bv);
    REQUIRE(fps.size() == 5);
    UpdateBitVectFromBinaryText(bv2, fps);
    REQUIRE(bv == bv2);

    bv.setBit(0);
    bv.setBit(1);
    bv.setBit(4);
    bv.setBit(17);
    bv.setBit(23);
    bv.setBit(32);

    fps = BitVectToBinaryText(bv);
    UpdateBitVectFromBinaryText(bv2, fps);
    REQUIRE(bv == bv2);
  }
}

TEST_CASE("test11SimilaritiesBV") {
  // similarity = 1.0
  ExplicitBitVect bv(10);
  bv.setBit(0);
  bv.setBit(1);
  bv.setBit(4);
  bv.setBit(6);
  bv.setBit(9);
  ExplicitBitVect bv2 = bv;

  REQUIRE_THAT(TanimotoSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(1., 1e-4));
  REQUIRE_THAT(DiceSimilarity(bv, bv2), Catch::Matchers::WithinAbs(1., 1e-4));
  REQUIRE_THAT(CosineSimilarity(bv, bv2), Catch::Matchers::WithinAbs(1., 1e-4));
  REQUIRE_THAT(KulczynskiSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(1., 1e-4));
  REQUIRE_THAT(SokalSimilarity(bv, bv2), Catch::Matchers::WithinAbs(1., 1e-4));
  REQUIRE_THAT(McConnaugheySimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(1., 1e-4));
  REQUIRE_THAT(BraunBlanquetSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(1., 1e-4));
  REQUIRE_THAT(RusselSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(0.5, 1e-4));
  REQUIRE_THAT(RogotGoldbergSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(1., 1e-4));
  REQUIRE_THAT(AllBitSimilarity(bv, bv2), Catch::Matchers::WithinAbs(1., 1e-4));

  // similarity = 0.0
  bv2 = ExplicitBitVect(10);
  bv2.setBit(2);
  bv2.setBit(3);
  bv2.setBit(5);
  bv2.setBit(7);
  bv2.setBit(8);

  REQUIRE_THAT(TanimotoSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(0., 1e-4));
  REQUIRE_THAT(DiceSimilarity(bv, bv2), Catch::Matchers::WithinAbs(0., 1e-4));
  REQUIRE_THAT(CosineSimilarity(bv, bv2), Catch::Matchers::WithinAbs(0., 1e-4));
  REQUIRE_THAT(KulczynskiSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(0., 1e-4));
  REQUIRE_THAT(SokalSimilarity(bv, bv2), Catch::Matchers::WithinAbs(0., 1e-4));
  REQUIRE_THAT(McConnaugheySimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(-1., 1e-4));
  REQUIRE_THAT(BraunBlanquetSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(0., 1e-4));
  REQUIRE_THAT(RusselSimilarity(bv, bv2), Catch::Matchers::WithinAbs(0., 1e-4));
  REQUIRE_THAT(RogotGoldbergSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(0., 1e-4));
  REQUIRE_THAT(AllBitSimilarity(bv, bv2), Catch::Matchers::WithinAbs(0., 1e-4));

  // similarity ~= 0.5
  bv.setBit(5);
  bv2 = ExplicitBitVect(10);
  bv2.setBit(0);
  bv2.setBit(2);
  bv2.setBit(4);
  bv2.setBit(5);
  bv2.setBit(8);
  bv2.setBit(9);

  REQUIRE_THAT(TanimotoSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(0.5, 1e-4));
  REQUIRE_THAT(DiceSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(0.6666, 1e-4));
  REQUIRE_THAT(CosineSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(0.6666, 1e-4));
  REQUIRE_THAT(KulczynskiSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(0.6666, 1e-4));
  REQUIRE_THAT(SokalSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(0.3333, 1e-4));
  REQUIRE_THAT(McConnaugheySimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(0.3333, 1e-4));
  REQUIRE_THAT(BraunBlanquetSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(0.6666, 1e-4));
  REQUIRE_THAT(RusselSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(0.4, 1e-4));
  REQUIRE_THAT(RogotGoldbergSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(0.5833, 1e-4));
  REQUIRE_THAT(AllBitSimilarity(bv, bv2),
               Catch::Matchers::WithinAbs(0.6, 1e-4));
}

TEST_CASE("test12SimilaritiesSparseBV") {
  // similarity = 1.0
  SparseBitVect sbv(10);
  sbv.setBit(0);
  sbv.setBit(1);
  sbv.setBit(4);
  sbv.setBit(6);
  sbv.setBit(9);
  SparseBitVect sbv2 = sbv;

  REQUIRE_THAT(TanimotoSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(1., 1e-4));
  REQUIRE_THAT(DiceSimilarity(sbv, sbv2), Catch::Matchers::WithinAbs(1., 1e-4));
  REQUIRE_THAT(CosineSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(1., 1e-4));
  REQUIRE_THAT(KulczynskiSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(1., 1e-4));
  REQUIRE_THAT(SokalSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(1., 1e-4));
  REQUIRE_THAT(McConnaugheySimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(1., 1e-4));
  REQUIRE_THAT(BraunBlanquetSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(1., 1e-4));
  REQUIRE_THAT(RusselSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0.5, 1e-4));
  REQUIRE_THAT(RogotGoldbergSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(1., 1e-4));
  REQUIRE_THAT(AllBitSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(1., 1e-4));

  // similarity = 0.0
  sbv2 = SparseBitVect(10);
  sbv2.setBit(2);
  sbv2.setBit(3);
  sbv2.setBit(5);
  sbv2.setBit(7);
  sbv2.setBit(8);

  REQUIRE_THAT(TanimotoSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0., 1e-4));
  REQUIRE_THAT(DiceSimilarity(sbv, sbv2), Catch::Matchers::WithinAbs(0., 1e-4));
  REQUIRE_THAT(CosineSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0., 1e-4));
  REQUIRE_THAT(KulczynskiSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0., 1e-4));
  REQUIRE_THAT(SokalSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0., 1e-4));
  REQUIRE_THAT(McConnaugheySimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(-1., 1e-4));
  REQUIRE_THAT(BraunBlanquetSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0., 1e-4));
  REQUIRE_THAT(RusselSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0., 1e-4));
  REQUIRE_THAT(RogotGoldbergSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0., 1e-4));
  REQUIRE_THAT(AllBitSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0., 1e-4));

  // similarity ~= 0.5
  sbv.setBit(5);
  sbv2 = SparseBitVect(10);
  sbv2.setBit(0);
  sbv2.setBit(2);
  sbv2.setBit(4);
  sbv2.setBit(5);
  sbv2.setBit(8);
  sbv2.setBit(9);

  REQUIRE_THAT(TanimotoSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0.5, 1e-4));
  REQUIRE_THAT(DiceSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0.6666, 1e-4));
  REQUIRE_THAT(CosineSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0.6666, 1e-4));
  REQUIRE_THAT(KulczynskiSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0.6666, 1e-4));
  REQUIRE_THAT(SokalSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0.3333, 1e-4));
  REQUIRE_THAT(McConnaugheySimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0.3333, 1e-4));
  REQUIRE_THAT(BraunBlanquetSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0.6666, 1e-4));
  REQUIRE_THAT(RusselSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0.4, 1e-4));
  REQUIRE_THAT(RogotGoldbergSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0.5833, 1e-4));
  REQUIRE_THAT(AllBitSimilarity(sbv, sbv2),
               Catch::Matchers::WithinAbs(0.6, 1e-4));
}

TEST_CASE("test13BitVectAllOnes") {
  {
    ExplicitBitVect bv(32, false);
    REQUIRE(bv.getNumOnBits() == 0);
    REQUIRE(!bv[0]);

    ExplicitBitVect bv2(32, true);
    REQUIRE(bv2.getNumOnBits() == 32);
    REQUIRE(bv2[0]);
  }
}

TEST_CASE("test14BitVectConcatenation") {
  {
    ExplicitBitVect bv(32, false);
    ExplicitBitVect bv2(32, true);
    ExplicitBitVect bv3 = bv + bv2;
    REQUIRE(bv3.getNumBits() == 64);
    REQUIRE(bv3.getNumOnBits() == 32);
    REQUIRE(bv3.getNumOffBits() == 32);
  }
}

TEST_CASE("test15BitmapOps") {
  {
    const unsigned char bv1[5] = {0x1, 0x1, 0x1, 0x1, 0x1};
    REQUIRE(CalcBitmapPopcount(bv1, 5) == 5);
  }
  {
    const unsigned char bv1[5] = {0x1, 0x1, 0x3, 0x1, 0x1};
    REQUIRE(CalcBitmapPopcount(bv1, 5) == 6);
  }
  {
    const unsigned char bv1[5] = {0x0, 0x0, 0x0, 0x0, 0x0};
    REQUIRE(CalcBitmapPopcount(bv1, 5) == 0);
  }
  {
    const unsigned char bv1[5] = {0x1, 0x1, 0x1, 0x1, 0x1};
    const unsigned char bv2[5] = {0x1, 0x1, 0x1, 0x1, 0x1};
    REQUIRE_THAT(CalcBitmapTanimoto(bv1, bv2, 5),
                 Catch::Matchers::WithinAbs(1.0, 1e-4));
  }
  {
    const unsigned char bv1[5] = {0x1, 0x1, 0x1, 0x1, 0x1};
    const unsigned char bv2[5] = {0x1, 0x1, 0x1, 0x0, 0x1};
    REQUIRE_THAT(CalcBitmapTanimoto(bv1, bv2, 5),
                 Catch::Matchers::WithinAbs(0.8, 1e-4));
  }
  {
    const unsigned char bv1[5] = {0x1, 0x1, 0x1, 0x1, 0x1};
    const unsigned char bv2[5] = {0x1, 0x1, 0x1, 0x1, 0x1};
    REQUIRE_THAT(CalcBitmapTversky(bv1, bv2, 5, 1., 1.),
                 Catch::Matchers::WithinAbs(1.0, 1e-4));
  }
  {
    const unsigned char bv1[5] = {0x1, 0x1, 0x1, 0x1, 0x1};
    const unsigned char bv2[5] = {0x1, 0x1, 0x1, 0x0, 0x1};
    REQUIRE_THAT(CalcBitmapTversky(bv1, bv2, 5, 1., 1.),
                 Catch::Matchers::WithinAbs(0.8, 1e-4));
  }
  {
    const unsigned char bv1[5] = {0x1, 0x1, 0x1, 0x1, 0x1};
    const unsigned char bv2[5] = {0x1, 0x1, 0x1, 0x1, 0x1};
    REQUIRE_THAT(CalcBitmapDice(bv1, bv2, 5),
                 Catch::Matchers::WithinAbs(1.0, 1e-4));
  }
  {
    const unsigned char bv1[5] = {0x1, 0x1, 0x1, 0x1, 0x1};
    const unsigned char bv2[5] = {0x1, 0x1, 0x1, 0x0, 0x1};
    REQUIRE_THAT(CalcBitmapDice(bv1, bv2, 5),
                 Catch::Matchers::WithinAbs(8. / 9, 1e-4));
  }
  {
    const unsigned char bv1[5] = {0x1, 0x1, 0x1, 0x1, 0x1};
    const unsigned char bv2[5] = {0x1, 0x1, 0x1, 0x1, 0x1};
    REQUIRE_THAT(CalcBitmapTversky(bv1, bv2, 5, 0.5, 0.5),
                 Catch::Matchers::WithinAbs(1.0, 1e-4));
  }
  {
    const unsigned char bv1[5] = {0x1, 0x1, 0x1, 0x1, 0x1};
    const unsigned char bv2[5] = {0x1, 0x1, 0x1, 0x0, 0x1};
    REQUIRE_THAT(CalcBitmapTversky(bv1, bv2, 5, 0.5, 0.5),
                 Catch::Matchers::WithinAbs(8. / 9, 1e-4));
  }

  {
    const unsigned char bv1[5] = {0x1, 0x0, 0x1, 0x0, 0x1};
    const unsigned char bv2[5] = {0x1, 0x1, 0x5, 0x1, 0x1};
    REQUIRE(CalcBitmapAllProbeBitsMatch(bv1, bv2, 5));
  }
  {
    const unsigned char bv1[5] = {0x1, 0x0, 0x1, 0x0, 0x1};
    const unsigned char bv2[5] = {0x1, 0x1, 0x2, 0x1, 0x1};
    REQUIRE(!CalcBitmapAllProbeBitsMatch(bv1, bv2, 5));
  }

  {
    const unsigned char bv1[5] = {0x1, 0x0, 0x1, 0x0, 0x1};
    const unsigned char bv2[5] = {0x1, 0x1, 0x5, 0x1, 0x1};
    REQUIRE_THAT(CalcBitmapTversky(bv1, bv2, 5, 1.0, 0.0),
                 Catch::Matchers::WithinAbs(1.0, 1e-4));
  }
}

TEST_CASE("test16BitVectProps") {
  ExplicitBitVect bv(32);
  for (int i = 0; i < 32; i += 2) {
    bv.setBit(i);
  }

  ExplicitBitVect bv2(bv.toString());
  REQUIRE(bv == bv2);

  Dict d;
  d.setVal<ExplicitBitVect>("exp", bv);
  const RDValue &value = d.getRawVal("exp");

  DataStructsExplicitBitVecPropHandler bv_handler;
  std::vector<CustomPropHandler *> handlers = {&bv_handler, bv_handler.clone()};
  for (auto handler : handlers) {
    REQUIRE(handler->canSerialize(value));
    RDValue bad_value = 1;
    REQUIRE(!handler->canSerialize(bad_value));
    std::stringstream ss;
    REQUIRE(handler->write(ss, value));
    RDValue newValue;
    REQUIRE(handler->read(ss, newValue));
    REQUIRE(from_rdvalue<ExplicitBitVect>(newValue) == bv);
    newValue.destroy();
  }
  delete handlers[1];
}

TEST_CASE("test17Github3994") {
  SparseIntVect<std::uint32_t> siv1(128);
  siv1.setVal(0, 3);
  siv1.setVal(100, 4);

  SparseIntVect<std::uint32_t> siv2(128);
  siv2.setVal(1, 12);
  siv2.setVal(17, 3);
  siv2.setVal(99, 4);

  auto siv3 = siv1;
  REQUIRE(siv3.getNonzeroElements().size() == 2);
  REQUIRE(siv3 == siv1);
  REQUIRE(siv3 != siv2);

  siv1 = siv2;
  REQUIRE(siv1.getNonzeroElements().size() == 3);
  REQUIRE(siv1 != siv3);
  REQUIRE(siv1 == siv2);
}

TEST_CASE("test14RealVect") {
  RealValueVect vect1(30);
  unsigned int i;
  for (i = 0; i < 15; ++i) {
    vect1.setVal(2 * i, 1.0);
  }

  REQUIRE(vect1.getLength() == 30);
  REQUIRE_THAT(vect1.getTotalVal(), Catch::Matchers::WithinAbs(15.0, 1e-6));
  for (i = 0; i < vect1.getLength(); ++i) {
    REQUIRE_THAT(vect1.getVal(i),
                 Catch::Matchers::WithinAbs((i + 1) % 2, 1e-6));
  }
  RealValueVect vect2(30);
  for (i = 0; i < vect2.getLength(); ++i) {
    vect2.setVal(i, double(1.0 / (i + 1.0)));
  }

  REQUIRE(vect2.getLength() == 30);
  for (i = 0; i < vect2.getLength(); ++i) {
    REQUIRE_THAT(vect2.getVal(i),
                 Catch::Matchers::WithinAbs(double(1.0 / (i + 1.0)), 1e-6));
  }

  // test copy constructor and operator[]
  RealValueVect vect3(vect2);
  REQUIRE(vect3.getLength() == 30);
  for (i = 0; i < vect3.getLength(); ++i) {
    REQUIRE_THAT(vect3[i],
                 Catch::Matchers::WithinAbs(double(1.0 / (i + 1.0)), 1e-6));
  }

  double Pi = 3.141592;
  RealValueVect vect4(60, Pi);
  REQUIRE_THAT(vect4.getTotalVal(), Catch::Matchers::WithinAbs(60 * Pi, 1e-6));
}

TEST_CASE("test15RealVectDists") {
  RealValueVect v1(30);
  RealValueVect v2(30);
  unsigned int i;
  for (i = 0; i < 15; ++i) {
    v1.setVal(2 * i, 1.0);
    v2.setVal(2 * i, 1.0);
  }
  REQUIRE_THAT(computeL1Norm(v1, v2), Catch::Matchers::WithinAbs(0, 1e-6));
  for (i = 0; i < 30; ++i) {
    v2.setVal(i, i % 2);
  }

  REQUIRE_THAT(computeL1Norm(v1, v2), Catch::Matchers::WithinAbs(30.0, 1e-6));

  for (i = 0; i < 30; ++i) {
    if (i % 3 == 0) {
      v2.setVal(i, 1.0);
    } else {
      v2.setVal(i, 0.0);
    }
  }

  REQUIRE_THAT(computeL1Norm(v1, v2), Catch::Matchers::WithinAbs(15.0, 1e-6));

  for (i = 0; i < 30; ++i) {
    v1.setVal(i, 0.0);
    v2.setVal(i, i / 10.0);
  }

  REQUIRE_THAT(computeL1Norm(v1, v2), Catch::Matchers::WithinAbs(43.5, 1e-6));
}

TEST_CASE("test16RealVectPickles") {
  RealValueVect v1(30);
  unsigned int i;
  for (i = 0; i < 15; ++i) {
    v1.setVal(2 * i, 1.1);
  }
  RealValueVect v2(v1.toString());
  CHECK_THAT(computeL1Norm(v1, v2), Catch::Matchers::WithinAbs(0.0, 0.001));
}

TEST_CASE("test17RealVectOps") {
  RealValueVect vect1(8);
  for (unsigned int i = 0; i < 4; ++i) {
    vect1.setVal(2 * i, 2.1);
  }
  REQUIRE(vect1.getLength() == 8);
  CHECK_THAT(vect1.getTotalVal(), Catch::Matchers::WithinAbs(8.4, 0.001));

  RealValueVect vect2(8);
  for (unsigned int i = 0; i < 4; ++i) {
    vect2.setVal(2 * i + 1, 2.2);
    vect2.setVal(2 * i, 1.1);
  }
  CHECK_THAT(vect2.getTotalVal(), Catch::Matchers::WithinAbs(13.2, 0.001));

  RealValueVect vect3 = vect1 & vect2;
  REQUIRE(vect3.getLength() == 8);
  CHECK_THAT(vect3.getTotalVal(), Catch::Matchers::WithinAbs(4.4, 0.001));

  RealValueVect vect4 = vect1 | vect2;
  REQUIRE(vect4.getLength() == 8);
  CHECK_THAT(vect4.getTotalVal(), Catch::Matchers::WithinAbs(17.2, 0.001));

  RealValueVect vect5 = vect1 + vect2;
  REQUIRE(vect5.getLength() == 8);
  CHECK_THAT(vect5.getTotalVal(), Catch::Matchers::WithinAbs(21.6, 0.001));

  vect5 = vect1 - vect2;
  CHECK_THAT(vect5.getTotalVal(), Catch::Matchers::WithinAbs(-4.8, 0.001));

  vect5 = vect2 - vect1;
  CHECK_THAT(vect5.getTotalVal(), Catch::Matchers::WithinAbs(4.8, 0.001));
}

TEST_CASE("old main") {
  try {
    throw IndexErrorException(3);
  } catch (IndexErrorException &) {
    BOOST_LOG(rdInfoLog) << "pass" << endl;
  }

  stringstream ss(ios_base::binary | ios_base::out | ios_base::in);
  int v1 = 4, v2 = 5, v3, v4;

  ss.write((const char *)&v1, sizeof(v1));
  ss.write((const char *)&v2, sizeof(v2));
  ss.seekp(0, ios_base::beg);
  RDKit::streamRead(ss, v3);
  RDKit::streamRead(ss, v4);

  TXTMSG("v3", v3);
  TXTMSG("v4", v4);

  SECTION("sparse") {
    SparseBitVect sparseFoo(10);
    Test(sparseFoo);
    TaniTest(sparseFoo);
    ProbeTest(sparseFoo);
  }
  SECTION("explicit") {
    ExplicitBitVect explicitFoo(10);
    Test(explicitFoo);
    TaniTest(explicitFoo);
  }
}
