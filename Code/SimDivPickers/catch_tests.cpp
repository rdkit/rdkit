//
//  Copyright (C) 2019
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"
#include <RDGeneral/types.h>
#include <RDGeneral/test.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include <SimDivPickers/LeaderPicker.h>

#include <iostream>
#include <fstream>

template <typename T>
class BVFunctor {
 public:
  BVFunctor(const T &obj) : d_obj(obj) {}
  ~BVFunctor() = default;
  double operator()(unsigned int i, unsigned int j) {
    double res = 1. - TanimotoSimilarity(*d_obj[i], *d_obj[j]);
    return res;
  }
  const T &d_obj;
};

TEST_CASE(
    "Leader Picker basics"
    "[LeaderPicker]") {
  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/SimDivPickers/Wrap/test_data/chembl_cyps.head.fps";
  std::ifstream inf(fName);
  std::string fpsText;
  std::getline(inf, fpsText);
  std::vector<std::unique_ptr<ExplicitBitVect>> fps;
  while (!inf.eof() && !fpsText.empty()) {
    fps.emplace_back(new ExplicitBitVect(fpsText.size() * 4));
    UpdateBitVectFromFPSText(*fps.back(), fpsText);
    std::getline(inf, fpsText);
  };
  REQUIRE(fps.size() == 1000);
  BVFunctor<std::vector<std::unique_ptr<ExplicitBitVect>>> bvf(fps);
  RDPickers::LeaderPicker pkr;
  SECTION("basics1") {
    double threshold = 0.8;
    auto res = pkr.lazyPick(bvf, fps.size(), 0, threshold);
    CHECK(res.size() == 146);
    for (unsigned i = 0; i < res.size(); ++i) {
      for (unsigned j = 0; j < i; ++j) {
        CHECK(bvf(res[i], res[j]) >= threshold);
      }
    }
  }
  SECTION("basics2") {
    double threshold = 0.9;
    auto res = pkr.lazyPick(bvf, fps.size(), 0, threshold);
    CHECK(res.size() == 14);
    for (unsigned i = 0; i < res.size(); ++i) {
      for (unsigned j = 0; j < i; ++j) {
        CHECK(bvf(res[i], res[j]) >= threshold);
      }
    }
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  SECTION("basics multithreaded") {
    double threshold = 0.8;
    RDKit::INT_VECT firstPicks;
    int nThreads = 0;  // use max available
    auto res =
        pkr.lazyPick(bvf, fps.size(), 0, firstPicks, threshold, nThreads);
    CHECK(res.size() == 146);
    for (unsigned i = 0; i < res.size(); ++i) {
      for (unsigned j = 0; j < i; ++j) {
        CHECK(bvf(res[i], res[j]) >= threshold);
      }
    }
  }
#endif
}
