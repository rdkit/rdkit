//
//  Copyright (C) 2017-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>
#include "MaxMinPicker.h"
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>

namespace {
double dist_on_line(unsigned int i, unsigned int j) {
  return std::fabs((double)i - (double)j);
}
}  // namespace
TEST_CASE("testGithub1421") {
  RDPickers::MaxMinPicker pkr;
  RDKit::INT_VECT picks;
  int poolSz = 1000;
  picks = pkr.lazyPick(dist_on_line, poolSz, 10, RDKit::INT_VECT(), 2748);
  for (auto pick : picks) {
    REQUIRE(pick < poolSz);
  }
}

TEST_CASE("testGithub2245") {
  const int MAX_ALLOWED_FAILURES = 3;
  int maxAllowedFailures;
  {
    RDPickers::MaxMinPicker pkr;
    int poolSz = 1000;
    auto picks1 = pkr.lazyPick(dist_on_line, poolSz, 10, RDKit::INT_VECT(), -1);
    for (maxAllowedFailures = MAX_ALLOWED_FAILURES; maxAllowedFailures;
         --maxAllowedFailures) {
      auto picks2 =
          pkr.lazyPick(dist_on_line, poolSz, 10, RDKit::INT_VECT(), -1);
      if (picks1 != picks2) {
        break;
      }
    }
    REQUIRE(maxAllowedFailures);
  }
  {
    RDPickers::MaxMinPicker pkr;
    int poolSz = 1000;
    auto picks1 = pkr.lazyPick(dist_on_line, poolSz, 10);
    for (maxAllowedFailures = MAX_ALLOWED_FAILURES; maxAllowedFailures;
         --maxAllowedFailures) {
      auto picks2 = pkr.lazyPick(dist_on_line, poolSz, 10);
      if (picks1 != picks2) {
        break;
      }
    }
    REQUIRE(maxAllowedFailures);
  }
  {
    RDPickers::MaxMinPicker pkr;
    int poolSz = 1000;
    auto picks1 =
        pkr.lazyPick(dist_on_line, poolSz, 10, RDKit::INT_VECT(), 0xf00d);
    auto picks2 =
        pkr.lazyPick(dist_on_line, poolSz, 10, RDKit::INT_VECT(), 0xf00d);
    REQUIRE(picks1 == picks2);
  }
}
