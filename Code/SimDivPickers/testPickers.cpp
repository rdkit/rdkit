//
//  Copyright (C) 2017 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include "MaxMinPicker.h"
#include <iostream>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>

namespace {
double dist_on_line(unsigned int i, unsigned int j) {
  return std::fabs((double)i - (double)j);
}
}  // namespace
void testGithub1421() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "Testing github issue 1421: MaxMinPicker picking non-existent element."
      << std::endl;
  RDPickers::MaxMinPicker pkr;
  RDKit::INT_VECT picks;
  int poolSz = 1000;
  picks = pkr.lazyPick(dist_on_line, poolSz, 10, RDKit::INT_VECT(), 2748);
  for (auto pick : picks) {
    TEST_ASSERT(pick < poolSz);
  }
  BOOST_LOG(rdErrorLog) << "Done" << std::endl;
}

void testGithub2245() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing github issue 2245: MinMax Diversity picker "
                           "seeding shows deterministic / non-random behaviour."
                        << std::endl;
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
    TEST_ASSERT(maxAllowedFailures);
  }
  {  // make sure the default is also random
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
    TEST_ASSERT(maxAllowedFailures);
  }
  {  // and we're still reproducible when we want to be
    RDPickers::MaxMinPicker pkr;
    int poolSz = 1000;
    auto picks1 =
        pkr.lazyPick(dist_on_line, poolSz, 10, RDKit::INT_VECT(), 0xf00d);
    auto picks2 =
        pkr.lazyPick(dist_on_line, poolSz, 10, RDKit::INT_VECT(), 0xf00d);
    TEST_ASSERT(picks1 == picks2);
  }
  BOOST_LOG(rdErrorLog) << "Done" << std::endl;
}

int main() {
  RDLog::InitLogs();
  testGithub1421();
  testGithub2245();
  return 0;
}
