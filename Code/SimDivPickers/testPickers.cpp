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
#include <boost/foreach.hpp>

namespace {
double dist_on_line(unsigned int i, unsigned int j) {
  return std::fabs((double)i - (double)j);
}
}
void testGithub1421() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "Testing github issue 1421: MaxMinPicker picking non-existent element."
      << std::endl;
  RDPickers::MaxMinPicker pkr;
  RDKit::INT_VECT picks;
  int poolSz = 1000;
  picks = pkr.lazyPick(dist_on_line, poolSz, 10, RDKit::INT_VECT(), 2748);
  BOOST_FOREACH (int pick, picks) { TEST_ASSERT(pick < poolSz); }
  BOOST_LOG(rdErrorLog) << "Done" << std::endl;
}

int main() {
  RDLog::InitLogs();
  testGithub1421();
  return 0;
}
