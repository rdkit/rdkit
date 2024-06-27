//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>

using namespace RDKit;

TEST_CASE(
    "Github #7162: Unexpected exact mass values are returned for radium and radon") {
  SECTION("radium") {
    CHECK(PeriodicTable::getTable()->getMostCommonIsotope(88) == 226);
    CHECK(PeriodicTable::getTable()->getMassForIsotope(88, 226) ==
          Catch::Approx(226.02540).epsilon(1e-4));
  }
  SECTION("radon") {
    CHECK(PeriodicTable::getTable()->getMostCommonIsotope(86) == 222);
    CHECK(PeriodicTable::getTable()->getMassForIsotope(86, 222) ==
          Catch::Approx(222.01757).epsilon(1e-4));
  }
}

TEST_CASE("periodic table row") {
  CHECK(PeriodicTable::getTable()->getRow(0) == 0);
  CHECK(PeriodicTable::getTable()->getRow(1) == 1);
  CHECK(PeriodicTable::getTable()->getRow(92) == 7);
}
