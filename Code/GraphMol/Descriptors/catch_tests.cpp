//
//  Copyright (C) 2021 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Descriptors/ConnectivityDescriptors.h>

using namespace RDKit;

TEST_CASE("Kier kappa2", "[2D]") {
  SECTION("values from https://doi.org/10.1002/qsar.19860050103") {
    std::vector<std::pair<std::string, double>> data = {
      // Table 5 from the paper
      {"c1ccccc15.Cl5", 1.987},
      {"c1ccccc15.F5", 1.735},
      {"c1ccccc15.[N+]5(=O)O", 2.259},
      {"c1ccccc15.C5(=O)C", 2.444},
#if 0
      // expected values from paper (differences are due to hybridization mismatches)
        {"c1ccccc15.N5(C)C", 2.646},
        {"c1ccccc15.C5(=O)N", 2.416},
        {"c1ccccc15.C5(=O)O", 2.416},
        {"c1ccccc15.S5(=O)(=O)C", 2.617},
        {"c1ccccc15.O5", 1.756},
#else
      {"c1ccccc15.N5(C)C", 2.53},
      {"c1ccccc15.C5(=O)N", 2.31},
      {"c1ccccc15.C5(=O)O", 2.31},
      {"c1ccccc15.S5(=O)(=O)C", 2.42},
      {"c1ccccc15.O5", 1.65},
#endif
    };
    for (const auto &pr : data) {
      std::unique_ptr<ROMol> m(SmilesToMol(pr.first));
      REQUIRE(m);
      auto k2 = Descriptors::calcKappa2(*m);
      CHECK(k2 == Approx(pr.second).epsilon(0.01));
    }
  }
}

TEST_CASE("Kier Phi", "[2D]") {
  SECTION("regression-test values from the paper") {
    std::vector<std::pair<std::string, double>> data = {
      // Table 1 from the paper
      {"CCCCCC", 5.00},
      {"CCCCCCC", 6.00},
      {"CCCCCCCC", 7.00},
      {"CCC(C)CC", 3.20},
      {"CC(C)C(C)C", 2.22},
      {"CC(C)(C)CC", 1.63},
      {"C1CCCC1", 0.92},
      {"C1CCCCC1", 1.54},
      {"C1CCCCCC1", 2.25},
      {"CCCCC=C", 4.53},
      {"C=CCCC=C", 4.07},
      {"C#CCCCC", 4.21},
      {"c1ccccc1", 0.91},
      // additional from Table 2
      {"C=CCC=CC", 4.09},
      {"CC=CC=CC", 4.09},
      {"C1=CCCCC1", 1.31},
      {"C1=CC=CCC1", 1.1},
      {"C1=CCC=CC1", 1.1},
      // additional from Table 3
      {"CCCCCCCCCC", 9.00},
      {"CC(C)CCCCC", 5.14},
      {"CCC(C)CCCC", 5.14},
      {"CC(C)CCCC", 4.17},
      {"CCC(C)CCC", 4.17},
      {"CCC(CC)CCC", 5.14},
      {"CCC(CC)CC", 4.17},
      {"CC(C)(C)CCC", 2.34},
      {"CC(C)C(C)CC", 3.06},
      {"CCC(C)(C)CC", 2.34},
      {"CC(C)(C)C(C)C", 1.85},
      // additional from table 4
      {"CCOCC", 3.93},
      {"CCC(=O)CC", 2.73},
      {"CCc1ccc(CC)cc1", 2.49},
      {"CCC(O)CC", 3.14},
      {"CCCC(Cl)(Cl)CCC", 4.69},
      {"CCC(F)C(F)CC", 3.75},
#if 0
      // expected values from paper (differences are due to hybridization mismatches)
        {"CCOC(=O)CC", 3.61},
        {"CCC(=O)Nc1ccc(CC)cc1", 3.65},
#else
      {"CCOC(=O)CC", 3.38},
      {"CCC(=O)Nc1ccc(CC)cc1", 3.50},
#endif

    };
    for (const auto &pr : data) {
      std::unique_ptr<ROMol> m(SmilesToMol(pr.first));
      REQUIRE(m);
      auto val = Descriptors::calcPhi(*m);
      CHECK(val == Approx(pr.second).epsilon(0.01));
    }
  }
}
