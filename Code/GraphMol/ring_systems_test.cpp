//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <GraphMol/test_fixtures.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("missing rings") {
  SECTION("case 1 - RDL") {
    auto m = "O=C=NC1=CC2C3=C(C=C1)C2=C(N=C=O)C=C3"_smiles;
    REQUIRE(m);
    MolOps::findRingFamilies(*m);
    auto arfs = m->getRingInfo()->atomRingFamilies();
    REQUIRE(arfs.size() == 3);
    CHECK(arfs[0] == std::vector<int>{5, 6, 7, 10});
    CHECK(arfs[1] == std::vector<int>{5, 6, 7, 10, 11, 15, 16});
    CHECK(arfs[2] == std::vector<int>{3, 4, 5, 6, 7, 8, 9, 10});
    auto arcs = m->getRingInfo()->atomRelevantCycles();
    REQUIRE(arcs.size() == 5);
    std::vector<std::vector<int>> expected{{5, 6, 7, 10},
                                           {5, 6, 16, 15, 11, 10},
                                           {6, 7, 10, 11, 15, 16},
                                           {3, 4, 5, 6, 7, 8, 9},
                                           {3, 4, 5, 10, 7, 8, 9}};
    CHECK(arcs == expected);
  }
  SECTION("case 1 - FindRings") {
    auto m = "O=C=NC1=CC2C3=C(C=C1)C2=C(N=C=O)C=C3"_smiles;
    REQUIRE(m);
    MolOps::symmetrizeSSSR(*m);
    auto arings = m->getRingInfo()->atomRings();
    REQUIRE(arings.size() == 5);
    std::vector<std::vector<int>> expected{{5, 6, 7, 10},
                                           {5, 6, 16, 15, 11, 10},
                                           {6, 7, 10, 11, 15, 16},
                                           {3, 4, 5, 6, 7, 8, 9},
                                           {3, 4, 5, 10, 7, 8, 9}};
    CHECK(arings == expected);
  }
}

TEST_CASE("performance") {
  std::string smi =
      "C1(CC3)CCC3CC(CC3)CCC3CC(CC3)CCC3C(CC3)CCC3CC(CC3)CCC3CC(CC3)CCC3C1";
  v2::SmilesParse::SmilesParserParams params;
  params.sanitize = false;
  auto m = v2::SmilesParse::MolFromSmiles(smi, params);
  REQUIRE(m);
  SECTION("case 1 - RDL") {
    std::chrono::high_resolution_clock::time_point t1 =
        std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1000; ++i) {
      MolOps::symmetrizeSSSR(*m);
    }
    std::chrono::high_resolution_clock::time_point t2 =
        std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "Time taken for 1000 calls to symmetrizeSSSR (new): "
              << time_span.count() << " seconds." << std::endl;
    std::cout << "mol has " << m->getRingInfo()->numRings() << " rings"
              << std::endl;
  }
  SECTION("case 2 - legacy") {
    std::chrono::high_resolution_clock::time_point t1 =
        std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1000; ++i) {
      bool includeDative = true;
      bool includeHBonds = true;
      bool legacyCalculation = true;
      MolOps::symmetrizeSSSR(*m, includeDative, includeHBonds,
                             legacyCalculation);
    }
    std::chrono::high_resolution_clock::time_point t2 =
        std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "Time taken for 1000 calls to symmetrizeSSSR (legacy): "
              << time_span.count() << " seconds." << std::endl;
    std::cout << "mol has " << m->getRingInfo()->atomRings().size() << " rings"
              << std::endl;
  }
}
