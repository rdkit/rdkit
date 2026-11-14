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
#include <RDStreams/streams.h>

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
    const auto arings = m->getRingInfo()->atomRings();
    REQUIRE(arings.size() == 5);
    std::vector<std::vector<int>> expected{{5, 6, 7, 10},
                                           {5, 6, 16, 15, 11, 10},
                                           {6, 7, 10, 11, 15, 16},
                                           {3, 4, 5, 6, 7, 8, 9},
                                           {3, 4, 5, 10, 7, 8, 9}};
    REQUIRE(arings.size() == expected.size());
    for (size_t i = 0; i < arings.size(); ++i) {
      CHECK(std::ranges::equal(arings[i], expected[i]));
    }
  }
}

TEST_CASE("environment variable") {
  std::string smi =
      "C1(CC3)CCC3CC(CC3)CCC3CC(CC3)CCC3C(CC3)CCC3CC(CC3)CCC3CC(CC3)CCC3C1";
  SECTION("legacy") {
    UseLegacyRingFindingFixture fix(true);
    auto m = v2::SmilesParse::MolFromSmiles(smi);
    REQUIRE(m);
    CHECK(m->getRingInfo()->numRings() == 24);
  }
  SECTION("RDL") {
    UseLegacyRingFindingFixture fix(false);
    auto m = v2::SmilesParse::MolFromSmiles(smi);
    REQUIRE(m);
    CHECK(m->getRingInfo()->numRings() == 70);
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
      MolOps::symmetrizeSSSR(*m, MolOps::SymmetrizeSSSRAlgorithm::LEGACY);
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

TEST_CASE("performance 2") {
  constexpr const char *smi =
      "COc1cc2ccc1n1c(=O)n(c3ccc(cc3OC)c3ccc(c(c3)OC)n3c(=O)n(c4ccc(cc4OC)"
      "c4ccc(c(c4)OC)n4c(=O)n(c5ccc(cc5OC)c5ccc(n6c(=O)n(c7ccc(c8ccc(n9c(=O)"
      "n(c%10ccc(c%11ccc(n%12c(=O)n(c%13ccc2cc%13OC)c(=O)n(c%12=O)c2ccc("
      "cc2OC)C)c(OC)c%11)cc%10OC)c(=O)n(c9=O)c2ccc(cc2OC)C)c(OC)c8)cc7OC)c(="
      "O)n(c6=O)c2ccc(cc2OC)C)c(c5)OC)c(=O)n(c4=O)c2ccc(cc2OC)C)c(=O)n(c3=O)"
      "c2ccc(cc2OC)C)c(=O)n(c1=O)c1ccc(cc1OC)C";
  v2::SmilesParse::SmilesParserParams params;
  params.sanitize = false;
  params.removeHs = false;
  auto m = v2::SmilesParse::MolFromSmiles(smi, params);
  REQUIRE(m);
  SECTION("case 1 - RDL") {
    std::chrono::high_resolution_clock::time_point t1 =
        std::chrono::high_resolution_clock::now();
    MolOps::symmetrizeSSSR(*m);
    std::chrono::high_resolution_clock::time_point t2 =
        std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "Time taken for 1 call to symmetrizeSSSR (new): "
              << time_span.count() << " seconds." << std::endl;
    std::cout << "mol has " << m->getRingInfo()->numRings() << " rings"
              << std::endl;
  }
  SECTION("case 2 - legacy") {
    std::chrono::high_resolution_clock::time_point t1 =
        std::chrono::high_resolution_clock::now();
    MolOps::symmetrizeSSSR(*m, MolOps::SymmetrizeSSSRAlgorithm::LEGACY);
    std::chrono::high_resolution_clock::time_point t2 =
        std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "Time taken for 1 call to symmetrizeSSSR (legacy): "
              << time_span.count() << " seconds." << std::endl;
    std::cout << "mol has " << m->getRingInfo()->atomRings().size() << " rings"
              << std::endl;
  }
}

#ifdef RDK_USE_BOOST_IOSTREAMS

TEST_CASE("bulk performance") {
  std::string rdbase = getenv("RDBASE");
  auto path = rdbase + "/Regress/Data/chembl_20_chiral.smi.gz";
  // auto path = rdbase + "/Regress/Data/znp.50k.smi.gz";
  std::unique_ptr<std::istream> strm = std::make_unique<gzstream>(path);
  double taccum_legacy = 0.0;
  double taccum_updated = 0.0;
  unsigned int nToDo = 10000;
  for (unsigned int i = 0; i < nToDo; ++i) {
    std::string line;
    std::getline(*strm, line);
    v2::SmilesParse::SmilesParserParams params;
    params.sanitize = false;
    params.removeHs = false;

    {
      auto m = v2::SmilesParse::MolFromSmiles(line, params);
      REQUIRE(m);
      std::chrono::high_resolution_clock::time_point t1 =
          std::chrono::high_resolution_clock::now();
      MolOps::symmetrizeSSSR(*m, MolOps::SymmetrizeSSSRAlgorithm::LEGACY);
      std::chrono::high_resolution_clock::time_point t2 =
          std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
      taccum_legacy += time_span.count();
    }
    {
      auto m = v2::SmilesParse::MolFromSmiles(line, params);
      REQUIRE(m);
      std::chrono::high_resolution_clock::time_point t1 =
          std::chrono::high_resolution_clock::now();
      MolOps::symmetrizeSSSR(*m);
      std::chrono::high_resolution_clock::time_point t2 =
          std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
      taccum_updated += time_span.count();
    }
  }
  std::cout << "Time taken for " << nToDo
            << " calls to symmetrizeSSSR (legacy): " << taccum_legacy
            << " seconds." << std::endl;
  std::cout << "Time taken for " << nToDo
            << " calls to symmetrizeSSSR (new): " << taccum_updated
            << " seconds." << std::endl;
}
#endif

TEST_CASE("max number of relevant cycles") {
  SECTION("clipping") {
    std::string unit = "-C1CCC(CC1)";
    std::string smiles = "C12CCC(CC1)";
    for (auto i = 0; i < 20; ++i) {
      smiles += unit;
    }
    smiles += "-C1CCC2(CC1)";
    v2::SmilesParse::SmilesParserParams params;
    params.sanitize = false;
    auto m = v2::SmilesParse::MolFromSmiles(smiles, params);
    REQUIRE(m);
    MolOps::findRingFamilies(*m);
    auto ri = m->getRingInfo();
    CHECK(ri->numRingFamilies() == 23);
    CHECK(ri->numRelevantCycles() == 4194326);
    std::vector<std::vector<int>> rings;
    auto nExpected = MAX_SYMMETRIZED_RINGS + ri->numRingFamilies();
    CHECK(MolOps::symmetrizeSSSR(*m, rings,
                                 MolOps::SymmetrizeSSSRAlgorithm::RDL) ==
          static_cast<int>(nExpected));
    CHECK(rings.size() == nExpected);
  }
  SECTION("max-out") {
    std::string unit = "-C1CCC(CC1)";
    std::string smiles = "C12CCC(CC1)";
    for (auto i = 0; i < 31; ++i) {
      smiles += unit;
    }
    smiles += "-C1CCC2(CC1)";
    v2::SmilesParse::SmilesParserParams params;
    params.sanitize = false;
    auto m = v2::SmilesParse::MolFromSmiles(smiles, params);
    REQUIRE(m);
    MolOps::findRingFamilies(*m);
    auto ri = m->getRingInfo();
    CHECK(ri->numRingFamilies() == 34);
    CHECK_THROWS_AS(ri->numRelevantCycles(), ValueErrorException);
    std::vector<std::vector<int>> rings;
    auto nExpected = MAX_SYMMETRIZED_RINGS + ri->numRingFamilies();
    CHECK(MolOps::symmetrizeSSSR(*m, rings,
                                 MolOps::SymmetrizeSSSRAlgorithm::RDL) ==
          static_cast<int>(nExpected));
    CHECK(rings.size() == nExpected);
  }
  SECTION("many individual rings") {
    std::string unit = "C1CC1.";
    std::string smiles = "";
    unsigned int nRings = MAX_SYMMETRIZED_RINGS + 1;
    for (auto i = 0u; i < nRings; ++i) {
      smiles += unit;
    }
    smiles += "C1CC1";
    v2::SmilesParse::SmilesParserParams params;
    params.sanitize = false;
    auto m = v2::SmilesParse::MolFromSmiles(smiles, params);
    REQUIRE(m);
    MolOps::findRingFamilies(*m);
    auto ri = m->getRingInfo();
    CHECK(ri->numRingFamilies() == nRings + 1);
    CHECK(ri->numRelevantCycles() == ri->numRingFamilies());
    std::vector<std::vector<int>> rings;
    CHECK(MolOps::symmetrizeSSSR(*m, rings,
                                 MolOps::SymmetrizeSSSRAlgorithm::RDL) ==
          static_cast<int>(nRings + 1));
    CHECK(rings.size() == nRings + 1);
  }
  SECTION("many individual rings + symmetry") {
    std::string unit = "C1CC1.";
    std::string smiles = "";
    unsigned int nRings = MAX_SYMMETRIZED_RINGS + 1;
    for (auto i = 0u; i < nRings; ++i) {
      smiles += unit;
    }
    smiles += "C1CC1";
    smiles += ".C12CCC(CC1)-C1CCC(CC1)-C1CCC2(CC1)";
    size_t nExtraRings = 4;
    size_t nExtraSymmetric = 7;
    v2::SmilesParse::SmilesParserParams params;
    params.sanitize = false;
    auto m = v2::SmilesParse::MolFromSmiles(smiles, params);
    REQUIRE(m);
    MolOps::findRingFamilies(*m);
    auto ri = m->getRingInfo();
    CHECK(ri->numRingFamilies() == nRings + 1 + nExtraRings);
    CHECK(ri->numRelevantCycles() == ri->numRingFamilies() + nExtraSymmetric);
    std::vector<std::vector<int>> rings;
    CHECK(MolOps::symmetrizeSSSR(*m, rings,
                                 MolOps::SymmetrizeSSSRAlgorithm::RDL) ==
          static_cast<int>(nRings + 1 + nExtraRings + nExtraSymmetric));
    CHECK(rings.size() == nRings + 1 + nExtraRings + nExtraSymmetric);
  }
  SECTION("SSSR rings are returned when clipping occurs") {
    std::string unit = "-C1CCC(CC1)";
    std::string smiles = "C12CCC(CC1)";
    for (auto i = 0; i < 20; ++i) {
      smiles += unit;
    }
    smiles += "-C1CCC2(CC1).C1CCCCCCC1";
    v2::SmilesParse::SmilesParserParams params;
    params.sanitize = false;
    auto m = v2::SmilesParse::MolFromSmiles(smiles, params);
    REQUIRE(m);
    std::vector<std::vector<int>> rings;
    CHECK(MolOps::findSSSR(*m, rings) == 24);
    CHECK(rings.size() == 24);
    for (auto i = 0; i < 22; ++i) {
      CHECK(rings[i].size() == 6);
    }
    CHECK(rings[22].size() == 8);
    MolOps::findRingFamilies(*m);
    auto ri = m->getRingInfo();
    auto nExpected = MAX_SYMMETRIZED_RINGS + ri->numRingFamilies();
    CHECK(MolOps::symmetrizeSSSR(*m, rings,
                                 MolOps::SymmetrizeSSSRAlgorithm::RDL) ==
          static_cast<int>(nExpected));
    CHECK(rings.size() == nExpected);
    for (auto i = 0; i < 22; ++i) {
      CHECK(rings[i].size() == 6);
    }
    CHECK(rings[22].size() == 8);
  }
}

TEST_CASE("rings from molecules with fragments are also ordered") {
  auto m = "C1CCCCC1.C1CC1.C1CCCC1.C1CCC1"_smiles;
  REQUIRE(m);

  auto checkit = [](const auto &rings) {
    REQUIRE(rings.size() == 4);
    CHECK(rings[0].size() == 3);
    CHECK(rings[1].size() == 4);
    CHECK(rings[2].size() == 5);
    CHECK(rings[3].size() == 6);
  };
  std::vector<std::vector<int>> rings;
  INFO("SSSR");
  CHECK(MolOps::findSSSR(*m, rings) == 4);
  checkit(rings);
  INFO("Symmetrized SSSR");
  CHECK(MolOps::symmetrizeSSSR(*m, rings) == 4);
  checkit(rings);
}