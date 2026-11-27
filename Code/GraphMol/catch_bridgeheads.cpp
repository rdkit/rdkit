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
#include <vector>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Bridgeheads.h>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("details") {
  SECTION("Find ring pairs") {
    std::vector<
        std::pair<std::string,
                  std::vector<std::pair<std::pair<unsigned int, unsigned int>,
                                        std::vector<unsigned int>>>>>
        testCases{
            {"C1CC2CCC1C2", {{{0, 1}, {5, 7}}}},
            {"C1CC1C1CC2CCC1C2", {{{1, 2}, {8, 11}}}},
            {"C1CCCC2C1CCC2", {}},
            {"C1CC2CCC1C23CCC3", {{{1, 2}, {5, 10}}}},
            {"C1CC2CCC1C2-C1CC2CCC1C2",
             {
                 {{0, 1}, {5, 15}},
                 {{2, 3}, {12, 16}},
             }},
        };
    for (const auto &[smiles, expected] : testCases) {
      INFO(smiles);
      auto mol = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(mol);
      auto ringPairs = RDKit::MolOps::details::findRingPairs(*mol);
      if (ringPairs != expected) {
        for (const auto &pair : ringPairs) {
          std::cerr << "Ring pair: " << pair.first.first << ","
                    << pair.first.second << std::endl;
          std::cerr << "Atoms: ";
          for (auto atom : pair.second) {
            std::cerr << atom << " ";
          }
          std::cerr << std::endl;
        }
      }
      CHECK(ringPairs == expected);
    }
  }
  SECTION("Find ring pairs, change minCommonBonds") {
    std::vector<
        std::pair<std::string,
                  std::vector<std::pair<std::pair<unsigned int, unsigned int>,
                                        std::vector<unsigned int>>>>>
        testCases{
            {"C1CC1C1CC2CCC1C2", {{{1, 2}, {8, 11}}}},
            {"C1CCCC2C1CCC2", {{{0, 1}, {4}}}},
        };
    for (const auto &[smiles, expected] : testCases) {
      INFO(smiles);
      auto mol = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(mol);
      auto ringPairs = RDKit::MolOps::details::findRingPairs(*mol, 1);
      CHECK(ringPairs == expected);
    }
  }
  SECTION("Find ring pairs, change considered atoms") {
    auto mol = "C1CC2CCC1C2-C1CC2CCC1C2"_smiles;
    REQUIRE(mol);
    {
      auto onlyWithAtoms = std::vector<unsigned int>{5};
      auto expected =
          std::vector<std::pair<std::pair<unsigned int, unsigned int>,
                                std::vector<unsigned int>>>{{{0, 1}, {5, 15}}};
      auto ringPairs =
          RDKit::MolOps::details::findRingPairs(*mol, 2, onlyWithAtoms);
      CHECK(ringPairs == expected);
    }
  }
}

TEST_CASE("Bridgeheads") {
  SECTION("Simple bridgehead detection") {
    std::vector<std::pair<std::string, std::vector<MolOps::Bridgehead>>>
        testCases = {
            // norbornane
            {"C1CC2CCC1C2", {{2, 5, {0, 1, 2, 6, 5}, {2, 3, 4, 5, 6}}}},
            // external atoms do not change the answer
            {"C1CC23CCC14C2.C3.C4", {{2, 5, {0, 1, 2, 6, 5}, {2, 3, 4, 5, 6}}}},
            // adamantane
            {
                "C1C3CC2CC(CC1C2)C3",
                {
                    {1, 7, {0, 1, 2, 3, 8, 7}, {0, 1, 9, 5, 6, 7}},
                    {1, 3, {0, 1, 2, 3, 8, 7}, {1, 2, 3, 4, 5, 9}},
                    {3, 7, {0, 1, 2, 3, 8, 7}, {3, 4, 5, 6, 7, 8}},
                    {1, 5, {0, 1, 9, 5, 6, 7}, {1, 2, 3, 4, 5, 9}},
                    {5, 7, {0, 1, 9, 5, 6, 7}, {3, 4, 5, 6, 7, 8}},
                    {3, 5, {1, 2, 3, 4, 5, 9}, {3, 4, 5, 6, 7, 8}},
                },
            },
            // norbornane with a spiro link on the bridge
            {"C1CC2CCC1C23CCC3", {{2, 5, {0, 1, 2, 6, 5}, {2, 3, 4, 5, 6}}}},
            // two spiro-linked norbornanes
            {"C12(C3CCC1CC3)C1CCC2CC1",
             {{1, 4, {0, 1, 2, 3, 4}, {0, 1, 6, 5, 4}},
              {7, 10, {0, 7, 8, 9, 10}, {0, 7, 12, 11, 10}}}},
            // cubane doesn't have any bridgeheads
            {"C12C3C4C1C5C2C3C45", {}},
        };
    for (const auto &[smiles, expected] : testCases) {
      INFO(smiles);
      auto mol = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(mol);
      auto bridgeheads = MolOps::findBridgeheads(*mol);
      CHECK(bridgeheads == expected);
    }
  }
  SECTION("test onlyWithAtoms parameter") {
    auto m = "C1C3CC2CC(CC1C2)C3"_smiles;
    REQUIRE(m);
    std::vector<unsigned int> onlyWithAtoms = {1, 5};
    auto bridgeheads = MolOps::findBridgeheads(*m, onlyWithAtoms);
    CHECK(bridgeheads.size() == 5);
    for (const auto &bh : bridgeheads) {
      CHECK((bh.atom1 == 1 || bh.atom1 == 5 || bh.atom2 == 1 || bh.atom2 == 5));
    }
  }
}

TEST_CASE("Bridgehead atoms") {
  SECTION("Simple bridgehead detection") {
    std::vector<std::pair<std::string, std::vector<unsigned int>>> testCases = {
        {"C1CC2CCC1C2", {2, 5}},               // norbornane
        {"C1C3CC2CC(CC1C2)C3", {1, 3, 5, 7}},  // adamantane
        {"C1CC2CCC1C23CCC3",
         {2, 5}},  // norbornane with a spiro link on the bridge
    };
    for (const auto &[smiles, expected] : testCases) {
      INFO(smiles);
      auto mol = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(mol);
      auto bridgeheads = MolOps::findBridgeheadAtoms(*mol);
      CHECK(bridgeheads == expected);
    }
  }
  SECTION("test onlyWithAtoms parameter") {
    // adamantane
    auto m = "C1C3CC2CC(CC1C2)C3"_smiles;
    REQUIRE(m);
    std::vector<unsigned int> onlyWithAtoms = {1, 5};
    auto bridgeheads = MolOps::findBridgeheadAtoms(*m, onlyWithAtoms);
    CHECK(bridgeheads == onlyWithAtoms);
  }
  SECTION("using SSSR also works") {
    v2::SmilesParse::SmilesParserParams params{.sanitize = false};
    auto m = v2::SmilesParse::MolFromSmiles("C1CC2CCC1CC2", params);
    REQUIRE(m);
    m->updatePropertyCache();
    MolOps::findSSSR(*m);
    auto bridgeheads = MolOps::findBridgeheadAtoms(*m);
    CHECK(bridgeheads == std::vector<unsigned int>{2, 5});
    MolOps::symmetrizeSSSR(*m);
    bridgeheads = MolOps::findBridgeheadAtoms(*m);
    CHECK(bridgeheads == std::vector<unsigned int>{2, 5});
  }
}