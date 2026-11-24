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
        testCases{{"C1CC2CCC1C2", {{{0, 1}, {5, 7}}}},
                  {"C1CC1C1CC2CCC1C2", {{{1, 2}, {8, 11}}}},
                  {"C1CCCC2C1CCC2", {}},
                  {"C1CC2CCC1C23CCC3", {{{1, 2}, {5, 10}}}}};
    for (const auto &[smiles, expected] : testCases) {
      INFO(smiles);
      auto mol = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(mol);
      auto ringPairs = RDKit::MolOps::details::findRingPairs(*mol);
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
    //          = {
    //     {"C1CC1C1CC2CCC1C2", {{{1, 2}, std::vector<unsigned int>{8, 11}}}},
    //     {"C1CCCC2C1CCC2", {{{0, 1}, std::vector<unsigned int>{4}}}},
    // };
    for (const auto &[smiles, expected] : testCases) {
      INFO(smiles);
      auto mol = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(mol);
      auto ringPairs = RDKit::MolOps::details::findRingPairs(*mol, 1);
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
        };
    for (const auto &[smiles, expected] : testCases) {
      INFO(smiles);
      auto mol = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(mol);
      auto bridgeheads = MolOps::findBridgeheads(*mol);
      std::cerr << "Bridgeheads, N=" << bridgeheads.size() << std::endl;
      for (const auto &bh : bridgeheads) {
        std::cerr << "\t" << bh << std::endl;
      }
      std::cerr << std::endl;
      CHECK(bridgeheads == expected);
      if (!bridgeheads.empty()) {
        std::cerr << "!!!!!" << bridgeheads[0] << std::endl;
      }
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
}