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
    std::vector<std::pair<std::string,
                          std::vector<std::pair<unsigned int, unsigned int>>>>
        testCases = {
            {"C1CC2CCC1C2", {{0, 1}}},
            {"C1CC1C1CC2CCC1C2", {{1, 2}}},
            {"C1CCCC2C1CCC2", {}},
        };
    for (const auto &[smiles, expected] : testCases) {
      INFO(smiles);
      auto mol = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(mol);
      auto ringPairs = RDKit::MolOps::details::findRingPairs(*mol);
      CHECK(ringPairs == expected);
    }
  }
  SECTION("Find ring pairs, change minCommonBonds") {
    std::vector<std::pair<std::string,
                          std::vector<std::pair<unsigned int, unsigned int>>>>
        testCases = {
            {"C1CC1C1CC2CCC1C2", {{1, 2}}},
            {"C1CCCC2C1CCC2", {{0, 1}}},
        };
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
            {"C1CC2CCC1C2", {{2, 5, {0, 1, 2, 6, 5}, {2, 3, 4, 5, 6}}}}};
    for (const auto &[smiles, expected] : testCases) {
      INFO(smiles);
      auto mol = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(mol);
      auto bridgeheads = MolOps::findBridgeheads(*mol);
      CHECK(bridgeheads == expected);
      std::cerr << "!!!!!" << bridgeheads[0] << std::endl;
    }
  }
}

TEST_CASE("Bridgehead atoms") {
  SECTION("Simple bridgehead detection") {
    std::vector<std::pair<std::string, std::vector<unsigned int>>> testCases = {
        {"C1CC2CCC1C2", {2, 5}}};
    for (const auto &[smiles, expected] : testCases) {
      INFO(smiles);
      auto mol = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(mol);
      auto bridgeheads = MolOps::findBridgeheadAtoms(*mol);
      CHECK(bridgeheads == expected);
    }
  }
}