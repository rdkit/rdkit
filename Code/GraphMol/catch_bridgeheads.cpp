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

TEST_CASE("Bridgeheads") {
  SECTION("Simple bridgehead detection") {
    std::vector<std::pair<std::string, std::vector<MolOps::Bridgehead>>>
        testCases = {
            {"C1CC2CCC1C2", {{2, 5, {2, 1, 0, 5, 6}, {2, 3, 4, 5, 6}}}}};
    for (const auto &[smiles, expected] : testCases) {
      auto mol = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(mol);
      auto bridgeheads = MolOps::findBridgeheads(*mol);
      CHECK(bridgeheads == expected);
    }
  }
}

TEST_CASE("Bridgehead atoms") {
  SECTION("Simple bridgehead detection") {
    std::vector<std::pair<std::string, std::vector<unsigned int>>> testCases = {
        {"C1CC2CCC1C2", {2, 5}}};
    for (const auto &[smiles, expected] : testCases) {
      auto mol = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(mol);
      auto bridgeheads = MolOps::findBridgeheadAtoms(*mol);
      CHECK(bridgeheads == expected);
    }
  }
}