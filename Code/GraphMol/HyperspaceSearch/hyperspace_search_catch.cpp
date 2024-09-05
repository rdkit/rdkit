//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include "HyperspaceSubstructureSearch.h"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <catch2/catch_all.hpp>

using namespace RDKit;
using namespace RDKit::HyperspaceSSSearch;

namespace RDKit::HyperspaceSSSearch::details {
std::vector<std::vector<std::unique_ptr<ROMol>>> splitMolecule(
    const ROMol &query, unsigned int maxBondSplits);
}  // namespace RDKit::HyperspaceSSSearch::details

using namespace RDKit::HyperspaceSSSearch::details;

std::string fName = getenv("RDBASE");
const std::string LIB_NAME =
    fName + "/Code/GraphMol/HyperspaceSearch/data/idorsia_toy_space_a.txt";

TEST_CASE("Test splits") {
  std::vector<std::string> smiles{"c1ccccc1CN1CCN(CC1)C(-O)c1ncc(F)cc1",
                                  "CC(C)OCc1nnc(N2CC(C)CC2)n1C1CCCC1"};
  std::vector<std::vector<size_t>> expCounts{{6, 60, 350}, {8, 58, 326}};
  for (size_t i = 0; i < smiles.size(); ++i) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles[i]);
    REQUIRE(mol);
    auto fragments = splitMolecule(*mol, 3);
    REQUIRE(fragments.size() == 3);
    CHECK(fragments[0].size() == expCounts[i][0]);
    CHECK(fragments[1].size() == expCounts[i][1]);
    CHECK(fragments[2].size() == expCounts[i][2]);
  }
}
TEST_CASE("Simple query") {
  auto query = "c1ccccc1C(=O)N1CCCC1"_smiles;
  auto results = SSSearch(*query, 1, LIB_NAME);
}
