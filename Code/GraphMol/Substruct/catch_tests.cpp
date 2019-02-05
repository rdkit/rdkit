//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <tuple>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;
typedef std::tuple<std::string, std::string, size_t> matchCase;

TEST_CASE("substructure parameters", "[substruct]") {
  SECTION("chirality") {
    auto mol1 = "CCC[C@@H]1CN(CCC)CCN1"_smiles;
    auto mol2 = "CCC[C@H]1CN(CCC)CCN1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);

    SubstructMatchParameters ps;
    // default is to ignore chirality:
    CHECK(SubstructMatch(*mol1, *mol2, ps).size() == 1);
    CHECK(SubstructMatch(*mol1, *mol1, ps).size() == 1);

    ps.useChirality = true;
    CHECK(SubstructMatch(*mol1, *mol2, ps).size() == 0);
    CHECK(SubstructMatch(*mol1, *mol1, ps).size() == 1);
  }
  SECTION("conjugated matching aromaticity 1") {
    auto mol1 = "C1=COC=C1"_smiles;
    REQUIRE(mol1);
    RWMol mol2(*mol1);
    MolOps::Kekulize(mol2);
    SubstructMatchParameters ps;
    CHECK(SubstructMatch(*mol1, mol2, ps).size() == 0);
    CHECK(SubstructMatch(mol2, *mol1, ps).size() == 0);

    ps.aromaticMatchesConjugated = true;
    CHECK(SubstructMatch(*mol1, mol2, ps).size() == 1);
    CHECK(SubstructMatch(mol2, *mol1, ps).size() == 1);
  }
  SECTION("conjugated matching aromaticity 2") {
    auto mol1 = "c1ccccc1"_smiles;
    REQUIRE(mol1);
    RWMol mol2(*mol1);
    MolOps::Kekulize(mol2);
    SubstructMatchParameters ps;
    CHECK(SubstructMatch(*mol1, mol2, ps).size() == 0);
    CHECK(SubstructMatch(mol2, *mol1, ps).size() == 0);

    ps.aromaticMatchesConjugated = true;
    CHECK(SubstructMatch(*mol1, mol2, ps).size() == 1);
    CHECK(SubstructMatch(mol2, *mol1, ps).size() == 1);
  }

  SECTION("conjugated matching aromaticity bulk") {
    std::vector<matchCase> examples;
    examples.push_back(std::make_tuple(std::string("c1ccccc1"), std::string("C1CCCCC1"), 0));
    examples.push_back(std::make_tuple(std::string("C1CCCCC1"), std::string("c1ccccc1"), 0));
    examples.push_back(std::make_tuple(std::string("O=C1C=CC(=O)C=C1"), std::string("c1ccccc1"), 1));
    SubstructMatchParameters ps;
    ps.aromaticMatchesConjugated = true;
    for (const auto &example : examples) {
      // std::cerr << "   " << std::get<0>(example) << " - "
      //           << std::get<1>(example) << std::endl;
      std::unique_ptr<RWMol> m1(SmilesToMol(std::get<0>(example)));
      REQUIRE(m1);
      std::unique_ptr<RWMol> m2(SmilesToMol(std::get<1>(example)));
      CHECK(SubstructMatch(*m1, *m2, ps).size() == std::get<2>(example));
    }
  }
  SECTION("looping") {
    auto mol1 = "CC(=O)C(=O)C(=O)"_smiles;
    auto mol2 = "C=O"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    for( auto match : SubstructMatch(*mol1,*mol2)){
      CHECK(match.size()==2);
    }
  }

}
