//
//  Copyright (c) 2023, Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Tests of the generalized substructure searching code
//

#include "catch.hpp"

#include <tuple>
#include <utility>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/GeneralizedSubstruct/XQMol.h>
#include <GraphMol/TautomerQuery/TautomerQuery.h>
#include <GraphMol/MolEnumerator/MolEnumerator.h>

using namespace RDKit;

TEST_CASE("tautomer basics") {
  auto mol = "Cc1n[nH]c(F)c1"_smiles;
  REQUIRE(mol);
  ExtendedQueryMol xqm =
      std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol));
  SECTION("substructure matching") {
    CHECK(SubstructMatch(*"CCc1n[nH]c(F)c1"_smiles, xqm).size() == 1);
    CHECK(SubstructMatch(*"CCc1[nH]nc(F)c1"_smiles, xqm).size() == 1);
    CHECK(SubstructMatch(*"CCc1[nH]ncc1"_smiles, xqm).empty());
  }
}

TEST_CASE("enumeration basics") {
  auto mol = "COCC |LN:1:1.3|"_smiles;
  REQUIRE(mol);
  ExtendedQueryMol xqm =
      std::make_unique<MolBundle>(MolEnumerator::enumerate(*mol));
  SECTION("substructure matching") {
    CHECK(SubstructMatch(*"COCC"_smiles, xqm).size() == 1);
    CHECK(SubstructMatch(*"COOCC"_smiles, xqm).size() == 1);
    CHECK(SubstructMatch(*"COOOCC"_smiles, xqm).size() == 1);
    CHECK(SubstructMatch(*"COOOOCC"_smiles, xqm).empty());
  }
}
