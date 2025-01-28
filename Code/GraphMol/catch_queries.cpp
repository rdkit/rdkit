//
//  Copyright (C) 2023 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

TEST_CASE(
    "github #6106: Dummy atoms should not be considered to be metals for M and MH queries") {
  const auto m = "C*[Fe]"_smiles;
  REQUIRE(m);

  SECTION("M") {
    std::unique_ptr<ATOM_OR_QUERY> q(makeMAtomQuery());
    REQUIRE(q);
    CHECK(!q->Match(m->getAtomWithIdx(0)));
    CHECK(!q->Match(m->getAtomWithIdx(1)));
    CHECK(q->Match(m->getAtomWithIdx(2)));
  }
  SECTION("MH") {
    std::unique_ptr<ATOM_OR_QUERY> q(makeMHAtomQuery());
    REQUIRE(q);
    CHECK(!q->Match(m->getAtomWithIdx(0)));
    CHECK(!q->Match(m->getAtomWithIdx(1)));
    CHECK(q->Match(m->getAtomWithIdx(2)));
  }
}
