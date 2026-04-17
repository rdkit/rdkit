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

TEST_CASE("range queries for atom ring membership") {
  auto m = "C1C2C1CC2C"_smiles;
  REQUIRE(m);
  SECTION("as range queries") {
    {
      std::unique_ptr<ATOM_RANGE_QUERY> q(makeAtomInRingOfSizeQuery(3, 4));
      REQUIRE(q);
      CHECK(q->Match(m->getAtomWithIdx(0)));
      CHECK(q->Match(m->getAtomWithIdx(1)));
      CHECK(q->Match(m->getAtomWithIdx(3)));
      CHECK(!q->Match(m->getAtomWithIdx(5)));
    }
    {
      std::unique_ptr<ATOM_RANGE_QUERY> q(makeAtomInRingOfSizeQuery(4, 6));
      REQUIRE(q);
      CHECK(!q->Match(m->getAtomWithIdx(0)));
      CHECK(q->Match(m->getAtomWithIdx(1)));
      CHECK(q->Match(m->getAtomWithIdx(3)));
      CHECK(!q->Match(m->getAtomWithIdx(5)));
    }
    {
      // this is how we use the queries in the SMARTS parser
      std::unique_ptr<ATOM_LESSEQUAL_QUERY> q(
          makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>(3, [](Atom const *at) {
            return queryAtomIsInRingOfSize(at, 3, -1);
          }));
      REQUIRE(q);
      CHECK(q->Match(m->getAtomWithIdx(0)));
      CHECK(q->Match(m->getAtomWithIdx(1)));
      CHECK(q->Match(m->getAtomWithIdx(3)));
      CHECK(!q->Match(m->getAtomWithIdx(5)));
    }
    {  // this is how we use the queries in the SMARTS parser
      std::unique_ptr<ATOM_GREATEREQUAL_QUERY> q(
          makeAtomSimpleQuery<ATOM_GREATEREQUAL_QUERY>(3, [](Atom const *at) {
            return queryAtomIsInRingOfSize(at, -1, 3);
          }));
      REQUIRE(q);
      CHECK(q->Match(m->getAtomWithIdx(0)));
      CHECK(q->Match(m->getAtomWithIdx(1)));
      CHECK(!q->Match(m->getAtomWithIdx(3)));
      CHECK(!q->Match(m->getAtomWithIdx(5)));
    }
  }
  SECTION("query function") {
    CHECK(queryAtomIsInRingOfSize(m->getAtomWithIdx(0), 3, -1) == 3);
    CHECK(queryAtomIsInRingOfSize(m->getAtomWithIdx(1), 3, -1) == 3);
    CHECK(queryAtomIsInRingOfSize(m->getAtomWithIdx(3), 3, -1) == 4);
    CHECK(queryAtomIsInRingOfSize(m->getAtomWithIdx(5), 3, -1) == -1);

    CHECK(queryAtomIsInRingOfSize(m->getAtomWithIdx(0), -1, 3) == 3);
    CHECK(queryAtomIsInRingOfSize(m->getAtomWithIdx(1), -1, 3) == 3);
    CHECK(queryAtomIsInRingOfSize(m->getAtomWithIdx(3), -1, 3) ==
          std::numeric_limits<int>::max());
    CHECK(queryAtomIsInRingOfSize(m->getAtomWithIdx(5), -1, 3) ==
          std::numeric_limits<int>::max());

    CHECK(queryAtomIsInRingOfSize(m->getAtomWithIdx(0), 3, 4) == 3);
    CHECK(queryAtomIsInRingOfSize(m->getAtomWithIdx(1), 3, 4) == 3);
    CHECK(queryAtomIsInRingOfSize(m->getAtomWithIdx(3), 3, 4) == 4);
    CHECK(queryAtomIsInRingOfSize(m->getAtomWithIdx(5), 3, 4) == -1);

    CHECK(queryAtomIsInRingOfSize(m->getAtomWithIdx(3), 0, 4) == 4);
    CHECK(queryAtomIsInRingOfSize(m->getAtomWithIdx(5), 0, 4) == -1);
  }
}