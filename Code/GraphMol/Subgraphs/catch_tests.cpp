//
//  Copyright (C) 2023-2026 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>
#include <RDGeneral/test.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>

#include <algorithm>

using namespace RDKit;

TEST_CASE("shortestPathsOnly") {
  SECTION("basics, findAllPathsOfLengthN") {
    auto m = "CC1CCC1"_smiles;
    REQUIRE(m);
    bool useBonds = true;
    bool useHs = false;
    int rootedAt = -1;
    bool onlyShortestPaths = true;
    auto ps = findAllPathsOfLengthN(*m, 4);
    CHECK(ps.size() == 3);

    ps = findAllPathsOfLengthN(*m, 4, useBonds, useHs, rootedAt,
                               onlyShortestPaths);
    CHECK(ps.empty());

    ps = findAllPathsOfLengthN(*m, 3);
    CHECK(ps.size() == 6);
    std::cerr << "---------" << std::endl;
    ps = findAllPathsOfLengthN(*m, 3, useBonds, useHs, rootedAt,
                               onlyShortestPaths);
    CHECK(ps.size() == 2);
  }
  SECTION("basics, findAllPathsOfLengthsMtoN") {
    auto m = "CC1CCC1"_smiles;
    REQUIRE(m);
    bool useBonds = true;
    bool useHs = false;
    int rootedAt = -1;
    bool onlyShortestPaths = true;
    auto ps = findAllPathsOfLengthsMtoN(*m, 3, 4);
    CHECK(ps.size() == 2);
    CHECK(ps[3].size() == 6);
    CHECK(ps[4].size() == 3);

    ps = findAllPathsOfLengthsMtoN(*m, 3, 4, useBonds, useHs, rootedAt,
                                   onlyShortestPaths);
    CHECK(ps.size() == 1);
    CHECK(ps[3].size() == 2);
  }
}

TEST_CASE("ignoreAtoms") {
  auto m = "CCCO"_smiles;
  REQUIRE(m);
  SECTION("findAllPathsOfLengthN") {
    bool useBonds = true;
    bool useHs = false;
    int rootedAt = -1;
    bool onlyShortestPaths = true;

    auto ps = findAllPathsOfLengthN(*m, 3, useBonds, useHs, rootedAt,
                                    onlyShortestPaths);
    CHECK(ps.size() == 1);

    // check that ignoreAtoms works:
    boost::dynamic_bitset<> ignoreAtoms(m->getNumAtoms());
    ignoreAtoms[0] = 1;
    auto ps2 = findAllPathsOfLengthN(*m, 3, useBonds, useHs, rootedAt,
                                     onlyShortestPaths, &ignoreAtoms);
    CHECK(ps2.empty());

    // make sure we ignore it even if we are rooted there.
    rootedAt = 1;
    auto ps3 = findAllPathsOfLengthN(*m, 3, useBonds, useHs, rootedAt,
                                     onlyShortestPaths, &ignoreAtoms);
    CHECK(ps3.empty());
  }
  SECTION("findAllPathsOfLengthsMtoN") {
    bool useBonds = true;
    bool useHs = false;
    int rootedAt = -1;
    bool onlyShortestPaths = true;

    auto ps = findAllPathsOfLengthsMtoN(*m, 3, 3, useBonds, useHs, rootedAt,
                                        onlyShortestPaths);
    CHECK(ps.size() == 1);

    // check that ignoreAtoms works:
    boost::dynamic_bitset<> ignoreAtoms(m->getNumAtoms());
    ignoreAtoms[0] = 1;
    auto ps2 = findAllPathsOfLengthsMtoN(*m, 3, 3, useBonds, useHs, rootedAt,
                                         onlyShortestPaths, &ignoreAtoms);
    CHECK(ps2.empty());

    // make sure we ignore it even if we are rooted there.
    rootedAt = 1;
    auto ps3 = findAllPathsOfLengthsMtoN(*m, 3, 3, useBonds, useHs, rootedAt,
                                         onlyShortestPaths, &ignoreAtoms);
    CHECK(ps3.empty());
  }
  SECTION("findAllSubgraphsOfLengthN") {
    bool useHs = false;
    int rootedAt = -1;

    auto ps = findAllSubgraphsOfLengthN(*m, 3, useHs, rootedAt);
    CHECK(ps.size() == 1);

    // check that ignoreAtoms works:
    boost::dynamic_bitset<> ignoreAtoms(m->getNumAtoms());
    ignoreAtoms[0] = 1;
    auto ps2 = findAllSubgraphsOfLengthN(*m, 3, useHs, rootedAt, &ignoreAtoms);
    CHECK(ps2.empty());

    // make sure we ignore it even if we are rooted there.
    rootedAt = 1;
    auto ps3 = findAllSubgraphsOfLengthN(*m, 3, useHs, rootedAt, &ignoreAtoms);
    CHECK(ps3.empty());
  }
  SECTION("findAllSubgraphsOfLengthsMtoN") {
    bool useHs = false;
    int rootedAt = -1;

    auto ps = findAllSubgraphsOfLengthsMtoN(*m, 3, 3, useHs, rootedAt);
    CHECK(ps.size() == 1);
    CHECK(ps[3].size() == 1);

    // check that ignoreAtoms works:
    boost::dynamic_bitset<> ignoreAtoms(m->getNumAtoms());
    ignoreAtoms[0] = 1;
    auto ps2 =
        findAllSubgraphsOfLengthsMtoN(*m, 3, 3, useHs, rootedAt, &ignoreAtoms);
    CHECK(ps2.size() == 1);
    CHECK(ps2[3].empty());

    // make sure we ignore it even if we are rooted there.
    rootedAt = 1;
    auto ps3 =
        findAllSubgraphsOfLengthsMtoN(*m, 3, 3, useHs, rootedAt, &ignoreAtoms);
    CHECK(ps3.size() == 1);
    CHECK(ps3[3].empty());
  }
}