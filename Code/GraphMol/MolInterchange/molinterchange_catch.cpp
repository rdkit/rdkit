//
//  Copyright (C) 2021 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include "MolInterchange.h"

using namespace RDKit;

TEST_CASE("queries to JSON", "[query]") {
  SECTION("SMARTS 1") {
    auto mol = "C[cH]:[#7+]"_smarts;
    REQUIRE(mol);
    auto json = MolInterchange::MolToJSONData(*mol);
    // no implicit Hs added with queries
    CHECK(json.find("{\"impHs\":3}") == std::string::npos);
    CHECK(json.find("\"descr\":\"AtomHCount\"") != std::string::npos);
    auto nmols = MolInterchange::JSONDataToMols(json);
    CHECK(nmols.size() == 1);
    CHECK(MolToSmarts(*nmols[0]) == MolToSmarts(*mol));
  }
  SECTION("Recursive SMARTS 1") {
    auto mol = "[C;R2]=,#!@[$(CO),r6]"_smarts;
    REQUIRE(mol);
    auto json = MolInterchange::MolToJSONData(*mol);
    CHECK(json.find("\"subquery\":") != std::string::npos);
    auto nmols = MolInterchange::JSONDataToMols(json);
    CHECK(nmols.size() == 1);
    CHECK(MolToSmarts(*nmols[0]) == MolToSmarts(*mol));
  }
  SECTION("mol blocks") {
    auto mol = R"CTAB(
  Mrv2102 04092105442D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Q -1.4583 0.4583 0 0
M  V30 2 C -0.1247 1.2283 0 0
M  V30 3 C 1.209 0.4583 0 0
M  V30 4 C 2.5427 1.2283 0 0
M  V30 5 C 1.209 -1.0817 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 4
M  V30 3 5 2 3
M  V30 4 1 3 5 TOPO=2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    auto json = MolInterchange::MolToJSONData(*mol);
    CHECK(json.find("\"type\":\"Q\"") != std::string::npos);
    auto ps = MolInterchange::JSONParseParameters();
    ps.useHCounts = false;
    auto nmols = MolInterchange::JSONDataToMols(json, ps);
    CHECK(nmols.size() == 1);
    auto mb = MolToV3KMolBlock(*nmols[0]);
    CHECK(mb.find(" Q ") != std::string::npos);
    CHECK(mb.find("TOPO=2") != std::string::npos);
  }
}