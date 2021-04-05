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

TEST_CASE("basic queries", "[query]") {
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
    auto mol = "[C;R2][$(CO),r6]"_smarts;
    REQUIRE(mol);
    auto json = MolInterchange::MolToJSONData(*mol);
    CHECK(json.find("\"subquery\":") != std::string::npos);
    auto nmols = MolInterchange::JSONDataToMols(json);
    CHECK(nmols.size() == 1);
    CHECK(MolToSmarts(*nmols[0]) == MolToSmarts(*mol));
  }
}