//
//  Copyright (C) 2020-2021 Brian P Kelley, Joann Prescott-Roy and other RDKit
//  contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Deprotect/Deprotect.h>
#include <boost/algorithm/string.hpp>

using namespace RDKit;
using namespace RDKit::Deprotect;
TEST_CASE("Deprotection basics", "[deprotect]") {
  const auto dps = getDeprotections();
  DeprotectData cp = dps[0];
  CHECK(cp == dps[0]);
  CHECK(dps[0] == cp);
  CHECK(dps[1] != dps[0]);
}

TEST_CASE("Standard deprotections", "[deprotect]") {
  SECTION("simple deprotections") {
    auto m = "N(C(=O)OC(C)(C)C)Cc1ccccc1NC(=O)OC(C)(C)C"_smiles;
    auto res = deprotect(*m);
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "NCc1ccccc1N");
    CHECK(res->getProp<int>("DEPROTECTION_COUNT") == 2);
    std::vector<std::string> expected{"Boc", "Boc"};
    CHECK(res->getProp<std::vector<std::string>>("DEPROTECTIONS") == expected);
  }
  SECTION("test deprotection examples") {
    for (auto &data : getDeprotections()) {
      std::vector<DeprotectData> vect = {data};
      std::vector<std::string> examples;
      boost::split(examples, data.example, boost::is_any_of(">"));
      std::unique_ptr<ROMol> start(SmilesToMol(examples[0]));
      std::unique_ptr<ROMol> end(SmilesToMol(examples[2]));
      // check the one
      auto res = deprotect(*start, vect);
      CHECK(MolToSmiles(*res) == MolToSmiles(*end));
      // check them all

      res = deprotect(*start);
      CHECK(MolToSmiles(*res) == MolToSmiles(*end));
    }
  }
}

TEST_CASE("Standard deprotections in place", "[deprotect]") {
  SECTION("simple deprotections") {
    auto m = "N(C(=O)OC(C)(C)C)Cc1ccccc1NC(=O)OC(C)(C)C"_smiles;
    CHECK(deprotectInPlace(*m));
    CHECK(MolToSmiles(*m) == "NCc1ccccc1N");
    CHECK(m->getProp<int>("DEPROTECTION_COUNT") == 2);
    std::vector<std::string> expected{"Boc", "Boc"};
    CHECK(m->getProp<std::vector<std::string>>("DEPROTECTIONS") == expected);
  }
  SECTION("test deprotection examples") {
    for (auto &data : getDeprotections()) {
      std::vector<DeprotectData> vect = {data};
      std::vector<std::string> examples;
      boost::split(examples, data.example, boost::is_any_of(">"));
      std::unique_ptr<RWMol> end(SmilesToMol(examples[2]));
      // check the one
      {
        std::unique_ptr<RWMol> start(SmilesToMol(examples[0]));
        deprotectInPlace(*start, vect);
        CHECK(MolToSmiles(*start) == MolToSmiles(*end));
      }

      // check them all
      {
        std::unique_ptr<RWMol> start(SmilesToMol(examples[0]));
        deprotectInPlace(*start);
        CHECK(MolToSmiles(*start) == MolToSmiles(*end));
      }
    }
  }
}
