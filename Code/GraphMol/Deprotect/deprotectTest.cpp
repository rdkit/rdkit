//
//  Copyright (C) 2020 Brian P Kelley, Joann Prescott-Roy
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Deprotect/Deprotect.h>
#include <boost/algorithm/string.hpp>

using namespace RDKit;
using namespace RDKit::Deprotect;
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
    for(auto &data : getDeprotections()) {
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
