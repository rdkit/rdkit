//
//  Copyright (C) 2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/test_fixtures.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;

TEST_CASE("Github #8424: direction on aromatic bonds in SMARTS") {
  SECTION("simplified") {
    auto m = "C/N=c1/[nH]cc(Br)nc1"_smiles;
    REQUIRE(m);
    auto smarts = MolToSmarts(*m);
    CHECK(smarts == "[#6]/[#7]=[#6]1/[#7]:[#6]:[#6](-[#35]):[#7]:[#6]:1");
  }
  SECTION("as reported") {
    auto m =
        "CN1C(=O)CN(c2cc(C3CC3)cn3cc(CC(=O)N/N=c4\\cnc(Br)c[nH]4)nc23)C1=O"_smiles;
    REQUIRE(m);
    auto smarts = MolToSmarts(*m);
    // should have slashes in both directions
    CHECK(smarts.find("/") != std::string::npos);
    CHECK(smarts.find("\\") != std::string::npos);
  }
}

TEST_CASE("repeated explicit H counts and charges") {
  SECTION("h counts") {
    std::vector<std::string> smartses = {
        "[N&H3&H0]",
        "[N&H0&H3]",
    };
    for (const auto &smarts : smartses) {
      INFO(smarts);
      auto m = v2::SmilesParse::MolFromSmarts(smarts);
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(0)->getNoImplicit());
      CHECK(m->getAtomWithIdx(0)->getNumExplicitHs() == 0);
    }
  }
  SECTION("charges") {
    std::vector<std::string> smartses = {"[N&+&+0]", "[N&+0&+2]"};
    for (const auto &smarts : smartses) {
      INFO(smarts);
      auto m = v2::SmilesParse::MolFromSmarts(smarts);
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(0)->getFormalCharge() == 0);
    }
  }
}

TEST_CASE("implicit Hs from SMILES should not make it into SMARTS") {
  SECTION("aromatic N") {
    auto m = "c1ccc[nH]1"_smiles;
    REQUIRE(m);
    auto smarts = MolToSmarts(*m);
    CHECK(smarts == "[#6]1:[#6]:[#6]:[#6]:[#7]:1");
  }
  SECTION("chirality") {
    // as of this writing, we still keep Hs in the SMARTS for chiral centers.
    //   it's inconsistent to do so, but we have explicitly punted on fixing
    //   this for now.
    auto m = "C[C@H](N)F"_smiles;
    REQUIRE(m);
    auto smarts = MolToSmarts(*m);
    CHECK(smarts == "[#6]-[#6@H](-[#7])-[#9]");
  }
}