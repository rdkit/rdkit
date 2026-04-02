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

void checkMatches(const std::string &smarts, const std::string &smiles,
                  unsigned int nMatches, unsigned int lenFirst,
                  bool addHs = false) {
  // utility function that will find the matches between a smarts and smiles
  // if they match the expected values
  //  smarts : smarts string
  //  smiles : smiles string
  //  nMatches : expected number of matches
  //  lenFirst : length of the first match
  //
  // Return the list of all matches just in case want to do additional testing
  INFO(smarts + " " + smiles);
  auto matcher = v2::SmilesParse::MolFromSmarts(smarts);
  REQUIRE(matcher);
  // we will at the same time test the serialization:
  std::string pickle;
  MolPickler::pickleMol(*matcher, pickle);
  ROMol matcher2(pickle);

  auto mol = v2::SmilesParse::MolFromSmiles(smiles);
  REQUIRE(mol);
  if (addHs) {
    MolOps::addHs(*mol);
  }
  MolOps::findSSSR(*mol);

  MatchVectType mV;
  auto matches = SubstructMatch(*mol, *matcher, mV);
  CHECK(matches);
  CHECK(mV.size() == lenFirst);
  std::vector<MatchVectType> mVV;
  auto uniquify = true;
  auto matchCount = SubstructMatch(*mol, *matcher, mVV, uniquify);
  CHECK(matchCount == nMatches);
  CHECK(mVV[0].size() == lenFirst);

  matches = SubstructMatch(*mol, matcher2, mV);
  CHECK(matches);
  CHECK(mV.size() == lenFirst);
  matchCount = SubstructMatch(*mol, matcher2, mVV, true);
  CHECK(matchCount == nMatches);
  CHECK(mVV[0].size() == lenFirst);
}

TEST_CASE("k SMARTS extensions") {
  SECTION("parsing and writing") {
    auto q = "[k4]"_smarts;
    REQUIRE(q);
    auto smarts = MolToSmarts(*q);
    CHECK(smarts == "[k4]");
  }
  SECTION("matching") {
    auto m = "C1CC2N1CCCCCC2"_smiles;
    REQUIRE(m);
    auto k4 = "[k4]"_smarts;
    REQUIRE(k4);
    std::vector<MatchVectType> matches;
    CHECK(SubstructMatch(*m, *k4, matches));
    CHECK(matches.size() == 4);
    auto q4 = "[r4]"_smarts;
    REQUIRE(q4);
    CHECK(SubstructMatch(*m, *q4, matches));
    CHECK(matches.size() == 4);
    auto k8 = "[k8]"_smarts;
    REQUIRE(k8);
    CHECK(SubstructMatch(*m, *k8, matches));
    CHECK(matches.size() == 8);
    auto q8 = "[r8]"_smarts;
    REQUIRE(q8);
    CHECK(SubstructMatch(*m, *q8, matches));
    CHECK(matches.size() == 6);
  }
  SECTION("ranges") {
    std::string smiles = "C1CC2N1CCCCCC2C";
    std::vector<std::pair<std::string, size_t>> smartses = {
        {"[r{4-5}]", 4},  {"[k{4-5}]", 4},  {"[k{4-}]", 10}, {"[r{4-}]", 10},
        {"[k{-5}]", 4},   {"[k{8-}]", 8},   {"[r{8-}]", 6},  {"[r{4-8}]", 10},
        {"[!r{4-5}]", 7}, {"[!k{4-5}]", 7}, {"[!k{4-}]", 1}, {"[!r{4-}]", 1},
        {"[!k{-5}]", 7},  {"[!k{8-}]", 3},  {"[!r{8-}]", 5}, {"[!r{4-8}]", 1},
        {"[k]", 10},      {"[r]", 10},      {"[!k]", 1},     {"[!r]", 1},
    };
    for (const auto &[sma, val] : smartses) {
      checkMatches(sma, smiles, val, 1);
    }
  }
}
