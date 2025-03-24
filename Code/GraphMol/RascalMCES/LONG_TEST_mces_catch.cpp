//
//  Copyright (C) 2025 David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <chrono>
#include <random>
#include <vector>

#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <catch2/catch_all.hpp>

#include <GraphMol/RascalMCES/RascalMCES.h>
#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/RascalMCES/RascalResult.h>
#include <GraphMol/RascalMCES/RascalDetails.h>

using namespace RDKit;
using namespace RDKit::RascalMCES;
using namespace RDKit::RascalMCES::details;

void check_smarts_ok(const RDKit::ROMol &mol1, const RDKit::ROMol &mol2,
                     const RascalResult &res) {
  auto qmol = v2::SmilesParse::MolFromSmarts(res.getSmarts());
  REQUIRE(qmol);
  RDKit::MatchVectType dont_care;
  CHECK(RDKit::SubstructMatch(mol1, *qmol, dont_care));
  CHECK(RDKit::SubstructMatch(mol2, *qmol, dont_care));
}

TEST_CASE("benchmarks") {
  // As well as timings, this also tests that the same result is produced from
  // random ordering of the input molecules i.e. the results aren't input order
  // dependent.
  std::vector<std::tuple<std::string, std::string, std::string, double, int,
                         unsigned int, int>>
      tests = {{"juglone", "Oc1cccc2C(=O)C=CC(=O)c12",
                "O1C(=O)C=Cc2cc(OC)c(O)cc12", 0.5, 100, 8, 10},
               {"methadone", "c1ccccc1C(C(=O)CC)(c1ccccc1)CC(C)N(C)C",
                "c1ccccc1C1(CCN(C)CC1)C(=O)OCC", 0.6, 100, 16, 20},
               {"symmetrical",
                "c1c(OC)c(OC)c(OC)cc1C(=O)OCCCOC(=O)c1cc(OC)c(OC)c(OC)c1",
                "c1c(OC)c(OC)c(OC)cc1C(=O)OCCOC(=O)c1cc(OC)c(OC)c(OC)c1", 0.7,
                100, 32, 50},
               {"testosterone", "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",
                "CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O", 0.6, 100, 16, 50}};
  std::vector<double> timings;
  std::random_device rd;
  std::mt19937 g(rd());
  // It's convenient if the same test is run each time.
  g.seed(1);
  RascalOptions opts;
  for (const auto &test : tests) {
    auto m1 = std::unique_ptr<ROMol>(RDKit::SmilesToMol(std::get<1>(test)));
    auto m2 = std::unique_ptr<ROMol>(RDKit::SmilesToMol(std::get<2>(test)));
    opts.similarityThreshold = std::get<3>(test);
    std::vector<unsigned int> atom1Inds(m1->getNumAtoms(), 0);
    std::iota(atom1Inds.begin(), atom1Inds.end(), 0);
    std::vector<unsigned int> atom2Inds(m2->getNumAtoms(), 0);
    std::iota(atom2Inds.begin(), atom2Inds.end(), 0);
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < std::get<4>(test); ++i) {
      std::shuffle(atom1Inds.begin(), atom1Inds.end(), g);
      std::unique_ptr<ROMol> randmol1(MolOps::renumberAtoms(*m1, atom1Inds));
      std::shuffle(atom2Inds.begin(), atom2Inds.end(), g);
      std::unique_ptr<ROMol> randmol2(MolOps::renumberAtoms(*m2, atom2Inds));

      auto res = rascalMCES(*randmol1, *randmol2, opts);
      REQUIRE(res.size() == 1);
      REQUIRE(std::get<5>(test) == res.front().getBondMatches().size());
      check_smarts_ok(*randmol1, *randmol2, res.front());
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    timings.push_back(
        std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() /
        std::get<4>(test));
  }

  // Make sure we haven't slowed it down a lot.
  //  for (size_t i = 0; i < tests.size(); ++i) {
  //    std::cout << "TIMING : " << std::get<0>(tests[i]) << " took " <<
  //    timings[i]
  //              << " milliseconds on average of " << std::get<4>(tests[i])
  //              << " runs\n";
  //  }
  for (size_t i = 0; i < tests.size(); ++i) {
    auto ref_time = std::get<6>(tests[i]);
#ifndef NDEBUG  // allow more time for debug builds
    ref_time *= 5;
#endif
    REQUIRE(timings[i] < ref_time);
  }
}

TEST_CASE("timeout") {
  auto m1 =
      "O[C@@H]1CC[C@H](C[C@H]1OC)C[C@@H](C)[C@@H]4CC(=O)[C@H](C)/C=C(\\C)[C@@H](O)[C@@H](OC)C(=O)[C@H](C)C[C@H](C)\\C=C\\C=C\\C=C(/C)[C@@H](OC)C[C@@H]2CC[C@@H](C)[C@@](O)(O2)C(=O)C(=O)N3CCCC[C@H]3C(=O)O4"_smiles;
  REQUIRE(m1);
  auto m2 =
      "CCC1C=C(CC(CC(C2C(CC(C(O2)(C(=O)C(=O)N3CCCCC3C(=O)OC(C(C(CC1=O)O)C)C(=CC4CCC(C(C4)OC)O)C)O)C)OC)OC)C)C"_smiles;
  REQUIRE(m2);

  {
    RascalOptions opts;
    opts.timeout = 10;
    opts.maxBondMatchPairs = 2000;
    auto res = rascalMCES(*m1, *m2, opts);
    REQUIRE(res.size() == 1);
    REQUIRE(res.front().getBondMatches().size() >= 39);
    check_smarts_ok(*m1, *m2, res.front());
  }
  {
    RascalOptions opts;
    opts.timeout = 70;
    opts.maxBondMatchPairs = 2000;
    auto res = rascalMCES(*m1, *m2, opts);
    REQUIRE(res.size() == 1);
    REQUIRE(res.front().getBondMatches().size() >= 44);
    check_smarts_ok(*m1, *m2, res.front());
  }
  {
    RascalOptions opts;
    opts.timeout = 120;
    opts.maxBondMatchPairs = 2000;
    auto res = rascalMCES(*m1, *m2, opts);
    REQUIRE(res.size() == 1);
    REQUIRE(res.front().getBondMatches().size() >= 44);
    check_smarts_ok(*m1, *m2, res.front());
  }
}

TEST_CASE("Duplicate Single Largest Frag") {
  {
    auto m1 = "c1ccc2c(c1)c(ncn2)CNCc3ccc(cc3)Cl"_smiles;
    REQUIRE(m1);
    auto m2 = "c1ccc2c(c1)c(ncn2)CCCc3ccc(cc3)Cl"_smiles;
    REQUIRE(m2);
    RascalOptions opts;
    opts.singleLargestFrag = true;
    opts.allBestMCESs = true;
    opts.timeout = -1;
    auto res = rascalMCES(*m1, *m2, opts);
    REQUIRE(res.size() == 1);
    CHECK(res.front().getAtomMatches().size() == 11);
    CHECK(res.front().getBondMatches().size() == 12);
  }
  {
    // Without the fix, there were 156 results sets, with different breaks
    // in the long chain giving different numbers of atoms in the largest
    // fragment.
    auto m1 = "CN(C)c1ccc(CC(=O)NCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1"_smiles;
    REQUIRE(m1);
    auto m2 = "CN(C)c1ccc(CC(=O)NCCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1"_smiles;
    REQUIRE(m2);
    RascalOptions opts;
    opts.singleLargestFrag = true;
    opts.allBestMCESs = true;
    opts.timeout = -1;
    auto res = rascalMCES(*m1, *m2, opts);
    REQUIRE(!res.empty());
    CHECK(res.size() == 4);
    CHECK(res.front().getAtomMatches().size() == 23);
    CHECK(res.front().getBondMatches().size() == 23);
    CHECK(res.back().getAtomMatches().size() == 23);
    CHECK(res.back().getBondMatches().size() == 23);
  }
}

TEST_CASE("Github8255 - incorrect MCES with singleLargestFrag=true") {
  RascalOptions opts;
  opts.singleLargestFrag = true;
  opts.allBestMCESs = true;
  opts.ignoreAtomAromaticity = true;
  opts.ringMatchesRingOnly = true;
  opts.completeAromaticRings = false;
  {
    auto m1 = "Cc1ccc(cn1)-c1ccc(cn1)-c1ccc2ccccc2n1"_smiles;
    REQUIRE(m1);
    auto m2 = "Cc1ccc(cn1)-c1ccc(cn1)-c1cc2ccccc2[nH]1"_smiles;
    REQUIRE(m2);
    auto res = rascalMCES(*m1, *m2, opts);
    CHECK(res.size() == 3);
    for (auto &r : res) {
      CHECK(r.getAtomMatches().size() == 22);
      CHECK(r.getBondMatches().size() == 24);
    }
  }
  {
    auto m1 = "Cc1ncc(cc1)c2ncc(cc2)c3ccc4c(c3)cn[nH]4"_smiles;
    REQUIRE(m1);
    auto m2 = "Cc1ncc(cc1)c2ncc(cc2)c3ccc4c(c3)nccn4"_smiles;
    REQUIRE(m2);
    auto res = rascalMCES(*m1, *m2, opts);
    CHECK(res.size() == 2);
    for (auto &r : res) {
      CHECK(r.getAtomMatches().size() == 20);
      CHECK(r.getBondMatches().size() == 22);
    }
  }
}

TEST_CASE("Github8360 - another incorrect MCES with singleLargestFrag=true") {
  RascalOptions opts;
  opts.singleLargestFrag = true;
  opts.allBestMCESs = true;
  opts.completeAromaticRings = false;
  {
    auto m2 = "c1ccc(cn1)-c1ccc2ccccc2n1"_smiles;
    REQUIRE(m2);
    auto m1 = "c1ccc(cn1)-c1cc2ccccc2[nH]1"_smiles;
    REQUIRE(m1);
    auto res = rascalMCES(*m1, *m2, opts);
    CHECK(res.size() == 3);
    for (auto &r : res) {
      CHECK(r.getAtomMatches().size() == 15);
      CHECK(r.getBondMatches().size() == 16);
    }
  }
  opts.singleLargestFrag = true;
  opts.allBestMCESs = true;
  opts.maxBondMatchPairs = 1200;
  {
    auto m1 = "C=CCCC=CCCC1CNCCC1"_smiles;
    REQUIRE(m1);
    auto m2 = "C=CCC=CCCCC1CNCCC1"_smiles;
    REQUIRE(m2);
    auto res = rascalMCES(*m1, *m2, opts);
    REQUIRE(!res.empty());
    CHECK(res.size() == 1);
    for (auto &r : res) {
      std::cout << r.getSmarts() << std::endl;
      CHECK(r.getAtomMatches().size() == 10);
      CHECK(r.getBondMatches().size() == 9);
    }
  }
}
