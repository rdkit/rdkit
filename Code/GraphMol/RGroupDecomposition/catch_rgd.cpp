//
//  Copyright (C) 2021 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"

#include <GraphMol/RDKitBase.h>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RGroupDecomposition/RGroupDecomp.h>
#include <GraphMol/RGroupDecomposition/RGroupUtils.h>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim_all.hpp>

using namespace RDKit;

template <typename T>
void initDataset(T &suppl, ROMOL_SPTR &core, std::vector<ROMOL_SPTR> &mols) {
  core.reset(suppl[0]);
  REQUIRE(core);
  for (unsigned int i = 1; i < suppl.length(); ++i) {
    mols.emplace_back(suppl[i]);
    REQUIRE(mols.back());
  }
}

std::string flatten_whitespace(const std::string &txt) {
  auto res = txt;
  boost::algorithm::trim_fill_if(res, "", boost::is_any_of(" \t\n"));
  return res;
}

std::string readReferenceData(const std::string &fname) {
  std::ifstream ins(fname);
  std::string res;
  ins.seekg(0, std::ios::end);
  res.reserve(ins.tellg());
  ins.seekg(0, std::ios::beg);
  res.assign((std::istreambuf_iterator<char>(ins)),
             std::istreambuf_iterator<char>());
  return res;
}
TEST_CASE("toJSONTests", "[unittests]") {
  std::string testDataDir =
      std::string(getenv("RDBASE")) +
      std::string("/Code/GraphMol/RGroupDecomposition/test_data/");
  std::string fName = testDataDir + "simple1.sdf";
  SDMolSupplier suppl(fName);
  std::vector<ROMOL_SPTR> cores(1);
  std::vector<ROMOL_SPTR> mols;
  initDataset(suppl, cores.front(), mols);
  SECTION("rows") {
    RGroupRows rows;
    auto n = RGroupDecompose(cores, mols, rows);
    CHECK(n == mols.size());
    CHECK(rows.size() == mols.size());
    std::string expected = R"JSON([
    {
        "Core": "Cc1cccc([*:2])c1[*:1]",
        "R1": "CO[*:1]",
        "R2": "[H][*:2]"
    },
    {
        "Core": "Cc1cccc([*:2])c1[*:1]",
        "R1": "CO[*:1]",
        "R2": "[H][*:2]"
    },
    {
        "Core": "Cc1cccc([*:2])c1[*:1]",
        "R1": "[H][*:1]",
        "R2": "CO[*:2]"
    }
])JSON";
    CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(expected));
  }
  SECTION("columns") {
    RGroupColumns cols;
    auto n = RGroupDecompose(cores, mols, cols);
    CHECK(n == mols.size());
    CHECK(cols.size() == mols.size());
    std::string expected = R"JSON([
  "Core": [
    "Cc1cccc([*:2])c1[*:1]",
    "Cc1cccc([*:2])c1[*:1]",
    "Cc1cccc([*:2])c1[*:1]"
  ],
  "R1": [
    "CO[*:1]",
    "CO[*:1]",
    "[H][*:1]"
  ],
  "R2": [
    "[H][*:2]",
    "[H][*:2]",
    "CO[*:2]"
  ]
]
)JSON";
    CHECK(flatten_whitespace(toJSON(cols)) == flatten_whitespace(expected));
  }
}
TEST_CASE("simple1") {
  std::string testDataDir =
      std::string(getenv("RDBASE")) +
      std::string("/Code/GraphMol/RGroupDecomposition/test_data/");
  std::string fName = testDataDir + "simple1.sdf";
  SDMolSupplier suppl(fName);
  std::vector<ROMOL_SPTR> cores(1);
  std::vector<ROMOL_SPTR> mols;
  initDataset(suppl, cores.front(), mols);
  SECTION("defaults") {
    RGroupRows rows;
    auto n = RGroupDecompose(cores, mols, rows);
    CHECK(n == mols.size());
    CHECK(rows.size() == mols.size());
    CHECK(flatten_whitespace(toJSON(rows)) ==
          flatten_whitespace(
              readReferenceData(testDataDir + "simple1.out1.json")));
  }
  SECTION("no symmetrization") {
    RGroupRows rows;
    RGroupDecompositionParameters ps;
    ps.matchingStrategy = RGroupMatching::NoSymmetrization;
    auto n = RGroupDecompose(cores, mols, rows, nullptr, ps);
    CHECK(n == mols.size());
    CHECK(rows.size() == mols.size());
    CHECK(flatten_whitespace(toJSON(rows)) ==
          flatten_whitespace(
              readReferenceData(testDataDir + "simple1.out2.json")));
  }
}

TEST_CASE("simple2 with specified R groups") {
  std::string testDataDir =
      std::string(getenv("RDBASE")) +
      std::string("/Code/GraphMol/RGroupDecomposition/test_data/");
  std::string fName = testDataDir + "simple2.sdf";
  SDMolSupplier suppl(fName);
  std::vector<ROMOL_SPTR> cores(1);
  std::vector<ROMOL_SPTR> mols;
  initDataset(suppl, cores.front(), mols);
  SECTION("defaults") {
    RGroupRows rows;
    auto n = RGroupDecompose(cores, mols, rows);
    CHECK(n == mols.size());
    CHECK(rows.size() == mols.size());
    CHECK(flatten_whitespace(toJSON(rows)) ==
          flatten_whitespace(
              readReferenceData(testDataDir + "simple2.out1.json")));
  }
  SECTION("only match at r groups") {
    RGroupRows rows;
    RGroupDecompositionParameters ps;
    ps.onlyMatchAtRGroups = true;
    std::vector<unsigned> unmatched;
    auto n = RGroupDecompose(cores, mols, rows, &unmatched, ps);
    CHECK(n == 2);
    CHECK(rows.size() == n);
    CHECK(unmatched.size() == mols.size() - n);
    CHECK(unmatched[0] == 2);
    CHECK(flatten_whitespace(toJSON(rows)) ==
          flatten_whitespace(
              readReferenceData(testDataDir + "simple2.out2.json")));
  }
}

TEST_CASE("jm7b00306 Snippet") {
  std::string testDataDir =
      std::string(getenv("RDBASE")) +
      std::string("/Code/GraphMol/RGroupDecomposition/test_data/");
  std::string fName = testDataDir + "jm7b00306.excerpt.sdf";
  SDMolSupplier suppl(fName);
  std::vector<ROMOL_SPTR> cores(1);
  std::vector<ROMOL_SPTR> mols;
  initDataset(suppl, cores.front(), mols);
  SECTION("defaults") {
    RGroupRows rows;
    std::vector<unsigned> unmatched;
    auto n = RGroupDecompose(cores, mols, rows, &unmatched);
    CHECK(n == mols.size() - 1);
    CHECK(rows.size() == n);
    // there is one structure in there that doesn't match the core
    CHECK(unmatched.size() == mols.size() - n);
    CHECK(unmatched[0] == 1);
    CHECK(flatten_whitespace(toJSON(rows)) ==
          flatten_whitespace(
              readReferenceData(testDataDir + "jm7b00306.excerpt.out1.json")));
  }
}

TEST_CASE("jm200186n Snippet") {
  std::string testDataDir =
      std::string(getenv("RDBASE")) +
      std::string("/Code/GraphMol/RGroupDecomposition/test_data/");
  std::string fName = testDataDir + "jm200186n.excerpt.sdf";
  SDMolSupplier suppl(fName);
  std::vector<ROMOL_SPTR> cores(1);
  std::vector<ROMOL_SPTR> mols;
  initDataset(suppl, cores.front(), mols);
  SECTION("defaults") {
    RGroupRows rows;
    std::vector<unsigned> unmatched;
    auto n = RGroupDecompose(cores, mols, rows, &unmatched);
    CHECK(n == mols.size() - 1);
    CHECK(rows.size() == n);
    // there is one structure in there that doesn't match the core
    CHECK(unmatched.size() == mols.size() - n);
    CHECK(unmatched[0] == 3);
    CHECK(flatten_whitespace(toJSON(rows)) ==
          flatten_whitespace(
              readReferenceData(testDataDir + "jm200186n.excerpt.out1.json")));
  }
}
