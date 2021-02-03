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
  boost::algorithm::trim_all_if(res, boost::is_any_of(" \t\n"));
  return res;
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
    CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON([
  {
    "Core":"Cc1cccc([*:2])c1[*:1]",
    "R1":"CO[*:1]",
    "R2":"[H][*:2]"
  },
  {
    "Core":"Cc1cccc([*:2])c1[*:1]",
    "R1":"CO[*:1]",
    "R2":"[H][*:2]"
  },
  {
    "Core":"Cc1cccc([*:2])c1[*:1]",
    "R1":"[H][*:1]",
    "R2":"CO[*:2]"
  }
  ])JSON"));
  }
  SECTION("no symmetrization") {
    RGroupRows rows;
    RGroupDecompositionParameters ps;
    ps.matchingStrategy = RGroupMatching::NoSymmetrization;
    auto n = RGroupDecompose(cores, mols, rows, nullptr, ps);
    CHECK(n == mols.size());
    CHECK(rows.size() == mols.size());
    CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON([
  {
  "Core":"Cc1c([*:3])ccc([*:2])c1[*:1]",
  "R1":"[H][*:1]",
  "R2":"[H][*:2]",
  "R3":"CO[*:3]"
  },
  {
  "Core":"Cc1c([*:3])ccc([*:2])c1[*:1]",
  "R1":"CO[*:1]",
  "R2":"[H][*:2]",
  "R3":"[H][*:3]"
  },
  {
  "Core":"Cc1c([*:3])ccc([*:2])c1[*:1]",
  "R1":"[H][*:1]",
  "R2":"CO[*:2]",
  "R3":"[H][*:3]"
  }
  ])JSON"));
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
    CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON([
  {
  "Core":"Cc1cc([*:1])ccc1[*:3]",
  "R1":"[H][*:1]",
  "R3":"CO[*:3]"
  },
  {
  "Core":"Cc1cc([*:1])ccc1[*:3]",
  "R1":"[H][*:1]",
  "R3":"CO[*:3]"
  },
  {
  "Core":"Cc1cc([*:1])ccc1[*:3]",
  "R1":"CO[*:1]",
  "R3":"[H][*:3]"
  }
  ])JSON"));
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
    CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON([
  {
  "Core":"Cc1ccccc1[*:3]",
  "R3":"CO[*:3]"
  },
  {
  "Core":"Cc1ccccc1[*:3]",
  "R3":"CO[*:3]"
  }
  ])JSON"));
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
    CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON([
  {
  "Core":"Cn1c(=O)n(C)c2cc([*:3])c(N3C(=O)c4cccc5c([*:1])c([*:2])cc(c45)C3=O)cc21",
  "R1":"OCCC[*:1]",
  "R3":"C[*:3]"
  },
  {
  "Core":"Cn1c(=O)n(C)c2cc([*:3])c(N3C(=O)c4cccc5c([*:1])c([*:2])cc(c45)C3=O)cc21",
  "R1":"[H][*:1]",
  "R3":"[H][*:3]"
  },
  {
  "Core":"Cn1c(=O)n(C)c2cc([*:3])c(N3C(=O)c4cccc5c([*:1])c([*:2])cc(c45)C3=O)cc21",
  "R1":"[H][*:1]",
  "R3":"C[*:3]"
  },
  {
  "Core":"Cn1c(=O)n(C)c2cc([*:3])c(N3C(=O)c4cccc5c([*:1])c([*:2])cc(c45)C3=O)cc21",
  "R1":"[H][*:1]",
  "R3":"CO[*:3]"
  },
  {
  "Core":"Cn1c(=O)n(C)c2cc([*:3])c(N3C(=O)c4cccc5c([*:1])c([*:2])cc(c45)C3=O)cc21",
  "R1":"[H][*:1]",
  "R3":"CN(C)[*:3]"
  }
  ])JSON"));
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
    CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON([
  {
  "Core":"c1c([*:3])cc([*:5])c(O[*:1])c1-c1cc([*:4])c(O[*:2])c([*:6])c1",
  "R2":"C[*:2]",
  "R3":"C=CC[*:3]",
  "R4":"CC(Br)C[*:4]",
  "R6":"[H][*:6]"
  },
  {
  "Core":"c1c([*:3])cc([*:5])c(O[*:1])c1-c1cc([*:4])c(O[*:2])c([*:6])c1",
  "R2":"C[*:2]",
  "R3":"CC(Br)C[*:3]",
  "R4":"C=CC[*:4]",
  "R6":"[H][*:6]"
  },
  {
  "Core":"c1c([*:3])cc([*:5])c(O[*:1])c1-c1cc([*:4])c(O[*:2])c([*:6])c1",
  "R2":"[H][*:2]",
  "R3":"C=CC[*:3]",
  "R4":"C=CC[*:4]",
  "R6":"[H][*:6]"
  },
  {
  "Core":"c1c([*:3])cc([*:5])c(O[*:1])c1-c1cc([*:4])c(O[*:2])c([*:6])c1",
  "R2":"C[*:2]",
  "R3":"C=CC[*:3]",
  "R4":"C=CC[*:4]",
  "R6":"[H][*:6]"
  },
  {
  "Core":"c1c([*:3])cc([*:5])c(O[*:1])c1-c1cc([*:4])c(O[*:2])c([*:6])c1",
  "R2":"C[*:2]",
  "R3":"C=CC[*:3]",
  "R4":"C=CC[*:4]",
  "R6":"F[*:6]"
  }
  ])JSON"));
  }
}
