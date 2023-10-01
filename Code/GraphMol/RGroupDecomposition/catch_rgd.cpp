//
//  Copyright (C) 2021 Greg Landrum and other RDKit contributors
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
  boost::algorithm::trim_fill_if(res, "", boost::is_any_of(" \t\r\n"));
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
        "Core": "Cc1cccc([*:1])c1[*:2]",
        "R1": "[H][*:1]",
        "R2": "CO[*:2]"
    },
    {
        "Core": "Cc1cccc([*:1])c1[*:2]",
        "R1": "[H][*:1]",
        "R2": "CO[*:2]"
    },
    {
        "Core": "Cc1cccc([*:1])c1[*:2]",
        "R1": "CO[*:1]",
        "R2": "[H][*:2]"
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
    "Cc1cccc([*:1])c1[*:2]",
    "Cc1cccc([*:1])c1[*:2]",
    "Cc1cccc([*:1])c1[*:2]"
  ],
  "R1": [
    "[H][*:1]",
    "[H][*:1]",
    "CO[*:1]"
  ],
  "R2": [
    "CO[*:2]",
    "CO[*:2]",
    "[H][*:2]"
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

TEST_CASE("simple3 with user labels on aromatic N") {
  std::string testDataDir =
      std::string(getenv("RDBASE")) +
      std::string("/Code/GraphMol/RGroupDecomposition/test_data/");
  std::string fName = testDataDir + "simple3.sdf";
  SDMolSupplier suppl(fName);
  std::vector<ROMOL_SPTR> cores(1);
  std::vector<ROMOL_SPTR> mols;
  initDataset(suppl, cores.front(), mols);
  SECTION("defaults (allH labels and R-groups are removed)") {
    RGroupRows rows;
    auto n = RGroupDecompose(cores, mols, rows);
    CHECK(n == mols.size());
    CHECK(rows.size() == mols.size());
    CHECK(flatten_whitespace(toJSON(rows)) ==
          flatten_whitespace(
              readReferenceData(testDataDir + "simple3.out1.json")));
  }
  SECTION("removeAllHydrogenRGroups = false (as defaults)") {
    RGroupRows rows;
    RGroupDecompositionParameters ps;
    ps.removeAllHydrogenRGroups = false;
    std::vector<unsigned> unmatched;
    auto n = RGroupDecompose(cores, mols, rows, &unmatched, ps);
    CHECK(n == mols.size());
    CHECK(rows.size() == mols.size());
    CHECK(unmatched.empty());
    CHECK(flatten_whitespace(toJSON(rows)) ==
          flatten_whitespace(
              readReferenceData(testDataDir + "simple3.out2.json")));
  }
  SECTION("removeAllHydrogenRGroupsAndLabels = false (allH labels retained)") {
    RGroupRows rows;
    RGroupDecompositionParameters ps;
    ps.removeAllHydrogenRGroupsAndLabels = false;
    std::vector<unsigned> unmatched;
    auto n = RGroupDecompose(cores, mols, rows, &unmatched, ps);
    CHECK(n == mols.size());
    CHECK(rows.size() == mols.size());
    CHECK(unmatched.empty());
    CHECK(flatten_whitespace(toJSON(rows)) ==
          flatten_whitespace(
              readReferenceData(testDataDir + "simple3.out3.json")));
  }
  SECTION(
      "removeAllHydrogenRGroupsAndLabels = false, removeAllHydrogenRGroups = "
      "false (allH labels and R-groups are retained)") {
    RGroupRows rows;
    RGroupDecompositionParameters ps;
    ps.removeAllHydrogenRGroups = false;
    ps.removeAllHydrogenRGroupsAndLabels = false;
    std::vector<unsigned> unmatched;
    auto n = RGroupDecompose(cores, mols, rows, &unmatched, ps);
    CHECK(n == mols.size());
    CHECK(rows.size() == mols.size());
    CHECK(unmatched.empty());
    CHECK(flatten_whitespace(toJSON(rows)) ==
          flatten_whitespace(
              readReferenceData(testDataDir + "simple3.out4.json")));
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

std::vector<ROMOL_SPTR> smisToMols(const std::vector<std::string> &smis) {
  std::vector<ROMOL_SPTR> mols;
  for (const auto &smi : smis) {
    auto m = SmilesToMol(smi);
    assert(m);
    mols.emplace_back(m);
  }
  return mols;
}

TEST_CASE("substructure parameters and RGD: chirality") {
  std::vector<std::string> smis = {"C1CN[C@H]1F", "C1CN[C@]1(O)F",
                                   "C1CN[C@@H]1F", "C1CN[CH]1F"};
  auto mols = smisToMols(smis);
  std::vector<std::string> csmis = {"C1CNC1[*:1]"};
  auto cores = smisToMols(csmis);
  std::vector<std::string> csmis2 = {"C1CN[C@H]1[*:1]"};
  auto chiral_cores = smisToMols(csmis2);
  SECTION("defaults") {
    RGroupRows rows;
    std::vector<unsigned> unmatched;
    RGroupDecompositionParameters params;
    params.allowMultipleRGroupsOnUnlabelled = true;
    {
      auto n = RGroupDecompose(cores, mols, rows, &unmatched, params);
      CHECK(n == mols.size());
      CHECK(rows.size() == n);
      CHECK(unmatched.empty());
      CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON(
[
  {
    "Core":"C1C[C@@]([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"[H][*:2]"
  },
  {
    "Core":"C1C[C@]([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"O[*:2]"
  },
  {
    "Core":"C1C[C@]([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"[H][*:2]"
  },
  {
    "Core":"C1CC([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"[H][*:2]"
  }
]
    )JSON"));
    }
    {
      auto n = RGroupDecompose(chiral_cores, mols, rows, &unmatched, params);
      CHECK(n == mols.size() - 1);
      CHECK(rows.size() == n);
      CHECK(unmatched.size() == 1);
      CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON(
[
  {
    "Core":"C1C[C@@]([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"[H][*:2]"
  },
  {
    "Core":"C1C[C@@]([*:1])([*:2])N1",
    "R1":"O[*:1]",
    "R2":"F[*:2]"
  },
  {
    "Core":"C1C[C@@]([*:1])([*:2])N1",
    "R1":"[H][*:1]",
    "R2":"F[*:2]"
  }
]
    )JSON"));
    }
  }

  SECTION("not using chirality") {
    // this time both cores return the same thing and stereo information is
    // removed from the chiral cores
    RGroupRows rows;
    std::vector<unsigned> unmatched;
    RGroupDecompositionParameters params;
    params.allowMultipleRGroupsOnUnlabelled = true;
    params.substructmatchParams.useChirality = false;
    {
      auto n = RGroupDecompose(cores, mols, rows, &unmatched, params);
      CHECK(n == mols.size());
      CHECK(rows.size() == n);
      CHECK(unmatched.empty());
      CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON(
[
  {
    "Core":"C1CC([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"[H][*:2]"
  },
  {
    "Core":"C1CC([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"O[*:2]"
  },
  {
    "Core":"C1CC([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"[H][*:2]"
  },
  {
    "Core":"C1CC([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"[H][*:2]"
  }
]
    )JSON"));
    }
    {
      auto n = RGroupDecompose(chiral_cores, mols, rows, &unmatched, params);
      CHECK(n == mols.size());
      CHECK(rows.size() == n);
      CHECK(unmatched.empty());
      CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON(
[
  {
    "Core":"C1CC([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"[H][*:2]"
  },
  {
    "Core":"C1CC([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"O[*:2]"
  },
  {
    "Core":"C1CC([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"[H][*:2]"
  },
  {
    "Core":"C1CC([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"[H][*:2]"
  }
]
    )JSON"));
    }
  }
}

TEST_CASE("substructure parameters and RGD: enhanced stereo") {
  std::vector<std::string> smis = {"F[C@H]1CCN1 |&1:1|", "C1CN[C@]1(O)F |&1:3|",
                                   "C1CN[C@@H]1F |&1:3|", "Cl[C@H]1CCN1 |o1:1|",
                                   "C1CN[CH]1F"};
  auto mols = smisToMols(smis);
  std::vector<std::string> csmis = {"C1CN[C@H]1[*:1] |&1:3|"};
  auto cores = smisToMols(csmis);
  std::vector<std::string> csmis2 = {"C1CN[C@H]1[*:1] |o1:3|"};
  auto cores2 = smisToMols(csmis2);
  SECTION("defaults: no enhanced stereo") {
    RGroupRows rows;
    std::vector<unsigned> unmatched;
    RGroupDecompositionParameters params;
    params.allowMultipleRGroupsOnUnlabelled = true;
    {
      auto n = RGroupDecompose(cores, mols, rows, &unmatched, params);
      CHECK(n == mols.size() - 1);
      CHECK(rows.size() == n);
      CHECK(unmatched.size() == mols.size() - n);
      // std::cerr << toJSON(rows) << std::endl;

      // the core output no longer is SMARTS as the core output is the portion
      // of the target that matches the core query.
      CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON(
[
  {
    "Core":"C1C[C@@]([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"[H][*:2]"
  },
  {
    "Core":"C1C[C@@]([*:1])([*:2])N1",
    "R1":"O[*:1]",
    "R2":"F[*:2]"
  },
  {
    "Core":"C1C[C@@]([*:1])([*:2])N1",
    "R1":"[H][*:1]",
    "R2":"F[*:2]"
  },
  {
    "Core":"C1C[C@@]([*:1])([*:2])N1",
    "R1":"Cl[*:1]",
    "R2":"[H][*:2]"
  }
]
    )JSON"));
    }
    {
      auto n = RGroupDecompose(cores2, mols, rows, &unmatched, params);
      CHECK(n == mols.size() - 1);
      CHECK(rows.size() == n);
      CHECK(unmatched.size() == 1);
      // std::cerr << toJSON(rows) << std::endl;

      // the core output no longer is SMARTS as the core output is the portion
      // of the target that matches the core query.
      CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON(
[
  {
    "Core":"C1C[C@@]([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"[H][*:2]"
  },
  {
    "Core":"C1C[C@@]([*:1])([*:2])N1",
    "R1":"O[*:1]",
    "R2":"F[*:2]"
  },
  {
    "Core":"C1C[C@@]([*:1])([*:2])N1",
    "R1":"[H][*:1]",
    "R2":"F[*:2]"
  },
  {
    "Core":"C1C[C@@]([*:1])([*:2])N1",
    "R1":"Cl[*:1]",
    "R2":"[H][*:2]"
  }
]
    )JSON"));
    }
  }

  SECTION("using enhanced stereo") {
    RGroupRows rows;
    std::vector<unsigned> unmatched;
    RGroupDecompositionParameters params;
    params.allowMultipleRGroupsOnUnlabelled = true;
    params.substructmatchParams.useEnhancedStereo = true;
    {
      auto n = RGroupDecompose(cores, mols, rows, &unmatched, params);
      CHECK(n == mols.size() - 2);
      CHECK(rows.size() == n);
      CHECK(unmatched.size() == mols.size() - n);
      // std::cerr << toJSON(rows) << std::endl;
      // the core output no longer is SMARTS as the core output is the portion
      // of the target that matches the core query.
      CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON(
[
  {
    "Core":"C1C[C@@]([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"[H][*:2]"
  },
  {
    "Core":"C1C[C@]([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"O[*:2]"
  },
  {
    "Core":"C1C[C@]([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"[H][*:2]"
  }
]
    )JSON"));
    }
    {
      auto n = RGroupDecompose(cores2, mols, rows, &unmatched, params);
      CHECK(n == mols.size() - 1);
      CHECK(rows.size() == n);
      CHECK(unmatched.size() == 1);
      // std::cerr << toJSON(rows) << std::endl;
      // the core output no longer is SMARTS as the core output is the portion
      // of the target that matches the core query.
      CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON(
[
  {
    "Core":"C1C[C@@]([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"[H][*:2]"
  },
  {
    "Core":"C1C[C@]([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"O[*:2]"
  },
  {
    "Core":"C1C[C@]([*:1])([*:2])N1",
    "R1":"F[*:1]",
    "R2":"[H][*:2]"
  },
  {
    "Core":"C1C[C@@]([*:1])([*:2])N1",
    "R1":"Cl[*:1]",
    "R2":"[H][*:2]"
  }
]
    )JSON"));
    }
  }
}

TEST_CASE("github4809: ring double bonds written as crossed bonds after RGD") {
  std::vector<std::string> smis = {"C1C=CCC2=C1C=CC=N2"};
  auto mols = smisToMols(smis);
  std::vector<std::string> csmis = {"c1ccnc([*:1])c1[*:2]"};
  auto cores = smisToMols(csmis);
  SECTION("basics") {
    RGroupRows rows;
    {
      auto n = RGroupDecompose(cores, mols, rows);
      CHECK(n == mols.size());
      CHECK(rows.size() == n);
      auto r1 = rows[0]["R1"];
      auto mb = MolToV3KMolBlock(*r1);
      CHECK(mb.find("CFG=2") == std::string::npos);
    }
  }
}

TEST_CASE("rgroupLabelling") {
  std::vector<std::string> smis = {"C1CN[C@H]1F", "C1CN[C@]1(O)F",
                                   "C1CN[C@@H]1F", "C1CN[CH]1F"};
  auto mols = smisToMols(smis);
  std::vector<std::string> csmis = {"C1CNC1[*:1]"};
  auto cores = smisToMols(csmis);
  SECTION("Isotope") {
    RGroupRows rows;
    std::vector<unsigned> unmatched;
    RGroupDecompositionParameters params;
    params.rgroupLabelling = RGroupLabelling::Isotope;
    params.allowMultipleRGroupsOnUnlabelled = true;
    {
      auto n = RGroupDecompose(cores, mols, rows, &unmatched, params);
      CHECK(n == mols.size());
      CHECK(rows.size() == n);
      CHECK(unmatched.empty());
      CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON(
[
  {
    "Core": "[1*][C@@]1([2*])CCN1",
    "R1":"[1*]F",
    "R2":"[2*][H]"
  },
  {
    "Core": "[1*][C@]1([2*])CCN1",
    "R1":"[1*]F",
    "R2":"[2*]O"
  },
  {
    "Core":"[1*][C@]1([2*])CCN1",
    "R1":"[1*]F",
    "R2":"[2*][H]"
  },
  {
    "Core":"[1*]C1([2*])CCN1",
    "R1":"[1*]F",
    "R2":"[2*][H]"
  }
]
    )JSON"));
    }
  }
  SECTION("RGroup") {
    RGroupRows rows;
    std::vector<unsigned> unmatched;
    RGroupDecompositionParameters params;
    params.rgroupLabelling = RGroupLabelling::MDLRGroup;
    params.allowMultipleRGroupsOnUnlabelled = true;
    {
      auto n = RGroupDecompose(cores, mols, rows, &unmatched, params);
      CHECK(n == mols.size());
      CHECK(rows.size() == n);
      CHECK(unmatched.empty());
      // in this case the labels don't show up in the output SMILES
      // Presumably the dummy atoms are no longer distinguishable without
      // the isotope labels as the smiles no longer contains chiralty.
      // Chirality is present in the core SMARTS
      CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON(
[
  {
    "Core":"*C1(*)CCN1",
    "R1":"*F",
    "R2":"*[H]"
  },
  {
    "Core":"*C1(*)CCN1",
    "R1":"*F",
    "R2":"*O"
  },
  {
    "Core":"*C1(*)CCN1",
    "R1":"*F",
    "R2":"*[H]"
  },
  {
    "Core":"*C1(*)CCN1",
    "R1":"*F",
    "R2":"*[H]"
  }
]
    )JSON"));
    }
  }
  SECTION("Isotope|Map") {
    RGroupRows rows;
    std::vector<unsigned> unmatched;
    RGroupDecompositionParameters params;
    params.allowMultipleRGroupsOnUnlabelled = true;
    params.rgroupLabelling =
        RGroupLabelling::Isotope | RGroupLabelling::AtomMap;
    {
      auto n = RGroupDecompose(cores, mols, rows, &unmatched, params);
      CHECK(n == mols.size());
      CHECK(rows.size() == n);
      CHECK(unmatched.empty());
      CHECK(flatten_whitespace(toJSON(rows)) == flatten_whitespace(R"JSON(

[
  {
    "Core":"C1C[C@@]([1*:1])([2*:2])N1",
    "R1":"F[1*:1]",
    "R2":"[H][2*:2]"
  },
  {
    "Core":"C1C[C@]([1*:1])([2*:2])N1",
    "R1":"F[1*:1]",
    "R2":"O[2*:2]"
  },
  {
    "Core":"C1C[C@]([1*:1])([2*:2])N1",
    "R1":"F[1*:1]",
    "R2":"[H][2*:2]"
  },
  {
    "Core":"C1CC([1*:1])([2*:2])N1",
    "R1":"F[1*:1]",
    "R2":"[H][2*:2]"
  }
]
    )JSON"));
    }
  }
}

TEST_CASE("MDL R labels from original core") {
  std::vector<std::string> smis = {"C1CN[C@H]1F", "C1CN[C@]1(O)F",
                                   "C1CN[C@@H]1F", "C1CN[CH]1F"};
  auto mols = smisToMols(smis);
  std::vector<std::string> csmis = {"[*]C1CCN1 |$_R1;;;;$|"};
  auto cores = smisToMols(csmis);
  SECTION("Map") {
    RGroupRows rows;
    std::vector<unsigned> unmatched;
    RGroupDecompositionParameters params;
    params.allowMultipleRGroupsOnUnlabelled = true;
    params.rgroupLabelling = RGroupLabelling::AtomMap;
    {
      auto n = RGroupDecompose(cores, mols, rows, &unmatched, params);
      CHECK(n == mols.size());
      CHECK(rows.size() == n);
      CHECK(unmatched.empty());
      CHECK(rows[0]["Core"]->getAtomWithIdx(4)->getAtomicNum() == 0);
      CHECK(!rows[0]["Core"]->getAtomWithIdx(4)->hasProp(
          common_properties::dummyLabel));
      CHECK(rows[0]["Core"]->getAtomWithIdx(5)->getAtomicNum() == 0);
      CHECK(!rows[0]["Core"]->getAtomWithIdx(5)->hasProp(
          common_properties::dummyLabel));
    }
  }
  SECTION("Map | MDL") {
    RGroupRows rows;
    std::vector<unsigned> unmatched;
    RGroupDecompositionParameters params;
    params.allowMultipleRGroupsOnUnlabelled = true;
    params.rgroupLabelling =
        RGroupLabelling::AtomMap | RGroupLabelling::MDLRGroup;
    {
      auto n = RGroupDecompose(cores, mols, rows, &unmatched, params);
      CHECK(n == mols.size());
      CHECK(rows.size() == n);
      CHECK(unmatched.empty());
      CHECK(rows[0]["Core"]->getAtomWithIdx(4)->getAtomicNum() == 0);
      CHECK(rows[0]["Core"]->getAtomWithIdx(4)->hasProp(
          common_properties::dummyLabel));
      CHECK(rows[0]["Core"]->getAtomWithIdx(5)->getAtomicNum() == 0);
      CHECK(rows[0]["Core"]->getAtomWithIdx(5)->hasProp(
          common_properties::dummyLabel));
    }
  }
}

TEST_CASE("Mol matches core") {
  auto core = "[*:1]c1[!#1]([*:2])cc([*:3])n([*:4])c(=O)1"_smarts;
  auto cmol = "Clc1c(C)cc(F)n(CC)c(=O)1"_smiles;
  auto nmol = "Clc1ncc(F)n(CC)c(=O)1"_smiles;
  auto smol = "Clc1ncc(F)n(CC)c(=S)1"_smiles;
  RGroupDecompositionParameters params;
  params.onlyMatchAtRGroups = true;
  RGroupDecomposition decomp(*core, params);
  CHECK(decomp.getMatchingCoreIdx(*cmol) == 0);
  CHECK(decomp.getMatchingCoreIdx(*nmol) == 0);
  CHECK(decomp.getMatchingCoreIdx(*smol) == -1);
  std::vector<MatchVectType> matches;
  CHECK(decomp.getMatchingCoreIdx(*cmol, &matches) == 0);
  CHECK(matches.size() == 1);
  CHECK(matches.front().size() == core->getNumAtoms());
  CHECK(decomp.getMatchingCoreIdx(*nmol, &matches) == 0);
  CHECK(matches.size() == 1);
  CHECK(matches.front().size() == core->getNumAtoms() - 1);
  CHECK(decomp.getMatchingCoreIdx(*smol, &matches) == -1);
  CHECK(matches.empty());
  MolOps::addHs(*cmol);
  MolOps::addHs(*nmol);
  MatchVectType match;
  CHECK(SubstructMatch(*cmol, *core, match));
  CHECK(match.size() == core->getNumAtoms());
  match.clear();
  CHECK(!SubstructMatch(*nmol, *core, match));
}
