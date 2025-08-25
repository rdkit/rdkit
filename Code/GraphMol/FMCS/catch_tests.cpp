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

#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/FMCS/TwoMolMCSS.h>
#include <GraphMol/FMCS/MatchTable.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <catch2/catch_all.hpp>

using namespace RDKit;

// Put together atom and bond match tables using default rules.
FMCS::MatchTable makeAtomMatchTable(const ROMol &mol1, const ROMol &mol2,
                                    const MCSAtomCompareParameters &p) {
  FMCS::MatchTable amt(mol1.getNumAtoms(), mol2.getNumAtoms());
  for (const auto a1 : mol1.atoms()) {
    for (const auto a2 : mol2.atoms()) {
      amt.set(a1->getIdx(), a2->getIdx(),
              MCSAtomCompareElements(p, mol1, a1->getIdx(), mol2, a2->getIdx(),
                                     nullptr));
    }
  }
  return amt;
}

FMCS::MatchTable makeBondMatchTable(const ROMol &mol1, const ROMol &mol2,
                                    const MCSBondCompareParameters &p) {
  FMCS::MatchTable bmt(mol1.getNumBonds(), mol2.getNumBonds());
  for (const auto b1 : mol1.bonds()) {
    for (const auto b2 : mol2.bonds()) {
      bmt.set(b1->getIdx(), b2->getIdx(),
              MCSBondCompareOrder(p, mol1, b1->getIdx(), mol2, b2->getIdx(),
                                  nullptr));
    }
  }
  return bmt;
}

void printCliques(
    const std::vector<std::vector<std::pair<unsigned int, unsigned int>>>
        &maxCliques) {
  for (const auto &clique : maxCliques) {
    std::cout << clique.size() << " : ";
    for (const auto c : clique) {
      std::cout << "(" << c.first << ", " << c.second << "), ";
    }
    std::cout << std::endl;
  }
}

TEST_CASE("BK Test") {
  auto mol1 = v2::SmilesParse::MolFromSmiles("CCC");
  REQUIRE(mol1);
  auto mol2 = v2::SmilesParse::MolFromSmiles("CC(C)C");
  REQUIRE(mol2);
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
  MCSParameters params;
  auto amt = makeAtomMatchTable(*mol1, *mol2, params.AtomCompareParameters);
  auto bmt = makeBondMatchTable(*mol1, *mol2, params.BondCompareParameters);
  TwoMolMCSS(*mol1, *mol2, 3, amt, bmt, false, 0, maxCliques);
  REQUIRE(maxCliques.size() == 6);
  CHECK(maxCliques.front().size() == 3);
  CHECK(maxCliques.front() ==
        std::vector<std::pair<unsigned int, unsigned int>>{
            {0, 0}, {1, 1}, {2, 2}});
  CHECK(maxCliques.back().size() == 3);
}

TEST_CASE("MedChemica Base Test") {
  auto mol1 = v2::SmilesParse::MolFromSmiles("c1ccc2c(c1)c(ncn2)Nc3ccc(cc3)Cl");
  REQUIRE(mol1);
  auto mol2 =
      v2::SmilesParse::MolFromSmiles("c1ccc2c(c1)c(ncn2)Nc3ccc(cc3)C#N");
  REQUIRE(mol2);
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
  unsigned int minMCSSSize(0.85 *
                           std::max(mol1->getNumAtoms(), mol2->getNumAtoms()));
  MCSParameters params;
  auto amt = makeAtomMatchTable(*mol1, *mol2, params.AtomCompareParameters);
  auto bmt = makeBondMatchTable(*mol1, *mol2, params.BondCompareParameters);
  TwoMolMCSS(*mol1, *mol2, minMCSSSize, amt, bmt, true, 0, maxCliques);
  CHECK(maxCliques.size() == 3);
  CHECK(maxCliques.front().size() == 17);
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> expCliques{
      {{0, 0},
       {1, 1},
       {2, 2},
       {3, 3},
       {4, 4},
       {5, 5},
       {6, 6},
       {7, 7},
       {8, 8},
       {9, 9},
       {10, 10},
       {11, 11},
       {12, 12},
       {13, 13},
       {14, 14},
       {15, 15},
       {16, 16}},
      {{0, 0},
       {1, 1},
       {2, 2},
       {3, 3},
       {4, 4},
       {5, 5},
       {6, 6},
       {7, 7},
       {8, 8},
       {9, 9},
       {10, 10},
       {11, 11},
       {12, 12},
       {13, 13},
       {14, 14},
       {15, 17},
       {16, 16}},
      {{0, 0},
       {1, 1},
       {2, 2},
       {3, 3},
       {4, 4},
       {5, 5},
       {6, 6},
       {7, 7},
       {8, 8},
       {9, 9},
       {10, 10},
       {11, 11},
       {12, 12},
       {13, 17},
       {14, 14},
       {15, 15},
       {16, 16}}};
  CHECK(maxCliques == expCliques);
  TwoMolMCSS(*mol1, *mol2, minMCSSSize, amt, bmt, false, 0, maxCliques);
  CHECK(maxCliques.size() == 30);
  CHECK(maxCliques.front().size() == 17);

  // Now compare new FMCS with old.
  std::vector<ROMOL_SPTR> mols;
  mols.push_back(ROMOL_SPTR(mol1.release()));
  mols.push_back(ROMOL_SPTR(mol2.release()));
  MCSResult res;
  params.MinMCSSSize = static_cast<unsigned int>(
      0.85 * std::max(mols[0]->getNumAtoms(), mols[1]->getNumAtoms()));
  {
    params.FastInitialSeed = true;
    res = findMCS(mols, &params);
    CHECK(res.NumAtoms == 17);
    CHECK(res.NumBonds == 19);
  }
  {
    params.FastInitialSeed = false;
    res = findMCS(mols, &params);
    CHECK(res.NumAtoms == 17);
    CHECK(res.NumBonds == 19);
  }
}

TEST_CASE("FMCS Slow 1") {
  // Spirocycle/open chain
  auto mol1 = v2::SmilesParse::MolFromSmiles(
      "CN1C[C@]2(CCN(C2=O)c2cncc3ccccc23)c2cc(Cl)ccc2C1=O");
  REQUIRE(mol1);
  auto mol2 = v2::SmilesParse::MolFromSmiles(
      "CN1C[C@@H](C(=O)Nc2cncc3ccccc23)c2cc(Cl)ccc2C1=O");
  REQUIRE(mol2);
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
  unsigned int minMCSSSize(0.85 *
                           std::max(mol1->getNumAtoms(), mol2->getNumAtoms()));
  MCSParameters params;
  auto amt = makeAtomMatchTable(*mol1, *mol2, params.AtomCompareParameters);
  auto bmt = makeBondMatchTable(*mol1, *mol2, params.BondCompareParameters);
  TwoMolMCSS(*mol1, *mol2, minMCSSSize, amt, bmt, true, 0, maxCliques);
  REQUIRE(maxCliques.size() == 7);
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> expCliques{
      {{0, 0},   {1, 1},   {2, 2},   {3, 3},   {6, 6},   {7, 4},   {8, 5},
       {9, 7},   {10, 8},  {11, 9},  {12, 10}, {13, 11}, {14, 12}, {15, 13},
       {16, 14}, {17, 15}, {18, 16}, {19, 17}, {20, 18}, {21, 19}, {22, 20},
       {23, 21}, {24, 22}, {25, 23}, {26, 24}, {27, 25}},
      {
          {0, 0},   {1, 1},   {2, 2},   {3, 3},   {4, 17},  {5, 18},  {6, 6},
          {7, 4},   {8, 5},   {9, 7},   {10, 8},  {11, 9},  {12, 10}, {13, 11},
          {14, 12}, {15, 13}, {16, 14}, {17, 15}, {18, 16}, {21, 19}, {22, 20},
          {23, 21}, {24, 22}, {25, 23}, {26, 24}, {27, 25},
      }};
  CHECK(maxCliques.front() == expCliques[0]);
  CHECK(maxCliques.back() == expCliques[1]);

  std::vector<ROMOL_SPTR> mols;
  mols.push_back(ROMOL_SPTR(mol1.release()));
  mols.push_back(ROMOL_SPTR(mol2.release()));
  MCSResult res;
  {
    params.FastInitialSeed = true;
    params.MinMCSSSize = static_cast<unsigned int>(
        0.85 * std::max(mols[0]->getNumAtoms(), mols[1]->getNumAtoms()));
    res = findMCS(mols, &params);
    CHECK(res.NumAtoms == 26);
    CHECK(res.NumBonds == 29);
  }
  {
    params.FastInitialSeed = false;
    res = findMCS(mols, &params);
    CHECK(res.NumAtoms == 26);
    CHECK(res.NumBonds == 29);
  }
}

TEST_CASE("FMCS Slow 2") {
  // spirocycle/spirocycle
  auto mol1 = v2::SmilesParse::MolFromSmiles(
      R"(CN1C[C@]2(CCN(C2=O)c2cncc3ccccc23)c2cc(Cl)ccc2C1=O)");
  REQUIRE(mol1);
  auto mol2 = v2::SmilesParse::MolFromSmiles(
      R"(CCN1C[C@]2(CCN(C2=O)c2cncc3ccccc23)c2cc(Cl)ccc2C1=O)");
  REQUIRE(mol2);
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
  unsigned int minMCSSSize(0.85 *
                           std::max(mol1->getNumAtoms(), mol2->getNumAtoms()));
  MCSParameters params;
  auto amt = makeAtomMatchTable(*mol1, *mol2, params.AtomCompareParameters);
  auto bmt = makeBondMatchTable(*mol1, *mol2, params.BondCompareParameters);
  TwoMolMCSS(*mol1, *mol2, minMCSSSize, amt, bmt, true, 0, maxCliques);
  REQUIRE(maxCliques.size() == 1);
  std::vector<std::pair<unsigned int, unsigned int>> expClique{
      {0, 1},   {1, 2},   {2, 3},   {3, 4},   {4, 5},   {5, 6},   {6, 7},
      {7, 8},   {8, 9},   {9, 10},  {10, 11}, {11, 12}, {12, 13}, {13, 14},
      {14, 15}, {15, 16}, {16, 17}, {17, 18}, {18, 19}, {19, 20}, {20, 21},
      {21, 22}, {22, 23}, {23, 24}, {24, 25}, {25, 26}, {26, 27}, {27, 28}};
  CHECK(maxCliques.front() == expClique);
  std::vector<ROMOL_SPTR> mols;
  mols.push_back(ROMOL_SPTR(mol1.release()));
  mols.push_back(ROMOL_SPTR(mol2.release()));
  MCSResult res;
  {
    params.FastInitialSeed = true;
    params.MinMCSSSize = static_cast<unsigned int>(
        0.85 * std::max(mols[0]->getNumAtoms(), mols[1]->getNumAtoms()));
    res = findMCS(mols, &params);
    CHECK(res.NumAtoms == 28);
    CHECK(res.NumBonds == 32);
  }
  {
    params.FastInitialSeed = false;
    res = findMCS(mols, &params);
    CHECK(res.NumAtoms == 28);
    CHECK(res.NumBonds == 32);
  }
}

TEST_CASE("FMCS Slow 3") {
  auto mol1 = v2::SmilesParse::MolFromSmiles(
      R"(CO\N=C(/C(=O)N[C@H]1[C@H]2SCC(C[n+]3cscc3C)=C(N2C1=O)C(O)=O)c1nc(N)sc1-c1ccccc1)");
  REQUIRE(mol1);
  auto mol2 = v2::SmilesParse::MolFromSmiles(
      R"(CO\N=C(/C(=O)N[C@H]1[C@H]2SCC(C[n+]3ccc4sccc4c3)=C(N2C1=O)C(O)=O)c1nc(N)sc1-c1ccccc1)");
  REQUIRE(mol2);
  MCSResult res;
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
  unsigned int minMCSSSize(0.85 *
                           std::max(mol1->getNumAtoms(), mol2->getNumAtoms()));
  MCSParameters params;
  auto amt = makeAtomMatchTable(*mol1, *mol2, params.AtomCompareParameters);
  auto bmt = makeBondMatchTable(*mol1, *mol2, params.BondCompareParameters);
  TwoMolMCSS(*mol1, *mol2, minMCSSSize, amt, bmt, true, 0, maxCliques);
  REQUIRE(maxCliques.size() == 4);
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> expCliques{
      {{0, 0},   {1, 1},   {2, 2},   {3, 3},   {4, 4},   {5, 5},
       {6, 6},   {7, 7},   {8, 8},   {9, 9},   {10, 10}, {11, 11},
       {12, 12}, {13, 13}, {14, 14}, {16, 20}, {17, 21}, {19, 22},
       {20, 23}, {21, 24}, {22, 25}, {23, 26}, {24, 27}, {25, 28},
       {26, 29}, {27, 30}, {28, 31}, {29, 32}, {30, 33}, {31, 34},
       {32, 35}, {33, 36}, {34, 37}, {35, 38}, {36, 39}, {37, 40}},
      {{0, 0},   {1, 1},   {2, 2},   {3, 3},   {4, 4},   {5, 5},
       {6, 6},   {7, 7},   {8, 8},   {9, 9},   {10, 10}, {11, 11},
       {12, 12}, {13, 13}, {14, 21}, {17, 14}, {18, 15}, {19, 22},
       {20, 23}, {21, 24}, {22, 25}, {23, 26}, {24, 27}, {25, 28},
       {26, 29}, {27, 30}, {28, 31}, {29, 32}, {30, 33}, {31, 34},
       {32, 35}, {33, 36}, {34, 37}, {35, 38}, {36, 39}, {37, 40}}};
  CHECK(maxCliques.front() == expCliques.front());
  CHECK(maxCliques.back() == expCliques.back());
  std::vector<ROMOL_SPTR> mols;
  mols.push_back(ROMOL_SPTR(mol1.release()));
  mols.push_back(ROMOL_SPTR(mol2.release()));
}

TEST_CASE("Split by Rascal") {
  // Rascal gives a 2-fragment MCES, 9 atoms and 3 atoms. This checks
  // that the new code gives the 10-atom single fragment MCSS.
  auto mol1 = v2::SmilesParse::MolFromSmiles("C=CCCC=CCCC1CNCCC1");
  REQUIRE(mol1);
  auto mol2 = v2::SmilesParse::MolFromSmiles("C=CCC=CCCCC1CNCCC1");
  REQUIRE(mol2);
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
  MCSParameters params;
  auto amt = makeAtomMatchTable(*mol1, *mol2, params.AtomCompareParameters);
  auto bmt = makeBondMatchTable(*mol1, *mol2, params.BondCompareParameters);
  TwoMolMCSS(*mol1, *mol2, 10, amt, bmt, true, 0, maxCliques);
  REQUIRE(maxCliques.size() == 1);
  CHECK(maxCliques.front() ==
        std::vector<std::pair<unsigned int, unsigned int>>{{2, 1},
                                                           {3, 2},
                                                           {4, 3},
                                                           {5, 4},
                                                           {6, 5},
                                                           {7, 6},
                                                           {8, 7},
                                                           {11, 12},
                                                           {12, 13},
                                                           {13, 8}});
}

TEST_CASE("Order Dependence") {
  // Once upon a time these gave different answers in the different orders.
  std::vector<ROMOL_SPTR> mols{
      R"(O=C(Cc1ccccc1)Nc1cccc(-c2nc3sccn3c2-c2ccnc(Nc3cccc(N4CCOCC4)c3)n2)c1)"_smiles,
      R"(CN(C)CCc1cccc(Nc2nccc(-c3c(-c4cccc(NC(=O)Cc5ccccc5)c4)nc4sccn34)n2)c1)"_smiles,
  };
  MCSParameters params;
  {
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
    auto amt =
        makeAtomMatchTable(*mols[0], *mols[1], params.AtomCompareParameters);
    auto bmt =
        makeBondMatchTable(*mols[0], *mols[1], params.BondCompareParameters);
    TwoMolMCSS(*mols[0], *mols[1], 11, amt, bmt, true, 0, maxCliques);
    REQUIRE(maxCliques.size() == 1);
    std::vector<std::pair<unsigned int, unsigned int>> expClique{
        {0, 25},  {1, 24},  {2, 26},  {3, 27},  {4, 28},  {5, 29},  {6, 30},
        {7, 31},  {8, 32},  {9, 23},  {10, 22}, {11, 21}, {12, 20}, {13, 19},
        {14, 18}, {15, 17}, {16, 34}, {17, 35}, {18, 36}, {19, 37}, {20, 38},
        {21, 39}, {22, 16}, {23, 15}, {24, 14}, {25, 13}, {26, 12}, {27, 11},
        {28, 10}, {29, 9},  {30, 41}, {31, 5},  {32, 4},  {33, 3},  {34, 1},
        {35, 0},  {39, 2},  {40, 8},  {41, 40}, {42, 33}};
    CHECK(maxCliques.front() == expClique);
  }
  {
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
    auto amt =
        makeAtomMatchTable(*mols[1], *mols[0], params.AtomCompareParameters);
    auto bmt =
        makeBondMatchTable(*mols[1], *mols[0], params.BondCompareParameters);
    TwoMolMCSS(*mols[1], *mols[0], 11, amt, bmt, true, 0, maxCliques);
    REQUIRE(maxCliques.size() == 1);
    std::vector<std::pair<unsigned int, unsigned int>> expClique{
        {0, 35},  {1, 34},  {2, 39},  {3, 33},  {4, 32},  {5, 31},  {8, 40},
        {9, 29},  {10, 28}, {11, 27}, {12, 26}, {13, 25}, {14, 24}, {15, 23},
        {16, 22}, {17, 15}, {18, 14}, {19, 13}, {20, 12}, {21, 11}, {22, 10},
        {23, 9},  {24, 1},  {25, 0},  {26, 2},  {27, 3},  {28, 4},  {29, 5},
        {30, 6},  {31, 7},  {32, 8},  {33, 42}, {34, 16}, {35, 17}, {36, 18},
        {37, 19}, {38, 20}, {39, 21}, {40, 41}, {41, 30}};
    CHECK(maxCliques.front() == expClique);
  }
}

TEST_CASE("Smaller than FMCS") {
  // The normal FMCS algorithm gives 29 atoms for this pair, TwoMolMCSS gives
  // 28.  This is correct, and due to the different ways the algorithms work.
  std::vector<ROMOL_SPTR> mols{
      R"(COc1cc2c(Nc3ccc(Cl)c(Cl)c3)ncnc2cc1OCc1nc(CN2CCOCC2)cs1)"_smiles,
      R"(CCN1CCC2N=C(COc3cc4ncnc(Nc5ccc(Cl)c(Cl)c5)c4cc3OC)SC2C1)"_smiles,
  };
  MCSParameters params;
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
  auto amt =
      makeAtomMatchTable(*mols[0], *mols[1], params.AtomCompareParameters);
  auto bmt =
      makeBondMatchTable(*mols[0], *mols[1], params.BondCompareParameters);
  TwoMolMCSS(*mols[0], *mols[1], 28, amt, bmt, true, 0, maxCliques);
  printCliques(maxCliques);
}
