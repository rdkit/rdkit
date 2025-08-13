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
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <catch2/catch_all.hpp>

using namespace RDKit;

TEST_CASE("BK Test") {
  auto mol1 = v2::SmilesParse::MolFromSmiles("CCC");
  REQUIRE(mol1);
  auto mol2 = v2::SmilesParse::MolFromSmiles("CC(C)C");
  REQUIRE(mol2);
  // auto mol1 = v2::SmilesParse::MolFromSmiles("NC1CC(Br)C1");
  // REQUIRE(mol1);
  // auto mol2 = v2::SmilesParse::MolFromSmiles("NC1CC(Cl)C1");
  // REQUIRE(mol2);
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
  unsigned int minMCSSSize(0.85 *
                           std::max(mol1->getNumAtoms(), mol2->getNumAtoms()));
  TwoMolMCSS(*mol1, *mol2, minMCSSSize, maxCliques);
  std::cout << "Number of max cliques: " << maxCliques.size() << std::endl;
  for (auto &clique : maxCliques) {
    std::sort(clique.begin(), clique.end());
  }
  std::sort(maxCliques.begin(), maxCliques.end());
  for (const auto &clique : maxCliques) {
    std::cout << clique.size() << " :";
    for (const auto &c : clique) {
      std::cout << " (" << c.first << ", " << c.second << ")";
    }
    std::cout << std::endl;
  }
  REQUIRE(!maxCliques.empty());
  CHECK(maxCliques.size() == 6);
  CHECK(maxCliques.front().size() == 3);
  auto smarts = makeSMARTSFromMCSS(*mol1, maxCliques.front());
  std::cout << smarts << std::endl;
  CHECK(smarts == "[#6]-[#6]-[#6]");
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
  TwoMolMCSS(*mol1, *mol2, minMCSSSize, maxCliques);
  std::cout << "Number of max cliques: " << maxCliques.size() << std::endl;
  std::sort(maxCliques.begin(), maxCliques.end());
  for (const auto &clique : maxCliques) {
    std::cout << clique.size() << " :";
    for (const auto &c : clique) {
      std::cout << " (" << c.first << ", " << c.second << ")";
    }
    std::cout << std::endl;
  }
  CHECK(!maxCliques.empty());
  CHECK(maxCliques.front().size() == 17);
#if 0
  std::vector<ROMOL_SPTR> mols;
  mols.push_back(ROMOL_SPTR(mol1.release()));
  mols.push_back(ROMOL_SPTR(mol2.release()));
  MCSResult res;
  MCSParameters params;
  params.Verbose = false;
  params.MinMCSSSize = static_cast<unsigned int>(
      0.85 * std::max(mols[0]->getNumAtoms(), mols[1]->getNumAtoms()));
  {
    params.FastInitialSeed = true;
    auto beg = std::chrono::high_resolution_clock::now();
    res = findMCS(mols, &params);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(end - beg);
    std::cout << "\nNew Elapsed Time: " << duration.count() << " microseconds"
              << std::endl;
    std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms
              << " atoms, " << res.NumBonds << " bonds\n";
    CHECK(res.NumAtoms == 17);
    CHECK(res.NumBonds == 19);
  }
  {
    params.FastInitialSeed = false;
    params.Verbose = true;
    auto beg = std::chrono::high_resolution_clock::now();
    res = findMCS(mols, &params);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(end - beg);
    std::cout << "\nOld Elapsed Time: " << duration.count() << " microseconds"
              << std::endl;
    std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms
              << " atoms, " << res.NumBonds << " bonds\n";
    CHECK(res.NumAtoms == 17);
    CHECK(res.NumBonds == 19);
  }
#endif
}

TEST_CASE("FMCS Slow 1") {
  // Spirocycle/open chain
  // Baseline time 1250-1300 ms.
  auto mol1 = v2::SmilesParse::MolFromSmiles(
      "CN1C[C@]2(CCN(C2=O)c2cncc3ccccc23)c2cc(Cl)ccc2C1=O");
  REQUIRE(mol1);
  auto mol2 = v2::SmilesParse::MolFromSmiles(
      "CN1C[C@@H](C(=O)Nc2cncc3ccccc23)c2cc(Cl)ccc2C1=O");
  REQUIRE(mol2);
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
  unsigned int minMCSSSize(0.85 *
                           std::max(mol1->getNumAtoms(), mol2->getNumAtoms()));
  TwoMolMCSS(*mol1, *mol2, minMCSSSize, maxCliques);
  std::cout << "Number of max cliques: " << maxCliques.size() << std::endl;
  for (const auto &clique : maxCliques) {
    std::cout << clique.size() << " :";
    for (const auto &c : clique) {
      std::cout << " (" << c.first << ", " << c.second << ")";
    }
    std::cout << std::endl;
  }
  REQUIRE(!maxCliques.empty());
  CHECK(maxCliques.front().size() == 26);
  auto smarts = makeSMARTSFromMCSS(*mol1, maxCliques.front());
  std::cout << smarts << std::endl;
  std::vector<ROMOL_SPTR> mols;
  mols.push_back(ROMOL_SPTR(mol1.release()));
  mols.push_back(ROMOL_SPTR(mol2.release()));
  MCSResult res;
  MCSParameters params;
  {
    params.Verbose = false;
    params.FastInitialSeed = true;
    params.MinMCSSSize = static_cast<unsigned int>(
        0.85 * std::max(mols[0]->getNumAtoms(), mols[1]->getNumAtoms()));
    auto beg = std::chrono::high_resolution_clock::now();
    res = findMCS(mols, &params);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(end - beg);
    std::cout << "\nNew Elapsed Time: " << duration.count() << " microseconds"
              << std::endl;
    std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms
              << " atoms, " << res.NumBonds << " bonds\n";
    CHECK(res.NumAtoms == 26);
    CHECK(res.NumBonds == 29);
  }
  {
    params.Verbose = true;
    params.FastInitialSeed = false;
    auto beg = std::chrono::high_resolution_clock::now();
    res = findMCS(mols, &params);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(end - beg);
    std::cout << "\nOld Elapsed Time: " << duration.count() << " microseconds"
              << std::endl;
    std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms
              << " atoms, " << res.NumBonds << " bonds\n";
    CHECK(res.NumAtoms == 26);
    CHECK(res.NumBonds == 29);
  }
}

TEST_CASE("FMCS Slow 2") {
  // spirocycle/spirocycle
  // Baseline time 373000-37400 ms.
  auto mol1 = v2::SmilesParse::MolFromSmiles(
      R"(CN1C[C@]2(CCN(C2=O)c2cncc3ccccc23)c2cc(Cl)ccc2C1=O)");
  REQUIRE(mol1);
  auto mol2 = v2::SmilesParse::MolFromSmiles(
      R"(CCN1C[C@]2(CCN(C2=O)c2cncc3ccccc23)c2cc(Cl)ccc2C1=O)");
  REQUIRE(mol2);
#if 1
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
  unsigned int minMCSSSize(0.85 *
                           std::max(mol1->getNumAtoms(), mol2->getNumAtoms()));
  TwoMolMCSS(*mol1, *mol2, minMCSSSize, maxCliques);
  std::cout << "Number of max cliques: " << maxCliques.size() << std::endl;
  for (const auto &clique : maxCliques) {
    std::cout << clique.size() << " :";
    for (const auto &c : clique) {
      std::cout << " (" << c.first << ", " << c.second << ")";
    }
    std::cout << std::endl;
  }
  REQUIRE(!maxCliques.empty());
  CHECK(maxCliques.front().size() == 28);
  auto smarts = makeSMARTSFromMCSS(*mol1, maxCliques.front());
  std::cout << smarts << std::endl;
#endif
#if 0
  std::vector<ROMOL_SPTR> mols;
  mols.push_back(ROMOL_SPTR(mol1.release()));
  mols.push_back(ROMOL_SPTR(mol2.release()));
  MCSResult res;
  MCSParameters params;
  {
    params.Verbose = false;
    params.FastInitialSeed = true;
    params.MinMCSSSize = static_cast<unsigned int>(
        0.85 * std::max(mols[0]->getNumAtoms(), mols[1]->getNumAtoms()));
    auto beg = std::chrono::high_resolution_clock::now();
    res = findMCS(mols, &params);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(end - beg);
    std::cout << "\nNew Elapsed Time: " << duration.count() << " microseconds"
              << std::endl;
    std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms
              << " atoms, " << res.NumBonds << " bonds\n";
    CHECK(res.NumAtoms == 28);
    CHECK(res.NumBonds == 32);
  }
#endif
#if 0
  {
    params.Verbose = true;
    params.FastInitialSeed = false;
    auto beg = std::chrono::high_resolution_clock::now();
    res = findMCS(mols, &params);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(end - beg);
    std::cout << "\nOld Elapsed Time: " << duration.count() << " microseconds"
              << std::endl;
    std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms
              << " atoms, " << res.NumBonds << " bonds\n";
    CHECK(res.NumAtoms == 28);
    CHECK(res.NumBonds == 32);
  }
#endif
}

TEST_CASE("FMCS Slow 3") {
  // Baseline time 130000000 ms.
  auto mol1 = v2::SmilesParse::MolFromSmiles(
      R"(CO\N=C(/C(=O)N[C@H]1[C@H]2SCC(C[n+]3cscc3C)=C(N2C1=O)C(O)=O)c1nc(N)sc1-c1ccccc1)");
  REQUIRE(mol1);
  auto mol2 = v2::SmilesParse::MolFromSmiles(
      R"(CO\N=C(/C(=O)N[C@H]1[C@H]2SCC(C[n+]3ccc4sccc4c3)=C(N2C1=O)C(O)=O)c1nc(N)sc1-c1ccccc1)");
  REQUIRE(mol2);
  std::vector<ROMOL_SPTR> mols;
  mols.push_back(ROMOL_SPTR(mol1.release()));
  mols.push_back(ROMOL_SPTR(mol2.release()));
  MCSResult res;
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
  unsigned int minMCSSSize(0.85 *
                           std::max(mol1->getNumAtoms(), mol2->getNumAtoms()));
  TwoMolMCSS(*mol1, *mol2, minMCSSSize, maxCliques);
  std::cout << "Number of max cliques: " << maxCliques.size() << std::endl;
  for (const auto &clique : maxCliques) {
    std::cout << clique.size() << " :";
    for (const auto &c : clique) {
      std::cout << " (" << c.first << ", " << c.second << ")";
    }
    std::cout << std::endl;
  }
  REQUIRE(!maxCliques.empty());
  CHECK(maxCliques.front().size() == 36);
  auto smarts = makeSMARTSFromMCSS(*mols[0], maxCliques.front());
  std::cout << smarts << std::endl;
  MCSParameters params;
  {
    params.Verbose = false;
    params.FastInitialSeed = true;
    params.MinMCSSSize = static_cast<unsigned int>(
        0.85 * std::max(mols[0]->getNumAtoms(), mols[1]->getNumAtoms()));
    auto beg = std::chrono::high_resolution_clock::now();
    res = findMCS(mols, &params);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(end - beg);
    std::cout << "\nNew Elapsed Time: " << duration.count() << " microseconds"
              << std::endl;
    std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms
              << " atoms, " << res.NumBonds << " bonds\n";
    CHECK(res.NumAtoms == 36);
    CHECK(res.NumBonds == 39);
  }
  {
    params.Verbose = true;
    params.FastInitialSeed = false;
    auto beg = std::chrono::high_resolution_clock::now();
    res = findMCS(mols, &params);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(end - beg);
    std::cout << "\nOld Elapsed Time: " << duration.count() << " microseconds"
              << std::endl;
    std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms
              << " atoms, " << res.NumBonds << " bonds\n";
    CHECK(res.NumAtoms == 36);
    CHECK(res.NumBonds == 39);
  }
}

TEST_CASE("Split by Rascal") {
  auto mol1 = v2::SmilesParse::MolFromSmiles("C=CCCC=CCCC1CNCCC1");
  REQUIRE(mol1);
  auto mol2 = v2::SmilesParse::MolFromSmiles("C=CCC=CCCCC1CNCCC1");
  REQUIRE(mol2);
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
  unsigned int minMCSSSize(0.85 *
                           std::max(mol1->getNumAtoms(), mol2->getNumAtoms()));
  TwoMolMCSS(*mol1, *mol2, minMCSSSize, maxCliques);
  std::cout << "Number of max cliques: " << maxCliques.size() << std::endl;
  for (const auto &clique : maxCliques) {
    std::cout << clique.size() << " :";
    for (const auto &c : clique) {
      std::cout << " (" << c.first << ", " << c.second << ")";
    }
    std::cout << std::endl;
  }
  REQUIRE(!maxCliques.empty());
  CHECK(maxCliques.size() == 1);
  CHECK(maxCliques.front().size() == 10);
  auto smarts = makeSMARTSFromMCSS(*mol1, maxCliques.front());
  std::cout << smarts << std::endl;
}