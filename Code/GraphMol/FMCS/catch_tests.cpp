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
#if 0
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
#endif
#if 1
  std::vector<ROMOL_SPTR> mols;
  mols.push_back(ROMOL_SPTR(mol1.release()));
  mols.push_back(ROMOL_SPTR(mol2.release()));
  MCSResult res;
  MCSParameters params;
  params.MinMCSSSize = static_cast<unsigned int>(
      0.85 * std::max(mols[0]->getNumAtoms(), mols[1]->getNumAtoms()));
  {
    params.Verbose = true;
    params.FastInitialSeed = true;
    auto beg = std::chrono::high_resolution_clock::now();
    res = findMCS(mols, &params);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(end - beg);
    std::cout << "\nMedChemica Base Test New Elapsed Time: " << duration.count()
              << " microseconds" << std::endl;
    std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms
              << " atoms, " << res.NumBonds << " bonds\n";
    CHECK(res.NumAtoms == 17);
    CHECK(res.NumBonds == 19);
  }
#endif
#if 1
  {
    params.FastInitialSeed = false;
    params.Verbose = true;
    auto beg = std::chrono::high_resolution_clock::now();
    res = findMCS(mols, &params);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(end - beg);
    std::cout << "\nMedChemica Base Test Old Elapsed Time: " << duration.count()
              << " microseconds" << std::endl;
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
  // auto smarts = makeSMARTSFromMCSS(*mol1, maxCliques.front());
  // std::cout << smarts << std::endl;

#if 1
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
    std::cout << "\nFMCS Slow 1 New Elapsed Time: " << duration.count()
              << " microseconds" << std::endl;
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
    std::cout << "\nFMCS Slow 1 Old Elapsed Time: " << duration.count()
              << " microseconds" << std::endl;
    std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms
              << " atoms, " << res.NumBonds << " bonds\n";
    CHECK(res.NumAtoms == 26);
    CHECK(res.NumBonds == 29);
  }
#endif
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
#if 1
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
    std::cout << "\nFMCS Slow 2 New Elapsed Time: " << duration.count()
              << " microseconds" << std::endl;
    std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms
              << " atoms, " << res.NumBonds << " bonds\n";
    CHECK(res.NumAtoms == 28);
    CHECK(res.NumBonds == 32);
  }
#endif
#if 1
  {
    params.Verbose = true;
    params.FastInitialSeed = false;
    auto beg = std::chrono::high_resolution_clock::now();
    res = findMCS(mols, &params);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(end - beg);
    std::cout << "\nFMCS Slow 2 Old Elapsed Time: " << duration.count()
              << " microseconds" << std::endl;
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
  auto smarts = makeSMARTSFromMCSS(*mol1, maxCliques.front());
  std::cout << smarts << std::endl;

  std::vector<ROMOL_SPTR> mols;
  mols.push_back(ROMOL_SPTR(mol1.release()));
  mols.push_back(ROMOL_SPTR(mol2.release()));
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
    std::cout << "\nFMCS Slow 3 New Elapsed Time: " << duration.count()
              << " microseconds" << std::endl;
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
    std::cout << "\nFMCS Slow 3 Old Elapsed Time: " << duration.count()
              << " microseconds" << std::endl;
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

TEST_CASE("Big Test") {
  std::vector<ROMOL_SPTR> mols{
      R"(C[S+](c1ccc(cc1)c2[nH]c(c(n2)c3ccc(cc3)F)c4ccncc4)[O-])"_smiles,
      R"(C[S+](c1ccc(cc1)c2[nH]c(c(n2)c3c(Cl)cc(cc3)F)c4ccncc4)[O-])"_smiles,
      R"(COc1ccc(cc1OC)Nc2c3cc(c(cc3ncn2)OC)OC)"_smiles,
      R"(c1cc(c(cc1C[C@@H](C(=O)O)N)O)O)"_smiles,
      R"(c1ccc2c(c1)c(nnc2Nc3ccc(cc3)Cl)Cc4ccncc4)"_smiles,
      R"(COc1cc2c(ccnc2cc1OC)Oc3cccc(c3)Br)"_smiles,
  };
  MCSParameters params;
  params.Verbose = false;
  {
    params.FastInitialSeed = true;
    std::chrono::microseconds totTime;
    for (size_t i = 0U; i < mols.size() - 1; ++i) {
      for (size_t j = i + 1; j < mols.size(); ++j) {
        params.MinMCSSSize = static_cast<unsigned int>(
            0.85 * std::max(mols[0]->getNumAtoms(), mols[1]->getNumAtoms()));
        auto beg = std::chrono::high_resolution_clock::now();
        std::cout << "\nDoing " << i << " vs " << j << " : "
                  << params.MinMCSSSize << std::endl;
        std::vector<ROMOL_SPTR> tmpMols{mols[i], mols[j]};
        auto res = findMCS(tmpMols, &params);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = duration_cast<std::chrono::microseconds>(end - beg);
        totTime += duration;
        std::cout << "New Elapsed Time: " << duration.count() << " microseconds"
                  << std::endl;
        std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms
                  << " atoms, " << res.NumBonds << " bonds\n";
      }
    }
    std::cout << "Total Time : " << totTime << std::endl;
  }
  {
    params.FastInitialSeed = false;
    std::chrono::microseconds totTime;
    for (size_t i = 0U; i < mols.size() - 1; ++i) {
      for (size_t j = i + 1; j < mols.size(); ++j) {
        params.MinMCSSSize = static_cast<unsigned int>(
            0.85 * std::max(mols[0]->getNumAtoms(), mols[1]->getNumAtoms()));
        auto beg = std::chrono::high_resolution_clock::now();
        std::cout << "\nDoing " << i << " vs " << j << " : "
                  << params.MinMCSSSize << std::endl;
        std::vector<ROMOL_SPTR> tmpMols{mols[i], mols[j]};
        auto res = findMCS(tmpMols, &params);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = duration_cast<std::chrono::microseconds>(end - beg);
        totTime += duration;
        std::cout << "Old Elapsed Time: " << duration.count() << " microseconds"
                  << std::endl;
        std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms
                  << " atoms, " << res.NumBonds << " bonds\n";
      }
    }
    std::cout << "Total Time : " << totTime << std::endl;
  }
}