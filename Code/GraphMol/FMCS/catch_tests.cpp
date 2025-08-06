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
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
  TwoMolMCSS(*mol1, *mol2, maxCliques);
  std::cout << "Number of max cliques: " << maxCliques.size() << std::endl;
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
  TwoMolMCSS(*mol1, *mol2, maxCliques);
  std::cout << "Number of max cliques: " << maxCliques.size() << std::endl;
  for (const auto &clique : maxCliques) {
    std::cout << clique.size() << " :";
    for (const auto &c : clique) {
      std::cout << " (" << c.first << ", " << c.second << ")";
    }
    std::cout << std::endl;
  }
  std::vector<ROMOL_SPTR> mols;
  mols.push_back(ROMOL_SPTR(mol1.release()));
  mols.push_back(ROMOL_SPTR(mol2.release()));
  MCSResult res;
  MCSParameters params;
  params.Verbose = false;
  // params.InitialSeed = smarts;
  auto beg = std::chrono::high_resolution_clock::now();
  res = findMCS(mols, &params);
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = duration_cast<std::chrono::microseconds>(end - beg);
  std::cout << "Elapsed Time: " << duration.count() / 2 << " microseconds"
            << std::endl;
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
}

TEST_CASE("FMCS Slow 1") {
  // Baseline time 1250-1300 ms.
  auto mol1 = v2::SmilesParse::MolFromSmiles("c1ccc2c(c1)c(ncn2)Nc3ccc(cc3)Cl");
  REQUIRE(mol1);
  auto mol2 =
      v2::SmilesParse::MolFromSmiles("c1ccc2c(c1)c(ncn2)Nc3ccc(cc3)C#N");
  REQUIRE(mol2);
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
  TwoMolMCSS(*mol1, *mol2, maxCliques);
  std::cout << "Number of max cliques: " << maxCliques.size() << std::endl;
  for (const auto &clique : maxCliques) {
    std::cout << clique.size() << " :";
    for (const auto &c : clique) {
      std::cout << " (" << c.first << ", " << c.second << ")";
    }
    std::cout << std::endl;
  }
  REQUIRE(!maxCliques.empty());
  auto smarts = makeSMARTSFromMCSS(*mol1, maxCliques.front());
  std::cout << smarts << std::endl;
#if 0
  std::vector<ROMOL_SPTR> mols;
  mols.push_back(ROMOL_SPTR(mol1.release()));
  mols.push_back(ROMOL_SPTR(mol2.release()));
  auto beg = std::chrono::high_resolution_clock::now();
  MCSResult res;
  MCSParameters params;
  params.Verbose = false;
  params.InitialSeed =
      "[#6]1:[#6]:[#6]:[#6]2:[#6](:[#6]:1):[#6](:[#7]:[#6]:[#7]:2)-[#7]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1";
  for (int i = 0; i < 1000; ++i) {
    res = findMCS(mols, &params);
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = duration_cast<std::chrono::microseconds>(end - beg);
  std::cout << "Elapsed Time: " << duration.count() / 1000 << " microseconds"
            << std::endl;
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  CHECK(
      res.SmartsString ==
      "[#6]1:[#6]:[#6]:[#6]2:[#6](:[#6]:1):[#6](:[#7]:[#6]:[#7]:2)-[#7]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1");
  CHECK(res.NumAtoms == 17);
  CHECK(res.NumBonds == 19);
#endif
}

TEST_CASE("FMCS Slow 2 - spirocycle") {
  // Baseline time 373000-37400 ms.
  auto mol1 = v2::SmilesParse::MolFromSmiles(
      R"(CN1CC2(CCN(C2=O)c2cncc3ccccc23)c2cc(Cl)ccc2C1=O)");
  REQUIRE(mol1);
  auto mol2 = v2::SmilesParse::MolFromSmiles(
      R"(CN1C(=O)c2ccc(Cl)cc2C11CCN(C1=O)c1cncc2ccccc12)");
  REQUIRE(mol2);
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
  TwoMolMCSS(*mol1, *mol2, maxCliques);
  std::cout << "Number of max cliques: " << maxCliques.size() << std::endl;
  for (const auto &clique : maxCliques) {
    std::cout << clique.size() << " :";
    for (const auto &c : clique) {
      std::cout << " (" << c.first << ", " << c.second << ")";
    }
    std::cout << std::endl;
  }
  REQUIRE(!maxCliques.empty());
  auto smarts = makeSMARTSFromMCSS(*mol1, maxCliques.front());
  std::cout << smarts << std::endl;
#if 0
  std::vector<ROMOL_SPTR> mols;
  mols.push_back(ROMOL_SPTR(mol1.release()));
  mols.push_back(ROMOL_SPTR(mol2.release()));
  auto beg = std::chrono::high_resolution_clock::now();
  MCSResult res;
  for (int i = 0; i < 10; ++i) {
    res = findMCS(mols);
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = duration_cast<std::chrono::microseconds>(end - beg);
  std::cout << "Elapsed Time: " << duration.count() / 10 << " microseconds"
            << std::endl;
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  CHECK(
      res.SmartsString ==
      "[#6]-[#7]-[#6](=[#8])-[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1-[#6]1-[#6]-[#6]-[#7](-[#6]-1=[#8])-[#6]1:[#6]:[#7]:[#6]:[#6]2:[#6]:1:[#6]:[#6]:[#6]:[#6]:2)-[#17]");
  CHECK(res.NumAtoms == 27);
  CHECK(res.NumBonds == 30);
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
#if 0
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
  TwoMolMCSS(*mols[0], *mols[1], maxCliques);
  std::cout << "Number of max cliques: " << maxCliques.size() << std::endl;
  for (const auto &clique : maxCliques) {
    std::cout << clique.size() << " :";
    for (const auto &c : clique) {
      std::cout << " (" << c.first << ", " << c.second << ")";
    }
    std::cout << std::endl;
  }
  REQUIRE(!maxCliques.empty());
  auto smarts = makeSMARTSFromMCSS(*mols[0], maxCliques.front());
  std::cout << smarts << std::endl;
#endif
  MCSParameters params;
  params.Verbose = false;
  // params.InitialSeed = smarts;
  auto beg = std::chrono::high_resolution_clock::now();
  res = findMCS(mols, &params);
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = duration_cast<std::chrono::microseconds>(end - beg);
  std::cout << "Elapsed Time: " << duration.count() / 2 << " microseconds"
            << std::endl;
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  CHECK(
      res.SmartsString ==
      "[#6]-[#8]-[#7]=[#6](-[#6](=[#8])-[#7]-[#6]1-[#6]2-[#16]-[#6]-[#6](=[#6](-[#7]-2-[#6]-1=[#8])-[#6](-[#8])=[#8])-[#6]-[#7](:[#6]):[#6]:[#6])-[#6]1:[#7]:[#6](:[#16]:[#6]:1-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#7]");
  CHECK(res.NumAtoms == 36);
  CHECK(res.NumBonds == 39);
}