//
//  Copyright (C) 2023 David Cosgrove
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

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <catch2/catch_all.hpp>

#include <GraphMol/RascalMCES/RascalMCES.h>
#include <GraphMol/RascalMCES/RascalClusterOptions.h>
#include <GraphMol/RascalMCES/RascalResult.h>

TEST_CASE("Small test", "[basics]") {
  std::string fName = getenv("RDBASE");
  fName += "/Contrib/Fastcluster/cdk2.smi";
  RDKit::SmilesMolSupplier suppl(fName, "\t", 1, 0, false);
  std::vector<std::shared_ptr<RDKit::ROMol>> mols;
  while (!suppl.atEnd()) {
    std::shared_ptr<RDKit::ROMol> mol(suppl.next());
    if (!mol) {
      continue;
    }
    mols.push_back(mol);
  }
  RDKit::RascalMCES::RascalClusterOptions clusOpts;
  auto clusters = RDKit::RascalMCES::rascalCluster(mols, clusOpts);
  REQUIRE(clusters.size() == 8);
  std::vector<size_t> expSizes{7, 7, 6, 2, 2, 2, 2, 20};
  for (size_t i = 0; i < 8; ++i) {
    REQUIRE(clusters[i].size() == expSizes[i]);
  }
}

TEST_CASE("BLSets subset", "[basics]") {
  std::string fName = getenv("RDBASE");
  fName += "/Code/GraphMol/RascalMCES/data/test_cluster1.smi";
  RDKit::SmilesMolSupplier suppl(fName, "\t", 1, 0, false);
  std::vector<std::shared_ptr<RDKit::ROMol>> mols;
  while (!suppl.atEnd()) {
    std::shared_ptr<RDKit::ROMol> mol(suppl.next());
    if (!mol) {
      continue;
    }
    mols.push_back(mol);
  }
  auto clusters = RDKit::RascalMCES::rascalCluster(mols);
  REQUIRE(clusters.size() == 12);
  std::vector<size_t> expSizes{8, 4, 4, 3, 3, 3, 2, 2, 2, 2, 2, 21};
  for (size_t i = 0; i < 12; ++i) {
    REQUIRE(clusters[i].size() == expSizes[i]);
  }
}

TEST_CASE("ChEMBL 1907596") {
  std::string fName = getenv("RDBASE");
  fName += "/Code/GraphMol/RascalMCES/data/chembl_1907596.smi";
  std::cout << fName << std::endl;
  RDKit::SmilesMolSupplier suppl(fName, "\t", 1, 0, false);
  std::vector<std::shared_ptr<RDKit::ROMol>> mols;
  while (!suppl.atEnd()) {
    std::shared_ptr<RDKit::ROMol> mol(suppl.next());
    if (!mol) {
      continue;
    }
    mols.push_back(mol);
  }
  RDKit::RascalMCES::RascalClusterOptions clusOpts;
  clusOpts.similarityCutoff = 0.7;
  auto clusters = RDKit::RascalMCES::rascalCluster(mols, clusOpts);
  REQUIRE(clusters.size() == 21);
  std::vector<size_t> expSizes{343, 71, 64, 33, 23, 11, 10, 6, 6, 5, 5,
                               4,   3,  3,  3,  3,  3,  2,  2, 2, 14};
  for (size_t i = 0; i < 21; ++i) {
    REQUIRE(clusters[i].size() == expSizes[i]);
  }
}

TEST_CASE("Small Butina test", "[basics]") {
  std::string fName = getenv("RDBASE");
  fName += "/Contrib/Fastcluster/cdk2.smi";
  RDKit::SmilesMolSupplier suppl(fName, "\t", 1, 0, false);
  std::vector<std::shared_ptr<RDKit::ROMol>> mols;
  while (!suppl.atEnd()) {
    std::shared_ptr<RDKit::ROMol> mol(suppl.next());
    if (!mol) {
      continue;
    }
    mols.push_back(mol);
  }
  RDKit::RascalMCES::RascalClusterOptions clusOpts;
  auto clusters = RDKit::RascalMCES::rascalButinaCluster(mols, clusOpts);
  unsigned int numMols = 0;
  for (const auto &cl : clusters) {
    numMols += cl.size();
  }
  REQUIRE(numMols == mols.size());
  REQUIRE(clusters.size() == 29);
  std::vector<size_t> expSizes{6, 6, 6, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  for (size_t i = 0; i < 29; ++i) {
    REQUIRE(clusters[i].size() == expSizes[i]);
  }
}

TEST_CASE("Small test, smaller number of threads", "[basics]") {
  // I'm not sure how to test whether this has had the desired effect,
  // but at least we'll know that it runs ok.
  std::string fName = getenv("RDBASE");
  fName += "/Contrib/Fastcluster/cdk2.smi";
  RDKit::SmilesMolSupplier suppl(fName, "\t", 1, 0, false);
  std::vector<std::shared_ptr<RDKit::ROMol>> mols;
  while (!suppl.atEnd()) {
    std::shared_ptr<RDKit::ROMol> mol(suppl.next());
    if (!mol) {
      continue;
    }
    mols.push_back(mol);
  }
  {
    RDKit::RascalMCES::RascalClusterOptions clusOpts;
    clusOpts.numThreads = 2;
    auto clusters = RDKit::RascalMCES::rascalCluster(mols, clusOpts);
    REQUIRE(clusters.size() == 8);
    std::vector<size_t> expSizes{7, 7, 6, 2, 2, 2, 2, 20};
    for (size_t i = 0; i < 8; ++i) {
      REQUIRE(clusters[i].size() == expSizes[i]);
    }
  }
  {
    RDKit::RascalMCES::RascalClusterOptions clusOpts;
    clusOpts.numThreads = -2;
    auto clusters = RDKit::RascalMCES::rascalCluster(mols, clusOpts);
    REQUIRE(clusters.size() == 8);
    std::vector<size_t> expSizes{7, 7, 6, 2, 2, 2, 2, 20};
    for (size_t i = 0; i < 8; ++i) {
      REQUIRE(clusters[i].size() == expSizes[i]);
    }
  }
}
