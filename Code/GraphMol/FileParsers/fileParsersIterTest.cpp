//
//   Copyright (C) 2002-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>

#include <vector>
#include <algorithm>
#include <execution>
#include <ranges>
#include <iostream>
#include <filesystem>

#include <GraphMol/RDKitBase.h>
#include "MolSupplier.h"
#include <GraphMol/SmilesParse/SmilesWrite.h>

static const std::string rdbase = getenv("RDBASE");

using namespace RDKit;

template <typename Supplier>
void iterTest(Supplier &reader, size_t len) {
  CHECK(reader.length() == len);
  std::vector<unsigned int> expected;
  expected.reserve(reader.length());
  for (unsigned int i = 0; i < reader.length(); ++i) {
    expected.push_back(reader[i]->getNumAtoms());
  }
  std::vector<unsigned int> actual;
  std::ranges::transform(
      reader, std::back_inserter(actual),
      [](const auto &mol) { return (size_t)mol->getNumAtoms(); });
  CHECK(actual == expected);
}

template <typename Supplier>
void forwardIterTest(Supplier &reader, size_t len) {
  std::vector<unsigned int> actual;
  std::ranges::transform(reader, std::back_inserter(actual),
                         [](const auto &mol) {
                           REQUIRE(mol);
                           return (size_t)mol->getNumAtoms();
                         });
  CHECK(actual.size() == len);
}

template <typename Supplier>
void cacheTest(Supplier &reader, size_t len) {
  reader.setCaching(true);
  CHECK(reader.length() == len);
  std::vector<std::shared_ptr<RWMol>> mols(reader.length());
  std::copy(reader.begin(), reader.end(), mols.begin());
  REQUIRE(mols.size() == reader.length());
  std::vector<size_t> expected;
  std::ranges::transform(mols, std::back_inserter(expected),
                         [](const auto &mol) {
                           REQUIRE(mol);
                           return (size_t)mol.get();
                         });
  std::vector<size_t> actual;
  std::ranges::transform(reader, std::back_inserter(actual),
                         [](const auto &mol) {
                           REQUIRE(mol);
                           return (size_t)mol.get();
                         });
  CHECK(actual == expected);
}

TEST_CASE("basic SDMolSupplier iteration") {
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  v2::FileParsers::SDMolSupplier reader(infile);
  SECTION("basics") { iterTest(reader, 16); }
  SECTION("with caching") { cacheTest(reader, 16); }
}

TEST_CASE("ForwardSDMolSupplier iteration") {
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  std::ifstream strm(infile);
  bool takeOwnership = false;
  v2::FileParsers::ForwardSDMolSupplier reader(&strm, takeOwnership);
  SECTION("basics") { forwardIterTest(reader, 16); }
  SECTION("error handling") {
    reader.next();
    CHECK_THROWS_AS(reader.begin(), ValueErrorException);
  }
  SECTION("pre-increment") {
    unsigned int i = 0;
    auto it = reader.begin();
    auto end = reader.end();
    while (it != end) {
      auto mol = *it;
      REQUIRE(mol);
      ++it;
      ++i;
    }
    CHECK(i == 16);
  }
  SECTION("post-increment") {
    unsigned int i = 0;
    auto it = reader.begin();
    auto end = reader.end();
    while (it != end) {
      auto mol = *it;
      REQUIRE(mol);
      it++;
      ++i;
    }
    CHECK(i == 16);
  }
}

TEST_CASE("cached SDMolSupplier error handling") {
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/good_bad_good_bad.sdf";
  SECTION("basics") {
    v2::FileParsers::SDMolSupplier reader(infile);
    reader.setCaching(true);
    CHECK(reader.length() == 4);
    std::vector<std::shared_ptr<RWMol>> mols(reader.length());
    std::copy(reader.begin(), reader.end(), mols.begin());
    REQUIRE(mols.size() == reader.length());
    CHECK(mols[0]);
    CHECK(!mols[1]);
    CHECK(mols[2]);
    CHECK(!mols[3]);
    std::vector<size_t> expected;
    std::ranges::transform(mols, std::back_inserter(expected),
                           [](const auto &mol) { return (size_t)mol.get(); });

    // now use the cached versions:
    std::vector<std::shared_ptr<RWMol>> nmols(reader.length());
    std::copy(reader.begin(), reader.end(), nmols.begin());
    REQUIRE(nmols.size() == reader.length());
    CHECK(nmols[0]);
    CHECK(!nmols[1]);
    CHECK(nmols[2]);
    CHECK(!nmols[3]);
    // confirm the caching
    std::vector<size_t> actual;
    std::ranges::transform(nmols, std::back_inserter(actual),
                           [](const auto &mol) { return (size_t)mol.get(); });
    CHECK(actual == expected);
  }
}

TEST_CASE("basic SmilesMolSupplier iteration") {
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/first_200.tpsa.csv";
  v2::FileParsers::SmilesMolSupplierParams params;
  params.delimiter = ',';
  params.smilesColumn = 0;
  params.nameColumn = -1;
  v2::FileParsers::SmilesMolSupplier reader(infile, params);
  SECTION("basics") { iterTest(reader, 200); }
  SECTION("with caching") { cacheTest(reader, 200); }
}

TEST_CASE("cached SmilesMolSupplier error handling") {
  v2::FileParsers::SmilesMolSupplierParams params;
  params.delimiter = ',';
  params.smilesColumn = 0;
  params.nameColumn = -1;
  std::string data = R"SMI(smiles,is_valid
CCO,1
CFC,0
CCN,1
c1cc1,0
c1cc,0
)SMI";
  SECTION("basics") {
    v2::FileParsers::SmilesMolSupplier reader;
    reader.setData(data, params);
    reader.setCaching(true);
    CHECK(reader.length() == 5);
    std::vector<std::shared_ptr<RWMol>> mols(reader.length());
    std::copy(reader.begin(), reader.end(), mols.begin());
    REQUIRE(mols.size() == reader.length());
    CHECK(mols[0]);
    CHECK(!mols[1]);
    CHECK(mols[2]);
    CHECK(!mols[3]);
    CHECK(!mols[4]);
    std::vector<size_t> expected;
    std::ranges::transform(mols, std::back_inserter(expected),
                           [](const auto &mol) { return (size_t)mol.get(); });

    // now use the cached versions:
    std::vector<std::shared_ptr<RWMol>> nmols(reader.length());
    std::copy(reader.begin(), reader.end(), nmols.begin());
    REQUIRE(nmols.size() == reader.length());
    CHECK(nmols[0]);
    CHECK(!nmols[1]);
    CHECK(nmols[2]);
    CHECK(!nmols[3]);
    CHECK(!nmols[4]);
    // confirm the caching
    std::vector<size_t> actual;
    std::ranges::transform(nmols, std::back_inserter(actual),
                           [](const auto &mol) { return (size_t)mol.get(); });
    CHECK(actual == expected);
  }
}

#ifdef RDK_BUILD_THREADSAFE_SSS
TEST_CASE("parallel reads") {
  auto *rdbase = std::getenv("RDBASE");
  REQUIRE(rdbase);
  SECTION("sdf") {
    auto path = std::filesystem::path(rdbase) /
                "Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";
    REQUIRE(std::filesystem::exists(path));
    v2::FileParsers::SDMolSupplier reader1(path.string());
    std::vector<unsigned int> nAts1(reader1.length());
    std::transform(reader1.begin(), reader1.end(), nAts1.begin(),
                   [](const auto mol) { return mol->getNumAtoms(); });
    // std::sort(nAts1.begin(), nAts1.end());
    auto start = std::chrono::high_resolution_clock::now();
    constexpr unsigned int numIters = 100;
    for (unsigned int iter = 0; iter < numIters; ++iter) {
      v2::FileParsers::SDMolSupplier reader2(path.string());
      reader2.setCaching(true);
      std::vector<unsigned int> nAts2(reader1.length());
      std::transform(std::execution::par, reader2.begin(), reader2.end(),
                     nAts2.begin(),
                     [](const auto mol) { return mol->getNumAtoms(); });
      REQUIRE(nAts1.size() == nAts2.size());
      REQUIRE(nAts1 == nAts2);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cerr << "Read of " << reader1.length() << "x" << numIters
              << " molecules took " << duration.count() << " ms" << std::endl;
  }
  SECTION("smiles") {
    auto path =
        std::filesystem::path(rdbase) / "Regress/Data/zinc.leads.500.q.smi";
    REQUIRE(std::filesystem::exists(path));
    v2::FileParsers::SmilesMolSupplierParams params;
    params.delimiter = '\t';
    params.smilesColumn = 0;
    params.nameColumn = 1;
    params.titleLine = false;

    v2::FileParsers::SmilesMolSupplier reader1(path.string(), params);
    std::vector<unsigned int> nAts1(reader1.length());
    std::transform(reader1.begin(), reader1.end(), nAts1.begin(),
                   [](const auto mol) { return mol->getNumAtoms(); });

    auto start = std::chrono::high_resolution_clock::now();
    constexpr unsigned int numIters = 100;
    for (unsigned int iter = 0; iter < numIters; ++iter) {
      v2::FileParsers::SmilesMolSupplier reader2(path.string(), params);
      reader2.setCaching(true);
      std::vector<unsigned int> nAts2(reader1.length());
      std::transform(std::execution::par, reader2.begin(), reader2.end(),
                     nAts2.begin(),
                     [](const auto mol) { return mol->getNumAtoms(); });
      REQUIRE(nAts1 == nAts2);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cerr << "Read of " << reader1.length() << "x" << numIters
              << " molecules took " << duration.count() << " ms" << std::endl;
  }
}
#endif
#if 1
TEST_CASE("benchmarking") {
  auto *rdbase = std::getenv("RDBASE");
  REQUIRE(rdbase);
  auto path =
      std::filesystem::path(rdbase) / "Regress/Data/zinc.leads.500.q.smi";
  REQUIRE(std::filesystem::exists(path));
  v2::FileParsers::SmilesMolSupplierParams params;
  params.delimiter = '\t';
  params.smilesColumn = 0;
  params.nameColumn = 1;
  params.titleLine = false;
  v2::FileParsers::SmilesMolSupplier reader(path.string(), params);
  reader.setCaching(true);

  SECTION("transform") {
    auto start = std::chrono::high_resolution_clock::now();
    // prime the cache:
    std::vector<unsigned int> nAts1;
    std::transform(reader.begin(), reader.end(), std::back_inserter(nAts1),
                   [](const auto mol) { return mol->getNumAtoms(); });

    double accum = 0.0;
    for (unsigned int iter = 0; iter < 1000; ++iter) {
      std::vector<unsigned int> nAts(reader.length());
      std::transform(std::execution::seq, reader.begin(), reader.end(),
                     nAts.begin(),
                     [](const auto mol) { return MolToSmiles(*mol).size(); });
      accum += nAts.size();
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cerr << "Base transform of " << reader.length() << " molecules took "
              << duration.count() << " ms" << std::endl;
    CHECK(accum > 0);
#if 1
    accum = 0.0;
    start = std::chrono::high_resolution_clock::now();
    for (unsigned int iter = 0; iter < 1000; ++iter) {
      std::vector<unsigned int> nAts(reader.length());
      std::transform(std::execution::par, reader.begin(), reader.end(),
                     nAts.begin(),
                     [](const auto mol) { return MolToSmiles(*mol).size(); });
      accum += nAts.size();
    }
    end = std::chrono::high_resolution_clock::now();
    duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cerr << "Parallel transform of " << reader.length()
              << " molecules took " << duration.count() << " ms" << std::endl;
    CHECK(accum > 0);
#endif
  }
}
#endif