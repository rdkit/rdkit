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
#include <ranges>

#include <GraphMol/RDKitBase.h>
#include "MolSupplier.h"

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
void cacheTest(Supplier &reader, size_t len) {
  reader.setCaching(true);
  CHECK(reader.length() == len);
  std::vector<std::shared_ptr<RWMol>> mols(reader.length());
  std::copy(reader.begin(), reader.end(), mols.begin());
  REQUIRE(mols.size() == reader.length());
  std::vector<size_t> expected;
  std::ranges::transform(mols, std::back_inserter(expected),
                         [](const auto &mol) { return (size_t)mol.get(); });
  std::vector<size_t> actual;
  std::ranges::transform(reader, std::back_inserter(actual),
                         [](const auto &mol) { return (size_t)mol.get(); });
  CHECK(actual == expected);
}

TEST_CASE("basic SDMolSupplier iteration") {
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  v2::FileParsers::SDMolSupplier reader(infile);
  SECTION("basics") { iterTest(reader, 16); }
  SECTION("with caching") { cacheTest(reader, 16); }
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
