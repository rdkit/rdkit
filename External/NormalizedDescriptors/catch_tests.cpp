//
//  Copyright (C) 2026 Brian Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Descriptors/Property.h>
#include <RDGeneral/Exceptions.h>
#include "NormalizedDescriptors.h"

using namespace RDKit::NormalizedDescriptors;
using Catch::Matchers::WithinAbs;

TEST_CASE("CumulativeHistogram lookup") {
  CumulativeHistogram table({0.0, 2.0, 4.0, 8.0}, {0.1, 0.5, 0.7, 0.9});
  SECTION("values at an edge take that bin") {
    CHECK(table.normalize(0.0) == 0.1);
    CHECK(table.normalize(2.0) == 0.5);
    CHECK(table.normalize(8.0) == 0.9);
  }
  SECTION("values between edges take the next bin, without interpolation") {
    CHECK(table.normalize(1.0) == 0.5);
    CHECK(table.normalize(3.0) == 0.7);
    CHECK(table(6.0) == 0.9);
  }
  SECTION("values outside the edges") {
    CHECK(table.normalize(-5.0) == 0.1);
    CHECK(table.normalize(8.5) == 1.0);
    CHECK(table.normalize(1e6) == 1.0);
  }
  SECTION("non-finite values") {
    CHECK(table.normalize(std::numeric_limits<double>::quiet_NaN()) == 0.0);
    CHECK(table.normalize(std::numeric_limits<double>::infinity()) == 1.0);
    CHECK(table.normalize(-std::numeric_limits<double>::infinity()) == 0.1);
  }
  SECTION("fractions are not clipped") {
    CumulativeHistogram bad({0.0, 1.0}, {-0.5, 1.5});
    CHECK(bad.normalize(0.0) == -0.5);
    CHECK(bad.normalize(1.0) == 1.5);
  }
  SECTION("single bin table") {
    CumulativeHistogram single({0.0}, {0.42});
    CHECK(single.normalize(-1.0) == 0.42);
    CHECK(single.normalize(0.0) == 0.42);
    CHECK(single.normalize(12.0) == 1.0);
  }
  SECTION("invalid tables") {
    CHECK_THROWS_AS(CumulativeHistogram({}, {}), ValueErrorException);
    CHECK_THROWS_AS(CumulativeHistogram({0.0, 1.0}, {0.0}),
                    ValueErrorException);
    CHECK_THROWS_AS(CumulativeHistogram({1.0, 0.0}, {0.0, 1.0}),
                    ValueErrorException);
  }
}

TEST_CASE("HistogramSet") {
  std::string text = R"TXT(# a comment
histogram foo 3
0 0.2
# comments are allowed between bins
5 0.5

10 1
histogram bar 2
-1 0.25
1 0.75
)TXT";
  HistogramSet tables;
  std::istringstream inStream(text);
  tables.loadFromStream(inStream);
  REQUIRE(tables.size() == 2);
  CHECK(tables.getNames() == std::vector<std::string>{"foo", "bar"});
  CHECK(tables.hasTable("foo"));
  CHECK(tables.normalize("foo", 2.5) == 0.5);
  CHECK(tables.normalize("bar", 0.0) == 0.75);

  SECTION("index-based access") {
    CHECK(tables.getTableIndex("foo") == 0);
    CHECK(tables.getTableIndex("bar") == 1);
    CHECK(tables.getTableIndex("baz") == -1);
    CHECK(&tables.getTable(1) == &tables.getTable("bar"));
    CHECK(tables.normalize(0, 2.5) == 0.5);
    CHECK_THROWS_AS(tables.getTable(2), IndexErrorException);
    CHECK_THROWS_AS(tables.normalize(2, 1.0), IndexErrorException);
    // replacing a table keeps its index
    tables.addTable("foo", CumulativeHistogram({0, 1}, {0.3, 0.6}));
    CHECK(tables.getTableIndex("foo") == 0);
    CHECK(tables.size() == 2);
    CHECK(tables.normalize(0, 0.5) == 0.6);
  }
  SECTION("missing descriptors normalize to 0") {
    CHECK(!tables.hasTable("baz"));
    CHECK(tables.normalize("baz", 1.0) == 0.0);
    CHECK_THROWS_AS(tables.getTable("baz"), KeyErrorException);
  }
  SECTION("round trip") {
    std::ostringstream outStream;
    tables.writeToStream(outStream);
    HistogramSet tables2;
    std::istringstream inStream2(outStream.str());
    tables2.loadFromStream(inStream2);
    REQUIRE(tables2.getNames() == tables.getNames());
    for (const auto &name : tables.getNames()) {
      CHECK(tables2.getTable(name).getEdges() ==
            tables.getTable(name).getEdges());
      CHECK(tables2.getTable(name).getFractions() ==
            tables.getTable(name).getFractions());
    }
  }
  SECTION("malformed input") {
    for (const auto &bad : {
             "foo 1\n0 0",                      // missing keyword
             "descriptor foo norm 0 1 1\n0 0",  // fitted-table header
             "histogram foo",                   // short header
             "histogram foo 0",                 // no bins
             "histogram foo 2\n0 0",            // too few bins
             "histogram foo 1\n0 0 1",          // extra field
             "histogram foo 1\na 0",            // bad number
             "histogram foo 1 7\n0 0",          // extra header field
             "histogram foo 2\n1 0\n0 1",       // unsorted edges
         }) {
      INFO(bad);
      HistogramSet badTables;
      std::istringstream badStream(bad);
      CHECK_THROWS_AS(badTables.loadFromStream(badStream), ValueErrorException);
    }
  }
}

TEST_CASE("default tables match descriptastorus") {
  const auto &tables = getDefaultTables();
  CHECK(tables.size() >= 201);
  // reference values from descriptastorus' RDKit2DHistogramNormalized, i.e.
  // bins[bisect(bins, (v,))][1] with the bins in descriptastorus' hists.py
  struct {
    const char *name;
    double value;
    double expected;
  } refs[] = {
      {"MolLogP", 2.5, 0.29897092796495756},
      {"MolLogP", -100.0, 1.000070004900343e-05},
      {"MolLogP", 100.0, 1.0},
      {"NumHDonors", 0.0, 0.5992519476363345},
      {"NumHDonors", 2.0, 0.9729981098676908},
      {"ExactMolWt", 350.1, 0.37919654375806305},
      {"TPSA", 75.3, 0.5698098866920684},
      {"fr_benzene", 1.0, 0.8641404898342884},
      {"qed", 0.6, 0.5122958607102497},
  };
  for (const auto &ref : refs) {
    INFO(ref.name << " " << ref.value);
    CHECK(tables.normalize(ref.name, ref.value) == ref.expected);
  }
}

namespace {
std::vector<std::string> splitTabs(const std::string &line) {
  std::vector<std::string> res;
  std::istringstream ss(line);
  std::string field;
  while (std::getline(ss, field, '\t')) {
    res.push_back(field);
  }
  return res;
}
}  // namespace

TEST_CASE("calcNormalizedDescriptors matches rdkit.Chem.Descriptors") {
  const auto &names = getNormalizedDescriptorNames();
  REQUIRE(names.size() == 217);
  CHECK(names.front() == "MaxAbsEStateIndex");
  CHECK(names.back() == "fr_urea");

  const char *rdbase = std::getenv("RDBASE");
  REQUIRE(rdbase);
  std::ifstream inStream(
      std::string(rdbase) +
      "/External/NormalizedDescriptors/test_data/reference_descriptors.tsv");
  REQUIRE(inStream);
  std::string line;
  REQUIRE(std::getline(inStream, line));
  auto header = splitTabs(line);
  REQUIRE(header.size() == names.size() + 1);
  for (size_t i = 0; i < names.size(); ++i) {
    CHECK(header[i + 1] == names[i]);
  }

  const auto &tables = getDefaultTables();
  NormalizedProperties props;
  REQUIRE(props.getDescriptorNames() == names);
  unsigned int nMols = 0;
  while (std::getline(inStream, line)) {
    auto fields = splitTabs(line);
    REQUIRE(fields.size() == names.size() + 1);
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(fields[0]));
    REQUIRE(mol);
    ++nMols;
    auto raw = calcDescriptorValues(*mol);
    auto normalized = calcNormalizedDescriptors(*mol);
    REQUIRE(raw.size() == names.size());
    REQUIRE(normalized.size() == names.size());
    CHECK(props.computeProperties(*mol) == normalized);
    CHECK(calcNormalizedDescriptors(*mol, tables) == normalized);
    for (size_t i = 0; i < names.size(); ++i) {
      INFO(fields[0] << " " << names[i]);
      double expected = std::stod(fields[i + 1]);
      double tol = 1e-4 * std::max(1.0, std::fabs(expected));
      CHECK_THAT(raw[i], WithinAbs(expected, tol));
      CHECK(normalized[i] == tables.normalize(names[i], raw[i]));
      CHECK(normalized[i] >= 0.0);
      CHECK(normalized[i] <= 1.0);
    }
  }
  CHECK(nMols == 10);
}

TEST_CASE("NormalizedProperties") {
  const auto &names = getNormalizedDescriptorNames();
  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol("c1ccccc1C(=O)O"));
  REQUIRE(mol);

  SECTION("descriptors come from the property registry") {
    CHECK(getPropertyName("ExactMolWt") == "exactmw");
    CHECK(getPropertyName("qed") == "qed");
    NormalizedProperties props;
    auto propNames = props.getPropertyNames();
    REQUIRE(propNames.size() == names.size());
    for (size_t i = 0; i < names.size(); ++i) {
      CHECK(propNames[i] == getPropertyName(names[i]));
    }
    // each one is registered on its own, so it can be used on its own
    auto available = RDKit::Descriptors::Properties::getAvailableProperties();
    auto raw = props.computeRawProperties(*mol);
    for (size_t i = 0; i < names.size(); ++i) {
      INFO(names[i]);
      CHECK(std::find(available.begin(), available.end(), propNames[i]) !=
            available.end());
      auto single =
          (*RDKit::Descriptors::Properties::getProperty(propNames[i]))(*mol);
      CHECK_THAT(single, WithinAbs(raw[i], 1e-12));
    }
  }
  SECTION("a subset, with other tables") {
    HistogramSet tables;
    tables.addTable("ExactMolWt", CumulativeHistogram({0, 200}, {0.25, 0.5}));
    NormalizedProperties props({"ExactMolWt", "NumHDonors"}, tables);
    CHECK(props.getDescriptorNames() ==
          std::vector<std::string>{"ExactMolWt", "NumHDonors"});
    CHECK(props.getPropertyNames() ==
          std::vector<std::string>{"exactmw", "NumHBD"});
    auto raw = props.computeRawProperties(*mol);
    CHECK_THAT(raw[0], WithinAbs(122.0368, 1e-4));
    CHECK(raw[1] == 1.0);
    // NumHDonors has no table, so it gets the failure value
    CHECK(props.computeProperties(*mol) == std::vector<double>{0.5, 0.0});
    props.computeProperties(*mol, true);
    CHECK(mol->getProp<double>("ExactMolWt") == 0.5);
  }
  SECTION("batches are normalized too") {
    std::unique_ptr<RDKit::ROMol> mol2(RDKit::SmilesToMol("CCO"));
    REQUIRE(mol2);
    NormalizedProperties props(getNormalizedDescriptorNames(),
                               getDefaultTables());
    std::vector<const RDKit::ROMol *> mols{mol.get(), nullptr, mol2.get()};
    for (int numThreads : {1, 2}) {
      auto rows = props.computePropertiesForMols(mols, numThreads);
      REQUIRE(rows.size() == 3);
      CHECK(rows[0] == props.computeProperties(*mol));
      CHECK(std::isnan(rows[1][0]));
      CHECK(rows[2] == props.computeProperties(*mol2));
    }
  }
  SECTION("unknown names") {
    CHECK_THROWS_AS(
        NormalizedProperties({"NotADescriptor"}, getDefaultTables()),
        KeyErrorException);
  }
}

TEST_CASE("Properties failures") {
  struct Failing : public RDKit::Descriptors::PropertyFunctor {
    Failing() : PropertyFunctor("AlwaysFails", "1.0.0") {}
    double operator()(const RDKit::ROMol &) const override {
      throw ValueErrorException("failed");
    }
  };
  RDKit::Descriptors::Properties::registerProperty(new Failing());
  std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol("CCO"));
  REQUIRE(mol);
  RDKit::Descriptors::Properties props({"AlwaysFails", "NumHeavyAtoms"});
  auto res = props.computeProperties(*mol);
  CHECK(std::isnan(res[0]));
  CHECK(res[1] == 3.0);

  NormalizedProperties nprops({"AlwaysFails", "NumHDonors"},
                              getDefaultTables());
  auto raw = nprops.computeRawProperties(*mol);
  CHECK(std::isnan(raw[0]));
  CHECK(nprops.computeProperties(*mol)[0] == 0.0);
}
