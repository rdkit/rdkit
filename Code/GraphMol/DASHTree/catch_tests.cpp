//
//  Copyright (C) 2026 Marc Lehner and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <GraphMol/DASHTree/DASHTree.h>
#include <GraphMol/DASHTree/MatchKey.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDLog.h>

using namespace RDKit;

namespace {
std::unique_ptr<ROMol> molWithHs(const std::string &smiles) {
  std::unique_ptr<ROMol> mol(SmilesToMol(smiles));
  REQUIRE(mol);
  return std::unique_ptr<ROMol>(MolOps::addHs(*mol));
}
}  // namespace

TEST_CASE("DASH atom features", "[DASHTree]") {
  // the classes are a property of the tree, read from the table it carries;
  // the synthetic tree has four of them
  const char *rdbase = std::getenv("RDBASE");
  REQUIRE(rdbase != nullptr);
  DASH::DASHTree tree(
      std::string(rdbase) +
          "/Code/GraphMol/DASHTree/test_data/dash_test_tree.dash",
      {"result"});
  REQUIRE(tree.getNumBranches() == 4);
  auto ethanol = molWithHs("CCO");
  CHECK(tree.getAtomFeatureIndex(ethanol->getAtomWithIdx(0)) == 0);
  CHECK(tree.getAtomFeature(2).atomicNum == 8);
  // an atom the tree does not cover must say so rather than misclassify
  CHECK(tree.getAtomFeatureIndex(
            molWithHs("C[Si](C)(C)C")->getAtomWithIdx(1)) == -1);
  CHECK(tree.getAtomFeatureIndex(molWithHs("c1ccccc1")->getAtomWithIdx(0)) ==
        -1);  // aromatic CH
  CHECK_THROWS_AS(tree.getAtomFeature(4), IndexErrorException);
}

TEST_CASE("DASH packed match key", "[DASHTree]") {
  // the two ways of building a key must agree, and a component out of range
  // must flag rather than alias onto some real node's key
  CHECK(DASH::packMatchKey(27, 2, 1) ==
        DASH::withConAtom(DASH::packPartialMatchKey(27, 1), 2));
  CHECK(DASH::packMatchKey(27, 0, 17) == DASH::noMatchKey);  // DATIVE
  CHECK(DASH::packMatchKey(DASH::maxKeyAtomType + 1, 0, 1) == DASH::noMatchKey);
  CHECK(DASH::withConAtom(DASH::packPartialMatchKey(27, 1),
                          DASH::maxKeyConAtom + 1) == DASH::noMatchKey);
  CHECK(DASH::packMatchKey(27, 0, 12) != DASH::packMatchKey(27, 0, 4));
  // and a key unpacks to what it was packed from, -1 included
  const std::uint16_t key = DASH::packMatchKey(27, 2, 4);
  CHECK(DASH::keyAtomType(key) == 27);
  CHECK(DASH::keyConAtom(key) == 2);
  CHECK(DASH::keyConType(key) == 4);
  CHECK(DASH::keyConAtom(DASH::packMatchKey(27, -1, -1)) == -1);
  CHECK(DASH::keyConType(DASH::packMatchKey(27, -1, -1)) == -1);
}

TEST_CASE("DASH integration on the synthetic tree", "[DASHTree]") {
  // test_data/dash_test_tree.dash is written by
  // rdkit/Chem/DASHTree/UnitTestDASHTree.py
  // and checked in: the four atom classes of ethanol and nothing else. Its
  // values are made up, so this checks the machinery -- the descent, the
  // property fallback, the stop flag, the normalisations and the node API --
  // rather than any chemistry.
  const char *rdbase = std::getenv("RDBASE");
  REQUIRE(rdbase != nullptr);
  DASH::DASHTree tree(
      std::string(rdbase) +
          "/Code/GraphMol/DASHTree/test_data/dash_test_tree.dash",
      {"result", "std"});
  CHECK(tree.getNumBranches() == 4);
  CHECK(tree.getNumNodes() == 18);
  auto mol = molWithHs("CCO");  // C0 H3 H4 H5, C1 H6 H7, O2 H8
  const std::uint32_t M = 0, P = 1, O = 2, H = 3;  // its four atom classes
  REQUIRE(tree.getAtomFeatureIndex(mol->getAtomWithIdx(2)) ==
          static_cast<int>(O));

  // every atom descends where the tree says, hydrogens through their neighbour
  const std::vector<std::pair<unsigned int, std::vector<std::uint32_t>>> paths =
      {{0, {M, 0, 1, 2, 3}}, {1, {P, 0, 1, 3}}, {2, {O, 0, 1}},
       {3, {H, 0, 1, 4}},    {6, {H, 0, 2}},    {8, {H, 0, 3, 5}}};
  std::vector<std::uint32_t> path;
  for (const auto &[atom, expected] : paths) {
    tree.getAtomNodePath(*mol, atom, path);
    CHECK(path == expected);
  }
  std::vector<unsigned int> atoms;
  tree.getMatchedSubstructure(*mol, 0, atoms);
  CHECK(atoms == std::vector<unsigned int>{0, 1, 2, 8});

  // the O descent stops at the flagged node; told to read the attentions
  // instead, it carries on. Atom 0's deepest node has no value, so its walk
  // falls back to the node above.
  DASH::DASHParams deeper;
  deeper.attentionThreshold = 100.0;
  tree.getAtomNodePath(*mol, 2, path, deeper);
  CHECK(path == std::vector<std::uint32_t>{O, 0, 1, 3});
  CHECK(tree.getAtomProperty(*mol, 2, "result") == -0.65);
  CHECK(tree.getAtomProperty(*mol, 2, "result", deeper) == -0.70);
  CHECK(tree.getAtomProperty(*mol, 0, "result") == -0.45);

  // the raw values, the deviations and the three normalisations
  const std::vector<double> raw = {-0.45, -0.18, -0.65, 0.12, 0.12,
                                   0.12,  0.08,  0.08,  0.31};
  const std::vector<double> dev = {0.03, 0.02, 0.05, 0.02, 0.02,
                                   0.02, 0.02, 0.02, 0.02};
  std::vector<double> charges, rawOut, devOut;
  std::vector<unsigned int> depths;
  tree.getPartialCharges(*mol, charges, rawOut, devOut, depths);
  CHECK(rawOut == raw);
  CHECK(devOut == dev);
  CHECK(depths == std::vector<unsigned int>{4, 3, 2, 3, 3, 3, 2, 2, 3});
  double deficit = 0.0, devTotal = 0.0;
  for (std::size_t i = 0; i < raw.size(); ++i) {
    deficit -= raw[i];
    devTotal += dev[i];
  }
  for (std::size_t i = 0; i < raw.size(); ++i) {
    CHECK_THAT(charges[i], Catch::Matchers::WithinAbs(
                               raw[i] + deficit * dev[i] / devTotal, 1e-12));
  }
  DASH::NormalizationOptions symmetric;
  symmetric.normalization = DASH::Normalization::SYMMETRIC;
  tree.getPartialCharges(*mol, charges, symmetric);
  for (std::size_t i = 0; i < raw.size(); ++i) {
    CHECK_THAT(charges[i],
               Catch::Matchers::WithinAbs(raw[i] + deficit / 9.0, 1e-12));
  }

  // charges are the normalised assignment with the formal charge as target;
  // any other target works the same way, and a non-positive default
  // deviation is refused rather than dividing by it
  std::vector<double> values;
  tree.getPartialCharges(*mol, charges);
  tree.getNormalizedMolProperty(*mol, 0.0, values);
  CHECK(values == charges);
  tree.getNormalizedMolProperty(*mol, 1.0, values);
  double total = 0.0;
  for (const auto value : values) {
    total += value;
  }
  CHECK_THAT(total, Catch::Matchers::WithinAbs(1.0, 1e-12));
  DASH::NormalizationOptions noStd;
  noStd.defaultStdValue = 0.0;
  CHECK_THROWS_AS(tree.getNormalizedMolProperty(*mol, 0.0, values, noStd),
                  ValueErrorException);

  // the node API reads the same records the descent used
  DASH::DASHTreeNode root = tree.getRoot(O);
  REQUIRE(root.getNumChildren() == 2);
  DASH::DASHTreeNode first = root.getChild(0);
  CHECK(first.getBranch() == O);
  CHECK(first.getAtomFeatureIndex() == static_cast<int>(H));
  CHECK(first.getConAtom() == 0);
  CHECK(first.getConType() == 1);
  CHECK(first.getAttention() == 11.0f);
  CHECK(first.isStopNode());
  CHECK(!root.getChild(1).isStopNode());
  CHECK(root.getChild(1).getValue("result") == -0.99);
  CHECK(first.getChild(0).getId() == 3);
  CHECK(std::isnan(tree.getNode(M, 3).getValue("result")));
  CHECK_THROWS_AS(root.getChild(2), IndexErrorException);
  CHECK_THROWS_AS(tree.getNode(4, 0), IndexErrorException);
  CHECK_THROWS_AS(root.getValue("no_such_property"), ValueErrorException);

  // a threaded batch agrees with the single call, and a null entry is empty
  std::vector<std::vector<double>> batch;
  tree.getPartialChargesBatch({mol.get(), nullptr, mol.get()}, batch,
                              DASH::NormalizationOptions(), 0);
  REQUIRE(batch.size() == 3);
  CHECK(batch[1].empty());
  tree.getPartialCharges(*mol, charges);
  CHECK(batch[0] == charges);
  CHECK(batch[2] == charges);

  // an atom outside the tree's classes is refused
  CHECK_THROWS_AS(tree.getPartialCharges(*molWithHs("C"), charges),
                  ValueErrorException);

  // an unsanitized molecule is still assigned, with a warning
  std::unique_ptr<RWMol> unsanitized(SmilesToMol("CCO", 0, false));
  unsanitized->updatePropertyCache();
  MolOps::addHs(*unsanitized);
  {
    RDLog::CaptureLog capture{rdWarningLog};
    tree.getPartialCharges(*mol, charges);
    CHECK(capture.messages().empty());
    tree.getPartialCharges(*unsanitized, charges);
    CHECK(capture.messages().find("sanitized") != std::string::npos);
  }

  // an int32 column's missing value falls back like a NaN, here to node 2
  DASH::DASHTree allColumns(tree.getFilename());
  CHECK(std::isnan(allColumns.getNode(M, 3).getValue("size")));
  CHECK(allColumns.getAtomProperty(*mol, 0, "size") == 2.0);

  // a tree moves, and the nodes taken from it move along
  DASH::DASHTree moved(std::move(tree));
  CHECK(moved.getNumBranches() == 4);
  CHECK(root.getChild(1).getValue("result") == -0.99);
}

TEST_CASE("DASH corrupt containers", "[DASHTree]") {
  // a tree file is input like any other: a link inside it is checked before it
  // is followed
  const char *rdbase = std::getenv("RDBASE");
  REQUIRE(rdbase != nullptr);
  std::ifstream in(std::string(rdbase) +
                       "/Code/GraphMol/DASHTree/test_data/dash_test_tree.dash",
                   std::ios::binary);
  const std::string good{std::istreambuf_iterator<char>(in), {}};
  const auto headerField = [&](std::size_t offset) {
    std::uint64_t value;
    std::memcpy(&value, good.data() + offset, sizeof(value));
    return value;
  };
  const std::uint64_t offBranchRoot = headerField(40);
  const std::uint64_t offNodeRecord = headerField(48);
  const std::uint64_t offPropertyDir = headerField(72);
  const std::string path = "dash_corrupt.dash";
  const auto load = [&](const std::string &bytes) {
    std::ofstream(path, std::ios::binary) << bytes;
    return std::make_unique<DASH::DASHTree>(path);
  };
  const auto patch = [&](std::uint64_t offset, std::uint32_t value) {
    std::string bytes = good;
    std::memcpy(&bytes[offset], &value, sizeof(value));
    return bytes;
  };

  // the header: the magic, the version, the byte order and the file size
  CHECK_THROWS_AS(load(patch(0, 0)), BadFileException);
  CHECK_THROWS_AS(load(patch(8, 0)), ValueErrorException);
  CHECK_THROWS_AS(load(patch(12, 0xEFBEADDE)), ValueErrorException);
  CHECK_THROWS_AS(load(good.substr(0, good.size() - 64)), ValueErrorException);

  // a name that fills its field has no terminator, and is read without one
  std::string bytes = good;
  bytes.replace(offPropertyDir, 32, 32, 'x');
  CHECK(load(bytes)->getAvailablePropertyNames()[0] == std::string(32, 'x'));

  // the first branch must start at node 0
  CHECK_THROWS_AS(load(patch(offBranchRoot, 1)), ValueErrorException);

  // the children of the root of branch 0 (node 0), pointing back at it or far
  // past the end of the file
  auto mol = molWithHs("CCO");
  std::vector<double> charges;
  for (const std::uint32_t firstChild : {0u, 0xFFFFFF00u}) {
    auto tree = load(patch(offNodeRecord + 4, firstChild));
    CHECK_THROWS_AS(tree->getPartialCharges(*mol, charges),
                    ValueErrorException);
    CHECK_THROWS_AS(tree->getRoot(0).getChild(0), ValueErrorException);
  }
  std::filesystem::remove(path);
}

