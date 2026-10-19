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

TEST_CASE("DASH packed match key", "[DASHTree]") {
  // a component out of range gives noMatchKey, not the key of a real node
  CHECK(DASH::packPartialMatchKey(27, 17) == DASH::noMatchKey);  // DATIVE
  CHECK(DASH::packPartialMatchKey(DASH::maxKeyAtomType + 1, 1) ==
        DASH::noMatchKey);
  CHECK(DASH::withConAtom(DASH::packPartialMatchKey(27, 1),
                          DASH::maxKeyConAtom + 1) == DASH::noMatchKey);
  // and a key unpacks to what it was packed from. A component it does not
  // have unpacks to -1.
  const std::uint16_t key =
      DASH::withConAtom(DASH::packPartialMatchKey(27, 4), 2);
  CHECK(DASH::keyAtomType(key) == 27);
  CHECK(DASH::keyConAtom(key) == 2);
  CHECK(DASH::keyConType(key) == 4);
  CHECK(DASH::keyConAtom(DASH::packPartialMatchKey(27, -1)) == -1);
  CHECK(DASH::keyConType(DASH::packPartialMatchKey(27, -1)) == -1);
}

TEST_CASE("DASH assignment", "[DASHTree]") {
  CHECK_THROWS_AS(DASH::DASHTree("no_such_file_at_all.dash"), BadFileException);

  const char *rdbase = std::getenv("RDBASE");
  REQUIRE(rdbase != nullptr);
  DASH::DASHTree tree(std::string(rdbase) + "/Data/DASHTree/initial_test.dash",
                      {"result", "std"});

  std::vector<std::unique_ptr<ROMol>> mols;
  std::vector<const ROMol *> ptrs;
  for (const auto *smiles : {"CCO", "CC(=O)[O-]", "CC(=O)Nc1ccc(O)cc1"}) {
    mols.push_back(molWithHs(smiles));
    ptrs.push_back(mols.back().get());
  }

  // the charges of a molecule must sum to its formal charge
  for (const auto *mol : ptrs) {
    const std::vector<double> charges = tree.getPartialCharges(*mol);
    REQUIRE(charges.size() == mol->getNumAtoms());
    double total = 0.0;
    for (const auto charge : charges) {
      total += charge;
    }
    CHECK_THAT(total, Catch::Matchers::WithinAbs(MolOps::getFormalCharge(*mol),
                                                 1e-12));
  }

  // and threading must not change an answer
  const std::vector<double> single = tree.getPartialCharges(*ptrs[1]);
  for (const int numThreads : {1, 0}) {  // 0 means every core
    const auto batch = tree.getPartialChargesBatch(
        ptrs, DASH::NormalizationOptions(), numThreads);
    REQUIRE(batch.size() == ptrs.size());
    CHECK(batch[1] == single);
  }
}

TEST_CASE("DASH integration on the synthetic tree", "[DASHTree]") {
  // rdkit/Chem/DASHTree/UnitTestDASHTree.py writes the checked-in
  // test_data/dash_test_tree.dash. It holds the four atom classes of ethanol
  // and nothing else. Its values are invented, so this case checks only the
  // machinery: the descent, the property fallback, the stop flag, the
  // normalizations and the node API.
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
  CHECK(tree.getAtomFeature(O).atomicNum == 8);
  // the tree gives -1 for an atom that it does not cover, not a wrong class
  CHECK(tree.getAtomFeatureIndex(
            molWithHs("C[Si](C)(C)C")->getAtomWithIdx(1)) == -1);

  // every atom descends where the tree says, hydrogens through their neighbour
  const std::vector<std::pair<unsigned int, std::vector<std::uint32_t>>> paths =
      {{0, {M, 0, 1, 2, 3}}, {1, {P, 0, 1, 3}}, {2, {O, 0, 1}},
       {3, {H, 0, 1, 4}},    {6, {H, 0, 2}},    {8, {H, 0, 3, 5}}};
  for (const auto &[atom, expected] : paths) {
    CHECK(tree.getAtomNodePath(*mol, atom) == expected);
  }
  CHECK(tree.getMatchedSubstructure(*mol, 0) ==
        std::vector<unsigned int>{0, 1, 2, 8});

  // the O descent stops at the flagged node. With a threshold of 100 it reads
  // the attentions instead and goes one node deeper. Atom 0's deepest node has
  // no value, so atom 0 gets the value of the node above.
  DASH::DASHParams deeper;
  deeper.attentionThreshold = 100.0;
  CHECK(tree.getAtomNodePath(*mol, 2, deeper) ==
        std::vector<std::uint32_t>{O, 0, 1, 3});
  CHECK(tree.getAtomProperty(*mol, 2, "result") == -0.65);
  CHECK(tree.getAtomProperty(*mol, 2, "result", deeper) == -0.70);
  CHECK(tree.getAtomProperty(*mol, 0, "result") == -0.45);

  // the raw values, the deviations and the three normalizations
  const std::vector<double> raw = {-0.45, -0.18, -0.65, 0.12, 0.12,
                                   0.12,  0.08,  0.08,  0.31};
  const std::vector<double> dev = {0.03, 0.02, 0.05, 0.02, 0.02,
                                   0.02, 0.02, 0.02, 0.02};
  std::vector<double> charges, rawOut, devOut;
  std::vector<unsigned int> depths;
  std::tie(charges, rawOut, devOut, depths) =
      tree.getPartialChargesDetails(*mol);
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
  charges = tree.getPartialCharges(*mol, symmetric);
  for (std::size_t i = 0; i < raw.size(); ++i) {
    CHECK_THAT(charges[i],
               Catch::Matchers::WithinAbs(raw[i] + deficit / 9.0, 1e-12));
  }

  // charges are the normalized assignment with the formal charge as target.
  // Any other target works the same way. The tree refuses a non-positive
  // default deviation, because the normalization divides by it.
  charges = tree.getPartialCharges(*mol);
  std::vector<double> values = tree.getNormalizedMolProperty(*mol, 0.0);
  CHECK(values == charges);
  values = tree.getNormalizedMolProperty(*mol, 1.0);
  double total = 0.0;
  for (const auto value : values) {
    total += value;
  }
  CHECK_THAT(total, Catch::Matchers::WithinAbs(1.0, 1e-12));
  DASH::NormalizationOptions noStd;
  noStd.defaultStdValue = 0.0;
  CHECK_THROWS_AS(tree.getNormalizedMolProperty(*mol, 0.0, noStd),
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
  CHECK_THROWS_AS(tree.getAtomFeature(4), IndexErrorException);
  CHECK_THROWS_AS(root.getValue("no_such_property"), ValueErrorException);

  // a threaded batch agrees with the single call, and a null entry is empty
  const auto batch = tree.getPartialChargesBatch(
      {mol.get(), nullptr, mol.get()}, DASH::NormalizationOptions(), 0);
  REQUIRE(batch.size() == 3);
  CHECK(batch[1].empty());
  charges = tree.getPartialCharges(*mol);
  CHECK(batch[0] == charges);
  CHECK(batch[2] == charges);

  // the tree refuses an atom outside its classes
  CHECK_THROWS_AS(tree.getPartialCharges(*molWithHs("C")), ValueErrorException);

  // the tree refuses an unsanitized molecule until its property cache is
  // updated and its hydrogens are added. It then assigns values to the
  // molecule, with a warning.
  std::unique_ptr<RWMol> unsanitized(SmilesToMol("CCO", 0, false));
  CHECK_THROWS_AS(tree.getPartialCharges(*unsanitized), ValueErrorException);
  unsanitized->updatePropertyCache();
  CHECK_THROWS_AS(tree.getPartialCharges(*unsanitized), ValueErrorException);
  MolOps::addHs(*unsanitized);
  {
    RDLog::CaptureLog capture{rdWarningLog};
    tree.getPartialCharges(*mol);
    CHECK(capture.messages().empty());
    tree.getPartialCharges(*unsanitized);
    CHECK(capture.messages().find("sanitized") != std::string::npos);
  }

  // an int32 column's missing value counts as a NaN, so atom 0 gets the value
  // of node 2
  DASH::DASHTree allColumns(tree.getFilename());
  CHECK(std::isnan(allColumns.getNode(M, 3).getValue("size")));
  CHECK(allColumns.getAtomProperty(*mol, 0, "size") == 2.0);

  // a tree moves, and the nodes taken from it stay valid
  DASH::DASHTree moved(std::move(tree));
  CHECK(moved.getNumBranches() == 4);
  CHECK(root.getChild(1).getValue("result") == -0.99);
}

TEST_CASE("DASH corrupt containers", "[DASHTree]") {
  // a tree file is input like any other. The tree checks a link in the file
  // before it follows the link.
  const char *rdbase = std::getenv("RDBASE");
  REQUIRE(rdbase != nullptr);
  std::ifstream in(std::string(rdbase) +
                       "/Code/GraphMol/DASHTree/test_data/dash_test_tree.dash",
                   std::ios::binary);
  REQUIRE(in);
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
  // an array past the end of the file, and an atom class the lookup cannot hold
  CHECK_THROWS_AS(load(patch(48, static_cast<std::uint32_t>(good.size()))),
                  ValueErrorException);
  CHECK_THROWS_AS(load(patch(headerField(96), 64)),
                  ValueErrorException);  // feature row 0 gets Z 64

  // a name that fills its field has no terminator, and is read without one
  std::string bytes = good;
  bytes.replace(offPropertyDir, 32, 32, 'x');
  CHECK(load(bytes)->getAvailablePropertyNames()[0] == std::string(32, 'x'));

  // the first branch must start at node 0
  CHECK_THROWS_AS(load(patch(offBranchRoot, 1)), ValueErrorException);

  // a patched link makes the children of the root of branch 0 (node 0) point
  // back at it or far past the end of the file
  auto mol = molWithHs("CCO");
  for (const std::uint32_t firstChild : {0u, 0xFFFFFF00u}) {
    auto tree = load(patch(offNodeRecord + 4, firstChild));
    CHECK_THROWS_AS(tree->getPartialCharges(*mol), ValueErrorException);
    CHECK_THROWS_AS(tree->getRoot(0).getChild(0), ValueErrorException);
  }
  std::filesystem::remove(path);
}

TEST_CASE("DASH atoms that fit the same node", "[DASHTree]") {
  // This case uses the initial_test tree in Data/DASHTree. Pyridine's meta
  // carbons look alike at every level. The two CH2 next to 5-aminoisoindoline's
  // ring NH look alike for five levels. The legacy rule takes the first in atom
  // order.
  const char *rdbase = std::getenv("RDBASE");
  REQUIRE(rdbase != nullptr);
  DASH::DASHTree tree(std::string(rdbase) + "/Data/DASHTree/initial_test.dash",
                      {"result", "std"});
  DASH::NormalizationOptions legacy;
  legacy.params.legacyTieBreaking = true;

  auto pyridine = molWithHs("c1ccncc1");
  std::vector<double> charges = tree.getPartialCharges(*pyridine);
  CHECK(charges[1] == charges[5]);
  charges = tree.getPartialCharges(*pyridine, legacy);
  CHECK(charges[1] != charges[5]);

  // the ring NH is atom 6 of one SMILES and atom 8 of the other
  auto first = molWithHs("Nc1ccc2CNCc2c1");
  auto second = molWithHs("Nc1ccc2c(c1)CNC2");
  CHECK(tree.getAtomNodePath(*first, 6) == tree.getAtomNodePath(*second, 8));
  CHECK(tree.getAtomNodePath(*first, 6, legacy.params) !=
        tree.getAtomNodePath(*second, 8, legacy.params));
}
