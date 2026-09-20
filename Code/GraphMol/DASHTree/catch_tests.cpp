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
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <GraphMol/DASHTree/AtomFeatures.h>
#include <GraphMol/DASHTree/DASHTree.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>

using namespace RDKit;

namespace {
std::unique_ptr<ROMol> molWithHs(const std::string &smiles) {
  std::unique_ptr<ROMol> mol(SmilesToMol(smiles));
  REQUIRE(mol);
  return std::unique_ptr<ROMol>(MolOps::addHs(*mol));
}
}  // namespace

TEST_CASE("DASH atom features", "[DASHTree]") {
  // the feature table is generated from the published tree, so check it
  // against molecules rather than against itself
  auto ethanol = molWithHs("CCO");
  CHECK(DASH::atomFeatureIndex(ethanol->getAtomWithIdx(0)) ==
        DASH::atomFeatureIndex(6, 4, 0, false, 3));
  CHECK(DASH::atomFeatureIndex(ethanol->getAtomWithIdx(2)) ==
        DASH::atomFeatureIndex(8, 2, 0, false, 1));
  // an atom the trees do not cover must say so rather than misclassify
  CHECK(DASH::atomFeatureIndex(14, 4, 0, false, 0) == -1);  // silicon
  CHECK(DASH::dashBondType(molWithHs("c1ccccc1")->getBondBetweenAtoms(0, 1)) ==
        4);
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

TEST_CASE("DASH assignment", "[DASHTree]") {
  CHECK_THROWS_AS(DASH::DASHTree("no_such_file_at_all.dash"), BadFileException);

  // A DASH tree is a couple of hundred megabytes of published data, far too
  // big for the repository. Point DASH_TREE_FILE at a container built by
  // tools/dash_convert.py to run the rest; without one it warns and passes.
  const char *treePath = std::getenv("DASH_TREE_FILE");
  if (!treePath || !*treePath) {
    WARN("set DASH_TREE_FILE to a container to run the DASH assignment tests");
    return;
  }
  const char *rdbase = std::getenv("RDBASE");
  REQUIRE(rdbase != nullptr);
  DASH::DASHTree tree(treePath, {"result", "std"});

  SDMolSupplier suppl(std::string(rdbase) +
                      "/Code/GraphMol/DASHTree/test_data/dash_test_mols.sdf");
  std::vector<std::unique_ptr<ROMol>> mols;
  std::vector<const ROMol *> ptrs;
  while (!suppl.atEnd()) {
    std::unique_ptr<ROMol> mol(suppl.next());
    REQUIRE(mol);
    mols.emplace_back(MolOps::addHs(*mol));
    ptrs.push_back(mols.back().get());
  }
  REQUIRE(mols.size() == 3);

  // the charges of a molecule have to add back up to what it started with
  for (const auto *mol : ptrs) {
    std::vector<double> charges;
    tree.getPartialCharges(*mol, charges);
    REQUIRE(charges.size() == mol->getNumAtoms());
    double total = 0.0;
    for (const auto charge : charges) {
      total += charge;
    }
    CHECK_THAT(total, Catch::Matchers::WithinAbs(MolOps::getFormalCharge(*mol),
                                                 1e-12));
  }

  // and threading may not change an answer
  std::vector<double> single;
  tree.getPartialCharges(*ptrs.back(), single);
  for (const int numThreads : {1, 0}) {  // 0 means every core
    std::vector<std::vector<double>> batch;
    tree.getPartialChargesBatch(ptrs, batch, DASH::ChargeOptions(), numThreads);
    REQUIRE(batch.size() == ptrs.size());
    CHECK(batch.back() == single);
  }

  // the nodes of a match are reachable by the ids the path reports, each is a
  // child of the one before, and the matched atoms line up with them
  const ROMol &ethanol = *ptrs.front();
  std::vector<std::uint32_t> path;
  tree.getAtomNodePath(ethanol, 0, path);
  REQUIRE(path.size() >= 3);
  DASH::DASHTreeNode node = tree.getRoot(path[0]);
  CHECK(node.getId() == path[1]);
  CHECK(node.getAtomFeatureIndex() == static_cast<int>(path[0]));
  CHECK(node.getConAtom() == -1);
  CHECK(node.getConType() == -1);
  for (std::size_t i = 2; i < path.size(); ++i) {
    bool isChild = false;
    for (unsigned int c = 0; c < node.getNumChildren() && !isChild; ++c) {
      isChild = node.getChild(c).getId() == path[i];
    }
    CHECK(isChild);
    node = tree.getNode(path[0], path[i]);
    CHECK(node.getBranch() == path[0]);
    CHECK(node.getConAtom() >= 0);
    CHECK(node.getConType() >= 1);
  }
  CHECK_THROWS_AS(node.getChild(node.getNumChildren()), ValueErrorException);
  CHECK_THROWS_AS(tree.getNode(tree.numBranches(), 0), ValueErrorException);
  CHECK_THROWS_AS(node.getValue("no_such_property"), ValueErrorException);
  // the deepest node on the path carrying a value is what the atom gets
  double deepest = std::numeric_limits<double>::quiet_NaN();
  for (std::size_t i = path.size(); i-- > 1 && std::isnan(deepest);) {
    deepest = tree.getNode(path[0], path[i]).getValue("result");
  }
  CHECK(deepest == tree.getAtomProperty(ethanol, 0, "result"));
  std::vector<unsigned int> atoms;
  tree.getMatchedSubstructure(ethanol, 0, atoms);
  REQUIRE(atoms.size() == path.size() - 1);
  CHECK(atoms.front() == 0);
}
