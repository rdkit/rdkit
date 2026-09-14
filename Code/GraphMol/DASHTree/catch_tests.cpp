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

#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include <GraphMol/DASHTree/AtomFeatures.h>
#include <GraphMol/DASHTree/DASHTree.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/BadFileException.h>

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
}
