//
//  Copyright (C) 2023 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>
#ifdef RDK_BUILD_THREADSAFE_SSS
#include <future>
#include <thread>
#endif

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

using namespace RDKit::v2;

TEST_CASE("v2 basics") {
  {
    auto mol = SmilesParse::MolFromSmiles("CCO[H]");
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 3);
  }
  {
    SmilesParse::SmilesParserParams ps;
    ps.removeHs = false;
    auto mol = SmilesParse::MolFromSmiles("CCO[H]", ps);
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 4);
  }
  {
    auto mol = SmilesParse::MolFromSmarts("[H]CC[R]");
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 4);
  }
  {
    SmilesParse::SmartsParserParams ps;
    ps.mergeHs = true;
    auto mol = SmilesParse::MolFromSmarts("[H]CC[R]", ps);
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 3);
  }
  {
    auto atm = SmilesParse::AtomFromSmiles("C");
    REQUIRE(atm);
  }
  {
    auto bnd = SmilesParse::BondFromSmiles("-");
    REQUIRE(bnd);
  }
  {
    auto atm = SmilesParse::AtomFromSmarts("[R]");
    REQUIRE(atm);
  }
  {
    auto bnd = SmilesParse::BondFromSmarts("@");
    REQUIRE(bnd);
  }
}
TEST_CASE("handling of aromatic Al in SMILES"){
  SECTION("basics") {
    auto mol = SmilesParse::MolFromSmiles("[Al+]1cccccccccc1");
    REQUIRE(mol);
    auto smi = RDKit::MolToSmiles(*mol);
    CHECK(smi.find("Al") != std::string::npos);
  }
}
