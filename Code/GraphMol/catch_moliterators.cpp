//
//  Copyright (C) 2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <algorithm>

using namespace RDKit;

TEST_CASE("mol.atoms()") {
  const auto m = "CC(C)CO"_smiles;
  REQUIRE(m);
  unsigned int ccount = 0;
  for (const auto atom : m->atoms()) {
    if (atom->getAtomicNum() == 6) {
      ++ccount;
    }
  }
  CHECK(ccount == 4);
  auto atoms = m->atoms();
  auto hasCarbon = std::any_of(atoms.begin(), atoms.end(), [](const auto atom) {
    return atom->getAtomicNum() == 6;
  });
  CHECK(hasCarbon);
  ccount = std::count_if(atoms.begin(), atoms.end(), [](const auto atom) {
    return atom->getAtomicNum() == 6;
  });
  CHECK(ccount == 4);
}

TEST_CASE("mol.bonds()") {
  const auto m = "OC(=O)C(=O)O"_smiles;
  REQUIRE(m);
  unsigned int doubleBondCount = 0;
  for (const auto bond : m->bonds()) {
    if (bond->getBondType() == Bond::DOUBLE) {
      ++doubleBondCount;
    }
  }
  CHECK(doubleBondCount == 2);
  auto bonds = m->bonds();
  auto hasDoubleBond = std::any_of(
      bonds.begin(), bonds.end(),
      [](const auto bond) { return bond->getBondType() == Bond::DOUBLE; });
  CHECK(hasDoubleBond);
  doubleBondCount = std::count_if(
      bonds.begin(), bonds.end(),
      [](const auto bond) { return bond->getBondType() == Bond::DOUBLE; });
  CHECK(doubleBondCount == 2);
}

TEST_CASE("mol.atomNeighbors()") {
  const auto m = "CC(C)CO"_smiles;
  REQUIRE(m);
  unsigned int count = 0;
  for (const auto atom : m->atomNeighbors(m->getAtomWithIdx(1))) {
    count += atom->getDegree();
  }
  CHECK(count == 4);
  for (auto atom : m->atomNeighbors(m->getAtomWithIdx(1))) {
    atom->setAtomicNum(7);
  }
  MolOps::sanitizeMol(*m);
  CHECK(MolToSmiles(*m) == "NC(N)NO");
}

TEST_CASE("mol.atomBonds()") {
  const auto m = "CC(=C)CO"_smiles;
  REQUIRE(m);
  double count = 0;
  for (const auto bond : m->atomBonds(m->getAtomWithIdx(1))) {
    count += bond->getBondTypeAsDouble();
  }
  CHECK(count == 4);
  for (auto bond : m->atomBonds(m->getAtomWithIdx(1))) {
    bond->setBondType(Bond::BondType::SINGLE);
  }
  MolOps::sanitizeMol(*m);
  CHECK(MolToSmiles(*m) == "CC(C)CO");
}