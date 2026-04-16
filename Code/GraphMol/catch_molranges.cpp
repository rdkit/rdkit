//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <algorithm>
#include <ranges>

using namespace RDKit;

TEST_CASE("ranges") {
  std::unique_ptr<RWMol> m{new RWMol()};
  REQUIRE(m);
  //  = "CC(=C)C=O"_smiles;
  m->addAtom(new Atom(6), true, true);
  m->addAtom(new Atom(6), true, true);
  m->addAtom(new Atom(6), true, true);
  m->addAtom(new Atom(6), true, true);
  m->addAtom(new Atom(8), true, true);
  m->addBond(0, 1, Bond::SINGLE);
  m->addBond(1, 2, Bond::DOUBLE);
  m->addBond(1, 3, Bond::SINGLE);
  m->addBond(3, 4, Bond::DOUBLE);
  SECTION("atoms and bonds") {
    auto atoms = m->atoms();
    auto bonds = m->bonds();
    CHECK(std::ranges::distance(atoms) == 5);
    CHECK(std::ranges::distance(bonds) == 4);
    {
      std::vector<unsigned int> atomDegrees;
      std::ranges::transform(atoms, std::back_inserter(atomDegrees),
                             [](const auto atom) { return atom->getDegree(); });
      CHECK(atomDegrees == std::vector<unsigned int>{1, 3, 1, 2, 1});
    }
    {
      std::vector<unsigned int> atomDegrees;
      std::ranges::transform(atoms | std::views::reverse,
                             std::back_inserter(atomDegrees),
                             [](const auto atom) { return atom->getDegree(); });
      CHECK(atomDegrees == std::vector<unsigned int>{1, 2, 1, 3, 1});
    }
    {
      std::vector<Bond::BondType> bondOrders;
      std::ranges::transform(
          bonds, std::back_inserter(bondOrders),
          [](const auto bond) { return bond->getBondType(); });
      CHECK(bondOrders ==
            std::vector<Bond::BondType>{Bond::SINGLE, Bond::DOUBLE,
                                        Bond::SINGLE, Bond::DOUBLE});
    }
    {
      std::vector<Bond::BondType> bondOrders;
      std::ranges::transform(
          bonds | std::views::reverse, std::back_inserter(bondOrders),
          [](const auto bond) { return bond->getBondType(); });
      CHECK(bondOrders ==
            std::vector<Bond::BondType>{Bond::DOUBLE, Bond::SINGLE,
                                        Bond::DOUBLE, Bond::SINGLE});
    }
  }
  SECTION("Neighbors") {
    auto neighbors = m->atomNeighbors(m->getAtomWithIdx(1));
    CHECK(std::ranges::distance(neighbors) == 3);
    std::vector<unsigned int> neighborIndices;
    std::ranges::transform(neighbors, std::back_inserter(neighborIndices),
                           [](const auto atom) { return atom->getIdx(); });
    CHECK(neighborIndices == std::vector<unsigned int>{0, 2, 3});
    auto abonds = m->atomBonds(m->getAtomWithIdx(1));
    CHECK(std::ranges::distance(abonds) == 3);
    std::vector<unsigned int> bondIndices;
    std::ranges::transform(abonds, std::back_inserter(bondIndices),
                           [](const auto bond) { return bond->getIdx(); });
    CHECK(bondIndices == std::vector<unsigned int>{0, 1, 2});
  }
}

TEST_CASE("algorithms") {
  std::unique_ptr<RWMol> m{new RWMol()};
  REQUIRE(m);
  //  = "COC(F)C=C"_smiles;
  m->addAtom(new Atom(6), true, true);
  m->addAtom(new Atom(8), true, true);
  m->addAtom(new Atom(6), true, true);
  m->addAtom(new Atom(9), true, true);
  m->addAtom(new Atom(6), true, true);
  m->addAtom(new Atom(6), true, true);
  m->addBond(0, 1, Bond::SINGLE);
  m->addBond(1, 2, Bond::SINGLE);
  m->addBond(2, 3, Bond::SINGLE);
  m->addBond(4, 5, Bond::DOUBLE);
  m->addBond(2, 4, Bond::SINGLE);
  SECTION("atom count_if, filter, and take") {
    auto atoms = m->atoms();
    auto numC = std::ranges::count_if(
        atoms, [](const auto atom) { return atom->getAtomicNum() == 6; });
    CHECK(numC == 4);
    std::vector<unsigned int> atomIndices;
    std::ranges::transform(atoms | std::views::filter([](const auto atom) {
                             return atom->getAtomicNum() == 6;
                           }),
                           std::back_inserter(atomIndices),
                           [](const auto atom) { return atom->getIdx(); });
    CHECK(atomIndices == std::vector<unsigned int>{0, 2, 4, 5});
    atomIndices.clear();
    std::ranges::transform(atoms | std::views::filter([](const auto atom) {
                             return atom->getAtomicNum() == 6;
                           }) | std::views::take(2),
                           std::back_inserter(atomIndices),
                           [](const auto atom) { return atom->getIdx(); });
    CHECK(atomIndices == std::vector<unsigned int>{0, 2});
  }
  SECTION("bond count_if, filter, and take") {
    auto bonds = m->bonds();
    auto numSingle = std::ranges::count_if(bonds, [](const auto bond) {
      return bond->getBondType() == Bond::SINGLE;
    });
    CHECK(numSingle == 4);
    std::vector<unsigned int> bondIndices;
    std::ranges::transform(bonds | std::views::filter([](const auto bond) {
                             return bond->getBondType() == Bond::SINGLE;
                           }),
                           std::back_inserter(bondIndices),
                           [](const auto bond) { return bond->getIdx(); });
    CHECK(bondIndices == std::vector<unsigned int>{0, 1, 2, 4});
    bondIndices.clear();
    std::ranges::transform(bonds | std::views::filter([](const auto bond) {
                             return bond->getBondType() == Bond::SINGLE;
                           }) | std::views::take(2),
                           std::back_inserter(bondIndices),
                           [](const auto bond) { return bond->getIdx(); });
    CHECK(bondIndices == std::vector<unsigned int>{0, 1});
  }
}
