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
    std::ranges::transform(bonds, std::back_inserter(bondOrders),
                           [](const auto bond) { return bond->getBondType(); });
    CHECK(bondOrders == std::vector<Bond::BondType>{Bond::SINGLE, Bond::DOUBLE,
                                                    Bond::SINGLE,
                                                    Bond::DOUBLE});
  }
  {
    std::vector<Bond::BondType> bondOrders;
    std::ranges::transform(bonds | std::views::reverse,
                           std::back_inserter(bondOrders),
                           [](const auto bond) { return bond->getBondType(); });
    CHECK(bondOrders == std::vector<Bond::BondType>{Bond::DOUBLE, Bond::SINGLE,
                                                    Bond::DOUBLE,
                                                    Bond::SINGLE});
  }
}
