//
// Copyright (C) 2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/DummyAtom.h>
#include <RDGeneral/types.h>

using namespace RDKit;

TEST_CASE("createDummyAtom") {
  auto atom = createDummyAtom();
  REQUIRE(atom != nullptr);

  SECTION("atomic number is 0") { CHECK(atom->getAtomicNum() == 0); }

  SECTION("has a null query") { CHECK(atom->hasQuery()); }

  SECTION("no implicit H suppression") { CHECK(!atom->getNoImplicit()); }
}

TEST_CASE("isAttachmentPointDummy") {
  SECTION("returns false for a carbon atom") {
    Atom carbon(6);
    CHECK(!isAttachmentPointDummy(carbon));
  }

  SECTION("returns false for a dummy atom without atomLabel") {
    RWMol mol;
    mol.addAtom(new Atom(6), false, true);
    auto *dummy = new Atom(0);
    mol.addAtom(dummy, false, true);
    mol.addBond(0, 1, Bond::SINGLE);
    CHECK(!isAttachmentPointDummy(*mol.getAtomWithIdx(1)));
  }

  SECTION("returns false for a dummy atom with wrong label prefix") {
    RWMol mol;
    mol.addAtom(new Atom(6), false, true);
    auto *dummy = new Atom(0);
    mol.addAtom(dummy, false, true);
    mol.addBond(0, 1, Bond::SINGLE);
    mol.getAtomWithIdx(1)->setProp(common_properties::atomLabel,
                                   std::string("_R1"));
    CHECK(!isAttachmentPointDummy(*mol.getAtomWithIdx(1)));
  }

  SECTION("returns false when degree is not 1") {
    RWMol mol;
    mol.addAtom(new Atom(6), false, true);
    mol.addAtom(new Atom(6), false, true);
    auto *dummy = new Atom(0);
    mol.addAtom(dummy, false, true);
    mol.addBond(0, 2, Bond::SINGLE);
    mol.addBond(1, 2, Bond::SINGLE);
    mol.getAtomWithIdx(2)->setProp(common_properties::atomLabel,
                                   std::string("_AP1"));
    CHECK(!isAttachmentPointDummy(*mol.getAtomWithIdx(2)));
  }

  SECTION("returns true for a degree-1 dummy with _AP label") {
    RWMol mol;
    mol.addAtom(new Atom(6), false, true);
    auto *dummy = new Atom(0);
    mol.addAtom(dummy, false, true);
    mol.addBond(0, 1, Bond::SINGLE);
    mol.getAtomWithIdx(1)->setProp(common_properties::atomLabel,
                                   std::string("_AP1"));
    CHECK(isAttachmentPointDummy(*mol.getAtomWithIdx(1)));
  }

  SECTION("matches any _AP prefix, not just _AP1") {
    RWMol mol;
    mol.addAtom(new Atom(6), false, true);
    auto *dummy = new Atom(0);
    mol.addAtom(dummy, false, true);
    mol.addBond(0, 1, Bond::SINGLE);
    mol.getAtomWithIdx(1)->setProp(common_properties::atomLabel,
                                   std::string("_AP99"));
    CHECK(isAttachmentPointDummy(*mol.getAtomWithIdx(1)));
  }
}

TEST_CASE("makeNewRGroup") {
  SECTION("throws for rGroupNum == 0") {
    CHECK_THROWS_AS(makeNewRGroup(0), std::invalid_argument);
  }

  SECTION("creates atom with correct properties for R1") {
    auto atom = makeNewRGroup(1);
    REQUIRE(atom != nullptr);
    CHECK(atom->getAtomicNum() == 0);
    CHECK(atom->getIsotope() == 1);

    std::string label;
    REQUIRE(atom->getPropIfPresent(common_properties::atomLabel, label));
    CHECK(label == "_R1");

    std::string dLabel;
    REQUIRE(atom->getPropIfPresent(common_properties::dummyLabel, dLabel));
    CHECK(dLabel == "R1");

    unsigned int rLabel = 0;
    REQUIRE(atom->getPropIfPresent(common_properties::_MolFileRLabel, rLabel));
    CHECK(rLabel == 1u);
  }

  SECTION("creates atom with correct properties for R10") {
    auto atom = makeNewRGroup(10);
    REQUIRE(atom != nullptr);
    CHECK(atom->getIsotope() == 10);

    std::string label;
    REQUIRE(atom->getPropIfPresent(common_properties::atomLabel, label));
    CHECK(label == "_R10");

    std::string dLabel;
    REQUIRE(atom->getPropIfPresent(common_properties::dummyLabel, dLabel));
    CHECK(dLabel == "R10");

    unsigned int rLabel = 0;
    REQUIRE(atom->getPropIfPresent(common_properties::_MolFileRLabel, rLabel));
    CHECK(rLabel == 10u);
  }
}

TEST_CASE("getRGroupNumber") {
  SECTION("returns nullopt for a non-dummy atom") {
    Atom carbon(6);
    CHECK(!getRGroupNumber(&carbon).has_value());
  }

  SECTION("returns nullopt for a dummy atom without _MolFileRLabel") {
    Atom dummy(0);
    CHECK(!getRGroupNumber(&dummy).has_value());
  }

  SECTION("returns correct number for an R-group atom") {
    auto atom = makeNewRGroup(5);
    auto result = getRGroupNumber(atom.get());
    REQUIRE(result.has_value());
    CHECK(*result == 5u);
  }

  SECTION("round-trips through makeNewRGroup for various indices") {
    for (unsigned int n : {1u, 2u, 10u, 20u, 99u}) {
      auto atom = makeNewRGroup(n);
      auto result = getRGroupNumber(atom.get());
      REQUIRE(result.has_value());
      CHECK(*result == n);
    }
  }
}
