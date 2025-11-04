//
//  Copyright (C) 2025 NVIDIA Corporation & Affiliates and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>

#include <GraphMol/Chirality.h>
#include <GraphMol/RDMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/StereoGroup.h>

using ::RDKit::StereoGroup;
using ::RDKit::StereoGroups;
using ::RDKit::Chirality::cleanupStereoGroups;

namespace {

StereoGroups makeGroups() {
  StereoGroups groups;
  groups.addGroup(RDKit::StereoGroupType::STEREO_ABSOLUTE,
                  std::vector<uint32_t>{}, std::vector<uint32_t>{});
  groups.addGroup(RDKit::StereoGroupType::STEREO_OR,
                  std::vector<uint32_t>{0, 1, 2}, std::vector<uint32_t>{1},
                  /*readId=*/5);
  groups.addGroup(RDKit::StereoGroupType::STEREO_AND,
                  std::vector<uint32_t>{2, 3}, std::vector<uint32_t>{0, 2});
  return groups;
}

std::vector<uint32_t> iteratorToIndexVector(
    const std::pair<const uint32_t *, const uint32_t *> &iter) {
  std::vector<uint32_t> res;
  for (const uint32_t *idx = iter.first; idx != iter.second; idx++) {
    res.push_back(*idx);
  }
  return res;
}

}  // namespace

TEST_CASE("Empty groups") {
  StereoGroups groups;
  CHECK_THROWS_AS(groups.getGroupType(0), Invar::Invariant);
  CHECK_THROWS_AS(groups.getReadID(0), Invar::Invariant);
  CHECK_THROWS_AS(groups.getWriteID(0), Invar::Invariant);
  CHECK_THROWS_AS(groups.getAtoms(0), Invar::Invariant);
  CHECK_THROWS_AS(groups.getBonds(0), Invar::Invariant);
  CHECK_THROWS_AS(groups.GroupsEq(0, 0), Invar::Invariant);
  // Check no crashes
  groups.assignStereoGroupIds();
  groups.forwardStereoGroupIds();
  groups.removeAtomFromGroups(4);
  groups.removeGroupsWithAtoms({4});
}

TEST_CASE("Basic get/set") {
  StereoGroups groups = makeGroups();

  CHECK(groups.getGroupType(0) == RDKit::StereoGroupType::STEREO_ABSOLUTE);
  CHECK(groups.getGroupType(1) == RDKit::StereoGroupType::STEREO_OR);
  CHECK(groups.getGroupType(2) == RDKit::StereoGroupType::STEREO_AND);

  CHECK(groups.getReadID(0) == 0);
  CHECK(groups.getReadID(1) == 5);
  CHECK(groups.getReadID(2) == 0);

  CHECK(groups.getWriteID(0) == 0);
  CHECK(groups.getWriteID(1) == 0);
  CHECK(groups.getWriteID(2) == 0);

  auto group0Atoms = iteratorToIndexVector(groups.getAtoms(0));
  auto group1Atoms = iteratorToIndexVector(groups.getAtoms(1));
  auto group2Atoms = iteratorToIndexVector(groups.getAtoms(2));

  auto group0Bonds = iteratorToIndexVector(groups.getBonds(0));
  auto group1Bonds = iteratorToIndexVector(groups.getBonds(1));
  auto group2Bonds = iteratorToIndexVector(groups.getBonds(2));

  CHECK(group0Atoms.empty());
  CHECK(group1Atoms == std::vector<uint32_t>{0, 1, 2});
  CHECK(group2Atoms == std::vector<uint32_t>{2, 3});

  CHECK(group0Bonds.empty());
  CHECK(group1Bonds == std::vector<uint32_t>{1});
  CHECK(group2Bonds == std::vector<uint32_t>{0, 2});
}

TEST_CASE("Group equality") {
  StereoGroups groups = makeGroups();
  // Add a group equal to 2
  groups.addGroup(RDKit::StereoGroupType::STEREO_AND,
                  std::vector<uint32_t>{2, 3}, std::vector<uint32_t>{0, 2});
  // Different only in atoms
  groups.addGroup(RDKit::StereoGroupType::STEREO_AND, std::vector<uint32_t>{2},
                  std::vector<uint32_t>{0, 2});
  // Different only in bonds
  groups.addGroup(RDKit::StereoGroupType::STEREO_AND,
                  std::vector<uint32_t>{2, 3}, std::vector<uint32_t>{0});
  CHECK(groups.GroupsEq(0, 0));
  CHECK_FALSE(groups.GroupsEq(0, 1));
  CHECK_FALSE(groups.GroupsEq(0, 2));
  CHECK_FALSE(groups.GroupsEq(1, 0));
  CHECK(groups.GroupsEq(1, 1));
  CHECK_FALSE(groups.GroupsEq(1, 2));
  CHECK_FALSE(groups.GroupsEq(2, 0));
  CHECK_FALSE(groups.GroupsEq(2, 1));
  CHECK(groups.GroupsEq(2, 2));
  CHECK(groups.GroupsEq(2, 3));
  CHECK(!groups.GroupsEq(2, 4));
  CHECK(!groups.GroupsEq(2, 5));
}

TEST_CASE("Remove atom from groups") {
  StereoGroups groups = makeGroups();
  groups.removeAtomFromGroups(2);

  // removeAtomFromGroups removes the empty group 0, so there are only 2 groups left
  REQUIRE(groups.getNumGroups() == 2);

  auto group0Atoms = iteratorToIndexVector(groups.getAtoms(0));
  auto group1Atoms = iteratorToIndexVector(groups.getAtoms(1));

  CHECK(group0Atoms == std::vector<uint32_t>{0, 1});
  CHECK(group1Atoms == std::vector<uint32_t>{3});
}

TEST_CASE("Remove groups with atoms") {
  StereoGroups groups = makeGroups();
  groups.removeGroupsWithAtoms({1});

  // Only group 1, containing atom 1, should be removed
  REQUIRE(groups.getNumGroups() == 2);

  auto group0Atoms = iteratorToIndexVector(groups.getAtoms(0));
  auto group1Atoms = iteratorToIndexVector(groups.getAtoms(1));
  auto group0Bonds = iteratorToIndexVector(groups.getBonds(0));
  auto group1Bonds = iteratorToIndexVector(groups.getBonds(1));

  CHECK(group0Atoms.empty());
  CHECK(group1Atoms == std::vector<uint32_t>{2, 3});
  CHECK(group0Bonds.empty());
  CHECK(group1Bonds == std::vector<uint32_t>{0, 2});

  CHECK(groups.readIds[0] == 0);
  CHECK(groups.readIds[1] == 0);
}

TEST_CASE("Assign IDs") {
  StereoGroups groups = makeGroups();

  std::vector<StereoGroup> legacyGroups;
  legacyGroups.push_back(
      StereoGroup(RDKit::StereoGroupType::STEREO_ABSOLUTE, {}, {}));
  legacyGroups.push_back(
      StereoGroup(RDKit::StereoGroupType::STEREO_OR, {}, {}, 5));
  legacyGroups.push_back(
      StereoGroup(RDKit::StereoGroupType::STEREO_AND, {}, {}));

  RDKit::assignStereoGroupIds(legacyGroups);
  groups.assignStereoGroupIds();

  SECTION("Not preassigned") {
    CHECK(legacyGroups[0].getWriteId() == groups.getWriteID(0));
    CHECK(legacyGroups[1].getWriteId() == groups.getWriteID(1));
    CHECK(legacyGroups[2].getWriteId() == groups.getWriteID(2));
  }

  SECTION("Preassigned, new group") {
    groups.addGroup(RDKit::StereoGroupType::STEREO_OR, {}, {});
    legacyGroups.push_back(
        StereoGroup(RDKit::StereoGroupType::STEREO_OR, {}, {}));

    RDKit::assignStereoGroupIds(legacyGroups);
    groups.assignStereoGroupIds();
    CHECK(legacyGroups[0].getWriteId() == groups.getWriteID(0));
    CHECK(legacyGroups[1].getWriteId() == groups.getWriteID(1));
    CHECK(legacyGroups[2].getWriteId() == groups.getWriteID(2));
    CHECK(legacyGroups[3].getWriteId() == groups.getWriteID(3));
  }

  SECTION("Duplicates") {
    groups.addGroup(RDKit::StereoGroupType::STEREO_OR, {}, {});
    legacyGroups.push_back(
        StereoGroup(RDKit::StereoGroupType::STEREO_OR, {}, {}));
    legacyGroups[3].setWriteId(1);
    groups.setWriteID(3, 1);
    RDKit::assignStereoGroupIds(legacyGroups);
    groups.assignStereoGroupIds();
    CHECK(legacyGroups[0].getWriteId() == groups.getWriteID(0));
    CHECK(legacyGroups[1].getWriteId() == groups.getWriteID(1));
    CHECK(legacyGroups[2].getWriteId() == groups.getWriteID(2));
  }
}

TEST_CASE("Forward IDs") {
  StereoGroups groups = makeGroups();
  CHECK(groups.getWriteID(0) == 0);
  CHECK(groups.getWriteID(1) == 0);
  CHECK(groups.getWriteID(2) == 0);
  groups.forwardStereoGroupIds();
  CHECK(groups.getWriteID(0) == 0);
  CHECK(groups.getWriteID(1) == 5);
  CHECK(groups.getWriteID(2) == 0);
}

TEST_CASE("Chirality cleanup") {
  RDKit::RDMol mol;
  RDKit::SmilesParseTemp temp;
  static const char *basicSmiles = "CCCC";

  RDKit::SmilesParserParams params;
  params.sanitize = false;
  params.removeHs = false;

  REQUIRE(RDKit::SmilesToMol(basicSmiles, params, mol, temp) == true);
  REQUIRE(mol.getNumAtoms() == 4);
  REQUIRE(mol.getNumBonds() == 3);
  mol.getAtom(0).setChiralTag(RDKit::AtomEnums::ChiralType::CHI_TETRAHEDRAL_CW);
  mol.getAtom(1).setChiralTag(RDKit::AtomEnums::ChiralType::CHI_TETRAHEDRAL_CW);
  mol.getAtom(2).setChiralTag(RDKit::AtomEnums::ChiralType::CHI_TETRAHEDRAL_CW);
  mol.getAtom(3).setChiralTag(RDKit::AtomEnums::ChiralType::CHI_TETRAHEDRAL_CW);
  mol.getBond(0).setStereo(RDKit::BondEnums::BondStereo::STEREOATROPCCW);
  mol.getBond(1).setStereo(RDKit::BondEnums::BondStereo::STEREOATROPCW);
  mol.getBond(2).setStereo(RDKit::BondEnums::BondStereo::STEREOATROPCW);

  SECTION("Unset") {
    cleanupStereoGroups(mol);
    CHECK(mol.getStereoGroups() == nullptr);
  }

  SECTION("Empty") {
    auto groups = std::make_unique<StereoGroups>();
    mol.setStereoGroups(std::move(groups));
    cleanupStereoGroups(mol);
    CHECK(mol.getStereoGroups()->getNumGroups() == 0);
  }

  SECTION("No removals") {
    auto groups = std::make_unique<StereoGroups>();
    groups->addGroup(RDKit::StereoGroupType::STEREO_OR, {0, 1, 2}, {1}, 5);
    groups->addGroup(RDKit::StereoGroupType::STEREO_AND, {2, 3}, {0, 2});
    mol.setStereoGroups(std::move(groups));
    cleanupStereoGroups(mol);
    const auto *res = mol.getStereoGroups();
    REQUIRE(res != nullptr);
    REQUIRE(res->getNumGroups() == 2);

    auto group0Atoms = iteratorToIndexVector(res->getAtoms(0));
    auto group1Atoms = iteratorToIndexVector(res->getAtoms(1));
    auto group0Bonds = iteratorToIndexVector(res->getBonds(0));
    auto group1Bonds = iteratorToIndexVector(res->getBonds(1));

    CHECK(group0Atoms == std::vector<uint32_t>{0, 1, 2});
    CHECK(group1Atoms == std::vector<uint32_t>{2, 3});
    CHECK(group0Bonds == std::vector<uint32_t>{1});
    CHECK(group1Bonds == std::vector<uint32_t>{0, 2});
  }

  SECTION("Don't remove already empty group") {
    auto groups = std::make_unique<StereoGroups>();
    groups->addGroup(RDKit::StereoGroupType::STEREO_OR, {}, {}, 5);
    groups->addGroup(RDKit::StereoGroupType::STEREO_AND, {2, 3}, {0, 2});
    mol.setStereoGroups(std::move(groups));
    cleanupStereoGroups(mol);
    const auto *res = mol.getStereoGroups();
    REQUIRE(res != nullptr);
    REQUIRE(res->getNumGroups() == 2);

    auto group0Atoms = iteratorToIndexVector(res->getAtoms(0));
    auto group0Bonds = iteratorToIndexVector(res->getBonds(0));
    auto group1Atoms = iteratorToIndexVector(res->getAtoms(1));
    auto group1Bonds = iteratorToIndexVector(res->getBonds(1));

    CHECK(group0Atoms == std::vector<uint32_t>{});
    CHECK(group0Bonds == std::vector<uint32_t>{});
    CHECK(group1Atoms == std::vector<uint32_t>{2, 3});
    CHECK(group1Bonds == std::vector<uint32_t>{0, 2});
  }

  SECTION("Atom removals, no empty groups") {
    auto groups = std::make_unique<StereoGroups>();
    mol.getAtom(2).setChiralTag(RDKit::AtomEnums::ChiralType::CHI_UNSPECIFIED);
    groups->addGroup(RDKit::StereoGroupType::STEREO_OR, {1, 2}, {}, 5);
    groups->addGroup(RDKit::StereoGroupType::STEREO_AND, {2, 3}, {0, 2});
    mol.setStereoGroups(std::move(groups));
    cleanupStereoGroups(mol);
    const auto *res = mol.getStereoGroups();
    REQUIRE(res != nullptr);
    REQUIRE(res->getNumGroups() == 2);

    auto group0Atoms = iteratorToIndexVector(res->getAtoms(0));
    auto group1Atoms = iteratorToIndexVector(res->getAtoms(1));
    auto group0Bonds = iteratorToIndexVector(res->getBonds(0));
    auto group1Bonds = iteratorToIndexVector(res->getBonds(1));

    CHECK(group0Atoms == std::vector<uint32_t>{1});
    CHECK(group1Atoms == std::vector<uint32_t>{3});
    CHECK(group0Bonds == std::vector<uint32_t>{});
    CHECK(group1Bonds == std::vector<uint32_t>{0, 2});
  }

  SECTION("Atom removals, empty groups removed") {
    auto groups = std::make_unique<StereoGroups>();
    mol.getAtom(2).setChiralTag(RDKit::AtomEnums::ChiralType::CHI_UNSPECIFIED);
    mol.getAtom(3).setChiralTag(RDKit::AtomEnums::ChiralType::CHI_UNSPECIFIED);
    groups->addGroup(RDKit::StereoGroupType::STEREO_OR, {1, 2}, {}, 5);
    groups->addGroup(RDKit::StereoGroupType::STEREO_AND, {2, 3}, {});
    groups->addGroup(RDKit::StereoGroupType::STEREO_AND, {0, 1}, {2});

    mol.setStereoGroups(std::move(groups));
    cleanupStereoGroups(mol);
    const auto *res = mol.getStereoGroups();
    REQUIRE(res != nullptr);
    REQUIRE(res->getNumGroups() == 2);

    auto group0Atoms = iteratorToIndexVector(res->getAtoms(0));
    auto group1Atoms = iteratorToIndexVector(res->getAtoms(1));
    auto group0Bonds = iteratorToIndexVector(res->getBonds(0));
    auto group1Bonds = iteratorToIndexVector(res->getBonds(1));

    CHECK(group0Atoms == std::vector<uint32_t>{1});
    CHECK(group1Atoms == std::vector<uint32_t>{0, 1});
    CHECK(group0Bonds == std::vector<uint32_t>{});
    CHECK(group1Bonds == std::vector<uint32_t>{2});
  }

  SECTION("Atom removals, group removed even if bonds in group") {
    auto groups = std::make_unique<StereoGroups>();
    mol.getAtom(2).setChiralTag(RDKit::AtomEnums::ChiralType::CHI_UNSPECIFIED);
    mol.getAtom(3).setChiralTag(RDKit::AtomEnums::ChiralType::CHI_UNSPECIFIED);
    groups->addGroup(RDKit::StereoGroupType::STEREO_OR, {1, 2}, {}, 5);
    groups->addGroup(RDKit::StereoGroupType::STEREO_AND, {2, 3},
                     {1});  // Note that the bond should not keep the group alive.

    mol.setStereoGroups(std::move(groups));
    cleanupStereoGroups(mol);
    const auto *res = mol.getStereoGroups();
    REQUIRE(res != nullptr);
    REQUIRE(res->getNumGroups() == 1);

    auto group0Atoms = iteratorToIndexVector(res->getAtoms(0));
    auto group0Bonds = iteratorToIndexVector(res->getBonds(0));

    CHECK(group0Atoms == std::vector<uint32_t>{1});
    CHECK(group0Bonds == std::vector<uint32_t>{});
  }
}