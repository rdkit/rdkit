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

#include <GraphMol/RDMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Bond.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/QueryAtom.h>

using namespace RDKit;

// TODO - consolidate rdmol/romol/rwmol test setup
std::unique_ptr<RDMol> parseSmiles(const char *smiles, bool sanitize = false,
                                   bool removeHs = false) {
  auto mol = std::make_unique<RDMol>();
  RDKit::SmilesParseTemp temp;
  SmilesParserParams params;
  params.sanitize = sanitize;
  params.removeHs = removeHs;

  REQUIRE(RDKit::SmilesToMol(smiles, params, *mol, temp) == true);
  return mol;
}

std::unique_ptr<RDMol> basicMol() {
  auto mol = parseSmiles("CCC=C");
  REQUIRE(mol->getNumAtoms() == 4);
  REQUIRE(mol->getNumBonds() == 3);
  REQUIRE(mol->getBond(2).getBondType() == BondEnums::BondType::DOUBLE);

  return mol;
}


TEST_CASE("Add atom") {
  RWMol mol;
  REQUIRE(mol.addAtom() == 0);
  REQUIRE(mol.getNumAtoms() == 1);
  REQUIRE(mol.getAtomWithIdx(0)->getAtomicNum() == 0);
  REQUIRE(mol.addAtom() == 1);
  REQUIRE(mol.getNumAtoms() == 2);
  REQUIRE(mol.getAtomWithIdx(1)->getAtomicNum() == 0);

  SECTION("Update label") {
    REQUIRE(mol.addAtom(true) == 2);
    CHECK(mol.getAtomWithBookmark(ci_RIGHTMOST_ATOM) == mol.getAtomWithIdx(2));
  }

  SECTION("Don't update label") {
    REQUIRE(mol.addAtom(false) == 2);
    CHECK(mol.getAtomWithBookmark(ci_RIGHTMOST_ATOM) == mol.getAtomWithIdx(1));
  }
}

TEST_CASE("Add atom failures") {
  SECTION("Null atom") {
    RWMol mol;
    CHECK_THROWS_AS(mol.addAtom(nullptr), Invar::Invariant);
  }

  SECTION("Take ownership of owned atom") {
    RWMol mol1;
    mol1.addAtom();
    auto *atom = mol1.getAtomWithIdx(0);
    REQUIRE(atom != nullptr);

    RWMol mol2;
    CHECK_THROWS_AS(mol2.addAtom(atom, false, /*takeOwnership=*/true),
                    Invar::Invariant);
  }

  SECTION("Try to add atom owned by self") {
    RWMol mol;
    mol.addAtom();
    auto *atom = mol.getAtomWithIdx(0);
    REQUIRE(atom != nullptr);

    CHECK_THROWS_AS(mol.addAtom(atom, false, /*takeOwnership=*/true),
                    Invar::Invariant);
  }
}

TEST_CASE("Add atom via Atom* API") {
  RWMol mol;
  // Add one in first to make sure we're indexing ok.
  REQUIRE(mol.addAtom() == 0);

  auto atom = std::make_unique<Atom>(6);
  atom->setHybridization(Atom::SP3);
  atom->setFormalCharge(1);
  atom->setIsAromatic(true);
  atom->setNumExplicitHs(2);
  atom->setNoImplicit(true);
  atom->setIsotope(13);
  atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
  atom->setNumRadicalElectrons(1);

  atom->setProp("test", 3);

  SECTION("No ownership change") {
    REQUIRE(mol.addAtom(atom.get(), false, /*takeOwnership=*/false) == 1);
    Atom *newAtomHandle = mol.getAtomWithIdx(1);

    // First check that we've created a new atom.
    CHECK(newAtomHandle != atom.get());

    // Next check properties have copied over.
    CHECK(newAtomHandle->getAtomicNum() == 6);
    CHECK(newAtomHandle->getHybridization() == Atom::SP3);
    CHECK(newAtomHandle->getFormalCharge() == 1);
    CHECK(newAtomHandle->getIsAromatic());
    CHECK(newAtomHandle->getNumExplicitHs() == 2);
    CHECK(newAtomHandle->getNoImplicit());
    CHECK(newAtomHandle->getIsotope() == 13);
    CHECK(newAtomHandle->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
    CHECK(newAtomHandle->getNumRadicalElectrons() == 1);

    // Finally make sure we don't depend on the old value.
    atom.reset();
    CHECK(newAtomHandle->getAtomicNum() == 6);
    CHECK(newAtomHandle->getProp<int>("test") == 3);
  }

  SECTION("Ownership change") {
    REQUIRE(mol.addAtom(atom.get(), false, /*takeOwnership=*/true) == 1);
    Atom *newAtomHandle = mol.getAtomWithIdx(1);

    // First check that we're using the same handle.
    CHECK(newAtomHandle == atom.release());

    // Next check properties have copied over.
    CHECK(newAtomHandle->getAtomicNum() == 6);
    CHECK(newAtomHandle->getHybridization() == Atom::SP3);
    CHECK(newAtomHandle->getFormalCharge() == 1);
    CHECK(newAtomHandle->getIsAromatic());
    CHECK(newAtomHandle->getNumExplicitHs() == 2);
    CHECK(newAtomHandle->getNoImplicit());
    CHECK(newAtomHandle->getIsotope() == 13);
    CHECK(newAtomHandle->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
    CHECK(newAtomHandle->getNumRadicalElectrons() == 1);

    CHECK(newAtomHandle->getProp<int>("test") == 3);
  }

  SECTION("Update label") {
    REQUIRE(mol.addAtom(atom.get(), true) == 1);
    CHECK(mol.getAtomWithBookmark(ci_RIGHTMOST_ATOM) == mol.getAtomWithIdx(1));
  }

  SECTION("Don't update label") {
    REQUIRE(mol.addAtom(atom.get(), false) == 1);
    CHECK(mol.getAtomWithBookmark(ci_RIGHTMOST_ATOM) != mol.getAtomWithIdx(1));
  }

  SECTION("Transferred monomer info") {
    atom->setMonomerInfo(new AtomPDBResidueInfo("dummy", /*serial=*/3));
    REQUIRE(mol.addAtom(atom.get(), false) == 1);
    auto monomerInfo = mol.getAtomWithIdx(1)->getMonomerInfo();
    REQUIRE(monomerInfo != nullptr);
    auto castInfo = dynamic_cast<AtomPDBResidueInfo *>(monomerInfo);
    REQUIRE(castInfo != nullptr);
    CHECK(castInfo->getSerialNumber() == 3);
  }
}

std::vector<int> getAtomNeighborIndices(const RWMol& mol, int atomIdx) {
  std::vector<int> neighbors;
  for (const auto& bondIter:  boost::make_iterator_range(mol.getAtomBonds(mol.getAtomWithIdx(atomIdx)))) {
    const auto &nbr = (mol)[bondIter];
    neighbors.push_back(nbr->getOtherAtomIdx(atomIdx));
  }
  return neighbors;
}

template<typename pointerT>
std::vector<int> stdListToVec(const std::list<pointerT*>& lst) {
  std::vector<int> res;
  for (const auto& element: lst) {
    res.push_back(element->getIdx());
  }
  return res;
}

TEST_CASE("Remove Atom") {
  auto molPtr = basicMol();
  RWMol &rwmol = molPtr->asRWMol();


  SECTION("Basic, no bond removal, no stereochem") {
    // Change up atom numbers for better testing.
    rwmol.getAtomWithIdx(1)->setAtomicNum(7);
    rwmol.getAtomWithIdx(2)->setAtomicNum(8);

    REQUIRE(rwmol.getNumAtoms() == 4);
    REQUIRE(rwmol.getNumBonds() == 3);
    REQUIRE(rwmol.getAtomWithIdx(0)->getAtomicNum() == 6);
    REQUIRE(rwmol.getAtomWithIdx(1)->getAtomicNum() == 7);
    REQUIRE(rwmol.getBondWithIdx(0)->getBeginAtomIdx() == 0);
    REQUIRE(rwmol.getBondWithIdx(0)->getEndAtomIdx() == 1);
    rwmol.removeAtom((unsigned int)0);

    CHECK(rwmol.getNumAtoms() == 3);
    CHECK(rwmol.getNumBonds() == 2);

    CHECK(rwmol.getAtomWithIdx(0)->getAtomicNum() == 7);
    CHECK(rwmol.getAtomWithIdx(1)->getAtomicNum() == 8);

    // Check that bond atom refs are updated.
    CHECK(rwmol.getBondWithIdx(0)->getBeginAtomIdx() == 0);
    CHECK(rwmol.getBondWithIdx(0)->getEndAtomIdx() == 1);

    // Check atom neighbor iteration unchanged except for indices
    std::vector<int> neighbors0 = getAtomNeighborIndices(rwmol, 0);
    std::vector<int> neighbors1 = getAtomNeighborIndices(rwmol, 1);
    std::vector<int> neighbors2 = getAtomNeighborIndices(rwmol, 2);
    CHECK_THAT(neighbors0, Catch::Matchers::UnorderedEquals(std::vector<int>{1}));
    CHECK_THAT(neighbors1, Catch::Matchers::UnorderedEquals(std::vector<int>{0, 2}));
    CHECK_THAT(neighbors2, Catch::Matchers::UnorderedEquals(std::vector<int>{1}));
  }

  SECTION("Removes associated bonds") {
    REQUIRE(rwmol.getNumAtoms() == 4);
    REQUIRE(rwmol.getNumBonds() == 3);

    rwmol.removeAtom(1);

    CHECK(rwmol.getNumAtoms() == 3);
    CHECK(rwmol.getNumBonds() == 1);

    std::vector<int> neighbors0 = getAtomNeighborIndices(rwmol, 0);
    std::vector<int> neighbors1 = getAtomNeighborIndices(rwmol, 1);
    std::vector<int> neighbors2 = getAtomNeighborIndices(rwmol, 2);
    CHECK_THAT(neighbors0, Catch::Matchers::UnorderedEquals(std::vector<int>{}));
    CHECK_THAT(neighbors1, Catch::Matchers::UnorderedEquals(std::vector<int>{2}));
    CHECK_THAT(neighbors2, Catch::Matchers::UnorderedEquals(std::vector<int>{1}));
  }

  SECTION("Stereo groups") {
    Atom* atom0 = rwmol.getAtomWithIdx(0);
    Atom* atom1 = rwmol.getAtomWithIdx(1);
    Atom* atom2 = rwmol.getAtomWithIdx(2);

    std::vector<StereoGroup> groups;
    std::vector<Atom *> atoms = {atom0, atom2, atom1};
    groups.push_back(StereoGroup(StereoGroupType::STEREO_OR, std::move(atoms), {}));
    std::vector<Atom *> atoms2 = {rwmol.getAtomWithIdx(1)};
    groups.push_back(StereoGroup(StereoGroupType::STEREO_AND, std::move(atoms2), {}));
    rwmol.setStereoGroups(std::move(groups));

    rwmol.removeAtom(1);

    const auto& resultGroups = rwmol.getStereoGroups();
    REQUIRE(resultGroups.size() == 1);
    const auto& group = resultGroups[0];

    CHECK(group.getGroupType() == StereoGroupType::STEREO_OR);
    CHECK_THAT(group.getAtoms(), Catch::Matchers::UnorderedEquals(std::vector<Atom *>({atom0, atom2})));
  }

  SECTION("Bookmarks") {
    rwmol.setAtomBookmark(rwmol.getAtomWithIdx(0), 0);
    rwmol.setAtomBookmark(rwmol.getAtomWithIdx(1), 0);
    rwmol.setAtomBookmark(rwmol.getAtomWithIdx(2), 0);
    rwmol.setAtomBookmark(rwmol.getAtomWithIdx(3), 0);
    rwmol.setAtomBookmark(rwmol.getAtomWithIdx(1), 1);
    rwmol.setAtomBookmark(rwmol.getAtomWithIdx(2), 1);
    rwmol.setAtomBookmark(rwmol.getAtomWithIdx(3), 2);

    rwmol.setBondBookmark(rwmol.getBondWithIdx(0), 0);
    rwmol.setBondBookmark(rwmol.getBondWithIdx(1), 0);
    rwmol.setBondBookmark(rwmol.getBondWithIdx(2), 0);
    rwmol.setBondBookmark(rwmol.getBondWithIdx(0), 1);
    rwmol.setBondBookmark(rwmol.getBondWithIdx(1), 2);
    rwmol.setBondBookmark(rwmol.getBondWithIdx(2), 2);

    // Before deletion, we have {bookmark: indices}
    // atoms:
    // 0: 0, 1, 2, 3
    // 1: 1, 2
    // 2: 3
    // bonds:
    // 0: 0, 1, 2
    // 1: 0, 1
    // 2: 2

    // We are deleting bonds 0 and 1, atom 1. Bond 2 now maps to 0. Atoms 2 and 3 now map to 1 and 2.
    rwmol.removeAtom(1);

    auto& atomBookmarks0 = rwmol.getAllAtomsWithBookmark(0);
    auto& atomBookmarks1 = rwmol.getAllAtomsWithBookmark(1);
    auto& atomBookmarks2 = rwmol.getAllAtomsWithBookmark(2);

    auto& bondBookmarks0 = rwmol.getAllBondsWithBookmark(0);
    CHECK_THROWS_AS(rwmol.getAllBondsWithBookmark(1), Invar::Invariant);
    auto& bondBookmarks2 = rwmol.getAllBondsWithBookmark(2);

    CHECK(stdListToVec(atomBookmarks0) == std::vector<int>{0, 1, 2});
    CHECK(stdListToVec(atomBookmarks1) == std::vector<int>{1});
    CHECK(stdListToVec(atomBookmarks2) == std::vector<int>{2});

    CHECK(stdListToVec(bondBookmarks0) == std::vector<int>{0});
    CHECK(stdListToVec(bondBookmarks2) == std::vector<int>{0});
  }

  SECTION("Properties") {
    for (int i = 0; i < 4; ++i) {
      rwmol.getAtomWithIdx(i)->setProp("test", i);
    }
    for (int i = 0; i < 3; ++i) {
      rwmol.getBondWithIdx(i)->setProp("test", i);
    }

    rwmol.removeAtom(1);
    CHECK(rwmol.getAtomWithIdx(0)->getProp<int>("test") == 0);
    CHECK(rwmol.getAtomWithIdx(1)->getProp<int>("test") == 2);
    CHECK(rwmol.getAtomWithIdx(2)->getProp<int>("test") == 3);
    CHECK(rwmol.getBondWithIdx(0)->getProp<int>("test") == 2);
  }

  SECTION("Does not clear computed props") {
    rwmol.getAtomWithIdx(0)->setProp("computed", 3, /*computed=*/ true);
    rwmol.setProp("computed", 3, /*computed=*/ true);
    rwmol.removeAtom(rwmol.getAtomWithIdx(1), /*clearComputedProps=*/false);
    CHECK(rwmol.getAtomWithIdx(0)->getProp<int>("computed") == 3);
    CHECK(rwmol.getProp<int>("computed") == 3);
  }

  SECTION("Clears computed props") {
    rwmol.getAtomWithIdx(0)->setProp("computed", 3, /*computed=*/ true);
    rwmol.setProp("computed", 3, /*computed=*/ true);
    rwmol.removeAtom(rwmol.getAtomWithIdx(1), /*clearComputedProps=*/true);
    CHECK(!rwmol.getAtomWithIdx(0)->hasProp("computed"));
    CHECK(!rwmol.hasProp("computed"));
  }

  SECTION("Conformers") {
    auto conf1 = std::make_unique<Conformer>(rwmol.getNumAtoms());
    auto conf2 = std::make_unique<Conformer>(rwmol.getNumAtoms());

    const int numAtoms = 4;
    for (int i = 0; i < numAtoms; i ++) {
      conf1->setAtomPos(i, RDGeom::Point3D(i, i, i));
      conf2->setAtomPos(i, RDGeom::Point3D(i + numAtoms, i + numAtoms, i + numAtoms));
    }

    uint32_t id1 = rwmol.addConformer(conf1.release(), true);
    uint32_t id2 = rwmol.addConformer(conf2.release(), true);

    rwmol.removeAtom(1);

    const Conformer* confResult = &rwmol.getConformer(id1);
    const Conformer* conf2Result = &rwmol.getConformer(id2);

    REQUIRE(rwmol.getNumConformers() == 2);
    REQUIRE(confResult != nullptr);
    REQUIRE(conf2Result != nullptr);

    CHECK(confResult->getNumAtoms() == 3);
    CHECK(conf2Result->getNumAtoms() == 3);

    CHECK(confResult->getAtomPos(0).x == 0);
    CHECK(confResult->getAtomPos(1).x == 2);
    CHECK(confResult->getAtomPos(2).x == 3);
    CHECK(conf2Result->getAtomPos(0).x == 4);
    CHECK(conf2Result->getAtomPos(1).x == 6);
    CHECK(conf2Result->getAtomPos(2).x == 7);
  }

  SECTION("Substance groups") {
    std::vector<SubstanceGroup> &sgs = getSubstanceGroups(rwmol);
    sgs.emplace_back(&rwmol, "REMOVE_HAS_ATOM");
    sgs.back().addAtomWithIdx(0);
    sgs.back().addAtomWithIdx(1);
    sgs.back().addAtomWithIdx(2);
    sgs.emplace_back(&rwmol, "REMOVE_HAS_REMOVED_BOND");
    sgs.back().addBondWithIdx(0);
    sgs.back().addBondWithIdx(2);
    sgs.emplace_back(&rwmol, "DO_NOT_REMOVE");
    sgs.back().addAtomWithIdx(0);
    sgs.back().addAtomWithIdx(2);
    sgs.back().addBondWithIdx(2);

    rwmol.removeAtom(1);
    std::vector<SubstanceGroup> &resultGroups = getSubstanceGroups(rwmol);
    REQUIRE(resultGroups.size() == 1);
    CHECK(resultGroups[0].getProp<std::string>("TYPE") == "DO_NOT_REMOVE");
    // Make sure atoms and bonds are shifted
    CHECK(resultGroups[0].getAtoms() == std::vector<unsigned int>{0, 1});
    CHECK(resultGroups[0].getBonds() == std::vector<unsigned int>{0});
  }
}

TEST_CASE("Replace atom") {
  RWMol mol;
  mol.addAtom();
  mol.addAtom();
  mol.getAtomWithIdx(1)->setAtomicNum(6);
  Atom* existingAtom = mol.getAtomWithIdx(1);

  auto atom = std::make_unique<Atom>(7);
  Atom* newAtom = atom.get();
  
  SECTION("Fails - null input") {
    CHECK_THROWS_AS(mol.replaceAtom(0, nullptr), Invar::Invariant);
  }

  SECTION("Fails - bad index") {
    CHECK_THROWS_AS(mol.replaceAtom(2, atom.get()), Invar::Invariant);
  }

  SECTION("Basics") {
    // replaceAtom does not take ownership, it makes a copy.
    mol.replaceAtom(1, newAtom);
    CHECK(mol.getAtomWithIdx(1) != newAtom);
    CHECK(mol.getAtomWithIdx(1)->getAtomicNum() == 7);
  }

  SECTION("Overwrites properties") {
    existingAtom->setProp("test", 3);
    existingAtom->setProp("should_be_deleted", 1);
    newAtom->setProp("otherprop", 5.0);
    newAtom->setProp("test", 4);
    mol.replaceAtom(1, newAtom);
    CHECK(mol.getAtomWithIdx(1)->getProp<int>("test") == 4);
    CHECK(mol.getAtomWithIdx(1)->getProp<double>("otherprop") == 5.0);
    CHECK(!mol.getAtomWithIdx(1)->hasProp("should_be_deleted"));
  }

  SECTION("Does not overwrite properties when preserved") {
    existingAtom->setProp("test", 3);
    newAtom->setProp("otherprop", 5.0);
    newAtom->setProp("test", 4);
    mol.replaceAtom(1, newAtom, false, /*preserveProps=*/true);
    CHECK(mol.getAtomWithIdx(1)->getProp<int>("test") == 3);
    CHECK(!mol.getAtomWithIdx(1)->hasProp("otherprop"));
  }

  SECTION("Updates bookmarks") {
    mol.setAtomBookmark(existingAtom, 3);
    mol.replaceAtom(1, newAtom);
    CHECK(mol.getAtomWithBookmark(3) == mol.getAtomWithIdx(1));
  }

  SECTION("Copies monomer info") {
    newAtom->setMonomerInfo(new AtomPDBResidueInfo("dummy", /*serial=*/3));
    mol.replaceAtom(1, newAtom);
    auto monomerInfo = mol.getAtomWithIdx(1)->getMonomerInfo();
    REQUIRE(monomerInfo != nullptr);
    auto castInfo = dynamic_cast<AtomPDBResidueInfo *>(monomerInfo);
    REQUIRE(castInfo != nullptr);
    CHECK(castInfo->getSerialNumber() == 3);
  }

  SECTION("Clears existing monomer info") {
    existingAtom->setMonomerInfo(new AtomPDBResidueInfo("dummy", /*serial=*/3));
    mol.replaceAtom(1, newAtom);
    CHECK(mol.getAtomWithIdx(1)->getMonomerInfo() == nullptr);
  }
}

TEST_CASE("Replace bond") {
  auto rdmol = basicMol();
  auto& rwmol = rdmol->asRWMol();
  Bond* existingBond = rwmol.getBondWithIdx(1);
  auto newBond = std::make_unique<Bond>();
  newBond->setBondType(Bond::BondType::TRIPLE);

  std::vector<SubstanceGroup> &sgs = getSubstanceGroups(rwmol);
  sgs.emplace_back(&rwmol, "TEST");
  sgs.back().addBondWithIdx(1);
  sgs.back().addBondWithIdx(2);
  sgs.emplace_back(&rwmol, "UNRELATED");
  sgs.back().addBondWithIdx(0);

  SECTION("Fails - null input") {
    CHECK_THROWS_AS(rwmol.replaceBond(0, nullptr), Invar::Invariant);
  }

  SECTION("Fails - bad index") {
    CHECK_THROWS_AS(rwmol.replaceBond(3, newBond.get()), Invar::Invariant);
  }

  SECTION("Basics") {
    // replacerwmol does not take ownership, it makes a copy.
    rwmol.replaceBond(1, newBond.get());
    CHECK(rwmol.getBondWithIdx(1) != newBond.get());
    CHECK(rwmol.getBondWithIdx(1)->getBondType() == Bond::BondType::TRIPLE);
  }

  SECTION("Overwrites properties") {
    existingBond->setProp("test", 3);
    existingBond->setProp("should_be_deleted", 1);
    newBond->setProp("otherprop", 5.0);
    newBond->setProp("test", 4);
    rwmol.replaceBond(1, newBond.get());
    CHECK(rwmol.getBondWithIdx(1)->getProp<int>("test") == 4);
    CHECK(rwmol.getBondWithIdx(1)->getProp<double>("otherprop") == 5.0);
    CHECK(!rwmol.getBondWithIdx(1)->hasProp("should_be_deleted"));
  }

  SECTION("Does not overwrite properties when preserved") {
    existingBond->setProp("test", 3);
    newBond->setProp("otherprop", 5.0);
    newBond->setProp("test", 4);
    rwmol.replaceBond(1, newBond.get(), /*preserveProps=*/true);
    CHECK(rwmol.getBondWithIdx(1)->getProp<int>("test") == 3);
    CHECK(!rwmol.getBondWithIdx(1)->hasProp("otherprop"));
  }

  SECTION("Updates bookmarks") {
    rwmol.setBondBookmark(existingBond, 3);
    rwmol.replaceBond(1, newBond.get());
    CHECK(rwmol.getBondWithBookmark(3) == rwmol.getBondWithIdx(1));
  }

  SECTION("Keeps substance groups") {
    rwmol.replaceBond(1, newBond.get(), false, /*keepSGroups=*/true);
    CHECK(sgs.size() == 2);
    CHECK(sgs[0].getProp<std::string>("TYPE") == "TEST");
  }

  SECTION("Removes substance groups referencing bond") {
    rwmol.replaceBond(1, newBond.get(), false, /*keepSGroups=*/false);
    CHECK(sgs.size() == 1);
    CHECK(sgs[0].getProp<std::string>("TYPE") == "UNRELATED");
  }

  SECTION("Explicit hydrogen handling - increase order") {
    // Set one atom to have fewer explicit Hs than order difference
    rwmol.getAtomWithIdx(1)->setNumExplicitHs(1);
    // Set one atom to have greater explict Hs than order difference
    // (not realistic but exercises code path)
    rwmol.getAtomWithIdx(2)->setNumExplicitHs(3);
    // Order difference of two
    REQUIRE(rwmol.getBondWithIdx(1)->getBondType() == Bond::BondType::SINGLE);
    rwmol.replaceBond(1, newBond.get());
    CHECK(rwmol.getAtomWithIdx(1)->getNumExplicitHs() == 0);
    CHECK(rwmol.getAtomWithIdx(2)->getNumExplicitHs() == 1);
  }

  SECTION("Explicit hydrogen handling - decrease order does nothing") {
    rwmol.getAtomWithIdx(1)->setNumExplicitHs(1);
    rwmol.getAtomWithIdx(2)->setNumExplicitHs(3);
    // Order difference of -1
    rwmol.getBondWithIdx(1)->setBondType(Bond::BondType::DOUBLE);
    newBond->setBondType(Bond::BondType::SINGLE);
    rwmol.replaceBond(1, newBond.get());
    CHECK(rwmol.getAtomWithIdx(1)->getNumExplicitHs() == 1);
    CHECK(rwmol.getAtomWithIdx(2)->getNumExplicitHs() == 3);
  }
}

// Expects a bond add between atoms 1 and 3
void checkBondConnections(const RWMol& mol, const Bond* bond) {
  REQUIRE(bond != nullptr);

  // Check neighbors of each atom
  CHECK_THAT(getAtomNeighborIndices(mol, 0),
             Catch::Matchers::UnorderedEquals(std::vector<int>({1})));
  CHECK_THAT(getAtomNeighborIndices(mol, 1),
             Catch::Matchers::UnorderedEquals(std::vector<int>({0, 2, 3})));
  CHECK_THAT(getAtomNeighborIndices(mol, 2),
             Catch::Matchers::UnorderedEquals(std::vector<int>({1, 3})));
  CHECK_THAT(getAtomNeighborIndices(mol, 3),
             Catch::Matchers::UnorderedEquals(std::vector<int>({1, 2})));
}

TEST_CASE("Add bond") {
  auto rdmol = basicMol();
  auto& rwmol = rdmol->asRWMol();
  REQUIRE(rwmol.getNumAtoms() == 4);
  REQUIRE(rwmol.getNumBonds() == 3);

  auto bondToAdd = std::make_unique<Bond>();
  REQUIRE(bondToAdd != nullptr);

  bondToAdd->setBondDir(Bond::BondDir::ENDDOWNRIGHT);
  bondToAdd->setBondType(Bond::BondType::DOUBLE);
  bondToAdd->setIsConjugated(true);
  bondToAdd->setStereo(Bond::BondStereo::STEREOCIS);
  bondToAdd->setProp("test", 3);

  SECTION("Index API fails, bad atom numbers") {
      CHECK_THROWS_AS(rwmol.addBond(3, 5, Bond::BondType::SINGLE), Invar::Invariant);
      CHECK_THROWS_AS(rwmol.addBond(5, 3, Bond::BondType::SINGLE), Invar::Invariant);
  }

  SECTION("Index API fails, same atoms for bond") {
      CHECK_THROWS_AS(rwmol.addBond(1, 1, Bond::BondType::SINGLE), Invar::Invariant);
  }

  SECTION("Index API fails, bond exists") {
      CHECK_THROWS_AS(rwmol.addBond(0, 1, Bond::BondType::SINGLE), Invar::Invariant);
  }

  SECTION("Index API succeeds") {
    uint32_t idx = rwmol.addBond(3, 1, Bond::BondType::TRIPLE);
    const Bond* bond = rwmol.getBondWithIdx(idx - 1);
    checkBondConnections(rwmol, bond);

    CHECK(bond->getBeginAtomIdx() == 3);
    CHECK(bond->getEndAtomIdx() == 1);
    CHECK(bond->getBondType() == Bond::BondType::TRIPLE);
  }

  SECTION("Index API with aromatic bond sets aromatic atoms") {
    uint32_t idx = rwmol.addBond(1, 3, Bond::BondType::AROMATIC);
    const Bond* bond = rwmol.getBondWithIdx(idx - 1);
    checkBondConnections(rwmol, bond);

    auto* atom1 = rwmol.getAtomWithIdx(1);
    auto* atom3 = rwmol.getAtomWithIdx(3);
    CHECK(bond->getIsAromatic());
    CHECK(atom1->getIsAromatic());
    CHECK(atom3->getIsAromatic());
  }

  SECTION("Atom* API fails, null atoms") {
    auto atom = std::make_unique<Atom>();
    CHECK_THROWS_AS(rwmol.addBond(atom.get(), nullptr, Bond::BondType::SINGLE), Invar::Invariant);
    CHECK_THROWS_AS(rwmol.addBond(nullptr, atom.get(), Bond::BondType::SINGLE), Invar::Invariant);
  }

  SECTION("Atom API succeeds") {
    auto* atom1 = rwmol.getAtomWithIdx(1);
    auto* atom2 = rwmol.getAtomWithIdx(3);
    REQUIRE(atom1 != nullptr);
    REQUIRE(atom2 != nullptr);

    uint32_t idx = rwmol.addBond(atom1, atom2, Bond::BondType::TRIPLE);
    // New number of bonds
    CHECK(idx == 4);
    auto* bond = rwmol.getBondWithIdx(idx - 1);
    checkBondConnections(rwmol, bond);
  }

  SECTION("Bond* API fails, null bond") {
    CHECK_THROWS_AS(rwmol.addBond(nullptr), Invar::Invariant);
  }

  SECTION("Bond API fails, takeownership when already owned by mol") {
    Bond* bond = rwmol.getBondWithIdx(0);
    REQUIRE(bond != nullptr);
    CHECK_THROWS_AS(rwmol.addBond(bond, /*takeOwnership=*/true), Invar::Invariant);
  }

  SECTION("Bond API fails, bad bond indices") {
    bondToAdd->setBeginAtomIdx(5);
    bondToAdd->setEndAtomIdx(3);
    CHECK_THROWS_AS(rwmol.addBond(bondToAdd.get(), /*takeOwnership=*/false), Invar::Invariant);
  }

  SECTION("Bond API fails, self indices") {
    bondToAdd->setBeginAtomIdx(1);
    bondToAdd->setEndAtomIdx(1);
    CHECK_THROWS_AS(rwmol.addBond(bondToAdd.get(), /*takeOwnership=*/false), Invar::Invariant);
  }

  SECTION("Bond API fail, existing bond in mol") {
    REQUIRE(rwmol.getBondBetweenAtoms(1, 2) != nullptr);
    bondToAdd->setBeginAtomIdx(1);
    bondToAdd->setEndAtomIdx(2);
    CHECK_THROWS_AS(rwmol.addBond(bondToAdd.get(), /*takeOwnership=*/false), Invar::Invariant);
  }

  SECTION("Bond API succeeds, don't take ownership") {
    bondToAdd->setBeginAtomIdx(1);
    bondToAdd->setEndAtomIdx(3);
    Bond* bondPtrOrig = bondToAdd.get();
    uint32_t idx = rwmol.addBond(bondToAdd.get(), /*takeOwnership=*/false);
    auto* bond = rwmol.getBondWithIdx(idx - 1);
    checkBondConnections(rwmol, bond);
    CHECK(bond != bondPtrOrig);

    CHECK(bond->getBondDir() == Bond::BondDir::ENDDOWNRIGHT);
    CHECK(bond->getBondType() == Bond::BondType::DOUBLE);
    CHECK(bond->getIsConjugated());
    CHECK(bond->getStereo() == Bond::BondStereo::STEREOCIS);
    CHECK(bond->getProp<int>("test") == 3);
  }

  SECTION("Bond API succeeds, does not update aromaticity") {
    bondToAdd->setBeginAtomIdx(1);
    bondToAdd->setEndAtomIdx(3);
    bondToAdd->setBondType(Bond::BondType::AROMATIC);
    uint32_t idx = rwmol.addBond(bondToAdd.get(), /*takeOwnership=*/false);
    auto* bond = rwmol.getBondWithIdx(idx - 1);
    checkBondConnections(rwmol, bond);
    auto* atom1 = rwmol.getAtomWithIdx(1);
    auto* atom3 = rwmol.getAtomWithIdx(3);
    CHECK(!bond->getIsAromatic());
    CHECK(!atom1->getIsAromatic());
    CHECK(!atom3->getIsAromatic());
  }

  SECTION("Bond API succeeds, take ownership") {
    bondToAdd->setBeginAtomIdx(1);
    bondToAdd->setEndAtomIdx(3);
    Bond* bondPtrOrig = bondToAdd.get();
    uint32_t idx = rwmol.addBond(bondToAdd.release(), /*takeOwnership=*/true);
    auto* bond = rwmol.getBondWithIdx(idx - 1);
    CHECK(bond == bondPtrOrig);
    checkBondConnections(rwmol, bond);

    CHECK(bond->getBondDir() == Bond::BondDir::ENDDOWNRIGHT);
    CHECK(bond->getBondType() == Bond::BondType::DOUBLE);
    CHECK(bond->getIsConjugated());
    CHECK(bond->getStereo() == Bond::BondStereo::STEREOCIS);
    CHECK(bond->getProp<int>("test") == 3);
  }
}

void addTestConformer(RWMol& mol, int startVal) {
  auto conf = new Conformer(mol.getNumAtoms());
  // Set each atom value to its index + startVal
  for (uint32_t i = 0; i < mol.getNumAtoms(); ++i) {
    conf->setAtomPos(i, RDGeom::Point3D(i + startVal, i + startVal, i + startVal));
  }
  mol.addConformer(conf, /*assignId=*/true);
}

TEST_CASE("InsertMol") {
  auto rdmol = basicMol();
  auto& rwmol = rdmol->asRWMol();
  rwmol.getAtomWithIdx(0)->setProp("test_prev_atom", 5.0);

  std::vector<SubstanceGroup> &sgs = getSubstanceGroups(rwmol);
  sgs.emplace_back(&rwmol, "TEST_MOL1");
  sgs.back().addAtomWithIdx(0);
  sgs.back().addBondWithIdx(1);

  std::vector<StereoGroup> groups;
  std::vector<Atom *> atoms = {rwmol.getAtomWithIdx(1)};
  groups.push_back(StereoGroup(StereoGroupType::STEREO_OR, std::move(atoms),
                               {}));
  rwmol.setStereoGroups(std::move(groups));

  RWMol mol2;
  auto atom1 = std::make_unique<Atom>(6);
  atom1->setProp("test_new_atom", 3);
  auto atom2 = std::make_unique<Atom>(7);
  auto atom3 = std::make_unique<Atom>(8);
  auto atom4 = std::make_unique<Atom>(9);

  mol2.addAtom(atom1.get());
  mol2.addAtom(atom2.get());
  mol2.addAtom(atom3.get());
  mol2.addAtom(atom4.get());
  mol2.addBond(0, 1, Bond::BondType::SINGLE);
  mol2.addBond(1, 2, Bond::BondType::DOUBLE);
  mol2.addBond(2, 3, Bond::BondType::SINGLE);
  mol2.getBondWithIdx(1)->setProp("test", 4);
  mol2.getBondWithIdx(1)->setStereoAtoms(0, 3);

  SECTION("Basics") {
    rwmol.insertMol(mol2);

    // Check that the new mol has been inserted
    REQUIRE(rwmol.getNumAtoms() == 8);
    REQUIRE(rwmol.getNumBonds() == 6);
    CHECK(rwmol.getAtomWithIdx(5)->getAtomicNum() == 7);
    CHECK(rwmol.getBondWithIdx(4)->getBondType() == Bond::BondType::DOUBLE);

    // Check bond neighbor indices
    CHECK_THAT(getAtomNeighborIndices(rwmol, 3),
             Catch::Matchers::UnorderedEquals(std::vector<int>({2})));
    CHECK_THAT(getAtomNeighborIndices(rwmol, 4),
               Catch::Matchers::UnorderedEquals(std::vector<int>({5})));
    CHECK_THAT(getAtomNeighborIndices(rwmol, 6),
               Catch::Matchers::UnorderedEquals(std::vector<int>({5, 7})));

    // Check properties
    CHECK(rwmol.getAtomWithIdx(0)->getProp<double>("test_prev_atom") == 5.0);
    CHECK(rwmol.getAtomWithIdx(4)->getProp<int>("test_new_atom") == 3);
    CHECK(rwmol.getBondWithIdx(4)->getProp<int>("test") == 4);

    CHECK(rwmol.getNumConformers() == 0);
  }

  SECTION("_ringStereoAtoms prop update") {
    mol2.getAtomWithIdx(3)->setProp(common_properties::_ringStereoAtoms, std::vector<int>({0, -1, 2}));
    rwmol.insertMol(mol2);

    std::vector<int> gotStereoAtoms;
    REQUIRE(rwmol.getAtomWithIdx(7)->getPropIfPresent(common_properties::_ringStereoAtoms, gotStereoAtoms));
    CHECK(gotStereoAtoms == std::vector<int>({4, -5, 6}));
  }

  SECTION("Stereo atoms") {
    mol2.getBondWithIdx(1)->setStereoAtoms(0, 3);
    rwmol.insertMol(mol2);
    const auto& result = rwmol.getBondWithIdx(4)->getStereoAtoms();
    CHECK(result == std::vector<int>{4, 7});
  }

  SECTION("Stereo groups") {
    std::vector<StereoGroup> groups2;
    std::vector<Atom *> atoms2 = {mol2.getAtomWithIdx(3)};
    groups2.push_back(StereoGroup(StereoGroupType::STEREO_AND, std::move(atoms2),
                                 {}));
    mol2.setStereoGroups(std::move(groups2));

    rwmol.insertMol(mol2);

    const auto& result = rwmol.getStereoGroups();
    REQUIRE(result.size() == 2);
    CHECK(result[0].getGroupType() == StereoGroupType::STEREO_OR);
    CHECK(result[1].getGroupType() == StereoGroupType::STEREO_AND);

    CHECK(result[0].getAtoms() == std::vector<Atom *>({rwmol.getAtomWithIdx(1)}));
    CHECK(result[1].getAtoms() == std::vector<Atom *>({rwmol.getAtomWithIdx(7)}));
  }

  SECTION("Substance groups") {
    std::vector<SubstanceGroup> &sgs2 = getSubstanceGroups(mol2);
    sgs2.emplace_back(&rwmol, "TEST_MOL2");
    sgs2.back().addAtomWithIdx(1);
    sgs2.back().addBondWithIdx(2);

    rwmol.insertMol(mol2);
    const std::vector<SubstanceGroup> &result = getSubstanceGroups(rwmol);
    REQUIRE(result.size() == 2);

    CHECK(result[0].getProp<std::string>("TYPE") == "TEST_MOL1");
    CHECK(result[1].getProp<std::string>("TYPE") == "TEST_MOL2");

    CHECK(result[0].getAtoms() == std::vector<unsigned int>{0});
    CHECK(result[0].getBonds() == std::vector<unsigned int>{1});
    CHECK(result[1].getAtoms() == std::vector<unsigned int>{5});
    CHECK(result[1].getBonds() == std::vector<unsigned int>{5});
  }

  SECTION("Conformers - only in original") {
    addTestConformer(rwmol, 0);
    rwmol.insertMol(mol2);

    REQUIRE(rwmol.getNumConformers() == 1);
    const auto& conf = rwmol.getConformer();
    CHECK(conf.getNumAtoms() == 8);
    CHECK(conf.getAtomPos(0).x == 0.0);
    CHECK(conf.getAtomPos(3).y == 3.0);
    // New atoms should be zeroed out.
    CHECK(conf.getAtomPos(4).x == 0.0);
    CHECK(conf.getAtomPos(7).z == 0.0);
  }

  SECTION("Conformers - only in mol to be merged") {
    addTestConformer(mol2, 4.0);
    rwmol.insertMol(mol2);

    REQUIRE(rwmol.getNumConformers() == 1);
    const auto& conf = rwmol.getConformer();
    CHECK(conf.getNumAtoms() == 8);
    CHECK(conf.getAtomPos(0).x == 0.0);
    CHECK(conf.getAtomPos(3).y == 0.0);

    CHECK(conf.getAtomPos(4).x == 4.0);
    CHECK(conf.getAtomPos(7).z == 7.0);
  }

  SECTION("Conformers - in both") {
    addTestConformer(rwmol, 0.0);
    addTestConformer(mol2, 4.0);

    rwmol.insertMol(mol2);
    REQUIRE(rwmol.getNumConformers() == 1);
    const auto& conf = rwmol.getConformer();
    CHECK(conf.getNumAtoms() == 8);
    CHECK(conf.getAtomPos(0).x == 0.0);
    CHECK(conf.getAtomPos(3).y == 3.0);

    CHECK(conf.getAtomPos(4).x == 4.0);
    CHECK(conf.getAtomPos(7).z == 7.0);
  }

  SECTION("Conformers in both, more in original") {
    addTestConformer(rwmol, 0.0);
    addTestConformer(rwmol, 10.0);
    addTestConformer(mol2, 4.0);

    rwmol.insertMol(mol2);
    REQUIRE(rwmol.getNumConformers() == 2);

    const auto& conf1 = rwmol.getConformer(0);
    CHECK(conf1.getNumAtoms() == 8);
    CHECK(conf1.getAtomPos(0).x == 0.0);
    CHECK(conf1.getAtomPos(3).y == 3.0);
    CHECK(conf1.getAtomPos(4).x == 0.0); // New atoms not copied in this case.

    const auto& conf2 = rwmol.getConformer(1);
    CHECK(conf2.getNumAtoms() == 8);
    CHECK(conf2.getAtomPos(0).x == 10.0);
    CHECK(conf2.getAtomPos(3).y == 13.0);
    CHECK(conf2.getAtomPos(4).z == 0.0);
  }

  SECTION("Conformers in both, more in mol to be merged") {
    addTestConformer(rwmol, 0.0);
    addTestConformer(mol2, 4.0);
    addTestConformer(mol2, 8.0);

    rwmol.insertMol(mol2);
    REQUIRE(rwmol.getNumConformers() == 1);
    const auto& conf = rwmol.getConformer();
    CHECK(conf.getNumAtoms() == 8);
    CHECK(conf.getAtomPos(0).x == 0.0);
    CHECK(conf.getAtomPos(3).y == 3.0);
    CHECK(conf.getAtomPos(4).x == 0.0); // New atoms not copied in this case.
  }
}

TEST_CASE("RWMol assignment operator") {
  SECTION("Basic molecule assignment") {
    auto mol1ptr = parseSmiles("CCC=C");
    RWMol& mol1 = mol1ptr->asRWMol();
    RWMol mol2;
    mol2 = mol1;

    CHECK(mol2.getNumAtoms() == 4);
    CHECK(mol2.getNumBonds() == 3);
    CHECK(mol2.getBondWithIdx(2)->getBondType() == Bond::BondType::DOUBLE);
    
    // Check atoms are properly copied
    for (unsigned int i = 0; i < mol2.getNumAtoms(); ++i) {
      CHECK(mol2.getAtomWithIdx(i)->getAtomicNum() == mol1.getAtomWithIdx(i)->getAtomicNum());
    }
    mol1ptr.reset(); // Check lifetimes.
    CHECK(mol2.getNumAtoms() == 4);
    CHECK(mol2.getNumBonds() == 3);
    CHECK(mol2.getBondWithIdx(2)->getBondType() == Bond::BondType::DOUBLE);
  }

  SECTION("Assignment with atom properties") {
    auto mol1ptr = parseSmiles("CC");
    RWMol& mol1 = mol1ptr->asRWMol();
    mol1.getAtomWithIdx(0)->setProp("test_prop", 42);
    mol1.getAtomWithIdx(1)->setProp("other_prop", "value");
    
    RWMol mol2;
    mol2 = mol1;
    mol1ptr.reset(); // Check lifetimes.

    CHECK(mol2.getAtomWithIdx(0)->getProp<int>("test_prop") == 42);
    CHECK(mol2.getAtomWithIdx(1)->getProp<std::string>("other_prop") == "value");
  }

  SECTION("Assignment with conformers") {
    auto mol1ptr = parseSmiles("CC");
    RWMol& mol1 = mol1ptr->asRWMol();
    addTestConformer(mol1, 1);
    addTestConformer(mol1, 2);
    
    RWMol mol2;
    mol2 = mol1;
    mol1ptr.reset(); // Check lifetimes.

    CHECK(mol2.getNumConformers() == 2);
    
    // Check conformer coordinates
    const auto conf1 = mol2.getConformer(0);
    CHECK(conf1.getAtomPos(0).x == 1.0);
    CHECK(conf1.getAtomPos(0).y == 1.0);
    CHECK(conf1.getAtomPos(0).z == 1.0);
    
    const auto conf2 = mol2.getConformer(1);
    CHECK(conf2.getAtomPos(0).x == 2.0);
    CHECK(conf2.getAtomPos(0).y == 2.0);
    CHECK(conf2.getAtomPos(0).z == 2.0);
  }

  SECTION("Assignment with stereo information") {
    auto mol1ptr = parseSmiles("C[C@H](F)Cl");
    RWMol& mol1 = mol1ptr->asRWMol();
    RWMol mol2;
    mol2 = mol1;

    CHECK(mol2.getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
    
    // Add stereo groups to source molecule
    std::vector<Atom*> atoms = {mol1.getAtomWithIdx(1)};
    std::vector<Bond*> bonds;
    std::vector<StereoGroup> groups;
    groups.emplace_back(StereoGroupType::STEREO_OR, atoms, bonds);
    mol1.setStereoGroups(std::move(groups));
    
    mol2 = mol1;
    const auto& resultGroups = mol2.getStereoGroups();
    CHECK(resultGroups.size() == 1);
    CHECK(resultGroups[0].getGroupType() == StereoGroupType::STEREO_OR);
    CHECK(resultGroups[0].getAtoms().size() == 1);
  }

  SECTION("Assignment with bookmarks") {
    auto mol1ptr = parseSmiles("CC");
    RWMol& mol1 = mol1ptr->asRWMol();
    mol1.setAtomBookmark(mol1.getAtomWithIdx(0), 1);
    mol1.setBondBookmark(mol1.getBondWithIdx(0), 2);
    
    RWMol mol2;
    mol2 = mol1;
    mol1ptr.reset(); // Check lifetimes.

    CHECK(mol2.hasAtomBookmark(1));
    CHECK(mol2.hasAtomBookmark(1));
    CHECK(mol2.getAtomWithBookmark(1)->getIdx() == 0);
    CHECK(mol2.getBondWithBookmark(2)->getIdx() == 0);
  }

// FIXME: It's a good idea to test this, but at least one build fails
// flagging self-assignment as invalid. Check if there's a less direct
// way to test this that the compiler won't identify.
#if 0
  SECTION("Self assignment") {
    auto mol1ptr = parseSmiles("CCC");
    RWMol& mol1 = mol1ptr->asRWMol();
    mol1 = mol1;

    CHECK(mol1.getNumAtoms() == 3);
    CHECK(mol1.getNumBonds() == 2);
  }
#endif

  SECTION("Assignment with queries") {
    auto mol1ptr = parseSmiles("CC");
    RWMol& mol1 = mol1ptr->asRWMol();
    auto query = std::make_unique<QueryAtom>(6);
    mol1.replaceAtom(0, query.get());
    
    RWMol mol2;
    mol2 = mol1;
    mol1ptr.reset(); // Check lifetimes.

    CHECK(mol2.getAtomWithIdx(0)->hasQuery());
    CHECK(!mol2.getAtomWithIdx(1)->hasQuery());
  }

  SECTION("Previously owned assignment") {
    auto mol1ptr = parseSmiles("CC");
    RWMol& mol1 = mol1ptr->asRWMol();
    auto mol2ptr = parseSmiles("CCC");
    RWMol& mol2 = mol2ptr->asRWMol();

    mol2 = mol1;
    mol1ptr.reset(); // Check lifetimes.

    CHECK(mol2.getNumAtoms() == 2);
    CHECK(mol2.getNumBonds() == 1);
  }

  SECTION("Both self-owned assignment") {
        auto mol1Ptr = std::make_unique<RWMol>();
        RWMol& mol1 = *mol1Ptr;
        RWMol mol2;
        auto atom = std::make_unique<Atom>(6);
        auto atom2 = std::make_unique<Atom>(7);
        mol1.addAtom(atom.get());
        mol1.addAtom(atom2.get());
        mol1.addBond(0, 1, Bond::BondType::SINGLE);
        mol2 = mol1;
        mol1Ptr.reset();

        CHECK(mol2.getNumAtoms() == 2);
        CHECK(mol2.getNumBonds() == 1);
  }
}

TEST_CASE("remove atom with MonomerInfo") {
  SECTION("Basics") {
    RWMol mol;
    auto atom1 = std::make_unique<Atom>(6);
    auto atom2 = std::make_unique<Atom>(7);
    atom1->setMonomerInfo(new AtomPDBResidueInfo("ALA1", /*serial=*/1));
    atom2->setMonomerInfo(new AtomPDBResidueInfo("ALA2", /*serial=*/2));
    mol.addAtom(atom1.get());
    mol.addAtom(atom2.get());
    mol.addBond(0, 1, Bond::BondType::SINGLE);

    mol.removeAtom(0u);
    REQUIRE(mol.getNumAtoms() == 1);
    REQUIRE(mol.getNumBonds() == 0);
    CHECK(mol.getAtomWithIdx(0)->getAtomicNum() == 7);
    CHECK(mol.getAtomWithIdx(0)->getMonomerInfo() != nullptr);
    CHECK(mol.getAtomWithIdx(0)->getMonomerInfo()->getName() == "ALA2");
  }
}