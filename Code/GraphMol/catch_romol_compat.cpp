//
//  Copyright (C) 2025 NVIDIA Corporation & Affiliates and other RDKit
//  contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <thread>

#include <catch2/catch_all.hpp>

#include <GraphMol/RDMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MonomerInfo.h>

using namespace RDKit;

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

struct MyCustomClass {
  static size_t copyCount;
  int value;
  MyCustomClass() : value(0) {}
  MyCustomClass(int initValue) : value(initValue) {}
  MyCustomClass(const MyCustomClass &that) : value(that.value) { ++copyCount; }
  MyCustomClass &operator=(const MyCustomClass &that) noexcept {
    value = that.value;
    ++copyCount;
    return *this;
  }
};
size_t MyCustomClass::copyCount = 0;

TEST_CASE("ROMol self initialization") {
  ROMol mol;
  CHECK(mol.getNumAtoms() == 0);
  CHECK(mol.getNumBonds() == 0);
  CHECK(mol.getNumHeavyAtoms() == 0);

  CHECK_THROWS_AS(mol.getAtomWithIdx(0), Invar::Invariant);
  CHECK_THROWS_AS(mol.getBondWithIdx(0), Invar::Invariant);
}

TEST_CASE("ROMol move constructor") {
  ROMol mol;
  RDMol *rdmol = &mol.asRDMol();
  auto *conf = new Conformer();
  mol.addConformer(conf, true);
  CHECK(&conf->getOwningMol() == &mol);
  mol.setProp("test", 1);
  ROMol mol2(std::move(mol));
  CHECK(&mol2.asRDMol() == rdmol);
  // This value is never used, but GCC can no longer deduce that it's
  // initialized, and maybe-uninitialized is currently treated as an error.
  int res = 0;
  CHECK(mol2.getPropIfPresent("test", res));
  CHECK(res == 1);
  CHECK(&conf->getOwningMol() == &mol2);
  CHECK(&rdmol->asROMol() == &mol2);
}

TEST_CASE("ROMol move assignment operator") {
  ROMol mol;
  RDMol *rdmol = &mol.asRDMol();
  auto *conf = new Conformer();
  mol.addConformer(conf, true);
  CHECK(&conf->getOwningMol() == &mol);
  mol.setProp("test", 1);
  ROMol mol2 = std::move(mol);
  CHECK(&mol2.asRDMol() == rdmol);
  // This value is never used, but GCC can no longer deduce that it's
  // initialized, and maybe-uninitialized is currently treated as an error.
  int res = 0;
  CHECK(mol2.getPropIfPresent("test", res));
  CHECK(res == 1);
  CHECK(&conf->getOwningMol() == &mol2);
  CHECK(&rdmol->asROMol() == &mol2);
}

TEST_CASE("ROMol move constructor - previous owner") {
  auto molPtr = basicMol();
  RDMol &mol = *molPtr;
  molPtr.release();  // Will be moved from.
  ROMol &romol = mol.asROMol();
  romol.setProp("test", 1);
  Atom *atom = romol.getAtomWithIdx(0);
  CHECK(&atom->getOwningMol() == &romol);
  ROMol mol2(std::move(romol));
  CHECK(&mol2.asRDMol() == &mol);
  CHECK(romol.getNumAtoms() == 0);
  CHECK(romol.getNumBonds() == 0);
  CHECK(mol2.getNumAtoms() == 4);
  CHECK(mol2.getNumBonds() == 3);
  CHECK(mol2.getAtomWithIdx(0) == atom);
  // This value is never used, but GCC can no longer deduce that it's
  // initialized, and maybe-uninitialized is currently treated as an error.
  int res = 0;
  CHECK(mol2.getPropIfPresent("test", res));
  CHECK(res == 1);
  CHECK(&atom->getOwningMol() == &mol2);
  CHECK(&mol2.asRDMol().asROMol() == &mol2);
}

TEST_CASE("ROMol copy constructor") {
  auto mol = std::make_unique<ROMol>();

  mol->setProp("prop", 3);
  ROMol mol2(*mol);
  mol.reset();
  // This value is never used, but GCC can no longer deduce that it's
  // initialized, and maybe-uninitialized is currently treated as an error.
  int res = 0;
  REQUIRE(mol2.getPropIfPresent<int>("prop", res));
  CHECK(res == 3);
}

template <typename MolT>
void createConformer(MolT &mol) {
  auto *conf = new Conformer(mol.getNumAtoms());
  conf->set3D(false);
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    conf->setAtomPos(i, RDGeom::Point3D(i, i, i));
  }
  mol.addConformer(conf, true);
}

void checkConformer(const ROMol &mol) {
  const Conformer &conf = mol.getConformer();
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    const RDGeom::Point3D &pos = conf.getAtomPos(i);
    CHECK(pos.x == i);
    CHECK(pos.y == i);
    CHECK(pos.z == i);
  }
}

TEST_CASE("ROMol copy constructor with owning RDMol") {
  auto rdmol = basicMol();
  auto &mol = rdmol->asROMol();

  createConformer(mol);
  mol.setProp("prop", 3);
  Atom *at = mol.getAtomWithIdx(1);
  at->setProp("atomprop", 5);
  mol.setAtomBookmark(at, 42);
  mol.setBondBookmark(mol.getBondWithIdx(0), 43);
  mol.getBondWithIdx(1)->setStereoAtoms(0, 3);

  // Stereo groups
  std::vector<StereoGroup> groups;
  std::vector<Atom *> atoms = {mol.getAtomWithIdx(1)};
  groups.push_back(
      StereoGroup(StereoGroupType::STEREO_OR, std::move(atoms), {}));
  mol.setStereoGroups(groups);

  // Substance groups
  std::vector<SubstanceGroup> &sgs = getSubstanceGroups(mol);
  sgs.emplace_back(&mol, "TEST");

  ROMol mol2(mol);

  mol.getConformer();
  mol2.getConformer();

  // This value is never used, but GCC can no longer deduce that it's
  // initialized, and maybe-uninitialized is currently treated as an error.
  int res = 0;
  REQUIRE(mol2.getPropIfPresent<int>("prop", res));
  CHECK(res == 3);

  REQUIRE(mol2.getAtomWithIdx(1)->getPropIfPresent<int>("atomprop", res));
  CHECK(res == 5);

  // Check that setting one doesn't affect the other now
  at->setProp("atomprop", 6);
  REQUIRE(mol2.getAtomWithIdx(1)->getPropIfPresent<int>("atomprop", res));
  CHECK(res == 5);

  CHECK(mol2.getAtomWithBookmark(42) == mol2.getAtomWithIdx(1));
  CHECK(mol2.getBondWithBookmark(43) == mol2.getBondWithIdx(0));

  checkConformer(mol);
  checkConformer(mol2);

  const std::vector<int> want = {0, 3};
  CHECK(mol2.getBondWithIdx(1)->getStereoAtoms() == want);

  // Get stereo groups
  const std::vector<StereoGroup> &resultStereoGroups = mol2.getStereoGroups();
  CHECK(resultStereoGroups.size() == 1);
  CHECK(resultStereoGroups[0].getGroupType() == StereoGroupType::STEREO_OR);
  CHECK(resultStereoGroups[0].getAtoms() ==
        std::vector<Atom *>{mol2.getAtomWithIdx(1)});

  // Get substance groups
  const std::vector<SubstanceGroup> &resultSubstanceGroups =
      getSubstanceGroups(mol2);
  CHECK(resultSubstanceGroups.size() == 1);
  CHECK(resultSubstanceGroups[0].getProp<std::string>("TYPE") == "TEST");
}

TEST_CASE("ROMol stable when RDMol is moved") {
  auto rdmol = basicMol();
  auto &mol = rdmol->asROMol();

  Atom *at = mol.getAtomWithIdx(1);
  at->setProp("atomprop", 5);

  RDMol rdmol2(std::move(*rdmol));
  // This value is never used, but GCC can no longer deduce that it's
  // initialized, and maybe-uninitialized is currently treated as an error.
  int res = 0;
  REQUIRE(mol.getAtomWithIdx(1)->getPropIfPresent<int>("atomprop", res));
  CHECK(res == 5);
}

TEST_CASE("ROMol Compat Bookmarks") {
  RDMol rdmol;
  RDKit::SmilesParseTemp temp;
  static const char *basicSmiles = "CCCC";

  SmilesParserParams params;
  params.sanitize = false;
  params.removeHs = false;

  REQUIRE(RDKit::SmilesToMol(basicSmiles, params, rdmol, temp) == true);
  REQUIRE(rdmol.getNumAtoms() == 4);
  REQUIRE(rdmol.getNumBonds() == 3);
  auto &mol = rdmol.asROMol();

  auto *atom0 = mol.getAtomWithIdx(0);
  auto *atom1 = mol.getAtomWithIdx(1);
  auto *atom2 = mol.getAtomWithIdx(2);
  auto *atom3 = mol.getAtomWithIdx(3);

  auto *bond0 = mol.getBondWithIdx(0);
  auto *bond1 = mol.getBondWithIdx(1);
  auto *bond2 = mol.getBondWithIdx(2);

  SECTION("Basic bookmarking") {
    mol.setAtomBookmark(atom1, 0);
    mol.setAtomBookmark(atom3, -6);
    mol.setBondBookmark(bond1, 0);

    CHECK(mol.getAtomWithBookmark(0) == atom1);
    CHECK(mol.getAtomWithBookmark(-6) == atom3);
    CHECK(mol.getBondWithBookmark(0) == bond1);
  }

  SECTION("Multiple atoms in bookmark") {
    mol.setAtomBookmark(atom1, 0);
    mol.setAtomBookmark(atom2, 0);
    mol.setAtomBookmark(atom3, 0);

    CHECK(mol.getAtomWithBookmark(0) == atom1);

    const std::list<Atom *> &res = mol.getAllAtomsWithBookmark(0);
    const std::list<Atom *> want = {atom1, atom2, atom3};
    CHECK(res == want);
  }

  SECTION("Multiple bonds in bookmark") {
    mol.setBondBookmark(bond2, 0);
    mol.setBondBookmark(bond1, 0);
    mol.setBondBookmark(bond0, 0);

    CHECK(mol.getBondWithBookmark(0) == bond2);

    const std::list<Bond *> &res = mol.getAllBondsWithBookmark(0);
    const std::list<Bond *> want = {bond2, bond1, bond0};
    CHECK(res == want);
  }

  SECTION("Replace atom bookmarks") {
    mol.setAtomBookmark(atom3, 0);
    mol.setAtomBookmark(atom2, 0);
    mol.setAtomBookmark(atom1, 0);

    mol.replaceAtomBookmark(atom0, 0);

    CHECK(mol.getAtomWithBookmark(0) == atom0);

    const std::list<Atom *> &res = mol.getAllAtomsWithBookmark(0);
    const std::list<Atom *> want = {atom0};
    CHECK(res == want);
  }

  SECTION("Multiple bookmarks same atom is fine") {
    mol.setAtomBookmark(atom1, 2);
    mol.setAtomBookmark(atom1, 1);
    mol.setAtomBookmark(atom1, 0);

    CHECK(mol.getAtomWithBookmark(0) == atom1);
    CHECK(mol.getAtomWithBookmark(1) == atom1);
    CHECK(mol.getAtomWithBookmark(2) == atom1);

    const std::list<Atom *> &res = mol.getAllAtomsWithBookmark(0);
    const std::list<Atom *> want = {atom1};
    CHECK(res == want);
  }

  SECTION("Missing bookmark throws") {
    mol.setAtomBookmark(atom0, 1);
    mol.setBondBookmark(bond0, 2);

    CHECK_THROWS_AS(mol.getAtomWithBookmark(0), Invar::Invariant);
    CHECK_THROWS_AS(mol.getBondWithBookmark(0), Invar::Invariant);

    CHECK_THROWS_AS(mol.getAllAtomsWithBookmark(0), Invar::Invariant);
    CHECK_THROWS_AS(mol.getAllBondsWithBookmark(0), Invar::Invariant);
  }

  SECTION("Clear atom bookmarks") {
    mol.setAtomBookmark(atom1, 0);
    mol.setAtomBookmark(atom1, 1);
    mol.setAtomBookmark(atom1, 2);

    mol.clearAtomBookmark(0);
    CHECK_THROWS_AS(mol.getAtomWithBookmark(0), Invar::Invariant);
    CHECK(mol.getAtomWithBookmark(1) == atom1);
    CHECK(mol.getAtomWithBookmark(2) == atom1);

    mol.clearAllAtomBookmarks();
    CHECK_THROWS_AS(mol.getAtomWithBookmark(0), Invar::Invariant);
    CHECK_THROWS_AS(mol.getAtomWithBookmark(1), Invar::Invariant);
    CHECK_THROWS_AS(mol.getAtomWithBookmark(2), Invar::Invariant);
  }

  SECTION("Clear bond bookmarks") {
    mol.setBondBookmark(bond1, 0);
    mol.setBondBookmark(bond1, 1);
    mol.setBondBookmark(bond1, 2);

    mol.clearBondBookmark(0);
    CHECK_THROWS_AS(mol.getBondWithBookmark(0), Invar::Invariant);
    CHECK(mol.getBondWithBookmark(1) == bond1);
    CHECK(mol.getBondWithBookmark(2) == bond1);

    mol.clearAllBondBookmarks();
    CHECK_THROWS_AS(mol.getBondWithBookmark(0), Invar::Invariant);
    CHECK_THROWS_AS(mol.getBondWithBookmark(1), Invar::Invariant);
    CHECK_THROWS_AS(mol.getBondWithBookmark(2), Invar::Invariant);
  }

  SECTION("Has bookmark") {
    mol.setAtomBookmark(atom1, 0);
    mol.setAtomBookmark(atom1, 1);
    mol.setAtomBookmark(atom1, 2);

    CHECK(mol.hasAtomBookmark(0));
    CHECK(mol.hasAtomBookmark(1));
    CHECK(mol.hasAtomBookmark(2));
    CHECK(!mol.hasAtomBookmark(3));

    mol.setBondBookmark(bond1, 0);
    mol.setBondBookmark(bond1, 1);
    mol.setBondBookmark(bond1, 2);

    CHECK(mol.hasBondBookmark(0));
    CHECK(mol.hasBondBookmark(1));
    CHECK(mol.hasBondBookmark(2));
    CHECK(!mol.hasBondBookmark(3));

    mol.clearAtomBookmark(0);
    CHECK(!mol.hasAtomBookmark(0));
  }

  SECTION("Full atom bookmark array compat") {
    mol.setAtomBookmark(atom2, 0);
    mol.setAtomBookmark(atom0, 0);
    mol.setAtomBookmark(atom1, 5);

    std::map<int, std::list<Atom *>> *atomBookmarks = mol.getAtomBookmarks();
    REQUIRE(atomBookmarks->size() == 2);
    REQUIRE(atomBookmarks->find(0) != atomBookmarks->end());
    std::list<Atom *> &res1 = atomBookmarks->at(0);
    std::list<Atom *> want1 = {atom2, atom0};
    CHECK(res1 == want1);
    REQUIRE(atomBookmarks->find(5) != atomBookmarks->end());
    std::list<Atom *> &res2 = atomBookmarks->at(5);
    std::list<Atom *> want2 = {atom1};
    CHECK(res2 == want2);
  }

  SECTION("Full bond bookmark array compat") {
    mol.setBondBookmark(bond2, 0);
    mol.setBondBookmark(bond0, 0);
    mol.setBondBookmark(bond1, 5);

    std::map<int, std::list<Bond *>> *bondBookmarks = mol.getBondBookmarks();
    REQUIRE(bondBookmarks->size() == 2);
    REQUIRE(bondBookmarks->find(0) != bondBookmarks->end());
    std::list<Bond *> &res1 = bondBookmarks->at(0);
    std::list<Bond *> want1 = {bond2, bond0};
    CHECK(res1 == want1);
    REQUIRE(bondBookmarks->find(5) != bondBookmarks->end());
    std::list<Bond *> &res2 = bondBookmarks->at(5);
    std::list<Bond *> want2 = {bond1};
    CHECK(res2 == want2);
  }

  SECTION("RDMol/ROMol bookmarks work") {
    rdmol.setAtomBookmark(2, 1);
    CHECK(mol.getAtomWithBookmark(1) == atom2);
  }
}

TEST_CASE("ROMol Compat Mol properties") {
  std::string signedIntToken("signed int");
  std::string unusedToken("unused");
  std::string unsignedIntToken("unsigned int");
  std::string floatToken("float");
  RDMol rdmol;
  auto &mol = rdmol.asROMol();

  SECTION("Self basic properties") {
    mol.setProp(signedIntToken, 3);
    mol.setProp(unsignedIntToken, 3u);
    mol.setProp(floatToken, 3.0f);

    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    REQUIRE(mol.getPropIfPresent<int>(signedIntToken, res));
    CHECK(res == 3);
    REQUIRE(!mol.getPropIfPresent<int>(unusedToken, res));
  }

  SECTION("Prop handle updates in place") {
    mol.setProp(signedIntToken, 3);

    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    REQUIRE(mol.getPropIfPresent<int>(signedIntToken, res));
    REQUIRE(res == 3);

    mol.setProp(signedIntToken, 4);
    REQUIRE(mol.getPropIfPresent<int>(signedIntToken, res));
    CHECK(res == 4);
  }

  SECTION("Clear Mol Props") {
    mol.setProp(floatToken, 3.0f);
    mol.setProp(unsignedIntToken, 3u);

    mol.clearProp(unsignedIntToken);
    mol.clearProp(floatToken);

    float resf;
    unsigned int resu;
    CHECK(!mol.getPropIfPresent<unsigned int>(unsignedIntToken, resu));
    CHECK(!mol.getPropIfPresent<float>(floatToken, resf));
  }

  SECTION("Computed props and clear") {
    mol.setProp(floatToken, 3.0f, true);
    mol.setProp(unsignedIntToken, 3u);

    mol.clearComputedProps();
    float resf;
    unsigned int resu;
    CHECK(mol.getPropIfPresent<unsigned int>(unsignedIntToken, resu));
    CHECK(!mol.getPropIfPresent<float>(floatToken, resf));
  }

  SECTION("nonPOD Data") {
    mol.setProp(signedIntToken, boost::shared_array<int>(new int[3]{1, 2, 3}));

    boost::shared_array<int> res;
    REQUIRE(mol.getPropIfPresent(signedIntToken, res));

    CHECK(res[1] == 2);
  }

  SECTION("Has prop") {
    mol.setProp(floatToken, 3.0f);
    mol.setProp(signedIntToken, 5);

    CHECK(mol.hasProp(floatToken));
    CHECK(mol.hasProp(signedIntToken));
    CHECK(!mol.hasProp(unsignedIntToken));
  }

  SECTION("clear() RDProp API") {
    mol.setProp(floatToken, 3.0f);
    mol.clear();
    CHECK(!mol.hasProp(floatToken));
  }

  SECTION("UpdateProps RDProp no preserve") {
    RDProps props;
    props.setProp(floatToken, 3.0f);
    props.setProp(signedIntToken, 5);
    std::string vector2Token("vector2");
    const std::vector<int> vector2Value{6, 7};
    props.setProp(vector2Token, vector2Value);

    mol.setProp(floatToken, 4.0f);
    std::string vectorToken("vector");
    const std::vector<int> vectorValue{4, 5};
    mol.setProp(vectorToken, vectorValue);

    // Update from self should do nothing
    mol.updateProps(mol, /*preserve=*/false);
    // These values are never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    float res = 0.0f;
    int intres = 0;
    std::vector<int> vectorRes;
    REQUIRE(mol.getPropIfPresent(floatToken, res));
    CHECK(res == 4.0f);
    CHECK_FALSE(mol.getPropIfPresent(signedIntToken, intres));
    REQUIRE(mol.getPropIfPresent(vectorToken, vectorRes));
    CHECK(vectorRes == vectorValue);
    CHECK_FALSE(mol.getPropIfPresent(vector2Token, vectorRes));

    RDMol rdmol2;
    rdmol2.addAtom();
    rdmol2.addAtom();
    rdmol2.addBond(0, 1, BondEnums::BondType::SINGLE);
    PropToken atomPropName("atomProp");
    PropToken bondPropName("bondProp");
    rdmol2.addAtomProp(atomPropName, 2);
    rdmol2.addBondProp(bondPropName, 3);
    CHECK(rdmol2.getBondPropIfPresent(bondPropName, 0, intres));

    // Atom and bond props should not become mol props
    mol.updateProps(rdmol2.asROMol(), /*preserve=*/true);

    CHECK_FALSE(mol.getPropIfPresent(atomPropName.getString(), intres));
    CHECK_FALSE(mol.getPropIfPresent(bondPropName.getString(), intres));
    CHECK(rdmol2.getAtomPropIfPresent(atomPropName, 0, intres));
    CHECK(rdmol2.getBondPropIfPresent(bondPropName, 0, intres));
    CHECK_FALSE(
        rdmol2.asROMol().getPropIfPresent(atomPropName.getString(), intres));
    CHECK_FALSE(
        rdmol2.asROMol().getPropIfPresent(bondPropName.getString(), intres));

    // Update from other with preserve should add and replace props, but not
    // remove
    mol.updateProps(props, /*preserve=*/true);

    REQUIRE(mol.getPropIfPresent(floatToken, res));
    CHECK(res == 3.0f);
    REQUIRE(mol.getPropIfPresent(signedIntToken, intres));
    CHECK(intres == 5);
    REQUIRE(mol.getPropIfPresent(vectorToken, vectorRes));
    CHECK(vectorRes == vectorValue);
    REQUIRE(mol.getPropIfPresent(vector2Token, vectorRes));
    CHECK(vectorRes == vector2Value);

    mol.setProp(floatToken, 4.0f);
    REQUIRE(mol.getPropIfPresent(floatToken, res));
    CHECK(res == 4.0f);
    mol.clearProp(signedIntToken);
    CHECK_FALSE(mol.getPropIfPresent(signedIntToken, intres));
    mol.clearProp(vector2Token);
    CHECK_FALSE(mol.getPropIfPresent(vector2Token, intres));

    // Update from other with preserve should add and replace props, and remove
    mol.updateProps(props, /*preserve=*/false);

    REQUIRE(mol.getPropIfPresent(floatToken, res));
    CHECK(res == 3.0f);
    REQUIRE(mol.getPropIfPresent(signedIntToken, intres));
    CHECK(intres == 5);
    CHECK_FALSE(mol.getPropIfPresent(vectorToken, vectorRes));
    REQUIRE(mol.getPropIfPresent(vector2Token, vectorRes));
    CHECK(vectorRes == vector2Value);
  }

  SECTION("Get Prop list") {
    const std::string privateToken = "_private";
    mol.setProp(floatToken, 3.0f, false);
    mol.setProp(signedIntToken, 5, true);
    mol.setProp(privateToken, 5, false);

    CHECK_THAT(mol.getPropList(true, true),
               Catch::Matchers::UnorderedEquals(std::vector<std::string>(
                   {floatToken, signedIntToken, privateToken,
                    std::string(detail::computedPropName)})));
    CHECK_THAT(mol.getPropList(true, false),
               Catch::Matchers::UnorderedEquals(
                   std::vector<std::string>({floatToken, privateToken})));
    CHECK_THAT(mol.getPropList(false, true),
               Catch::Matchers::UnorderedEquals(std::vector<std::string>({
                   floatToken,
                   signedIntToken,
               })));
    CHECK_THAT(mol.getPropList(false, false),
               Catch::Matchers::UnorderedEquals(
                   std::vector<std::string>({floatToken})));
  }

  SECTION("Custom Data Copying") {
    mol.setProp(signedIntToken, MyCustomClass(7));
    size_t copyCount1 = MyCustomClass::copyCount;
    MyCustomClass copy = mol.getProp<MyCustomClass>(signedIntToken);
    size_t copyCount2 = MyCustomClass::copyCount;
    const MyCustomClass &a = mol.getProp<const MyCustomClass &>(signedIntToken);
    size_t copyCount3 = MyCustomClass::copyCount;
    const MyCustomClass &b = mol.getProp<const MyCustomClass &>(signedIntToken);
    size_t copyCount4 = MyCustomClass::copyCount;

    CHECK(copy.value == 7);
    CHECK(a.value == 7);
    CHECK(b.value == 7);
    CHECK(copyCount2 > copyCount1);
    CHECK(copyCount3 == copyCount2);
    CHECK(copyCount4 == copyCount3);
  }
}

template <typename T>
void addPropToAllAtoms(ROMol &mol, const std::string &token, const T &val,
                       bool computed = false) {
  for (auto atom : mol.atoms()) {
    atom->setProp<T>(token, val, computed);
  }
}

template <typename T>
void addPropToAllBonds(ROMol &mol, const std::string &token, const T &val,
                       bool computed = false) {
  for (auto bond : mol.bonds()) {
    bond->setProp<T>(token, val, computed);
  }
}

TEST_CASE("ROMol compat Atom/bond properties") {
  std::string signedIntToken("signed int");
  std::string unusedToken("unused");
  std::string unsignedIntToken("unsigned int");
  std::string floatToken("float");

  RDMol rdmol;
  RDKit::SmilesParseTemp temp;
  static const char *basicSmiles = "CCCC";

  SmilesParserParams params;
  params.sanitize = false;
  params.removeHs = false;

  REQUIRE(RDKit::SmilesToMol(basicSmiles, params, rdmol, temp) == true);
  REQUIRE(rdmol.getNumAtoms() == 4);
  REQUIRE(rdmol.getNumBonds() == 3);
  auto &mol = rdmol.asROMol();

  Atom *atom0 = mol.getAtomWithIdx(0);
  Atom *atom1 = mol.getAtomWithIdx(1);
  Atom *atom2 = mol.getAtomWithIdx(2);

  REQUIRE(atom0 != nullptr);
  REQUIRE(atom1 != nullptr);
  REQUIRE(atom2 != nullptr);

  Bond *bond0 = mol.getBondWithIdx(0);
  Bond *bond1 = mol.getBondWithIdx(1);

  REQUIRE(bond0 != nullptr);
  REQUIRE(bond1 != nullptr);

  SECTION("Atom properties access and override") {
    addPropToAllAtoms(mol, signedIntToken, 3);
    REQUIRE(atom0 != nullptr);
    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    REQUIRE(atom0->getPropIfPresent<int>(signedIntToken, res));
    CHECK(res == 3);
    atom0->setProp(signedIntToken, 4);
    REQUIRE(atom0->getPropIfPresent<int>(signedIntToken, res));
    CHECK(res == 4);
    REQUIRE(atom2->getPropIfPresent<int>(signedIntToken, res));
    CHECK(res == 3);
  }

  SECTION("Atom non-POD") {
    addPropToAllAtoms(mol, signedIntToken, std::vector<int>({1, 2, 3}));
    std::vector<int> res;
    REQUIRE(atom1->getPropIfPresent(signedIntToken, res));
    CHECK(res[1] == 2);
  }

  SECTION("Bond properties access and override") {
    addPropToAllBonds(mol, signedIntToken, 3);
    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    REQUIRE(bond1->getPropIfPresent<int>(signedIntToken, res));
    CHECK(res == 3);

    bond1->setProp(signedIntToken, 4);
    REQUIRE(bond1->getPropIfPresent<int>(signedIntToken, res));
    CHECK(res == 4);

    Bond *bond2 = mol.getBondWithIdx(2);
    REQUIRE(bond2 != nullptr);
    REQUIRE(bond2->getPropIfPresent<int>(signedIntToken, res));
    CHECK(res == 3);
  }

  SECTION("Bond non-POD") {
    addPropToAllBonds(mol, signedIntToken,
                      boost::shared_array<int>(new int[3]{1, 2, 3}));
    Bond *bond = mol.getBondWithIdx(0);
    boost::shared_array<int> res;
    REQUIRE(bond->getPropIfPresent(signedIntToken, res));
    CHECK(res[1] == 2);
  }

  SECTION("Clear atom properties") {
    addPropToAllAtoms(mol, floatToken, 3.0f);
    addPropToAllAtoms(mol, signedIntToken, 5);

    atom0->clearProp(floatToken);

    // Check that we've cleared atom 0
    float res;
    int intres;
    CHECK(!atom0->getPropIfPresent<float>(floatToken, res));

    // Check that we haven't cleared atom 1
    REQUIRE(atom1->getPropIfPresent<float>(floatToken, res));

    // Check that we haven't cleared the int prop
    REQUIRE(atom0->getPropIfPresent<int>(signedIntToken, intres));
  }

  SECTION("Clear bond properties") {
    addPropToAllBonds(mol, floatToken, 3.0f);
    addPropToAllBonds(mol, signedIntToken, 5);

    bond0->clearProp(floatToken);

    // Check that we've cleared atom 0
    float res;
    int intres;
    CHECK(!bond0->getPropIfPresent<float>(floatToken, res));

    // Check that we haven't cleared atom 1
    REQUIRE(bond1->getPropIfPresent<float>(floatToken, res));

    // Check that we haven't cleared the int prop
    REQUIRE(bond0->getPropIfPresent<int>(signedIntToken, intres));
  }

  SECTION("Clear computed also clears atoms and bonds") {
    addPropToAllAtoms(mol, floatToken, 3.0f, true);
    addPropToAllBonds(mol, signedIntToken, 5, true);
    addPropToAllBonds(mol, floatToken, 4.0f, false);

    mol.clearComputedProps();

    float res;
    int intres;
    CHECK(!atom0->getPropIfPresent<float>(floatToken, res));
    CHECK(!bond0->getPropIfPresent<int>(signedIntToken, intres));
    // Non-computed should remain.
    CHECK(bond0->getPropIfPresent<float>(floatToken, res));
  }

  SECTION("Has prop") {
    addPropToAllAtoms(mol, floatToken, 3.0f, true);
    addPropToAllBonds(mol, signedIntToken, 5, true);

    CHECK(atom0->hasProp(floatToken));
    CHECK(bond0->hasProp(signedIntToken));
    CHECK(!atom0->hasProp(signedIntToken));
  }

  SECTION("clear() RDProp API") {
    addPropToAllAtoms(mol, floatToken, 3.0f, true);
    atom0->clear();
    CHECK(!atom0->hasProp(floatToken));
  }

  SECTION("Store different types with same key") {
    atom1->setProp(floatToken, 5u);
    atom0->setProp(floatToken, -3.0f);

    bond0->setProp(floatToken, 4.0f);
    std::vector<int> storeVal = {1, 2, 3};
    bond1->setProp(floatToken, storeVal);

    // These values are never used, but GCC can no longer deduce that they're
    // initialized, and maybe-uninitialized is currently treated as an error.
    float resf = 0.0f;
    unsigned int resu = 0;
    REQUIRE(atom0->getPropIfPresent<float>(floatToken, resf));
    CHECK(resf == -3.0f);
    REQUIRE(atom1->getPropIfPresent<unsigned int>(floatToken, resu));
    CHECK(resu == 5u);

    std::vector<int> resv;
    REQUIRE(bond1->getPropIfPresent(std::string(floatToken), resv));
    CHECK(resv == storeVal);
    REQUIRE(bond0->getPropIfPresent<float>(floatToken, resf));
    CHECK(resf == 4.0f);
  }

  SECTION("Clear all") {
    atom0->setProp(floatToken, 3.0f);
    atom1->setProp(signedIntToken, 5);
    bond0->setProp(floatToken, 4.0f);
    bond1->setProp(signedIntToken, 6);
    mol.setProp(unsignedIntToken, 7);

    mol.clear();
    CHECK(!atom0->hasProp(floatToken));
    CHECK(!atom1->hasProp(signedIntToken));
    CHECK(!bond0->hasProp(floatToken));
    CHECK(!bond1->hasProp(signedIntToken));
    CHECK(!mol.hasProp(unsignedIntToken));
  }

  SECTION("Custom Atom Data Copying") {
    std::string testPropName("test");
    atom0->clearProp(testPropName);
    atom1->clearProp(testPropName);
    atom0->setProp(testPropName, MyCustomClass(11));
    atom1->setProp(testPropName, MyCustomClass(17));
    size_t copyCount1 = MyCustomClass::copyCount;
    MyCustomClass copy = atom0->getProp<MyCustomClass>(testPropName);
    size_t copyCount2 = MyCustomClass::copyCount;
    const MyCustomClass &a0 =
        atom0->getProp<const MyCustomClass &>(testPropName);
    size_t copyCount3 = MyCustomClass::copyCount;
    const MyCustomClass &b0 =
        atom0->getProp<const MyCustomClass &>(testPropName);
    size_t copyCount4 = MyCustomClass::copyCount;

    CHECK(copy.value == 11);
    CHECK(a0.value == 11);
    CHECK(b0.value == 11);
    CHECK(copyCount2 > copyCount1);
    CHECK(copyCount3 == copyCount2);
    CHECK(copyCount4 == copyCount3);

    copy = atom1->getProp<MyCustomClass>(testPropName);
    size_t copyCount5 = MyCustomClass::copyCount;
    const MyCustomClass &a1 =
        atom1->getProp<const MyCustomClass &>(testPropName);
    size_t copyCount6 = MyCustomClass::copyCount;
    const MyCustomClass &b1 =
        atom1->getProp<const MyCustomClass &>(testPropName);
    size_t copyCount7 = MyCustomClass::copyCount;

    CHECK(copy.value == 17);
    CHECK(a1.value == 17);
    CHECK(b1.value == 17);
    CHECK(copyCount5 > copyCount4);
    CHECK(copyCount6 == copyCount5);
    CHECK(copyCount7 == copyCount6);
  }

  SECTION("Custom Bond Data Copying") {
    std::string testPropName("testb");
    bond0->setProp(testPropName, MyCustomClass(12));
    bond1->setProp(testPropName, MyCustomClass(18));
    size_t copyCount1 = MyCustomClass::copyCount;
    MyCustomClass copy = bond0->getProp<MyCustomClass>(testPropName);

    size_t copyCount2 = MyCustomClass::copyCount;
    const MyCustomClass &a0 =
        bond0->getProp<const MyCustomClass &>(testPropName);
    size_t copyCount3 = MyCustomClass::copyCount;
    const MyCustomClass &b0 =
        bond0->getProp<const MyCustomClass &>(testPropName);
    size_t copyCount4 = MyCustomClass::copyCount;

    CHECK(copy.value == 12);
    CHECK(a0.value == 12);
    CHECK(b0.value == 12);
    CHECK(copyCount2 > copyCount1);
    CHECK(copyCount3 == copyCount2);
    CHECK(copyCount4 == copyCount3);

    copy = bond1->getProp<MyCustomClass>(testPropName);
    size_t copyCount5 = MyCustomClass::copyCount;
    const MyCustomClass &a1 =
        bond1->getProp<const MyCustomClass &>(testPropName);
    size_t copyCount6 = MyCustomClass::copyCount;
    const MyCustomClass &b1 =
        bond1->getProp<const MyCustomClass &>(testPropName);
    size_t copyCount7 = MyCustomClass::copyCount;

    CHECK(copy.value == 18);
    CHECK(a1.value == 18);
    CHECK(b1.value == 18);
    CHECK(copyCount5 > copyCount4);
    CHECK(copyCount6 == copyCount5);
    CHECK(copyCount7 == copyCount6);
  }
}

template <typename T1, typename T2>
void setPropsForUpdatePropsTest(T1 &prop1, T2 &prop2) {
  // basics
  prop1.clear();
  prop2.clear();
  prop1.setProp("first", 3);
  prop1.setProp("second", "value");
  prop2.setProp("first", 5.0);
  prop2.setProp("third", 7);

  // computed behavior
  // First, one that's only in source
  prop1.setProp("firstComputed", 4, true);
  // Next, one that's computed in both.
  prop1.setProp("computedInBoth", "value1", true);
  prop2.setProp("computedInBoth", "value2", true);
  // Next, computed only in source
  prop1.setProp("computedOnlyIn1", 1, true);
  prop2.setProp("computedOnlyIn1", 2, false);
  // Finally, computed only in dest
  prop1.setProp("computedOnlyIn2", 1, false);
  prop2.setProp("computedOnlyIn2", 2, true);
}

template <typename T>
void checkUpdatePropsValuesNoPreserve(const T &prop2,
                                      const std::string &message) {
  INFO(message);
  CHECK(prop2.template getProp<int>("first") == 3);
  CHECK(prop2.template getProp<std::string>("second") == "value");
  CHECK(!prop2.hasProp("third"));

  CHECK(prop2.template getProp<int>("firstComputed") == 4);
  CHECK(prop2.template getProp<std::string>("computedInBoth") == "value1");
  CHECK(prop2.template getProp<int>("computedOnlyIn1") == 1);
  CHECK(prop2.template getProp<int>("computedOnlyIn2") == 1);

  // Get a list of the computed props
  std::vector<std::string> computedProps = prop2.getPropList(false, true);
  std::vector<std::string> nonComputedProps = prop2.getPropList(false, false);

  for (const auto &prop : nonComputedProps) {
    computedProps.erase(
        std::remove(computedProps.begin(), computedProps.end(), prop),
        computedProps.end());
  }
  // Handle the private __computedProps implementation of RDProps
  computedProps.erase(std::remove(computedProps.begin(), computedProps.end(),
                                  std::string("__computedProps")),
                      computedProps.end());

  std::vector<std::string> wantComputedProps = {
      "firstComputed", "computedInBoth", "computedOnlyIn1"};
  CHECK_THAT(computedProps,
             Catch::Matchers::UnorderedEquals(wantComputedProps));
}

template <typename T>
void checkUpdatePropsValuesPreserve(const T &prop2,
                                    const std::string &message) {
  INFO(message);
  CHECK(prop2.template getProp<int>("first") == 3);
  CHECK(prop2.template getProp<std::string>("second") == "value");
  CHECK(prop2.template getProp<int>("third") == 7);

  CHECK(prop2.template getProp<int>("firstComputed") == 4);
  CHECK(prop2.template getProp<std::string>("computedInBoth") == "value1");
  CHECK(prop2.template getProp<int>("computedOnlyIn1") == 1);
  CHECK(prop2.template getProp<int>("computedOnlyIn2") == 1);

  // Get computed props
  std::vector<std::string> computedProps = prop2.getPropList(false, true);
  std::vector<std::string> nonComputedProps = prop2.getPropList(false, false);

  for (const auto &prop : nonComputedProps) {
    computedProps.erase(
        std::remove(computedProps.begin(), computedProps.end(), prop),
        computedProps.end());
  }
  // Handle the private __computedProps implementation of RDProps
  computedProps.erase(std::remove(computedProps.begin(), computedProps.end(),
                                  std::string("__computedProps")),
                      computedProps.end());

  std::vector<std::string> wantComputedProps = {
      "firstComputed", "computedInBoth", "computedOnlyIn1"};
  CHECK_THAT(computedProps,
             Catch::Matchers::UnorderedEquals(wantComputedProps));
}

template <typename T1, typename T2>
void runFullUpdatePropsSequence(T1 &prop1, T2 &prop2,
                                const std::string &message) {
  setPropsForUpdatePropsTest(prop1, prop2);
  prop2.updateProps(prop1, /*preserve=*/false);
  checkUpdatePropsValuesNoPreserve(prop2, message + ": no preserve");

  setPropsForUpdatePropsTest(prop1, prop2);
  prop2.updateProps(prop1, /*preserve=*/true);
  checkUpdatePropsValuesPreserve(prop2, message + ": preserve");
}

TEST_CASE("Atom/bond updateprops") {
  auto mol = basicMol();
  ROMol &romol = mol->asROMol();
  auto mol2 = basicMol();
  ROMol &romol2 = mol2->asROMol();
  Atom soloAtom(6);
  Atom soloAtom2(8);

  Bond soloBond;
  Bond soloBond2;

  SECTION("Check RDProps behavior assumptions") {
    RDProps prop1, prop2;
    runFullUpdatePropsSequence(prop1, prop2, "RDProps");
  }

  SECTION("Solo Atom to solo Atom") {
    runFullUpdatePropsSequence(soloAtom, soloAtom2, "Solo atom to solo atom");
  }

  SECTION("Mol Atom to mol Atom") {
    runFullUpdatePropsSequence(*romol.getAtomWithIdx(1),
                               *romol2.getAtomWithIdx(2),
                               "mol atom to mol atom");
  }

  SECTION("Mol Atom to non-mol Atom") {
    runFullUpdatePropsSequence(*romol.getAtomWithIdx(1), soloAtom,
                               "mol atom to non-mol atom");
  }

  SECTION("Non-Atom to Atom") {
    Conformer conf;
    runFullUpdatePropsSequence(conf, *romol.getAtomWithIdx(1),
                               "non-atom to atom");
  }
  // Note: Bond to generic non-Bond not supported, as the inherited APIs take
  // in an RDProps.

  SECTION("Solo Bond to solo Bond") {
    runFullUpdatePropsSequence(soloBond, soloBond2, "solo bond to solo bond");
  }

  SECTION("Mol Bond to mol Bond") {
    runFullUpdatePropsSequence(*romol.getBondWithIdx(1),
                               *romol2.getBondWithIdx(2),
                               "mol bond to mol bond");
  }

  SECTION("Mol Bond to non-mol Bond") {
    runFullUpdatePropsSequence(*romol.getBondWithIdx(1), soloBond,
                               "mol bond to nonmol bond");
  }

  SECTION("Non-Bond to Bond") {
    Conformer conf;
    runFullUpdatePropsSequence(conf, *romol.getBondWithIdx(1),
                               "nonbond to bond");
  }

  SECTION("Mol to Mol") {
    runFullUpdatePropsSequence(romol, romol2, "mol to mol");
  }

  SECTION("Non-mol to Mol") {
    Conformer conf;
    runFullUpdatePropsSequence(conf, romol, "non mol to mol");
  }
}

TEST_CASE("Atom Monomomer info") {
  auto mol = basicMol();
  ROMol &romol = mol->asROMol();
  Atom *atom = romol.getAtomWithIdx(0);
  Atom *atom2 = romol.getAtomWithIdx(1);
  const Atom *atom3 =
      romol.getAtomWithIdx(2);  // Used to check const overload of null entry
  Atom *atom4 = romol.getAtomWithIdx(
      3);  // Used to check non-const overload of null entry

  SECTION("Has None") {
    CHECK(atom->getMonomerInfo() == nullptr);
    CHECK(atom2->getMonomerInfo() == nullptr);
  }

  SECTION("Set/get") {
    auto *baseClassInfo = new AtomMonomerInfo(AtomMonomerInfo::OTHER);
    atom->setMonomerInfo(baseClassInfo);
    auto *PDBInfo = new AtomPDBResidueInfo("dummy", /*serial=*/3);
    atom2->setMonomerInfo(PDBInfo);

    CHECK(atom4->getMonomerInfo() == nullptr);
    AtomMonomerInfo *res1 = atom->getMonomerInfo();
    const AtomMonomerInfo *res2 = atom2->getMonomerInfo();

    CHECK(res1->getMonomerType() == AtomMonomerInfo::OTHER);
    CHECK(res2->getMonomerType() == AtomMonomerInfo::PDBRESIDUE);

    auto *res2PDB = dynamic_cast<const AtomPDBResidueInfo *>(res2);
    CHECK(res2PDB != nullptr);
    CHECK(res2PDB->getSerialNumber() == 3);

    // Const overload
    const Atom *atomConst = atom;
    const AtomMonomerInfo *res1Const = atomConst->getMonomerInfo();
    CHECK(res1Const->getMonomerType() == AtomMonomerInfo::OTHER);
    CHECK(atom3->getMonomerInfo() == nullptr);
  }

  SECTION("Override") {
    auto *baseClassInfo = new AtomMonomerInfo(AtomMonomerInfo::OTHER);
    atom->setMonomerInfo(baseClassInfo);
    auto *PDBInfo = new AtomPDBResidueInfo("dummy", /*serial=*/3);
    atom->setMonomerInfo(PDBInfo);

    AtomMonomerInfo *res1 = atom->getMonomerInfo();
    CHECK(res1->getMonomerType() == RDKit::AtomMonomerInfo::PDBRESIDUE);
  }
}

TEST_CASE("Atom constructors/ownership") {
  SECTION("Default") {
    Atom atom;
    CHECK(!atom.hasOwningMol());

    Atom atom2(6);
    CHECK(!atom2.hasOwningMol());
    CHECK(atom2.getAtomicNum() == 6);

    Atom atom3("C");
    CHECK(!atom3.hasOwningMol());
    CHECK(atom3.getAtomicNum() == 6);
  }

  SECTION("Copy, no previous owner") {
    Atom atom(6);
    atom.setProp("test", 3, /*computed=*/true);
    atom.setProp("complexTest", std::vector<int>({1, 2, 3}));
    atom.setMonomerInfo(new AtomMonomerInfo(AtomMonomerInfo::OTHER));

    Atom atom2(atom);
    CHECK(!atom2.hasOwningMol());
    CHECK(atom2.getAtomicNum() == 6);

    Atom atom3 = atom;
    CHECK(!atom3.hasOwningMol());
    CHECK(atom3.getAtomicNum() == 6);

    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    REQUIRE(atom3.getPropIfPresent<int>("test", res));
    CHECK(res == 3);
    std::vector<int> resv;
    REQUIRE(atom3.getPropIfPresent("complexTest", resv));
    CHECK(resv == std::vector<int>({1, 2, 3}));

    AtomMonomerInfo *info = atom3.getMonomerInfo();
    REQUIRE(info != nullptr);
    CHECK(info->getMonomerType() == AtomMonomerInfo::OTHER);
  }

  SECTION("Copy, previous owner") {
    auto molPtr = basicMol();
    RDMol &mol = *molPtr;
    auto &romol = mol.asROMol();

    Atom *atom = romol.getAtomWithIdx(0);
    atom->setProp("test", 3, /*computed=*/true);
    atom->setProp("complexTest", std::vector<int>({1, 2, 3}));
    atom->setMonomerInfo(new AtomMonomerInfo(AtomMonomerInfo::OTHER));

    REQUIRE(atom != nullptr);
    REQUIRE(atom->getAtomicNum() == 6);
    REQUIRE(&atom->getOwningMol() == &romol);

    Atom atom2(*atom);
    CHECK(!atom2.hasOwningMol());
    CHECK(atom2.getAtomicNum() == 6);

    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    REQUIRE(atom2.getPropIfPresent<int>("test", res));
    CHECK(res == 3);
    std::vector<int> resv;
    REQUIRE(atom2.getPropIfPresent("complexTest", resv));
    CHECK(resv == std::vector<int>({1, 2, 3}));

    AtomMonomerInfo *info = atom2.getMonomerInfo();
    REQUIRE(info != nullptr);
    CHECK(info->getMonomerType() == AtomMonomerInfo::OTHER);
  }

  SECTION("move, no previous owner") {
    Atom atom(6);
    atom.setProp("test", 3, /*computed=*/true);
    atom.setProp("complexTest", std::vector<int>({1, 2, 3}));
    atom.setMonomerInfo(new AtomMonomerInfo(AtomMonomerInfo::OTHER));
    Atom atom2 = std::move(atom);
    CHECK(!atom2.hasOwningMol());
    CHECK(atom2.getAtomicNum() == 6);

    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    REQUIRE(atom2.getPropIfPresent<int>("test", res));
    CHECK(res == 3);
    std::vector<int> resv;
    REQUIRE(atom2.getPropIfPresent("complexTest", resv));
    CHECK(resv == std::vector<int>({1, 2, 3}));

    AtomMonomerInfo *info = atom2.getMonomerInfo();
    REQUIRE(info != nullptr);
    CHECK(info->getMonomerType() == AtomMonomerInfo::OTHER);
  }

  SECTION("Move, previous owner") {
    auto molPtr = basicMol();
    RDMol &mol = *molPtr;
    auto &romol = mol.asROMol();

    Atom *atom = romol.getAtomWithIdx(0);
    atom->setProp("test", 3, /*computed=*/true);
    atom->setProp("complexTest", std::vector<int>({1, 2, 3}));
    atom->setMonomerInfo(new AtomMonomerInfo(AtomMonomerInfo::OTHER));

    REQUIRE(atom != nullptr);
    REQUIRE(atom->getAtomicNum() == 6);
    REQUIRE(&atom->getOwningMol() == &romol);

    Atom atom2(std::move(*atom));
    CHECK(atom2.hasOwningMol());
    CHECK(atom2.getAtomicNum() == 6);

    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    REQUIRE(atom2.getPropIfPresent<int>("test", res));
    CHECK(res == 3);
    std::vector<int> resv;
    REQUIRE(atom2.getPropIfPresent("complexTest", resv));
    CHECK(resv == std::vector<int>({1, 2, 3}));

    AtomMonomerInfo *info = atom2.getMonomerInfo();
    REQUIRE(info != nullptr);
    CHECK(info->getMonomerType() == AtomMonomerInfo::OTHER);
  }
}

TEST_CASE("Bond constructors/ownership") {
  SECTION("Default") {
    Bond bond;
    CHECK(!bond.hasOwningMol());

    Bond bond2(Bond::BondType::SINGLE);
    CHECK(!bond2.hasOwningMol());
    CHECK(bond2.getBondType() == Bond::BondType::SINGLE);
  }

  SECTION("Constructors, no previous owner") {
    Bond bond(Bond::BondType::SINGLE);
    bond.setProp("test", 3, /*computed=*/true);
    bond.setProp("complexTest", std::vector<int>({1, 2, 3}));
    Bond bond2(bond);
    CHECK(!bond2.hasOwningMol());

    Bond bond3(std::move(bond));
    CHECK(!bond3.hasOwningMol());

    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    REQUIRE(bond3.getPropIfPresent<int>("test", res));
    CHECK(res == 3);
    std::vector<int> resv;
    REQUIRE(bond3.getPropIfPresent("complexTest", resv));
    CHECK(resv == std::vector<int>({1, 2, 3}));
  }

  SECTION("Copy assignment operator, previous owner") {
    auto molPtr = basicMol();
    RDMol &mol = *molPtr;
    auto &romol = mol.asROMol();

    Bond *bond = romol.getBondWithIdx(0);
    bond->setProp("test", 3, /*computed=*/true);
    bond->setProp("complexTest", std::vector<int>({1, 2, 3}));
    REQUIRE(bond != nullptr);
    REQUIRE(bond->getBondType() == Bond::BondType::SINGLE);
    REQUIRE(&bond->getOwningMol() == &romol);

    Bond bond2 = *bond;
    CHECK(!bond2.hasOwningMol());
    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    REQUIRE(bond2.getPropIfPresent<int>("test", res));
    CHECK(res == 3);
    std::vector<int> resv;
    REQUIRE(bond2.getPropIfPresent("complexTest", resv));
    CHECK(resv == std::vector<int>({1, 2, 3}));
  }

  SECTION("Move assign operator, no previous owner") {
    Bond bond(Bond::BondType::SINGLE);
    bond.setProp("test", 3, /*computed=*/true);
    bond.setProp("complexTest", std::vector<int>({1, 2, 3}));
    Bond bond2 = std::move(bond);
    CHECK(!bond2.hasOwningMol());
    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    REQUIRE(bond2.getPropIfPresent<int>("test", res));
    CHECK(res == 3);
    std::vector<int> resv;
    REQUIRE(bond2.getPropIfPresent("complexTest", resv));
    CHECK(resv == std::vector<int>({1, 2, 3}));
  }

  SECTION("Move constructor, previous owner") {
    auto molPtr = basicMol();
    RDMol &mol = *molPtr;
    auto &romol = mol.asROMol();

    Bond *bond = romol.getBondWithIdx(0);
    bond->setProp("test", 3, /*computed=*/true);
    bond->setProp("complexTest", std::vector<int>({1, 2, 3}));
    REQUIRE(bond != nullptr);
    REQUIRE(bond->getBondType() == Bond::BondType::SINGLE);
    REQUIRE(&bond->getOwningMol() == &romol);

    Bond bond2(std::move(*bond));
    CHECK(bond2.hasOwningMol());
    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    REQUIRE(bond2.getPropIfPresent<int>("test", res));
    CHECK(res == 3);
    std::vector<int> resv;
    REQUIRE(bond2.getPropIfPresent("complexTest", resv));
    CHECK(resv == std::vector<int>({1, 2, 3}));
  }
}

TEST_CASE("Atom map number") {
  auto molPtr = basicMol();
  auto &mol = molPtr->asROMol();
  Atom *atom0 = mol.getAtomWithIdx(0);
  Atom *atom1 = mol.getAtomWithIdx(1);

  SECTION("Set and get atom map number") {
    atom0->setAtomMapNum(1);
    atom1->setAtomMapNum(2);
    CHECK(atom0->getAtomMapNum() == 1);
    CHECK(atom1->getAtomMapNum() == 2);
  }

  SECTION("Strict/non-strict") {
    atom0->setAtomMapNum(-10, false);
    CHECK(atom0->getAtomMapNum() == -10);
    CHECK_THROWS_AS(atom1->setAtomMapNum(-10, true), Invar::Invariant);
  }

  SECTION("Default") { CHECK(atom0->getAtomMapNum() == 0); }

  SECTION("Set as 0 clears") {
    atom0->setAtomMapNum(1);
    atom0->setAtomMapNum(0);
    atom1->setAtomMapNum(0);
    CHECK(atom0->getAtomMapNum() == 0);
    CHECK(atom1->getAtomMapNum() == 0);
  }
}

// Test cases taken from comment above function implementation.
TEST_CASE("GetPerturbationOrder") {
  const std::string fourBondedCarbon = "C(N)(N)(N)N";
  auto mol = parseSmiles(fourBondedCarbon.c_str());
  auto &romol = mol->asROMol();
  Atom *atom0 = romol.getAtomWithIdx(0);

  INT_LIST inputOrder;
  inputOrder.push_back(1);
  inputOrder.push_back(0);
  inputOrder.push_back(2);
  inputOrder.push_back(3);
  CHECK(atom0->getPerturbationOrder(inputOrder) == 1);

  inputOrder.clear();
  inputOrder.push_back(1);
  inputOrder.push_back(2);
  inputOrder.push_back(3);
  inputOrder.push_back(0);
  CHECK(atom0->getPerturbationOrder(inputOrder) == 3);

  inputOrder.clear();
  inputOrder.push_back(1);
  inputOrder.push_back(2);
  inputOrder.push_back(0);
  inputOrder.push_back(3);
  CHECK(atom0->getPerturbationOrder(inputOrder) == 2);
}

TEST_CASE("Atom valence conversions from RDMol") {
  auto molPtr = basicMol();
  auto &mol = molPtr->asROMol();
  Atom *atom0 = mol.getAtomWithIdx(0);
  Atom *atom1 = mol.getAtomWithIdx(1);
  Atom *atom2 = mol.getAtomWithIdx(2);
  Atom *atom3 = mol.getAtomWithIdx(3);

  SECTION("Neither implicit nor explicit set") {
    CHECK_THROWS_AS(atom0->getValence(Atom::ValenceType::EXPLICIT),
                    Invar::Invariant);
    CHECK_THROWS_AS(atom0->getValence(Atom::ValenceType::IMPLICIT),
                    Invar::Invariant);
  }

  SECTION("Explicit set") {
    atom0->calcExplicitValence(false);
    CHECK(atom0->getValence(Atom::ValenceType::EXPLICIT) == 1);
    CHECK_THROWS_AS(atom0->getValence(Atom::ValenceType::IMPLICIT),
                    Invar::Invariant);
  }

  SECTION("No implict set") {
    atom0->setNoImplicit(true);
    CHECK(atom0->getValence(Atom::ValenceType::IMPLICIT) == 0);
  }

  SECTION("Correct after calculation") {
    atom0->calcImplicitValence(false);
    atom1->calcImplicitValence(false);
    atom2->calcImplicitValence(false);
    atom3->calcImplicitValence(false);

    CHECK(atom0->getValence(Atom::ValenceType::EXPLICIT) == 1);
    CHECK(atom1->getValence(Atom::ValenceType::EXPLICIT) == 2);
    // double bond
    CHECK(atom2->getValence(Atom::ValenceType::EXPLICIT) == 3);
    CHECK(atom3->getValence(Atom::ValenceType::EXPLICIT) == 2);

    CHECK(atom0->getValence(Atom::ValenceType::IMPLICIT) == 3);
    CHECK(atom1->getValence(Atom::ValenceType::IMPLICIT) == 2);
    CHECK(atom2->getValence(Atom::ValenceType::IMPLICIT) == 1);
    CHECK(atom3->getValence(Atom::ValenceType::IMPLICIT) == 2);
  }

  SECTION("Atom with no bonds but explicit hydrogens") {
    auto singleAtomMol = parseSmiles("C", false, false);
    auto &romol = singleAtomMol->asROMol();
    auto *atom = romol.getAtomWithIdx(0);
    atom->setNumExplicitHs(3);
    atom->calcExplicitValence(false);
    CHECK(atom->getValence(Atom::ValenceType::EXPLICIT) == 3);
  }
}

TEST_CASE("Bond getOtherAtom and getOtherAtomIdx") {
  auto molPtr = basicMol();
  auto &mol = molPtr->asROMol();
  Bond *bond1 = mol.getBondWithIdx(1);
  REQUIRE(bond1 != nullptr);
  REQUIRE(bond1->getBeginAtomIdx() == 1);
  REQUIRE(bond1->getEndAtomIdx() == 2);

  SECTION("getOtherAtom") {
    Atom *atom1 = bond1->getBeginAtom();
    Atom *atom2 = bond1->getEndAtom();
    Atom *other1 = bond1->getOtherAtom(atom1);
    Atom *other2 = bond1->getOtherAtom(atom2);
    CHECK(other1 == atom2);
    CHECK(other2 == atom1);
  }

  SECTION("getOtherAtomFails") {
    CHECK_THROWS_AS(bond1->getOtherAtom(nullptr), Invar::Invariant);

    // getOtherAtom needs an owning molecule
    std::unique_ptr<Bond> isolatedBond(new Bond(RDKit::Bond::BondType::SINGLE));
    isolatedBond->setBeginAtomIdx(bond1->getBeginAtomIdx());
    isolatedBond->setEndAtomIdx(bond1->getEndAtomIdx());
    CHECK_THROWS_AS(isolatedBond->getOtherAtom(bond1->getBeginAtom()),
                    Invar::Invariant);
  }

  SECTION("getOtherAtomIdx") {
    CHECK(bond1->getOtherAtomIdx(1) == 2);
    CHECK(bond1->getOtherAtomIdx(2) == 1);
  }

  SECTION("GetOtherAtomIdxFails") {
    CHECK_THROWS_AS(bond1->getOtherAtomIdx(3), Invar::Invariant);
  }
}

TEST_CASE("Update property cache") {
  auto molPtr = basicMol();
  auto &mol = molPtr->asROMol();
  Atom *atom0 = mol.getAtomWithIdx(0);
  Atom *atom1 = mol.getAtomWithIdx(1);

  SECTION("Needs update - explicit") {
    CHECK(atom0->needsUpdatePropertyCache() == true);
  }

  SECTION("Needs update - implicit") {
    atom1->calcExplicitValence(false);
    CHECK(atom1->needsUpdatePropertyCache() == true);
  }

  SECTION("Does not need - noimplict") {
    atom0->setNoImplicit(true);
    atom0->calcExplicitValence(false);
    CHECK(atom0->needsUpdatePropertyCache() == false);
  }

  SECTION("Update cache") {
    CHECK(atom0->needsUpdatePropertyCache() == true);
    mol.updatePropertyCache();
    CHECK(atom0->needsUpdatePropertyCache() == false);
  }
}

TEST_CASE("Bond Stereo Atoms") {
  RDMol rdmol;
  RDKit::SmilesParseTemp temp;
  static const char *basicSmiles = "CC1C(C)C1";

  SmilesParserParams params;
  params.sanitize = false;
  params.removeHs = false;

  REQUIRE(RDKit::SmilesToMol(basicSmiles, params, rdmol, temp) == true);
  REQUIRE(rdmol.getNumAtoms() == 5);
  REQUIRE(rdmol.getNumBonds() == 5);

  SECTION("Sync various combinations") {
    // Start with only RDMol
    REQUIRE(rdmol.getBondIndexBetweenAtoms(1, 2) == 1);
    REQUIRE(rdmol.getBondIndexBetweenAtoms(1, 4) < 5);
    REQUIRE(rdmol.getBondIndexBetweenAtoms(0, 1) < 5);
    REQUIRE(rdmol.getBondIndexBetweenAtoms(2, 3) < 5);
    REQUIRE(rdmol.getBondIndexBetweenAtoms(2, 4) < 5);
    {
      const atomindex_t *newStereoAtoms = rdmol.getBondStereoAtoms(1);
      REQUIRE(newStereoAtoms != nullptr);
      CHECK(newStereoAtoms[0] == atomindex_t(-1));
      CHECK(newStereoAtoms[1] == atomindex_t(-1));
    }
    rdmol.setBondStereoAtoms(1, 4, 4);
    {
      const atomindex_t *newStereoAtoms = rdmol.getBondStereoAtoms(1);
      REQUIRE(newStereoAtoms != nullptr);
      CHECK(newStereoAtoms[0] == 4);
      CHECK(newStereoAtoms[1] == 4);
    }
    rdmol.clearBondStereoAtoms(1);
    {
      const atomindex_t *newStereoAtoms = rdmol.getBondStereoAtoms(1);
      REQUIRE(newStereoAtoms != nullptr);
      CHECK(newStereoAtoms[0] == atomindex_t(-1));
      CHECK(newStereoAtoms[1] == atomindex_t(-1));
    }

    auto &mol = rdmol.asROMol();

    Bond *bond1 = mol.getBondWithIdx(1);
    REQUIRE(mol.getBondBetweenAtoms(1, 2) == bond1);
    REQUIRE(mol.getBondBetweenAtoms(1, 4) != nullptr);
    REQUIRE(mol.getBondBetweenAtoms(0, 1) != nullptr);
    REQUIRE(mol.getBondBetweenAtoms(2, 3) != nullptr);
    REQUIRE(mol.getBondBetweenAtoms(2, 4) != nullptr);
    {
      const RDKit::INT_VECT &stereoAtoms = bond1->getStereoAtoms();
      REQUIRE(stereoAtoms.size() == 0);
      const atomindex_t *newStereoAtoms = rdmol.getBondStereoAtoms(1);
      REQUIRE(newStereoAtoms != nullptr);
      CHECK(newStereoAtoms[0] == atomindex_t(-1));
      CHECK(newStereoAtoms[1] == atomindex_t(-1));
    }
    bond1->setStereoAtoms(0, 4);
    {
      const RDKit::INT_VECT &stereoAtoms = bond1->getStereoAtoms();
      REQUIRE(stereoAtoms.size() == 2);
      CHECK(stereoAtoms[0] == 0);
      CHECK(stereoAtoms[1] == 4);
      const atomindex_t *newStereoAtoms = rdmol.getBondStereoAtoms(1);
      REQUIRE(newStereoAtoms != nullptr);
      CHECK(newStereoAtoms[0] == 0);
      CHECK(newStereoAtoms[1] == 4);
    }
    rdmol.setBondStereoAtoms(1, 4, 3);
    {
      const RDKit::INT_VECT &stereoAtoms = bond1->getStereoAtoms();
      REQUIRE(stereoAtoms.size() == 2);
      CHECK(stereoAtoms[0] == 4);
      CHECK(stereoAtoms[1] == 3);
      const atomindex_t *newStereoAtoms = rdmol.getBondStereoAtoms(1);
      REQUIRE(newStereoAtoms != nullptr);
      CHECK(newStereoAtoms[0] == 4);
      CHECK(newStereoAtoms[1] == 3);
    }
    bond1->setStereoAtoms(0, 3);
    {
      const RDKit::INT_VECT &stereoAtoms = bond1->getStereoAtoms();
      REQUIRE(stereoAtoms.size() == 2);
      CHECK(stereoAtoms[0] == 0);
      CHECK(stereoAtoms[1] == 3);
      const atomindex_t *newStereoAtoms = rdmol.getBondStereoAtoms(1);
      REQUIRE(newStereoAtoms != nullptr);
      CHECK(newStereoAtoms[0] == 0);
      CHECK(newStereoAtoms[1] == 3);
    }
  }
}

TEST_CASE("Stereo Groups") {
  RDMol rdmol;
  RDKit::SmilesParseTemp temp;
  static const char *basicSmiles = "CCCCCC";

  SmilesParserParams params;
  params.sanitize = false;
  params.removeHs = false;

  REQUIRE(RDKit::SmilesToMol(basicSmiles, params, rdmol, temp) == true);
  REQUIRE(rdmol.getNumAtoms() == 6);
  REQUIRE(rdmol.getNumBonds() == 5);
  auto &mol = rdmol.asROMol();

  [[maybe_unused]] Atom *atom0 = mol.getAtomWithIdx(0);
  [[maybe_unused]] Atom *atom1 = mol.getAtomWithIdx(1);
  [[maybe_unused]] Atom *atom2 = mol.getAtomWithIdx(2);
  [[maybe_unused]] Atom *atom3 = mol.getAtomWithIdx(3);
  [[maybe_unused]] Atom *atom4 = mol.getAtomWithIdx(4);
  [[maybe_unused]] Atom *atom5 = mol.getAtomWithIdx(5);

  [[maybe_unused]] Bond *bond0 = mol.getBondWithIdx(0);
  [[maybe_unused]] Bond *bond1 = mol.getBondWithIdx(1);
  [[maybe_unused]] Bond *bond2 = mol.getBondWithIdx(2);
  [[maybe_unused]] Bond *bond3 = mol.getBondWithIdx(3);
  [[maybe_unused]] Bond *bond4 = mol.getBondWithIdx(4);

  SECTION("Empty") { CHECK(mol.getStereoGroups().size() == 0); }

  SECTION("Set/get") {
    std::vector<StereoGroup> groups;
    std::vector<Atom *> atoms = {atom0, atom2, atom1};
    std::vector<Bond *> bonds = {bond4, bond3};
    groups.push_back(StereoGroup(StereoGroupType::STEREO_OR, atoms, bonds));

    std::vector<Atom *> atoms2 = {atom5, atom4};
    std::vector<Bond *> bonds2 = {bond2};

    groups.push_back(StereoGroup(StereoGroupType::STEREO_AND, atoms2, bonds2));

    mol.setStereoGroups(std::move(groups));

    const std::vector<StereoGroup> &res = mol.getStereoGroups();
    REQUIRE(res.size() == 2);
    CHECK(res[0].getGroupType() == StereoGroupType::STEREO_OR);
    CHECK(res[1].getGroupType() == StereoGroupType::STEREO_AND);

    CHECK(res[0].getAtoms().size() == 3);
    CHECK(res[1].getAtoms().size() == 2);

    CHECK_THAT(res[0].getAtoms(), Catch::Matchers::UnorderedEquals(atoms));
    CHECK_THAT(res[1].getAtoms(), Catch::Matchers::UnorderedEquals(atoms2));

    CHECK_THAT(res[0].getBonds(), Catch::Matchers::UnorderedEquals(bonds));
    CHECK_THAT(res[1].getBonds(), Catch::Matchers::UnorderedEquals(bonds2));
  }

  SECTION("Set ROMol / get RDMol") {
    std::vector<StereoGroup> groups;
    std::vector<Atom *> atoms = {atom0, atom2, atom1};
    std::vector<Bond *> bonds = {bond4, bond3};
    groups.push_back(StereoGroup(StereoGroupType::STEREO_OR, std::move(atoms),
                                 std::move(bonds)));

    std::vector<Atom *> atoms2 = {atom5, atom4};
    std::vector<Bond *> bonds2 = {bond2};

    groups.push_back(StereoGroup(StereoGroupType::STEREO_AND, std::move(atoms2),
                                 std::move(bonds2)));

    mol.setStereoGroups(std::move(groups));

    const StereoGroups *res = rdmol.getStereoGroups();

    REQUIRE(res != nullptr);

    const auto &types = res->stereoTypes;
    const auto &atomBegins = res->atomBegins;
    const auto &bondBegins = res->bondBegins;
    const auto &atomIndices = res->atoms;
    const auto &bondIndices = res->bonds;

    CHECK(types == std::vector<StereoGroupType>({StereoGroupType::STEREO_OR,
                                                 StereoGroupType::STEREO_AND}));

    CHECK(atomBegins == std::vector<uint32_t>({0, 3, 5}));
    CHECK(bondBegins == std::vector<uint32_t>({0, 2, 3}));

    std::vector<uint32_t> atomSection1(3);
    std::vector<uint32_t> atomSection2(2);
    std::vector<uint32_t> bondSection1(2);
    std::vector<uint32_t> bondSection2(1);

    std::copy(atomIndices.begin(), atomIndices.begin() + 3,
              atomSection1.begin());
    std::copy(atomIndices.begin() + 3, atomIndices.begin() + 5,
              atomSection2.begin());
    std::copy(bondIndices.begin(), bondIndices.begin() + 2,
              bondSection1.begin());
    std::copy(bondIndices.begin() + 2, bondIndices.begin() + 3,
              bondSection2.begin());

    CHECK_THAT(atomSection1, Catch::Matchers::UnorderedEquals(
                                 std::vector<uint32_t>({0, 2, 1})));
    CHECK_THAT(atomSection2,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>({5, 4})));
    CHECK_THAT(bondSection1,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>({4, 3})));
    CHECK_THAT(bondSection2,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>({2})));
  }

  SECTION("Set RDMol / get ROMol") {
    auto groups = std::make_unique<StereoGroups>();
    groups->addGroup(StereoGroupType::STEREO_OR, {0, 2, 1}, {4, 3}, 1);
    groups->addGroup(StereoGroupType::STEREO_AND, {5, 4}, {2}, 2);
    rdmol.setStereoGroups(std::move(groups));

    const std::vector<StereoGroup> &res = mol.getStereoGroups();
    REQUIRE(res.size() == 2);
    CHECK(res[0].getGroupType() == StereoGroupType::STEREO_OR);
    CHECK(res[1].getGroupType() == StereoGroupType::STEREO_AND);

    CHECK(res[0].getAtoms().size() == 3);
    CHECK(res[1].getAtoms().size() == 2);

    CHECK(res[0].getReadId() == 1);
    CHECK(res[1].getReadId() == 2);

    CHECK_THAT(res[0].getAtoms(),
               Catch::Matchers::UnorderedEquals(
                   std::vector<Atom *>({atom0, atom2, atom1})));
    CHECK_THAT(res[1].getAtoms(), Catch::Matchers::UnorderedEquals(
                                      std::vector<Atom *>({atom5, atom4})));

    CHECK_THAT(res[0].getBonds(), Catch::Matchers::UnorderedEquals(
                                      std::vector<Bond *>({bond3, bond4})));
    CHECK_THAT(res[1].getBonds(),
               Catch::Matchers::UnorderedEquals(std::vector<Bond *>({bond2})));
  }

  SECTION("Threadsafety") {
    auto groups = std::make_unique<StereoGroups>();
    groups->addGroup(StereoGroupType::STEREO_OR, {0, 2, 1}, {4, 3}, 1);
    groups->addGroup(StereoGroupType::STEREO_AND, {5, 4}, {2}, 2);
    rdmol.setStereoGroups(std::move(groups));

    const size_t numThreads = 8;
    std::vector<std::thread> threads;
    std::atomic_int atomicInt(0);
    std::atomic<const std::vector<StereoGroup> *> atomicPtr(nullptr);
    std::atomic_int successCount(0);

    for (size_t i = 0; i < numThreads; ++i) {
      threads.emplace_back([&atomicInt, &atomicPtr, numThreads = numThreads,
                            &mol, &successCount]() {
        atomicInt.fetch_add(1);
        // Just spin as a barrier, to increase the chances of all threads
        // being in sync
        while (uint32_t(atomicInt.load()) != numThreads) {
        }
        // Have all threads call getStereoGroups at the same time
        const std::vector<StereoGroup> &res = mol.getStereoGroups();
        const std::vector<StereoGroup> *value = nullptr;
        if (!atomicPtr.compare_exchange_strong(value, &res)) {
          // All threads should have received the same address without
          // crashing
          CHECK(&res == value);
          if (&res == value) {
            ++successCount;
          }
        } else {
          ++successCount;
        }
      });
    }
    for (auto &thread : threads) {
      thread.join();
    }
    CHECK(successCount.load() == numThreads);
  }
}

TEST_CASE("Substance groups") {
  auto mol = basicMol();
  auto &romol = mol->asROMol();

  SECTION("No substance groups") { CHECK(getSubstanceGroups(romol).empty()); }

  SECTION("Set from RDMol") {
    std::vector<SubstanceGroup> &sgs = mol->getSubstanceGroups();
    sgs.emplace_back(&mol->asROMol(), "TEST");

    const std::vector<SubstanceGroup> &sgs2 = mol->getSubstanceGroups();
    CHECK(sgs2.size() == 1);
    CHECK(sgs2[0].getProp<std::string>("TYPE") == "TEST");
    CHECK(&sgs2[0].getOwningMol() == &mol->asROMol());
  }

  SECTION("Set from ROMol") {
    std::vector<SubstanceGroup> &sgs = getSubstanceGroups(romol);
    sgs.emplace_back(&romol, "TEST");

    const std::vector<SubstanceGroup> &sgs2 = mol->getSubstanceGroups();
    CHECK(sgs2.size() == 1);
    CHECK(sgs2[0].getProp<std::string>("TYPE") == "TEST");
    CHECK(&sgs2[0].getOwningMol() == &mol->asROMol());
  }

  SECTION("Clear") {
    std::vector<SubstanceGroup> &sgs = getSubstanceGroups(romol);
    sgs.emplace_back(&romol, "TEST");
    sgs.clear();
    CHECK(getSubstanceGroups(romol).empty());
  }
}

TEST_CASE("Conformers ROMol") {
  auto mol = basicMol();
  auto &romol = mol->asROMol();
  auto conf = std::make_unique<Conformer>(romol.getNumAtoms());
  for (size_t i = 0; i < romol.getNumAtoms(); ++i) {
    conf->setAtomPos(i, RDGeom::Point3D(i, i, i));
  }

  SECTION("No conformers") { CHECK(romol.getNumConformers() == 0); }

  SECTION("Add conformer to empty mol") {
    ROMol emptyMol;
    auto emptyConf = std::make_unique<Conformer>();
    emptyMol.addConformer(emptyConf.release());
    REQUIRE(emptyMol.getNumConformers() == 1);
    CHECK(emptyMol.getConformer().getNumAtoms() == 0);
    CHECK(&emptyMol.getConformer().getOwningMol() == &emptyMol);
  }

  SECTION("Add conformer fails - nullptr input") {
    CHECK_THROWS_AS(romol.addConformer(nullptr), Invar::Invariant);
  }

  SECTION("Add conformer fails - wrong number of atoms") {
    auto badConf = std::make_unique<Conformer>(romol.getNumAtoms() - 1);
    CHECK_THROWS_AS(romol.addConformer(badConf.get()), Invar::Invariant);
  }

  SECTION("Add conformer") {
    romol.addConformer(conf.release());
    CHECK(romol.getNumConformers() == 1);
    const auto &resultConf = romol.getConformer();
    CHECK(resultConf.getNumAtoms() == 4);
    CHECK(&resultConf.getOwningMol() == &romol);
    CHECK(resultConf.getId() == 0);

    // Check positions
    for (size_t i = 0; i < romol.getNumAtoms(); ++i) {
      CHECK(resultConf.getAtomPos(i).x == i);
    }
  }

  SECTION("Add conformer with nonlinear ID") {
    conf->setId(3);
    romol.addConformer(conf.release());
    CHECK(romol.getNumConformers() == 1);
    const auto &resultConf = romol.getConformer();
    CHECK(resultConf.getId() == 3);
  }

  SECTION("Add conformer with assignId") {
    conf->setId(3);
    romol.addConformer(conf.release(), /*assignID=*/true);
    CHECK(romol.getNumConformers() == 1);
    const auto &resultConf = romol.getConformer();
    CHECK(resultConf.getId() == 0);
  }

  SECTION("Add 2D conformer") {
    conf->set3D(false);
    romol.addConformer(conf.release());
    CHECK(romol.getNumConformers() == 1);
    const auto &resultConf = romol.getConformer();
    CHECK(resultConf.is3D() == false);
  }

  SECTION("Remove conformer") {
    auto secondConf = std::make_unique<Conformer>(romol.getNumAtoms());
    for (size_t i = 0; i < romol.getNumAtoms(); ++i) {
      secondConf->setAtomPos(
          i, RDGeom::Point3D(i + romol.getNumAtoms(), i + romol.getNumAtoms(),
                             i + romol.getNumAtoms()));
    }
    secondConf->set3D(false);
    romol.addConformer(conf.release());
    romol.addConformer(secondConf.release(), /*assignId=*/true);
    CHECK(romol.getNumConformers() == 2);
    romol.removeConformer(0);
    CHECK(romol.getNumConformers() == 1);

    auto &resultConf = romol.getConformer();
    CHECK(resultConf.getNumAtoms() == 4);
    CHECK(resultConf.getId() == 1);
    CHECK(resultConf.is3D() == false);
    // Check positions
    for (size_t i = 0; i < romol.getNumAtoms(); ++i) {
      CHECK(resultConf.getAtomPos(i).x == i + romol.getNumAtoms());
    }
  }

  SECTION("Remove conformer no-op null ID") {
    romol.addConformer(conf.release());
    CHECK(romol.getNumConformers() == 1);
    romol.removeConformer(1);
    CHECK(romol.getNumConformers() == 1);
  }

  SECTION("Clear conformers") {
    auto secondConf = std::make_unique<Conformer>(romol.getNumAtoms());
    romol.addConformer(conf.release());
    romol.addConformer(secondConf.release());
    CHECK(romol.getNumConformers() == 2);
    romol.clearConformers();
    CHECK(romol.getNumConformers() == 0);
  }
}

TEST_CASE("Conformers RDMol/ROMol interop") {
  SECTION("RDMol setup, ROMol access") {
    RDMol mol;
    mol.addAtom();
    mol.addAtom();

    auto [confId, confPos] = mol.addConformer(3, true);
    std::iota(confPos, confPos + 6, 0);

    auto [confId2, confPos2] = mol.addConformer(-1, false);
    std::iota(confPos2, confPos2 + 6, 6);

    auto &romol = mol.asROMol();
    CHECK(romol.getNumConformers() == 2);
    const auto &conf = romol.getConformer(confId);
    CHECK(conf.getNumAtoms() == 2);
    CHECK(conf.getAtomPos(0).x == 0);
    CHECK(conf.getAtomPos(1).x == 3);
    CHECK(conf.getId() == confId);
    CHECK(conf.is3D() == true);

    const auto &conf2 = romol.getConformer(confId2);
    CHECK(conf2.getNumAtoms() == 2);
    CHECK(conf2.getAtomPos(0).x == 6);
    CHECK(conf2.getAtomPos(1).x == 9);
    CHECK(conf2.getId() == confId2);
    CHECK(conf2.is3D() == false);

    mol.removeConformer(confId);
    CHECK(romol.getNumConformers() == 1);
    CHECK(romol.getConformer().getId() == confId2);
  }

  SECTION("RDMol setup, add to existing ROMol") {
    // Identical setup to above, except the RDMol APIs are called after
    // the compat data is set up, exercising the in-call updates rather than
    // the initial setup of compat data.
    RDMol mol;
    auto &romol = mol.asROMol();

    mol.addAtom();
    mol.addAtom();

    std::vector<double> pos1(6, 0);
    std::iota(pos1.begin(), pos1.end(), 0);
    std::vector<double> pos2(6, 0);
    std::iota(pos2.begin(), pos2.end(), 6);

    auto confId = mol.addConformer(pos1.data(), 3, true);
    auto confId2 = mol.addConformer(pos2.data(), -1, false);

    CHECK(romol.getNumConformers() == 2);
    const auto &conf = romol.getConformer(confId);
    CHECK(conf.getNumAtoms() == 2);
    CHECK(conf.getAtomPos(0).x == 0);
    CHECK(conf.getAtomPos(1).x == 3);
    CHECK(conf.getId() == confId);
    CHECK(conf.is3D() == true);

    const auto &conf2 = romol.getConformer(confId2);
    CHECK(conf2.getNumAtoms() == 2);
    CHECK(conf2.getAtomPos(0).x == 6);
    CHECK(conf2.getAtomPos(1).x == 9);
    CHECK(conf2.getId() == confId2);
    CHECK(conf2.is3D() == false);

    mol.removeConformer(confId);
    CHECK(romol.getNumConformers() == 1);
    CHECK(romol.getConformer().getId() == confId2);
  }

  SECTION("ROMol setup, RDMol access") {
    RDMol mol;
    mol.addAtom();
    mol.addAtom();
    ROMol &romol = mol.asROMol();

    auto conf1 = std::make_unique<Conformer>(2);
    conf1->setAtomPos(0, RDGeom::Point3D(0, 1, 2));
    conf1->setAtomPos(1, RDGeom::Point3D(3, 4, 5));
    auto conf2 = std::make_unique<Conformer>(2);
    conf2->setAtomPos(0, RDGeom::Point3D(6, 7, 8));
    conf2->setAtomPos(1, RDGeom::Point3D(9, 10, 11));

    conf1->setId(3);
    conf2->set3D(false);

    romol.addConformer(conf1.release());
    romol.addConformer(conf2.release(), true);

    CHECK(mol.getNumConformers() == 2);
    const double *conf1Pos = mol.getConformerPositions(3);
    const double *conf2Pos = mol.getConformerPositions(4);
    REQUIRE(conf1Pos != nullptr);
    REQUIRE(conf2Pos != nullptr);

    CHECK(conf1Pos[0] == 0);
    CHECK(conf1Pos[3] == 3);
    CHECK(conf2Pos[0] == 6);
    CHECK(conf2Pos[3] == 9);

    CHECK(mol.is3DConformer(3) == true);
    CHECK(mol.is3DConformer(4) == false);

    romol.removeConformer(3);
    CHECK(mol.getNumConformers() == 1);
    double *pos = mol.getConformerPositions(4);
    REQUIRE(pos != nullptr);
    CHECK(pos[0] == 6);
    CHECK(pos[3] == 9);
    CHECK(mol.is3DConformer(4) == false);
  }

  SECTION("Syncing updates") {
    RDMol mol;
    mol.addAtom();
    mol.addAtom();
    ROMol &romol = mol.asROMol();
    auto conf = std::make_unique<Conformer>(2);
    conf->setAtomPos(0, RDGeom::Point3D(0, 1, 2));
    romol.addConformer(conf.release());

    double *pos = mol.getConformerPositions(0);
    pos[0]++;
    pos[1]++;
    pos[2]++;

    Conformer &confHandle = romol.getConformer();
    CHECK(confHandle.getAtomPos(0).x == 1);
    CHECK(confHandle.getAtomPos(0).y == 2);
    CHECK(confHandle.getAtomPos(0).z == 3);

    confHandle.getAtomPos(0).x++;
    confHandle.getAtomPos(0).y++;
    confHandle.getAtomPos(0).z++;

    const double *nextCheckedPos = mol.getConformerPositions(0);
    CHECK(nextCheckedPos[0] == 2);
    CHECK(nextCheckedPos[1] == 3);
    CHECK(nextCheckedPos[2] == 4);
  }

  SECTION("Multithreaded access") {
    RDMol mol;
    mol.addAtom();
    mol.addAtom();
    ROMol &romol = mol.asROMol();
    auto conf = std::make_unique<Conformer>(2);
    conf->setAtomPos(0, RDGeom::Point3D(0, 1, 2));
    romol.addConformer(conf.release());
    auto &romolPos = romol.getConformer();
    // Update the compat position to make sure they're not equal to the
    // RDMol positions.
    auto &atomPos = romolPos.getAtomPos(0);
    atomPos.x++;
    atomPos.y++;
    atomPos.z++;
    REQUIRE(atomPos.x == 1);
    REQUIRE(atomPos.y == 2);
    REQUIRE(atomPos.z == 3);

    auto checkPos = [&mol]() {
      const double *rdmolPos = mol.getConformerPositions(0);
      CHECK(rdmolPos[0] == 1);
      CHECK(rdmolPos[1] == 2);
      CHECK(rdmolPos[2] == 3);
    };
    checkPos();
    std::thread checkThread1(checkPos);
    std::thread checkThread2(checkPos);
    checkThread1.join();
    checkThread2.join();
  }
}

TEST_CASE("RingInfo RDMol/ROMol interop") {
  SECTION("Test of sync bug") {
    // Two pentagonal rings with a shared edge
    auto mol = parseSmiles("C1CCC2C1CCC2");
    MolOps::findSSSR(mol->asROMol());

    REQUIRE(mol->asROMol().getRingInfo()->numRings() == 2);

    // There was a bug in non-const RDMol::getRingInfo where it incorrectly
    // marked the new and compat data structures as in sync, instead of
    // the new data structure being modified.
    RingInfoCache &ringInfoNew = mol->getRingInfo();

    CHECK(ringInfoNew.numRings() == 2);
    CHECK((ringInfoNew.ringBegins[0] == 0 && ringInfoNew.ringBegins[1] == 5 &&
           ringInfoNew.ringBegins[2] == 10));
    CHECK(ringInfoNew.findRingType == FIND_RING_TYPE_SSSR);

    // Modify the new data structure
    ringInfoNew.reset();
    ringInfoNew.findRingType = FIND_RING_TYPE_FAST;
    ringInfoNew.isInit = true;

    // Ensure that the compat data structure receives the modification.
    CHECK(mol->asROMol().getRingInfo()->numRings() == 0);
    CHECK(mol->asROMol().getRingInfo()->getRingType() == FIND_RING_TYPE_FAST);
  }

  SECTION("Threadsafety") {
    auto mol = parseSmiles("C1CCC2C1CCC2");
    MolOps::findSSSR(mol->asROMol());

    REQUIRE(mol->asROMol().getRingInfo()->numRings() == 2);

    const size_t numThreads = 8;
    std::vector<std::thread> threads;
    std::atomic_int atomicInt(0);
    std::atomic<const uint32_t *> atomicPtr(nullptr);
    std::atomic_int successCount(0);

    for (size_t i = 0; i < numThreads; ++i) {
      // Some compilers require numThreads to be captured, but Clang will
      // error that the numThreads lambda capture isn't necessary.
      // Specifying the capture as "numThreads = numThreads" forces copying,
      // which works around the error in Clang and still captures it for other
      // compilers.
      threads.emplace_back([&atomicInt, &atomicPtr, numThreads = numThreads,
                            &mol, &successCount]() {
        atomicInt.fetch_add(1);
        // Just spin as a barrier, to increase the chances of all threads
        // being in sync
        while (uint32_t(atomicInt.load()) != numThreads) {
        }
        // Have all threads call getRingInfo at the same time
        const RingInfoCache &res = mol->getRingInfo();
        // The ring info is a member variable, so they'd get the same
        // address of res regardless of whether the results clobbered each
        // other, so check one of the vectors.
        const uint32_t *atomsInRings = res.atomsInRings.data();
        const uint32_t *value = nullptr;
        if (!atomicPtr.compare_exchange_strong(value, atomsInRings)) {
          // All threads should have received the same address without
          // crashing
          CHECK(atomsInRings == value);
          if (atomsInRings == value) {
            ++successCount;
          }
        } else {
          ++successCount;
        }
      });
    }
    for (auto &thread : threads) {
      thread.join();
    }
    CHECK(successCount.load() == numThreads);
  }
}