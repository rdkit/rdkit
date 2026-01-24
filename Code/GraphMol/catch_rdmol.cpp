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
#include <catch2/catch_all.hpp>
#include <catch2/catch_session.hpp>

#include <GraphMol/RDMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

namespace {

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

static PropToken unusedToken("unused");
static PropToken unsignedIntToken("unsigned int");
static PropToken floatToken("float");
static PropToken complexToken("nonPOD");
static PropToken signedIntToken("signed int");

}  // namespace

void checkPropsAndGroups(const RDMol &mol) {
  // Check mol props
  // These values are never used, but GCC can no longer deduce that they're
  // initialized, and maybe-uninitialized is currently treated as an error.
  float resf = 0.0f;
  unsigned int resu = 0;
  CHECK(mol.getMolPropIfPresent(floatToken, resf));
  CHECK(resf == 2.0f);

  const int *resi = mol.getAtomPropArrayIfPresent<int>(signedIntToken);
  REQUIRE(resi != nullptr);
  CHECK(resi[1] == 3);

  CHECK(mol.getBondPropIfPresent(unsignedIntToken, 1, resu));
  CHECK(resu == 4u);

  const auto *stereoGroups = mol.getStereoGroups();
  REQUIRE(stereoGroups != nullptr);
  CHECK(stereoGroups->getGroupType(0) == StereoGroupType::STEREO_ABSOLUTE);
  const auto &[begAtom, endAtom] = stereoGroups->getAtoms(0);
  CHECK(*begAtom == 1);
  CHECK(*(begAtom + 1) == 2);

  const auto &substanceGroups = mol.getSubstanceGroups();
  REQUIRE(substanceGroups.size() == 1);
  CHECK(substanceGroups[0].getProp<std::string>("TYPE") == "TEST");
  CHECK(&substanceGroups[0].getOwningMol() == &mol.asROMol());
}

void checkConformers(const RDMol &mol) {
  REQUIRE(mol.getNumConformers() == 2);

  const auto *conf1 = mol.getConformerPositions(0);
  const auto *conf2 = mol.getConformerPositions(5);

  REQUIRE(conf1 != nullptr);
  REQUIRE(conf2 != nullptr);

  for (int i = 0; i < 12; ++i) {
    CHECK(conf1[i] == 1.0);
    CHECK(conf2[i] == 2.0);
  }

  CHECK(mol.is3DConformer(0));
  CHECK(!mol.is3DConformer(5));
}

TEST_CASE("RDMol Constructors") {
  auto molPtr = basicMol();
  RDMol &mol = *molPtr;
  mol.setAtomBookmark(1, 0);

  AtomData &atom1 = mol.getAtom(1);
  atom1.setAtomicNum(4);

  mol.setMolProp(floatToken, 2.0f);
  mol.addAtomProp(signedIntToken, 3);
  mol.setSingleBondProp(unsignedIntToken, 1, 4u);

  // Substance groups
  std::vector<SubstanceGroup> &sgs = mol.getSubstanceGroups();
  sgs.emplace_back(&mol.asROMol(), "TEST");

  // Stereo groups
  auto stereoGroups = std::make_unique<StereoGroups>();
  stereoGroups->addGroup(StereoGroupType::STEREO_ABSOLUTE, {1, 2}, {1});
  mol.setStereoGroups(std::move(stereoGroups));

  // Conformers
  auto [conf1Id, conf1Pos] = mol.addConformer();
  std::fill(conf1Pos, conf1Pos + 12, 1.0);
  auto [conf2Id, conf2Pos] = mol.addConformer(5, false);
  std::fill(conf2Pos, conf2Pos + 12, 2.0);

  SECTION("Default") {
    RDMol mol2;
    CHECK(mol2.getNumAtoms() == 0);
    CHECK(mol2.getNumBonds() == 0);
  }

  SECTION("Copy constructor, no params") {
    RDMol mol2(mol);
    CHECK(mol2.getNumAtoms() == 4);
    CHECK(mol2.getNumBonds() == 3);

    CHECK(mol2.hasAtomBookmark(0));
    mol.clearAllAtomBookmarks();
    CHECK(mol2.hasAtomBookmark(0));

    AtomData &mol2Atom1 = mol2.getAtom(1);
    // First check that the set before copy propagated.
    CHECK(mol2Atom1.getAtomicNum() == 4);
    atom1.setAtomicNum(8);
    // Now check that the copies aren't linked.
    CHECK(mol2Atom1.getAtomicNum() == 4);

    checkPropsAndGroups(mol2);
    checkConformers(mol2);
  }

  SECTION("Copy constructor, quick") {
    RDMol mol2(mol, true);
    CHECK(mol2.getNumAtoms() == 4);
    CHECK(mol2.getNumBonds() == 3);

    // Bookmarks and properties should not have been copied.
    CHECK(!mol2.hasAtomBookmark(0));
    CHECK(!mol2.hasProp(floatToken));
    CHECK(mol2.getAtomPropArrayIfPresent<int>(signedIntToken) == nullptr);
  }

  SECTION("Copy constructor, not quick with ID") {
    RDMol mol2(mol, false, 5);
    CHECK(mol2.getNumAtoms() == 4);
    CHECK(mol2.getNumBonds() == 3);

    CHECK(mol2.getNumConformers() == 1);
    const double *confPos = mol2.getConformerPositions(5);
    REQUIRE(confPos != nullptr);
    CHECK(confPos[0] == 2.0);
    CHECK(confPos[11] == 2.0);
  }

  SECTION("copy constructor, not quick with non-existent ID") {
    RDMol mol2(mol, false, 6);
    CHECK(mol2.getNumAtoms() == 4);
    CHECK(mol2.getNumBonds() == 3);

    CHECK(mol2.getNumConformers() == 0);
  }

  SECTION("Copy assignment operator") {
    RDMol mol2;
    mol2 = mol;
    CHECK(mol2.getNumAtoms() == 4);
    CHECK(mol2.getNumBonds() == 3);

    CHECK(mol2.hasAtomBookmark(0));
    mol.clearAllAtomBookmarks();
    CHECK(mol2.hasAtomBookmark(0));

    AtomData &mol2Atom1 = mol2.getAtom(1);
    // First check that the set before copy propagated.
    CHECK(mol2Atom1.getAtomicNum() == 4);
    atom1.setAtomicNum(8);
    // Now check that the copies aren't linked.
    CHECK(mol2Atom1.getAtomicNum() == 4);

    checkPropsAndGroups(mol2);
    checkConformers(mol2);
  }

  SECTION("Move constructor") {
    RDMol mol2(std::move(mol));

    CHECK(mol2.getNumAtoms() == 4);
    CHECK(mol2.getNumBonds() == 3);
    CHECK(mol2.hasAtomBookmark(0));

    AtomData &mol2Atom = mol2.getAtom(1);
    CHECK(mol2Atom.getAtomicNum() == 4);

    checkPropsAndGroups(mol2);
    checkConformers(mol2);
  }

  SECTION("Move assignment operator") {
    RDMol mol2;
    mol2 = std::move(mol);

    CHECK(mol2.getNumAtoms() == 4);
    CHECK(mol2.getNumBonds() == 3);
    CHECK(mol2.hasAtomBookmark(0));

    AtomData &mol2Atom = mol2.getAtom(1);
    CHECK(mol2Atom.getAtomicNum() == 4);

    checkPropsAndGroups(mol2);
    checkConformers(mol2);
  }

  SECTION("Copy constructor with compat data copies back") {
    auto &romol = mol.asROMol();
    romol.setAtomBookmark(romol.getAtomWithIdx(1), 5);
    romol.setBondBookmark(romol.getBondWithIdx(2), 6);
    romol.getBondWithIdx(2)->setStereoAtoms(1, 2);

    RDMol mol2(mol);

    CHECK(mol2.getAtomWithBookmark(0) == 1);
    CHECK(mol2.getAtomWithBookmark(5) == 1);
    CHECK(mol2.getBondWithBookmark(6) == 2);

    const auto *stereoAtoms = mol2.getBondStereoAtoms(2);
    REQUIRE(stereoAtoms != nullptr);
    CHECK(stereoAtoms[0] == 1);
    CHECK(stereoAtoms[1] == 2);
  }
}

TEST_CASE("RDMol Bookmarks") {
  RDMol mol;
  RDKit::SmilesParseTemp temp;
  static const char *basicSmiles = "CCCC";

  SmilesParserParams params;
  params.sanitize = false;
  params.removeHs = false;

  REQUIRE(RDKit::SmilesToMol(basicSmiles, params, mol, temp) == true);
  REQUIRE(mol.getNumAtoms() == 4);
  REQUIRE(mol.getNumBonds() == 3);

  SECTION("Basic bookmarking") {
    mol.setAtomBookmark(1, 0);
    mol.setAtomBookmark(3, -6);
    mol.setBondBookmark(1, 0);

    CHECK(mol.getAtomWithBookmark(0) == 1);
    CHECK(mol.getAtomWithBookmark(-6) == 3);
    CHECK(mol.getBondWithBookmark(0) == 1);
  }

  SECTION("Multiple atoms in bookmark") {
    mol.setAtomBookmark(1, 0);
    mol.setAtomBookmark(2, 0);
    mol.setAtomBookmark(3, 0);

    CHECK(mol.getAtomWithBookmark(0) == 1);

    std::vector<int> res = mol.getAllAtomsWithBookmarks(0);
    std::vector<int> want = {1, 2, 3};
    CHECK(res == want);
  }

  SECTION("Multiple bonds in bookmark") {
    mol.setBondBookmark(3, 0);
    mol.setBondBookmark(2, 0);
    mol.setBondBookmark(1, 0);

    CHECK(mol.getBondWithBookmark(0) == 3);

    std::vector<int> res = mol.getAllBondsWithBookmarks(0);
    std::vector<int> want = {3, 2, 1};
    CHECK(res == want);
  }

  SECTION("Replace atom bookmarks") {
    mol.setAtomBookmark(3, 1);
    mol.setAtomBookmark(0, 0);
    mol.setAtomBookmark(2, 0);

    mol.replaceAtomBookmark(1, 0);

    CHECK(mol.getAtomWithBookmark(0) == 1);

    std::vector<int> res = mol.getAllAtomsWithBookmarks(0);
    std::vector<int> want = {1};
    CHECK(res == want);
  }

  SECTION("Multiple bookmarks same atom is fine") {
    mol.setAtomBookmark(1, 0);
    mol.setAtomBookmark(1, 1);
    mol.setAtomBookmark(1, 2);

    CHECK(mol.getAtomWithBookmark(0) == 1);
    CHECK(mol.getAtomWithBookmark(1) == 1);
    CHECK(mol.getAtomWithBookmark(2) == 1);

    std::vector<int> res = mol.getAllAtomsWithBookmarks(0);
    std::vector<int> want = {1};
    CHECK(res == want);
  }

  SECTION("Missing bookmark throws") {
    mol.setAtomBookmark(1, 1);
    mol.setBondBookmark(1, 2);

    CHECK_THROWS_AS(mol.getAtomWithBookmark(0), Invar::Invariant);
    CHECK_THROWS_AS(mol.getBondWithBookmark(0), Invar::Invariant);

    CHECK_THROWS_AS(mol.getAllAtomsWithBookmarks(0), Invar::Invariant);
    CHECK_THROWS_AS(mol.getAllBondsWithBookmarks(0), Invar::Invariant);
  }

  SECTION("Clear atom bookmarks") {
    mol.setAtomBookmark(1, 0);
    mol.setAtomBookmark(1, 1);
    mol.setAtomBookmark(1, 2);

    mol.clearAtomBookmark(0);
    CHECK_THROWS_AS(mol.getAtomWithBookmark(0), Invar::Invariant);
    CHECK(mol.getAtomWithBookmark(1) == 1);
    CHECK(mol.getAtomWithBookmark(2) == 1);

    mol.clearAllAtomBookmarks();
    CHECK_THROWS_AS(mol.getAtomWithBookmark(0), Invar::Invariant);
    CHECK_THROWS_AS(mol.getAtomWithBookmark(1), Invar::Invariant);
    CHECK_THROWS_AS(mol.getAtomWithBookmark(2), Invar::Invariant);
  }

  SECTION("Clear bond bookmarks") {
    mol.setBondBookmark(1, 0);
    mol.setBondBookmark(1, 1);
    mol.setBondBookmark(1, 2);

    mol.clearBondBookmark(0);
    CHECK_THROWS_AS(mol.getBondWithBookmark(0), Invar::Invariant);
    CHECK(mol.getBondWithBookmark(1) == 1);
    CHECK(mol.getBondWithBookmark(2) == 1);

    mol.clearAllBondBookmarks();
    CHECK_THROWS_AS(mol.getBondWithBookmark(0), Invar::Invariant);
    CHECK_THROWS_AS(mol.getBondWithBookmark(1), Invar::Invariant);
    CHECK_THROWS_AS(mol.getBondWithBookmark(2), Invar::Invariant);
  }

  SECTION("Clear individual atom idx bookmark") {
    mol.setAtomBookmark(0, 1);
    mol.setAtomBookmark(1, 1);
    mol.setAtomBookmark(2, 1);

    mol.clearAtomBookmark(1, 1);
    std::vector<int> res = mol.getAllAtomsWithBookmarks(1);
    std::vector<int> want = {0, 2};
    CHECK(res == want);
  }

  SECTION("Clear individual bond idx bookmark") {
    mol.setBondBookmark(0, 1);
    mol.setBondBookmark(1, 1);
    mol.setBondBookmark(2, 1);

    mol.clearBondBookmark(1, 1);
    std::vector<int> res = mol.getAllBondsWithBookmarks(1);
    std::vector<int> want = {0, 2};
    CHECK(res == want);
  }

  SECTION("Has bookmark") {
    mol.setAtomBookmark(1, 0);
    mol.setAtomBookmark(1, 1);
    mol.setAtomBookmark(1, 2);

    CHECK(mol.hasAtomBookmark(0));
    CHECK(mol.hasAtomBookmark(1));
    CHECK(mol.hasAtomBookmark(2));
    CHECK(!mol.hasAtomBookmark(3));

    mol.setBondBookmark(1, 0);
    mol.setBondBookmark(1, 1);
    mol.setBondBookmark(1, 2);

    CHECK(mol.hasBondBookmark(0));
    CHECK(mol.hasBondBookmark(1));
    CHECK(mol.hasBondBookmark(2));
    CHECK(!mol.hasBondBookmark(3));

    mol.clearAtomBookmark(0);
    CHECK(!mol.hasAtomBookmark(0));
  }
}

TEST_CASE("RDMol Mol properties") {
  RDMol mol;

  SECTION("Self basic properties") {
    mol.addMolProp(signedIntToken, 3);
    mol.addMolProp(unsignedIntToken, 3u);
    mol.addMolProp(floatToken, 3.0f);

    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    REQUIRE(mol.getMolPropIfPresent<int>(signedIntToken, res));
    CHECK(res == 3);

    REQUIRE(!mol.getMolPropIfPresent<int>(unusedToken, res));
  }

  SECTION("Prop handle updates in place") {
    mol.addMolProp(signedIntToken, 3);
    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    REQUIRE(mol.getMolPropIfPresent<int>(signedIntToken, res));
    REQUIRE(res == 3);

    mol.setMolProp(signedIntToken, 4);
    REQUIRE(mol.getMolPropIfPresent<int>(signedIntToken, res));
    CHECK(res == 4);
  }

  SECTION("Clear Mol Props") {
    mol.addMolProp(floatToken, 3.0f);
    mol.addMolProp(unsignedIntToken, 3u);

    mol.clearMolPropIfPresent(unsignedIntToken);
    mol.clearMolPropIfPresent(floatToken);

    float resf;
    unsigned int resu;
    CHECK(!mol.getMolPropIfPresent<unsigned int>(unsignedIntToken, resu));
    CHECK(!mol.getMolPropIfPresent<float>(floatToken, resf));
  }

  SECTION("Computed props and clear") {
    mol.addMolProp(floatToken, 3.0f, true);
    mol.addMolProp(unsignedIntToken, 3u);

    mol.clearComputedProps();
    float resf;
    unsigned int resu;
    CHECK(mol.getMolPropIfPresent<unsigned int>(unsignedIntToken, resu));
    CHECK(!mol.getMolPropIfPresent<float>(floatToken, resf));
  }

  SECTION("Mol prop nonstandard data - nonPOD") {
    mol.addMolProp(signedIntToken,
                   boost::shared_array<int>(new int[3]{1, 2, 3}));

    boost::shared_array<int> res;
    REQUIRE(mol.getMolPropIfPresent(signedIntToken, res));

    CHECK(res[1] == 2);
  }

  SECTION("String to int conversion") {
    mol.addMolProp(signedIntToken, std::string("3"));
    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    REQUIRE(mol.getMolPropIfPresent<int>(signedIntToken, res));
    CHECK(res == 3);
  }

  SECTION("char to int conversion") {
    mol.addMolProp(signedIntToken, "30");
    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    REQUIRE(mol.getMolPropIfPresent<int>(signedIntToken, res));
    CHECK(res == 30);
  }
}

TEST_CASE("RDMol Atom/bond properties") {
  RDMol mol;
  RDKit::SmilesParseTemp temp;
  static const char *basicSmiles = "CCCC";

  SmilesParserParams params;
  params.sanitize = false;
  params.removeHs = false;

  REQUIRE(RDKit::SmilesToMol(basicSmiles, params, mol, temp) == true);
  REQUIRE(mol.getNumAtoms() == 4);
  REQUIRE(mol.getNumBonds() == 3);

  SECTION("Atom properties access and override") {
    int *atomProps = mol.addAtomProp(signedIntToken, 3);
    REQUIRE(atomProps != nullptr);
    CHECK(*atomProps == 3);
    CHECK(*(atomProps + 3) == 3);

    int *atomProps2 = mol.getAtomPropArrayIfPresent<int>(signedIntToken);
    CHECK(atomProps == atomProps2);
    REQUIRE(mol.addAtomProp(signedIntToken, 4));
    int *atomProps3 = mol.getAtomPropArrayIfPresent<int>(signedIntToken);
    CHECK(*atomProps3 == 3);
    CHECK(*(atomProps + 3) == 3);
  }

  SECTION("Atom single access/set from existing property") {
    REQUIRE(mol.addAtomProp(signedIntToken, 3) != nullptr);
    mol.setSingleAtomProp(signedIntToken, 1, 4);
    int res1, res2;
    REQUIRE(mol.getAtomPropIfPresent<int>(signedIntToken, 0, res1));
    REQUIRE(mol.getAtomPropIfPresent<int>(signedIntToken, 1, res2));
    CHECK(res1 == 3);
    CHECK(res2 == 4);
    CHECK(!mol.getAtomPropIfPresent<int>(unusedToken, 2, res1));
  }

  SECTION("Atom single access/set new property") {
    mol.setSingleAtomProp<int>(signedIntToken, 1, 4);
    // Not present, the above should only set the second atom present bit.
    int res1, res2;
    REQUIRE(!mol.getAtomPropIfPresent<int>(signedIntToken, 0, res1));
    REQUIRE(mol.getAtomPropIfPresent<int>(signedIntToken, 1, res2));
    CHECK(res2 == 4);
  }

  SECTION("Atom array props basic") {
    const std::vector<int> defaultVals = {1, 2, 3};
    mol.addAtomProp<std::vector<int>>(signedIntToken, defaultVals);

    std::vector<int> res;
    REQUIRE(mol.getAtomPropIfPresent<std::vector<int>>(signedIntToken, 1, res));
    REQUIRE(mol.getAtomPropArrayIfPresent<int>(unsignedIntToken) == nullptr);

    CHECK(res == defaultVals);
  }

  SECTION("Atom prop remove all atoms") {
    mol.setSingleAtomProp(signedIntToken, 1, 7);
    mol.setSingleAtomProp(signedIntToken, 2, 11);
    int *prop = mol.getAtomPropArrayIfPresent<int>(signedIntToken);
    REQUIRE(prop != nullptr);
    CHECK(prop[1] == 7);
    CHECK(prop[2] == 11);

    mol.clearSingleAtomProp(signedIntToken, 1);
    mol.clearSingleAtomProp(signedIntToken, 2);

    prop = mol.getAtomPropArrayIfPresent<int>(signedIntToken);
    CHECK(prop == nullptr);
  }

  SECTION("Bond properties access and override") {
    int *BondProps = mol.addBondProp(signedIntToken, 3);
    REQUIRE(BondProps != nullptr);
    CHECK(*BondProps == 3);
    CHECK(*(BondProps + 1) == 3);

    int *BondProps2 = mol.getBondPropArrayIfPresent<int>(signedIntToken);
    CHECK(BondProps == BondProps2);
    REQUIRE(mol.addBondProp(signedIntToken, 4));
    int *BondProps3 = mol.getBondPropArrayIfPresent<int>(signedIntToken);
    CHECK(*BondProps3 == 3);
    CHECK(*(BondProps + 1) == 3);
  }

  SECTION("Bond single access/set from existing property") {
    REQUIRE(mol.addBondProp(signedIntToken, 3) != nullptr);
    mol.setSingleBondProp(signedIntToken, 1, 4);
    int res1, res2;
    REQUIRE(mol.getBondPropIfPresent<int>(signedIntToken, 0, res1));
    REQUIRE(mol.getBondPropIfPresent<int>(signedIntToken, 1, res2));
    CHECK(res1 == 3);
    CHECK(res2 == 4);
    CHECK(!mol.getBondPropIfPresent<int>(unusedToken, 2, res1));
  }

  SECTION("Bond single access/set new property") {
    mol.setSingleBondProp<int>(signedIntToken, 1, 4);
    // Not present, the above should only set the second Bond present bit.
    int res1, res2;
    REQUIRE(!mol.getBondPropIfPresent<int>(signedIntToken, 0, res1));
    REQUIRE(mol.getBondPropIfPresent<int>(signedIntToken, 1, res2));
    CHECK(res2 == 4);
  }

  SECTION("Bond array props basic") {
    const std::vector<int> defaultVals = {1, 2, 3};
    mol.addBondProp<std::vector<int>>(signedIntToken, defaultVals);

    std::vector<int> res;
    REQUIRE(mol.getBondPropIfPresent<std::vector<int>>(signedIntToken, 1, res));

    REQUIRE(mol.getBondPropArrayIfPresent<int>(unsignedIntToken) == nullptr);
    CHECK(res == defaultVals);
  }

  SECTION("Bond prop remove all bonds") {
    mol.setSingleBondProp(signedIntToken, 1, 7);
    mol.setSingleBondProp(signedIntToken, 2, 11);
    int *prop = mol.getBondPropArrayIfPresent<int>(signedIntToken);
    REQUIRE(prop != nullptr);
    CHECK(prop[1] == 7);
    CHECK(prop[2] == 11);

    mol.clearSingleBondProp(signedIntToken, 1);
    mol.clearSingleBondProp(signedIntToken, 2);

    prop = mol.getBondPropArrayIfPresent<int>(signedIntToken);
    CHECK(prop == nullptr);
  }

  SECTION("Get array prop of different type returns null") {
    mol.setSingleAtomProp<int>(signedIntToken, 1, 4);
    mol.setSingleBondProp<int>(signedIntToken, 1, 4);
    REQUIRE(mol.getBondPropArrayIfPresent<float>(signedIntToken) == nullptr);
    REQUIRE(mol.getAtomPropArrayIfPresent<float>(signedIntToken) == nullptr);
  }

  SECTION("Atom/bond/mol do not clash") {
    mol.addMolProp(floatToken, 3.0f);
    mol.addAtomProp(floatToken, 4.0f);
    mol.addBondProp(floatToken, 5.0f);

    // These values are never used, but GCC can no longer deduce that they're
    // initialized, and maybe-uninitialized is currently treated as an error.
    float molRes = 0.0f;
    float AtomRes = 0.0f;
    float bondRes = 0.0f;

    REQUIRE(mol.getMolPropIfPresent<float>(floatToken, molRes));
    REQUIRE(mol.getAtomPropIfPresent<float>(floatToken, 0, AtomRes));
    REQUIRE(mol.getBondPropIfPresent<float>(floatToken, 0, bondRes));

    CHECK(molRes == 3.0f);
    CHECK(AtomRes == 4.0f);
    CHECK(bondRes == 5.0f);
  }

  SECTION("Overrides existing prop of different type") {
    mol.addBondProp(floatToken, 3.0f);
    mol.addAtomProp(floatToken, 4.0f);

    mol.addBondProp(floatToken, 6u);
    mol.addAtomProp(floatToken, 6u);

    // These values are never used, but GCC can no longer deduce that they're
    // initialized, and maybe-uninitialized is currently treated as an error.
    unsigned int atomRes = 0;
    unsigned int bondRes = 0;
    REQUIRE(mol.getAtomPropIfPresent<unsigned int>(floatToken, 0, atomRes));
    REQUIRE(mol.getBondPropIfPresent<unsigned int>(floatToken, 0, bondRes));
    CHECK(atomRes == 6u);
    CHECK(bondRes == 6u);
  }

  SECTION("Clear Full properties") {
    REQUIRE(mol.addAtomProp(floatToken, 3.0f) != nullptr);
    REQUIRE(mol.addAtomProp(signedIntToken, 5) != nullptr);
    REQUIRE(mol.addBondProp(floatToken, 3.0f) != nullptr);

    mol.clearAtomPropIfPresent(floatToken);
    mol.clearBondPropIfPresent(floatToken);
    CHECK(!mol.getAtomPropArrayIfPresent<float>(floatToken));
    CHECK(!mol.getBondPropArrayIfPresent<int>(floatToken));

    int *res = mol.getAtomPropArrayIfPresent<int>(signedIntToken);
    REQUIRE(res != nullptr);
    CHECK(*res == 5);
  }

  SECTION("Clear atom properties") {
    REQUIRE(mol.addAtomProp(floatToken, 3.0f) != nullptr);
    mol.clearSingleAtomProp(floatToken, 1);
    float res;
    CHECK(mol.getAtomPropIfPresent<float>(floatToken, 0, res));
    CHECK(!mol.getAtomPropIfPresent<float>(floatToken, 1, res));
    CHECK(mol.getAtomPropIfPresent<float>(floatToken, 2, res));
  }

  SECTION("Clear bond properties") {
    REQUIRE(mol.addBondProp(floatToken, 3.0f) != nullptr);
    mol.clearSingleBondProp(floatToken, 1);
    float res;
    CHECK(mol.getBondPropIfPresent<float>(floatToken, 0, res));
    CHECK(!mol.getBondPropIfPresent<float>(floatToken, 1, res));
    CHECK(mol.getBondPropIfPresent<float>(floatToken, 2, res));
  }

  SECTION("Non scalar types non array accessible") {
    REQUIRE(mol.addAtomProp<std::vector<int>>(signedIntToken, {1, 2, 3}) ==
            nullptr);
    std::vector<int> res;
    CHECK(mol.getAtomPropIfPresent<std::vector<int>>(signedIntToken, 0, res));
    CHECK(res == std::vector<int>{1, 2, 3});
    CHECK(mol.getAtomPropArrayIfPresent<std::vector<int>>(signedIntToken) ==
          nullptr);
  }

  SECTION("Request scalar from nonscalar fails") {
    REQUIRE(mol.addAtomProp<std::vector<int>>(signedIntToken, {1, 2, 3}) ==
            nullptr);
    int res;
    CHECK_THROWS_AS(mol.getAtomPropIfPresent<int>(signedIntToken, 0, res),
                    std::bad_any_cast);
  }

  SECTION("Request non-scalar from scalar fails") {
    REQUIRE(mol.addAtomProp(signedIntToken, 3) != nullptr);
    std::vector<int> res;
    CHECK_THROWS_AS(
        mol.getAtomPropIfPresent<std::vector<int>>(signedIntToken, 0, res),
        std::bad_any_cast);
  }
}

static PropToken testToken("test");

template <typename initT, typename resT>
void checkCast(RDMol &mol) {
  std::string message =
      "Testing with type: " + std::string(typeid(initT).name()) + "->" +
      std::string(typeid(resT).name());
  INFO(message);

  // Check basic positive, negative, and overflow values
  constexpr std::array<int64_t, 3> testVals = {
      5, -10, int64_t(std::numeric_limits<int>::max()) + 1};
  for (auto &testVal : testVals) {
    auto typedVal = initT(testVal);
    mol.clearAtomPropIfPresent(testToken);
    mol.addAtomProp<initT>(testToken, typedVal);
    RDValue caster = typedVal;
    resT wantVal;
    bool expectThrowAnyCast = false;
    bool expectThrowBoostNumeric = false;
    try {
      wantVal = rdvalue_cast<resT>(caster);
    } catch (std::bad_any_cast &e) {
      expectThrowAnyCast = true;
    } catch (boost::bad_numeric_cast &e) {
      expectThrowBoostNumeric = true;
    }

    if (expectThrowAnyCast) {
      CHECK_THROWS_AS(mol.getAtomPropIfPresent<resT>(testToken, 0, wantVal),
                      std::bad_any_cast);
    } else if (expectThrowBoostNumeric) {
      CHECK_THROWS_AS(mol.getAtomPropIfPresent<resT>(testToken, 0, wantVal),
                      boost::bad_numeric_cast);
    } else {
      // This value is never used, but GCC can no longer deduce that it's
      // initialized, and maybe-uninitialized is currently treated as an error.
      resT gotVal = resT();
      REQUIRE(mol.getAtomPropIfPresent<resT>(testToken, 0, gotVal));
      CHECK(gotVal == wantVal);
    }
    RDValue::cleanup_rdvalue(caster);
  }
}

TEMPLATE_TEST_CASE("Property default type conversions", "", bool, char, int16_t,
                   int, unsigned int, int64_t, uint64_t, float, double) {
  RDMol mol;
  RDKit::SmilesParseTemp temp;
  static const char *basicSmiles = "CCCC";

  SmilesParserParams params;
  params.sanitize = false;
  params.removeHs = false;

  REQUIRE(RDKit::SmilesToMol(basicSmiles, params, mol, temp) == true);
  REQUIRE(mol.getNumAtoms() == 4);
  REQUIRE(mol.getNumBonds() == 3);
  SECTION("Cast to char") { checkCast<TestType, char>(mol); }
  SECTION("Cast to bool") { checkCast<TestType, bool>(mol); }
  SECTION("Cast to int16") { checkCast<TestType, int16_t>(mol); }
  SECTION("Cast to int") { checkCast<TestType, int>(mol); }
  SECTION("Cast to unsigned int") { checkCast<TestType, unsigned int>(mol); }
  SECTION("Cast to int64_t") { checkCast<TestType, int64_t>(mol); }
  SECTION("Cast to uint64_t") { checkCast<TestType, uint64_t>(mol); }
  SECTION("Cast to float") { checkCast<TestType, float>(mol); }
  SECTION("Cast to double") { checkCast<TestType, double>(mol); }
}

TEST_CASE("Bond stereo atoms") {
  auto molPtr = basicMol();
  RDMol &mol = *molPtr;

  SECTION("Basic set") {
    mol.setBondStereoAtoms(1, 1, 2);
    const auto *stereoAtoms = mol.getBondStereoAtoms(1);
    REQUIRE(stereoAtoms != nullptr);
    CHECK(stereoAtoms[0] == 1);
    CHECK(stereoAtoms[1] == 2);
    CHECK(mol.hasBondStereoAtoms(1));

    CHECK(!mol.hasBondStereoAtoms(2));
    auto *stereoAtoms2 = mol.getBondStereoAtoms(2);
    CHECK(stereoAtoms2 != nullptr);
    CHECK(stereoAtoms2[0] == atomindex_t(-1));
  }

  SECTION("Clear") {
    mol.setBondStereoAtoms(1, 1, 2);
    CHECK(mol.hasBondStereoAtoms(1));
    mol.clearBondStereoAtoms(1);
    const auto *stereoAtoms = mol.getBondStereoAtoms(1);
    CHECK(stereoAtoms != nullptr);
    CHECK(stereoAtoms[0] == atomindex_t(-1));
    CHECK(!mol.hasBondStereoAtoms(1));
  }

  SECTION("Get, populated from ROMol API") {
    ROMol &romol = mol.asROMol();
    romol.getBondWithIdx(1)->setStereoAtoms(0, 3);
    const auto *stereoAtoms = mol.getBondStereoAtoms(1);
    REQUIRE(stereoAtoms != nullptr);
    CHECK(stereoAtoms[0] == 0);
    CHECK(stereoAtoms[1] == 3);
    CHECK(mol.hasBondStereoAtoms(1));

    CHECK(!mol.hasBondStereoAtoms(2));
    auto *stereoAtoms2 = mol.getBondStereoAtoms(2);
    CHECK(stereoAtoms2 != nullptr);
    CHECK(stereoAtoms2[0] == atomindex_t(-1));
  }
}

TEST_CASE("RDMol add atom") {
  RDMol mol;
  RDKit::SmilesParseTemp temp;
  static const char *basicSmiles = "CNPO";
  SmilesParserParams params;
  params.sanitize = false;
  params.removeHs = false;

  REQUIRE(RDKit::SmilesToMol(basicSmiles, params, mol, temp) == true);
  REQUIRE(mol.getNumAtoms() == 4);
  REQUIRE(mol.getNumBonds() == 3);

  mol.setSingleAtomProp(signedIntToken, 1, 3);
  mol.setSingleBondProp(signedIntToken, 1,
                        4);  // Checking filter out of non-atom props
  REQUIRE(mol.hasAtomProp(signedIntToken, 1));

  std::vector<double> conf1(3 * mol.getNumAtoms(), 1.0);
  std::vector<double> conf2(3 * mol.getNumAtoms(), 2.0);
  mol.addConformer(conf1.data());
  mol.addConformer(conf2.data());

  AtomData &newAtom = mol.addAtom();
  SECTION("Default values and set") {
    CHECK(mol.getNumAtoms() == 5);
    CHECK(mol.getNumBonds() == 3);
    CHECK(newAtom.getAtomicNum() == 0);
    newAtom.setAtomicNum(6);

    AtomData &secondHandle = mol.getAtom(4);
    CHECK(secondHandle.getAtomicNum() == 6);

    const double *confPos = mol.getConformerPositions(0);
    REQUIRE(confPos != nullptr);
    const double *confPos2 = mol.getConformerPositions(1);
    REQUIRE(confPos2 != nullptr);
    for (size_t i = 0; i < 12; i++) {
      CHECK(confPos[i] == 1.0);
      CHECK(confPos2[i] == 2.0);
    }
    for (size_t i = 12; i < 15; i++) {
      CHECK(confPos[12] == 0.0);
      CHECK(confPos2[12] == 0.0);
    }
  }

  SECTION("Properties") {
    // Check that previous properties were not displaced
    CHECK(mol.hasAtomProp(signedIntToken, 1));
    CHECK(!mol.hasAtomProp(signedIntToken, 2));
    mol.setSingleAtomProp(signedIntToken, 4, 5);
    CHECK(mol.hasAtomProp(signedIntToken, 4));
  }
}

std::vector<uint32_t> iteratorToIndexVector(
    const std::pair<const uint32_t *, const uint32_t *> &iter) {
  std::vector<uint32_t> res;
  for (const uint32_t *idx = iter.first; idx != iter.second; idx++) {
    res.push_back(*idx);
  }
  return res;
}

TEST_CASE("Add Bond") {
  auto molPtr = basicMol();
  RDMol &mol = *molPtr;

  SECTION("Bond to self throws") {
    CHECK_THROWS_AS(mol.addBond(1, 1, BondEnums::BondType::SINGLE),
                    Invar::Invariant);
  }

  SECTION("Bond to out of bounds throws") {
    CHECK_THROWS_AS(mol.addBond(1, 5, BondEnums::BondType::SINGLE),
                    Invar::Invariant);
  }

  SECTION("Basic set") {
    BondData &newBond = mol.addBond(1, 3, BondEnums::BondType::DOUBLE);
    CHECK(mol.getNumBonds() == 4);
    CHECK(!newBond.getIsAromatic());
    CHECK(newBond.getBondType() == BondEnums::BondType::DOUBLE);
    newBond.setBondDir(BondEnums::BondDir::ENDUPRIGHT);

    BondData &sameBondNewRef = mol.getBond(3);
    CHECK(sameBondNewRef.getBondDir() == BondEnums::BondDir::ENDUPRIGHT);
  }

  SECTION("Iterate over bond atoms") {
    mol.addBond(1, 3, BondEnums::BondType::DOUBLE);

    auto atom1BondIterator = mol.getAtomBonds(1);
    auto atom2BondIterator = mol.getAtomBonds(2);
    auto atom3BondIterator = mol.getAtomBonds(3);

    std::vector<uint32_t> atom1BondIdxs =
        iteratorToIndexVector(atom1BondIterator);
    std::vector<uint32_t> atom2BondIdxs =
        iteratorToIndexVector(atom2BondIterator);
    std::vector<uint32_t> atom3BondIdxs =
        iteratorToIndexVector(atom3BondIterator);

    CHECK_THAT(atom1BondIdxs, Catch::Matchers::UnorderedEquals(
                                  std::vector<uint32_t>{0, 1, 3}));
    CHECK_THAT(atom2BondIdxs,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{1, 2}));
    CHECK_THAT(atom3BondIdxs,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{2, 3}));

    auto neighbors1 = mol.getAtomNeighbors(1);
    auto neighbors2 = mol.getAtomNeighbors(2);
    auto neighbors3 = mol.getAtomNeighbors(3);

    std::vector<uint32_t> neighbors1Idxs = iteratorToIndexVector(neighbors1);
    std::vector<uint32_t> neighbors2Idxs = iteratorToIndexVector(neighbors2);
    std::vector<uint32_t> neighbors3Idxs = iteratorToIndexVector(neighbors3);

    CHECK_THAT(neighbors1Idxs, Catch::Matchers::UnorderedEquals(
                                   std::vector<uint32_t>{0, 2, 3}));
    CHECK_THAT(neighbors2Idxs,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{1, 3}));
    CHECK_THAT(neighbors3Idxs,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{1, 2}));
  }

  SECTION("Set reverse order") {
    BondData &newBond = mol.addBond(3, 1, BondEnums::BondType::DOUBLE);
    CHECK(mol.getNumBonds() == 4);
    CHECK(!newBond.getIsAromatic());
    CHECK(newBond.getBondType() == BondEnums::BondType::DOUBLE);
    CHECK(newBond.getBeginAtomIdx() == 3);
    CHECK(newBond.getEndAtomIdx() == 1);
  }

  SECTION("Set aromatic") {
    AtomData &atom1 = mol.getAtom(1);
    AtomData &atom3 = mol.getAtom(3);
    CHECK(atom1.getIsAromatic() == false);
    CHECK(atom3.getIsAromatic() == false);
    BondData &newBond = mol.addBond(1, 3, BondEnums::BondType::AROMATIC);
    CHECK(newBond.getIsAromatic());
    CHECK(atom1.getIsAromatic());
    CHECK(atom3.getIsAromatic());
  }

  SECTION("Set aromatic, noAromaticity") {
    AtomData &atom1 = mol.getAtom(1);
    AtomData &atom3 = mol.getAtom(3);
    CHECK(atom1.getIsAromatic() == false);
    CHECK(atom3.getIsAromatic() == false);
    BondData &newBond = mol.addBond(1, 3, BondEnums::BondType::AROMATIC, false);
    CHECK(newBond.getIsAromatic());
    CHECK(!atom1.getIsAromatic());
    CHECK(!atom3.getIsAromatic());
  }

  SECTION("Properties") {
    mol.setSingleAtomProp(signedIntToken, 1,
                          5);  // Checking filtering out of non-bond type props
    mol.setSingleBondProp(signedIntToken, 1, 3);
    mol.setSingleBondProp(signedIntToken, 2, 4);
    mol.setSingleBondProp(complexToken, 0, std::vector<int>({1, 2, 3}));

    CHECK(mol.hasBondProp(signedIntToken, 1));
    CHECK(mol.hasBondProp(signedIntToken, 2));
    CHECK(!mol.hasBondProp(signedIntToken, 0));

    mol.addBond(1, 3, BondEnums::BondType::DOUBLE);
    CHECK(!mol.hasBondProp(signedIntToken, 3));
    CHECK(mol.hasBondProp(signedIntToken, 1));
    CHECK(mol.hasBondProp(signedIntToken, 2));
    CHECK(mol.hasBondProp(complexToken, 0));
    CHECK(!mol.hasBondProp(complexToken, 3));

    mol.setSingleBondProp(signedIntToken, 3, 5);
    CHECK(mol.hasBondProp(signedIntToken, 3));
  }
}

TEST_CASE("Atom valence") {
  SECTION("Neither implicit nor explicit set") {
    RDMol mol;
    mol.addAtom().setAtomicNum(6);
    CHECK_THROWS_AS(mol.getAtom(0).getExplicitValence(), Invar::Invariant);
    CHECK_THROWS_AS(mol.getAtom(0).getImplicitValence(), Invar::Invariant);
    CHECK(mol.getAtom(0).needsUpdatePropertyCache());
  }

  SECTION("Explicit valence set") {
    RDMol mol;
    mol.addAtom().setAtomicNum(6);
    mol.calcAtomExplicitValence(0, false);
    CHECK(mol.getAtom(0).getExplicitValence() == 0);
    CHECK_THROWS_AS(mol.getAtom(0).getImplicitValence(), Invar::Invariant);

    // Still needs implicit.
    CHECK(mol.getAtom(0).needsUpdatePropertyCache());
  }

  SECTION("No implicit set") {
    RDMol mol;
    AtomData &atom = mol.addAtom();
    atom.setAtomicNum(6);
    atom.setNoImplicit(true);
    CHECK(atom.getImplicitValence() == 0);

    // Check that no implicit means no need to update implicit.
    CHECK(atom.needsUpdatePropertyCache());
    mol.calcAtomExplicitValence(0, false);
    CHECK(!atom.needsUpdatePropertyCache());
  }

  SECTION("Correct after calculation") {
    auto molPtr = basicMol();
    RDMol &mol = *molPtr;

    mol.calcAtomImplicitValence(0, false);
    mol.calcAtomImplicitValence(1, false);
    mol.calcAtomImplicitValence(2, false);
    mol.calcAtomImplicitValence(3, false);

    CHECK(!mol.getAtom(0).needsUpdatePropertyCache());

    CHECK(mol.getAtom(0).getExplicitValence() == 1);
    CHECK(mol.getAtom(1).getExplicitValence() == 2);
    CHECK(mol.getAtom(2).getExplicitValence() == 3);  // Double bond to 3
    CHECK(mol.getAtom(3).getExplicitValence() == 2);

    CHECK(mol.getAtom(0).getImplicitValence() == 3);
    CHECK(mol.getAtom(1).getImplicitValence() == 2);
    CHECK(mol.getAtom(2).getImplicitValence() == 1);
    CHECK(mol.getAtom(3).getImplicitValence() == 2);
  }
}

TEST_CASE("Conformers") {
  auto molPtr = basicMol();
  RDMol &mol = *molPtr;
  const int numDoublesInConf = 3 * mol.getNumAtoms();

  std::vector<double> conf0(numDoublesInConf, 1.0);
  std::vector<double> conf1(numDoublesInConf, 2.0);
  std::vector<double> conf2(numDoublesInConf, 3.0);

  SECTION("No conformers") {
    CHECK(mol.getNumConformers() == 0);
    CHECK_THROWS_AS(mol.getConformerPositions(), ConformerException);
    CHECK_THROWS_AS(mol.getConformerPositions(3), ConformerException);
  }

  SECTION("Access missing conformer") {
    mol.addConformer(conf1.data());
    CHECK_THROWS_AS(mol.getConformerPositions(15), ConformerException);
    CHECK_THROWS_AS(mol.is3DConformer(15), ConformerException);
  }

  SECTION("Add in order, and get") {
    uint32_t id1 = mol.addConformer(conf1.data());
    uint32_t id2 = mol.addConformer(conf0.data());
    uint32_t id3 = mol.addConformer(conf2.data());

    // Check that reserving is OK
    mol.allocateConformers(10);

    CHECK(mol.getNumConformers() == 3);
    CHECK(id1 == 0);
    CHECK(id2 == 1);
    CHECK(id3 == 2);

    double *res1 = mol.getConformerPositions(id1);
    CHECK(res1[0] == 2.0);
    CHECK(res1[numDoublesInConf - 1] == 2.0);
    double *res2 = mol.getConformerPositions(id2);
    CHECK(res2[0] == 1.0);
    CHECK(res2[numDoublesInConf - 1] == 1.0);
    double *res3 = mol.getConformerPositions(id3);
    CHECK(res3[0] == 3.0);
    CHECK(res3[numDoublesInConf - 1] == 3.0);

    CHECK(mol.is3DConformer(id1));
    CHECK(mol.is3DConformer(id2));
    CHECK(mol.is3DConformer(id3));
  }

  SECTION("Get default conformer") {
    mol.addConformer(conf1.data());
    mol.addConformer(conf0.data());
    double *res1 = mol.getConformerPositions();
    CHECK(res1[0] == 2.0);
  }

  SECTION("Add in order, then add one non3D") {
    uint32_t id1 = mol.addConformer(conf1.data());
    uint32_t id2 = mol.addConformer(conf0.data());
    uint32_t id3 = mol.addConformer(conf2.data(), -1, false);

    CHECK(mol.getNumConformers() == 3);
    double *res1 = mol.getConformerPositions(id1);
    CHECK(res1[0] == 2.0);
    CHECK(res1[numDoublesInConf - 1] == 2.0);
    double *res2 = mol.getConformerPositions(id2);
    CHECK(res2[0] == 1.0);
    CHECK(res2[numDoublesInConf - 1] == 1.0);
    const double *res3 = mol.getConformerPositions(id3);
    CHECK(res3[0] == 3.0);
    CHECK(res3[numDoublesInConf - 1] == 3.0);

    CHECK(mol.is3DConformer(id1));
    CHECK(mol.is3DConformer(id2));
    CHECK(!mol.is3DConformer(id3));
  }

  SECTION("Add in order, then add one with a different ID") {
    uint32_t id1 = mol.addConformer(conf1.data());
    uint32_t id2 = mol.addConformer(conf0.data());
    uint32_t id3 = mol.addConformer(conf2.data(), 42);

    CHECK(id3 == 42);

    CHECK(mol.getNumConformers() == 3);
    double *res1 = mol.getConformerPositions(id1);
    CHECK(res1[0] == 2.0);
    CHECK(res1[numDoublesInConf - 1] == 2.0);
    double *res2 = mol.getConformerPositions(id2);
    CHECK(res2[0] == 1.0);
    CHECK(res2[numDoublesInConf - 1] == 1.0);
    const double *res3 = mol.getConformerPositions(id3);
    CHECK(res3[0] == 3.0);
    CHECK(res3[numDoublesInConf - 1] == 3.0);

    CHECK(mol.is3DConformer(id1));
    CHECK(mol.is3DConformer(id2));
    CHECK(mol.is3DConformer(id3));
  }

  SECTION("Add with an ID, then add others with default") {
    uint32_t id1 = mol.addConformer(conf1.data(), 42);
    uint32_t id2 = mol.addConformer(conf0.data(), -1, false);
    uint32_t id3 = mol.addConformer(conf2.data());

    CHECK(id1 == 42);

    CHECK(mol.getNumConformers() == 3);
    double *res1 = mol.getConformerPositions(id1);
    CHECK(res1[0] == 2.0);
    CHECK(res1[numDoublesInConf - 1] == 2.0);
    double *res2 = mol.getConformerPositions(id2);
    CHECK(res2[0] == 1.0);
    CHECK(res2[numDoublesInConf - 1] == 1.0);
    const double *res3 = mol.getConformerPositions(id3);
    CHECK(res3[0] == 3.0);
    CHECK(res3[numDoublesInConf - 1] == 3.0);

    CHECK(mol.is3DConformer(id1));
    CHECK(!mol.is3DConformer(id2));
    CHECK(mol.is3DConformer(id3));
  }

  SECTION("Clear all conformers") {
    mol.addConformer(conf1.data(), 42);
    mol.addConformer(conf0.data(), -1, false);
    mol.addConformer(conf2.data());

    CHECK(mol.getNumConformers() == 3);
    mol.clearConformers();
    CHECK(mol.getNumConformers() == 0);
  }

  SECTION("Clear one conformer, previously set IDs") {
    uint32_t id1 = mol.addConformer(conf1.data(), 42);
    uint32_t id2 = mol.addConformer(conf0.data(), -1, false);
    uint32_t id3 = mol.addConformer(conf2.data());

    CHECK(mol.getNumConformers() == 3);
    mol.removeConformer(id1);
    CHECK(mol.getNumConformers() == 2);

    CHECK_THROWS_AS(mol.getConformerPositions(id1), ConformerException);
    double *res2 = mol.getConformerPositions(id2);
    CHECK(res2[0] == 1.0);
    CHECK(res2[numDoublesInConf - 1] == 1.0);
    double *res3 = mol.getConformerPositions(id3);
    CHECK(res3[0] == 3.0);
    CHECK(res3[numDoublesInConf - 1] == 3.0);

    CHECK_THROWS_AS(mol.is3DConformer(id1), ConformerException);
    CHECK(!mol.is3DConformer(id2));
    CHECK(mol.is3DConformer(id3));
  }

  SECTION("Clear one conformer, no previously set IDs") {
    uint32_t id1 = mol.addConformer(conf1.data());
    uint32_t id2 = mol.addConformer(conf0.data());
    uint32_t id3 = mol.addConformer(conf2.data());

    CHECK(mol.getNumConformers() == 3);
    mol.removeConformer(id1);
    CHECK(mol.getNumConformers() == 2);

    CHECK_THROWS_AS(mol.getConformerPositions(id1), ConformerException);
    double *res2 = mol.getConformerPositions(id2);
    CHECK(res2[0] == 1.0);
    CHECK(res2[numDoublesInConf - 1] == 1.0);
    double *res3 = mol.getConformerPositions(id3);
    CHECK(res3[0] == 3.0);
    CHECK(res3[numDoublesInConf - 1] == 3.0);

    CHECK_THROWS_AS(mol.is3DConformer(id1), ConformerException);
    CHECK(mol.is3DConformer(id2));
    CHECK(mol.is3DConformer(id3));
  }

  SECTION("Remove nonexistent conformer does nothing") {
    mol.addConformer(conf1.data());
    mol.addConformer(conf0.data());

    mol.removeConformer(50);
    CHECK(mol.getNumConformers() == 2);
  }

  SECTION("Valid behavior with no atoms") {
    RDMol emptyMol;
    CHECK(emptyMol.getNumConformers() == 0);
    CHECK_THROWS_AS(emptyMol.getConformerPositions(), ConformerException);
    CHECK_THROWS_AS(emptyMol.getConformerPositions(3), ConformerException);
    CHECK(emptyMol.addConformer(nullptr) == 0);
    CHECK(emptyMol.addConformer(nullptr, 5) == 5);
    CHECK(emptyMol.addConformer(nullptr) == 6);
    CHECK(emptyMol.getNumConformers() == 3);
    CHECK_NOTHROW(emptyMol.getConformerPositions());
    CHECK_NOTHROW(emptyMol.getConformerPositions(5));
    CHECK(emptyMol.is3DConformer(0) == true);
    emptyMol.removeConformer(5);
    CHECK(emptyMol.getNumConformers() == 2);
    emptyMol.clearConformers();
    CHECK(emptyMol.getNumConformers() == 0);
  }

  SECTION("Multiple conformers same ID valid") {
    uint32_t id1 = mol.addConformer(conf1.data(), 3, false);
    uint32_t id2 = mol.addConformer(conf0.data(), 3, true);
    CHECK(mol.getNumConformers() == 2);
    CHECK(id1 == id2);
    double *res1 = mol.getConformerPositions(id1);
    CHECK(res1[0] == 2.0);
    CHECK(!mol.is3DConformer(id1));
  }
}

TEST_CASE("Remove bond") {
  auto molPtr = basicMol();
  RDMol &mol = *molPtr;

  SECTION("Nullop") {
    mol.removeBond(5);
    CHECK(mol.getNumBonds() == 3);
  }

  SECTION("Basic remove") {
    int atomIdx1, atomIdx2;
    {
      BondData &bondToRemove = mol.getBond(1);
      atomIdx1 = bondToRemove.getBeginAtomIdx();
      atomIdx2 = bondToRemove.getEndAtomIdx();
    }

    uint32_t degree1 = mol.getAtomDegree(atomIdx1);
    uint32_t degree2 = mol.getAtomDegree(atomIdx2);

    mol.removeBond(1);
    CHECK(mol.getNumBonds() == 2);
    // Check that the index has shifted down.
    CHECK(mol.getBond(0).getBondType() == BondEnums::BondType::SINGLE);
    CHECK(mol.getBond(1).getBondType() == BondEnums::BondType::DOUBLE);

    CHECK(mol.getAtomDegree(atomIdx1) == degree1 - 1);
    CHECK(mol.getAtomDegree(atomIdx2) == degree2 - 1);

    auto neighbors0 = mol.getAtomNeighbors(0);
    auto neighbors1 = mol.getAtomNeighbors(1);
    auto neighbors2 = mol.getAtomNeighbors(2);
    auto neighbors3 = mol.getAtomNeighbors(3);

    std::vector<uint32_t> neighbors0Idxs = iteratorToIndexVector(neighbors0);
    std::vector<uint32_t> neighbors1Idxs = iteratorToIndexVector(neighbors1);
    std::vector<uint32_t> neighbors2Idxs = iteratorToIndexVector(neighbors2);
    std::vector<uint32_t> neighbors3Idxs = iteratorToIndexVector(neighbors3);

    CHECK_THAT(neighbors0Idxs,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{1}));
    CHECK_THAT(neighbors1Idxs,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{0}));
    CHECK_THAT(neighbors2Idxs,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{3}));
    CHECK_THAT(neighbors3Idxs,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{2}));
  }

  SECTION("Properties") {
    mol.setSingleBondProp(signedIntToken, 1, 3);
    mol.setSingleBondProp(signedIntToken, 2, 4);

    mol.removeBond(1);

    // Check that we shifted down, and index 1 now points to what was previously
    // index 2.
    CHECK(mol.getNumBonds() == 2);
    // This value is never used, but GCC can no longer deduce that it's
    // initialized, and maybe-uninitialized is currently treated as an error.
    int res = 0;
    CHECK(mol.getBondPropIfPresent(signedIntToken, 1, res));
    CHECK(res == 4);
  }

  SECTION("Bookmarks") {
    mol.setBondBookmark(0, 4);
    mol.setBondBookmark(0, 3);
    mol.setBondBookmark(1, 0);
    mol.setBondBookmark(1, 3);
    mol.setBondBookmark(1, 4);
    mol.setBondBookmark(2, 0);
    mol.setBondBookmark(2, 3);
    mol.removeBond(1);
    // The bookmark for bond 0 should be the same.
    CHECK(mol.getBondWithBookmark(4) == 0);
    // The bookmark for bond 2 should now be 1.
    CHECK(mol.getBondWithBookmark(0) == 1);

    std::vector<int> res = mol.getAllBondsWithBookmarks(3);
    // The previous bond 1 bookmark is gone, but bond 2 shifts to 1 now.
    CHECK(res == std::vector<int>{0, 1});

    // Only atom 0 should still be here.
    res = mol.getAllBondsWithBookmarks(4);
    CHECK(res == std::vector<int>{0});
  }

  SECTION("Substance groups") {
    auto &sgs = mol.getSubstanceGroups();
    REQUIRE(sgs.size() == 0);

    ROMol &romol = mol.asROMol();

    SubstanceGroup &sg = sgs.emplace_back(&romol, "A");
    sg.addAtomWithIdx(0);
    sg.addBondWithIdx(1);
    sg.addBondWithIdx(2);
    sg.setProp("index", 0);
    SubstanceGroup &sg2 = sgs.emplace_back(&romol, "B");
    sg2.addAtomWithIdx(1);
    sg2.addBondWithIdx(0);
    sg2.setProp("PARENT", 0);
    SubstanceGroup &sg3 = sgs.emplace_back(&romol, "C");
    sg3.addAtomWithIdx(2);
    sg3.addBondWithIdx(2);

    mol.removeBond(1);

    CHECK(sgs.size() == 1);
    // SG 0 should be removed. Since it is the parent of 1, so is 1. 2 Remains
    // shifted down.
    CHECK(sgs[0].getAtoms() == std::vector<unsigned int>{2});
    CHECK(sgs[0].getBonds() == std::vector<unsigned int>{1});
  }
}

TEST_CASE("Remove bond with stereo") {
  // RDKit::SmilesParseTemp temp;

  // The middle bond is trans, so the stereo atoms are 0 and 3, on the far ends
  // of the double bond respectively. The other 2 ends of the 4 bonds that
  // determine the cis/trans stereo are hydrogens.
  static const char *basicSmiles = "C/C=C/C";

  RDKit::v2::SmilesParse::SmilesParserParams params;
  params.sanitize = true;
  params.removeHs = true;

  // REQUIRE(RDKit::SmilesToMol(basicSmiles, params, mol, temp) == true);
  std::unique_ptr<RWMol> molOwner =
      RDKit::v2::SmilesParse::MolFromSmiles(basicSmiles, params);
  RDMol &mol = molOwner->asRDMol();
  REQUIRE(mol.getNumAtoms() == 4);
  REQUIRE(mol.getNumBonds() == 3);

  auto &doubleBond = mol.getBond(1);
  doubleBond.setStereo(BondEnums::BondStereo::STEREOTRANS);
  mol.setBondStereoAtoms(1, 0, 3);
  REQUIRE(mol.getBond(1).getStereo() == BondEnums::BondStereo::STEREOTRANS);
  REQUIRE(mol.getBondStereoAtoms(1)[0] == 0);
  REQUIRE(mol.getBondStereoAtoms(1)[1] == 3);
  REQUIRE(mol.hasBondStereoAtoms(1));

  mol.removeBond(0);
  CHECK(mol.getNumBonds() == 2);
  CHECK(mol.getBond(1).getStereo() == BondEnums::BondStereo::STEREONONE);
  CHECK(!mol.hasBondStereoAtoms(1));
}

TEST_CASE("Remove Atom") {
  auto molPtr = basicMol();
  RDMol &mol = *molPtr;

  SECTION("Basic, no bond removal, no stereochem") {
    // Change up atom numbers for better testing.
    mol.getAtom(1).setAtomicNum(7);
    mol.getAtom(2).setAtomicNum(8);

    REQUIRE(mol.getNumAtoms() == 4);
    REQUIRE(mol.getNumBonds() == 3);
    REQUIRE(mol.getAtom(0).getAtomicNum() == 6);
    REQUIRE(mol.getAtom(1).getAtomicNum() == 7);
    REQUIRE(mol.getBond(0).getBeginAtomIdx() == 0);
    REQUIRE(mol.getBond(0).getEndAtomIdx() == 1);
    mol.removeAtom(0);

    CHECK(mol.getNumAtoms() == 3);
    CHECK(mol.getNumBonds() == 2);

    CHECK(mol.getAtom(0).getAtomicNum() == 7);
    CHECK(mol.getAtom(1).getAtomicNum() == 8);

    // Check that bond atom refs are updated.
    CHECK(mol.getBond(0).getBeginAtomIdx() == 0);
    CHECK(mol.getBond(0).getEndAtomIdx() == 1);

    // Check atom neighbor iteration unchanged except for indices
    auto neighbors0 = mol.getAtomNeighbors(0);
    auto neighbors1 = mol.getAtomNeighbors(1);
    auto neighbors2 = mol.getAtomNeighbors(2);

    std::vector<uint32_t> neighbors0Idxs = iteratorToIndexVector(neighbors0);
    std::vector<uint32_t> neighbors1Idxs = iteratorToIndexVector(neighbors1);
    std::vector<uint32_t> neighbors2Idxs = iteratorToIndexVector(neighbors2);

    CHECK_THAT(neighbors0Idxs,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{1}));
    CHECK_THAT(neighbors1Idxs,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{0, 2}));
    CHECK_THAT(neighbors2Idxs,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{1}));
  }

  SECTION("Removes associated bonds") {
    REQUIRE(mol.getNumAtoms() == 4);
    REQUIRE(mol.getNumBonds() == 3);

    mol.removeAtom(1);
    // C is now disconnected, O and F are still connected.

    CHECK(mol.getNumAtoms() == 3);
    CHECK(mol.getNumBonds() == 1);

    auto neighbors0 = mol.getAtomNeighbors(0);
    auto neighbors1 = mol.getAtomNeighbors(1);
    auto neighbors2 = mol.getAtomNeighbors(2);

    std::vector<uint32_t> neighbors0Idxs = iteratorToIndexVector(neighbors0);
    std::vector<uint32_t> neighbors1Idxs = iteratorToIndexVector(neighbors1);
    std::vector<uint32_t> neighbors2Idxs = iteratorToIndexVector(neighbors2);

    CHECK_THAT(neighbors0Idxs,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{}));
    CHECK_THAT(neighbors1Idxs,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{2}));
    CHECK_THAT(neighbors2Idxs,
               Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{1}));
  }

  SECTION("Stereo groups") {
    auto groups = std::make_unique<StereoGroups>();
    groups->addGroup(RDKit::StereoGroupType::STEREO_ABSOLUTE,
                     std::vector<uint32_t>{1}, std::vector<uint32_t>{});
    groups->addGroup(RDKit::StereoGroupType::STEREO_OR,
                     std::vector<uint32_t>{0, 1, 2}, std::vector<uint32_t>{},
                     /*readId=*/5);
    groups->addGroup(RDKit::StereoGroupType::STEREO_AND,
                     std::vector<uint32_t>{2, 3}, std::vector<uint32_t>{});
    mol.setStereoGroups(std::move(groups));

    mol.removeAtom(1);
    StereoGroups *res = mol.getStereoGroups();
    REQUIRE(res != nullptr);
    REQUIRE(res->getNumGroups() == 2);
    auto group1Atoms = iteratorToIndexVector(res->getAtoms(0));
    auto group2Atoms = iteratorToIndexVector(res->getAtoms(1));

    CHECK(group1Atoms == std::vector<uint32_t>{0, 1});
    CHECK(group2Atoms == std::vector<uint32_t>{1, 2});
  }

  SECTION("Bookmarks") {
    mol.setAtomBookmark(0, 0);
    mol.setAtomBookmark(1, 0);
    mol.setAtomBookmark(2, 0);
    mol.setAtomBookmark(3, 0);
    mol.setAtomBookmark(1, 1);
    mol.setAtomBookmark(2, 1);
    mol.setAtomBookmark(3, 2);

    mol.setBondBookmark(0, 0);
    mol.setBondBookmark(1, 0);
    mol.setBondBookmark(2, 0);
    mol.setBondBookmark(0, 1);
    mol.setBondBookmark(1, 1);
    mol.setBondBookmark(2, 2);

    // Before deletion, we have {bookmark: indices}
    // atoms:
    // 0: 0, 1, 2, 3
    // 1: 1, 2
    // 2: 3
    // bonds:
    // 0: 0, 1, 2
    // 1: 0, 1
    // 2: 2

    // We are deleting bonds 0 and 1, atom 1. Bond 2 now maps to 0. Atoms 2 and
    // 3 now map to 1 and 2.
    mol.removeAtom(1);

    std::vector<int> atomBookmarks0 = mol.getAllAtomsWithBookmarks(0);
    std::vector<int> atomBookmarks1 = mol.getAllAtomsWithBookmarks(1);
    std::vector<int> atomBookmarks2 = mol.getAllAtomsWithBookmarks(2);

    std::vector<int> bondBookmarks0 = mol.getAllBondsWithBookmarks(0);
    std::vector<int> bondBookmarks1 = mol.getAllBondsWithBookmarks(1);
    std::vector<int> bondBookmarks2 = mol.getAllBondsWithBookmarks(2);

    CHECK(atomBookmarks0 == std::vector<int>{0, 1, 2});
    CHECK(atomBookmarks1 == std::vector<int>{1});
    CHECK(atomBookmarks2 == std::vector<int>{2});

    CHECK(bondBookmarks0 == std::vector<int>{0});
    CHECK(bondBookmarks1 == std::vector<int>{});
    CHECK(bondBookmarks2 == std::vector<int>{0});
  }

  SECTION("Properties") {
    {
      int *atomProps = mol.addAtomProp(signedIntToken, 0);
      std::iota(atomProps, atomProps + mol.getNumAtoms(), 0);
      int *bondProps = mol.addBondProp(signedIntToken, 0);
      std::iota(bondProps, bondProps + mol.getNumBonds(), 0);
    }
    mol.removeAtom(1);
    int *atomProps = mol.getAtomPropArrayIfPresent<int>(signedIntToken);
    int *bondProps = mol.getBondPropArrayIfPresent<int>(signedIntToken);
    CHECK(atomProps[0] == 0);
    CHECK(atomProps[1] == 2);
    CHECK(atomProps[2] == 3);
    CHECK(bondProps[0] == 2);
  }

  SECTION("Conformers") {
    const int numDoublesInConf = 3 * mol.getNumAtoms();

    std::vector<double> conf0(numDoublesInConf, 1.0);
    std::vector<double> conf1(numDoublesInConf, 2.0);
    std::vector<double> conf2(numDoublesInConf, 3.0);

    uint32_t id1 = mol.addConformer(conf0.data());
    uint32_t id2 = mol.addConformer(conf1.data(), 3);
    uint32_t id3 = mol.addConformer(conf2.data(), -1, false);

    mol.removeAtom(1);
    const int newNumDoublesInConf = 3 * mol.getNumAtoms();

    CHECK(mol.getNumConformers() == 3);
    double *res1 = mol.getConformerPositions(id1);
    CHECK(res1[0] == 1.0);
    CHECK(res1[newNumDoublesInConf - 1] == 1.0);
    double *res2 = mol.getConformerPositions(id2);
    CHECK(res2[0] == 2.0);
    CHECK(res2[newNumDoublesInConf - 1] == 2.0);
    const double *res3 = mol.getConformerPositions(id3);
    CHECK(res3[0] == 3.0);
    CHECK(res3[newNumDoublesInConf - 1] == 3.0);

    CHECK(mol.is3DConformer(id1));
    CHECK(mol.is3DConformer(id2));
    CHECK(!mol.is3DConformer(id3));
  }
}

void checkBasicEdit(const RDMol &mol) {
  CHECK(mol.getNumAtoms() == 3);
  CHECK(mol.getNumBonds() == 1);

  // Check that the atom and bond indices have shifted down.
  CHECK(mol.getAtom(0).getAtomicNum() == 6);
  CHECK(mol.getAtom(1).getAtomicNum() == 6);
  CHECK(mol.getAtom(2).getAtomicNum() == 7);

  CHECK(mol.getBond(0).getBondType() == BondEnums::BondType::DOUBLE);

  // check neighbors
  auto neighbors0 = mol.getAtomNeighbors(0);
  auto neighbors1 = mol.getAtomNeighbors(1);
  auto neighbors2 = mol.getAtomNeighbors(2);

  std::vector<uint32_t> neighbors0Idxs = iteratorToIndexVector(neighbors0);
  std::vector<uint32_t> neighbors1Idxs = iteratorToIndexVector(neighbors1);
  std::vector<uint32_t> neighbors2Idxs = iteratorToIndexVector(neighbors2);

  CHECK_THAT(neighbors0Idxs,
             Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{}));
  CHECK_THAT(neighbors1Idxs,
             Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{2}));
  CHECK_THAT(neighbors2Idxs,
             Catch::Matchers::UnorderedEquals(std::vector<uint32_t>{1}));
}

TEST_CASE("Batch edits") {
  auto molPtr = parseSmiles("OCC=NC");
  RDMol &mol = *molPtr;
  mol.setAtomBookmark(0, 0);
  mol.setAtomBookmark(1, 1);
  mol.setBondBookmark(0, 0);
  mol.setBondBookmark(3, 0);
  mol.setSingleAtomProp(signedIntToken, 0, 1);
  mol.setSingleAtomProp(signedIntToken, 1, 2);
  mol.setSingleAtomProp(signedIntToken, 4, 3);
  mol.setSingleBondProp(signedIntToken, 0, 3);
  mol.setSingleBondProp(signedIntToken, 2, 4);

  SECTION("Empty") {
    mol.beginBatchEdit();
    mol.commitBatchEdit();
    CHECK(mol.getNumAtoms() == 5);
    CHECK(mol.getNumBonds() == 4);
  }

  SECTION("Edit empty mol") {
    RDMol empty;
    empty.beginBatchEdit();
    empty.commitBatchEdit();
    CHECK(empty.getNumAtoms() == 0);
    CHECK(empty.getNumBonds() == 0);
  }

  SECTION("Basics") {
    mol.beginBatchEdit();
    mol.removeAtom(0);
    mol.removeBond(1);
    mol.removeAtom(
        4);  // Note that if incorrectly set, this would be past nAtoms.
    // Include that atoms and bonds aren't deleted until end of edit.

    CHECK(mol.getNumAtoms() == 5);
    CHECK(mol.getNumBonds() == 4);
    mol.commitBatchEdit();

    checkBasicEdit(mol);

    // Include that after finishing, we can begin again.
    mol.beginBatchEdit();
  }

  SECTION("commit does nothing if no begin") {
    mol.commitBatchEdit();
    CHECK(mol.getNumAtoms() == 5);
    CHECK(mol.getNumBonds() == 4);
  }

  SECTION("Can't begin twice") {
    mol.beginBatchEdit();
    CHECK_THROWS_AS(mol.beginBatchEdit(), ValueErrorException);
  }

  SECTION("Remove all atoms and bonds") {
    mol.beginBatchEdit();
    mol.removeAtom(0);
    mol.removeAtom(1);
    mol.removeAtom(2);
    mol.removeAtom(3);
    mol.removeAtom(4);
    mol.removeBond(0);
    mol.removeBond(1);
    mol.removeBond(2);
    mol.removeBond(3);
    mol.commitBatchEdit();

    CHECK(mol.getNumAtoms() == 0);
    CHECK(mol.getNumBonds() == 0);
  }

  SECTION("Adding atoms/bonds during edit") {
    mol.beginBatchEdit();
    mol.addAtom().setAtomicNum(5);
    mol.addBond(2, 4, BondEnums::BondType::AROMATIC);
    mol.removeAtom(3);
    mol.removeBond(0);
    mol.addAtom().setAtomicNum(9);
    mol.removeAtom(5);

    mol.commitBatchEdit();

    CHECK(mol.getNumAtoms() == 5);
    CHECK(mol.getNumBonds() == 2);

    // OCC=NC 9
    CHECK(mol.getAtom(0).getAtomicNum() == 8);
    CHECK(mol.getAtom(1).getAtomicNum() == 6);
    CHECK(mol.getAtom(2).getAtomicNum() == 6);
    CHECK(mol.getAtom(3).getAtomicNum() == 6);
    CHECK(mol.getAtom(4).getAtomicNum() == 9);

    CHECK(mol.getBond(0).getBondType() == BondEnums::BondType::SINGLE);
    CHECK(mol.getBond(1).getBondType() == BondEnums::BondType::AROMATIC);
  }

  SECTION("Bookmarks updated") {
    mol.beginBatchEdit();
    mol.removeAtom(0);
    mol.removeBond(1);
    mol.removeAtom(4);
    mol.commitBatchEdit();

    CHECK(!mol.hasAtomBookmark(0));
    CHECK(mol.getAtomWithBookmark(1) == 0);
    CHECK(!mol.hasBondBookmark(0));
  }

  SECTION("Properties updated") {
    mol.beginBatchEdit();
    mol.removeAtom(0);
    mol.removeBond(1);
    mol.removeAtom(4);
    mol.commitBatchEdit();

    int res;
    CHECK(mol.getAtomPropIfPresent(signedIntToken, 0, res));
    CHECK(res == 2);
    CHECK(mol.getBondPropIfPresent(signedIntToken, 0, res));
    CHECK(res == 4);
  }

  SECTION("Copy during edit") {
    mol.beginBatchEdit();
    mol.removeAtom(0);
    mol.removeBond(1);
    mol.removeAtom(4);

    CHECK(mol.getNumAtoms() == 5);
    CHECK(mol.getNumBonds() == 4);

    RDMol mol2(mol);

    mol2.commitBatchEdit();
    checkBasicEdit(mol2);
  }

  SECTION("Move during edit") {
    mol.beginBatchEdit();
    mol.removeAtom(0);
    mol.removeBond(1);
    mol.removeAtom(4);

    CHECK(mol.getNumAtoms() == 5);
    CHECK(mol.getNumBonds() == 4);

    RDMol mol2(std::move(mol));

    mol2.commitBatchEdit();
    checkBasicEdit(mol2);
  }

  SECTION("Rollback") {
    mol.beginBatchEdit();
    mol.removeAtom(0);
    mol.removeBond(1);
    mol.removeAtom(4);

    mol.rollbackBatchEdit();
    mol.commitBatchEdit();

    CHECK(mol.getNumAtoms() == 5);
    CHECK(mol.getNumBonds() == 4);
  }

  SECTION("Double delete is nullop") {
    mol.beginBatchEdit();
    mol.removeAtom(0);
    mol.removeAtom(0);
    mol.commitBatchEdit();
    CHECK(mol.getNumAtoms() == 4);
    CHECK(mol.getNumBonds() == 3);
  }

  SECTION("Delete end atom is safe after previous deletion") {
    mol.beginBatchEdit();
    mol.removeAtom(3);
    mol.removeAtom(4);
    mol.commitBatchEdit();
    CHECK(mol.getNumAtoms() == 3);
    CHECK(mol.getNumBonds() == 2);
  }
}

TEST_CASE("benchmarking construction from smiles") {
  // from hanoittest.cpp
  std::vector<std::string> smileses = {
      "C[C@@H]1CCC[C@H](C)[C@H]1C",
      "N[C@@]1(C[C@H]([18F])C1)C(=O)O",
      "CC12CCCC1C1CCC3CC(O)CCC3(C)C1CC2",
      "CC(C)CCCC[C@@H]1C[C@H](/C=C/[C@]2(C)CC[C@H](O)CC2)[C@@H](O)[C@H]1O",
      "C[C@@]12CCC[C@H]1[C@@H]1CC[C@H]3C[C@@H](O)CC[C@]3(C)[C@H]1CC2",
      "CCCN[C@H]1CC[C@H](NC)CC1",
      "O=S(=O)(NC[C@H]1CC[C@H](CNCc2ccc3ccccc3c2)CC1)c1ccc2ccccc2c1",
      "CC(C)[C@H]1CC[C@H](C(=O)N[C@H](Cc2ccccc2)C(=O)O)CC1",
      "O=[N+]([O-])c1ccccc1S(=O)(=O)NC[C@H]1CC[C@H](CNCC2Cc3ccccc3CC2)CC1",
      "Oc1ccc2c(Cc3ccc(OCCN4CCCCC4)cc3)c([C@H]3CC[C@H](O)CC3)sc2c1",
      "O=C(c1ccc(OCCN2CCCCC2)cc1)c1c2ccc(O)cc2sc1[C@H]1CC[C@H](O)CC1",
      "N#Cc1ccc2c(c1)CCN(CC[C@@H]1CC[C@@H](NC(=O)c3ccnc4ccccc34)CC1)C2",
      "COCCOC[C@H](CC1(C(=O)N[C@H]2CC[C@@H](C(=O)O)CC2)CCCC1)C(=O)O",
      "c1ccc(CN[C@H]2CC[C@H](Nc3ccc4[nH]ncc4c3)CC2)cc1",
      "[CH]12=[CH]3[CH]4=[CH]5[CH-]16.[Fe]23456",
      "[CH]12=[CH]3[CH]4=[CH]5[CH]6=[CH]17.[Fe]234567",
      "[cH]12[cH]3[cH]4[cH]5[cH]6[cH]17.[Fe]234567",
      "[cH]12[cH]3[cH]4[cH]5[nH]16.[Fe]23456",
      "[CH]12=[CH]3[CH]4=[CH]5[NH]16.[Fe]23456",
  };
  std::uint64_t total = 0;
  for (auto i = 0u; i < 100; ++i) {
    for (const auto &smiles : smileses) {
      auto molPtr = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(molPtr);
      total += molPtr->getNumAtoms();
    }
  }
  REQUIRE(total > 100);
}

// Explicit main needed as Catch returns odd exit codes by default and messes
// with gcovr.
int main(int argc, char *argv[]) { return Catch::Session().run(argc, argv); }
