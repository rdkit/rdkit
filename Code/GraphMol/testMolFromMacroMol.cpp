//
// Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/Atom.h>
#include <GraphMol/MacroMol.h>
#include <GraphMol/MacroMolTemplate.h>
#include <GraphMol/MacroMolTemplateLibrary.h>
#include <GraphMol/MolFromMacroMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <catch2/catch_all.hpp>

#include <memory>
#include <string>
#include <vector>

using namespace RDKit;

namespace {

constexpr bool doIsomericSmiles = true;
constexpr bool doKekule = false;
constexpr int rootedAtAtom = 0;
constexpr bool canonical = false;

std::unique_ptr<MacroMolTemplate> makeAlanineTemplate() {
  v2::SmilesParse::SmilesParserParams params;
  params.removeHs = false;
  const std::string smiles = "C[C@H](N[H])C(=O)O";
  auto alanine = v2::SmilesParse::MolFromSmiles(smiles, params);
  MacroMolTemplateBuilder builder(*alanine, MonomerClass::AminoAcid, "ALA",
                                  "A", smiles);
  builder.setMainGroup({0, 1, 2, 4, 5})
      .addLeavingGroup({{3}, 2, 3, 1})
      .addLeavingGroup({{6}, 4, 6, 2});
  return builder.build();
}

std::unique_ptr<MacroMolTemplate> makeGlycineTemplate() {
  v2::SmilesParse::SmilesParserParams params;
  params.removeHs = false;
  const std::string smiles = "O=C(CN[H])O";
  auto glycine = v2::SmilesParse::MolFromSmiles(smiles, params);
  MacroMolTemplateBuilder builder(*glycine, MonomerClass::AminoAcid, "GLY",
                                  "G", smiles);
  builder.setMainGroup({0, 1, 2, 3})
      .addLeavingGroup({{4}, 3, 4, 1})
      .addLeavingGroup({{5}, 1, 5, 2});
  return builder.build();
}

const SubstanceGroup &getOnlyMainSgroup(const ROMol &mol) {
  const auto &sgroups = getSubstanceGroups(mol);
  REQUIRE(sgroups.size() == 1);
  return sgroups.front();
}

}  // namespace

TEST_CASE("MolFromMacroMol converts a single amino acid macro atom",
          "[MolFromMacroMol]") {
  MacroMolTemplateLibrary templates;
  templates.addTemplate(makeAlanineTemplate());

  MacroMol macroMol;
  macroMol.addMacroAtom("A", MonomerClass::AminoAcid);

  auto mol = MolFromMacroMol(macroMol, templates);

  REQUIRE(mol);
  CHECK(mol->getNumAtoms() == 7);
  CHECK(mol->getNumBonds() == 6);
  CHECK(MolToSmiles(*mol, doIsomericSmiles, doKekule, rootedAtAtom,
                    canonical) ==
        "C[C@H](N[H])C(=O)O");
  const auto &mainSgroup = getOnlyMainSgroup(*mol);
  CHECK(mainSgroup.getProp<std::string>("CLASS") == "AminoAcid");
  CHECK(mainSgroup.getProp<std::string>("LABEL") == "ALA");
  CHECK(mainSgroup.getAtoms() ==
        std::vector<unsigned int>({0, 1, 2, 4, 5}));
  REQUIRE(mainSgroup.getAttachPoints().size() == 2);
  CHECK(mainSgroup.getAttachPoints()[0].lvIdx == 3);
  CHECK(mainSgroup.getAttachPoints()[1].lvIdx == 6);
}

TEST_CASE("MolFromMacroMol converts connected amino acid macro atoms",
          "[MolFromMacroMol]") {
  MacroMolTemplateLibrary templates;
  templates.addTemplate(makeAlanineTemplate());
  templates.addTemplate(makeGlycineTemplate());

  MacroMol macroMol;
  const auto alanineMacroAtom = macroMol.addMacroAtom("A", MonomerClass::AminoAcid);
  const auto glycineMacroAtom = macroMol.addMacroAtom("G", MonomerClass::AminoAcid);
  macroMol.addMacroBond(alanineMacroAtom, glycineMacroAtom, 2, 1);

  auto mol = MolFromMacroMol(macroMol, templates);

  REQUIRE(mol);
  CHECK(mol->getNumAtoms() == 11);
  CHECK(mol->getNumBonds() == 10);
  CHECK(MolToSmiles(*mol, doIsomericSmiles, doKekule, rootedAtAtom,
                    canonical) ==
        "C[C@H](N[H])C(=O)NCC(=O)O");
  const auto &sgroups = getSubstanceGroups(*mol);
  REQUIRE(sgroups.size() == 2);
  CHECK(sgroups[0].getProp<std::string>("LABEL") == "ALA");
  CHECK(sgroups[1].getProp<std::string>("LABEL") == "GLY");
  REQUIRE(sgroups[0].getAttachPoints().size() == 2);
  REQUIRE(sgroups[1].getAttachPoints().size() == 2);
  CHECK(sgroups[0].getAttachPoints()[0].lvIdx >= 0);
  CHECK(sgroups[0].getAttachPoints()[1].lvIdx == -1);
  CHECK(sgroups[1].getAttachPoints()[0].lvIdx == -1);
  CHECK(sgroups[1].getAttachPoints()[1].lvIdx >= 0);
}

TEST_CASE("MolFromMacroMol converts mixed macro and atomistic atoms",
          "[MolFromMacroMol]") {
  MacroMolTemplateLibrary templates;
  templates.addTemplate(makeAlanineTemplate());

  MacroMol macroMol;
  const auto macroAtom = macroMol.addMacroAtom("A", MonomerClass::AminoAcid);
  const auto atomisticAtom = macroMol.addAtom(new Atom(6), false, true);
  macroMol.addMacroAtomToAtomBond(macroAtom, atomisticAtom, 1);

  auto mol = MolFromMacroMol(macroMol, templates);

  REQUIRE(mol);
  CHECK(mol->getNumAtoms() == 7);
  CHECK(mol->getNumBonds() == 6);
  CHECK(MolToSmiles(*mol, doIsomericSmiles, doKekule, rootedAtAtom,
                    canonical) ==
        "C[C@H](NC)C(=O)O");
  const auto &mainSgroup = getOnlyMainSgroup(*mol);
  CHECK(mainSgroup.getProp<std::string>("LABEL") == "ALA");
  CHECK(mainSgroup.getAttachPoints()[0].lvIdx == -1);
  CHECK(mainSgroup.getAttachPoints()[1].lvIdx >= 0);
}

TEST_CASE("MolFromMacroMol removes complete multi-atom leaving groups",
          "[MolFromMacroMol]") {
  const std::string smiles = "COC";
  auto parsed = v2::SmilesParse::MolFromSmiles(smiles);
  MacroMolTemplateBuilder builder(*parsed, MonomerClass::Chemical, "METHYL",
                                  "Me", smiles);
  builder.setMainGroup({0}).addLeavingGroup({{1, 2}, 0, 1, 1});

  MacroMolTemplateLibrary templates;
  templates.addTemplate(builder.build());

  SECTION("unused group remains") {
    MacroMol macroMol;
    macroMol.addMacroAtom("Me", MonomerClass::Chemical);
    auto mol = MolFromMacroMol(macroMol, templates);
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 3);
    CHECK(mol->getNumBonds() == 2);
    const auto &mainSgroup = getOnlyMainSgroup(*mol);
    CHECK(mainSgroup.getAttachPoints()[0].lvIdx == 1);
  }

  SECTION("used group is removed in full") {
    MacroMol macroMol;
    const auto macroAtom =
        macroMol.addMacroAtom("Me", MonomerClass::Chemical);
    const auto regularAtom = macroMol.addAtom(new Atom(6), false, true);
    macroMol.addMacroAtomToAtomBond(macroAtom, regularAtom, 1);

    auto mol = MolFromMacroMol(macroMol, templates);
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 2);
    CHECK(mol->getNumBonds() == 1);
    const auto &mainSgroup = getOnlyMainSgroup(*mol);
    CHECK(mainSgroup.getProp<std::string>("CLASS") == "Chemical");
    CHECK(mainSgroup.getProp<std::string>("LABEL") == "METHYL");
    CHECK(mainSgroup.getAttachPoints()[0].lvIdx == -1);
  }
}

TEST_CASE("MolFromMacroMol preserves chirality at an attachment atom",
          "[MolFromMacroMol]") {
  // The leaving atom (Cl) is bonded directly to the stereocenter, so the
  // stereocenter's bond to it is replaced by the inter-monomer bond.
  const std::string smiles = "[Cl][C@](F)(Br)I";
  auto parsed = v2::SmilesParse::MolFromSmiles(smiles);
  MacroMolTemplateBuilder builder(*parsed, MonomerClass::Chemical, "CFBRI",
                                  "X", smiles);
  builder.setMainGroup({1, 2, 3, 4}).addLeavingGroup({{0}, 1, 0, 1});

  MacroMolTemplateLibrary templates;
  templates.addTemplate(builder.build());

  MacroMol macroMol;
  const auto macroAtom = macroMol.addMacroAtom("X", MonomerClass::Chemical);
  const auto regularAtom = macroMol.addAtom(new Atom(6), false, true);
  macroMol.addMacroAtomToAtomBond(macroAtom, regularAtom, 1);

  auto mol = MolFromMacroMol(macroMol, templates);
  REQUIRE(mol);
  MolOps::sanitizeMol(*mol);

  // The new carbon takes the place of the Cl, so the configuration should
  // match the template with Cl replaced by C.
  auto expected = v2::SmilesParse::MolFromSmiles("C[C@](F)(Br)I");
  CHECK(MolToSmiles(*mol) == MolToSmiles(*expected));
}

TEST_CASE("MolFromMacroMol preserves aromaticity of template bonds",
          "[MolFromMacroMol]") {
  const std::string smiles = "c1ccc(cc1)O";
  auto parsed = v2::SmilesParse::MolFromSmiles(smiles);
  MacroMolTemplateBuilder builder(*parsed, MonomerClass::Chemical, "PHENYL",
                                  "Ph", smiles);
  builder.setMainGroup({0, 1, 2, 3, 4, 5}).addLeavingGroup({{6}, 3, 6, 1});

  MacroMolTemplateLibrary templates;
  templates.addTemplate(builder.build());

  MacroMol macroMol;
  const auto macroAtom = macroMol.addMacroAtom("Ph", MonomerClass::Chemical);
  const auto regularAtom = macroMol.addAtom(new Atom(6), false, true);
  macroMol.addMacroAtomToAtomBond(macroAtom, regularAtom, 1);

  auto mol = MolFromMacroMol(macroMol, templates);
  REQUIRE(mol);
  for (const auto *bond : mol->bonds()) {
    if (bond->getBondType() == Bond::AROMATIC) {
      CHECK(bond->getIsAromatic());
    }
  }
  REQUIRE_NOTHROW(MolOps::sanitizeMol(*mol));
  CHECK(MolToSmiles(*mol) == "Cc1ccccc1");
}

TEST_CASE("MolFromMacroMol handles a stereocenter in a leaving group",
          "[MolFromMacroMol]") {
  // The leaving group {1, 2, 3} contains the stereocenter at atom 1.
  const std::string smiles = "C[C@H](F)Cl";
  auto parsed = v2::SmilesParse::MolFromSmiles(smiles);
  MacroMolTemplateBuilder builder(*parsed, MonomerClass::Chemical, "CHFCL",
                                  "Y", smiles);
  builder.setMainGroup({0}).addLeavingGroup({{1, 2, 3}, 0, 1, 1});

  MacroMolTemplateLibrary templates;
  templates.addTemplate(builder.build());

  SECTION("unused group keeps its configuration") {
    MacroMol macroMol;
    macroMol.addMacroAtom("Y", MonomerClass::Chemical);
    auto mol = MolFromMacroMol(macroMol, templates);
    REQUIRE(mol);
    MolOps::sanitizeMol(*mol);
    CHECK(MolToSmiles(*mol) == MolToSmiles(*parsed));
  }

  SECTION("used group is removed with its stereocenter") {
    MacroMol macroMol;
    const auto macroAtom = macroMol.addMacroAtom("Y", MonomerClass::Chemical);
    const auto regularAtom = macroMol.addAtom(new Atom(8), false, true);
    macroMol.addMacroAtomToAtomBond(macroAtom, regularAtom, 1);

    auto mol = MolFromMacroMol(macroMol, templates);
    REQUIRE(mol);
    MolOps::sanitizeMol(*mol);
    CHECK(MolToSmiles(*mol) == "CO");
  }
}

TEST_CASE("MolFromMacroMol preserves double bond stereo",
          "[MolFromMacroMol]") {
  const auto convert = [](const std::string &smiles,
                          std::vector<unsigned int> &&mainAtoms,
                          MacroMolLeavingGroup &&leavingGroup, bool bond) {
    auto parsed = v2::SmilesParse::MolFromSmiles(smiles);
    MacroMolTemplateBuilder builder(*parsed, MonomerClass::Chemical, "DB",
                                    "Z", smiles);
    builder.setMainGroup(std::move(mainAtoms))
        .addLeavingGroup(std::move(leavingGroup));
    MacroMolTemplateLibrary templates;
    templates.addTemplate(builder.build());

    MacroMol macroMol;
    const auto macroAtom = macroMol.addMacroAtom("Z", MonomerClass::Chemical);
    if (bond) {
      const auto regularAtom = macroMol.addAtom(new Atom(7), false, true);
      macroMol.addMacroAtomToAtomBond(macroAtom, regularAtom, 1);
    }
    auto mol = MolFromMacroMol(macroMol, templates);
    REQUIRE(mol);
    // the stereo must be stored on the bond, not just implied by bond
    // directions
    for (const auto *bond : mol->bonds()) {
      if (bond->getBondType() == Bond::DOUBLE &&
          bond->getBeginAtom()->getAtomicNum() == 6 &&
          bond->getEndAtom()->getAtomicNum() == 6) {
        CHECK(bond->getStereo() > Bond::STEREOANY);
        CHECK(bond->getStereoAtoms().size() == 2);
      }
    }
    MolOps::sanitizeMol(*mol);
    return MolToSmiles(*mol);
  };

  SECTION("double bond away from the leaving group, unused") {
    CHECK(convert("C/C=C/C(=O)O", {0, 1, 2, 3, 4}, {{5}, 3, 5, 1}, false) ==
          "C/C=C/C(=O)O");
  }
  SECTION("double bond away from the leaving group, used") {
    CHECK(convert("C/C=C/C(=O)O", {0, 1, 2, 3, 4}, {{5}, 3, 5, 1}, true) ==
          "C/C=C/C(N)=O");
  }
  SECTION("leaving atom is a stereo atom, unused") {
    CHECK(convert("O/C=C/C", {1, 2, 3}, {{0}, 1, 0, 1}, false) == "C/C=C/O");
  }
  SECTION("leaving atom is a stereo atom, used") {
    CHECK(convert("O/C=C/C", {1, 2, 3}, {{0}, 1, 0, 1}, true) == "C/C=C/N");
  }
  SECTION("double bond between ordinary atoms") {
    MacroMol macroMol;
    for (unsigned int i = 0; i < 4; ++i) {
      macroMol.addAtom(new Atom(6), false, true);
    }
    macroMol.addBond(0u, 1u, Bond::SINGLE);
    macroMol.addBond(1u, 2u, Bond::DOUBLE);
    macroMol.addBond(2u, 3u, Bond::SINGLE);
    auto *doubleBond = macroMol.getBondWithIdx(1);
    doubleBond->setStereoAtoms(0, 3);
    doubleBond->setStereo(Bond::STEREOCIS);

    auto mol = MolFromMacroMol(macroMol, MacroMolTemplateLibrary());
    REQUIRE(mol);
    const auto *resultBond = mol->getBondBetweenAtoms(1, 2);
    REQUIRE(resultBond);
    CHECK(resultBond->getStereo() == Bond::STEREOCIS);
    CHECK(resultBond->getStereoAtoms() == std::vector<int>({0, 3}));
    MolOps::sanitizeMol(*mol);
    CHECK(MolToSmiles(*mol) == "C/C=C\\C");
  }
}

TEST_CASE("MolFromMacroMol works with the built-in templates",
          "[MolFromMacroMol]") {
  const auto &templates = getGlobalMacroMolTemplateLibrary();
  v2::SmilesParse::SmilesParserParams params;
  params.removeHs = false;

  MacroMol macroMol;
  const auto ala = macroMol.addMacroAtom("A", MonomerClass::AminoAcid);
  const auto gly = macroMol.addMacroAtom("G", MonomerClass::AminoAcid);
  macroMol.addMacroBond(ala, gly, 2, 1);

  SECTION("linear dipeptide") {
    auto mol = MolFromMacroMol(macroMol, templates);
    REQUIRE(mol);
    MolOps::sanitizeMol(*mol);
    // L-Ala-Gly; the unused leaving groups (N-terminal H and C-terminal OH)
    // remain
    auto expected = v2::SmilesParse::MolFromSmiles(
        "[H]N[C@@H](C)C(=O)NCC(=O)O", params);
    CHECK(MolToSmiles(*mol) == MolToSmiles(*expected));
    const auto &sgroups = getSubstanceGroups(*mol);
    REQUIRE(sgroups.size() == 2);
    CHECK(sgroups[0].getProp<std::string>("LABEL") == "ALA");
    CHECK(sgroups[1].getProp<std::string>("LABEL") == "GLY");
  }

  SECTION("cyclic dipeptide") {
    // a second macro bond between the same macro atoms
    macroMol.addMacroBond(gly, ala, 2, 1);
    auto mol = MolFromMacroMol(macroMol, templates);
    REQUIRE(mol);
    MolOps::sanitizeMol(*mol);
    auto expected = v2::SmilesParse::MolFromSmiles("N1[C@@H](C)C(=O)NCC1=O");
    CHECK(MolToSmiles(*mol) == MolToSmiles(*expected));
  }
}

TEST_CASE("MolFromMacroMol checks macro bond types", "[MolFromMacroMol]") {
  // Converts a template bonded at attachment point 1 to an ordinary carbon
  const auto convert = [](const std::string &smiles,
                          std::vector<unsigned int> &&mainAtoms,
                          MacroMolLeavingGroup &&leavingGroup,
                          Bond::BondType bondType) {
    auto parsed = v2::SmilesParse::MolFromSmiles(smiles);
    MacroMolTemplateBuilder builder(*parsed, MonomerClass::Chemical, "BT",
                                    "B", smiles);
    builder.setMainGroup(std::move(mainAtoms))
        .addLeavingGroup(std::move(leavingGroup));
    MacroMolTemplateLibrary templates;
    templates.addTemplate(builder.build());

    MacroMol macroMol;
    const auto macroAtom = macroMol.addMacroAtom("B", MonomerClass::Chemical);
    const auto regularAtom = macroMol.addAtom(new Atom(6), false, true);
    macroMol.addMacroAtomToAtomBond(macroAtom, regularAtom, 1, bondType);
    return MolFromMacroMol(macroMol, templates);
  };
  using Catch::Matchers::ContainsSubstring;

  SECTION("double-bonded leaving atom") {
    auto mol = convert("CC=O", {0, 1}, {{2}, 1, 2, 1}, Bond::DOUBLE);
    REQUIRE(mol);
    MolOps::sanitizeMol(*mol);
    CHECK(MolToSmiles(*mol) == "C=CC");

    CHECK_THROWS_WITH(
        convert("CC=O", {0, 1}, {{2}, 1, 2, 1}, Bond::SINGLE),
        ContainsSubstring("must have the same bond type as the bond to the "
                          "leaving atom"));
  }

  SECTION("single-bonded leaving atom") {
    CHECK_THROWS_WITH(
        convert("CCO", {0, 1}, {{2}, 1, 2, 1}, Bond::DOUBLE),
        ContainsSubstring("must have the same bond type as the bond to the "
                          "leaving atom"));
  }

  SECTION("attachment point without a leaving atom") {
    for (const auto bondType : {Bond::HYDROGEN, Bond::ZERO}) {
      auto mol = convert("CO", {0, 1}, {{}, 1, -1, 1}, bondType);
      REQUIRE(mol);
      CHECK(mol->getNumAtoms() == 3);
      CHECK(mol->getBondBetweenAtoms(1, 2)->getBondType() == bondType);
      // the hydrogen on the oxygen is not displaced
      CHECK(mol->getAtomWithIdx(1)->getTotalNumHs() == 1);
    }

    CHECK_THROWS_WITH(
        convert("CO", {0, 1}, {{}, 1, -1, 1}, Bond::SINGLE),
        ContainsSubstring("must be a hydrogen or zero-order bond"));
  }
}
