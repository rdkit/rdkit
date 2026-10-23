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
#include <GraphMol/FileParsers/MolFromMacroMol.h>
#include <GraphMol/MacroMol.h>
#include <GraphMol/MacroMolTemplate.h>
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
  SmilesParserParams params;
  params.removeHs = false;
  const std::string smiles = "C[C@H](N[H:1])C(=O)[OH:2]";
  auto alanine = std::unique_ptr<RWMol>(SmilesToMol(smiles, params));
  MacroMolTemplateBuilder builder(*alanine, MonomerClass::AminoAcid, "ALA",
                                  "A", smiles);
  builder.setMainGroup({0, 1, 2, 4, 5})
      .addLeavingGroup({{3}, 2, 3, 1})
      .addLeavingGroup({{6}, 4, 6, 2});
  return builder.build();
}

std::unique_ptr<MacroMolTemplate> makeGlycineTemplate() {
  SmilesParserParams params;
  params.removeHs = false;
  const std::string smiles = "O=C(CN[H:1])[OH:2]";
  auto glycine = std::unique_ptr<RWMol>(SmilesToMol(smiles, params));
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
        "C[C@H](N[H:1])C(=O)[OH:2]");
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
        "C[C@H](N[H:1])C(=O)NCC(=O)[OH:2]");
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
        "C[C@H](NC)C(=O)[OH:2]");
  const auto &mainSgroup = getOnlyMainSgroup(*mol);
  CHECK(mainSgroup.getProp<std::string>("LABEL") == "ALA");
  CHECK(mainSgroup.getAttachPoints()[0].lvIdx == -1);
  CHECK(mainSgroup.getAttachPoints()[1].lvIdx >= 0);
}

TEST_CASE("MolFromMacroMol removes complete multi-atom leaving groups",
          "[MolFromMacroMol]") {
  const std::string smiles = "COC";
  auto parsed = std::unique_ptr<RWMol>(SmilesToMol(smiles));
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
