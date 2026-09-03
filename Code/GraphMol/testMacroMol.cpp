//
// Copyright (C) 2026 Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/RDLog.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/MacroMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <catch2/catch_all.hpp>

#include <memory>
#include <string>
#include <type_traits>
#include <utility>

using namespace RDKit;

namespace {

std::unique_ptr<MacroMolTemplate> makeMacroMolTemplate(
    const std::string &name, const std::string &symbol,
    MonomerClass monomerClass = MonomerClass::AminoAcid) {
  RWMol mol;
  mol.addAtom(new Atom(6), false, true);
  MacroMolTemplateBuilder builder(mol, monomerClass, name, symbol, "");
  builder.setMainGroup({0});
  return builder.build();
}

}  // namespace

static_assert(std::is_same_v<
              decltype(std::declval<MacroMol &>().getLocalTemplateLibrary()),
              MacroMolTemplateLibrary &>);
static_assert(std::is_same_v<
              decltype(std::declval<const MacroMol &>()
                           .getLocalTemplateLibrary()),
              const MacroMolTemplateLibrary &>);
static_assert(std::is_nothrow_move_constructible_v<MacroMol>);
static_assert(std::is_nothrow_move_assignable_v<MacroMol>);

TEST_CASE("testBuildMacroMol") {
  // Build a simple MacroMol with three macro atoms and two bonds, and check
  // that the MacroMol has the expected number of atoms and bonds, and that the
  // "sequence" of template names in the macro atoms is as expected.
  auto macro_mol = std::make_unique<MacroMol>();
  auto macro_atom_1 = macro_mol->addMacroAtom("A", MonomerClass::AminoAcid);
  auto macro_atom_2 = macro_mol->addMacroAtom("C", MonomerClass::AminoAcid);
  auto macro_atom_3 = macro_mol->addMacroAtom("D", MonomerClass::AminoAcid);
  auto num_bonds_1 = macro_mol->addMacroBond(macro_atom_1, macro_atom_2, 2, 1,
                                              Bond::BondType::SINGLE);
  auto num_bonds_2 = macro_mol->addMacroBond(macro_atom_2, macro_atom_3, 2, 1,
                                              Bond::BondType::SINGLE);
  auto bond_idx_1 = num_bonds_1 - 1;
  auto bond_idx_2 = num_bonds_2 - 1;
  auto bond_1 = macro_mol->getBondWithIdx(bond_idx_1);
  auto bond_2 = macro_mol->getBondWithIdx(bond_idx_2);
  REQUIRE(bond_1);
  REQUIRE(bond_2);
  CHECK(macro_mol->getNumAtoms() == 3);
  CHECK(macro_mol->getNumBonds() == 2);
  CHECK(num_bonds_1 == 1);
  CHECK(num_bonds_2 == 2);
  CHECK(bond_1->getIdx() == bond_idx_1);
  CHECK(bond_1->getBeginAtomIdx() == macro_atom_1);
  CHECK(bond_1->getEndAtomIdx() == macro_atom_2);
  CHECK(bond_2->getIdx() == bond_idx_2);
  CHECK(bond_2->getBeginAtomIdx() == macro_atom_2);
  CHECK(bond_2->getEndAtomIdx() == macro_atom_3);

  std::string sequence;
  for (const auto &atom : macro_mol->atoms()) {
    const auto *info = atom->getMacroAtomInfo();
    REQUIRE(info);
    sequence += info->getSymbol();
    CHECK(info->getMonomerClass() == MonomerClass::AminoAcid);
  }
  CHECK(sequence == "ACD");
  for (auto bond : {bond_1, bond_2}) {
    CHECK(bond->getBondType() == Bond::UNSPECIFIED);
    const auto *info = bond->getMacroBondInfo();
    REQUIRE(info);
    CHECK(info->getNumBonds() == 1);
    CHECK(info->getBond(0).beginAttachPt == 2);
    CHECK(info->getBond(0).endAttachPt == 1);
    CHECK(info->getBond(0).bondType ==
          static_cast<unsigned int>(Bond::BondType::SINGLE));
  }
}

TEST_CASE("testMultipleConnectionsSameMacroAtoms") {
  // Build a MacroMol with two macro atoms and two bonds between them. This
  // results in one graph bond with two macro bond records.
  auto macro_mol = std::make_unique<MacroMol>();
  auto macro_atom_1 = macro_mol->addMacroAtom("C", MonomerClass::AminoAcid);
  auto macro_atom_2 = macro_mol->addMacroAtom("C", MonomerClass::AminoAcid);
  auto num_bonds_1 = macro_mol->addMacroBond(macro_atom_1, macro_atom_2, 2, 1);
  auto num_bonds_2 = macro_mol->addMacroBond(macro_atom_1, macro_atom_2, 3, 3,
                                            Bond::BondType::DOUBLE);
  CHECK(num_bonds_1 == num_bonds_2);
  auto bond_idx = num_bonds_1 - 1;
  CHECK(macro_mol->getNumBonds() == 1);
  auto graph_bond = macro_mol->getBondWithIdx(bond_idx);
  CHECK(graph_bond->getBeginAtomIdx() == macro_atom_1);
  CHECK(graph_bond->getEndAtomIdx() == macro_atom_2);
  CHECK(graph_bond->getBondType() == Bond::BondType::UNSPECIFIED);
  const auto *bond_info = graph_bond->getMacroBondInfo();
  REQUIRE(bond_info);
  REQUIRE(bond_info->getNumBonds() == 2);
  CHECK(bond_info->getBond(0).beginAttachPt == 2);
  CHECK(bond_info->getBond(0).endAttachPt == 1);
  CHECK(bond_info->getBond(0).bondType ==
        static_cast<unsigned int>(Bond::BondType::SINGLE));
  CHECK(bond_info->getBond(1).beginAttachPt == 3);
  CHECK(bond_info->getBond(1).endAttachPt == 3);
  CHECK(bond_info->getBond(1).bondType ==
        static_cast<unsigned int>(Bond::BondType::DOUBLE));
}

TEST_CASE("testAddAtomToMacroAtomBond") {
  // Build a MacroMol with one macro atom and one regular atom, and add a bond
  // between them. Check that the MacroMol has the expected number of atoms and
  // bonds, and that the properties of the macro atom and the bond are as
  // expected.
  auto atomMacroAtom = std::make_unique<MacroMol>();
  auto atom = new Atom(6);
  auto atom_idx = atomMacroAtom->addAtom(atom, false, true);
  auto macro_atom_idx =
      atomMacroAtom->addMacroAtom("C", MonomerClass::AminoAcid);
  auto num_bonds =
      atomMacroAtom->addAtomToMacroAtomBond(atom_idx, macro_atom_idx, 1);
  auto bond_idx = num_bonds - 1;
  CHECK(atomMacroAtom->getNumAtoms() == 2);
  CHECK(atomMacroAtom->getNumBonds() == 1);
  CHECK(num_bonds == atomMacroAtom->getNumBonds());
  auto macro_atom = atomMacroAtom->getAtomWithIdx(macro_atom_idx);
  const auto *macro_info = macro_atom->getMacroAtomInfo();
  REQUIRE(macro_info);
  CHECK(macro_info->getSymbol() == "C");
  CHECK(macro_info->getMonomerClass() == MonomerClass::AminoAcid);
  auto bond = atomMacroAtom->getBondWithIdx(bond_idx);
  CHECK(bond->getBeginAtomIdx() == atom_idx);
  CHECK(bond->getEndAtomIdx() == macro_atom_idx);
  CHECK(bond->getBondType() == Bond::BondType::UNSPECIFIED);
  const auto *bond_info = bond->getMacroBondInfo();
  REQUIRE(bond_info);
  REQUIRE(bond_info->getNumBonds() == 1);
  CHECK(bond_info->getBond(0).beginAttachPt == -1);
  CHECK(bond_info->getBond(0).endAttachPt == 1);
  CHECK(bond_info->getBond(0).bondType ==
        static_cast<unsigned int>(Bond::BondType::SINGLE));
}

TEST_CASE("testAddMacroAtomToAtomBond") {
  // Build a MacroMol with one macro atom and one regular atom, and add a bond
  // between them. Check that the MacroMol has the expected number of atoms and
  // bonds, and that the properties of the macro atom and the bond are as
  // expected.
  auto macroAtomAtom = std::make_unique<MacroMol>();
  auto macro_atom_idx =
      macroAtomAtom->addMacroAtom("C", MonomerClass::AminoAcid);
  auto atom = new Atom(6);
  auto atom_idx = macroAtomAtom->addAtom(atom, false, true);
  auto num_bonds =
      macroAtomAtom->addMacroAtomToAtomBond(macro_atom_idx, atom_idx, 1);
  auto bond_idx = num_bonds - 1;
  CHECK(macroAtomAtom->getNumAtoms() == 2);
  CHECK(macroAtomAtom->getNumBonds() == 1);
  CHECK(num_bonds == macroAtomAtom->getNumBonds());
  auto macro_atom = macroAtomAtom->getAtomWithIdx(macro_atom_idx);
  const auto *macro_info = macro_atom->getMacroAtomInfo();
  REQUIRE(macro_info);
  CHECK(macro_info->getSymbol() == "C");
  CHECK(macro_info->getMonomerClass() == MonomerClass::AminoAcid);
  auto bond = macroAtomAtom->getBondWithIdx(bond_idx);
  CHECK(bond->getBeginAtomIdx() == macro_atom_idx);
  CHECK(bond->getEndAtomIdx() == atom_idx);
  CHECK(bond->getBondType() == Bond::BondType::UNSPECIFIED);
  const auto *bond_info = bond->getMacroBondInfo();
  REQUIRE(bond_info);
  REQUIRE(bond_info->getNumBonds() == 1);
  CHECK(bond_info->getBond(0).beginAttachPt == 1);
  CHECK(bond_info->getBond(0).endAttachPt == -1);
  CHECK(bond_info->getBond(0).bondType ==
        static_cast<unsigned int>(Bond::BondType::SINGLE));
}

TEST_CASE("testAddBond") {
  // Build a MacroMol with two macro atoms and attempt to add a regular bond
  // between them.
  auto macro_mol = std::make_unique<MacroMol>();
  auto macro_atom_1 = macro_mol->addMacroAtom("A", MonomerClass::AminoAcid);
  auto macro_atom_2 = macro_mol->addMacroAtom("C", MonomerClass::AminoAcid);
  CHECK_THROWS_AS(macro_mol->addBond(macro_atom_1, macro_atom_2),
                  Invar::Invariant);
  auto atom = new Atom(6);
  auto atom_idx = macro_mol->addAtom(atom, false, true);
  CHECK_THROWS_AS(macro_mol->addBond(atom_idx, macro_atom_1), Invar::Invariant);

  auto macro_atom = macro_mol->getAtomWithIdx(macro_atom_1);
  auto regular_atom = macro_mol->getAtomWithIdx(atom_idx);
  CHECK_THROWS_AS(macro_mol->addBond(macro_atom, regular_atom),
                  Invar::Invariant);
  CHECK_THROWS_AS(macro_mol->addBond(regular_atom, macro_atom),
                  Invar::Invariant);

  auto macro_bond = std::make_unique<Bond>(Bond::SINGLE);
  macro_bond->setBeginAtomIdx(macro_atom_1);
  macro_bond->setEndAtomIdx(macro_atom_2);
  CHECK_THROWS_AS(macro_mol->addBond(macro_bond.get()), Invar::Invariant);

  auto second_atom = new Atom(6);
  auto second_atom_idx = macro_mol->addAtom(second_atom, false, true);
  auto second_regular_atom = macro_mol->getAtomWithIdx(second_atom_idx);
  CHECK(macro_mol->addBond(regular_atom, second_regular_atom) == 1);
  CHECK(macro_mol->getNumBonds() == 1);
}

TEST_CASE("MacroMol owns a local template library") {
  SECTION("default construction") {
    MacroMol macroMol;
    CHECK(macroMol.getLocalTemplateLibrary().getByName(
              MonomerClass::AminoAcid, "ALA") == nullptr);

    auto alanine = makeMacroMolTemplate("ALA", "A");
    const auto *alaninePtr = alanine.get();
    macroMol.getLocalTemplateLibrary().addTemplate(std::move(alanine));

    CHECK(alanine == nullptr);
    CHECK(macroMol.getLocalTemplateLibrary().getBySymbol(
              MonomerClass::AminoAcid, "A") == alaninePtr);
    CHECK(macroMol.getLocalTemplateLibrary().getByName(
              MonomerClass::AminoAcid, "ALA") == alaninePtr);
  }

  SECTION("ownership transfer") {
    auto localTemplateLibrary =
        std::make_unique<MacroMolTemplateLibrary>();
    auto alanine = makeMacroMolTemplate("ALA", "A");
    const auto *alaninePtr = alanine.get();
    localTemplateLibrary->addTemplate(std::move(alanine));

    MacroMol macroMol(std::move(localTemplateLibrary));

    CHECK(localTemplateLibrary == nullptr);
    CHECK(macroMol.getLocalTemplateLibrary().getBySymbol(
              MonomerClass::AminoAcid, "A") == alaninePtr);
  }

  SECTION("replacement ownership transfer") {
    MacroMol macroMol;
    auto localTemplateLibrary =
        std::make_unique<MacroMolTemplateLibrary>();
    auto alanine = makeMacroMolTemplate("ALA", "A");
    const auto *alaninePtr = alanine.get();
    localTemplateLibrary->addTemplate(std::move(alanine));

    macroMol.setLocalTemplateLibrary(std::move(localTemplateLibrary));

    CHECK(localTemplateLibrary == nullptr);
    CHECK(macroMol.getLocalTemplateLibrary().getBySymbol(
              MonomerClass::AminoAcid, "A") == alaninePtr);
  }

  SECTION("null library") {
    CHECK_THROWS_AS(
        MacroMol(std::unique_ptr<MacroMolTemplateLibrary>{}),
        Invar::Invariant);
    MacroMol macroMol;
    CHECK_THROWS_AS(
        macroMol.setLocalTemplateLibrary(
            std::unique_ptr<MacroMolTemplateLibrary>{}),
        Invar::Invariant);
  }
}

TEST_CASE("MacroMol validates references against only its local library") {
  SECTION("empty and ordinary-atom-only molecules are valid") {
    MacroMol macroMol;
    CHECK(macroMol.checkLocalTemplateReferences());
    macroMol.addAtom(new Atom(6), false, true);
    CHECK(macroMol.checkLocalTemplateReferences());
  }

  SECTION("symbol lookup") {
    MacroMol macroMol;
    macroMol.getLocalTemplateLibrary().addTemplate(
        makeMacroMolTemplate("ALA", "A"));
    macroMol.addMacroAtom("A", MonomerClass::AminoAcid);
    CHECK(macroMol.checkLocalTemplateReferences());
  }

  SECTION("template-name lookup") {
    MacroMol macroMol;
    macroMol.getLocalTemplateLibrary().addTemplate(
        makeMacroMolTemplate("ALA", "A"));
    macroMol.addMacroAtom("ALA", MonomerClass::AminoAcid);
    CHECK(macroMol.checkLocalTemplateReferences());
  }

  SECTION("missing template") {
    MacroMol macroMol;
    macroMol.addMacroAtom("A", MonomerClass::AminoAcid);
    CHECK_FALSE(macroMol.checkLocalTemplateReferences());
  }

  SECTION("monomer class is part of the lookup key") {
    MacroMol macroMol;
    macroMol.getLocalTemplateLibrary().addTemplate(
        makeMacroMolTemplate("ALA", "A"));
    macroMol.addMacroAtom("A", MonomerClass::NucleicAcid);
    CHECK_FALSE(macroMol.checkLocalTemplateReferences());
  }

  SECTION("one unresolved macro atom invalidates the molecule") {
    MacroMol macroMol;
    macroMol.getLocalTemplateLibrary().addTemplate(
        makeMacroMolTemplate("ALA", "A"));
    macroMol.addMacroAtom("A", MonomerClass::AminoAcid);
    macroMol.addMacroAtom("C", MonomerClass::AminoAcid);
    CHECK_FALSE(macroMol.checkLocalTemplateReferences());
  }
}

TEST_CASE("MacroMol local-template library rejects null templates") {
  MacroMol macroMol;
  CHECK_THROWS_AS(macroMol.getLocalTemplateLibrary().addTemplate(nullptr),
                  Invar::Invariant);
}

TEST_CASE("MacroMol copies have independent local template libraries") {
  MacroMol source;
  source.getLocalTemplateLibrary().addTemplate(
      makeMacroMolTemplate("ALA", "A"));
  source.addMacroAtom("A", MonomerClass::AminoAcid);
  const auto *sourceTemplate = source.getLocalTemplateLibrary().getBySymbol(
      MonomerClass::AminoAcid, "A");
  REQUIRE(sourceTemplate);

  MacroMol copy(source);
  CHECK(&copy.getLocalTemplateLibrary() !=
        &source.getLocalTemplateLibrary());
  const auto *copiedTemplate = copy.getLocalTemplateLibrary().getBySymbol(
      MonomerClass::AminoAcid, "A");
  REQUIRE(copiedTemplate);
  CHECK(copiedTemplate != sourceTemplate);
  CHECK(copiedTemplate->getName() == sourceTemplate->getName());
  CHECK(copy.checkLocalTemplateReferences());

  MacroMol assigned;
  assigned.getLocalTemplateLibrary().addTemplate(
      makeMacroMolTemplate("CYS", "C"));
  assigned = source;
  CHECK(&assigned.getLocalTemplateLibrary() !=
        &source.getLocalTemplateLibrary());
  const auto *assignedTemplate = assigned.getLocalTemplateLibrary().getBySymbol(
      MonomerClass::AminoAcid, "A");
  REQUIRE(assignedTemplate);
  CHECK(assignedTemplate != sourceTemplate);
  CHECK(assignedTemplate != copiedTemplate);
  CHECK(assigned.getLocalTemplateLibrary().getBySymbol(
            MonomerClass::AminoAcid, "C") == nullptr);
  CHECK(assigned.checkLocalTemplateReferences());
}

TEST_CASE("MacroMol moves its local template library") {
  MacroMol source;
  source.getLocalTemplateLibrary().addTemplate(
      makeMacroMolTemplate("ALA", "A"));
  source.addMacroAtom("A", MonomerClass::AminoAcid);
  const auto *sourceLibrary = &source.getLocalTemplateLibrary();
  const auto *sourceTemplate = sourceLibrary->getBySymbol(
      MonomerClass::AminoAcid, "A");
  REQUIRE(sourceTemplate);

  MacroMol moved(std::move(source));
  CHECK(&moved.getLocalTemplateLibrary() == sourceLibrary);
  CHECK(moved.getLocalTemplateLibrary().getBySymbol(
            MonomerClass::AminoAcid, "A") == sourceTemplate);
  CHECK(moved.checkLocalTemplateReferences());

  MacroMol assigned;
  assigned.getLocalTemplateLibrary().addTemplate(
      makeMacroMolTemplate("CYS", "C"));
  assigned = std::move(moved);
  CHECK(&assigned.getLocalTemplateLibrary() == sourceLibrary);
  CHECK(assigned.getLocalTemplateLibrary().getBySymbol(
            MonomerClass::AminoAcid, "A") == sourceTemplate);
  CHECK(assigned.checkLocalTemplateReferences());
}
