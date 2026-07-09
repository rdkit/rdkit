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

using namespace RDKit;

TEST_CASE("testBuildMacroMol") {
  // Build a simple MacroMol with three macro atoms and two bonds, and check
  // that the MacroMol has the expected number of atoms and bonds, and that the
  // "sequence" of template names in the macro atoms is as expected.
  auto macro_mol = std::make_unique<MacroMol>();
  auto macro_atom_1 = macro_mol->addMacroAtom(MonomerClass::AA, "A");
  auto macro_atom_2 = macro_mol->addMacroAtom(MonomerClass::AA, "C");
  auto macro_atom_3 = macro_mol->addMacroAtom(MonomerClass::AA, "D");
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
  CHECK(bond_1->getBondType() == Bond::SINGLE);
  CHECK(bond_2->getIdx() == bond_idx_2);
  CHECK(bond_2->getBeginAtomIdx() == macro_atom_2);
  CHECK(bond_2->getEndAtomIdx() == macro_atom_3);
  CHECK(bond_2->getBondType() == Bond::SINGLE);
  std::string sequence;
  for (const auto &atom : macro_mol->atoms()) {
    const auto *info = atom->getMacroAtomInfo();
    REQUIRE(info);
    sequence += info->getSymbol();
    CHECK(info->getMonomerClass() == MonomerClass::AA);
  }
  CHECK(sequence == "ACD");

  ROMol copied(*macro_mol);
  const auto *copiedInfo =
      copied.getAtomWithIdx(macro_atom_1)->getMacroAtomInfo();
  REQUIRE(copiedInfo);
  CHECK(copiedInfo->getSymbol() == "A");
  CHECK(copiedInfo->getMonomerClass() == MonomerClass::AA);
}

TEST_CASE("testMultipleConnectionsSameMacroAtoms") {
  // Build a MacroMol with two macro atoms and two bonds between them. This
  // results in an exception because the same macro atoms cannot be connected by
  // multiple bonds. Eventually, we will want to allow this, but for now we
  // raise an exception.
  auto macro_mol = std::make_unique<MacroMol>();
  auto macro_atom_1 = macro_mol->addMacroAtom(MonomerClass::AA, "C");
  auto macro_atom_2 = macro_mol->addMacroAtom(MonomerClass::AA, "C");
  macro_mol->addMacroBond(macro_atom_1, macro_atom_2, 2, 1);
  // Raises an exception because the same macro atoms cannot be connected by
  // multiple bonds
  CHECK_THROWS_AS(macro_mol->addMacroBond(macro_atom_1, macro_atom_2, 3, 3),
                  Invar::Invariant);
}

TEST_CASE("testAddAtomToMacroAtomBond") {
  // Build a MacroMol with one macro atom and one regular atom, and add a bond
  // between them. Check that the MacroMol has the expected number of atoms and
  // bonds, and that the properties of the macro atom and the bond are as
  // expected.
  auto atomMacroAtom = std::make_unique<MacroMol>();
  auto atom = new Atom(6);
  auto atom_idx = atomMacroAtom->addAtom(atom, false, true);
  auto macro_atom_idx = atomMacroAtom->addMacroAtom(MonomerClass::AA, "C");
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
  CHECK(macro_info->getMonomerClass() == MonomerClass::AA);
  auto bond = atomMacroAtom->getBondWithIdx(bond_idx);
  CHECK(bond->getBeginAtomIdx() == atom_idx);
  CHECK(bond->getEndAtomIdx() == macro_atom_idx);
}

TEST_CASE("testAddMacroAtomToAtomBond") {
  // Build a MacroMol with one macro atom and one regular atom, and add a bond
  // between them. Check that the MacroMol has the expected number of atoms and
  // bonds, and that the properties of the macro atom and the bond are as
  // expected.
  auto macroAtomAtom = std::make_unique<MacroMol>();
  auto macro_atom_idx = macroAtomAtom->addMacroAtom(MonomerClass::AA, "C");
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
  CHECK(macro_info->getMonomerClass() == MonomerClass::AA);
  auto bond = macroAtomAtom->getBondWithIdx(bond_idx);
  CHECK(bond->getBeginAtomIdx() == macro_atom_idx);
  CHECK(bond->getEndAtomIdx() == atom_idx);
}

TEST_CASE("testAddBond") {
  // Build a MacroMol with two macro atoms and attempt to add a regular bond
  // between them.
  auto macro_mol = std::make_unique<MacroMol>();
  auto macro_atom_1 = macro_mol->addMacroAtom(MonomerClass::AA, "A");
  auto macro_atom_2 = macro_mol->addMacroAtom(MonomerClass::AA, "C");
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
