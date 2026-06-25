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
  macro_mol->addMacroBond(macro_atom_1, macro_atom_2, 2, 1);
  macro_mol->addMacroBond(macro_atom_2, macro_atom_3, 2, 1);
  CHECK(macro_mol->getNumAtoms() == 3);
  CHECK(macro_mol->getNumBonds() == 2);
  std::string sequence;
  for (const auto &atom : macro_mol->atoms()) {
    std::string templateName =
        atom->getProp<std::string>(common_properties::dummyLabel);
    sequence += templateName;
    CHECK(atom->getProp<std::string>(common_properties::molAtomClass) == "AA");
  }
  CHECK(sequence == "ACD");
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

TEST_CASE("testSmilesCanonicalization") {
  // Build a MacroMol with two macro atoms with the same template name and
  // class name, and check that the smiles of the two macro atoms are the same
  auto macro_mol_1 = std::make_unique<MacroMol>();
  auto m1_atom_1 = macro_mol_1->addMacroAtom(MonomerClass::AA, "A");
  auto m1_atom_2 = macro_mol_1->addMacroAtom(MonomerClass::AA, "C");
  auto m1_atom_3 = macro_mol_1->addMacroAtom(MonomerClass::AA, "D");
  macro_mol_1->addMacroBond(m1_atom_1, m1_atom_2, 2, 1);
  macro_mol_1->addMacroBond(m1_atom_2, m1_atom_3, 2, 1);

  auto macro_mol_2 = std::make_unique<MacroMol>();
  auto m2_atom_1 = macro_mol_2->addMacroAtom(MonomerClass::AA, "C");
  auto m2_atom_2 = macro_mol_2->addMacroAtom(MonomerClass::AA, "D");
  auto m2_atom_3 = macro_mol_2->addMacroAtom(MonomerClass::AA, "A");
  macro_mol_2->addMacroBond(m2_atom_1, m2_atom_2, 2, 1);
  macro_mol_2->addMacroBond(m2_atom_3, m2_atom_1, 2, 1);

  std::string m1_smiles = MolToSmiles(*macro_mol_1);
  std::string m2_smiles = MolToSmiles(*macro_mol_2);
  CHECK(m1_smiles == m2_smiles);
}
