//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>

#include <GraphMol/RWMol.h>
#include <GraphMol/MonomerMol/MonomerLibrary.h>
#include <GraphMol/MonomerMol/MonomerMol.h>
#include <GraphMol/MonomerMol/Conversions.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("MonomerLibrary") {
  SECTION("GlobalLibraryBasics") {
    // The global library is available by default
    CHECK(MonomerLibrary::isUsingGlobalLibrary() == true);

    // Get the global library
    auto& globalLib = MonomerLibrary::getGlobalLibrary();

    // Check that built-in monomers are available
    CHECK(globalLib.hasMonomer("A", "PEPTIDE"));
    CHECK(globalLib.hasMonomer("G", "PEPTIDE"));
    CHECK(globalLib.hasMonomer("C", "PEPTIDE"));

    // Check that non-existent monomers return false
    CHECK(globalLib.hasMonomer("XYZ", "PEPTIDE") == false);
    CHECK(globalLib.hasMonomer("A", "RNA") == false);

    // Get monomer data
    auto alanineData = globalLib.getMonomerData("A", "PEPTIDE");
    REQUIRE(alanineData.has_value());
    CHECK(alanineData->find("C[C@H]") != std::string::npos);  // Alanine SMILES

    // Get monomer info from PDB code
    auto alaInfo = globalLib.getMonomerInfo("ALA");
    REQUIRE(alaInfo.has_value());
    CHECK(std::get<0>(*alaInfo) == "A");  // symbol
    CHECK(std::get<2>(*alaInfo) == "PEPTIDE");  // class

    // Get PDB code from symbol
    auto pdbCode = globalLib.getPdbCode("A", "PEPTIDE");
    REQUIRE(pdbCode.has_value());
    CHECK(*pdbCode == "ALA");
  }

  SECTION("GlobalLibraryWithMonomerMol") {
    // When no custom library is set, MonomerMol uses the global library
    MonomerMol mol;
    CHECK(mol.hasCustomLibrary() == false);

    // getMonomerLibrary() returns the global library
    auto& lib = mol.getMonomerLibrary();
    CHECK(lib.hasMonomer("A", "PEPTIDE"));

    // Build a simple peptide using the global library
    mol.addMonomer("A", 1, "PEPTIDE", "PEPTIDE1");
    mol.addMonomer("G");
    mol.addMonomer("C");
    mol.addConnection(0, 1, BACKBONE_LINKAGE);
    mol.addConnection(1, 2, BACKBONE_LINKAGE);

    CHECK(mol.getNumAtoms() == 3);
    CHECK(mol.getNumBonds() == 2);

    // Convert to atomistic - uses global library
    auto atomistic = toAtomistic(mol);
    CHECK(atomistic->getNumAtoms() > 0);
  }

  SECTION("CustomLibraryViaConstructor") {
    // Create a custom library with built-ins plus a non-standard monomer
    auto customLib = std::make_shared<MonomerLibrary>(true);  // load built-ins

    // The custom library starts with built-in definitions
    CHECK(customLib->hasMonomer("A", "PEPTIDE"));

    // Add a custom monomer (fictitious amino acid "X")
    customLib->addMonomerFromSmiles(
        "CC(C)(C)[C@H](N[H:1])C(=O)[OH:2]",  // tert-butyl glycine
        "X",
        "PEPTIDE",
        "TBG"
    );
    CHECK(customLib->hasMonomer("X", "PEPTIDE"));

    // Create MonomerMol with the custom library
    MonomerMol mol(customLib);
    CHECK(mol.hasCustomLibrary() == true);

    // The custom monomer is accessible
    auto& lib = mol.getMonomerLibrary();
    CHECK(lib.hasMonomer("X", "PEPTIDE"));

    // Build a peptide with the custom monomer
    mol.addMonomer("A", 1, "PEPTIDE", "PEPTIDE1");
    mol.addMonomer("X");  // custom monomer
    mol.addMonomer("G");
    mol.addConnection(0, 1, BACKBONE_LINKAGE);
    mol.addConnection(1, 2, BACKBONE_LINKAGE);

    CHECK(mol.getNumAtoms() == 3);

    // Convert to atomistic - uses custom library
    auto atomistic = toAtomistic(mol);
    CHECK(atomistic->getNumAtoms() > 0);

    // Verify the custom monomer was expanded correctly
    std::string smiles = MolToSmiles(*atomistic);
    CHECK(smiles.find("C(C)(C)") != std::string::npos);  // tert-butyl group
  }

  SECTION("CustomLibraryViaSetMonomerLibrary") {
    // Create a MonomerMol first (uses global by default)
    MonomerMol mol;
    CHECK(mol.hasCustomLibrary() == false);

    // Create and set a custom library with built-ins
    auto customLib = std::make_shared<MonomerLibrary>(true);  // load built-ins
    customLib->addMonomerFromSmiles(
        "NCCC[C@H](N[H:1])C(=O)[OH:2]",  // ornithine-like
        "Z",
        "PEPTIDE",
        "ORN"
    );

    mol.setMonomerLibrary(customLib);
    CHECK(mol.hasCustomLibrary() == true);

    // Now the custom monomer is available
    CHECK(mol.getMonomerLibrary().hasMonomer("Z", "PEPTIDE"));

    // Build peptide with custom monomer
    mol.addMonomer("Z", 1, "PEPTIDE", "PEPTIDE1");
    mol.addMonomer("A");
    mol.addConnection(0, 1, BACKBONE_LINKAGE);

    auto atomistic = toAtomistic(mol);
    CHECK(atomistic->getNumAtoms() > 0);
  }

  SECTION("ClearCustomLibrary") {
    auto customLib = std::make_shared<MonomerLibrary>();
    customLib->addMonomerFromSmiles("CC", "Y", "PEPTIDE");

    MonomerMol mol(customLib);
    CHECK(mol.hasCustomLibrary() == true);

    // Clear by setting nullptr - reverts to global
    mol.setMonomerLibrary(nullptr);
    CHECK(mol.hasCustomLibrary() == false);

    // Now uses global library again
    auto& lib = mol.getMonomerLibrary();
    CHECK(&lib == &MonomerLibrary::getGlobalLibrary());
  }

  SECTION("SharedLibraryBetweenMolecules") {
    // Multiple MonomerMols can share the same custom library
    auto sharedLib = std::make_shared<MonomerLibrary>();
    sharedLib->addMonomerFromSmiles("CC", "Q1", "PEPTIDE", "QQ1");

    MonomerMol mol1(sharedLib);
    MonomerMol mol2(sharedLib);

    CHECK(mol1.hasCustomLibrary() == true);
    CHECK(mol2.hasCustomLibrary() == true);

    // Both reference the same library
    CHECK(&mol1.getMonomerLibrary() == &mol2.getMonomerLibrary());

    // Adding to the shared library affects both
    sharedLib->addMonomerFromSmiles("CCC", "Q2", "PEPTIDE", "QQ2");
    CHECK(mol1.getMonomerLibrary().hasMonomer("Q2", "PEPTIDE"));
    CHECK(mol2.getMonomerLibrary().hasMonomer("Q2", "PEPTIDE"));
  }

  SECTION("CopyPreservesLibrary") {
    auto customLib = std::make_shared<MonomerLibrary>();
    customLib->addMonomerFromSmiles("CC", "W1", "PEPTIDE");

    MonomerMol mol1(customLib);
    mol1.addMonomer("A", 1, "PEPTIDE", "PEPTIDE1");

    // Copy constructor preserves the library
    MonomerMol mol2(mol1);
    CHECK(mol2.hasCustomLibrary() == true);
    CHECK(&mol2.getMonomerLibrary() == &mol1.getMonomerLibrary());

    // Copy assignment also preserves
    MonomerMol mol3;
    mol3 = mol1;
    CHECK(mol3.hasCustomLibrary() == true);
  }

  SECTION("MovePreservesLibrary") {
    auto customLib = std::make_shared<MonomerLibrary>();
    customLib->addMonomerFromSmiles("CC", "M1", "PEPTIDE");

    MonomerMol mol1(customLib);
    mol1.addMonomer("A", 1, "PEPTIDE", "PEPTIDE1");

    // Move constructor transfers the library
    MonomerMol mol2(std::move(mol1));
    CHECK(mol2.hasCustomLibrary() == true);
    CHECK(mol2.getMonomerLibrary().hasMonomer("M1", "PEPTIDE"));
  }

  SECTION("GetMonomer") {
    auto& lib = MonomerLibrary::getGlobalLibrary();

    // Get a parsed molecule for alanine
    auto alaMol = lib.getMonomer("A", "PEPTIDE");
    REQUIRE(alaMol != nullptr);
    CHECK(alaMol->getNumAtoms() > 0);

    // Second call should return the same cached molecule
    auto alaMol2 = lib.getMonomer("A", "PEPTIDE");
    CHECK(alaMol.get() == alaMol2.get());

    // Non-existent monomer returns nullptr
    auto notFound = lib.getMonomer("NOTEXIST", "PEPTIDE");
    CHECK(notFound == nullptr);
  }

  SECTION("AddMonomerWithMol") {
    auto customLib = std::make_shared<MonomerLibrary>();

    // Pre-parse a molecule
    auto mol = std::shared_ptr<ROMol>(SmilesToMol("CC(N)C(=O)O"));

    // Add with pre-parsed mol (no original data needed)
    customLib->addMonomer(mol, "TEST", "PEPTIDE", "TST");

    // getMol returns the same pre-parsed molecule
    auto retrieved = customLib->getMonomer("TEST", "PEPTIDE");
    CHECK(retrieved.get() == mol.get());
  }

  SECTION("AddMonomerFromSDF") {
    auto customLib = std::make_shared<MonomerLibrary>();

    // Simple alanine-like structure in SDF format
    std::string sdfData = R"(
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2500    1.2990    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.2500   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  2  0
M  END
)";

    customLib->addMonomerFromSDF(sdfData, "SDF1", "PEPTIDE", "SD1");

    CHECK(customLib->hasMonomer("SDF1", "PEPTIDE"));

    // getMol should parse the SDF and return a molecule
    auto mol = customLib->getMonomer("SDF1", "PEPTIDE");
    REQUIRE(mol != nullptr);
    CHECK(mol->getNumAtoms() == 4);
    CHECK(mol->getNumBonds() == 3);
  }

  SECTION("GlobalLibraryConfiguration") {
    // Save original state
    bool originalState = MonomerLibrary::isUsingGlobalLibrary();

    // Can toggle global library mode
    MonomerLibrary::useGlobalLibrary(false);
    CHECK(MonomerLibrary::isUsingGlobalLibrary() == false);

    MonomerLibrary::useGlobalLibrary(true);
    CHECK(MonomerLibrary::isUsingGlobalLibrary() == true);

    // Restore original state
    MonomerLibrary::useGlobalLibrary(originalState);
  }

  SECTION("EmptyLibrary") {
    // Create an empty library (default behavior)
    auto emptyLib = std::make_shared<MonomerLibrary>();

    // Should not have any built-in monomers
    CHECK(emptyLib->hasMonomer("A", "PEPTIDE") == false);
    CHECK(emptyLib->hasMonomer("G", "PEPTIDE") == false);
    CHECK(emptyLib->hasMonomer("C", "PEPTIDE") == false);

    // Can still add custom monomers
    emptyLib->addMonomerFromSmiles(
        "CC[C@H](N[H:1])C(=O)[OH:2]",
        "CUSTOM",
        "PEPTIDE",
        "CUS"
    );
    CHECK(emptyLib->hasMonomer("CUSTOM", "PEPTIDE"));

    // Use with MonomerMol
    MonomerMol mol(emptyLib);
    mol.addMonomer("CUSTOM", 1, "PEPTIDE", "PEPTIDE1");
    CHECK(mol.getNumAtoms() == 1);

    // Convert to atomistic works with custom monomer
    auto atomistic = toAtomistic(mol);
    CHECK(atomistic->getNumAtoms() > 0);
  }

  SECTION("LibraryWithBuiltins") {
    // Explicitly create library with built-ins
    auto libWithBuiltins = std::make_shared<MonomerLibrary>(true);

    // Should have built-in monomers
    CHECK(libWithBuiltins->hasMonomer("A", "PEPTIDE"));
    CHECK(libWithBuiltins->hasMonomer("G", "PEPTIDE"));
  }
}
