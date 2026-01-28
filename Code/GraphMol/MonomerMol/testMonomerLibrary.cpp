//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "MonomerLibrary.h"
#include "MonomerMol.h"

using namespace RDKit;

TEST_CASE("MonomerLibrary basic operations", "[MonomerLibrary]") {
    SECTION("empty library") {
        MonomerLibrary lib;
        CHECK(lib.empty());
        CHECK(lib.size() == 0);
        CHECK_FALSE(lib.hasMonomer("A", "PEPTIDE"));
    }

    SECTION("add monomer from SMILES") {
        MonomerLibrary lib;

        // Add a simple monomer
        lib.addMonomerFromSmiles("X", "CUSTOM", "[*:1]CC[*:2]", "XXX");

        CHECK_FALSE(lib.empty());
        CHECK(lib.size() == 1);
        CHECK(lib.hasMonomer("X", "CUSTOM"));
        CHECK_FALSE(lib.hasMonomer("X", "PEPTIDE"));  // Different type

        // Check we can retrieve the SMILES
        auto smiles = lib.getMonomerSmiles("X", "CUSTOM");
        REQUIRE(smiles.has_value());
        CHECK_FALSE(smiles->empty());

        // Check we can retrieve the molecule
        auto mol = lib.getMonomerMol("X", "CUSTOM");
        REQUIRE(mol != nullptr);
        CHECK(mol->getNumAtoms() > 0);

        // Check PDB code lookup
        auto pdb = lib.getPdbCode("X", "CUSTOM");
        REQUIRE(pdb.has_value());
        CHECK(*pdb == "XXX");
    }

    SECTION("add monomer from RDKit molecule") {
        MonomerLibrary lib;

        // Create a molecule directly
        auto mol = std::shared_ptr<ROMol>(SmilesToMol("[*:1]CCC[*:2]"));
        REQUIRE(mol != nullptr);

        lib.addMonomerFromMol("Y", "CUSTOM", mol, "YYY");

        CHECK(lib.hasMonomer("Y", "CUSTOM"));

        // Retrieve and verify
        auto retrieved = lib.getMonomerMol("Y", "CUSTOM");
        REQUIRE(retrieved != nullptr);
        CHECK(retrieved->getNumAtoms() == mol->getNumAtoms());
    }

    SECTION("duplicate identical definition is allowed") {
        MonomerLibrary lib;

        lib.addMonomerFromSmiles("A", "TEST", "[*:1]C[*:2]");

        // Adding identical definition should not throw
        CHECK_NOTHROW(lib.addMonomerFromSmiles("A", "TEST", "[*:1]C[*:2]"));
        CHECK(lib.size() == 1);  // Still only one monomer
    }

    SECTION("duplicate conflicting definition throws") {
        MonomerLibrary lib;

        lib.addMonomerFromSmiles("A", "TEST", "[*:1]C[*:2]");

        // Adding conflicting definition should throw
        CHECK_THROWS_AS(
            lib.addMonomerFromSmiles("A", "TEST", "[*:1]CC[*:2]"),
            std::runtime_error);
    }

    SECTION("PDB code aliases") {
        MonomerLibrary lib;

        lib.addMonomerFromSmiles("H", "PEPTIDE", "[*:1]Cc1cnc[nH]1[*:2]", "HIS");

        // Add aliases for protonation states
        lib.addPdbAlias("HID", "H", "PEPTIDE");
        lib.addPdbAlias("HIE", "H", "PEPTIDE");

        // All aliases should resolve to the same monomer
        auto his = lib.getHelmInfo("HIS");
        auto hid = lib.getHelmInfo("HID");
        auto hie = lib.getHelmInfo("HIE");

        REQUIRE(his.has_value());
        REQUIRE(hid.has_value());
        REQUIRE(hie.has_value());

        CHECK(std::get<0>(*his) == "H");
        CHECK(std::get<0>(*hid) == "H");
        CHECK(std::get<0>(*hie) == "H");
    }

    SECTION("clear library") {
        MonomerLibrary lib;

        lib.addMonomerFromSmiles("A", "TEST", "[*:1]C[*:2]");
        lib.addMonomerFromSmiles("B", "TEST", "[*:1]CC[*:2]");
        CHECK(lib.size() == 2);

        lib.clear();
        CHECK(lib.empty());
        CHECK(lib.size() == 0);
    }
}

TEST_CASE("MonomerLibrary default amino acids", "[MonomerLibrary]") {
    auto lib = MonomerLibrary::createWithDefaults();

    SECTION("has standard amino acids") {
        // Check a few standard amino acids
        CHECK(lib->hasMonomer("A", "PEPTIDE"));  // Alanine
        CHECK(lib->hasMonomer("G", "PEPTIDE"));  // Glycine
        CHECK(lib->hasMonomer("W", "PEPTIDE"));  // Tryptophan
        CHECK(lib->hasMonomer("Y", "PEPTIDE"));  // Tyrosine

        // Check special amino acids
        CHECK(lib->hasMonomer("U", "PEPTIDE"));  // Selenocysteine
        CHECK(lib->hasMonomer("O", "PEPTIDE"));  // Pyrrolysine
    }

    SECTION("PDB code lookups work") {
        auto ala = lib->getHelmInfo("ALA");
        REQUIRE(ala.has_value());
        CHECK(std::get<0>(*ala) == "A");
        CHECK(std::get<2>(*ala) == "PEPTIDE");

        // Protonation variants
        auto hid = lib->getHelmInfo("HID");
        REQUIRE(hid.has_value());
        CHECK(std::get<0>(*hid) == "H");
    }

    SECTION("can retrieve molecules") {
        auto alaMol = lib->getMonomerMol("A", "PEPTIDE");
        REQUIRE(alaMol != nullptr);
        CHECK(alaMol->getNumAtoms() > 0);

        auto glyMol = lib->getMonomerMol("G", "PEPTIDE");
        REQUIRE(glyMol != nullptr);
        // Glycine is smaller than alanine
        CHECK(glyMol->getNumAtoms() < alaMol->getNumAtoms());
    }
}

TEST_CASE("Global MonomerLibrary", "[MonomerLibrary]") {
    // Reset to defaults first
    MonomerLibrary::setGlobalLibrary(nullptr);

    SECTION("global library is lazily initialized") {
        auto lib1 = MonomerLibrary::getGlobalLibrary();
        REQUIRE(lib1 != nullptr);
        CHECK_FALSE(lib1->empty());

        // Same instance returned
        auto lib2 = MonomerLibrary::getGlobalLibrary();
        CHECK(lib1.get() == lib2.get());
    }

    SECTION("can set custom global library") {
        auto customLib = std::make_shared<MonomerLibrary>();
        customLib->addMonomerFromSmiles("CUSTOM", "TEST", "[*:1]CCCC[*:2]");

        MonomerLibrary::setGlobalLibrary(customLib);

        auto retrieved = MonomerLibrary::getGlobalLibrary();
        CHECK(retrieved.get() == customLib.get());
        CHECK(retrieved->hasMonomer("CUSTOM", "TEST"));

        // Standard amino acids should NOT be present
        CHECK_FALSE(retrieved->hasMonomer("A", "PEPTIDE"));

        // Reset for other tests
        MonomerLibrary::setGlobalLibrary(nullptr);
    }

    SECTION("nullptr resets to defaults") {
        auto customLib = std::make_shared<MonomerLibrary>();
        MonomerLibrary::setGlobalLibrary(customLib);

        MonomerLibrary::setGlobalLibrary(nullptr);

        auto lib = MonomerLibrary::getGlobalLibrary();
        CHECK(lib.get() != customLib.get());
        CHECK(lib->hasMonomer("A", "PEPTIDE"));  // Defaults restored
    }
}

TEST_CASE("Shared library across multiple MonomerMols", "[MonomerLibrary]") {
    SECTION("multiple molecules share the same library") {
        // Create a shared library with custom definitions
        auto sharedLib = std::make_shared<MonomerLibrary>();
        sharedLib->addMonomerFromSmiles("X", "CUSTOM", "[*:1]C(=O)N[*:2]", "XAA");
        sharedLib->addMonomerFromSmiles("Y", "CUSTOM", "[*:1]C(=O)O[*:2]", "YAA");

        // Create multiple MonomerMols using the same library
        MonomerMol mol1;
        mol1.setMonomerLibrary(sharedLib);

        MonomerMol mol2;
        mol2.setMonomerLibrary(sharedLib);

        MonomerMol mol3;
        mol3.setMonomerLibrary(sharedLib);

        // All should share the same library instance
        CHECK(mol1.getMonomerLibrary().get() == sharedLib.get());
        CHECK(mol2.getMonomerLibrary().get() == sharedLib.get());
        CHECK(mol3.getMonomerLibrary().get() == sharedLib.get());

        // Modifications to the shared library are visible to all
        sharedLib->addMonomerFromSmiles("Z", "CUSTOM", "[*:1]CC(=O)[*:2]", "ZAA");

        CHECK(mol1.getEffectiveLibrary()->hasMonomer("Z", "CUSTOM"));
        CHECK(mol2.getEffectiveLibrary()->hasMonomer("Z", "CUSTOM"));
        CHECK(mol3.getEffectiveLibrary()->hasMonomer("Z", "CUSTOM"));
    }

    SECTION("copied MonomerMol shares library with original") {
        auto lib = std::make_shared<MonomerLibrary>();
        lib->addMonomerFromSmiles("A", "TEST", "[*:1]C[*:2]");

        MonomerMol original;
        original.setMonomerLibrary(lib);

        // Copy the molecule
        MonomerMol copy(original);

        // Both should share the same library
        CHECK(copy.getMonomerLibrary().get() == lib.get());
        CHECK(copy.getMonomerLibrary().get() == original.getMonomerLibrary().get());
    }
}

TEST_CASE("Per-molecule library (non-shared)", "[MonomerLibrary]") {
    SECTION("each molecule has its own library") {
        // Create separate libraries for each molecule
        auto lib1 = std::make_shared<MonomerLibrary>();
        lib1->addMonomerFromSmiles("X", "TYPE1", "[*:1]CC[*:2]");

        auto lib2 = std::make_shared<MonomerLibrary>();
        lib2->addMonomerFromSmiles("X", "TYPE2", "[*:1]CCC[*:2]");  // Different!

        MonomerMol mol1;
        mol1.setMonomerLibrary(lib1);

        MonomerMol mol2;
        mol2.setMonomerLibrary(lib2);

        // Each molecule has its own definitions
        CHECK(mol1.getEffectiveLibrary()->hasMonomer("X", "TYPE1"));
        CHECK_FALSE(mol1.getEffectiveLibrary()->hasMonomer("X", "TYPE2"));

        CHECK(mol2.getEffectiveLibrary()->hasMonomer("X", "TYPE2"));
        CHECK_FALSE(mol2.getEffectiveLibrary()->hasMonomer("X", "TYPE1"));

        // Libraries are different instances
        CHECK(mol1.getMonomerLibrary().get() != mol2.getMonomerLibrary().get());
    }

    SECTION("molecule without library falls back to global") {
        // Reset global to defaults
        MonomerLibrary::setGlobalLibrary(nullptr);

        MonomerMol mol;
        // No library set

        CHECK(mol.getMonomerLibrary() == nullptr);

        // getEffectiveLibrary should return global
        auto effective = mol.getEffectiveLibrary();
        REQUIRE(effective != nullptr);
        CHECK(effective->hasMonomer("A", "PEPTIDE"));  // From global defaults
    }
}

TEST_CASE("MonomerLibrary with different monomer types", "[MonomerLibrary]") {
    MonomerLibrary lib;

    SECTION("same symbol different types") {
        // Add same symbol for different monomer types
        lib.addMonomerFromSmiles("A", "PEPTIDE", "[*:1]C[C@H](N)C(=O)[*:2]");
        lib.addMonomerFromSmiles("A", "RNA", "[*:1]C1OC(n2cnc3c(N)ncnc23)C(O)C1[*:2]");

        CHECK(lib.hasMonomer("A", "PEPTIDE"));
        CHECK(lib.hasMonomer("A", "RNA"));

        // They should be different molecules
        auto peptideA = lib.getMonomerMol("A", "PEPTIDE");
        auto rnaA = lib.getMonomerMol("A", "RNA");

        REQUIRE(peptideA != nullptr);
        REQUIRE(rnaA != nullptr);

        // RNA adenosine is larger than alanine
        CHECK(rnaA->getNumAtoms() > peptideA->getNumAtoms());
    }

    SECTION("custom monomer types") {
        // You can use any string as monomer type
        lib.addMonomerFromSmiles("PEG", "LINKER", "[*:1]OCCOCCOCCO[*:2]");
        lib.addMonomerFromSmiles("CLICK", "LINKER",
                                  "[*:1]c1cn(nn1)C[*:2]");

        CHECK(lib.hasMonomer("PEG", "LINKER"));
        CHECK(lib.hasMonomer("CLICK", "LINKER"));
    }
}

TEST_CASE("MonomerLibrary source data preservation", "[MonomerLibrary]") {
    SECTION("SMILES source is preserved") {
        MonomerLibrary lib;
        std::string originalSmiles = "[*:1]CC(=O)N[*:2]";
        lib.addMonomerFromSmiles("X", "TEST", originalSmiles);

        auto sourceData = lib.getMonomerSourceData("X", "TEST");
        REQUIRE(sourceData.has_value());
        CHECK(*sourceData == originalSmiles);

        auto def = lib.getMonomerDef("X", "TEST");
        REQUIRE(def != nullptr);
        CHECK(def->sourceFormat == "SMILES");
    }

    SECTION("mol source has empty source data") {
        MonomerLibrary lib;
        auto mol = std::shared_ptr<ROMol>(SmilesToMol("[*:1]CC[*:2]"));
        lib.addMonomerFromMol("X", "TEST", mol);

        auto sourceData = lib.getMonomerSourceData("X", "TEST");
        REQUIRE(sourceData.has_value());
        CHECK(sourceData->empty());

        auto def = lib.getMonomerDef("X", "TEST");
        REQUIRE(def != nullptr);
        CHECK(def->sourceFormat == "MOL");
    }
}

TEST_CASE("MonomerLibrary error handling", "[MonomerLibrary]") {
    SECTION("invalid SMILES throws") {
        MonomerLibrary lib;
        CHECK_THROWS_AS(
            lib.addMonomerFromSmiles("X", "TEST", "not a valid smiles!!!"),
            std::runtime_error);
    }

    SECTION("null molecule throws") {
        MonomerLibrary lib;
        CHECK_THROWS_AS(
            lib.addMonomerFromMol("X", "TEST", nullptr),
            std::runtime_error);
    }

    SECTION("PDB alias for non-existent monomer throws") {
        MonomerLibrary lib;
        CHECK_THROWS_AS(
            lib.addPdbAlias("XXX", "NONEXISTENT", "PEPTIDE"),
            std::runtime_error);
    }

    SECTION("getMonomerMol returns nullptr for non-existent") {
        MonomerLibrary lib;
        auto mol = lib.getMonomerMol("NONEXISTENT", "TYPE");
        CHECK(mol == nullptr);
    }

    SECTION("getMonomerSmiles returns nullopt for non-existent") {
        MonomerLibrary lib;
        auto smiles = lib.getMonomerSmiles("NONEXISTENT", "TYPE");
        CHECK_FALSE(smiles.has_value());
    }
}

TEST_CASE("MonomerLibrary use case: building peptides", "[MonomerLibrary]") {
    // This test demonstrates a realistic use case

    SECTION("build a simple peptide using library") {
        // Get the global library with standard amino acids
        auto lib = MonomerLibrary::getGlobalLibrary();

        // Create a MonomerMol and attach the library
        MonomerMol peptide;
        peptide.setMonomerLibrary(lib);

        // Add monomers
        auto ala = peptide.addMonomer("A", 1, "PEPTIDE1");
        auto gly = peptide.addMonomer("G", 2, "PEPTIDE1");
        auto val = peptide.addMonomer("V", 3, "PEPTIDE1");

        // Connect them
        peptide.addConnection(ala, gly, ConnectionType::FORWARD);
        peptide.addConnection(gly, val, ConnectionType::FORWARD);

        CHECK(peptide.getNumAtoms() == 3);
        CHECK(peptide.getNumBonds() == 2);

        // The library should be accessible
        CHECK(peptide.getEffectiveLibrary()->hasMonomer("A", "PEPTIDE"));
    }
}

TEST_CASE("MonomerLibrary use case: custom monomers for drug-peptide conjugate",
          "[MonomerLibrary]") {
    // Create a library with both standard amino acids and custom linkers
    auto lib = std::make_shared<MonomerLibrary>();

    // Add some standard amino acids manually (or copy from defaults)
    lib->addMonomerFromSmiles("A", "PEPTIDE",
                               "C[C@H](N[*:1])C(=O)[*:2]", "ALA");
    lib->addMonomerFromSmiles("K", "PEPTIDE",
                               "NCCCC[C@H](N[*:1])C(=O)[*:2]", "LYS");

    // Add a custom PEG linker
    lib->addMonomerFromSmiles("PEG4", "LINKER",
                               "[*:1]OCCOCCOCCOCC[*:2]");

    // Add a custom drug-like moiety
    lib->addMonomerFromSmiles("DRUG", "PAYLOAD",
                               "[*:1]c1ccc(F)c(N)c1");

    CHECK(lib->size() == 4);
    CHECK(lib->hasMonomer("A", "PEPTIDE"));
    CHECK(lib->hasMonomer("K", "PEPTIDE"));
    CHECK(lib->hasMonomer("PEG4", "LINKER"));
    CHECK(lib->hasMonomer("DRUG", "PAYLOAD"));

    // Create a MonomerMol using this custom library
    MonomerMol conjugate;
    conjugate.setMonomerLibrary(lib);

    // The effective library should have our custom monomers
    auto effective = conjugate.getEffectiveLibrary();
    CHECK(effective->hasMonomer("PEG4", "LINKER"));
    CHECK(effective->hasMonomer("DRUG", "PAYLOAD"));
}
