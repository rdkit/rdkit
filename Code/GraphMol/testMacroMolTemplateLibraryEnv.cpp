// Copyright (C) 2026 Schrödinger and other RDKit contributors
//
// Tests for the RDKIT_MONOMER_LIBRARY environment variable.
//
// These tests are kept separate from testMacroMolTemplateLibrary.cpp because
// the global MacroMol template library is loaded only once per process, on
// first use, and the result is kept for the lifetime of the process. Once the
// tests in that file have loaded the default templates, setting
// RDKIT_MONOMER_LIBRARY has no further effect, so these tests must run in a
// process of their own that sets the variable before the library is first
// used.

#include <GraphMol/MacroMolTemplateLibrary.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/utils.h>
#include <catch2/catch_all.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace RDKit;
using Catch::Matchers::ContainsSubstring;

TEST_CASE("RDKIT_MONOMER_LIBRARY replaces the built-in MacroMol templates") {
  const auto *rdbase = std::getenv("RDBASE");
  REQUIRE(rdbase);
  const auto customFile = std::filesystem::path(rdbase) / "Code" / "GraphMol" /
                          "test_data" / "custom_monomers.json";

  // A failed load is not cached, so the error cases must come before the
  // successful one.
  const auto missingFile =
      std::filesystem::temp_directory_path() / "rdkit_no_such_monomers.json";
  std::filesystem::remove(missingFile);
  setenv("RDKIT_MONOMER_LIBRARY", missingFile.string().c_str(), 1);
  CHECK_THROWS_WITH(getGlobalMacroMolTemplateLibrary(),
                    ContainsSubstring("Could not open monomer library") &&
                        ContainsSubstring(missingFile.string()));

  // Writes a temporary monomer library and checks that loading it fails
  const auto checkLoadError = [](const std::string &json,
                                 const std::string &expectedMessage) {
    const auto file =
        std::filesystem::temp_directory_path() / "rdkit_bad_monomers.json";
    {
      std::ofstream out(file);
      out << json;
    }
    setenv("RDKIT_MONOMER_LIBRARY", file.string().c_str(), 1);
    CHECK_THROWS_WITH(getGlobalMacroMolTemplateLibrary(),
                      ContainsSubstring("Error loading monomer library") &&
                          ContainsSubstring(file.string()) &&
                          ContainsSubstring(expectedMessage));
    std::filesystem::remove(file);
  };
  const std::string glycine = R"({"SYMBOL": "G", "POLYMER_TYPE": "PEPTIDE",
      "CORE_SMILES": "*NCC(*)=O |$_R1;;;;_R2;$|",
      "RGROUPS": {"R1": "*[H] |$_R1;$|", "R2": "*O |$_R2;$|"},
      "MONOMER_TYPE": "Backbone", "PDBCODE": "GLY"})";

  // malformed JSON
  checkLoadError(R"([{"SYMBOL": "Sar")", "");
  // an R-group missing from RGROUPS
  checkLoadError(R"([{"SYMBOL": "X", "POLYMER_TYPE": "PEPTIDE",
                     "CORE_SMILES": "*NCC(*)=O |$_R1;;;;_R2;$|",
                     "RGROUPS": {"R1": "*[H] |$_R1;$|"},
                     "MONOMER_TYPE": "Backbone", "PDBCODE": "XXX"}])",
                 "do not match its RGROUPS");
  // fields that are an object instead of a string
  for (const auto *symbol : {R"({"x": "y"})", "{}"}) {
    checkLoadError(R"([{"SYMBOL": )" + std::string(symbol) + R"(,
                       "POLYMER_TYPE": "PEPTIDE",
                       "CORE_SMILES": "*NCC(*)=O |$_R1;;;;_R2;$|",
                       "RGROUPS": {"R1": "*[H] |$_R1;$|", "R2": "*O |$_R2;$|"},
                       "MONOMER_TYPE": "Backbone", "PDBCODE": "XXX"}])",
                   "Missing or invalid 'SYMBOL'");
  }
  // duplicate monomers
  checkLoadError("[" + glycine + ", " + glycine + "]",
                 "already contains an entry");

  setenv("RDKIT_MONOMER_LIBRARY", customFile.string().c_str(), 1);
  const auto &library = getGlobalMacroMolTemplateLibrary();
  CHECK(library.getTemplates().size() == 5);
  CHECK(library.getBySymbol(MonomerClass::AminoAcid, "Sar"));
  CHECK(!library.getBySymbol(MonomerClass::AminoAcid, "A"));

  // a leaving group with more than one atom
  const auto *ester = library.getByName(MonomerClass::AminoAcid, "BAM");
  REQUIRE(ester);
  CHECK(ester->getMol().getNumAtoms() == 8);
  const auto &leavingGroups = ester->getLeavingGroups();
  REQUIRE(leavingGroups.size() == 2);
  const std::vector<unsigned int> expectedEsterAtoms{5, 7};
  CHECK(leavingGroups[1].atomIdxs == expectedEsterAtoms);
  CHECK(leavingGroups[1].attachAtomIdx == 4);
  CHECK(leavingGroups[1].leavingAtomIdx == 5);
  CHECK(leavingGroups[1].attachPoint == 2);

  // R-groups with stereo; the chirality must not depend on where the dummy
  // atom is written in the R-group SMILES
  for (const auto *name : {"HE1", "HE2"}) {
    const auto *templ = library.getByName(MonomerClass::AminoAcid, name);
    REQUIRE(templ);
    CHECK(MolToSmiles(templ->getMol()) == "[H]NCC(=O)[C@@H](C)O");
  }
  const auto *propenyl = library.getByName(MonomerClass::AminoAcid, "PRP");
  REQUIRE(propenyl);
  CHECK(MolToSmiles(propenyl->getMol()) == "[H]NCC(=O)/C=C/C");

  MacroMolTemplateLibrary local;
  addBuiltinMacroMolTemplates(local);
  CHECK(local.getTemplates().size() == 5);
}
