//
// Copyright (C) 2026 Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/MacroMolTemplate.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <catch2/catch_all.hpp>

#include <memory>

using namespace RDKit;

TEST_CASE("testBuildMacroMolEntry") {
  // Build a simple MacroMolEntry and check that the public fields are
  // populated as expected.
  MacroMolEntry macroMolEntry;
  macroMolEntry.monomerClass = "AA";
  macroMolEntry.templateName = "ALA";
  macroMolEntry.symbol = "A";
  macroMolEntry.original_data = "N[C@@H](C)C(=O)O";
  macroMolEntry.molTemplate = std::make_shared<MacroMolTemplate>();

  CHECK(macroMolEntry.monomerClass == "AA");
  CHECK(macroMolEntry.templateName == "ALA");
  CHECK(macroMolEntry.symbol == "A");
  CHECK(macroMolEntry.original_data == "N[C@@H](C)C(=O)O");
  CHECK(macroMolEntry.molTemplate);
}

TEST_CASE("testMacroMolTemplateLibraryLookup") {
  // Build a MacroMolTemplateLibrary with two amino-acid entries and check that
  // lookup by template name and symbol returns the expected entry.
  MacroMolTemplateLibrary templateLibrary;
  auto alanineEntry = std::make_shared<MacroMolEntry>();
  alanineEntry->monomerClass = "AA";
  alanineEntry->templateName = "ALA";
  alanineEntry->symbol = "A";
  alanineEntry->molTemplate = std::make_shared<MacroMolTemplate>();

  auto cysteineEntry = std::make_shared<MacroMolEntry>();
  cysteineEntry->monomerClass = "AA";
  cysteineEntry->templateName = "CYS";
  cysteineEntry->symbol = "C";
  cysteineEntry->molTemplate = std::make_shared<MacroMolTemplate>();

  templateLibrary.addEntry(alanineEntry);
  templateLibrary.addEntry(cysteineEntry);

  CHECK(templateLibrary.getByTemplateName("AA", "ALA") == alanineEntry);
  CHECK(templateLibrary.getBySymbol("AA", "A") == alanineEntry);
  CHECK(templateLibrary.getByTemplateName("AA", "CYS") == cysteineEntry);
  CHECK(templateLibrary.getBySymbol("AA", "C") == cysteineEntry);
}

TEST_CASE("testMacroMolTemplateLibrarySeparatesMonomerClasses") {
  // Build a MacroMolTemplateLibrary with amino-acid and nucleic-acid entries
  // sharing the same template name and symbol. Check that the monomer class is
  // part of the lookup key.
  MacroMolTemplateLibrary templateLibrary;
  auto aaEntry = std::make_shared<MacroMolEntry>();
  aaEntry->monomerClass = "AA";
  aaEntry->templateName = "ALA";
  aaEntry->symbol = "A";
  aaEntry->molTemplate = std::make_shared<MacroMolTemplate>();

  auto naEntry = std::make_shared<MacroMolEntry>();
  naEntry->monomerClass = "NA";
  naEntry->templateName = "ADE";
  naEntry->symbol = "A";
  naEntry->molTemplate = std::make_shared<MacroMolTemplate>();

  templateLibrary.addEntry(aaEntry);
  templateLibrary.addEntry(naEntry);

  CHECK(templateLibrary.getByTemplateName("AA", "ALA") == aaEntry);
  CHECK(templateLibrary.getBySymbol("AA", "A") == aaEntry);
  CHECK(templateLibrary.getByTemplateName("NA", "ADE") == naEntry);
  CHECK(templateLibrary.getBySymbol("NA", "A") == naEntry);
}

TEST_CASE("testMacroMolTemplateMainAndLeavingGroups") {
  // Build an alanine template with explicit peptide-bond leaving groups:
  //   [H] on the N side and [OH] on the C-terminal O side. Atom indices:
  //   0:H(N)  1:N  2:C(alpha)  3:C(methyl)  4:C(carbonyl)  5:O(=O)  6:O(H)
  SmilesParserParams params;
  params.removeHs = false;
  std::unique_ptr<RWMol> parsed(SmilesToMol("[H]N[C@@H](C)C(=O)O", params));
  REQUIRE(parsed);
  MacroMolTemplate macroMolTemplate(*parsed);
  const auto &constMacroMolTemplate = macroMolTemplate;

  REQUIRE(constMacroMolTemplate.getMainSgroup() == nullptr);

  CHECK(macroMolTemplate.getNumAtoms() == 7);
  CHECK(macroMolTemplate.getAtomWithIdx(0)->getSymbol() == "H");
  CHECK(macroMolTemplate.getAtomWithIdx(1)->getSymbol() == "N");
  CHECK(macroMolTemplate.getAtomWithIdx(2)->getSymbol() == "C");
  CHECK(macroMolTemplate.getAtomWithIdx(3)->getSymbol() == "C");
  CHECK(macroMolTemplate.getAtomWithIdx(4)->getSymbol() == "C");
  CHECK(macroMolTemplate.getAtomWithIdx(5)->getSymbol() == "O");
  CHECK(macroMolTemplate.getAtomWithIdx(6)->getSymbol() == "O");

  macroMolTemplate.setMainGroup({1, 2, 3, 4, 5}, "AA");
  // The amino nitrogen (1) attaches to the leaving hydrogen (0).
  macroMolTemplate.addLeavingGroup({0}, 1, 0, 1);
  // The carbonyl carbon (4) attaches to the leaving hydroxyl oxygen (6).
  macroMolTemplate.addLeavingGroup({6}, 4, 6, 2);

  const auto *mainSgroup = constMacroMolTemplate.getMainSgroup();
  REQUIRE(mainSgroup != nullptr);
  CHECK(mainSgroup->getProp<std::string>("TYPE") == "SUP");
  CHECK(mainSgroup->getProp<std::string>("CLASS") == "AA");
  CHECK(mainSgroup->getAtoms() == std::vector<unsigned int>({1, 2, 3, 4, 5}));

  auto leavingGroups = macroMolTemplate.getLeavingGroups();
  REQUIRE(leavingGroups.size() == 2);
  CHECK(leavingGroups[0]->getProp<std::string>("TYPE") == "SUP");
  CHECK(leavingGroups[0]->getProp<std::string>("CLASS") == "LGRP");
  CHECK(leavingGroups[0]->getAtoms() == std::vector<unsigned int>({0}));
  CHECK(leavingGroups[1]->getProp<std::string>("TYPE") == "SUP");
  CHECK(leavingGroups[1]->getProp<std::string>("CLASS") == "LGRP");
  CHECK(leavingGroups[1]->getAtoms() == std::vector<unsigned int>({6}));

  const auto &attachPoints = mainSgroup->getAttachPoints();
  REQUIRE(attachPoints.size() == 2);
  CHECK(attachPoints[0].aIdx == 1);
  CHECK(attachPoints[0].lvIdx == 0);
  CHECK(attachPoints[0].id == "1");
  CHECK(attachPoints[1].aIdx == 4);
  CHECK(attachPoints[1].lvIdx == 6);
  CHECK(attachPoints[1].id == "2");
}

TEST_CASE("testMacroMolTemplateLibraryMissingTemplate") {
  // Build a MacroMolTemplateLibrary and check that missing lookups return an
  // empty shared pointer.
  MacroMolTemplateLibrary templateLibrary;

  CHECK(!templateLibrary.getByTemplateName("AA", "ALA"));
  CHECK(!templateLibrary.getBySymbol("AA", "A"));
}
