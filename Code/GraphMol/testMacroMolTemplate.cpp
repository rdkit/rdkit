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
#include <catch2/catch_all.hpp>

#include <memory>

using namespace RDKit;

TEST_CASE("testBuildMacroMolTemplate") {
  // Build a simple MacroMolTemplate and check that the public fields are
  // populated as expected.
  MacroMolTemplate macroMolTemplate;
  macroMolTemplate.monomerClass = "AA";
  macroMolTemplate.templateName = "ALA";
  macroMolTemplate.symbol = "A";
  macroMolTemplate.original_data = "N[C@@H](C)C(=O)O";
  macroMolTemplate.mol = std::make_shared<RWMol>();

  CHECK(macroMolTemplate.monomerClass == "AA");
  CHECK(macroMolTemplate.templateName == "ALA");
  CHECK(macroMolTemplate.symbol == "A");
  CHECK(macroMolTemplate.original_data == "N[C@@H](C)C(=O)O");
  CHECK(macroMolTemplate.mol);
}

TEST_CASE("testMacroMolTemplateLibraryLookup") {
  // Build a MacroMolTemplateLibrary with two amino-acid templates and check
  // that lookup by template name and symbol returns the expected template.
  MacroMolTemplateLibrary templateLibrary;
  auto alanineTemplate = std::make_shared<MacroMolTemplate>();
  alanineTemplate->monomerClass = "AA";
  alanineTemplate->templateName = "ALA";
  alanineTemplate->symbol = "A";
  alanineTemplate->mol = std::make_shared<RWMol>();

  auto cysteineTemplate = std::make_shared<MacroMolTemplate>();
  cysteineTemplate->monomerClass = "AA";
  cysteineTemplate->templateName = "CYS";
  cysteineTemplate->symbol = "C";
  cysteineTemplate->mol = std::make_shared<RWMol>();

  templateLibrary.addTemplate(alanineTemplate);
  templateLibrary.addTemplate(cysteineTemplate);

  CHECK(templateLibrary.getByTemplateName("AA", "ALA") == alanineTemplate);
  CHECK(templateLibrary.getBySymbol("AA", "A") == alanineTemplate);
  CHECK(templateLibrary.getByTemplateName("AA", "CYS") == cysteineTemplate);
  CHECK(templateLibrary.getBySymbol("AA", "C") == cysteineTemplate);
}

TEST_CASE("testMacroMolTemplateLibrarySeparatesMonomerClasses") {
  // Build a MacroMolTemplateLibrary with amino-acid and nucleic-acid templates
  // sharing the same template name and symbol. Check that the monomer class is
  // part of the lookup key.
  MacroMolTemplateLibrary templateLibrary;
  auto aaTemplate = std::make_shared<MacroMolTemplate>();
  aaTemplate->monomerClass = "AA";
  aaTemplate->templateName = "ALA";
  aaTemplate->symbol = "A";
  aaTemplate->mol = std::make_shared<RWMol>();

  auto naTemplate = std::make_shared<MacroMolTemplate>();
  naTemplate->monomerClass = "NA";
  naTemplate->templateName = "ADE";
  naTemplate->symbol = "A";
  naTemplate->mol = std::make_shared<RWMol>();

  templateLibrary.addTemplate(aaTemplate);
  templateLibrary.addTemplate(naTemplate);

  CHECK(templateLibrary.getByTemplateName("AA", "ALA") == aaTemplate);
  CHECK(templateLibrary.getBySymbol("AA", "A") == aaTemplate);
  CHECK(templateLibrary.getByTemplateName("NA", "ADE") == naTemplate);
  CHECK(templateLibrary.getBySymbol("NA", "A") == naTemplate);
}

TEST_CASE("testMacroMolTemplateLibraryMissingTemplate") {
  // Build a MacroMolTemplateLibrary and check that missing lookups return an
  // empty shared pointer.
  MacroMolTemplateLibrary templateLibrary;

  CHECK(!templateLibrary.getByTemplateName("AA", "ALA"));
  CHECK(!templateLibrary.getBySymbol("AA", "A"));
}
