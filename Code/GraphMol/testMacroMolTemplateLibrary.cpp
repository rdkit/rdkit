// Copyright (C) 2026 Schrödinger and other RDKit contributors

#include <GraphMol/MacroMol.h>
#include <GraphMol/MacroMolTemplateLibrary.h>
#include <GraphMol/SubstanceGroup.h>
#include <catch2/catch_all.hpp>

#include <algorithm>
#include <set>
#include <string>
#include <vector>

using namespace RDKit;

TEST_CASE("global MacroMol template library contains valid built-ins") {
  const auto &library = getGlobalMacroMolTemplateLibrary();
  REQUIRE(library.getTemplates().size() == 22);

  std::set<std::string> symbols;
  std::set<std::string> names;
  for (const auto *templ : library.getTemplates()) {
    REQUIRE(templ);
    CHECK(symbols.insert(templ->getSymbol()).second);
    CHECK(names.insert(templ->getName()).second);
    CHECK(templ->getMol().getNumAtoms() > 0);
    CHECK(!templ->getMainAtomIdxs().empty());
    CHECK(getSubstanceGroups(templ->getMol()).size() ==
          templ->getLeavingGroups().size() + 1);
    CHECK(library.getBySymbol(templ->getMonomerClass(), templ->getSymbol()) ==
          templ);
    CHECK(library.getByName(templ->getMonomerClass(), templ->getName()) ==
          templ);
  }
}

TEST_CASE("built-in MacroMol templates can populate a local library") {
  MacroMolTemplateLibrary local;
  addBuiltinMacroMolTemplates(local);
  REQUIRE(local.getTemplates().size() == 22);
}
