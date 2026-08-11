// Copyright (C) 2026 Schrödinger and other RDKit contributors

#include <GraphMol/FileParsers/MolToMacroMol.h>
#include <GraphMol/MacroMolTemplate.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SubstanceGroup.h>
#include <catch2/catch_all.hpp>

#include <memory>
#include <string>
#include <vector>

using namespace RDKit;

namespace {

std::unique_ptr<MacroMolTemplate> makeTemplate(
    const std::string &name, const std::string &symbol,
    const std::string &smiles, const std::vector<unsigned int> &mainAtoms,
    const std::vector<MacroMolLeavingGroup> &leavingGroups = {}) {
  auto parsed = std::unique_ptr<RWMol>(SmilesToMol(smiles));
  MacroMolTemplateBuilder builder(*parsed, MonomerClass::Other, name, symbol,
                                  smiles);
  builder.setMainGroup(mainAtoms);
  for (const auto &leavingGroup : leavingGroups) {
    builder.addLeavingGroup(leavingGroup);
  }
  return std::move(builder).build();
}

}  // namespace

TEST_CASE("MolToMacroMol annotates complete atomistic molecules") {
  MacroMolTemplateLibrary templates;
  templates.addTemplate(
      makeTemplate("ETHYL", "Et", "CCO", {0, 1}, {{{2}, 1, 2, 1}}));

  auto mol = std::unique_ptr<RWMol>(SmilesToMol("CCC"));
  auto grouped = MolToMacroMol(*mol, templates);

  REQUIRE(grouped);
  CHECK(grouped->getNumAtoms() == mol->getNumAtoms());
  CHECK(grouped->getNumBonds() == mol->getNumBonds());
  const auto &groups = getSubstanceGroups(*grouped);
  REQUIRE(groups.size() == 1);
  CHECK(groups[0].getProp<std::string>("CLASS") == "Other");
  CHECK(groups[0].getProp<std::string>("LABEL") == "Et");
  CHECK(groups[0].getAtoms() == std::vector<unsigned int>({0, 1}));
  REQUIRE(groups[0].getAttachPoints().size() == 1);
  CHECK(groups[0].getAttachPoints()[0].aIdx == 1);
  CHECK(groups[0].getAttachPoints()[0].lvIdx == -1);
  CHECK(groups[0].getAttachPoints()[0].id == "1");
}

TEST_CASE("MolToMacroMol keeps deterministic largest-first grouping") {
  MacroMolTemplateLibrary templates;
  templates.addTemplate(makeTemplate("SMALL", "S", "CC", {0, 1}));
  templates.addTemplate(makeTemplate("LARGE", "L", "CCC", {0, 1, 2}));

  auto mol = std::unique_ptr<RWMol>(SmilesToMol("CCC"));
  auto grouped = MolToMacroMol(*mol, templates);

  REQUIRE(grouped);
  const auto &groups = getSubstanceGroups(*grouped);
  REQUIRE(groups.size() == 1);
  CHECK(groups[0].getProp<std::string>("LABEL") == "L");
}

TEST_CASE("CollapseMacroMol collapses a grouped unit and preserves neighbors") {
  MacroMolTemplateLibrary templates;
  templates.addTemplate(
      makeTemplate("ETHYL", "Et", "CCO", {0, 1}, {{{2}, 1, 2, 1}}));
  auto mol = std::unique_ptr<RWMol>(SmilesToMol("CCC"));
  auto grouped = MolToMacroMol(*mol, templates);

  auto collapsed = CollapseMacroMol(*grouped, templates);
  REQUIRE(collapsed);
  CHECK(collapsed->getNumAtoms() == 2);
  CHECK(collapsed->getNumBonds() == 1);
  const auto *macroInfo = collapsed->getAtomWithIdx(0)->getMacroAtomInfo();
  REQUIRE(macroInfo);
  CHECK(macroInfo->getSymbol() == "Et");
  const auto *bondInfo = collapsed->getBondWithIdx(0)->getMacroBondInfo();
  REQUIRE(bondInfo);
  CHECK(bondInfo->getBond(0).beginAttachPt == 1);
  CHECK(bondInfo->getBond(0).endAttachPt == -1);
}

TEST_CASE("CollapseMacroMol creates macro bonds between grouped units") {
  MacroMolTemplateLibrary templates;
  templates.addTemplate(makeTemplate("UNIT", "X", "NCON", {1, 2},
                                      {{{0}, 1, 0, 1}, {{3}, 2, 3, 2}}));
  auto mol = std::unique_ptr<RWMol>(SmilesToMol("COCO"));
  auto grouped = MolToMacroMol(*mol, templates);
  REQUIRE(getSubstanceGroups(*grouped).size() == 2);

  auto collapsed = CollapseMacroMol(*grouped, templates);
  REQUIRE(collapsed);
  CHECK(collapsed->getNumAtoms() == 2);
  REQUIRE(collapsed->getNumBonds() == 1);
  const auto *info = collapsed->getBondWithIdx(0)->getMacroBondInfo();
  REQUIRE(info);
  CHECK(info->getBond(0).beginAttachPt == 2);
  CHECK(info->getBond(0).endAttachPt == 1);
}

TEST_CASE("MolToMacroMol leaves unmatched molecules atomistic") {
  MacroMolTemplateLibrary templates;
  auto mol = std::unique_ptr<RWMol>(SmilesToMol("C/C=C/C"));
  mol->getBondBetweenAtoms(1, 2)->setProp<int>("customTag", 7);
  auto grouped = MolToMacroMol(*mol, templates);

  REQUIRE(grouped);
  CHECK(getSubstanceGroups(*grouped).empty());
  CHECK(grouped->getNumAtoms() == 4);
  CHECK(grouped->getNumBonds() == 3);
  CHECK(grouped->getBondBetweenAtoms(1, 2)->getStereo() ==
        mol->getBondBetweenAtoms(1, 2)->getStereo());
  CHECK(grouped->getBondBetweenAtoms(1, 2)->getProp<int>("customTag") == 7);
}
