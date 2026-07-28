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
#include <GraphMol/Atom.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <catch2/catch_all.hpp>

#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace RDKit;

namespace {

std::unique_ptr<MacroMolTemplate> makeTemplate(
    const std::string &templateName, const std::string &symbol,
    const std::string &smiles, const std::vector<unsigned int> &mainGroupAtoms,
    MonomerClass monomerClass = MonomerClass::AminoAcid) {
  std::unique_ptr<RWMol> mol(SmilesToMol(smiles));
  if (!mol) {
    throw ValueErrorException("could not parse template SMILES");
  }
  auto result = std::make_unique<MacroMolTemplate>(
      *mol, monomerClass, templateName, symbol, smiles);
  result->setMainGroup(mainGroupAtoms);
  return result;
}

}  // namespace

static_assert(std::is_base_of_v<RWMol, MacroMolTemplate>);
static_assert(!std::is_default_constructible_v<MacroMolTemplate>);
static_assert(std::is_same_v<
              decltype(std::declval<const MacroMolTemplate &>()
                           .getTemplateName()),
              const std::string &>);
static_assert(std::is_same_v<
              decltype(std::declval<const MacroMolTemplateLibrary &>()
                           .entries()),
              const std::vector<const MacroMolTemplate *> &>);
static_assert(std::is_same_v<
              decltype(std::declval<const MacroMolTemplateLibrary &>()
                           .getByTemplateName(MonomerClass::AminoAcid, "")),
              const MacroMolTemplate *>);

TEST_CASE("MacroMolTemplate owns its molecule and metadata") {
  RWMol mol;
  MacroMolTemplate macroMolTemplate(
      mol, MonomerClass::AminoAcid, "ALA", "A", "N[C@@H](C)C(=O)O");

  CHECK(macroMolTemplate.getNumAtoms() == 0);
  CHECK(macroMolTemplate.getMonomerClass() == MonomerClass::AminoAcid);
  CHECK(macroMolTemplate.getTemplateName() == "ALA");
  CHECK(macroMolTemplate.getSymbol() == "A");
  CHECK(macroMolTemplate.getOriginalData() == "N[C@@H](C)C(=O)O");
}

TEST_CASE("MacroMolTemplateLibrary lookups return const templates") {
  MacroMolTemplateLibrary templateLibrary;
  auto alanine = makeTemplate("ALA", "A", "CC", {0, 1});
  auto cysteine = makeTemplate("CYS", "C", "CS", {0, 1});
  const auto *alaninePtr = alanine.get();
  const auto *cysteinePtr = cysteine.get();

  templateLibrary.addTemplate(std::move(alanine));
  templateLibrary.addTemplate(std::move(cysteine));

  CHECK(templateLibrary.getByTemplateName(MonomerClass::AminoAcid, "ALA") ==
        alaninePtr);
  CHECK(templateLibrary.getBySymbol(MonomerClass::AminoAcid, "A") ==
        alaninePtr);
  CHECK(templateLibrary.getByTemplateName(MonomerClass::AminoAcid, "CYS") ==
        cysteinePtr);
  CHECK(templateLibrary.getBySymbol(MonomerClass::AminoAcid, "C") ==
        cysteinePtr);
}

TEST_CASE("MacroMolTemplateLibrary separates monomer classes") {
  MacroMolTemplateLibrary templateLibrary;
  auto aminoAcid =
      makeTemplate("ALA", "A", "C", {0}, MonomerClass::AminoAcid);
  auto nucleicAcid =
      makeTemplate("ADE", "A", "N", {0}, MonomerClass::NucleicAcid);
  const auto *aminoAcidPtr = aminoAcid.get();
  const auto *nucleicAcidPtr = nucleicAcid.get();

  templateLibrary.addTemplate(std::move(aminoAcid));
  templateLibrary.addTemplate(std::move(nucleicAcid));

  CHECK(templateLibrary.getByTemplateName(MonomerClass::AminoAcid, "ALA") ==
        aminoAcidPtr);
  CHECK(templateLibrary.getBySymbol(MonomerClass::AminoAcid, "A") ==
        aminoAcidPtr);
  CHECK(templateLibrary.getByTemplateName(MonomerClass::NucleicAcid, "ADE") ==
        nucleicAcidPtr);
  CHECK(templateLibrary.getBySymbol(MonomerClass::NucleicAcid, "A") ==
        nucleicAcidPtr);
}

TEST_CASE("MacroMolTemplate main and leaving groups") {
  // Build an alanine template with explicit peptide-bond leaving groups:
  //   [H] on the N side and [OH] on the C-terminal O side. Atom indices:
  //   0:H(N)  1:N  2:C(alpha)  3:C(methyl)  4:C(carbonyl)  5:O(=O)  6:O(H)
  SmilesParserParams params;
  params.removeHs = false;
  std::unique_ptr<RWMol> parsed(SmilesToMol("[H]N[C@@H](C)C(=O)O", params));
  REQUIRE(parsed);
  MacroMolTemplate macroMolTemplate(*parsed, MonomerClass::AminoAcid, "ALA",
                                    "A", "[H]N[C@@H](C)C(=O)O");
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

  macroMolTemplate.setMainGroup({1, 2, 3, 4, 5});
  // The amino nitrogen (1) attaches to the leaving hydrogen (0).
  macroMolTemplate.addLeavingGroup({0}, 1, 0, 1);
  // The carbonyl carbon (4) attaches to the leaving hydroxyl oxygen (6).
  macroMolTemplate.addLeavingGroup({6}, 4, 6, 2);

  const auto *mainSgroup = constMacroMolTemplate.getMainSgroup();
  REQUIRE(mainSgroup != nullptr);
  CHECK(mainSgroup->getProp<std::string>("TYPE") == "SUP");
  CHECK(mainSgroup->getProp<std::string>("CLASS") == "AminoAcid");
  CHECK(mainSgroup->getAtoms() == std::vector<unsigned int>({1, 2, 3, 4, 5}));

  auto leavingGroups = constMacroMolTemplate.getLeavingGroups();
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

TEST_CASE("MacroMolTemplateLibrary orders templates by main-group size") {
  MacroMolTemplateLibrary templateLibrary;
  auto small = makeTemplate("SMALL", "S", "CC", {0});
  auto firstLarge = makeTemplate("FIRST_LARGE", "L1", "CCC", {0, 1, 2});
  auto secondLarge = makeTemplate("SECOND_LARGE", "L2", "CCN", {0, 1, 2});
  const auto *smallPtr = small.get();
  const auto *firstLargePtr = firstLarge.get();
  const auto *secondLargePtr = secondLarge.get();

  templateLibrary.addTemplate(std::move(small));
  templateLibrary.addTemplate(std::move(firstLarge));
  templateLibrary.addTemplate(std::move(secondLarge));

  const auto &entries = templateLibrary.entries();
  REQUIRE(entries.size() == 3);
  CHECK(entries[0] == firstLargePtr);
  CHECK(entries[1] == secondLargePtr);
  CHECK(entries[2] == smallPtr);
}

TEST_CASE("MacroMolTemplateLibrary rejects invalid templates") {
  SECTION("null template") {
    MacroMolTemplateLibrary templateLibrary;
    CHECK_THROWS_AS(
        templateLibrary.addTemplate(std::unique_ptr<MacroMolTemplate>{}),
        ValueErrorException);
  }

  SECTION("missing main group") {
    MacroMolTemplateLibrary templateLibrary;
    RWMol mol;
    auto macroMolTemplate = std::make_unique<MacroMolTemplate>(
        mol, MonomerClass::AminoAcid, "ALA", "A", "");
    CHECK_THROWS_AS(templateLibrary.addTemplate(std::move(macroMolTemplate)),
                    ValueErrorException);
  }

  SECTION("multiple main groups") {
    MacroMolTemplateLibrary templateLibrary;
    auto macroMolTemplate = makeTemplate("ALA", "A", "CC", {0});
    SubstanceGroup secondMainGroup(macroMolTemplate.get(), "SUP");
    secondMainGroup.setProp("CLASS", std::string("AminoAcid"));
    secondMainGroup.setAtoms({1});
    addSubstanceGroup(*macroMolTemplate, secondMainGroup);

    CHECK_THROWS_AS(templateLibrary.addTemplate(std::move(macroMolTemplate)),
                    ValueErrorException);
  }
}

TEST_CASE("MacroMolTemplateLibrary rejects duplicate lookup keys") {
  SECTION("duplicate template name") {
    MacroMolTemplateLibrary templateLibrary;
    templateLibrary.addTemplate(makeTemplate("ALA", "A", "C", {0}));

    CHECK_THROWS_AS(
        templateLibrary.addTemplate(makeTemplate("ALA", "X", "N", {0})),
        ValueErrorException);
  }

  SECTION("duplicate symbol") {
    MacroMolTemplateLibrary templateLibrary;
    templateLibrary.addTemplate(makeTemplate("ALA", "A", "C", {0}));

    CHECK_THROWS_AS(
        templateLibrary.addTemplate(makeTemplate("OTHER", "A", "N", {0})),
        ValueErrorException);
  }
}

TEST_CASE("MacroMolTemplateLibrary missing lookups return nullptr") {
  MacroMolTemplateLibrary templateLibrary;

  CHECK(templateLibrary.getByTemplateName(MonomerClass::AminoAcid, "ALA") ==
        nullptr);
  CHECK(templateLibrary.getBySymbol(MonomerClass::AminoAcid, "A") == nullptr);
}
