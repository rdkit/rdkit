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
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileWriters.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <catch2/catch_all.hpp>

#include <array>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace RDKit;

namespace {

std::unique_ptr<MacroMolTemplate> makeTemplate(
    const std::string &name, const std::string &symbol,
    const std::string &smiles, std::vector<unsigned int> mainGroupAtoms,
    std::vector<MacroMolLeavingGroup> leavingGroups = {},
    MonomerClass monomerClass = MonomerClass::AminoAcid) {
  std::unique_ptr<RWMol> mol(SmilesToMol(smiles));
  if (!mol) {
    throw ValueErrorException("could not parse template SMILES");
  }
  MacroMolTemplateBuilder builder(*mol, monomerClass, name, symbol, smiles);
  builder.setMainGroup(std::move(mainGroupAtoms));
  for (auto &leavingGroup : leavingGroups) {
    builder.addLeavingGroup(std::move(leavingGroup));
  }
  return std::move(builder).build();
}

std::unique_ptr<MacroMolTemplate> makeAnnotatedAlanineTemplate() {
  SmilesParserParams params;
  params.removeHs = false;
  std::unique_ptr<RWMol> parsed(SmilesToMol("[H]N[C@@H](C)C(=O)O", params));
  if (!parsed) {
    throw ValueErrorException("could not parse alanine template SMILES");
  }

  MacroMolTemplateBuilder builder(*parsed, MonomerClass::AminoAcid, "ALA",
                                  "A", "[H]N[C@@H](C)C(=O)O");
  return std::move(builder.setMainGroup({1, 2, 3, 4, 5})
                       .addLeavingGroup({{0}, 1, 0, 1})
                       .addLeavingGroup({{6}, 4, 6, 2}))
      .build();
}

void checkSerializedAlanineTemplate(const ROMol &mol) {
  const auto &sgroups = getSubstanceGroups(mol);
  REQUIRE(sgroups.size() == 3);

  const auto &mainSgroup = sgroups[0];
  CHECK(mainSgroup.getProp<std::string>("TYPE") == "SUP");
  CHECK(mainSgroup.getProp<std::string>("CLASS") == "AminoAcid");
  CHECK(mainSgroup.getProp<std::string>("LABEL") == "ALA");
  CHECK(mainSgroup.getAtoms() ==
        std::vector<unsigned int>({1, 2, 3, 4, 5}));
  const auto &attachPoints = mainSgroup.getAttachPoints();
  REQUIRE(attachPoints.size() == 2);
  const SubstanceGroup::AttachPoint firstAttachPoint{1, 0, "1"};
  const SubstanceGroup::AttachPoint secondAttachPoint{4, 6, "2"};
  CHECK(attachPoints[0] == firstAttachPoint);
  CHECK(attachPoints[1] == secondAttachPoint);

  CHECK(sgroups[1].getProp<std::string>("CLASS") == "LGRP");
  CHECK(sgroups[1].getAtoms() == std::vector<unsigned int>{0});
  CHECK(sgroups[2].getProp<std::string>("CLASS") == "LGRP");
  CHECK(sgroups[2].getAtoms() == std::vector<unsigned int>{6});
}

}  // namespace

static_assert(!std::is_base_of_v<ROMol, MacroMolTemplate>);
static_assert(!std::is_base_of_v<RWMol, MacroMolTemplate>);
static_assert(!std::is_default_constructible_v<MacroMolLeavingGroup>);
static_assert(std::is_constructible_v<MacroMolLeavingGroup,
                                      std::vector<unsigned int>, unsigned int,
                                      unsigned int, int>);
static_assert(!std::is_default_constructible_v<MacroMolTemplate>);
static_assert(std::is_copy_constructible_v<MacroMolTemplate>);
static_assert(std::is_move_constructible_v<MacroMolTemplate>);
static_assert(!std::is_copy_assignable_v<MacroMolTemplate>);
static_assert(std::is_same_v<
              decltype(std::declval<const MacroMolTemplate &>().getMol()),
              const ROMol &>);
static_assert(std::is_same_v<
              decltype(std::declval<const MacroMolTemplate &>()
                           .getMainSgroup()),
              const SubstanceGroup &>);
static_assert(std::is_same_v<
              decltype(std::declval<const MacroMolTemplateLibrary &>()
                           .entries()),
              const std::vector<const MacroMolTemplate *> &>);

TEST_CASE("MacroMolTemplate owns a logically read-only molecule and metadata") {
  auto templ = makeTemplate("ALA", "A", "C", {0});

  CHECK(templ->getMol().getNumAtoms() == 1);
  CHECK(templ->getMonomerClass() == MonomerClass::AminoAcid);
  CHECK(templ->getName() == "ALA");
  CHECK(templ->getSymbol() == "A");
  CHECK(templ->getOriginalData() == "C");
  CHECK(templ->getMainAtomIdxs() == std::vector<unsigned int>{0});
  CHECK(templ->getLeavingGroups().empty());

  MacroMolTemplate copied(*templ);
  CHECK(copied.getMol().getNumAtoms() == 1);
  CHECK(&copied.getMainSgroup().getOwningMol() == &copied.getMol());
}

TEST_CASE("MacroMolTemplate mirrors typed main and leaving groups") {
  auto templ = makeAnnotatedAlanineTemplate();

  CHECK(templ->getMol().getNumAtoms() == 7);
  CHECK(templ->getMainAtomIdxs() ==
        std::vector<unsigned int>({1, 2, 3, 4, 5}));
  const auto &leavingGroups = templ->getLeavingGroups();
  REQUIRE(leavingGroups.size() == 2);
  CHECK(leavingGroups[0].atomIdxs == std::vector<unsigned int>{0});
  CHECK(leavingGroups[0].attachAtomIdx == 1);
  CHECK(leavingGroups[0].leavingAtomIdx == 0);
  CHECK(leavingGroups[0].attachPoint == 1);
  CHECK(leavingGroups[1].atomIdxs == std::vector<unsigned int>{6});

  const auto &mainSgroup = templ->getMainSgroup();
  CHECK(mainSgroup.getProp<std::string>("TYPE") == "SUP");
  CHECK(mainSgroup.getProp<std::string>("CLASS") == "AminoAcid");
  CHECK(mainSgroup.getProp<std::string>("LABEL") == "ALA");
  CHECK(mainSgroup.getAtoms() == templ->getMainAtomIdxs());
  const auto &attachPoints = mainSgroup.getAttachPoints();
  REQUIRE(attachPoints.size() == 2);
  const SubstanceGroup::AttachPoint firstAttachPoint{1, 0, "1"};
  const SubstanceGroup::AttachPoint secondAttachPoint{4, 6, "2"};
  CHECK(attachPoints[0] == firstAttachPoint);
  CHECK(attachPoints[1] == secondAttachPoint);

  const auto &sgroups = getSubstanceGroups(templ->getMol());
  REQUIRE(sgroups.size() == 3);
  CHECK(sgroups[1].getProp<std::string>("CLASS") == "LGRP");
  CHECK(sgroups[1].getAtoms() == leavingGroups[0].atomIdxs);
  CHECK(sgroups[2].getProp<std::string>("CLASS") == "LGRP");
  CHECK(sgroups[2].getAtoms() == leavingGroups[1].atomIdxs);
}

TEST_CASE("Monomer classes are valid generic SGroup classes") {
  const std::array monomerClasses{
      MonomerClass::AminoAcid, MonomerClass::NucleicAcid,
      MonomerClass::Chemical, MonomerClass::Other};

  for (const auto monomerClass : monomerClasses) {
    const std::string className = monomerClassToString(monomerClass);
    CHECK(SubstanceGroupChecks::isValidClass(className));

    auto templ = makeTemplate(className, "X", "C", {0}, {}, monomerClass);
    const auto molBlock = MolToV3KMolBlock(templ->getMol());
    std::unique_ptr<RWMol> roundTripped(
        MolBlockToMol(molBlock, false, false, true));
    REQUIRE(roundTripped);
    const auto &sgroups = getSubstanceGroups(*roundTripped);
    REQUIRE(sgroups.size() == 1);
    CHECK(sgroups[0].getProp<std::string>("CLASS") == className);
  }

  CHECK_FALSE(SubstanceGroupChecks::isValidClass("NotAMonomerClass"));
}

TEST_CASE("MacroMolTemplate SGroups survive generic MOL and SDF round trips") {
  auto templ = makeAnnotatedAlanineTemplate();

  SECTION("V2000 MOL") {
    const auto molBlock = MolToV2KMolBlock(templ->getMol());
    std::unique_ptr<RWMol> roundTripped(
        MolBlockToMol(molBlock, false, false, true));
    REQUIRE(roundTripped);
    checkSerializedAlanineTemplate(*roundTripped);
  }

  SECTION("V3000 MOL") {
    const auto molBlock = MolToV3KMolBlock(templ->getMol());
    std::unique_ptr<RWMol> roundTripped(
        MolBlockToMol(molBlock, false, false, true));
    REQUIRE(roundTripped);
    checkSerializedAlanineTemplate(*roundTripped);
  }

  SECTION("SDF") {
    v2::FileParsers::MolFileParserParams params;
    params.sanitize = false;
    params.removeHs = false;
    params.strictParsing = true;
    v2::FileParsers::SDMolSupplier supplier;
    supplier.setData(SDWriter::getText(templ->getMol()), params);
    auto roundTripped = supplier.next();
    REQUIRE(roundTripped);
    checkSerializedAlanineTemplate(*roundTripped);
  }
}

TEST_CASE("MacroMolTemplateBuilder validates completed definitions") {
  SECTION("main group is required") {
    RWMol mol;
    MacroMolTemplateBuilder builder(mol, MonomerClass::Other, "X", "X", "");
    CHECK_THROWS_AS(std::move(builder).build(), ValueErrorException);
  }
  SECTION("metadata is required") {
    auto mol = std::unique_ptr<RWMol>(SmilesToMol("C"));
    MacroMolTemplateBuilder builder(*mol, MonomerClass::Other, "", "X", "C");
    builder.setMainGroup({0});
    CHECK_THROWS_AS(std::move(builder).build(), ValueErrorException);
  }
  SECTION("main indices must be unique and in range") {
    auto mol = std::unique_ptr<RWMol>(SmilesToMol("CC"));
    MacroMolTemplateBuilder duplicate(*mol, MonomerClass::Other, "X", "X",
                                      "CC");
    duplicate.setMainGroup({0, 0});
    CHECK_THROWS_AS(std::move(duplicate).build(), ValueErrorException);

    MacroMolTemplateBuilder outOfRange(*mol, MonomerClass::Other, "X", "X",
                                       "CC");
    outOfRange.setMainGroup({0, 2});
    CHECK_THROWS_AS(std::move(outOfRange).build(), ValueErrorException);
  }
  SECTION("groups must form a disjoint atom partition") {
    auto mol = std::unique_ptr<RWMol>(SmilesToMol("CCC"));
    MacroMolTemplateBuilder gap(*mol, MonomerClass::Other, "X", "X", "CCC");
    gap.setMainGroup({0, 1});
    CHECK_THROWS_AS(std::move(gap).build(), ValueErrorException);

    MacroMolTemplateBuilder overlap(*mol, MonomerClass::Other, "X", "X",
                                    "CCC");
    overlap.setMainGroup({0, 1}).addLeavingGroup({{1, 2}, 1, 2, 1});
    CHECK_THROWS_AS(std::move(overlap).build(), ValueErrorException);
  }
  SECTION("attachment ids must be positive and unique") {
    auto mol = std::unique_ptr<RWMol>(SmilesToMol("CCC"));
    MacroMolTemplateBuilder duplicate(*mol, MonomerClass::Other, "X", "X",
                                      "CCC");
    duplicate.setMainGroup({1})
        .addLeavingGroup({{0}, 1, 0, 1})
        .addLeavingGroup({{2}, 1, 2, 1});
    CHECK_THROWS_AS(std::move(duplicate).build(), ValueErrorException);
  }
  SECTION("leaving groups must be connected") {
    auto mol = std::unique_ptr<RWMol>(SmilesToMol("CCC.C"));
    MacroMolTemplateBuilder builder(*mol, MonomerClass::Other, "X", "X",
                                    "CCC.C");
    builder.setMainGroup({0}).addLeavingGroup({{1, 2, 3}, 0, 1, 1});
    CHECK_THROWS_AS(std::move(builder).build(), ValueErrorException);
  }
  SECTION("leaving groups have exactly one declared boundary bond") {
    auto mol = std::unique_ptr<RWMol>(SmilesToMol("C1CC1"));
    MacroMolTemplateBuilder builder(*mol, MonomerClass::Other, "X", "X",
                                    "C1CC1");
    builder.setMainGroup({0}).addLeavingGroup({{1, 2}, 0, 1, 1});
    CHECK_THROWS_AS(std::move(builder).build(), ValueErrorException);
  }
}

TEST_CASE("MacroMolTemplateLibrary owns, orders, and looks up templates") {
  MacroMolTemplateLibrary library;
  auto small = makeTemplate("SMALL", "S", "C", {0});
  auto firstLarge = makeTemplate("FIRST_LARGE", "L1", "CCC", {0, 1, 2});
  auto secondLarge = makeTemplate("SECOND_LARGE", "L2", "CCN", {0, 1, 2});
  const auto *smallPtr = small.get();
  const auto *firstLargePtr = firstLarge.get();
  const auto *secondLargePtr = secondLarge.get();

  library.addTemplate(std::move(small));
  library.addTemplate(std::move(firstLarge));
  library.addTemplate(std::move(secondLarge));

  const std::vector<const MacroMolTemplate *> expectedEntries{
      firstLargePtr, secondLargePtr, smallPtr};
  CHECK(library.entries() == expectedEntries);
  CHECK(library.getByName(MonomerClass::AminoAcid, "SMALL") == smallPtr);
  CHECK(library.getBySymbol(MonomerClass::AminoAcid, "S") == smallPtr);
  CHECK(library.getByName(MonomerClass::AminoAcid, "missing") == nullptr);
}

TEST_CASE("MacroMolTemplateLibrary separates classes and rejects duplicates") {
  MacroMolTemplateLibrary library;
  library.addTemplate(
      makeTemplate("ALA", "A", "C", {0}, {}, MonomerClass::AminoAcid));
  library.addTemplate(
      makeTemplate("ADE", "A", "N", {0}, {}, MonomerClass::NucleicAcid));

  CHECK(library.getBySymbol(MonomerClass::AminoAcid, "A")->getName() ==
        "ALA");
  CHECK(library.getBySymbol(MonomerClass::NucleicAcid, "A")->getName() ==
        "ADE");
  CHECK_THROWS_AS(
      library.addTemplate(
          makeTemplate("ALA", "X", "N", {0}, {}, MonomerClass::AminoAcid)),
      ValueErrorException);
  CHECK_THROWS_AS(
      library.addTemplate(
          makeTemplate("OTHER", "A", "N", {0}, {}, MonomerClass::AminoAcid)),
      ValueErrorException);
  CHECK_THROWS_AS(
      library.addTemplate(std::unique_ptr<MacroMolTemplate>{}),
      ValueErrorException);
}
