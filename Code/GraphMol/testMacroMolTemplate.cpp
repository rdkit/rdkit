//
// Copyright (C) 2026 Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/Atom.h>
#include <GraphMol/MacroMolTemplate.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Invariant.h>
#include <catch2/catch_all.hpp>

#include <memory>
#include <string>
#include <type_traits>
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
  return builder.build();
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
  builder.setMainGroup({1, 2, 3, 4, 5})
      .addLeavingGroup({{0}, 1, 0, 1})
      .addLeavingGroup({{6}, 4, 6, 2});
  return builder.build();
}

}  // namespace

// MacroMolTemplate owns an ROMol instead of extending the molecule hierarchy.
static_assert(!std::is_base_of_v<ROMol, MacroMolTemplate>);

TEST_CASE("MacroMolTemplate owns a logically read-only molecule and metadata") {
  auto templ = makeTemplate("ALA", "A", "C", {0});

  CHECK(templ->getMol().getNumAtoms() == 1);
  CHECK(templ->getMonomerClass() == MonomerClass::AminoAcid);
  CHECK(templ->getName() == "ALA");
  CHECK(templ->getSymbol() == "A");
  CHECK(templ->getOriginalData() == "C");
  CHECK(templ->getSubclass().empty());
  CHECK(templ->getMainAtomIdxs() == std::vector<unsigned int>{0});
  CHECK(templ->getLeavingGroups().empty());

  // Copying must repair the SGroup's owning-molecule back-reference.
  MacroMolTemplate copied(*templ);
  CHECK(copied.getMol().getNumAtoms() == 1);
  CHECK(&copied.getMainSgroup().getOwningMol() == &copied.getMol());
}

TEST_CASE("MacroMolTemplate preserves its SCSR subclass") {
  RWMol mol;
  mol.addAtom(new Atom(6), false, true);
  MacroMolTemplateBuilder builder(mol, MonomerClass::AminoAcid, "ALA", "A",
                                  "C");
  builder.setMainGroup({0}).setSubclass("AA");

  auto templ = builder.build();

  CHECK(templ->getSubclass() == "AA");
  MacroMolTemplate copied(*templ);
  CHECK(copied.getSubclass() == "AA");
}

TEST_CASE("MacroMolTemplateBuilder provides mutable molecule access") {
  RWMol mol;
  mol.addAtom(new Atom(6), false, true);
  MacroMolTemplateBuilder builder(mol, MonomerClass::Chemical, "CARBON", "C",
                                  "C");

  builder.getMol().getAtomWithIdx(0)->setIsotope(13);
  builder.setMainGroup({0});
  auto templ = builder.build();

  CHECK(templ->getMol().getAtomWithIdx(0)->getIsotope() == 13);
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

TEST_CASE("MacroMolTemplateBuilder validates completed definitions") {
  SECTION("main group can only be set once") {
    auto mol = std::unique_ptr<RWMol>(SmilesToMol("C"));
    MacroMolTemplateBuilder builder(*mol, MonomerClass::Other, "X", "X",
                                    "C");
    builder.setMainGroup({0});
    CHECK_THROWS_AS(builder.setMainGroup({0}), Invar::Invariant);
  }
  SECTION("main group is required") {
    RWMol mol;
    MacroMolTemplateBuilder builder(mol, MonomerClass::Other, "X", "X", "");
    CHECK_THROWS_AS(builder.build(), ValueErrorException);
  }
  SECTION("metadata is required") {
    auto mol = std::unique_ptr<RWMol>(SmilesToMol("C"));
    MacroMolTemplateBuilder builder(*mol, MonomerClass::Other, "", "X", "C");
    builder.setMainGroup({0});
    CHECK_THROWS_AS(builder.build(), ValueErrorException);
  }
  SECTION("main indices must be unique and in range") {
    auto mol = std::unique_ptr<RWMol>(SmilesToMol("CC"));
    MacroMolTemplateBuilder duplicate(*mol, MonomerClass::Other, "X", "X",
                                      "CC");
    duplicate.setMainGroup({0, 0});
    CHECK_THROWS_AS(duplicate.build(), ValueErrorException);

    MacroMolTemplateBuilder outOfRange(*mol, MonomerClass::Other, "X", "X",
                                       "CC");
    outOfRange.setMainGroup({0, 2});
    CHECK_THROWS_AS(outOfRange.build(), ValueErrorException);
  }
  SECTION("groups must form a disjoint atom partition") {
    auto mol = std::unique_ptr<RWMol>(SmilesToMol("CCC"));
    MacroMolTemplateBuilder gap(*mol, MonomerClass::Other, "X", "X", "CCC");
    gap.setMainGroup({0, 1});
    CHECK_THROWS_AS(gap.build(), ValueErrorException);

    MacroMolTemplateBuilder overlap(*mol, MonomerClass::Other, "X", "X",
                                    "CCC");
    overlap.setMainGroup({0, 1}).addLeavingGroup({{1, 2}, 1, 2, 1});
    CHECK_THROWS_AS(overlap.build(), ValueErrorException);
  }
  SECTION("attachment ids must be positive and unique") {
    auto mol = std::unique_ptr<RWMol>(SmilesToMol("CCC"));
    MacroMolTemplateBuilder duplicate(*mol, MonomerClass::Other, "X", "X",
                                      "CCC");
    duplicate.setMainGroup({1})
        .addLeavingGroup({{0}, 1, 0, 1})
        .addLeavingGroup({{2}, 1, 2, 1});
    CHECK_THROWS_AS(duplicate.build(), ValueErrorException);
  }
  SECTION("the template molecule must be connected") {
    auto mol = std::unique_ptr<RWMol>(SmilesToMol("CCC.C"));
    MacroMolTemplateBuilder builder(*mol, MonomerClass::Other, "X", "X",
                                    "CCC.C");
    builder.setMainGroup({0}).addLeavingGroup({{1, 2, 3}, 0, 1, 1});
    CHECK_THROWS_AS(builder.build(), ValueErrorException);
  }
  SECTION("leaving groups have exactly one declared boundary bond") {
    auto mol = std::unique_ptr<RWMol>(SmilesToMol("C1CC1"));
    MacroMolTemplateBuilder builder(*mol, MonomerClass::Other, "X", "X",
                                    "C1CC1");
    builder.setMainGroup({0}).addLeavingGroup({{1, 2}, 0, 1, 1});
    CHECK_THROWS_AS(builder.build(), ValueErrorException);
  }
}

TEST_CASE("MacroMolTemplateBuilder can build without being consumed") {
  auto mol = std::unique_ptr<RWMol>(SmilesToMol("C"));
  MacroMolTemplateBuilder builder(*mol, MonomerClass::Other, "X", "X", "C");
  builder.setMainGroup({0});

  const auto first = builder.build();
  const auto second = builder.build();
  CHECK(first->getName() == second->getName());
}

TEST_CASE("MacroMolTemplateLibrary owns and looks up templates") {
  MacroMolTemplateLibrary library;
  auto small = makeTemplate("SMALL", "S", "C", {0});
  const auto *smallPtr = small.get();

  library.addTemplate(std::move(small));

  CHECK(library.getByName(MonomerClass::AminoAcid, "SMALL") == smallPtr);
  CHECK(library.getBySymbol(MonomerClass::AminoAcid, "S") == smallPtr);
  CHECK(library.getByName(MonomerClass::AminoAcid, "missing") == nullptr);
}

TEST_CASE("MacroMolTemplateLibrary separates classes and rejects duplicates") {
  // Name and symbol keys are scoped by MonomerClass.
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
      Invar::Invariant);
}
