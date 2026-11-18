// Copyright (C) 2026 Schrödinger and other RDKit contributors

#include <GraphMol/MacroMol.h>
#include <GraphMol/MacroMolTemplateLibrary.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SubstanceGroup.h>
#include <catch2/catch_all.hpp>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace RDKit;

TEST_CASE("global MacroMol template library contains valid built-ins") {
  const auto &library = getGlobalMacroMolTemplateLibrary();
  REQUIRE(library.getTemplates().size() == 70);

  std::set<std::pair<MonomerClass, std::string>> symbols;
  std::set<std::pair<MonomerClass, std::string>> names;
  for (const auto *templ : library.getTemplates()) {
    REQUIRE(templ);
    const auto monomerClass = templ->getMonomerClass();
    CHECK(symbols.insert({monomerClass, templ->getSymbol()}).second);
    CHECK(names.insert({monomerClass, templ->getName()}).second);
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

TEST_CASE("built-in alanine template preserves its definition") {
  const auto *templ = getGlobalMacroMolTemplateLibrary().getBySymbol(
      MonomerClass::AminoAcid, "A");
  REQUIRE(templ);

  CHECK(templ->getMonomerClass() == MonomerClass::AminoAcid);
  CHECK(templ->getName() == "ALA");
  CHECK(templ->getSymbol() == "A");
  CHECK(templ->getSubclass() == "Backbone");
  // the original data is the JSON definition of the monomer
  boost::property_tree::ptree definition;
  std::istringstream originalData(templ->getOriginalData());
  boost::property_tree::read_json(originalData, definition);
  CHECK(definition.get<std::string>("SYMBOL") == "A");
  CHECK(definition.get<std::string>("PDBCODE") == "ALA");
  CHECK(definition.get<std::string>("CORE_SMILES") ==
        "*N[C@@H](C)C(*)=O |$_R1;;;;;_R2;$,atomProp:1.pdbName. N  :2.pdbName. "
        "CA :3.pdbName. CB :4.pdbName. C  :6.pdbName. O  |");
  CHECK(definition.get<std::string>("RGROUPS.R2") ==
        "*O |$_R2;$,atomProp:1.pdbName. OXT|");

  // Atoms follow CORE_SMILES, with each R-group atom taking the place of its
  // dummy atom
  const auto &mol = templ->getMol();
  REQUIRE(mol.getNumAtoms() == 7);
  CHECK(MolToSmiles(mol) == "[H]N[C@@H](C)C(=O)O");
  CHECK(mol.getAtomWithIdx(0)->getProp<std::string>("pdbName") == " H  ");
  CHECK(mol.getAtomWithIdx(2)->getProp<std::string>("pdbName") == " CA ");
  CHECK(mol.getAtomWithIdx(5)->getProp<std::string>("pdbName") == " OXT");
  const std::vector<unsigned int> expectedMainAtoms{1, 2, 3, 4, 6};
  CHECK(templ->getMainAtomIdxs() == expectedMainAtoms);

  const auto &leavingGroups = templ->getLeavingGroups();
  REQUIRE(leavingGroups.size() == 2);
  const std::vector<unsigned int> expectedFirstLeavingAtoms{0};
  const std::vector<unsigned int> expectedSecondLeavingAtoms{5};
  CHECK(leavingGroups[0].atomIdxs == expectedFirstLeavingAtoms);
  CHECK(leavingGroups[0].attachAtomIdx == 1);
  CHECK(leavingGroups[0].leavingAtomIdx == 0);
  CHECK(leavingGroups[0].attachPoint == 1);
  CHECK(leavingGroups[1].atomIdxs == expectedSecondLeavingAtoms);
  CHECK(leavingGroups[1].attachAtomIdx == 4);
  CHECK(leavingGroups[1].leavingAtomIdx == 5);
  CHECK(leavingGroups[1].attachPoint == 2);
}

TEST_CASE("built-in MacroMol templates can populate a local library") {
  MacroMolTemplateLibrary local;
  addBuiltinMacroMolTemplates(local);
  REQUIRE(local.getTemplates().size() == 70);
}
