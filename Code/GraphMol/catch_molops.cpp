//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/test_fixtures.h>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("github #7556: chiral sulfur in conjugated rings") {
  SECTION("as reported") {
    auto m = "CC1=CC(Cl)=CC2=C1N=[S@](C)N=C2N"_smiles;
    REQUIRE(m);
    CHECK(!m->getBondBetweenAtoms(8, 9)->getIsConjugated());
    CHECK(!m->getBondBetweenAtoms(9, 11)->getIsConjugated());
    REQUIRE(m->getAtomWithIdx(9)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  }
}

TEST_CASE("getAvgMolWt") {
  SECTION("basics") {
    auto mol = "C"_smiles;
    REQUIRE(mol);
    auto amw = MolOps::getAvgMolWt(*mol);
    CHECK_THAT(amw, Catch::Matchers::WithinAbs(16.043, 0.001));
    amw = MolOps::getAvgMolWt(*mol, true);
    CHECK_THAT(amw, Catch::Matchers::WithinAbs(12.011, 0.001));
    MolOps::addHs(*mol);
    amw = MolOps::getAvgMolWt(*mol);
    CHECK_THAT(amw, Catch::Matchers::WithinAbs(16.043, 0.001));
    amw = MolOps::getAvgMolWt(*mol, true);
    CHECK_THAT(amw, Catch::Matchers::WithinAbs(12.011, 0.001));
  }
  SECTION("Hs in SMILES") {
    auto mol = "[CH4]"_smiles;
    REQUIRE(mol);
    auto amw = MolOps::getAvgMolWt(*mol);
    CHECK_THAT(amw, Catch::Matchers::WithinAbs(16.043, 0.001));
    amw = MolOps::getAvgMolWt(*mol, true);
    CHECK_THAT(amw, Catch::Matchers::WithinAbs(12.011, 0.001));
  }
  SECTION("isotopes") {
    auto mol = "C[2H]"_smiles;
    REQUIRE(mol);
    auto amw = MolOps::getAvgMolWt(*mol);
    CHECK_THAT(amw, Catch::Matchers::WithinAbs(17.0, 0.1));
  }
}

TEST_CASE("getExactMolWt") {
  SECTION("basics") {
    auto mol = "C"_smiles;
    REQUIRE(mol);
    auto mw = MolOps::getExactMolWt(*mol);
    CHECK_THAT(mw, Catch::Matchers::WithinAbs(16.031, 0.001));
    mw = MolOps::getExactMolWt(*mol, true);
    CHECK_THAT(mw, Catch::Matchers::WithinAbs(12.000, 0.001));
    MolOps::addHs(*mol);
    mw = MolOps::getExactMolWt(*mol);
    CHECK_THAT(mw, Catch::Matchers::WithinAbs(16.031, 0.001));
    mw = MolOps::getExactMolWt(*mol, true);
    CHECK_THAT(mw, Catch::Matchers::WithinAbs(12.000, 0.001));
  }
  SECTION("Hs in SMILES") {
    auto mol = "[CH4]"_smiles;
    REQUIRE(mol);
    auto mw = MolOps::getExactMolWt(*mol);
    CHECK_THAT(mw, Catch::Matchers::WithinAbs(16.031, 0.001));
    mw = MolOps::getExactMolWt(*mol, true);
    CHECK_THAT(mw, Catch::Matchers::WithinAbs(12.000, 0.001));
  }
  SECTION("isotopes") {
    auto mol = "C[2H]"_smiles;
    REQUIRE(mol);
    auto mw = MolOps::getExactMolWt(*mol);
    CHECK_THAT(mw, Catch::Matchers::WithinAbs(17.037, 0.001));
    mw = MolOps::getExactMolWt(*mol, true);
    CHECK_THAT(mw, Catch::Matchers::WithinAbs(12.000, 0.001));
  }
  SECTION("Cl") {
    auto mol = "Cl"_smiles;
    REQUIRE(mol);
    auto mw = MolOps::getExactMolWt(*mol);
    CHECK_THAT(mw, Catch::Matchers::WithinAbs(34.9688 + 1.0078, 0.001));
    mw = MolOps::getAvgMolWt(*mol);
    CHECK_THAT(mw, Catch::Matchers::WithinAbs(35.453 + 1.008, 0.001));

    mol = "[35ClH]"_smiles;
    REQUIRE(mol);
    mw = MolOps::getExactMolWt(*mol);
    CHECK_THAT(mw, Catch::Matchers::WithinAbs(34.9688 + 1.0078, 0.001));
    mw = MolOps::getAvgMolWt(*mol);
    CHECK_THAT(mw, Catch::Matchers::WithinAbs(34.9688 + 1.008, 0.001));

    mol = "[36ClH]"_smiles;
    REQUIRE(mol);
    mw = MolOps::getExactMolWt(*mol);
    CHECK_THAT(mw, Catch::Matchers::WithinAbs(35.9683 + 1.0078, 0.001));
    mw = MolOps::getAvgMolWt(*mol);
    CHECK_THAT(mw, Catch::Matchers::WithinAbs(35.9683 + 1.008, 0.001));
  }
}

TEST_CASE("getMolFormula") {
  SECTION("basics") {
    auto mol = "C"_smiles;
    REQUIRE(mol);
    auto formula = MolOps::getMolFormula(*mol);
    CHECK(formula == "CH4");
    MolOps::addHs(*mol);
    formula = MolOps::getMolFormula(*mol);
    CHECK(formula == "CH4");

    mol = "[CH4]"_smiles;
    REQUIRE(mol);
    formula = MolOps::getMolFormula(*mol);
    CHECK(formula == "CH4");

    mol = "CO"_smiles;
    REQUIRE(mol);
    formula = MolOps::getMolFormula(*mol);
    CHECK(formula == "CH4O");
    MolOps::addHs(*mol);
    formula = MolOps::getMolFormula(*mol);
    CHECK(formula == "CH4O");

    mol = "C(=O)N"_smiles;
    REQUIRE(mol);
    formula = MolOps::getMolFormula(*mol);
    CHECK(formula == "CH3NO");

    mol = "C(=O)=O"_smiles;
    REQUIRE(mol);
    formula = MolOps::getMolFormula(*mol);
    CHECK(formula == "CO2");

    mol = "C(=O)[O-]"_smiles;
    REQUIRE(mol);
    formula = MolOps::getMolFormula(*mol);
    CHECK(formula == "CHO2-");

    mol = "C([O-])[O-]"_smiles;
    REQUIRE(mol);
    formula = MolOps::getMolFormula(*mol);
    CHECK(formula == "CH2O2-2");

    mol = "C([NH3+])[O-]"_smiles;
    REQUIRE(mol);
    formula = MolOps::getMolFormula(*mol);
    CHECK(formula == "CH5NO");

    mol = "C([NH3+])O"_smiles;
    REQUIRE(mol);
    formula = MolOps::getMolFormula(*mol);
    CHECK(formula == "CH6NO+");

    mol = "C([NH3+])[NH3+]"_smiles;
    REQUIRE(mol);
    formula = MolOps::getMolFormula(*mol);
    CHECK(formula == "CH8N2+2");
  }
  SECTION("H isotopes") {
    auto mol = "[2H]C([3H])O"_smiles;
    REQUIRE(mol);
    auto formula = MolOps::getMolFormula(*mol);
    CHECK(formula == "CH4O");
    formula = MolOps::getMolFormula(*mol, true);
    CHECK(formula == "CH2DTO");
    formula = MolOps::getMolFormula(*mol, true, false);
    CHECK(formula == "CH2[2H][3H]O");
  }

  SECTION("isotopes") {
    auto mol = "[13CH3]C[13CH2]C"_smiles;
    REQUIRE(mol);
    auto formula = MolOps::getMolFormula(*mol);
    CHECK(formula == "C4H10");
    formula = MolOps::getMolFormula(*mol, true);
    CHECK(formula == "C2[13C]2H10");
    formula = MolOps::getMolFormula(*mol, true, false);
    CHECK(formula == "C2[13C]2H10");

    mol = "[13CH3]C([2H])O"_smiles;
    REQUIRE(mol);
    formula = MolOps::getMolFormula(*mol);
    CHECK(formula == "C2H6O");
    formula = MolOps::getMolFormula(*mol, true);
    CHECK(formula == "C[13C]H5DO");
    formula = MolOps::getMolFormula(*mol, true, false);
    CHECK(formula == "C[13C]H5[2H]O");

    mol = "[13CH3]C[13CH2]CB(O)O[2H]"_smiles;
    REQUIRE(mol);
    formula = MolOps::getMolFormula(*mol);
    CHECK(formula == "C4H11BO2");
    formula = MolOps::getMolFormula(*mol, true);
    CHECK(formula == "C2[13C]2H10DBO2");
    formula = MolOps::getMolFormula(*mol, true, false);
    CHECK(formula == "C2[13C]2H10[2H]BO2");
  }
}
