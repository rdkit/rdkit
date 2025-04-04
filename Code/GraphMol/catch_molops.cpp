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

#include <GraphMol/FileParsers/FileParsers.h>
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

TEST_CASE(
    "github #8121: symmetric ring finding not returning correct results for molecules with fragments") {
  auto twoCubanes = "C12C3C4C1C5C2C3C45.C12C3C4C1C5C2C3C45"_smiles;
  REQUIRE(twoCubanes);
  auto rinfo = twoCubanes->getRingInfo();
  CHECK(rinfo->numRings() == 12);
}

TEST_CASE("check division by zero in setTerminalAtomCoords") {
  SECTION("degree 4") {
    auto m = R"CTAB(
  Mrv2311 11162401483D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 F -0.7971 -0.9945 1.275 0
M  V30 3 F 0.7971 0.9945 -1.275 0
M  V30 4 F -1.2069 -0.1568 -0.8768 0
M  V30 5 Cl 1.1223 1.1513 0.8154 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);

    CHECK_NOTHROW(MolOps::setTerminalAtomCoords(*m, 4, 0));
  }
  SECTION("degree 2, aligned 2nd neighbors") {
    // This looks like a weird mol, but it's an intermediate
    // state in AddHs.
    auto mb = R"CTAB(
  Mrv1908 06032010402D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -19.8317 16.5 0 1 CHG=-1
M  V30 2 N -18.2917 16.5 0 2 CHG=1
M  V30 3 N -16.7517 16.5 0 3
M  V30 4 H 0 0 0 3
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 3 2 3
M  V30 3 1 1 4
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB";

    v2::FileParsers::MolFileParserParams p;
    p.removeHs = false;

    auto m = v2::FileParsers::MolFromMolBlock(mb, p);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 4);

    CHECK_NOTHROW(MolOps::setTerminalAtomCoords(*m, 3, 0));
  }
}

TEST_CASE("cleanuporganometallics and carbon") {
  SECTION("basics") {
    std::vector<std::pair<std::string,
                          std::vector<std::pair<unsigned int, unsigned int>>>>
        data = {
            {"C=[CH2][Fe]", {{1, 2}}},  // this is a silly example, but it's the
                                        // simplest one I could think of
            {"[CH2]1=[CH2]2.[Fe]12", {{1, 2}, {0, 2}}},
            {"[CH2]1=[CH]2-[CH2-]3.[Fe]123", {{1, 3}, {0, 3}, {2, 3}}},
            {"[CH]12=[CH]3[CH]4=[CH]5[CH-]16.[Fe]23456",
             {{0, 5}, {1, 5}, {2, 5}, {3, 5}, {4, 5}}},
            {"[Fe]12.[CH2]1=[CH2]2",
             {{1, 0}, {2, 0}}},  // reverse the original atom order
            {"[CH]12=[CH]3[CH]4=[CH]5[CH]6=[CH]17.[Fe]234567",
             {{0, 6}, {1, 6}, {2, 6}, {3, 6}, {4, 6}, {5, 6}}},
            {"[cH]12[cH]3[cH]4[cH]5[cH]6[cH]17.[Fe]234567",
             {{0, 6}, {1, 6}, {2, 6}, {3, 6}, {4, 6}, {5, 6}}},
            {"[cH]12[cH]3[cH]4[cH]5[nH]16.[Fe]23456",
             {{0, 5}, {1, 5}, {2, 5}, {3, 5}, {4, 5}}},
            {"[CH]12=[CH]3[CH]4=[CH]5[NH]16.[Fe]23456",
             {{0, 5}, {1, 5}, {2, 5}, {3, 5}, {4, 5}}},
        };
    for (const auto &pr : data) {
      INFO(pr.first);
      auto mol = v2::SmilesParse::MolFromSmiles(pr.first);
      REQUIRE(mol);
      for (const auto &pair : pr.second) {
        REQUIRE(mol->getBondBetweenAtoms(pair.first, pair.second) != nullptr);
        CHECK(mol->getBondBetweenAtoms(pair.first, pair.second)
                  ->getBeginAtomIdx() == pair.first);
        CHECK(
            mol->getBondBetweenAtoms(pair.first, pair.second)->getBondType() ==
            Bond::BondType::DATIVE);
      }
    }
  }
  SECTION("no dative bonds") {
    std::vector<std::string> smileses = {
        "C=[CH][Fe]",
        "[CH]1=[CH]2.[Fe]12",
        "[CH]1=[C]2-[CH-]3.[Fe]123",
        "[C]12=[C]3[C]4=[C]5[C-]16.[Fe]23456",
    };
    for (const auto &smiles : smileses) {
      INFO(smiles);
      auto mol = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(mol);
      for (const auto bond : mol->bonds()) {
        CHECK(bond->getBondType() != Bond::BondType::DATIVE);
      }
    }
  }
  SECTION("github #8312") {
    std::string mb = R"CTAB(
  ChemDraw03012503262D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 21 24 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.822529 0.405221 0.000000 0 VAL=5
M  V30 2 C -1.822529 -0.419779 0.000000 0 VAL=5
M  V30 3 C -1.239167 -1.003142 0.000000 0
M  V30 4 C -0.414167 -1.003142 0.000000 0
M  V30 5 C 0.169196 -0.419779 0.000000 0 VAL=5
M  V30 6 C 0.169196 0.405221 0.000000 0 VAL=5
M  V30 7 C -0.414167 0.988584 0.000000 0
M  V30 8 C -1.239167 0.988584 0.000000 0
M  V30 9 Pt 1.335541 -0.028481 0.000000 0 VAL=6
M  V30 10 C 1.697726 0.712766 0.000000 0
M  V30 11 C 1.796387 -0.712766 0.000000 0
M  V30 12 H 2.520757 0.769729 0.000000 0
M  V30 13 H 1.236879 1.397051 0.000000 0
M  V30 14 H 2.059910 1.454013 0.000000 0
M  V30 15 H 2.619419 -0.655803 0.000000 0
M  V30 16 H 1.434203 -1.454013 0.000000 0
M  V30 17 H 2.257234 -1.397051 0.000000 0
M  V30 18 H -2.619419 0.618747 0.000000 0
M  V30 19 H -2.619419 -0.633305 0.000000 0
M  V30 20 H 0.382722 1.202110 0.000000 0
M  V30 21 H 0.581696 -1.134250 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 1 7 8
M  V30 8 1 8 1
M  V30 9 1 9 10
M  V30 10 1 9 11
M  V30 11 1 10 12
M  V30 12 1 10 13
M  V30 13 1 10 14
M  V30 14 1 11 15
M  V30 15 1 11 16
M  V30 16 1 11 17
M  V30 17 1 9 6
M  V30 18 1 9 5
M  V30 19 1 9 1
M  V30 20 1 9 2
M  V30 21 1 1 18
M  V30 22 1 2 19
M  V30 23 1 6 20
M  V30 24 1 5 21
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB";
    auto mol = v2::FileParsers::MolFromMolBlock(mb);
    REQUIRE(mol);
    REQUIRE(mol->getBondBetweenAtoms(5, 8));
    CHECK(mol->getBondBetweenAtoms(5, 8)->getBondType() ==
          Bond::BondType::DATIVE);
    REQUIRE(mol->getBondBetweenAtoms(4, 8));
    CHECK(mol->getBondBetweenAtoms(4, 8)->getBondType() ==
          Bond::BondType::DATIVE);
    REQUIRE(mol->getBondBetweenAtoms(0, 8));
    CHECK(mol->getBondBetweenAtoms(0, 8)->getBondType() ==
          Bond::BondType::DATIVE);
    REQUIRE(mol->getBondBetweenAtoms(1, 8));
    CHECK(mol->getBondBetweenAtoms(1, 8)->getBondType() ==
          Bond::BondType::DATIVE);
    REQUIRE(mol->getBondBetweenAtoms(9, 8));
    CHECK(mol->getBondBetweenAtoms(9, 8)->getBondType() ==
          Bond::BondType::SINGLE);
    REQUIRE(mol->getBondBetweenAtoms(10, 8));
    CHECK(mol->getBondBetweenAtoms(10, 8)->getBondType() ==
          Bond::BondType::SINGLE);
  }
}