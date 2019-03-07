//
//
//  Copyright (C) 2018-2019 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;

TEST_CASE("Github #1972", "[SMILES,bug]") {
  SECTION("basics") {
    std::vector<std::vector<std::string>> smiles = {
        {"[C@@]1(Cl)(F)(I).Br1", "[C@@](Br)(Cl)(F)(I)"},
        {"[C@@](Cl)(F)(I)1.Br1", "[C@@](Cl)(F)(I)Br"},
        {"[C@@](Cl)1(F)(I).Br1", "[C@@](Cl)(Br)(F)(I)"},
        {"[C@@](Cl)(F)1(I).Br1", "[C@@](Cl)(F)(Br)(I)"}};
    for (const auto &pr : smiles) {
      std::unique_ptr<ROMol> m1(SmilesToMol(pr[0]));
      std::unique_ptr<ROMol> m2(SmilesToMol(pr[1]));
      REQUIRE(m1);
      REQUIRE(m2);
      auto csmi1 = MolToSmiles(*m1);
      auto csmi2 = MolToSmiles(*m2);
      CHECK(csmi1 == csmi2);
    }
  }
  SECTION("further examples") {
    std::vector<std::vector<std::string>> smiles = {
        {"[C@@]1(Cl)2(I).Br1.F2", "[C@@](Br)(Cl)(F)(I)"},
        {"[C@@](Cl)2(I)1.Br1.F2", "[C@@](Cl)(F)(I)Br"},
        {"[C@@]12(Cl)(I).Br1.F2", "[C@@](Br)(F)(Cl)(I)"},
        {"[C@@]21(Cl)(I).Br1.F2", "[C@@](F)(Br)(Cl)(I)"},
        {"[C@@](Cl)12(I).Br1.F2", "[C@@](Cl)(Br)(F)(I)"},
        {"[C@@](Cl)21(I).Br1.F2", "[C@@](Cl)(F)(Br)(I)"},
        {"[C@@](Cl)(I)21.Br1.F2", "[C@@](Cl)(I)(F)(Br)"},
        {"[C@@](Cl)(I)12.Br1.F2", "[C@@](Cl)(I)(Br)(F)"}};
    for (const auto &pr : smiles) {
      std::unique_ptr<ROMol> m1(SmilesToMol(pr[0]));
      std::unique_ptr<ROMol> m2(SmilesToMol(pr[1]));
      REQUIRE(m1);
      REQUIRE(m2);
      auto csmi1 = MolToSmiles(*m1);
      auto csmi2 = MolToSmiles(*m2);
      CHECK(csmi1 == csmi2);
    }
  }
}

TEST_CASE("Github #2029", "[SMILES,bug]") {
  SECTION("wedging") {
    std::unique_ptr<ROMol> m1(SmilesToMol("CN[C@H](Cl)C(=O)O"));
    REQUIRE(m1);
    m1->getBondWithIdx(1)->setBondDir(Bond::BEGINWEDGE);
    bool doKekule = false, allBondsExplicit = false;
    CHECK("" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(1), -1, doKekule,
                                           allBondsExplicit));
    allBondsExplicit = true;
    CHECK("-" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(1), -1, doKekule,
                                            allBondsExplicit));
  }
  SECTION("direction") {
    std::unique_ptr<ROMol> m1(SmilesToMol("C/C=C/C"));
    REQUIRE(m1);
    bool doKekule = false, allBondsExplicit = false;
    CHECK("" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(0), -1, doKekule,
                                           allBondsExplicit));
    CHECK("" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(2), -1, doKekule,
                                           allBondsExplicit));
    allBondsExplicit = true;
    CHECK("/" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(0), -1, doKekule,
                                            allBondsExplicit));
    CHECK("/" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(2), -1, doKekule,
                                            allBondsExplicit));
  }
  SECTION("aromatic double bonds") {
    std::unique_ptr<RWMol> m1(SmilesToMol("c1ccccc1"));
    REQUIRE(m1);
    bool markAtomsBonds = false;
    MolOps::Kekulize(*m1, markAtomsBonds);
    bool doKekule = false, allBondsExplicit = false;
    CHECK("" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(0), -1, doKekule,
                                           allBondsExplicit));
    CHECK("" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(1), -1, doKekule,
                                           allBondsExplicit));
    allBondsExplicit = true;
    CHECK("=" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(0), -1, doKekule,
                                            allBondsExplicit));
    CHECK("-" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(1), -1, doKekule,
                                            allBondsExplicit));
  }
}

TEST_CASE("Smiles literals", "[SMILES]") {
  auto mol = "c1ccccc1"_smiles;
  REQUIRE(mol);
  CHECK(6 == mol->getNumAtoms());
  auto fail1 = "c1ccccc"_smiles;
  REQUIRE(!fail1);
  auto fail2 = "c1cccn1"_smiles;
  REQUIRE(!fail2);
}

TEST_CASE("Smarts literals", "[Smarts]") {
  auto mol = "c1ccc[c,n]c1"_smarts;
  REQUIRE(mol);
  CHECK(6 == mol->getNumAtoms());
  auto fail1 = "c1ccccc"_smarts;
  REQUIRE(!fail1);
  auto mol2 = "c1cccn1"_smarts;
  REQUIRE(mol2);
}

TEST_CASE(
    "github #2197 and #2237: handling of aromatic main group atoms in SMARTS",
    "[Smarts]") {
  std::vector<std::string> smarts = {
      "[si]1ccccc1",
      "[as]1ccccc1",
      "[se]1ccccc1",
      "[te]1ccccc1",

  };
  SECTION("#2197") {
    for (const auto sma : smarts) {
      std::unique_ptr<ROMol> mol(SmartsToMol(sma));
      REQUIRE(mol);
      CHECK(6 == mol->getNumAtoms());
      REQUIRE(mol->getAtomWithIdx(0)->hasQuery());
      REQUIRE(static_cast<QueryAtom *>(mol->getAtomWithIdx(0))
                  ->getQuery()
                  ->getDescription() == "AtomType");
    }
  }
  SECTION("#2237") {
    for (const auto sma : smarts) {
      std::unique_ptr<ROMol> mol(SmartsToMol(sma));
      REQUIRE(mol);
      REQUIRE(MolToSmarts(*mol) == sma);
    }
  }
}

TEST_CASE("github #2257: writing cxsmiles", "[smiles,cxsmiles]") {
  SECTION("basics") {
    auto mol = "OCC"_smiles;
    REQUIRE(mol);
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "CCO");
  }
  SECTION("atom labels") {
    auto mol = "CCC |$R1;;R2$|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(0)->getProp<std::string>(
              common_properties::atomLabel) == "R1");
    CHECK(mol->getAtomWithIdx(2)->getProp<std::string>(
              common_properties::atomLabel) == "R2");
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "CCC |$R1;;R2$|");
  }
  SECTION("atom ordering") {
    auto mol = "OC(F)C |$R1;;R2;R3$|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(0)->getProp<std::string>(
              common_properties::atomLabel) == "R1");
    CHECK(mol->getAtomWithIdx(2)->getProp<std::string>(
              common_properties::atomLabel) == "R2");
    CHECK(mol->getAtomWithIdx(3)->getProp<std::string>(
              common_properties::atomLabel) == "R3");
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "CC(O)F |$R3;;R1;R2$|");
  }
  SECTION("atom values") {
    auto mol = "COCC |$_AV:;bar;;foo$|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(3)->getProp<std::string>(
              common_properties::molFileValue) == "foo");
    CHECK(mol->getAtomWithIdx(1)->getProp<std::string>(
              common_properties::molFileValue) == "bar");
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "CCOC |$_AV:foo;;bar;$|");
  }

  SECTION("radicals") {
    auto mol = "[Fe]N([O])[O] |^1:2,3|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(1)->getNumRadicalElectrons() == 0);
    CHECK(mol->getAtomWithIdx(2)->getNumRadicalElectrons() == 1);
    CHECK(mol->getAtomWithIdx(3)->getNumRadicalElectrons() == 1);

    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "[O]N([O])[Fe] |^1:0,2|");
  }
  SECTION("radicals2") {
    auto mol = "[CH]C[CH2] |^1:2,^2:0|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(1)->getNumRadicalElectrons() == 0);
    CHECK(mol->getAtomWithIdx(2)->getNumRadicalElectrons() == 1);
    CHECK(mol->getAtomWithIdx(0)->getNumRadicalElectrons() == 2);

    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "[CH]C[CH2] |^1:2,^2:0|");
  }
  SECTION("coordinates") {
    auto mol = "OC |(0,.75,;0,-.75,)|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getNumConformers() == 1);

    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "CO |(0,-0.75,;0,0.75,)|");
  }
  SECTION("coordinates3d") {
    auto mol = "OC |(0,.75,0.1;0,-.75,-0.1)|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getNumConformers() == 1);

    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "CO |(0,-0.75,-0.1;0,0.75,0.1)|");
  }
  SECTION("atom props") {
    auto mol = "N1CC1C |atomProp:0.p2.v2:0.p1.v1:1.p2.v2:1.p1.v1;2;3|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 4);
    CHECK(mol->getAtomWithIdx(0)->hasProp("p1"));
    CHECK(mol->getAtomWithIdx(0)->getProp<std::string>("p1") == "v1");
    CHECK(mol->getAtomWithIdx(0)->hasProp("p2"));
    CHECK(mol->getAtomWithIdx(0)->getProp<std::string>("p2") == "v2");
    CHECK(mol->getAtomWithIdx(1)->hasProp("p2"));
    CHECK(mol->getAtomWithIdx(1)->getProp<std::string>("p2") == "v2");
    CHECK(mol->getAtomWithIdx(1)->hasProp("p1"));
    CHECK(mol->getAtomWithIdx(1)->getProp<std::string>("p1") == "v1;2;3");

    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "CC1CN1 |atomProp:2.p2.v2:2.p1.v1;2;3:3.p2.v2:3.p1.v1|");
  }
  SECTION("atom props and values") {
    //"CN |$_AV:atomv0;atomv1$,atomProp:0.p2.v2:1.p2.v1|";
    auto mol = "CN |atomProp:0.p2.v2:1.p1.v1,$_AV:val1;val2$|"_smiles;
    REQUIRE(mol);
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "CN |$_AV:val1;val2$,atomProp:0.p2.v2:1.p1.v1|");
  }
  SECTION("enhanced stereo 1") {
    auto mol = "C[C@H](F)[C@H](C)[C@@H](C)Br |a:1,o1:4,5|"_smiles;
    REQUIRE(mol);
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "C[C@H](F)[C@H](C)[C@@H](C)Br |a:1,o1:4,5|");
  }

  SECTION("enhanced stereo 2") {
    auto mol = "C[C@H](O)[C@H](CC)F |o1:1,3|"_smiles;
    REQUIRE(mol);
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "CC[C@H](F)[C@H](C)O |o1:2,4|");
  }

  SECTION("enhanced stereo 3") {
    auto mol = "C[C@@H]1N[C@H](C)[C@@H]([C@H](C)[C@@H]1C)C1[C@@H](C)O[C@@H](C)[C@@H](C)[C@H]1C |a:5,o1:1,8,o2:14,16,&1:11,18,&2:3,6,r|"_smiles;
    REQUIRE(mol);
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "C[C@@H]1N[C@H](C)[C@H](C2[C@@H](C)O[C@@H](C)[C@@H](C)[C@H]2C)[C@H](C)[C@@H]1C |a:5,o1:1,18,o2:10,12,&1:3,16,&2:7,14|");
  }

  SECTION("enhanced stereo 4") {
    auto mol = "C[C@@H]1CCO[C@H](C)C1 |a:1,5,r|"_smiles;
    REQUIRE(mol);
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "C[C@@H]1CCO[C@H](C)C1 |a:1,5|");
  }

  SECTION("mol fragments1") {
    auto mol = "Cl.OC |(1,0,0;0,.75,0.1;0,-.75,-0.1)|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getNumConformers() == 1);

    std::vector<int> atomsToUse = {1, 2};
    auto smi = MolFragmentToCXSmiles(*mol, atomsToUse);
    CHECK(smi == "CO |(0,-0.75,-0.1;0,0.75,0.1)|");
  }
  SECTION("mol fragments2") {
    auto mol = "Cl.N1CC1C |atomProp:1.p2.v1:1.p1.v1:2.p2.v2:2.p1.v2|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 5);
    CHECK(!mol->getAtomWithIdx(0)->hasProp("p1"));
    CHECK(mol->getAtomWithIdx(1)->hasProp("p1"));
    CHECK(mol->getAtomWithIdx(1)->getProp<std::string>("p1") == "v1");

    std::vector<int> atomsToUse = {1, 2, 3, 4};
    auto smi = MolFragmentToCXSmiles(*mol, atomsToUse);
    CHECK(smi == "CC1CN1 |atomProp:2.p2.v2:2.p1.v2:3.p2.v1:3.p1.v1|");
  }

  SECTION("mol fragments3") {
    auto mol = "Cl.[CH]C[CH2] |^1:3,^2:1|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(2)->getNumRadicalElectrons() == 0);
    CHECK(mol->getAtomWithIdx(3)->getNumRadicalElectrons() == 1);
    CHECK(mol->getAtomWithIdx(1)->getNumRadicalElectrons() == 2);

    std::vector<int> atomsToUse = {1, 2, 3};
    auto smi = MolFragmentToCXSmiles(*mol, atomsToUse);
    CHECK(smi == "[CH]C[CH2] |^1:2,^2:0|");
  }
}

TEST_CASE("Github #2148", "[bug, Smiles, Smarts]") {
  SECTION("SMILES") {
    auto mol = "C(=C\\F)\\4.O=C1C=4CCc2ccccc21"_smiles;
    REQUIRE(mol);
    REQUIRE(mol->getBondBetweenAtoms(0, 5));
    CHECK(mol->getBondBetweenAtoms(0, 5)->getBondType() == Bond::DOUBLE);
  }
  SECTION("SMILES edges") {
    auto m1 = "C/C=C/C"_smiles;
    REQUIRE(m1);
    CHECK(m1->getBondBetweenAtoms(2, 1)->getBondType() == Bond::DOUBLE);
    CHECK(m1->getBondBetweenAtoms(2, 1)->getStereo() != Bond::STEREONONE);

    {
      std::vector<std::string> smis = {"C1=C/C.C/1", "C/1=C/C.C1",
                                       "C-1=C/C.C/1", "C/1=C/C.C-1"};
      for (auto smi : smis) {
        std::unique_ptr<RWMol> mol(SmilesToMol(smi));
        REQUIRE(mol);
        CHECK(mol->getBondBetweenAtoms(0, 3)->getBondType() == Bond::SINGLE);
        CHECK(mol->getBondBetweenAtoms(0, 3)->getBondDir() != Bond::NONE);
        CHECK(mol->getBondBetweenAtoms(0, 1)->getBondType() == Bond::DOUBLE);
        CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() != Bond::STEREONONE);
      }
    }
  }

  SECTION("Writing SMILES") {
    auto mol = "C/C=c1/ncc(=C)cc1"_smiles;
    REQUIRE(mol);
    REQUIRE(mol->getBondBetweenAtoms(1, 2));
    CHECK(mol->getBondBetweenAtoms(1, 2)->getBondType() == Bond::DOUBLE);
    CHECK(mol->getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOE);
    auto smi = MolToSmiles(*mol);
    CHECK(smi == "C=c1cc/c(=C\\C)nc1");
  }
}

TEST_CASE("Github #2298", "[bug, Smarts, substructure]") {
  SubstructMatchParameters ps;
  ps.useQueryQueryMatches = true;
  SECTION("basics") {
    auto m1 = "[#6]"_smarts;
    REQUIRE(m1);
    CHECK(SubstructMatch(*m1, *m1, ps).size() == 1);
    auto m2 = "[C]"_smarts;
    REQUIRE(m2);
    CHECK(SubstructMatch(*m2, *m2, ps).size() == 1);
    auto m3 = "[C]"_smarts;
    REQUIRE(m3);
    CHECK(SubstructMatch(*m3, *m3, ps).size() == 1);
  }
  SECTION("a bit more complex") {
    auto m1 = "[CH0+2]"_smarts;
    REQUIRE(m1);
    CHECK(SubstructMatch(*m1, *m1, ps).size() == 1);
  }
}

TEST_CASE("dative ring closures", "[bug, smiles]") {
  SECTION("first closure1") {
    auto m1 = "N->1CCN->[Pt]1"_smiles;
    REQUIRE(m1);
    REQUIRE(m1->getBondBetweenAtoms(0, 4));
    CHECK(m1->getBondBetweenAtoms(0, 4)->getBondType() == Bond::DATIVE);
    CHECK(m1->getBondBetweenAtoms(0, 4)->getBeginAtomIdx() == 0);
  }
  SECTION("first closure2") {
    auto m1 = "[Pt]<-1CCCN1"_smiles;
    REQUIRE(m1);
    REQUIRE(m1->getBondBetweenAtoms(0, 4));
    CHECK(m1->getBondBetweenAtoms(0, 4)->getBondType() == Bond::DATIVE);
    CHECK(m1->getBondBetweenAtoms(0, 4)->getBeginAtomIdx() == 4);
  }
  SECTION("second closure1") {
    auto m1 = "N1CCN->[Pt]<-1"_smiles;
    REQUIRE(m1);
    REQUIRE(m1->getBondBetweenAtoms(0, 4));
    CHECK(m1->getBondBetweenAtoms(0, 4)->getBondType() == Bond::DATIVE);
    CHECK(m1->getBondBetweenAtoms(0, 4)->getBeginAtomIdx() == 0);
  }
  SECTION("second closure2") {
    auto m1 = "[Pt]1CCCN->1"_smiles;
    REQUIRE(m1);
    REQUIRE(m1->getBondBetweenAtoms(0, 4));
    CHECK(m1->getBondBetweenAtoms(0, 4)->getBondType() == Bond::DATIVE);
    CHECK(m1->getBondBetweenAtoms(0, 4)->getBeginAtomIdx() == 4);
  }
  SECTION("branch1") {
    auto m1 = "N(->[Pt])C"_smiles;
    REQUIRE(m1);
    REQUIRE(m1->getBondBetweenAtoms(0, 1));
    CHECK(m1->getBondBetweenAtoms(0, 1)->getBondType() == Bond::DATIVE);
    CHECK(m1->getBondBetweenAtoms(0, 1)->getBeginAtomIdx() == 0);
  }
  SECTION("branch2") {
    auto m1 = "N(->[Pt])C"_smiles;
    REQUIRE(m1);
    REQUIRE(m1->getBondBetweenAtoms(0, 1));
    CHECK(m1->getBondBetweenAtoms(0, 1)->getBondType() == Bond::DATIVE);
    CHECK(m1->getBondBetweenAtoms(0, 1)->getBeginAtomIdx() == 0);
  }
}