//
//  Copyright (C) 2018-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>
#ifdef RDK_BUILD_THREADSAFE_SSS
#include <future>
#include <thread>
#endif

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/test_fixtures.h>
#include <GraphMol/Canon.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/SmilesParse/CanonicalizeStereoGroups.h>
#include <GraphMol/MolOps.h>

using namespace RDKit;

TEST_CASE("Github #1972", "[SMILES][bug]") {
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

TEST_CASE("Github #2029", "[SMILES][bug]") {
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
    for (const auto &sma : smarts) {
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
    for (const auto &sma : smarts) {
      std::unique_ptr<ROMol> mol(SmartsToMol(sma));
      REQUIRE(mol);
      REQUIRE(MolToSmarts(*mol) == sma);
    }
  }
}

TEST_CASE("github #2257: writing cxsmiles", "[smiles][cxsmiles]") {
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
    CHECK(smi == "[O][N]([O])[Fe] |^1:0,2|");
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
    auto mol = "C[C@H](F)[C@H](C)[C@@H](C)Br |a:1,o1:3,5|"_smiles;
    REQUIRE(mol);
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "C[C@H](F)[C@@H](C)[C@H](C)Br |a:1,o1:3,5|");
  }

  SECTION("enhanced stereo 2") {
    auto mol = "C[C@H](O)[C@H](CC)F |o1:1,3|"_smiles;
    REQUIRE(mol);
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "CC[C@H](F)[C@H](C)O |o1:2,4|");
  }

  SECTION("enhanced stereo 3") {
    auto mol =
        "C[C@@H]1N[C@H](C)[C@@H]([C@H](C)[C@@H]1C)C1[C@@H](C)O[C@@H](C)[C@@H](C)[C@H]1C |a:5,o1:1,8,o2:14,16,&1:11,18,&2:3,6,r|"_smiles;
    REQUIRE(mol);
    auto smi = MolToCXSmiles(*mol);
    CHECK(
        smi ==
        "C[C@H]1N[C@H](C)[C@H](C2[C@H](C)O[C@H](C)[C@H](C)[C@@H]2C)[C@H](C)[C@H]1C |a:5,o1:1,18,o2:10,12,&1:3,16,&2:7,14|");
  }

  SECTION("enhanced stereo 4") {
    auto mol = "C[C@@H]1CCO[C@H](C)C1 |a:1,5,r|"_smiles;
    REQUIRE(mol);
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "C[C@@H]1CCO[C@H](C)C1 |a:1,5|");
  }

  SECTION("enhanced stereo with other properties") {
    auto mol = "CC[C@H](C)O |atomProp:3.p2.v2,o1:2|"_smiles;
    REQUIRE(mol);
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "CC[C@H](C)O |atomProp:3.p2.v2,o1:2|");
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

TEST_CASE("Github #2148", "[bug][Smiles][Smarts]") {
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
    const auto useLegacy = GENERATE(true, false);
    CAPTURE(useLegacy);
    UseLegacyStereoPerceptionFixture fx(useLegacy);
    auto mol = "C/C=c1/ncc(=C)cc1"_smiles;
    REQUIRE(mol);
    REQUIRE(mol->getBondBetweenAtoms(1, 2));
    CHECK(mol->getBondBetweenAtoms(1, 2)->getBondType() == Bond::DOUBLE);
    if (useLegacy) {
      CHECK(mol->getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOE);
    } else {
      CHECK(mol->getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOTRANS);
      CHECK(mol->getBondBetweenAtoms(1, 2)->getStereoAtoms() ==
            std::vector<int>{0, 3});
    }
    auto smi = MolToSmiles(*mol);
    CHECK(smi == "C=c1cc/c(=C\\C)nc1");
  }
}

TEST_CASE("Github #2298", "[bug][Smarts][substructure]") {
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

TEST_CASE("dative ring closures", "[bug][smiles]") {
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

TEST_CASE("github#2450: getAtomSmarts() fails for free atoms", "[bug]") {
  SECTION("original report") {
    std::unique_ptr<QueryAtom> qat(new QueryAtom());
    qat->setQuery(makeAtomNumQuery(6));
    auto smarts = SmartsWrite::GetAtomSmarts(qat.get());
    CHECK(smarts == "[#6]");
  }
  SECTION("query bonds") {
    std::unique_ptr<QueryBond> qbnd(new QueryBond(Bond::AROMATIC));
    auto smarts = SmartsWrite::GetBondSmarts(qbnd.get());
    CHECK(smarts == ":");
  }
  SECTION("SMILES works too") {
    std::unique_ptr<Bond> bnd(new Bond(Bond::AROMATIC));
    auto smiles = SmilesWrite::GetBondSmiles(bnd.get());
    CHECK(smiles == ":");
  }
}

TEST_CASE("MolFragmentToSmarts", "[Smarts]") {
  SECTION("BasicFragment") {
    auto m = "CCCCCN"_smiles;
    std::vector<int> indices = {3, 4, 5};
    const auto smarts = MolFragmentToSmarts(*m, indices);
    CHECK(smarts == "[#6]-[#6]-[#7]");
  }
  SECTION("FragmentWithParity1") {
    auto m = "C[C@H](F)CCCN"_smiles;
    std::vector<int> indices = {0, 1, 2, 3};
    const auto smarts = MolFragmentToSmarts(*m, indices);
    CHECK(smarts == "[#6]-[#6@H](-[#9])-[#6]");
  }
  SECTION("FragmentWithParity2") {
    auto m = "C[C@](F)(Cl)CCCN"_smiles;
    std::vector<int> indices = {0, 1, 2, 4};
    const auto smarts = MolFragmentToSmarts(*m, indices);
    CHECK(smarts == "[#6]-[#6@@](-[#9])-[#6]");
  }
  SECTION("FragmentLosingParity") {
    auto m = "C[C@H](F)CCCN"_smiles;
    std::vector<int> indices = {0, 1, 2};
    const auto smarts = MolFragmentToSmarts(*m, indices);
    CHECK(smarts == "[#6]-[#6@H]-[#9]");
  }
  SECTION("FragmentWithSpecifiedBonds") {
    auto m = "C1CC1O"_smiles;
    std::vector<int> atomIndices = {0, 1, 2};
    std::vector<int> bondIndices = {0};
    const auto smarts = MolFragmentToSmarts(*m, atomIndices, &bondIndices);
    CHECK(smarts == "[#6]-[#6].[#6]");
  }
  SECTION("SmartsFragmentFromQueryMol") {
    auto m = "CCCC[C,N]N"_smarts;
    std::vector<int> indices = {3, 4, 5};
    const auto smarts = MolFragmentToSmarts(*m, indices);
    CHECK(smarts == "C[C,N]N");
  }
}

TEST_CASE("github #2667: MolToCXSmiles generates error for empty molecule",
          "[bug][cxsmiles]") {
  SECTION("basics") {
    auto mol = ""_smiles;
    REQUIRE(mol);
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "");
  }
}

TEST_CASE("github #2604: support range-based charge queries from SMARTS",
          "[ranges][smarts]") {
  SECTION("positive") {
    auto query = "[N+{0-1}]"_smarts;
    REQUIRE(query);
    {
      auto m1 = "CN"_smiles;
      REQUIRE(m1);
      CHECK(SubstructMatch(*m1, *query).size() == 1);
    }
    {
      auto m1 = "C[NH3+]"_smiles;
      REQUIRE(m1);
      CHECK(SubstructMatch(*m1, *query).size() == 1);
    }
    {
      auto m1 = "C[NH2+2]"_smiles;
      REQUIRE(m1);
      CHECK(SubstructMatch(*m1, *query).empty());
    }
    {
      auto m1 = "C[NH-]"_smiles;
      REQUIRE(m1);
      CHECK(SubstructMatch(*m1, *query).empty());
    }
  }
  SECTION("negative") {
    auto query = "[N-{0-1}]"_smarts;
    REQUIRE(query);
    {
      auto m1 = "CN"_smiles;
      REQUIRE(m1);
      CHECK(SubstructMatch(*m1, *query).size() == 1);
    }
    {
      auto m1 = "C[NH-]"_smiles;
      REQUIRE(m1);
      CHECK(SubstructMatch(*m1, *query).size() == 1);
    }
    {
      auto m1 = "C[N-2]"_smiles;
      REQUIRE(m1);
      CHECK(SubstructMatch(*m1, *query).empty());
    }
    {
      auto m1 = "C[NH3+]"_smiles;
      REQUIRE(m1);
      CHECK(SubstructMatch(*m1, *query).empty());
    }
  }
}

TEST_CASE("_smarts fails gracefully", "[smarts]") {
  SECTION("empty") {
    auto mol = ""_smarts;
    REQUIRE(mol);
  }
  SECTION("syntax error") {
    auto mol = "C1C"_smarts;
    REQUIRE(!mol);
  }
}

TEST_CASE(
    "github #2801: MolToSmarts may generate invalid SMARTS for bond queries",
    "[bug][smarts]") {
  SECTION("original_report") {
    auto q1 = "*~CCC"_smarts;
    REQUIRE(q1);
    Bond *qb = q1->getBondBetweenAtoms(0, 1);
    BOND_EQUALS_QUERY *bq1 = makeBondOrderEqualsQuery(qb->getBondType());
    qb->setQuery(bq1);
    BOND_EQUALS_QUERY *bq2 = makeBondIsInRingQuery();
    bq2->setNegation(true);
    qb->expandQuery(bq2, Queries::COMPOSITE_AND, true);
    std::string smarts = MolToSmarts(*q1);
    CHECK(smarts == "*!@CCC");
    std::unique_ptr<RWMol> q2(SmartsToMol(smarts));
    REQUIRE(q2);
  }
  SECTION("composite_or") {
    auto q1 = "*~CCC"_smarts;
    REQUIRE(q1);
    Bond *qb = q1->getBondBetweenAtoms(0, 1);
    BOND_EQUALS_QUERY *bq1 = makeBondOrderEqualsQuery(qb->getBondType());
    qb->setQuery(bq1);
    BOND_EQUALS_QUERY *bq2 = makeBondIsInRingQuery();
    bq2->setNegation(true);
    qb->expandQuery(bq2, Queries::COMPOSITE_OR, true);
    // this used to yield *,!@CCC
    std::string smarts = MolToSmarts(*q1);
    CHECK(smarts == "*!@CCC");
    std::unique_ptr<RWMol> q2(SmartsToMol(smarts));
    REQUIRE(q2);
  }
  SECTION("composite_lowand") {
    auto q1 = "*~CCC"_smarts;
    REQUIRE(q1);
    Bond *qb = q1->getBondBetweenAtoms(0, 1);
    BOND_EQUALS_QUERY *bq1 = makeBondOrderEqualsQuery(qb->getBondType());
    qb->setQuery(bq1);
    BOND_EQUALS_QUERY *bq2 = makeBondOrderEqualsQuery(qb->getBondType());
    qb->expandQuery(bq2, Queries::COMPOSITE_OR, true);
    BOND_EQUALS_QUERY *bq3 = makeBondIsInRingQuery();
    bq3->setNegation(true);
    qb->expandQuery(bq3, Queries::COMPOSITE_AND, true);
    std::string smarts = MolToSmarts(*q1);
    CHECK(smarts == "*!@CCC");
    std::unique_ptr<RWMol> q2(SmartsToMol(smarts));
    REQUIRE(q2);
  }
}

TEST_CASE("large rings", "[smarts]") {
  auto query = "[r24]"_smarts;
  auto m_r24 = "C1CCCCCCCCCCCCCCCCCCCCCCC1"_smiles;
  auto m_r23 = "C1CCCCCCCCCCCCCCCCCCCCCC1"_smiles;

  CHECK(SubstructMatch(*m_r23, *query).empty());
  CHECK(SubstructMatch(*m_r24, *query).size() == 24);
}

TEST_CASE("random smiles vectors", "[smiles]") {
  auto m = "C1OCC1N(CO)(Cc1ccccc1NCCl)"_smiles;
  REQUIRE(m);
  SECTION("basics") {
    std::vector<std::string> tgt = {
        "c1cc(CN(C2COC2)CO)c(cc1)NCCl", "N(CCl)c1c(CN(C2COC2)CO)cccc1",
        "N(CCl)c1ccccc1CN(C1COC1)CO", "OCN(Cc1ccccc1NCCl)C1COC1",
        "C(N(C1COC1)Cc1c(cccc1)NCCl)O"};
    unsigned int randomSeed = 0xf00d;
    auto smiV = MolToRandomSmilesVect(*m, 5, randomSeed);
    CHECK(smiV == tgt);
  }
  SECTION("options1") {
    std::vector<std::string> tgt = {
        "C1-C=C(-C-N(-C2-C-O-C-2)-C-O)-C(=C-C=1)-N-C-Cl",
        "N(-C-Cl)-C1-C(-C-N(-C2-C-O-C-2)-C-O)=C-C=C-C=1",
        "N(-C-Cl)-C1=C-C=C-C=C-1-C-N(-C1-C-O-C-1)-C-O",
        "O-C-N(-C-C1=C-C=C-C=C-1-N-C-Cl)-C1-C-O-C-1",
        "C(-N(-C1-C-O-C-1)-C-C1-C(=C-C=C-C=1)-N-C-Cl)-O"};
    RWMol nm(*m);
    MolOps::Kekulize(nm, true);
    unsigned int randomSeed = 0xf00d;
    bool isomericSmiles = true;
    bool kekuleSmiles = true;
    bool allBondsExplicit = true;
    bool allHsExplicit = false;
    auto smiV =
        MolToRandomSmilesVect(nm, 5, randomSeed, isomericSmiles, kekuleSmiles,
                              allBondsExplicit, allHsExplicit);
    CHECK(smiV == tgt);
  }
  SECTION("options2") {
    std::vector<std::string> tgt = {
        "[cH]1[cH][c]([CH2][N]([CH]2[CH2][O][CH2]2)[CH2][OH])[c]([cH][cH]1)[NH]"
        "[CH2][Cl]",
        "[NH]([CH2][Cl])[c]1[c]([CH2][N]([CH]2[CH2][O][CH2]2)[CH2][OH])[cH][cH]"
        "[cH][cH]1",
        "[NH]([CH2][Cl])[c]1[cH][cH][cH][cH][c]1[CH2][N]([CH]1[CH2][O][CH2]1)["
        "CH2][OH]",
        "[OH][CH2][N]([CH2][c]1[cH][cH][cH][cH][c]1[NH][CH2][Cl])[CH]1[CH2][O]["
        "CH2]1",
        "[CH2]([N]([CH]1[CH2][O][CH2]1)[CH2][c]1[c]([cH][cH][cH][cH]1)[NH][CH2]"
        "[Cl])[OH]"};
    RWMol nm(*m);
    MolOps::Kekulize(nm, false);
    unsigned int randomSeed = 0xf00d;
    bool isomericSmiles = true;
    bool kekuleSmiles = false;
    bool allBondsExplicit = false;
    bool allHsExplicit = true;
    auto smiV =
        MolToRandomSmilesVect(nm, 5, randomSeed, isomericSmiles, kekuleSmiles,
                              allBondsExplicit, allHsExplicit);
    CHECK(smiV == tgt);
  }
}

TEST_CASE(
    "github #3197: Molecule constructed from CXSMILES cannot be translated to "
    "SMARTS",
    "[smarts][bug]") {
  auto m = "C* |$;M_p$|"_smiles;
  REQUIRE(m);
  SECTION("smarts writing") {
    auto smarts = MolToSmarts(*m);
    // this will change if/when the definition of the query changes, just have
    // to update then
    CHECK(
        smarts ==
        "[#6]-[!#0&!#2&!#5&!#6&!#7&!#8&!#9&!#10&!#14&!#15&!#16&!#17&!#18&!#33&!#"
        "34&!#35&!#36&!#52&!#53&!#54&!#85&!#86&!#1]");
  }
  SECTION("serialization") {
    std::string pkl;
    MolPickler::pickleMol(*m, pkl, PicklerOps::PropertyPickleOptions::AllProps);
    ROMol cpy(pkl);
    auto osmi = MolToCXSmiles(*m);
    CHECK(osmi == "*C |$M_p;$|");
    auto smi = MolToCXSmiles(cpy);
    CHECK(smi == osmi);
    QueryAtom *oa1 = static_cast<QueryAtom *>(m->getAtomWithIdx(1));
    QueryAtom *a1 = static_cast<QueryAtom *>(m->getAtomWithIdx(1));
    REQUIRE(oa1->hasQuery());
    REQUIRE(a1->hasQuery());
    size_t osz =
        oa1->getQuery()->endChildren() - oa1->getQuery()->beginChildren();
    size_t sz = a1->getQuery()->endChildren() - a1->getQuery()->beginChildren();
    // we don't need to test the exact size (since that may change), but let's
    // at least be sure it's not unreasonable:
    CHECK(osz > 0);
    CHECK(osz < 200);
    CHECK(osz == sz);
  }
}

TEST_CASE("d primitive in SMARTS", "[smarts][extension]") {
  SmilesParserParams ps;
  ps.removeHs = false;
  std::unique_ptr<ROMol> m(SmilesToMol("[H]OCO[2H]", ps));
  REQUIRE(m);
  CHECK(m->getNumAtoms() == 5);
  SECTION("basics") {
    auto q = "[d2]"_smarts;
    REQUIRE(q);
    CHECK(SubstructMatch(*m, *q).size() == 2);
  }
  SECTION("comparison to D") {
    auto q = "[D2]"_smarts;
    REQUIRE(q);
    CHECK(SubstructMatch(*m, *q).size() == 3);
  }
}

TEST_CASE(
    "github #3342: unspecified branch bonds in SMARTS don't have aromaticity "
    "set",
    "[smarts][bug]") {
  SECTION("as reported") {
    auto m = "c1(ccccc1)"_smarts;
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(0, 1));
    CHECK(m->getBondBetweenAtoms(0, 1)->getBondType() ==
          Bond::BondType::AROMATIC);
    CHECK(m->getBondBetweenAtoms(0, 1)->getIsAromatic());
  }
}

TEST_CASE("github #3320: incorrect bond properties from CXSMILES",
          "[cxsmiles][bug]") {
  SECTION("as reported") {
    auto m = "[Cl-][Pt++]1([Cl-])NCCN1C1CCCCC1 |C:6.6,3.2,0.0,2.1|"_smiles;
    REQUIRE(m);
    std::vector<std::pair<unsigned, unsigned>> bonds = {
        {0, 1}, {3, 1}, {2, 1}, {6, 1}};
    for (const auto &pr : bonds) {
      auto bnd = m->getBondBetweenAtoms(pr.first, pr.second);
      REQUIRE(bnd);
      CHECK(bnd->getBondType() == Bond::BondType::DATIVE);
      CHECK(bnd->getBeginAtomIdx() == pr.first);
    }
  }
  SECTION("as reported") {
    auto m = "[Cl-][Pt++]1([Cl-])NCC3C2CCCCC2.N13 |C:12.12,3.2,0.0,2.1|"_smiles;
    REQUIRE(m);
    std::vector<std::pair<unsigned, unsigned>> bonds = {
        {0, 1}, {3, 1}, {2, 1}, {12, 1}};
    for (const auto &pr : bonds) {
      auto bnd = m->getBondBetweenAtoms(pr.first, pr.second);
      REQUIRE(bnd);
      CHECK(bnd->getBondType() == Bond::BondType::DATIVE);
      CHECK(bnd->getBeginAtomIdx() == pr.first);
    }
  }
}

TEST_CASE("github #3774: MolToSmarts inverts direction of dative bond",
          "[smarts][bug]") {
  SECTION("as reported") {
    {
      auto m = "N->[Cu+]"_smiles;
      REQUIRE(m);
      CHECK(MolToSmarts(*m) == "[#7]->[Cu+]");
      CHECK(MolToSmiles(*m) == "[NH3]->[Cu+]");
    }
    {
      auto m = "N<-[Cu+]"_smiles;
      REQUIRE(m);
      CHECK(MolToSmarts(*m) == "[#7]<-[Cu+]");
      CHECK(MolToSmiles(*m) == "[NH2]<-[Cu+]");
    }
  }
  SECTION("from smarts") {
    {
      auto m = "N->[Cu+]"_smarts;
      REQUIRE(m);
      CHECK(MolToSmarts(*m) == "N->[#29&+]");
    }
    {
      auto m = "N<-[Cu+]"_smarts;
      REQUIRE(m);
      CHECK(MolToSmarts(*m) == "N<-[#29&+]");
    }
  }
}

TEST_CASE("Hydrogen bonds", "[smiles]") {
  SECTION("basics") {
    auto m = "CC1O[H]O=C(C)C1 |H:4.3|"_smiles;
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(3, 4));
    CHECK(m->getBondBetweenAtoms(3, 4)->getBondType() ==
          Bond::BondType::HYDROGEN);
  }
}

TEST_CASE("Github #2788: doKekule=true should kekulize the molecule",
          "[smiles]") {
  SECTION("basics1") {
    auto m = "c1ccccc1"_smiles;
    REQUIRE(m);
    bool doIsomeric = true;
    bool doKekule = true;
    CHECK(MolToSmiles(*m, doIsomeric, doKekule) == "C1=CC=CC=C1");
  }
  SECTION("basics2") {
    auto m = "c1cc[nH]c1"_smiles;
    REQUIRE(m);
    bool doIsomeric = true;
    bool doKekule = true;
    CHECK(MolToSmiles(*m, doIsomeric, doKekule) == "C1=CNC=C1");
  }

  SECTION("can thrown exceptions") {
    int debugParse = 0;
    bool sanitize = false;
    std::unique_ptr<RWMol> m{SmilesToMol("c1ccnc1", debugParse, sanitize)};
    REQUIRE(m);
    bool doIsomeric = true;
    bool doKekule = false;
    {
      RWMol tm(*m);
      CHECK(MolToSmiles(tm, doIsomeric, doKekule) == "c1ccnc1");
    }
    doKekule = true;
    {
      RWMol tm(*m);
      CHECK_THROWS_AS(MolToSmiles(tm, doIsomeric, doKekule), KekulizeException);
    }
  }
}

TEST_CASE("bogus recursive SMARTS", "[smarts]") {
  std::string sma = "C)foo";
  CHECK(SmartsToMol(sma) == nullptr);
}

TEST_CASE(
    "Github #3998 MolFragmentToSmiles failing in Kekulization with "
    "kekuleSmiles=true") {
  auto mol = "Cc1ccccc1"_smiles;
  REQUIRE(mol);
  SECTION("normal") {
    std::vector<int> ats{0};
    std::string smi = MolFragmentToSmiles(*mol, ats);
    CHECK(smi == "C");
  }
  SECTION("kekulized") {
    std::vector<int> ats{0};
    bool doIsomericSmiles = true;
    bool doKekule = true;
    std::string smi = MolFragmentToSmiles(*mol, ats, nullptr, nullptr, nullptr,
                                          doIsomericSmiles, doKekule);
    CHECK(smi == "C");
  }
  SECTION("including ring parts") {
    std::vector<int> ats{0, 1, 2};
    bool doIsomericSmiles = true;
    bool doKekule = true;
    std::string smi = MolFragmentToSmiles(*mol, ats, nullptr, nullptr, nullptr,
                                          doIsomericSmiles, doKekule);
    CHECK(smi == "C:CC");
  }
}

TEST_CASE("Github #4319 add CXSMARTS support") {
  // note: the CXSMARTS support uses exactly the same code as the CXSMILES
  // parser/writer. We aren't testing that here since it's tested already with
  // the CXSMILES tests. The goal here is just to make sure that it's being
  // called by default and that we can control its behavior with the
  // SmartsParseParams structure
  SECTION("defaults") {
    auto mol = "CCC |$foo;;bar$|"_smarts;
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 3);
    CHECK(mol->getAtomWithIdx(0)->getProp<std::string>(
              common_properties::atomLabel) == "foo");
    CHECK(mol->getAtomWithIdx(2)->getProp<std::string>(
              common_properties::atomLabel) == "bar");
    CHECK(!mol->getAtomWithIdx(1)->hasProp(common_properties::atomLabel));
  }
  SECTION("params") {
    std::string sma = "CCC |$foo;;bar$|";
    SmartsParserParams ps;
    const std::unique_ptr<RWMol> mol(SmartsToMol(sma, ps));
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 3);
    CHECK(mol->getAtomWithIdx(0)->getProp<std::string>(
              common_properties::atomLabel) == "foo");
    CHECK(mol->getAtomWithIdx(2)->getProp<std::string>(
              common_properties::atomLabel) == "bar");
    CHECK(!mol->getAtomWithIdx(1)->hasProp(common_properties::atomLabel));
  }
  SECTION("no cxsmarts") {
    std::string sma = "CCC |$foo;;bar$|";
    SmartsParserParams ps;
    ps.allowCXSMILES = false;
    ps.parseName = false;
    const std::unique_ptr<RWMol> mol(SmartsToMol(sma, ps));
    REQUIRE(!mol);
  }
  SECTION("name") {
    std::string sma = "CCC foobar";
    SmartsParserParams ps;
    ps.parseName = true;
    const std::unique_ptr<RWMol> mol(SmartsToMol(sma, ps));
    REQUIRE(mol);
    REQUIRE(mol->getProp<std::string>(common_properties::_Name) == "foobar");
  }
  SECTION("writer") {
    auto mol = "CCC |$foo;;bar$|"_smarts;
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 3);
    CHECK(MolToSmarts(*mol) == "CCC");
    CHECK(MolToCXSmarts(*mol) == "CCC |$foo;;bar$|");
  }
  SECTION("writer, check reordering") {
    auto mol = "CC1.OC1 |$foo;;;bar$|"_smarts;
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 4);
    CHECK(MolToSmarts(*mol) == "CCCO");
    CHECK(MolToCXSmarts(*mol) == "CCCO |$foo;;bar;$|");
  }

  SECTION("parser, confirm enhanced stereo working") {
    auto mol = "[#6][C@]([#8])(F)Cl |&1:1|"_smarts;
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 5);
    CHECK(MolToSmarts(*mol) == "[#6][C@]([#8])(F)Cl");
    CHECK(MolToCXSmarts(*mol) == "[#6][C@]([#8])(F)Cl |&1:1|");

    {
      auto smol = "C[C@](O)(F)Cl |&1:1|"_smiles;
      REQUIRE(smol);
      SubstructMatchParameters sssparams;
      sssparams.useEnhancedStereo = true;
      sssparams.useChirality = true;
      CHECK(SubstructMatch(*smol, *mol, sssparams).size() == 1);
    }
    {
      auto smol = "C[C@](O)(F)Cl |o1:1|"_smiles;
      REQUIRE(smol);
      SubstructMatchParameters sssparams;
      sssparams.useEnhancedStereo = true;
      sssparams.useChirality = true;
      CHECK(SubstructMatch(*smol, *mol, sssparams).empty());
    }
  }

  SECTION("CXSMARTS parsing bug") {
    {  // no cxsmarts
      auto mol = "C[C@H]([F,Cl,Br])[C@H](C)[C@@H](C)Br"_smarts;
      REQUIRE(mol);
      CHECK(mol->getAtomWithIdx(2)->getQuery()->getDescription() == "AtomOr");
      CHECK(MolToSmarts(*mol) == "C[C@&H1]([F,Cl,Br])[C@&H1](C)[C@@&H1](C)Br");
      CHECK(MolToCXSmarts(*mol) ==
            "C[C@&H1]([F,Cl,Br])[C@&H1](C)[C@@&H1](C)Br");
    }
    {  // make sure that doesn't break anything
      auto mol = "C[C@H]([F,Cl,Br])[C@H](C)[C@@H](C)Br |a:1,o1:4,5|"_smarts;
      REQUIRE(mol);
      CHECK(mol->getAtomWithIdx(2)->getQuery()->getDescription() == "AtomOr");
      CHECK(MolToSmarts(*mol) == "C[C@&H1]([F,Cl,Br])[C@&H1](C)[C@@&H1](C)Br");
      CHECK(MolToCXSmarts(*mol) ==
            "C[C@&H1]([F,Cl,Br])[C@&H1](C)[C@@&H1](C)Br |a:1,o1:4,5|");
    }
  }
}

TEST_CASE("Github #4233: data groups in CXSMILES neither parsed nor written") {
  SECTION("basics") {
    {
      auto mol = "C/C=C/C |SgD:2,1:FIELD:info::::|"_smiles;
      REQUIRE(mol);
      const auto &sgs = getSubstanceGroups(*mol);
      REQUIRE(sgs.size() == 1);
      CHECK(sgs[0].getAtoms() == std::vector<unsigned int>{2, 1});
      CHECK(sgs[0].getProp<std::string>("TYPE") == "DAT");
      CHECK(sgs[0].getProp<std::string>("FIELDNAME") == "FIELD");
      CHECK(sgs[0].getProp<std::vector<std::string>>("DATAFIELDS") ==
            std::vector<std::string>{"info"});
      CHECK(MolToCXSmiles(*mol) == "C/C=C/C |SgD:2,1:FIELD:info::::|");
    }
    {
      auto mol = "C/C=C/C |SgD:2,1:FIELD:foo:like:info:tag:|"_smiles;
      REQUIRE(mol);
      const auto &sgs = getSubstanceGroups(*mol);
      REQUIRE(sgs.size() == 1);
      CHECK(sgs[0].getAtoms() == std::vector<unsigned int>{2, 1});
      CHECK(sgs[0].getProp<std::string>("TYPE") == "DAT");
      CHECK(sgs[0].getProp<std::string>("FIELDNAME") == "FIELD");
      CHECK(sgs[0].getProp<std::vector<std::string>>("DATAFIELDS") ==
            std::vector<std::string>{"foo"});
      CHECK(sgs[0].getProp<std::string>("QUERYOP") == "like");
      CHECK(sgs[0].getProp<std::string>("FIELDINFO") == "info");
      CHECK(sgs[0].getProp<std::string>("FIELDTAG") == "tag");
      CHECK(MolToCXSmiles(*mol) ==
            "C/C=C/C |SgD:2,1:FIELD:foo:like:info:tag:|");
    }
    {  // data on a dummy atom
      auto mol = "CC-* |$;;star_e$,SgD:2:querydata:val::::|"_smiles;
      REQUIRE(mol);
      const auto &sgs = getSubstanceGroups(*mol);
      REQUIRE(sgs.size() == 1);
      CHECK(sgs[0].getAtoms().size() == 1);
      CHECK(sgs[0].getAtoms()[0] == 2);
      CHECK(sgs[0].getProp<std::string>("TYPE") == "DAT");
      CHECK(sgs[0].getProp<std::vector<std::string>>("DATAFIELDS") ==
            std::vector<std::string>{"val"});
      CHECK(MolToCXSmiles(*mol) == "*CC |$star_e;;$,SgD:0:querydata:val::::|");
    }
    {
      auto mol =
          "C1=CC=C(C=C1)C1=CC=CC=C1 |c:0,2,4,9,11,t:7,SgD:8,9,11,10,7,6:PieceName:Ring1::::,SgD:1,2,3,4,5,0:PieceName:Ring2::::|"_smiles;
      REQUIRE(mol);
      const auto &sgs = getSubstanceGroups(*mol);
      REQUIRE(sgs.size() == 2);
      CHECK(sgs[0].getAtoms() == std::vector<unsigned int>{8, 9, 11, 10, 7, 6});
      CHECK(sgs[0].getProp<std::string>("TYPE") == "DAT");
      CHECK(sgs[0].getProp<std::string>("FIELDNAME") == "PieceName");
      CHECK(sgs[0].getProp<std::vector<std::string>>("DATAFIELDS") ==
            std::vector<std::string>{"Ring1"});
      CHECK(sgs[1].getAtoms() == std::vector<unsigned int>{1, 2, 3, 4, 5, 0});
      CHECK(sgs[1].getProp<std::string>("TYPE") == "DAT");
      CHECK(sgs[1].getProp<std::string>("FIELDNAME") == "PieceName");
      CHECK(sgs[1].getProp<std::vector<std::string>>("DATAFIELDS") ==
            std::vector<std::string>{"Ring2"});

      CHECK(MolToCXSmiles(*mol) ==
            "c1ccc(-c2ccccc2)cc1 "
            "|SgD:6,7,9,8,5,4:PieceName:Ring1::::,SgD:1,2,3,10,11,0:PieceName:"
            "Ring2::::|");

      // make sure we can round-trip that result through mol blocks
      auto mb = MolToV3KMolBlock(*mol);
      std::unique_ptr<RWMol> mol2(MolBlockToMol(mb));
      REQUIRE(mol2);
      const auto &sgs2 = getSubstanceGroups(*mol2);
      REQUIRE(sgs2.size() == 2);
      CHECK(sgs2[0].getAtoms() ==
            std::vector<unsigned int>{8, 9, 11, 10, 7, 6});
      CHECK(sgs2[0].getProp<std::string>("TYPE") == "DAT");
      CHECK(sgs2[0].getProp<std::string>("FIELDNAME") == "PieceName");
      CHECK(sgs2[0].getProp<std::vector<std::string>>("DATAFIELDS") ==
            std::vector<std::string>{"Ring1"});
      CHECK(sgs2[1].getAtoms() == std::vector<unsigned int>{1, 2, 3, 4, 5, 0});
      CHECK(sgs2[1].getProp<std::string>("TYPE") == "DAT");
      CHECK(sgs2[1].getProp<std::string>("FIELDNAME") == "PieceName");
      CHECK(sgs2[1].getProp<std::vector<std::string>>("DATAFIELDS") ==
            std::vector<std::string>{"Ring2"});
    }
  }
}

TEST_CASE("polymer SGroups") {
  SECTION("basics") {
    {
      auto mol =
          "*CC(*)C(*)N* |$star_e;;;star_e;;star_e;;star_e$,Sg:n:6,1,2,4::ht:6,0,:4,2,|"_smiles;
      REQUIRE(mol);
      const auto &sgs = getSubstanceGroups(*mol);
      REQUIRE(sgs.size() == 1);
      CHECK(sgs[0].getAtoms() == std::vector<unsigned int>{6, 1, 2, 4});
      CHECK(sgs[0].getBonds() == std::vector<unsigned int>{6, 0, 4, 2});
      CHECK(sgs[0].getProp<std::vector<unsigned int>>("XBHEAD") ==
            std::vector<unsigned int>{6, 0});
      CHECK(sgs[0].getProp<std::vector<unsigned int>>("XBCORR") ==
            std::vector<unsigned int>{6, 4, 0, 2});
      CHECK(sgs[0].getProp<std::string>("TYPE") == "SRU");
      CHECK(sgs[0].getProp<std::string>("CONNECT") == "HT");
      CHECK(sgs[0].getProp<unsigned int>("index") == 1);

      auto smi = MolToCXSmiles(*mol);
      CHECK(smi ==
            "*CC(*)C(*)N* "
            "|$star_e;;;star_e;;star_e;;star_e$,Sg:n:6,1,2,4::ht:6,0:4,2:|");
      // auto mb = MolToV3KMolBlock(*mol);
      // std::cerr << mb << std::endl;
    }
    {
      auto mol =
          "*CC(*)C(*)N* |$star_e;;;star_e;;star_e;;star_e$,Sg:n:6,1,2,4::hh&#44;f:6,0,:4,2,|"_smiles;
      REQUIRE(mol);
      const auto &sgs = getSubstanceGroups(*mol);
      REQUIRE(sgs.size() == 1);
      CHECK(sgs[0].getAtoms() == std::vector<unsigned int>{6, 1, 2, 4});
      CHECK(sgs[0].getBonds() == std::vector<unsigned int>{6, 0, 4, 2});
      CHECK(sgs[0].getProp<std::vector<unsigned int>>("XBHEAD") ==
            std::vector<unsigned int>{6, 0});
      CHECK(sgs[0].getProp<std::vector<unsigned int>>("XBCORR") ==
            std::vector<unsigned int>{6, 2, 0, 4});
      CHECK(sgs[0].getProp<std::string>("TYPE") == "SRU");
      CHECK(sgs[0].getProp<std::string>("CONNECT") == "HH");
      CHECK(sgs[0].getProp<unsigned int>("index") == 1);

      auto smi = MolToCXSmiles(*mol);
      CHECK(smi ==
            "*CC(*)C(*)N* "
            "|$star_e;;;star_e;;star_e;;star_e$,Sg:n:6,1,2,4::hh:6,0:2,4:|");
    }
    {  // minimal
      auto mol = "*-CCO-* |$star_e;;;;star_e$,Sg:n:1,2,3|"_smiles;
      REQUIRE(mol);
      const auto &sgs = getSubstanceGroups(*mol);
      REQUIRE(sgs.size() == 1);
      CHECK(sgs[0].getAtoms() == std::vector<unsigned int>{1, 2, 3});
      CHECK(sgs[0].getBonds() == std::vector<unsigned int>{0, 3});
      CHECK(sgs[0].getProp<std::vector<unsigned int>>("XBCORR") ==
            sgs[0].getBonds());
      CHECK(sgs[0].getProp<std::vector<unsigned int>>("XBHEAD") ==
            std::vector<unsigned int>{0});
      CHECK(sgs[0].getProp<std::string>("TYPE") == "SRU");
      CHECK(sgs[0].getProp<std::string>("CONNECT") == "EU");
      CHECK(sgs[0].getProp<unsigned int>("index") == 1);

      auto smi = MolToCXSmiles(*mol);
      CHECK(smi == "*CCO* |$star_e;;;;star_e$,Sg:n:1,2,3::eu:::|");
    }
    {  // single atom
      auto mol = "*-C-* |$star_e;;star_e$,Sg:n:1::ht|"_smiles;
      REQUIRE(mol);
      const auto &sgs = getSubstanceGroups(*mol);
      REQUIRE(sgs.size() == 1);
      CHECK(sgs[0].getAtoms() == std::vector<unsigned int>{1});
      CHECK(sgs[0].getBonds() == std::vector<unsigned int>{0, 1});
      CHECK(sgs[0].getProp<std::vector<unsigned int>>("XBCORR") ==
            sgs[0].getBonds());
      CHECK(sgs[0].getProp<std::vector<unsigned int>>("XBHEAD") ==
            std::vector<unsigned int>{0});
      CHECK(sgs[0].getProp<std::string>("TYPE") == "SRU");
      CHECK(sgs[0].getProp<std::string>("CONNECT") == "HT");
      CHECK(sgs[0].getProp<unsigned int>("index") == 1);

      auto smi = MolToCXSmiles(*mol);
      CHECK(smi == "*C* |$star_e;;star_e$,Sg:n:1::ht:::|");
    }
  }
  SECTION("multiple s groups") {
    {
      auto mol =
          "*-NCCO-* |$star_e;;;;;star_e$,Sg:n:1,2::ht,Sg:any:3,4::hh|"_smiles;
      REQUIRE(mol);
      const auto &sgs = getSubstanceGroups(*mol);
      REQUIRE(sgs.size() == 2);
      CHECK(sgs[0].getAtoms() == std::vector<unsigned int>{1, 2});
      CHECK(sgs[0].getBonds() == std::vector<unsigned int>{0, 2});
      CHECK(sgs[0].getProp<std::vector<unsigned int>>("XBCORR") ==
            sgs[0].getBonds());
      CHECK(sgs[0].getProp<std::vector<unsigned int>>("XBHEAD") ==
            std::vector<unsigned int>{0});
      CHECK(sgs[0].getProp<std::string>("TYPE") == "SRU");
      CHECK(sgs[0].getProp<std::string>("CONNECT") == "HT");
      CHECK(sgs[0].getProp<unsigned int>("index") == 1);

      CHECK(sgs[1].getAtoms() == std::vector<unsigned int>{3, 4});
      CHECK(sgs[1].getBonds() == std::vector<unsigned int>{2, 4});
      CHECK(sgs[1].getProp<std::vector<unsigned int>>("XBCORR") ==
            sgs[1].getBonds());
      CHECK(sgs[1].getProp<std::vector<unsigned int>>("XBHEAD") ==
            std::vector<unsigned int>{2});
      CHECK(sgs[1].getProp<std::string>("TYPE") == "ANY");
      CHECK(sgs[1].getProp<std::string>("CONNECT") == "HH");
      CHECK(sgs[1].getProp<unsigned int>("index") == 2);

      auto smi = MolToCXSmiles(*mol);
      CHECK(smi ==
            "*NCCO* |$star_e;;;;;star_e$,,,Sg:n:1,2::ht:::,Sg:any:3,4::hh:::|");
    }

    {  // multiple s groups + data
      auto mol =
          "CCNCCO-* "
          "|$;;;;;;star_e$,SgD:1:atomdata:val::::,Sg:n:2,3::ht,Sg:any:4,5::hh|"_smiles;
      REQUIRE(mol);
      const auto &sgs = getSubstanceGroups(*mol);
      REQUIRE(sgs.size() == 3);
      CHECK(sgs[0].getAtoms() == std::vector<unsigned int>{1});
      CHECK(sgs[0].getProp<std::string>("TYPE") == "DAT");
      CHECK(sgs[0].getProp<std::string>("FIELDNAME") == "atomdata");
      CHECK(sgs[0].getProp<std::vector<std::string>>("DATAFIELDS") ==
            std::vector<std::string>{"val"});
      CHECK(sgs[0].getProp<unsigned int>("index") == 1);

      CHECK(sgs[1].getAtoms() == std::vector<unsigned int>{2, 3});
      CHECK(sgs[1].getBonds() == std::vector<unsigned int>{1, 3});
      CHECK(sgs[1].getProp<std::vector<unsigned int>>("XBCORR") ==
            sgs[1].getBonds());
      CHECK(sgs[1].getProp<std::vector<unsigned int>>("XBHEAD") ==
            std::vector<unsigned int>{1});
      CHECK(sgs[1].getProp<std::string>("TYPE") == "SRU");
      CHECK(sgs[1].getProp<std::string>("CONNECT") == "HT");
      CHECK(sgs[1].getProp<unsigned int>("index") == 2);

      CHECK(sgs[2].getAtoms() == std::vector<unsigned int>{4, 5});
      CHECK(sgs[2].getBonds() == std::vector<unsigned int>{3, 5});
      CHECK(sgs[2].getProp<std::vector<unsigned int>>("XBCORR") ==
            sgs[2].getBonds());
      CHECK(sgs[2].getProp<std::vector<unsigned int>>("XBHEAD") ==
            std::vector<unsigned int>{3});
      CHECK(sgs[2].getProp<std::string>("TYPE") == "ANY");
      CHECK(sgs[2].getProp<std::string>("CONNECT") == "HH");
      CHECK(sgs[2].getProp<unsigned int>("index") == 3);

      auto smi = MolToCXSmiles(*mol);
      CHECK(smi ==
            "*OCCNCC "
            "|$star_e;;;;;;$,SgD:5:atomdata:val::::,,,,Sg:n:4,3::ht:::,Sg:any:"
            "2,1::hh:::|");
    }
  }
}

TEST_CASE("SGroup hierarchy") {
  SECTION("basics") {
    auto mol =
        "*-CNC(C-*)O-* "
        "|$star_e;;;;;star_e;;star_e$,Sg:any:2,1::ht,Sg:any:4,3,2,1,0,6::ht,"
        "SgH:1:0|"_smiles;
    REQUIRE(mol);
    const auto &sgs = getSubstanceGroups(*mol);
    REQUIRE(sgs.size() == 2);
    CHECK(sgs[0].getAtoms() == std::vector<unsigned int>{2, 1});
    CHECK(sgs[0].getProp<std::string>("TYPE") == "ANY");
    CHECK(sgs[0].getProp<unsigned int>("PARENT") == 2);
    CHECK(sgs[0].getProp<unsigned int>("index") == 1);
    CHECK(sgs[1].getAtoms() == std::vector<unsigned int>{4, 3, 2, 1, 0, 6});
    CHECK(sgs[1].getProp<std::string>("TYPE") == "ANY");
    CHECK(sgs[1].getProp<unsigned int>("index") == 2);
    CHECK(!sgs[1].hasProp("PARENT"));
    CHECK(MolToCXSmiles(*mol) ==
          "*CNC(C*)O* "
          "|$star_e;;;;;star_e;;star_e$,,,Sg:any:2,1::ht:::,Sg:any:4,3,2,1,0,6:"
          ":ht:::,SgH:1:0|");
  }
  SECTION("nested") {
    auto mol =
        "*-CNC(CC(-*)C-*)O-* "
        "|$star_e;;;;;;star_e;;star_e;;star_e$,"
        "SgD:4:internal data:val::::,SgD:7:atom value:value2::::,"
        "Sg:n:7::ht,Sg:n:2::ht,Sg:any:5,7,8,4,3,2,1,0,9::ht,"
        "SgH:4:2.3.0,2:1|"_smiles;
    REQUIRE(mol);
    const auto &sgs = getSubstanceGroups(*mol);
    REQUIRE(sgs.size() == 5);
    CHECK(sgs[0].getProp<unsigned int>("PARENT") == 5);
    CHECK(sgs[2].getProp<unsigned int>("PARENT") == 5);
    CHECK(sgs[3].getProp<unsigned int>("PARENT") == 5);
    CHECK(sgs[1].getProp<unsigned int>("PARENT") == 3);
    CHECK(
        MolToCXSmiles(*mol) ==
        "*CNC(CC(*)C*)O* |$star_e;;;;;;star_e;;star_e;;star_e$,SgD:4:internal "
        "data:val::::,SgD:7:atom "
        "value:value2::::,,,,,,Sg:n:7::ht:::,Sg:n:2::ht:::,Sg:any:5,7,8,4,3,2,"
        "1,0,9::ht:::,SgH:2:1,4:0.2.3|");
  }
}

TEST_CASE("Linknode writing") {
  SECTION("single") {
    auto mol = "OC1CCC(F)C1 |LN:1:1.3.2.6|"_smiles;
    REQUIRE(mol);
    std::string lns;
    CHECK(mol->getPropIfPresent(common_properties::molFileLinkNodes, lns));
    CHECK(lns == "1 3 2 2 3 2 7");
    CHECK(MolToCXSmiles(*mol) == "OC1CCC(F)C1 |LN:1:1.3.2.6|");
  }
  SECTION("multiple") {
    auto mol = "FC1CCC(O)C1 |LN:1:1.3.2.6,4:1.4.3.6|"_smiles;
    REQUIRE(mol);
    std::string lns;
    CHECK(mol->getPropIfPresent(common_properties::molFileLinkNodes, lns));
    CHECK(lns == "1 3 2 2 3 2 7|1 4 2 5 4 5 7");
    CHECK(MolToCXSmiles(*mol) == "OC1CCC(F)C1 |LN:4:1.3.3.6,1:1.4.2.6|");
  }
  SECTION("two-coordinate") {
    auto mol = "C1OCCC1C |LN:0:1.5,1:1.3|"_smiles;
    REQUIRE(mol);
    CHECK(MolToCXSmiles(*mol) == "CC1CCOC1 |LN:5:1.5,4:1.3|");
  }
}

TEST_CASE("smilesBondOutputOrder") {
  SECTION("basics") {
    auto m = "OCCN.CCO"_smiles;
    REQUIRE(m);
    CHECK(MolToSmiles(*m) == "CCO.NCCO");
    std::vector<unsigned int> order;
    m->getProp(common_properties::_smilesAtomOutputOrder, order);
    CHECK(order == std::vector<unsigned int>{4, 5, 6, 3, 2, 1, 0});
    m->getProp(common_properties::_smilesBondOutputOrder, order);
    CHECK(order == std::vector<unsigned int>{3, 4, 2, 1, 0});
  }
  SECTION("Github #5585: incorrect ordering for ring closures") {
    auto m = "OC1CCCN1"_smiles;
    REQUIRE(m);
    CHECK(MolToSmiles(*m) == "OC1CCCN1");
    std::vector<unsigned int> order;
    m->getProp(common_properties::_smilesAtomOutputOrder, order);
    CHECK(order == std::vector<unsigned int>{0, 1, 2, 3, 4, 5});
    m->getProp(common_properties::_smilesBondOutputOrder, order);
    CHECK(order == std::vector<unsigned int>{0, 1, 2, 3, 4, 5});
  }
}

TEST_CASE("Github #4320: Support toggling components of CXSMILES output") {
  SECTION("sgroups") {
    auto mol =
        "*-CNC(CC(-*)C-*)O-* "
        "|$star_e;;;;;;star_e;;star_e;;star_e$,"
        "SgD:4:internal data:val::::,SgD:7:atom value:value2::::,"
        "Sg:n:7::ht,Sg:n:2::ht,Sg:any:5,7,8,4,3,2,1,0,9::ht,"
        "SgH:4:2.3.0,2:1|"_smiles;
    SmilesWriteParams ps;
    {
      auto cxsmi = MolToCXSmiles(*mol, ps, SmilesWrite::CXSmilesFields::CX_ALL);
      CHECK(cxsmi ==
            "*CNC(CC(*)C*)O* "
            "|$star_e;;;;;;star_e;;star_e;;star_e$,SgD:4:internal "
            "data:val::::,SgD:7:atom "
            "value:value2::::,,,,,,Sg:n:7::ht:::,Sg:n:2::ht:::,Sg:any:5,7,8,4,"
            "3,2,"
            "1,0,9::ht:::,SgH:2:1,4:0.2.3|");
      CHECK(std::unique_ptr<ROMol>(SmilesToMol(cxsmi)));
    }
    {
      auto cxsmi =
          MolToCXSmiles(*mol, ps, SmilesWrite::CXSmilesFields::CX_NONE);
      CHECK(cxsmi == "*CNC(CC(*)C*)O*");
      CHECK(std::unique_ptr<ROMol>(SmilesToMol(cxsmi)));
    }
    {
      auto cxsmi =
          MolToCXSmiles(*mol, ps,
                        SmilesWrite::CXSmilesFields::CX_ALL ^
                            SmilesWrite::CXSmilesFields::CX_ATOM_LABELS);
      CHECK(
          cxsmi ==
          "*CNC(CC(*)C*)O* |SgD:4:internal data:val::::,SgD:7:atom "
          "value:value2::::,,,,,,Sg:n:7::ht:::,Sg:n:2::ht:::,Sg:any:5,7,8,4,3,"
          "2,1,0,9::ht:::,SgH:2:1,4:0.2.3|");
      CHECK(std::unique_ptr<ROMol>(SmilesToMol(cxsmi)));
    }
    {
      auto cxsmi = MolToCXSmiles(*mol, ps,
                                 SmilesWrite::CXSmilesFields::CX_ALL ^
                                     SmilesWrite::CXSmilesFields::CX_SGROUPS);
      CHECK(cxsmi ==
            "*CNC(CC(*)C*)O* "
            "|$star_e;;;;;;star_e;;star_e;;star_e$,,,Sg:n:7::ht:::,Sg:n:2::ht::"
            ":,Sg:any:5,7,8,4,3,2,1,0,9::ht:::,SgH:2:0.1|");
      CHECK(std::unique_ptr<ROMol>(SmilesToMol(cxsmi)));
    }
    {
      auto cxsmi = MolToCXSmiles(*mol, ps,
                                 SmilesWrite::CXSmilesFields::CX_ALL ^
                                     SmilesWrite::CXSmilesFields::CX_POLYMER);
      CHECK(cxsmi ==
            "*CNC(CC(*)C*)O* "
            "|$star_e;;;;;;star_e;;star_e;;star_e$,SgD:4:internal "
            "data:val::::,SgD:7:atom value:value2::::,,,|");
      CHECK(std::unique_ptr<ROMol>(SmilesToMol(cxsmi)));
    }
  }
  SECTION("coordinates") {
    auto m = "CC |(0,.75,;0,-.75,)|"_smiles;
    REQUIRE(m);
    SmilesWriteParams ps;
    auto cxsmi = MolToCXSmiles(*m, ps,
                               SmilesWrite::CXSmilesFields::CX_ALL ^
                                   SmilesWrite::CXSmilesFields::CX_COORDS);
    CHECK(cxsmi == "CC");
    CHECK(std::unique_ptr<ROMol>(SmilesToMol(cxsmi)));
  }
  SECTION("enhanced stereo") {
    auto m = "C[C@H](F)Cl |o1:1|"_smiles;
    REQUIRE(m);
    SmilesWriteParams ps;
    auto cxsmi =
        MolToCXSmiles(*m, ps,
                      SmilesWrite::CXSmilesFields::CX_ALL ^
                          SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO);
    CHECK(cxsmi == "C[C@H](F)Cl");
    CHECK(std::unique_ptr<ROMol>(SmilesToMol(cxsmi)));
  }
  SECTION("link nodes") {
    auto m = "OC1CCC(F)C1 |LN:1:1.3.2.6|"_smiles;
    REQUIRE(m);
    SmilesWriteParams ps;
    auto cxsmi = MolToCXSmiles(*m, ps,
                               SmilesWrite::CXSmilesFields::CX_ALL ^
                                   SmilesWrite::CXSmilesFields::CX_LINKNODES);
    CHECK(cxsmi == "OC1CCC(F)C1");
    CHECK(std::unique_ptr<ROMol>(SmilesToMol(cxsmi)));
  }
  SECTION("radicals") {
    auto m = "[O][C][O] |^1:0,2|"_smiles;
    REQUIRE(m);
    SmilesWriteParams ps;
    auto cxsmi = MolToCXSmiles(*m, ps,
                               SmilesWrite::CXSmilesFields::CX_ALL ^
                                   SmilesWrite::CXSmilesFields::CX_RADICALS);
    CHECK(cxsmi == "[O][C][O]");
    CHECK(std::unique_ptr<ROMol>(SmilesToMol(cxsmi)));
  }
  SECTION("values") {
    auto m = "COCC |$_AV:;bar;;foo$|"_smiles;
    REQUIRE(m);
    SmilesWriteParams ps;
    auto cxsmi =
        MolToCXSmiles(*m, ps,
                      SmilesWrite::CXSmilesFields::CX_ALL ^
                          SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES);
    CHECK(cxsmi == "CCOC");
  }
}

TEST_CASE("non-tetrahedral chirality") {
  SECTION("allowed values") {
    {
      auto m = "C[Fe@SP3](Cl)(F)N"_smiles;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_SQUAREPLANAR);
      CHECK(m->getAtomWithIdx(1)->getProp<unsigned int>(
                common_properties::_chiralPermutation) == 3);
    }
    {
      auto m = "C[Fe@SP3](Cl)(F)N"_smiles;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_SQUAREPLANAR);
      CHECK(m->getAtomWithIdx(1)->getProp<unsigned int>(
                common_properties::_chiralPermutation) == 3);
    }
  }
}
TEST_CASE(
    "Github #4503: MolFromSmiles and MolFromSmarts incorrectly accepting input "
    "with spaces") {
  SECTION("SMILES defaults") {
    std::unique_ptr<RWMol> m{SmilesToMol("NON sense extra")};
    CHECK(m);
    CHECK(m->hasProp("_Name"));
    CHECK(m->getProp<std::string>("_Name") == "sense extra");

    CHECK_THROWS_AS(SmilesToMol("NON |sense|"), SmilesParseException);
  }
  SECTION("SMARTS defaults") {
    std::unique_ptr<RWMol> m{SmilesToMol("NON sense extra")};
    CHECK(m);
    CHECK(m->hasProp("_Name"));
    CHECK(m->getProp<std::string>("_Name") == "sense extra");

    CHECK_THROWS_AS(SmartsToMol("NON |sense|"), SmilesParseException);
  }
  SECTION("SMILES no names") {
    SmilesParserParams ps;
    ps.parseName = false;
    CHECK_THROWS_AS(SmilesToMol("NON sense", ps), SmilesParseException);
  }
  SECTION("SMARTS no names") {
    SmartsParserParams ps;
    ps.parseName = false;
    CHECK_THROWS_AS(SmartsToMol("NON sense", ps), SmilesParseException);
  }
  SECTION("SMILES names without parsing CXSMILES") {
    SmilesParserParams ps;
    ps.allowCXSMILES = false;
    {
      std::unique_ptr<RWMol> m{SmilesToMol("NON sense extra", ps)};
      CHECK(m);
      CHECK(m->hasProp("_Name"));
      CHECK(m->getProp<std::string>("_Name") == "sense extra");
    }
    {
      std::unique_ptr<RWMol> m{SmilesToMol("NON |$N1;O2;N3$| sense extra", ps)};
      CHECK(m);
      CHECK(m->hasProp("_Name"));
      CHECK(m->getProp<std::string>("_Name") == "|$N1;O2;N3$| sense extra");
    }
  }

  SECTION("SMILES not strict") {
    SmilesParserParams ps;
    ps.strictCXSMILES = false;
    {
      std::unique_ptr<RWMol> m{SmilesToMol("NON sense", ps)};
      CHECK(m);
      CHECK(m->hasProp("_Name"));
    }
    {
      std::unique_ptr<RWMol> m{SmilesToMol("NON |sense|", ps)};
      CHECK(m);
      CHECK(!m->hasProp("_Name"));
    }
    ps.parseName = false;
    {
      std::unique_ptr<RWMol> m{SmilesToMol("NON sense", ps)};
      CHECK(m);
      CHECK(!m->hasProp("_Name"));
    }
    {
      std::unique_ptr<RWMol> m{SmilesToMol("NON |sense|", ps)};
      CHECK(m);
      CHECK(!m->hasProp("_Name"));
    }
  }
  SECTION("SMARTS not strict") {
    SmartsParserParams ps;
    ps.strictCXSMILES = false;
    {
      std::unique_ptr<RWMol> m{SmartsToMol("NON sense", ps)};
      CHECK(m);
      CHECK(m->hasProp("_Name"));
    }
    {
      std::unique_ptr<RWMol> m{SmartsToMol("NON |sense|", ps)};
      CHECK(m);
      CHECK(!m->hasProp("_Name"));
    }
    ps.parseName = false;
    {
      std::unique_ptr<RWMol> m{SmartsToMol("NON sense", ps)};
      CHECK(m);
      CHECK(!m->hasProp("_Name"));
    }
    {
      std::unique_ptr<RWMol> m{SmartsToMol("NON |sense|", ps)};
      CHECK(m);
      CHECK(!m->hasProp("_Name"));
    }
  }
  SECTION("SMARTS CXExtensions + names") {
    SmartsParserParams ps;
    ps.strictCXSMILES = false;
    {  // CXSMILES + name
      std::unique_ptr<RWMol> m{SmartsToMol("NON |$_AV:bar;;foo$| name", ps)};
      CHECK(m);
      CHECK(m->hasProp("_Name"));
      CHECK(m->getProp<std::string>("_Name") == "name");
    }
    {  // CXSMILES fails, so we don't read the name
      std::unique_ptr<RWMol> m{SmartsToMol("NON |sense| name", ps)};
      CHECK(m);
      CHECK(!m->hasProp("_Name"));
    }
    ps.parseName = false;
    {  // CXSMILES, skip the name
      std::unique_ptr<RWMol> m{SmartsToMol("NON |$_AV:bar;;foo$| name", ps)};
      CHECK(m);
      CHECK(!m->hasProp("_Name"));
    }
    ps.parseName = true;
    ps.allowCXSMILES = false;
    ps.strictCXSMILES = true;
    {  // CXSMILES + name, but not parsing the CXSMILES, so we read it as the
       // name
      std::unique_ptr<RWMol> m{SmartsToMol("NON |$_AV:bar;;foo$| name", ps)};
      CHECK(m);
      CHECK(m->hasProp("_Name"));
      CHECK(m->getProp<std::string>("_Name") == "|$_AV:bar;;foo$| name");
    }
  }
  SECTION("SMILES CXExtensions + names") {
    SmilesParserParams ps;
    ps.strictCXSMILES = false;
    {  // CXSMILES + name
      std::unique_ptr<RWMol> m{SmilesToMol("NON |$_AV:bar;;foo$| name", ps)};
      CHECK(m->hasProp("_Name"));
      CHECK(m->getProp<std::string>("_Name") == "name");
    }
    {  // CXSMILES fails, so we don't read the name
      std::unique_ptr<RWMol> m{SmilesToMol("NON |sense| name", ps)};
      CHECK(m);
      CHECK(!m->hasProp("_Name"));
    }
    ps.parseName = false;
    {  // CXSMILES, skip the name
      std::unique_ptr<RWMol> m{SmilesToMol("NON |$_AV:bar;;foo$| name", ps)};
      CHECK(m);
      CHECK(!m->hasProp("_Name"));
    }
    ps.parseName = true;
    ps.allowCXSMILES = false;
    ps.strictCXSMILES = true;
    {  // CXSMILES + name, but not parsing the CXSMILES, so we read it as the
       // name
      std::unique_ptr<RWMol> m{SmilesToMol("NON |$_AV:bar;;foo$| name", ps)};
      CHECK(m);
      CHECK(m->hasProp("_Name"));
      CHECK(m->getProp<std::string>("_Name") == "|$_AV:bar;;foo$| name");
    }
  }
}

TEST_CASE("Github #4582: double bonds and ring closures") {
  const auto useLegacy = GENERATE(true, false);
  CAPTURE(useLegacy);
  UseLegacyStereoPerceptionFixture fxn(useLegacy);
  auto mol = R"CTAB(CHEMBL409450
     RDKit          2D

 22 25  0  0  0  0  0  0  0  0999 V2000
   -1.1669    1.3591    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6820    0.6916    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9516    1.1041    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9516    0.2791    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5647    1.6562    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3931    2.4631    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6085    2.7181    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9954    2.1661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1669    0.0242    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9120   -0.7604    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3969   -1.4279    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9120   -2.0953    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1274   -1.8404    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1274   -1.0154    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5871   -2.2529    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3016   -1.8404    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3016   -1.0154    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5871   -0.6029    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2219   -1.4279    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1430    0.6916    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.5555    1.4061    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0160   -2.2529    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  2  0
  1  8  1  0
  9  2  1  0
  2 20  2  0
  3  4  1  0
  3  5  1  0
  9  4  1  0
  5  6  2  0
  6  7  1  0
  8  7  2  0
  9 10  2  0
 10 11  1  0
 10 14  1  0
 11 12  1  0
 11 19  2  0
 13 12  1  0
 13 14  2  0
 13 15  1  0
 14 18  1  0
 15 16  2  0
 16 17  1  0
 22 16  1  0
 17 18  2  0
 20 21  1  0
M  END)CTAB"_ctab;
  REQUIRE(mol);
  auto dbond = mol->getBondBetweenAtoms(1, 19);
  REQUIRE(dbond);
  CHECK(dbond->getBondType() == Bond::BondType::DOUBLE);
  if (useLegacy) {
    CHECK(dbond->getStereo() == Bond::BondStereo::STEREOE);
    CHECK(dbond->getStereoAtoms() == std::vector<int>{8, 20});
  } else {
    CHECK(dbond->getStereo() == Bond::BondStereo::STEREOCIS);
    CHECK(dbond->getStereoAtoms() == std::vector<int>{0, 20});
  }
  auto csmiles = MolToSmiles(*mol);
  CHECK(csmiles == "O=C1Nc2cc(Br)ccc2/C1=C1/Nc2ccccc2/C1=N\\O");

#if 0
  // FIX: this test fails. The SMILES generated for this random permutation
  // has the stereo for one of the double bonds wrong.
  SECTION("specific random output") {
    SmilesWriteParams ps;
    ps.doRandom = true;
    getRandomGenerator(51)();
    auto rsmiles = MolToSmiles(*mol, ps);
    std::unique_ptr<RWMol> mol2(SmilesToMol(rsmiles));
    REQUIRE(mol2);
    auto smiles = MolToSmiles(*mol2);
    if (smiles != csmiles) {
      std::cerr << smiles << std::endl;
    }
    CHECK(smiles == csmiles);
  }
#endif

  SECTION("bulk random output order") {
    auto csmiles = MolToSmiles(*mol);
    CHECK(csmiles == "O=C1Nc2cc(Br)ccc2/C1=C1/Nc2ccccc2/C1=N\\O");
    SmilesWriteParams ps;
    ps.doRandom = true;
    for (auto i = 0u; i < 100; ++i) {
      if (i == 50) {
        // we know this one fails (it's explicitly tested above)
        continue;
      }
      getRandomGenerator(i + 1)();
      auto rsmiles = MolToSmiles(*mol, ps);
      std::unique_ptr<RWMol> mol2(SmilesToMol(rsmiles));
      REQUIRE(mol2);
      auto smiles = MolToSmiles(*mol2);
      if (smiles != csmiles) {
        std::cerr << ">>> " << i << " " << rsmiles << std::endl;
      }
      CHECK(smiles == csmiles);
    }
  }
}

TEST_CASE("Github #4582 continued: double bonds and ring closures") {
  SECTION("basics") {
    auto mol = R"CTAB(
  Mrv2108 10072106112D

 11 11  0  0  0  0            999 V2000
    1.2097   -1.0517    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3651   -1.1116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0571    0.9536    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2670   -0.5667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7001    0.4032    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3269    0.2561    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2180    0.8882    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7601   -0.4305    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5763   -1.7848    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1262   -2.4676    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4920   -3.1990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  8  1  1  0  0  0  0
  1  2  2  0  0  0  0
  5  3  1  0  0  0  0
  7  3  1  0  0  0  0
  4  2  1  0  0  0  0
  5  8  1  0  0  0  0
  6  7  1  0  0  0  0
  1  9  1  0  0  0  0
  9 10  2  0  0  0  0
 10 11  1  0  0  0  0
  6  4  1  0  0  0  0
M  END)CTAB"_ctab;
    REQUIRE(mol);
    auto csmiles = MolToSmiles(*mol);
    CHECK(csmiles == "C/N=C/C1=C/CCCCCC1");
  }
  SECTION("github #3967 part 1") {
    // duplicating a test from elsewhere
    auto mol = "C=c1s/c2n(c1=O)CCCCCCC\\N=2"_smiles;
    REQUIRE(mol);
    auto smi = MolToSmiles(*mol);
    CHECK(smi == "C=c1s/c2n(c1=O)CCCCCCC\\N=2");
  }
  SECTION("github #3967 part 2") {
    // duplicating a test from elsewhere
    auto mol = R"SMI(C1=CC/C=C2C3=C/CC=CC=CC\3C\2C=C1)SMI"_smiles;
    REQUIRE(mol);
    auto smi = MolToSmiles(*mol);
    CHECK(smi == R"SMI(C1=CC/C=C2C3=C\CC=CC=CC/3C\2C=C1)SMI");
  }
  SECTION("CHEMBL3623347") {
    auto mol = R"CTAB(CHEMBL3623347
     RDKit          2D

 44 47  0  0  0  0  0  0  0  0999 V2000
   -2.0000    1.0700    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1400   -0.4700    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6400   -0.8300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4400    0.4800    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4200    1.6600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2200   -1.9300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6700   -2.0400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9400    1.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4900   -1.0400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1200    0.7400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5400    0.1400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5400    1.3300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6000    0.4700    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4000    1.6300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4000    2.9600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8800    2.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6604    2.1462    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6965    1.6501    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2764    1.7277    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9496    2.6967    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4143   -0.1105    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6960    2.8278    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2128   -2.2139    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2997   -3.4050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7586   -4.5137    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1099   -3.2485    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7200   -1.1200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2300   -0.7900    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1740   -2.2308    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8927   -3.2754    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0667   -4.5284    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.7380   -5.8707    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9120   -7.1238    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5833   -8.4661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7573   -9.7192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4287  -11.0615    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9267  -11.1538    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4623  -12.2277    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5892  -10.1533    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6027  -12.3146    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.2740  -13.6569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4719  -13.7282    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6136  -14.6588    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5388   -1.7411    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  6
  3  2  1  0
  3  4  1  0
  4  5  1  0
  5  1  1  0
 28  6  1  0
  6  7  2  0
 10  8  1  0
  8 14  1  0
 13  9  1  0
  9  7  1  0
 10 28  1  0
 27 11  1  0
 11 12  1  0
 12 10  1  0
 13 14  1  0
 14 15  1  0
 15 16  1  0
 16  1  1  0
  1 13  1  0
 12 17  1  1
 12 18  1  0
 10 19  1  1
 14 20  1  1
 13 21  1  6
  5 22  1  6
  3 23  1  1
 23 24  2  0
 24 25  1  0
 24 26  1  0
 27 28  1  0
 27 29  2  0
  6 30  1  0
 30 31  2  0
 31 32  1  0
 32 33  1  0
 33 34  1  0
 34 35  1  0
 35 36  1  0
 36 37  1  0
 37 38  2  0
 37 39  1  0
 36 40  1  1
 40 41  1  0
 41 42  1  0
 41 43  2  0
 28 44  1  1
M  END)CTAB"_ctab;
    REQUIRE(mol);
    auto csmiles = MolToSmiles(*mol);
    // clang-format off
    CHECK(
        csmiles ==
        "CC(=O)N[C@@H](CCCC/N=C/C1=C/C[C@@H]2[C@](C)(CC[C@@]23O[C@@H](C=C(C)C)C[C@@H]3C)C[C@H]2[C@@H]1C(=O)C[C@@]2(C)O)C(=O)O");
  }
  // clang-format on
}

TEST_CASE(
    "unspecified bonds at beginning of branch in SMARTS not properly tagged") {
  auto m = "C(C)C"_smarts;
  REQUIRE(m);
  CHECK(m->getBondWithIdx(0)->getQuery()->getDescription() ==
        "SingleOrAromaticBond");
  CHECK(m->getBondWithIdx(1)->getQuery()->getDescription() ==
        "SingleOrAromaticBond");
}

TEST_CASE("Github #4878: cannot parse coordinate bonds from CXSMARTS",
          "[cxsmiles]") {
  SECTION("basics") {
    auto m = "[Fe]OC |C:1.0|"_smarts;
    REQUIRE(m);
    CHECK(m->getBondWithIdx(0)->getBondType() == Bond::BondType::DATIVE);
  }
  SECTION("branches") {
    {
      auto m = "[Fe](O)OC |C:2.1|"_smarts;
      REQUIRE(m);
      CHECK(m->getBondWithIdx(1)->getBondType() == Bond::BondType::DATIVE);
    }
    {
      auto m = "[Fe](OC)O |C:1.0|"_smarts;
      REQUIRE(m);
      CHECK(m->getBondWithIdx(0)->getBondType() == Bond::BondType::DATIVE);
    }
  }
  SECTION("rings") {
    {
      auto m = "O1CC[Fe]1 |C:0.3|"_smarts;
      REQUIRE(m);
      CHECK(m->getBondBetweenAtoms(0, 3)->getBondType() ==
            Bond::BondType::DATIVE);
    }
  }
}

TEST_CASE("Github #4981: Invalid SMARTS for negated single-atoms", "[smarts]") {
  SECTION("basics") {
    auto mol = "[!C]"_smarts;
    REQUIRE(mol);
    CHECK(SmartsWrite::GetAtomSmarts(mol->getAtomWithIdx(0)) == "[!C]");
  }
  SECTION("as reported") {
    auto mol = "[NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]"_smarts;
    REQUIRE(mol);
    CHECK(MolToSmarts(*mol).find("N!C") == std::string::npos);
  }
}

TEST_CASE("Pol and Mod atoms in CXSMILES", "[cxsmiles]") {
  SECTION("Pol basics") {
    auto mol = "CC* |$;;Pol_p$|"_smiles;
    REQUIRE(mol);
    std::string val;
    CHECK(mol->getAtomWithIdx(2)->getPropIfPresent(
        common_properties::dummyLabel, val));
    CHECK(val == "Pol");
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "*CC |$Pol_p;;$|");
  }
  SECTION("Mod basics") {
    auto mol = "CC* |$;;Mod_p$|"_smiles;
    REQUIRE(mol);
    std::string val;
    CHECK(mol->getAtomWithIdx(2)->getPropIfPresent(
        common_properties::dummyLabel, val));
    CHECK(val == "Mod");
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "*CC |$Mod_p;;$|");
  }
}

TEST_CASE("empty atom label block", "[cxsmiles]") {
  SECTION("basics") {
    auto m = R"CTAB(
  MJ201100

  8  8  0  0  0  0  0  0  0  0999 V2000
   -1.0491    1.5839    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7635    1.1714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7635    0.3463    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0491   -0.0661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3346    0.3463    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3346    1.1714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3798    1.5839    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
    0.3798   -0.0661    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  6  7  1  0  0  2  0
  5  8  1  0  0  2  0
M  RGP  2   7   1   8   2
M  END)CTAB"_ctab;
    REQUIRE(m);
    m->clearConformers();
    auto smi = MolToCXSmiles(*m);
    CHECK(smi.find("$") == std::string::npos);
  }
}

TEST_CASE("Github #5372: errors with fragments and doRandom=True") {
  SECTION("basics") {
    auto m = "C.C"_smiles;
    REQUIRE(m);
    SmilesWriteParams ps;
    ps.doRandom = true;
    auto smi = MolToSmiles(*m, ps);
    CHECK(smi == "C.C");
  }
}

TEST_CASE("github #5466 writing floating point atom props cxsmiles",
          "[smiles][cxsmiles]") {
  SECTION("simple") {
    auto mol = "C"_smiles;
    REQUIRE(mol);

    mol->getAtomWithIdx(0)->setProp<std::string>("foo", "7.6");
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "C |atomProp:0.foo.7&#46;6|");

    // 7.5 is exactly representable in IEEE so this helps :)
    mol->getAtomWithIdx(0)->setProp<double>("foo", 7.5);
    smi = MolToCXSmiles(*mol);
    CHECK(smi == "C |atomProp:0.foo.7&#46;5|");
  }
  SECTION("label with .") {
    auto mol = "C"_smiles;
    REQUIRE(mol);

    mol->getAtomWithIdx(0)->setProp<int>("foo.foo", 7);
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "C |atomProp:0.foo&#46;foo.7|");
  }
}

TEST_CASE("CXSMILES and 3d conformers") {
  SECTION("detect 3d and perceive stereo") {
    auto m =
        "CC(O)(F)Cl |(0.856693,-0.0843716,-0.11692;-0.539611,0.360041,0.18126;-1.47941,-0.628083,-0.170155;-0.665621,0.560167,1.55524;-0.999017,1.83232,-0.693661)|"_smiles;
    REQUIRE(m);
    REQUIRE(m->getNumConformers() == 1);
    CHECK(m->getConformer().is3D());
    CHECK(m->getAtomWithIdx(1)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
  SECTION("detect 2d") {
    auto m =
        "CC(O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,)|"_smiles;
    REQUIRE(m);
    REQUIRE(m->getNumConformers() == 1);
    CHECK(!m->getConformer().is3D());
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
  SECTION("detect 2d when z coords are all zero") {
    auto m =
        "CC(O)Cl |(-3.9163,5.4767,0;-3.9163,3.9367,0;-2.5826,3.1667,0;-5.25,3.1667,0)|"_smiles;
    REQUIRE(m);
    REQUIRE(m->getNumConformers() == 1);
    CHECK(!m->getConformer().is3D());
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
}

TEST_CASE("wiggly and wedged bonds in CXSMILES") {
  SECTION("basic reading/writing") {
    auto m = "CC(O)F |w:1.2|"_smiles;
    REQUIRE(m);
    unsigned int bondcfg = 0;
    CHECK(m->getBondWithIdx(2)->getPropIfPresent("_MolFileBondCfg", bondcfg));
    CHECK(bondcfg == 2);
    // make sure we end up with a wiggly bond in output mol blocks:
    Chirality::reapplyMolBlockWedging(*m);
    auto mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("CFG=2") != std::string::npos);
    // make sure we end up with the wiggly bond in the output CXSMILES:
    auto cxsmi = MolToCXSmiles(*m, SmilesWriteParams(),
                               SmilesWrite::CXSmilesFields::CX_ALL,
                               RestoreBondDirOptionTrue);
    CHECK(cxsmi == "CC(O)F |w:1.2|");
    // but we can turn that off
    SmilesWriteParams ps;
    cxsmi =
        MolToCXSmiles(*m, ps, SmilesWrite::CX_ALL ^ SmilesWrite::CX_BOND_CFG);
    CHECK(cxsmi == "CC(O)F");
  }

  SECTION("CXSMILES wiggly bond over-rides atomic stereo") {
    auto m = "C[C@H](O)F |w:1.2|"_smiles;
    REQUIRE(m);
    unsigned int bondcfg = 0;
    CHECK(m->getBondWithIdx(2)->getPropIfPresent("_MolFileBondCfg", bondcfg));
    CHECK(bondcfg == 2);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);

    // make sure we end up with a wiggly bond in output mol blocks:
    Chirality::reapplyMolBlockWedging(*m);
    auto mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("CFG=2") != std::string::npos);
    // make sure we end up with the wiggly bond in the output CXSMILES:
    auto cxsmi = MolToCXSmiles(*m, SmilesWriteParams(),
                               SmilesWrite::CXSmilesFields::CX_ALL,
                               RestoreBondDirOptionTrue);
    CHECK(cxsmi == "CC(O)F |w:1.2|");
  }
  SECTION("make sure order gets reversed when needed") {
    auto m = "CC(O)Cl |w:1.0|"_smiles;
    REQUIRE(m);
    CHECK(m->getBondWithIdx(0)->getBeginAtomIdx() == 1);
    unsigned int bondcfg = 0;
    CHECK(m->getBondWithIdx(0)->getPropIfPresent("_MolFileBondCfg", bondcfg));
    CHECK(bondcfg == 2);
  }

  SECTION("make sure wedging gets applied when coordinates are there") {
    auto m =
        "CC(O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,),wU:1.0|"_smiles;
    REQUIRE(m);
    unsigned int bondcfg = 0;
    CHECK(m->getBondWithIdx(0)->getPropIfPresent("_MolFileBondCfg", bondcfg));
    CHECK(bondcfg == 1);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
    m = "CC(O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,),wD:1.0|"_smiles;
    REQUIRE(m);
    CHECK(m->getBondWithIdx(0)->getPropIfPresent("_MolFileBondCfg", bondcfg));
    CHECK(bondcfg == 3);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    invertMolBlockWedgingInfo(*m);
    CHECK(m->getBondWithIdx(0)->getPropIfPresent("_MolFileBondCfg", bondcfg));
    CHECK(bondcfg == 1);
    Chirality::reapplyMolBlockWedging(*m);
    MolOps::assignChiralTypesFromBondDirs(*m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
    invertMolBlockWedgingInfo(*m);
    CHECK(m->getBondWithIdx(0)->getPropIfPresent("_MolFileBondCfg", bondcfg));
    CHECK(bondcfg == 3);
    Chirality::reapplyMolBlockWedging(*m);
    MolOps::assignChiralTypesFromBondDirs(*m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    Chirality::clearMolBlockWedgingInfo(*m);
    m->getAtomWithIdx(1)->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(!m->getBondWithIdx(0)->getPropIfPresent("_MolFileBondCfg", bondcfg));
    Chirality::reapplyMolBlockWedging(*m);
    MolOps::assignChiralTypesFromBondDirs(*m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
  }

  SECTION("writing examples") {
    auto m = "OC(C)Cl"_smiles;
    REQUIRE(m);
    {
      ROMol nm(*m);
      nm.getBondWithIdx(1)->setBondDir(Bond::BondDir::UNKNOWN);
      auto cxsmi = MolToCXSmiles(nm, SmilesWriteParams(),
                                 SmilesWrite::CXSmilesFields::CX_ALL,
                                 RestoreBondDirOptionTrue);
      CHECK(cxsmi == "CC(O)Cl |w:1.0|");
    }
    {
      ROMol nm(*m);
      nm.getBondWithIdx(1)->setProp(common_properties::_MolFileBondCfg, 2);
      auto cxsmi = MolToCXSmiles(nm, SmilesWriteParams(),
                                 SmilesWrite::CXSmilesFields::CX_ALL,
                                 RestoreBondDirOptionClear);
      CHECK(cxsmi == "CC(O)Cl |w:1.0|");
    }
  }

  SECTION("writing wedges and dashes") {
    auto m = R"CTAB(
  RDKit             2D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Cl -5.25 3.1667 0 0
M  V30 2 C -3.9163 3.9367 0 0
M  V30 3 O -2.5826 3.1667 0 0
M  V30 4 C -3.9163 5.4767 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    auto cxsmi = MolToCXSmiles(*m);
    CHECK(cxsmi.find("wU:1.0") != std::string::npos);
    SmilesWriteParams ps;
    cxsmi =
        MolToCXSmiles(*m, ps, SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS);
    CHECK(cxsmi.find("wU:1.0") == std::string::npos);
    // change the bond dir. This also tests that the wedging overrides the
    // CFG property
    m->getBondWithIdx(2)->setBondDir(Bond::BondDir::BEGINDASH);
    cxsmi = MolToCXSmiles(*m, SmilesWriteParams(),
                          SmilesWrite::CXSmilesFields::CX_ALL,
                          RestoreBondDirOptionTrue);
    CHECK(cxsmi.find("wU:1.0") != std::string::npos);
    cxsmi =
        MolToCXSmiles(*m, ps, SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS);
    CHECK(cxsmi.find("wU:1.0") == std::string::npos);
    m->getBondWithIdx(2)->setBondDir(Bond::BondDir::UNKNOWN);
    cxsmi = MolToCXSmiles(*m, SmilesWriteParams(),
                          SmilesWrite::CXSmilesFields::CX_ALL,
                          RestoreBondDirOptionTrue);
    CHECK(cxsmi.find("wU:1.0") != std::string::npos);
    // wiggly bonds get written even if we don't output coords:
    cxsmi =
        MolToCXSmiles(*m, ps, SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS,
                      RestoreBondDirOptionClear);
    CHECK(cxsmi.find("w:1.0") == std::string::npos);
  }

  SECTION("double bond stereo") {
    auto m_from_ctab = R"CTAB(
  Mrv2211 09152216122D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.125 -12.0417 0 0
M  V30 2 C -5.7913 -11.2717 0 0
M  V30 3 C -4.4576 -12.0417 0 0
M  V30 4 C -8.4587 -11.2717 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3 CFG=2
M  V30 3 1 1 4
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m_from_ctab);
    auto m = "CC=CC |w:2.2|"_smiles;
    REQUIRE(m);
    unsigned int bondcfg = 0;
    CHECK(m->getBondWithIdx(2)->getPropIfPresent("_MolFileBondCfg", bondcfg));
    CHECK(bondcfg == 2);
  }

  SECTION("some invalid cxsmiles") {
    std::vector<std::string> smileses = {"CC(O)Cl |wU:0.0,wU:1.0|",
                                         "CC(O)Cl |wUP:1.0|", "CC(O)Cl |wU:1|"};
    for (const auto &smi : smileses) {
      INFO(smi);
      CHECK_THROWS_AS(SmilesToMol(smi), SmilesParseException);
    }
  }
}

TEST_CASE("ring bond stereochemistry in CXSMILES") {
  UseLegacyStereoPerceptionFixture useLegacy(false);
  SECTION("basic reading") {
    std::vector<std::pair<std::string, Bond::BondStereo>> tests = {
        {"C1CCCC/C=C/CCC1 |t:5|", Bond::BondStereo::STEREOTRANS},
        {"C1CCCCC=CCCC1 |t:5|", Bond::BondStereo::STEREOTRANS},
        {"C1CCCCC=CCCC1 |c:5|", Bond::BondStereo::STEREOCIS},
        {"C1CCCC/C=C/CCC1 |c:5|", Bond::BondStereo::STEREOCIS},
        {"C1CCCCC=CCCC1 |ctu:5|", Bond::BondStereo::STEREOANY},
        {"C1CCCC/C=C/CCC1 |ctu:5|", Bond::BondStereo::STEREOANY},
    };
    for (const auto &[smi, val] : tests) {
      std::unique_ptr<RWMol> m{SmilesToMol(smi)};
      INFO(smi);
      REQUIRE(m);
      CHECK(m->getBondWithIdx(5)->getBondType() == Bond::BondType::DOUBLE);
      CHECK(m->getBondWithIdx(5)->getStereo() == val);
    }
  }
  SECTION("read multiple values") {
    std::unique_ptr<RWMol> m{SmilesToMol("C1CCCCC=CCCC=1 |t:5,9|")};
    REQUIRE(m);
    CHECK(m->getBondWithIdx(5)->getBondType() == Bond::BondType::DOUBLE);
    CHECK(m->getBondWithIdx(5)->getStereo() == Bond::BondStereo::STEREOTRANS);
    CHECK(m->getBondWithIdx(9)->getBondType() == Bond::BondType::DOUBLE);
    CHECK(m->getBondWithIdx(9)->getStereo() == Bond::BondStereo::STEREOTRANS);
  }
  SECTION("skip small rings") {
    std::vector<std::pair<std::string, Bond::BondStereo>> tests = {
        {"C1=CCCCC1 |ctu:0|", Bond::BondStereo::STEREONONE},
        {"C1=CCCCC1 |c:0|", Bond::BondStereo::STEREONONE}};
    for (const auto &[smi, val] : tests) {
      std::unique_ptr<RWMol> m{SmilesToMol(smi)};
      INFO(smi);
      REQUIRE(m);
      CHECK(m->getBondWithIdx(0)->getBondType() == Bond::BondType::DOUBLE);
      CHECK(m->getBondWithIdx(0)->getStereo() == val);
    }
  }
  SECTION("basic writing") {
    std::vector<std::pair<std::string, std::string>> tests = {
        {"C1CCCC/C=C/CCC1 |t:5|", "C1=C/CCCCCCCC/1 |t:0|"},
        {"C1CCCCC=CCCC1 |t:5|", "C1=C/CCCCCCCC/1 |t:0|"},
        {"C1CCCCC=CCCC1 |c:5|", "C1=C\\CCCCCCCC/1 |c:0|"},
        {"C1CCCC/C=C/CCC1 |c:5|", "C1=C\\CCCCCCCC/1 |c:0|"},
        {"C1=CCCCCCCCC1 |ctu:0|", "C1=CCCCCCCCC1 |ctu:0|"},
        {"C=CCCCCCCCC |ctu:0|", "C=CCCCCCCCC"}
        // we don't write the markers for non-ring bonds
    };
    for (const auto &[smi, val] : tests) {
      std::unique_ptr<RWMol> m{SmilesToMol(smi)};
      INFO(smi);
      REQUIRE(m);
      SmilesWriteParams ops;
      auto cxsmi = MolToCXSmiles(
          *m, ops, SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS);
      CHECK(cxsmi == val);
    }
  }
}

TEST_CASE(
    "Github #5722: check w/c/t/ctu CX labels use bond positions from SMILES",
    "[SMILES][bug]") {
  SECTION("'w:' label") {
    auto m = "CC1CN1C=CC1CC1 |w:4.5|"_smiles;
    REQUIRE(m);
    auto b = m->getBondWithIdx(4);
    REQUIRE(b->getBondType() == Bond::BondType::DOUBLE);
    CHECK(b->getBondDir() == Bond::BondDir::UNKNOWN);
  }
  SECTION("'ctu:' label") {
    auto m = "CC1CN1C=CC1CC1 |ctu:5|"_smiles;
    REQUIRE(m);
    auto b = m->getBondWithIdx(4);
    REQUIRE(b->getBondType() == Bond::BondType::DOUBLE);
    CHECK(b->getStereo() == Bond::STEREOANY);
  }
  SECTION("'c:' label") {
    UseLegacyStereoPerceptionFixture useLegacy(false);

    auto m = "CC1CN1C=CC1CC1 |c:5|"_smiles;
    REQUIRE(m);
    auto b = m->getBondWithIdx(4);
    REQUIRE(b->getBondType() == Bond::BondType::DOUBLE);
    CHECK(b->getStereo() == Bond::STEREOCIS);
  }
  SECTION("'t:' label") {
    UseLegacyStereoPerceptionFixture useLegacy(false);

    auto m = "CC1CN1C=CC1CC1 |t:5|"_smiles;
    REQUIRE(m);
    auto b = m->getBondWithIdx(4);
    REQUIRE(b->getBondType() == Bond::BondType::DOUBLE);
    CHECK(b->getStereo() == Bond::STEREOTRANS);
  }
}

TEST_CASE("Github #5683: SMARTS bond ordering should be the same as SMILES") {
  SECTION("basics") {
    auto m = "C(OC)C"_smarts;
    REQUIRE(m);
    CHECK(m->getBondBetweenAtoms(0, 1)->getIdx() == 0);
    CHECK(m->getBondBetweenAtoms(1, 2)->getIdx() == 1);
    CHECK(m->getBondBetweenAtoms(0, 3)->getIdx() == 2);
  }
  SECTION("as reported: SMARTS which should generate an exception") {
    auto sma = "O=C(C=CCC1CCCCC1)N1N=Cc2ccccc2C1c1ccccc1 |w:3.1|";
    CHECK_THROWS_AS(SmartsToMol(sma), SmilesParseException);
  }
}

TEST_CASE("SMARTS for molecule with zero order bond should match itself") {
  auto m = "CN(<-[Li+])C"_smiles;
  m->getBondBetweenAtoms(1, 2)->setBondType(Bond::BondType::ZERO);
  const auto sma = MolToSmarts(*m);
  std::unique_ptr<ROMol> q{SmartsToMol(sma)};
  SubstructMatchParameters ps;
  ps.maxMatches = 1;
  CHECK(SubstructMatch(*m, *q, ps).size() == 1);
}

TEST_CASE("smilesSymbol in SMARTS", "[smarts][smilesSymbol]") {
  // Probably just the first case is going to be useful to people
  SECTION("smilesSymbol without queries") {
    auto m = "CC*"_smiles;
    REQUIRE(m);
    m->getAtomWithIdx(2)->setProp(common_properties::smilesSymbol, "Xa");
    CHECK(MolToSmarts(*m) == "[#6]-[#6]-[Xa]");
  }
  SECTION("smilesSymbol with query on other atoms") {
    auto m = "[#6]cc"_smarts;
    REQUIRE(m);
    m->getAtomWithIdx(0)->setProp(common_properties::smilesSymbol, "Xa");
    CHECK(MolToSmarts(*m) == "[Xa]cc");
  }
  SECTION("smilesSymbol with aromaticity query") {
    auto m = "ccc"_smarts;
    REQUIRE(m);
    m->getAtomWithIdx(0)->setProp(common_properties::smilesSymbol, "Xa");
    CHECK(MolToSmarts(*m) == "[Xa;a]cc");
  }
  SECTION("smilesSymbol with aliphatic query") {
    auto m = "CCC"_smarts;
    REQUIRE(m);
    m->getAtomWithIdx(0)->setProp(common_properties::smilesSymbol, "Xa");
    CHECK(MolToSmarts(*m) == "[Xa;A]CC");
  }
  SECTION("smilesSymbol with additional queries") {
    auto m = "[#6]C[#6]"_smarts;
    REQUIRE(m);
    auto atom = m->getAtomWithIdx(1);

    atom->expandQuery(makeAtomExplicitDegreeQuery(3), Queries::COMPOSITE_AND);
    atom->expandQuery(makeAtomTotalValenceQuery(4), Queries::COMPOSITE_AND);
    atom->expandQuery(makeAtomFormalChargeQuery(2), Queries::COMPOSITE_AND);

    atom->setProp(common_properties::smilesSymbol, "Xa");

    CHECK(MolToSmarts(*m) == "[#6][Xa;A&D3&v4&+2][#6]");
  }
  SECTION("degree query") {
    auto m = "[X3]-C"_smarts;
    REQUIRE(m);
    m->getAtomWithIdx(0)->setProp(common_properties::smilesSymbol, "Xa");
    CHECK(MolToSmarts(*m) == "[Xa;X3]-C");
  }
  SECTION("atom list") {
    auto m = "[C,N,O]C"_smarts;
    REQUIRE(m);
    m->getAtomWithIdx(0)->setProp(common_properties::smilesSymbol, "Xa");
    CHECK(MolToSmarts(*m) == "[Xa;C,N,O]C");
  }
}

TEST_CASE("Atropisomer output in CXSMILES", "[SMILES]") {
  SECTION("'WithChiralAtom' label") {
    auto mol =
        "CC1=C(N2C=CC=C2[C@H](C)Cl)C(C)CCC1 |(2.679,0.4142,;1.3509,1.181,;0.0229,0.4141,;0.0229,-1.1195,;1.2645,-2.0302,;0.7901,-3.4813,;-0.7446,-3.4813,;-1.219,-2.0302,;-2.679,-1.5609,;-3.0039,-0.0556,;-3.8202,-2.595,;-1.3054,1.1809,;-2.6335,0.4141,;-1.3054,2.7145,;0.0229,3.4813,;1.3509,2.7146,),wD:2.11,wU:8.10,&1:8|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getNumConformers() == 1);

    RDKit::SmilesWriteParams ps;
    ps.canonical = false;
    unsigned int flags = SmilesWrite::CXSmilesFields::CX_COORDS |
                         SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES |
                         SmilesWrite::CXSmilesFields::CX_ATOM_PROPS |
                         SmilesWrite::CXSmilesFields::CX_BOND_CFG |
                         SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO;

    auto smi = MolToCXSmiles(*mol, ps, flags,
                             RestoreBondDirOption::RestoreBondDirOptionTrue);

    CHECK(
        smi ==
        "CC1=C(n2cccc2[C@H](C)Cl)C(C)CCC1 |(2.679,0.4142,;1.3509,1.181,;0.0229,0.4141,;0.0229,-1.1195,;1.2645,-2.0302,;0.7901,-3.4813,;-0.7446,-3.4813,;-1.219,-2.0302,;-2.679,-1.5609,;-3.0039,-0.0556,;-3.8202,-2.595,;-1.3054,1.1809,;-2.6335,0.4141,;-1.3054,2.7145,;0.0229,3.4813,;1.3509,2.7146,),wD:2.11,wU:8.10,&1:8|");

    flags = SmilesWrite::CXSmilesFields::CX_COORDS |
            SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES |
            SmilesWrite::CXSmilesFields::CX_ATOM_PROPS |
            SmilesWrite::CXSmilesFields::CX_BOND_ATROPISOMER |
            SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO;

    smi = MolToCXSmiles(*mol, ps, flags,
                        RestoreBondDirOption::RestoreBondDirOptionTrue);

    CHECK(
        smi ==
        "CC1=C(n2cccc2[C@H](C)Cl)C(C)CCC1 |(2.679,0.4142,;1.3509,1.181,;0.0229,0.4141,;0.0229,-1.1195,;1.2645,-2.0302,;0.7901,-3.4813,;-0.7446,-3.4813,;-1.219,-2.0302,;-2.679,-1.5609,;-3.0039,-0.0556,;-3.8202,-2.595,;-1.3054,1.1809,;-2.6335,0.4141,;-1.3054,2.7145,;0.0229,3.4813,;1.3509,2.7146,),wD:2.11,&1:8|");

    flags = SmilesWrite::CXSmilesFields::CX_COORDS |
            SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES |
            SmilesWrite::CXSmilesFields::CX_ATOM_PROPS |
            SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO;

    smi = MolToCXSmiles(*mol, ps, flags,
                        RestoreBondDirOption::RestoreBondDirOptionTrue);

    CHECK(
        smi ==
        "CC1=C(n2cccc2[C@H](C)Cl)C(C)CCC1 |(2.679,0.4142,;1.3509,1.181,;0.0229,0.4141,;0.0229,-1.1195,;1.2645,-2.0302,;0.7901,-3.4813,;-0.7446,-3.4813,;-1.219,-2.0302,;-2.679,-1.5609,;-3.0039,-0.0556,;-3.8202,-2.595,;-1.3054,1.1809,;-2.6335,0.4141,;-1.3054,2.7145,;0.0229,3.4813,;1.3509,2.7146,),&1:8|");
  }
}

TEST_CASE("Dative  bond in cxsmiles double double def", "[bug][cxsmiles]") {
  SECTION("basics") {
    SmilesParserParams smilesParserParams;
    smilesParserParams.sanitize = true;
    smilesParserParams.allowCXSMILES = true;

    std::unique_ptr<RWMol> smilesMol(
        SmilesToMol("C1CCC2=[N]1[Fe](\\[O]=C(\\C)/C=C/C1CCCC1)[N]1=C(CCC1)C2",
                    smilesParserParams));
    RDKit::Chirality::reapplyMolBlockWedging(*smilesMol);
    {
      SmilesWriteParams ps;
      ps.canonical = true;

      std::string smilesOut = MolToSmiles(*smilesMol, ps);

      CHECK(smilesOut ==
            "CC(/C=C/C1CCCC1)=[O]->[Fe]1<-[N]2=C(CC3=[N]->1CCC3)CCC2");
    }
  }
}

TEST_CASE("Fieldname not found in SuperatomSgroup in CXSmiles",
          "[bug][cxsmiles]") {
  SECTION("basics") {
    SmilesParserParams smilesParserParams;
    smilesParserParams.sanitize = true;
    smilesParserParams.allowCXSMILES = true;

    std::unique_ptr<RWMol> smilesMol(
        SmilesToMol("CC |SgD:0:::|", smilesParserParams));
    RDKit::Chirality::reapplyMolBlockWedging(*smilesMol);
    {
      SmilesWriteParams ps;
      ps.canonical = true;

      std::string smilesOut = MolToCXSmiles(*smilesMol, ps);

      CHECK(smilesOut == "CC |SgD:0::::::|");
    }
  }
}
TEST_CASE("ensure unused features are not used") {
  SECTION("isotopes") {
    auto mol1 = "FOCN[15F]"_smiles;
    REQUIRE(mol1);
    auto mol2 = "[15F]OCNF"_smiles;
    REQUIRE(mol2);
    std::vector<unsigned int> ranks;
    SmilesWriteParams ps;
    ps.doIsomericSmiles = true;
    auto smiles = MolToSmiles(*mol1, ps);
    CHECK(smiles == "FOCN[15F]");
    smiles = MolToSmiles(*mol2, ps);
    CHECK(smiles == "FNCO[15F]");

    ps.doIsomericSmiles = false;
    smiles = MolToSmiles(*mol1, ps);
    CHECK(smiles == "FNCOF");
    smiles = MolToSmiles(*mol2, ps);
    CHECK(smiles == "FNCOF");
  }

  SECTION("chirality") {
    auto mol1 = "FC(Cl)OCN[C@H](F)Cl"_smiles;
    REQUIRE(mol1);
    auto mol2 = "FC(Cl)NCO[C@H](F)Cl"_smiles;
    REQUIRE(mol2);
    std::vector<unsigned int> ranks;
    SmilesWriteParams ps;
    ps.doIsomericSmiles = true;
    auto smiles = MolToSmiles(*mol1, ps);
    CHECK(smiles == "FC(Cl)OCN[C@H](F)Cl");
    smiles = MolToSmiles(*mol2, ps);
    CHECK(smiles == "FC(Cl)NCO[C@H](F)Cl");

    ps.doIsomericSmiles = false;
    smiles = MolToSmiles(*mol1, ps);
    CHECK(smiles == "FC(Cl)NCOC(F)Cl");
    smiles = MolToSmiles(*mol2, ps);
    CHECK(smiles == "FC(Cl)NCOC(F)Cl");
  }

  SECTION("ring stereo") {
    auto mol1 = "CC1CCC(CC1)NO[C@@H]1CC[C@H](C)CC1"_smiles;
    REQUIRE(mol1);
    auto mol2 = "CC1CCC(CC1)ON[C@@H]1CC[C@H](C)CC1"_smiles;
    REQUIRE(mol2);
    std::vector<unsigned int> ranks;
    SmilesWriteParams ps;
    ps.doIsomericSmiles = true;
    auto smiles = MolToSmiles(*mol1, ps);
    CHECK(smiles == "CC1CCC(NO[C@H]2CC[C@@H](C)CC2)CC1");
    smiles = MolToSmiles(*mol2, ps);
    CHECK(smiles == "CC1CCC(ON[C@H]2CC[C@@H](C)CC2)CC1");

    ps.doIsomericSmiles = false;
    smiles = MolToSmiles(*mol1, ps);
    CHECK(smiles == "CC1CCC(NOC2CCC(C)CC2)CC1");
    smiles = MolToSmiles(*mol2, ps);
    CHECK(smiles == "CC1CCC(NOC2CCC(C)CC2)CC1");
  }

  SECTION("enhanced stereo") {
    // if we aren't doing CXSMILES then the enhanced stereo shouldn't enter into
    // consideration in canonicalization
    std::unique_ptr<ROMol> mol1 = "F[C@H](Cl)NCO[C@H](F)Cl |&1:6|"_smiles;
    REQUIRE(mol1);
    std::unique_ptr<ROMol> mol2 = "F[C@H](Cl)OCN[C@H](F)Cl |&1:6|"_smiles;
    REQUIRE(mol2);
    SmilesWriteParams ps;
    ps.doIsomericSmiles = true;
    auto smiles = MolToSmiles(*mol1, ps);
    CHECK(smiles == "F[C@H](Cl)NCO[C@H](F)Cl");
    smiles = MolToSmiles(*mol2, ps);
    CHECK(smiles == "F[C@H](Cl)NCO[C@H](F)Cl");

    smiles = MolToCXSmiles(*mol1, ps);
    CHECK(smiles == "F[C@H](Cl)NCO[C@H](F)Cl |&1:6|");
    smiles = MolToCXSmiles(*mol2, ps);
    CHECK(smiles == "F[C@H](Cl)OCN[C@H](F)Cl |&1:6|");

    ps.doIsomericSmiles = false;
    smiles = MolToSmiles(*mol1, ps);
    CHECK(smiles == "FC(Cl)NCOC(F)Cl");
    smiles = MolToSmiles(*mol2, ps);
    CHECK(smiles == "FC(Cl)NCOC(F)Cl");
  }

  SECTION("problematic cases") {
    auto mol1 =
        "[C@H]1CC(C[C@@H](C)[C@@H](C)O)C[C@@H](C)C1 |o1:1,11,o2:5,7|"_smiles;
    REQUIRE(mol1);
    Canon::canonicalizeEnhancedStereo(*mol1);
    auto smiles = MolToSmiles(*mol1);
    CHECK(smiles == "C[C@H]1C[CH]CC(C[C@@H](C)[C@@H](C)O)C1");
  }
}
TEST_CASE("enhanced stereo canonicalized") {
  SECTION("case takes forever") {
    std::unique_ptr<RDKit::ROMol> m1 =
        "C[C@H]1CC[C@]2(NC1)O[C@H]1C[C@H]3[C@H]4CC=C5C[C@H](O[C@H]6O[C@H](CO)[C@H](O[C@H]7O[C@H](C)[C@H](O)[C@H](O)[C@H]7O)[C@H](O)[C@H]6O[C@H]6O[C@H](C)[C@H](O)[C@H](O)[C@H]6O)CC[C@]5(C)[C@H]4CC[C@]3(C)[C@H]1[C@H]2C "
        "|&1:1,&2:4,&3:8,&4:10,&5:11,&6:16,&7:18,&8:20,&9:23,&10:25,&11:27,&12:29,&13:31,&14:33,&15:35,&16:37,&17:39,&18:41,&19:43,&20:45,&21:47,&22:51,&23:53,&24:56,&25:58,&26:59|"_smiles;
    REQUIRE(m1);
    std::unique_ptr<RDKit::ROMol> m2 =
        "C[C@@H]1CC[C@]2(NC1)O[C@H]1C[C@H]3[C@H]4CC=C5C[C@H](O[C@H]6O[C@H](CO)[C@H](O[C@H]7O[C@H](C)[C@H](O)[C@H](O)[C@H]7O)[C@H](O)[C@H]6O[C@H]6O[C@H](C)[C@H](O)[C@H](O)[C@H]6O)CC[C@]5(C)[C@H]4CC[C@]3(C)[C@H]1[C@H]2C "
        "|&1:1,&2:4,&3:8,&4:10,&5:11,&6:16,&7:18,&8:20,&9:23,&10:25,&11:27,&12:29,&13:31,&14:33,&15:35,&16:37,&17:39,&18:41,&19:43,&20:45,&21:47,&22:51,&23:53,&24:56,&25:58,&26:59|"_smiles;
    REQUIRE(m2);
    RDKit::canonicalizeStereoGroups(m1);
    RDKit::canonicalizeStereoGroups(m2);

    CHECK(MolToCXSmiles(*m1) == MolToCXSmiles(*m2));
  }

  SECTION("case takes long but not forever") {
    std::unique_ptr<RDKit::ROMol> m1 =
        "C[C@H]1CC[C@]2(NC1)O[C@H]1C[C@H]3[C@H]4CC=C5C[C@H](O[C@H]6O[C@H](CO)[C@H](O[C@H]7O[C@H](C)[C@H](O)[C@H](O)[C@H]7O)[C@H](O)[C@H]6O[C@H]6O[C@H](C)[C@H](O)[C@H](O)[C@H]6O)CC[C@]5(C)[C@H]4CC[C@]3(C)[C@H]1[C@H]2C "
        "|&1:1,&2:4,&3:8,&4:10,&5:11,&6:16,&7:18,&8:20,&9:23,&10:25,&11:27,&12:29|"_smiles;
    REQUIRE(m1);
    std::unique_ptr<RDKit::ROMol> m2 =
        "C[C@@H]1CC[C@]2(NC1)O[C@H]1C[C@H]3[C@H]4CC=C5C[C@H](O[C@H]6O[C@H](CO)[C@H](O[C@H]7O[C@H](C)[C@H](O)[C@H](O)[C@H]7O)[C@H](O)[C@H]6O[C@H]6O[C@H](C)[C@H](O)[C@H](O)[C@H]6O)CC[C@]5(C)[C@H]4CC[C@]3(C)[C@H]1[C@H]2C "
        "|&1:1,&2:4,&3:8,&4:10,&5:11,&6:16,&7:18,&8:20,&9:23,&10:25,&11:27,&12:29|"_smiles;
    REQUIRE(m2);
    RDKit::canonicalizeStereoGroups(m1);
    RDKit::canonicalizeStereoGroups(m2);

    CHECK(MolToCXSmiles(*m1) == MolToCXSmiles(*m2));
  }

  SECTION("basic") {
    std::unique_ptr<ROMol> mol1 = "F[C@H](Cl)NCO[C@H](F)Cl |&1:6|"_smiles;
    REQUIRE(mol1);
    std::unique_ptr<ROMol> mol2 = "F[C@H](Cl)OCN[C@H](F)Cl |&1:6|"_smiles;
    REQUIRE(mol2);

    SmilesWriteParams ps;
    ps.doIsomericSmiles = true;

    RDKit::canonicalizeStereoGroups(mol1);
    RDKit::canonicalizeStereoGroups(mol2);

    auto smiles = MolToCXSmiles(*mol1, ps);
    CHECK(smiles == "F[C@H](Cl)NCO[C@H](F)Cl |a:1,&1:6|");
    smiles = MolToCXSmiles(*mol2, ps);
    CHECK(smiles == "F[C@H](Cl)OCN[C@H](F)Cl |a:1,&1:6|");
  }

  SECTION("trimethylcyclohexane") {
    std::unique_ptr<ROMol> mol1 =
        "C[C@H]1C[C@@H](C)C[C@@H](C)C1 |o1:1,o2:6,o3:3|"_smiles;
    std::unique_ptr<ROMol> mol2 = "C[C@H]1C[C@@H](C)C[C@@H](C)C1 |o1:1|"_smiles;
    std::unique_ptr<ROMol> mol3 =
        "C[C@@H]1C[C@H](C)C[C@H](C)C1 |o1:1,o2:6,o3:3|"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    REQUIRE(mol3);

    SmilesWriteParams ps;
    RDKit::canonicalizeStereoGroups(mol1);
    auto smiles = MolToCXSmiles(*mol1, ps);
    CHECK(smiles == "C[C@H]1C[C@H](C)C[C@H](C)C1 |a:1,o1:3,6|");

    RDKit::canonicalizeStereoGroups(mol2);
    smiles = MolToCXSmiles(*mol2, ps);
    CHECK(smiles == "C[C@H]1C[C@H](C)C[C@H](C)C1 |a:1,o1:3,6|");

    RDKit::canonicalizeStereoGroups(mol3);
    smiles = MolToCXSmiles(*mol3, ps);
    CHECK(smiles == "C[C@H]1C[C@H](C)C[C@H](C)C1 |a:1,o1:3,6|");
  }

  SECTION("multiFrag test") {
    std::unique_ptr<ROMol> mol1 =
        "C[C@H]1C[C@@H](C)C[C@@H](C)C1.F[C@H](Cl)NCO[C@H](F)Cl |o1:6,o2:3,o3:1,o4:15|"_smiles;
    REQUIRE(mol1);

    SmilesWriteParams ps;
    RDKit::canonicalizeStereoGroups(mol1);
    auto smiles = MolToCXSmiles(*mol1, ps);
    CHECK(
        smiles ==
        "C[C@H]1C[C@H](C)C[C@H](C)C1.F[C@H](Cl)NCO[C@H](F)Cl |a:1,10,o1:3,6,o2:15|");
  }
}

std::unique_ptr<ROMol> getSmartsRootedAtAtom(const ROMol &mol, int root_idx) {
  bool doIsomericSmarts = true;
  auto smarts = MolToSmarts(mol, doIsomericSmarts, root_idx);
  return std::unique_ptr<ROMol>{SmartsToMol(smarts)};
}

TEST_CASE("Test rootedAtAtom argument", "[smarts]") {
  SubstructMatchParameters sssparams;
  sssparams.useChirality = GENERATE(false, true);
  CAPTURE(sssparams.useChirality);

  SECTION("fully substituted chiral center in linear mol") {
    auto mol1 = "C[C@](O)(F)CCCl"_smiles;
    auto mol2 = "C[C@](F)(O)CCCl"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    REQUIRE(mol1->getNumAtoms() == 7);
    REQUIRE(mol2->getNumAtoms() == 7);

    auto root_idx = GENERATE(range(-1, 6));
    CAPTURE(root_idx);
    auto qmol1 = getSmartsRootedAtAtom(*mol1, root_idx);
    auto qmol2 = getSmartsRootedAtAtom(*mol2, root_idx);

    CHECK(SubstructMatch(*mol1, *qmol1, sssparams).size() == 1);
    CHECK(SubstructMatch(*mol2, *qmol2, sssparams).size() == 1);
    CHECK(SubstructMatch(*mol2, *qmol1, sssparams).size() ==
          static_cast<size_t>(static_cast<size_t>(!sssparams.useChirality)));
    CHECK(SubstructMatch(*mol1, *qmol2, sssparams).size() ==
          static_cast<size_t>(static_cast<size_t>(!sssparams.useChirality)));
  }

  SECTION("chiral center w/ implicit H in linear mol") {
    auto mol1 = "C[C@H](F)CCCl"_smiles;
    auto mol2 = "C[C@@H](F)CCCl"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    REQUIRE(mol1->getNumAtoms() == 6);
    REQUIRE(mol2->getNumAtoms() == 6);

    auto root_idx = GENERATE(range(-1, 5));
    CAPTURE(root_idx);
    auto qmol1 = getSmartsRootedAtAtom(*mol1, root_idx);
    auto qmol2 = getSmartsRootedAtAtom(*mol2, root_idx);

    CHECK(SubstructMatch(*mol1, *qmol1, sssparams).size() == 1);
    CHECK(SubstructMatch(*mol2, *qmol2, sssparams).size() == 1);
    CHECK(SubstructMatch(*mol2, *qmol1, sssparams).size() ==
          static_cast<size_t>(static_cast<size_t>(!sssparams.useChirality)));
    CHECK(SubstructMatch(*mol1, *qmol2, sssparams).size() ==
          static_cast<size_t>(static_cast<size_t>(!sssparams.useChirality)));
  }

  SECTION("fully substituted, asymmetric chiral atoms (2) in ring") {
    auto mol1 = "C[C@](F)1CC[C@](N)(F)CC1"_smiles;
    auto mol2 = "F[C@](C)1CC[C@](N)(F)CC1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    REQUIRE(mol1->getNumAtoms() == 10);
    REQUIRE(mol2->getNumAtoms() == 10);

    auto root_idx = GENERATE(range(-1, 9));
    CAPTURE(root_idx);
    auto qmol1 = getSmartsRootedAtAtom(*mol1, root_idx);
    auto qmol2 = getSmartsRootedAtAtom(*mol2, root_idx);

    CHECK(SubstructMatch(*mol1, *qmol1, sssparams).size() == 1);
    CHECK(SubstructMatch(*mol2, *qmol2, sssparams).size() == 1);
    CHECK(SubstructMatch(*mol2, *qmol1, sssparams).size() ==
          static_cast<size_t>(!sssparams.useChirality));
    CHECK(SubstructMatch(*mol1, *qmol2, sssparams).size() ==
          static_cast<size_t>(!sssparams.useChirality));
  }

  SECTION("partially substituted, asymmetric chiral atoms (2) in ring") {
    auto mol1 = "C[C@H]1CC[C@@H](N)CC1"_smiles;
    auto mol2 = "C[C@H]1CC[C@H](N)CC1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    REQUIRE(mol1->getNumAtoms() == 8);
    REQUIRE(mol2->getNumAtoms() == 8);

    auto root_idx = GENERATE(range(-1, 7));
    CAPTURE(root_idx);
    auto qmol1 = getSmartsRootedAtAtom(*mol1, root_idx);
    auto qmol2 = getSmartsRootedAtAtom(*mol2, root_idx);

    CHECK(SubstructMatch(*mol1, *qmol1, sssparams).size() == 1);
    CHECK(SubstructMatch(*mol2, *qmol2, sssparams).size() == 1);
    CHECK(SubstructMatch(*mol2, *qmol1, sssparams).size() ==
          static_cast<size_t>(!sssparams.useChirality));
    CHECK(SubstructMatch(*mol1, *qmol2, sssparams).size() ==
          static_cast<size_t>(!sssparams.useChirality));
  }

  SECTION("Issue #7572: prevent rootAtAtom if more than one fragment exists.") {
    auto mol1 =
        "C(C)Oc1c(Cl)c(C(=O)[O-])c(Cl)c(OCC)c1C(=O)[O-].O=C([O-])C#CC#CC(=O)[O-].[Zn][Zn]"_smiles;
    REQUIRE(mol1);
    auto ps = SmilesWriteParams();
    ps.rootedAtAtom = 20;

    // Fixed by #8328. Multiple fragments are supported now.
    CHECK(
        MolToSmiles(*mol1, ps) ==
        "CCOc1c(Cl)c(C(=O)[O-])c(Cl)c(OCC)c1C(=O)[O-].O=C([O-])C#CC#CC(=O)[O-].[Zn][Zn]");
  }
}

TEST_CASE(
    "Github #5499: STEREOANY bonds lead to non-stable SMILES/SMARTS strings") {
  SECTION("as reported") {
    std::string mb = R"CTAB(7643724
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 14 15 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 2.776307 0.296081 0.000000 0
M  V30 2 N 2.708144 -1.320684 0.000000 0
M  V30 3 N 2.976836 -2.283790 0.000000 0
M  V30 4 N 4.328520 -0.579106 0.000000 0
M  V30 5 C 1.812993 0.027197 0.000000 0
M  V30 6 C 1.770944 -0.971719 0.000000 0
M  V30 7 C 3.329205 -0.537040 0.000000 0
M  V30 8 C 3.124872 1.233298 0.000000 0
M  V30 9 C 0.968792 0.563284 0.000000 0
M  V30 10 C 0.884677 -1.434947 0.000000 0
M  V30 11 C 0.082334 0.100263 0.000000 0
M  V30 12 C 0.040268 -0.899052 0.000000 0
M  V30 13 C 2.487521 2.003850 0.000000 0
M  V30 14 C 2.836086 2.941066 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 5
M  V30 2 1 1 7
M  V30 3 1 1 8
M  V30 4 1 2 3
M  V30 5 1 2 6
M  V30 6 1 2 7
M  V30 7 2 4 7 CFG=2
M  V30 8 2 5 6
M  V30 9 1 5 9
M  V30 10 1 6 10
M  V30 11 1 8 13
M  V30 12 2 9 11
M  V30 13 2 10 12
M  V30 14 1 11 12
M  V30 15 2 13 14
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB";
    auto params = v2::FileParsers::MolFileParserParams();
    params.sanitize = false;
    auto m = v2::FileParsers::MolFromMolBlock(mb, params);
    REQUIRE(m);
    auto osmi = MolToSmiles(*m);
    auto smiparams = v2::SmilesParse::SmilesParserParams();
    smiparams.sanitize = false;
    auto m2 = v2::SmilesParse::MolFromSmiles(osmi, smiparams);
    REQUIRE(m2);
    CHECK(MolToSmiles(*m2) == osmi);
  }
}

TEST_CASE("leaks on unclosed rings") {
  // These are Ok except for the unmatched ring closures
  SECTION("Ok SMARTS") {
    auto m = "C1.C2.C3.C4.C5"_smarts;
    REQUIRE(!m);
  }
  SECTION("Ok SMILES") {
    auto m = "C1.C2.C3.C4.C5"_smiles;
    REQUIRE(!m);
  }
  // These are bogus, but fail parsing AFTER we've seen
  // the unmatched ring closure
  SECTION("Bad SMARTS") {
    auto m = "C1C)foo"_smarts;
    REQUIRE(!m);
  }
  SECTION("Bad SMILES") {
    auto m = "C1C)foo"_smiles;
    REQUIRE(!m);
  }
}

TEST_CASE("Github #7295") {
  SECTION("basics") {
    auto m = "C[C@@H]1CC[C@@H](C(=O)O)CC1.Cl"_smiles;
    REQUIRE(m);
    auto smi1 = MolToSmiles(*m);
    m->clearComputedProps();
    auto smi2 = MolToSmiles(*m);
    CHECK(smi1 == smi2);
    auto m2(*m);
    bool cleanIt = true;
    MolOps::assignStereochemistry(m2, cleanIt);
    auto smi3 = MolToSmiles(m2);
    CHECK(smi1 == smi3);
  }
}

TEST_CASE("simpleSmiles") {
  SECTION("basics") {
    auto m = "CCN(CCO)CCCCC[C@H]1CC[C@H](N(C)C(=O)Oc2ccc(Cl)cc2)CC1.Cl"_smiles;
    REQUIRE(m);
    CHECK(MolToSmiles(*m) ==
          "CCN(CCO)CCCCC[C@H]1CC[C@H](N(C)C(=O)Oc2ccc(Cl)cc2)CC1.Cl");
  }
}

TEST_CASE("CX_BOND_ATROPISOMER option requires ring info", "[bug][cxsmiles]") {
  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase +
      "/Code/GraphMol/FileParsers/test_data/atropisomers/RP-6306_atrop1.sdf";

  auto m = v2::FileParsers::MolFromMolFile(fName);
  REQUIRE(m);

  auto atropBond = m->getBondWithIdx(3);
  REQUIRE(atropBond->getStereo() == Bond::STEREOATROPCW);

  // Clear ring info to check that atropisomer wedging doesn't fail
  // if the info is not present
  bool includeRingInfo = true;
  m->clearComputedProps(includeRingInfo);

  auto ps = SmilesWriteParams();
  auto flags = SmilesWrite::CXSmilesFields::CX_BOND_ATROPISOMER;

  // This will fail if there's no ring information
  auto smi = MolToCXSmiles(*m, ps, flags);
  CHECK(smi == "Cc1cc2c(C(N)=O)c(N)n(-c3c(C)ccc(O)c3C)c2nc1C |wD:10.9|");
}

TEST_CASE("Github #7372: SMILES output option to disable dative bonds") {
  SECTION("basics") {
    auto m = "[NH3]->[Fe]-[NH2]"_smiles;
    REQUIRE(m);
    auto smi = MolToSmiles(*m);
    CHECK(smi == "[NH2][Fe]<-[NH3]");
    SmilesWriteParams ps;
    ps.includeDativeBonds = false;
    auto newSmi = MolToSmiles(*m, ps);
    CHECK(newSmi == "[NH2][Fe][NH3]");
    // ensure that representation round trips:
    auto m2 = v2::SmilesParse::MolFromSmiles(newSmi);
    REQUIRE(m2);
    CHECK(MolToSmiles(*m2) == smi);
  }
  SECTION("SMARTS basics") {
    auto m = "[NH3]->[Fe]-[NH2]"_smiles;
    REQUIRE(m);
    auto smi = MolToSmarts(*m);
    CHECK(smi == "[#7]->[Fe]-[#7]");
    SmilesWriteParams ps;
    ps.includeDativeBonds = false;
    auto newSmi = MolToSmarts(*m, ps);
    CHECK(newSmi == "[#7]-[Fe]-[#7]");
  }
}

void strip_atom_properties(RWMol *molecule) {
  for (auto atom : molecule->atoms()) {
    for (auto property : atom->getPropList(false, false)) {
      atom->clearProp(property);
    }
  }

  // return molecule;
}

TEST_CASE("Remove CX extension from SMILES", "[cxsmiles]") {
  SECTION("basics") {
    std::unique_ptr<RWMol> molecule(RDKit::SmilesToMol(
        "N[C@@H]([O-])C1=[CH:1]C(=[13CH]C(=C1)N(=O)=O)C(N)[O-] |$_AV:;bar;;foo;;;;;;;;;;;$,c:5,7,t:3|",
        0, false));
    REQUIRE(molecule);
    strip_atom_properties(molecule.get());
    std::string stripped_smiles = RDKit::MolToCXSmiles(*molecule);

    CHECK(stripped_smiles ==
          "NC([O-])C1=[13CH]C(N(=O)=O)=CC([C@@H](N)[O-])=C1");
  }
}

TEST_CASE("Canonicalization of meso structures") {
  SECTION("basics") {
    std::vector<std::pair<std::vector<std::string>, std::string>> data = {
        {{"N[C@H]1CC[C@@H](O)CC1", "N[C@@H]1CC[C@H](O)CC1"},
         "N[C@H]1CC[C@@H](O)CC1"},
        {{"C[C@@H](Cl)C[C@H](C)Cl", "Cl[C@H](C)C[C@H](C)Cl",
          "C[C@@H](Cl)C[C@@H](Cl)C", "C[C@H](Cl)C[C@@H](C)Cl"},
         "C[C@H](Cl)C[C@@H](C)Cl"},
    };
    for (const auto &[smileses, expected] : data) {
      for (const auto &smi : smileses) {
        auto m = v2::SmilesParse::MolFromSmiles(smi);
        REQUIRE(m);
        auto osmi = MolToSmiles(*m);
        INFO(smi);
        CHECK(osmi == expected);
      }
    }
  }
}

TEST_CASE("Ignore atom map numbers") {
  SmilesWriteParams params;
  auto m1 = "[NH2:1]c1ccccc1"_smiles;
  CHECK(MolToSmiles(*m1, params) == "c1ccc([NH2:1])cc1");
  params.ignoreAtomMapNumbers = true;
  CHECK(MolToSmiles(*m1, params) == "[NH2:1]c1ccccc1");
  auto m2 = "Nc1ccccc1"_smiles;
  m1->getAtomWithIdx(0)->setAtomMapNum(0);
  CHECK(MolToSmiles(*m1, params) == MolToSmiles(*m2, params));
  CHECK(MolToSmiles(*m1, true, false, -1, true, false, false, false, true) ==
        MolToSmiles(*m2, true, false, -1, true, false, false, false, true));
}

TEST_CASE("Github #7340", "[Reaction][CX][CXSmiles]") {
  SECTION("Test getCXExtensions with a Vector") {
    // Create the MOL_SPTR_VECT to hold the molecular pointers
    const auto mols = {
        "CCO* |$;;;_R1$(0,0,0;1.5,0,0;1.5,1.5,0;0,1.5,0)|"_smiles,
        "C1CCCCC1 |$;label2;$|"_smiles,
        "CC(=O)O |$;label1;$|"_smiles,
        "*-C-* |$star_e;;star_e$,Sg:n:1::ht|"_smiles,
    };

    std::vector<ROMol *> mol_vect;
    mol_vect.reserve(mols.size());
    for (const auto &mol : mols) {
      mol_vect.push_back(mol.get());
    }

    // Write to smiles to populate atom and bond output order properties
    for (const auto &entry : mol_vect) {
      MolToSmiles(*entry);
    }

    std::string cxExt = SmilesWrite::getCXExtensions(
        mol_vect, RDKit::SmilesWrite::CXSmilesFields::CX_ALL);

    CHECK(
        cxExt ==
        "|(0,1.5,;1.5,1.5,;1.5,0,;0,0,;0,0,;0,0,;0,0,;0,0,;0,0,;0,0,;0,0,;0,0,;0,0,;0,0,;0,0,;0,0,;0,0,),$_R1;;;;;label2;;;;;;label1;;;star_e;;star_e$,Sg:n:15::ht:::|");
  }

  SECTION("Expects an error") {
    const auto mols = {
        "CCO* |$;;;_R1$(0,0,0;1.5,0,0;1.5,1.5,0;0,1.5,0)|"_smiles,
        "C1CCCCC1 |$;label2;$|"_smiles,
        "CC(=O)O |$;label1;$|"_smiles,
        "*-C-* |$star_e;;star_e$,Sg:n:1::ht|"_smiles,
    };

    std::vector<ROMol *> mol_vect;
    mol_vect.reserve(mols.size());
    for (const auto &mol : mols) {
      mol_vect.push_back(mol.get());
    }

    CHECK_THROWS_AS(SmilesWrite::getCXExtensions(
                        mol_vect, RDKit::SmilesWrite::CXSmilesFields::CX_ALL),
                    ValueErrorException);
  }
}

TEST_CASE("trimethylcyclohexane") {
  SECTION("Basic") {
    UseLegacyStereoPerceptionFixture useLegacy(false);

    auto smi = "C[C@H]1C[C@@H](C)C[C@@H](C)C1";
    RDKit::v2::SmilesParse::SmilesParserParams smilesParserParams;
    auto m1 = RDKit::v2::SmilesParse::MolFromSmiles(smi, smilesParserParams);
    auto smiOut = RDKit::MolToCXSmiles(*m1);
    CHECK(smiOut == "C[C@@H]1C[C@H](C)C[C@H](C)C1");
  }
  SECTION("WithEnhancedStereo") {
    UseLegacyStereoPerceptionFixture useLegacy(false);

    auto smi = "C[C@H]1C[C@@H](C)C[C@@H](C)C1 |o1:1,o2:6,o3:3|";
    RDKit::v2::SmilesParse::SmilesParserParams smilesParserParams;
    auto m1 = RDKit::v2::SmilesParse::MolFromSmiles(smi, smilesParserParams);
    auto smiOut = RDKit::MolToCXSmiles(*m1);
    CHECK(smiOut == "C[C@H]1C[C@H](C)C[C@H](C)C1 |o1:1,o2:3,o3:6|");
  }
}

TEST_CASE("DnaTestError", "DnaTestError") {
  SECTION("basics") {
    auto smi =
        "[H]OC[C@H]1O[C@@H]2[C@H](O)[C@@H]1OP(=O)(O)OC[C@H]1O[C@@H]3[C@H](O)[C@@H]1OP(=O)(O)OC[C@H]1O[C@@H]4[C@H](O)[C@@H]1OP(=O)(O)OC[C@H]1O[C@@H]5[C@H](O)[C@@H]1OP(=O)(O)OC[C@H]1O[C@@H]6[C@H](O)[C@@H]1OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O[H])N1C=CC7=N8~N9C(=NC%10=C(N=CN%10[C@@H]%10O[C@H](COP(=O)(O)O[C@H]%11[C@@H](O)[C@@H](O[C@@H]%11COP(=O)(O)O[C@H]%11[C@@H](O)[C@@H](O[C@@H]%11COP(=O)(O)O[C@H]%11[C@@H](O)[C@@H](O[C@@H]%11COP(=O)(O)O[C@H]%11[C@@H](O)[C@@H](O[C@@H]%11COP(=O)(O)O[C@H]%11[C@@H](O)[C@@H](O[C@@H]%11CO[H])N%11C=NC%12=C%11N=C%11~O=C%13N2C=CC2=O~NC%12=N%11~N2%13)N2C=NC%11=C2N=C2N~O=C%12N3C=CC3=N%12~N2C%11=O~N3)N2C=CC3=N%11~N%12C(=NC%13=C(N=CN%134)C%12=O~N3)N~O=C2%11)N2C=C(C)C3=O~NC4=N%11~N3C2=O~C%11=NC2=C4N=CN25)N2C=NC3=C2N=C2~O=C4N6C=CC5=O~NC3=N2~N54)[C@@H](O[H])[C@H]%10O)C9=O~N7)N~O=C18 |(4.67218,-35.5914,;5.34598,-36.5843,;6.84288,-36.4764,;7.68418,-37.7162,;7.17218,-39.1261,;8.35478,-40.0487,;9.59778,-39.2091,;10.7257,-39.6187,;9.18338,-37.7674,;10.1024,-36.5843,;17.9626,-38.0598,;18.5623,-39.0992,;18.5623,-37.0204,;25.2033,-36.5843,;26.7002,-36.4764,;27.5415,-37.7162,;27.0295,-39.1261,;28.2121,-40.0487,;29.4551,-39.2091,;30.583,-39.6187,;29.0407,-37.7674,;29.9597,-36.5843,;37.8199,-38.0598,;38.4196,-39.0992,;38.4196,-37.0204,;45.0606,-36.5843,;46.5575,-36.4764,;47.3988,-37.7162,;46.8868,-39.1261,;48.0694,-40.0487,;49.3124,-39.2091,;50.4403,-39.6187,;48.898,-37.7674,;49.817,-36.5843,;57.6772,-38.0598,;58.2769,-39.0992,;58.2769,-37.0204,;64.9179,-36.5843,;66.4148,-36.4764,;67.2561,-37.7162,;66.7441,-39.1261,;67.9267,-40.0487,;69.1697,-39.2091,;70.2976,-39.6187,;68.7553,-37.7674,;69.6743,-36.5843,;77.5345,-38.0598,;78.1342,-39.0992,;78.1342,-37.0204,;84.7752,-36.5843,;86.2721,-36.4764,;87.1134,-37.7162,;86.6014,-39.1261,;87.784,-40.0487,;89.027,-39.2091,;90.1549,-39.6187,;88.6126,-37.7674,;89.5316,-36.5843,;97.3918,-38.0598,;97.9915,-39.0992,;97.9915,-37.0204,;104.632,-36.5843,;106.129,-36.4764,;106.971,-37.7162,;106.459,-39.1261,;107.641,-40.0487,;108.884,-39.2091,;110.012,-39.6187,;108.47,-37.7674,;109.389,-36.5843,;110.578,-36.7473,;106.21,-27.8713,;106.96,-26.5723,;108.46,-26.5724,;109.21,-27.8714,;108.46,-29.1704,;109.007,-16.4468,;107.581,-15.9834,;106.466,-16.9872,;106.778,-18.4544,;108.205,-18.9178,;108.205,-20.417,;106.778,-20.8806,;105.897,-19.6672,;107.641,-10.2628,;106.459,-9.34016,;106.971,-7.93026,;106.129,-6.69046,;104.632,-6.79836,;97.3918,-8.27388,;97.9915,-9.31328,;97.9915,-7.23447,;89.5316,-6.79836,;88.6126,-7.98146,;89.027,-9.42316,;90.1549,-9.83276,;87.784,-10.2628,;86.6014,-9.34016,;87.1134,-7.93026,;86.2721,-6.69046,;84.7752,-6.79836,;77.5345,-8.27388,;78.1342,-9.31328,;78.1342,-7.23447,;69.6743,-6.79836,;68.7553,-7.98146,;69.1697,-9.42316,;70.2976,-9.83276,;67.9267,-10.2628,;66.7441,-9.34016,;67.2561,-7.93026,;66.4148,-6.69046,;64.9179,-6.79836,;57.6772,-8.27388,;58.2769,-9.31328,;58.2769,-7.23447,;49.817,-6.79836,;48.898,-7.98146,;49.3124,-9.42316,;50.4403,-9.83276,;48.0694,-10.2628,;46.8868,-9.34016,;47.3988,-7.93026,;46.5575,-6.69046,;45.0606,-6.79836,;37.8199,-8.27388,;38.4196,-9.31328,;38.4196,-7.23447,;29.9597,-6.79836,;29.0407,-7.98146,;29.4551,-9.42316,;30.583,-9.83276,;28.2121,-10.2628,;27.0295,-9.34016,;27.5415,-7.93026,;26.7002,-6.69046,;25.2033,-6.79836,;17.9626,-8.27388,;18.5623,-9.31328,;18.5623,-7.23447,;10.1024,-6.79836,;9.18338,-7.98146,;9.59778,-9.42316,;10.7257,-9.83276,;8.35478,-10.2628,;7.17218,-9.34016,;7.68418,-7.93026,;6.84288,-6.69046,;5.34598,-6.79836,;4.67218,-5.80546,;6.58934,-19.3588,;7.47114,-20.5722,;8.89774,-20.1086,;8.89724,-18.6094,;7.47074,-18.146,;7.15874,-16.6788,;8.27334,-15.675,;7.07393,-26.0527,;7.67393,-27.0919,;6.92393,-28.3909,;7.67383,-29.69,;9.17383,-29.69,;9.92383,-28.391,;11.1239,-28.391,;11.1532,-17.9763,;10.0118,-17.6056,;9.69994,-16.1384,;9.17383,-27.092,;26.4675,-19.6672,;27.3493,-20.8806,;28.7759,-20.417,;28.7754,-18.9178,;27.3489,-18.4544,;27.0369,-16.9872,;28.1515,-15.9834,;27.9019,-14.8096,;26.9311,-30.2096,;27.5312,-29.1704,;26.7813,-27.8713,;27.5313,-26.5723,;29.0312,-26.5724,;29.7812,-27.8714,;29.0312,-29.1704,;29.5781,-16.4468,;29.89,-17.914,;31.0314,-18.2847,;30.9813,-27.8714,;46.6386,-17.9427,;47.3886,-16.6437,;48.8885,-16.6438,;49.6385,-17.9428,;48.8885,-19.2418,;49.4354,-26.3755,;48.0088,-25.9121,;46.8942,-26.9159,;47.2062,-28.3831,;48.6327,-28.8465,;48.6332,-30.3457,;47.2066,-30.8093,;46.3248,-29.5959,;49.7473,-27.8427,;50.8887,-28.2134,;50.8386,-17.9428,;47.7592,-24.7383,;46.7884,-20.281,;47.3885,-19.2418,;66.3458,-18.2025,;67.0957,-19.5016,;68.5957,-19.5016,;69.1957,-20.5409,;69.3457,-18.2026,;70.5458,-18.2026,;70.7251,-27.9049,;69.5837,-27.5342,;69.2718,-26.067,;68.5957,-16.9036,;67.0958,-16.9035,;66.4958,-15.8643,;67.8452,-25.6036,;66.7306,-26.6074,;67.0426,-28.0746,;68.4691,-28.538,;68.4696,-30.0372,;67.043,-30.5008,;66.1612,-29.2874,;86.0185,-19.3588,;86.9003,-20.5722,;88.3269,-20.1086,;88.3264,-18.6094,;86.8999,-18.146,;86.5879,-16.6788,;87.7025,-15.675,;86.5031,-26.0527,;87.1031,-27.0919,;86.3531,-28.3909,;87.103,-29.69,;88.603,-29.69,;89.353,-28.391,;90.5531,-28.391,;90.5824,-17.9763,;89.441,-17.6056,;89.1291,-16.1384,;88.603,-27.092,;108.47,-7.98146,;109.389,-6.79836,;110.578,-6.96136,;108.884,-9.42316,;110.012,-9.83276,;109.319,-17.914,;110.461,-18.2847,;110.41,-27.8714,;107.331,-14.8096,;106.36,-30.2096,;106.96,-29.1704,),wD:3.2,15.15,27.28,39.41,51.54,63.67,65.76,84.90,86.93,96.103,98.106,108.116,110.119,120.129,122.132,132.142,134.145,144.155,146.158,wU:5.4,6.6,8.9,17.17,18.19,20.22,29.30,30.32,32.35,41.43,42.45,44.48,53.56,54.58,56.61,66.71,68.74,93.99,94.101,105.112,106.114,117.125,118.127,129.138,130.140,141.151,142.153,243.285,246.289,H:76.81,156.169,164.177,166.182,175.194,183.203,185.208,191.215,202.228,204.231,211.241,215.244,217.249,231.268,239.276,241.281,249.293,251.296,a:3,5,6,8,15,17,18,20,27,29,30,32,39,41,42,44,51,53,54,56,63,65,66,68,84,86,93,94,96,98,105,106,108,110,117,118,120,122,129,130,132,134,141,142,144,146,243,246|";
    RDKit::v2::SmilesParse::SmilesParserParams sp;
    sp.sanitize = false;
    sp.removeHs = false;
    sp.strictCXSMILES = false;
    sp.removeHs = false;

    auto mol = MolFromSmiles(smi, sp);
    RDKit::Chirality::removeNonExplicit3DChirality(*(mol.get()));

    CHECK(mol);
    CHECK(mol->getNumAtoms() == 254);
    CHECK(mol->getNumBonds() == 300);

    // do the operations to the mol

    unsigned int operationThatFailed = 0;
    const unsigned int sanitizeOps = RDKit::MolOps::SANITIZE_ALL ^
                                     RDKit::MolOps::SANITIZE_CLEANUP ^
                                     RDKit::MolOps::SANITIZE_PROPERTIES;

    RDKit::MolOps::sanitizeMol(*mol, operationThatFailed, sanitizeOps);

    RDKit::MolOps::RemoveHsParameters removeParams;
    removeParams.removeHydrides = true;
    removeParams.removeInSGroups = true;
    removeParams.removeMapped = true;
    removeParams.removeNonimplicit = true;
    removeParams.removeWithWedgedBond = true;

    RDKit::MolOps::removeHs(*mol, removeParams, false);

    // 2nd sanitize pass

    RDKit::MolOps::sanitizeMol(*mol, operationThatFailed, sanitizeOps);

    std::string expectedOutSmi =
        "Cc1cn2[C@@H]3O[C@@H]4COP(=O)(O)O[C@H]5[C@@H](O)[C@@H]6n7ccc8NO=c9c%10ncn([C@@H]%11O[C@H](COP(=O)(O)O[C@H]%12[C@@H](O)[C@@H]%13n%14ccc%15NO=c%16c%17ncn([C@@H]%18O[C@H](COP(=O)(O)O[C@H]%19[C@@H](O)[C@H](n%20cnc%21c%22NO=c%23ccn([C@@H]%24O[C@H](CO)[C@@H](OP(=O)(O)OC[C@H]%12O%13)[C@H]%24O)c%12=Oc(nc%21%20)n%22-[nH]%23%12)O[C@@H]%19CO)[C@@H](OP(=O)(O)OC[C@H]5O6)[C@H]%18O)c%17nc5NO=c%14n%15-[nH]5%16)[C@@H](OP(=O)(O)OC[C@H]5O[C@@H](n6cnc%12c%13NO=c1[nH]1-n%13c(O=c21)nc%126)[C@H](O)[C@@H]5OP(=O)(O)OC[C@H]1O[C@H]2n5ccc6=ONc%12c%13ncn([C@@H]%14O[C@H](COP(=O)(O)O[C@H]4[C@H]3O)[C@@H](OP(=O)(O)OC[C@H]3O[C@@H](n4cnc%15c%16=ONc%17ccn([C@@H]%18O[C@H](COP(=O)(O)O[C@H]1[C@H]2O)[C@@H](O)[C@H]%18O)c1=ONc(nc%154)[nH]%16-n%171)[C@H](O)[C@@H]3O)[C@H]%14O)c%13nc1O=c5[nH]6-n%121)[C@H]%11O)c%10nc1NO=c7n8-[nH]19 |(69.1957,-20.5409,;68.5957,-19.5016,;67.0957,-19.5016,;66.3458,-18.2025,;67.9267,-10.2628,;66.7441,-9.34016,;67.2561,-7.93026,;66.4148,-6.69046,;64.9179,-6.79836,;57.6772,-8.27388,;58.2769,-9.31328,;58.2769,-7.23447,;49.817,-6.79836,;48.898,-7.98146,;49.3124,-9.42316,;50.4403,-9.83276,;48.0694,-10.2628,;46.6386,-17.9427,;47.3886,-16.6437,;48.8885,-16.6438,;49.6385,-17.9428,;50.8386,-17.9428,;50.8887,-28.2134,;49.7473,-27.8427,;48.6327,-28.8465,;48.6332,-30.3457,;47.2066,-30.8093,;46.3248,-29.5959,;48.0694,-40.0487,;46.8868,-39.1261,;47.3988,-37.7162,;46.5575,-36.4764,;45.0606,-36.5843,;37.8199,-38.0598,;38.4196,-39.0992,;38.4196,-37.0204,;29.9597,-36.5843,;29.0407,-37.7674,;29.4551,-39.2091,;30.583,-39.6187,;28.2121,-40.0487,;26.7813,-27.8713,;27.5313,-26.5723,;29.0312,-26.5724,;29.7812,-27.8714,;30.9813,-27.8714,;31.0314,-18.2847,;29.89,-17.914,;28.7754,-18.9178,;28.7759,-20.417,;27.3493,-20.8806,;26.4675,-19.6672,;28.2121,-10.2628,;27.0295,-9.34016,;27.5415,-7.93026,;26.7002,-6.69046,;25.2033,-6.79836,;17.9626,-8.27388,;18.5623,-9.31328,;18.5623,-7.23447,;10.1024,-6.79836,;9.18338,-7.98146,;9.59778,-9.42316,;10.7257,-9.83276,;8.35478,-10.2628,;6.58934,-19.3588,;7.47114,-20.5722,;8.89774,-20.1086,;8.89724,-18.6094,;10.0118,-17.6056,;11.1532,-17.9763,;11.1239,-28.391,;9.92383,-28.391,;9.17383,-29.69,;7.67383,-29.69,;6.92393,-28.3909,;8.35478,-40.0487,;7.17218,-39.1261,;7.68418,-37.7162,;6.84288,-36.4764,;5.34598,-36.5843,;9.18338,-37.7674,;10.1024,-36.5843,;17.9626,-38.0598,;18.5623,-39.0992,;18.5623,-37.0204,;25.2033,-36.5843,;26.7002,-36.4764,;27.5415,-37.7162,;27.0295,-39.1261,;9.59778,-39.2091,;10.7257,-39.6187,;7.67393,-27.0919,;7.07393,-26.0527,;8.27334,-15.675,;7.15874,-16.6788,;7.47074,-18.146,;9.69994,-16.1384,;9.17383,-27.092,;7.17218,-9.34016,;7.68418,-7.93026,;6.84288,-6.69046,;5.34598,-6.79836,;29.0407,-7.98146,;29.9597,-6.79836,;37.8199,-8.27388,;38.4196,-9.31328,;38.4196,-7.23447,;45.0606,-6.79836,;46.5575,-6.69046,;47.3988,-7.93026,;46.8868,-9.34016,;29.4551,-9.42316,;30.583,-9.83276,;27.3489,-18.4544,;27.0369,-16.9872,;28.1515,-15.9834,;27.9019,-14.8096,;26.9311,-30.2096,;27.5312,-29.1704,;29.0312,-29.1704,;29.5781,-16.4468,;48.898,-37.7674,;49.817,-36.5843,;57.6772,-38.0598,;58.2769,-39.0992,;58.2769,-37.0204,;64.9179,-36.5843,;66.4148,-36.4764,;67.2561,-37.7162,;66.7441,-39.1261,;67.9267,-40.0487,;66.1612,-29.2874,;67.043,-30.5008,;68.4696,-30.0372,;68.4691,-28.538,;69.5837,-27.5342,;70.7251,-27.9049,;70.5458,-18.2026,;69.3457,-18.2026,;68.5957,-16.9036,;69.2718,-26.067,;67.8452,-25.6036,;66.4958,-15.8643,;67.0958,-16.9035,;66.7306,-26.6074,;67.0426,-28.0746,;69.1697,-39.2091,;70.2976,-39.6187,;68.7553,-37.7674,;69.6743,-36.5843,;77.5345,-38.0598,;78.1342,-39.0992,;78.1342,-37.0204,;84.7752,-36.5843,;86.2721,-36.4764,;87.1134,-37.7162,;86.6014,-39.1261,;87.784,-40.0487,;86.3531,-28.3909,;87.103,-29.69,;88.603,-29.69,;89.353,-28.391,;90.5531,-28.391,;90.5824,-17.9763,;89.441,-17.6056,;88.3264,-18.6094,;88.3269,-20.1086,;86.9003,-20.5722,;86.0185,-19.3588,;87.784,-10.2628,;86.6014,-9.34016,;87.1134,-7.93026,;86.2721,-6.69046,;84.7752,-6.79836,;77.5345,-8.27388,;78.1342,-9.31328,;78.1342,-7.23447,;69.6743,-6.79836,;68.7553,-7.98146,;69.1697,-9.42316,;70.2976,-9.83276,;88.6126,-7.98146,;89.5316,-6.79836,;97.3918,-8.27388,;97.9915,-9.31328,;97.9915,-7.23447,;104.632,-6.79836,;106.129,-6.69046,;106.971,-7.93026,;106.459,-9.34016,;107.641,-10.2628,;105.897,-19.6672,;106.778,-20.8806,;108.205,-20.417,;108.205,-18.9178,;109.319,-17.914,;110.461,-18.2847,;110.41,-27.8714,;109.21,-27.8714,;108.46,-26.5724,;106.96,-26.5723,;106.21,-27.8713,;107.641,-40.0487,;106.459,-39.1261,;106.971,-37.7162,;106.129,-36.4764,;104.632,-36.5843,;97.3918,-38.0598,;97.9915,-39.0992,;97.9915,-37.0204,;89.5316,-36.5843,;88.6126,-37.7674,;89.027,-39.2091,;90.1549,-39.6187,;108.47,-37.7674,;109.389,-36.5843,;108.884,-39.2091,;110.012,-39.6187,;106.96,-29.1704,;106.36,-30.2096,;107.331,-14.8096,;107.581,-15.9834,;106.466,-16.9872,;106.778,-18.4544,;109.007,-16.4468,;108.46,-29.1704,;108.884,-9.42316,;110.012,-9.83276,;108.47,-7.98146,;109.389,-6.79836,;89.027,-9.42316,;90.1549,-9.83276,;86.8999,-18.146,;86.5879,-16.6788,;87.7025,-15.675,;86.5031,-26.0527,;87.1031,-27.0919,;88.603,-27.092,;89.1291,-16.1384,;49.3124,-39.2091,;50.4403,-39.6187,;47.2062,-28.3831,;46.8942,-26.9159,;48.0088,-25.9121,;47.7592,-24.7383,;46.7884,-20.281,;47.3885,-19.2418,;48.8885,-19.2418,;49.4354,-26.3755,),wD:4.4,6.6,30.30,52.52,54.54,78.78,88.87,64.106,100.109,110.118,16.121,129.145,156.179,170.194,172.196,189.214,191.217,203.228,205.231,wU:13.12,14.14,28.28,37.36,38.38,61.60,62.62,76.76,81.81,40.90,90.93,103.112,112.124,122.139,131.147,147.170,149.173,158.181,179.202,180.206,182.208,212.237,213.241,215.243,217.246,227.261,229.264,231.267,240.282,H:21.21,46.45,70.70,94.96,97.103,117.130,121.135,138.154,140.158,143.161,164.187,197.223,221.249,225.257,235.272,239.277,246.288,248.293,a:4,6,13,14,16,28,30,37,38,40,52,54,61,62,64,76,78,81,88,90,100,103,110,112,122,129,131,147,149,156,158,170,172,179,180,182,189,191,203,205,212,213,215,217,227,229,231,240|";
    std::string outMolSmi = "";
    RDKit::SmilesWriteParams swp;

    swp.doIsomericSmiles = true;
    swp.doKekule = false;
    swp.canonical = true;
    swp.cleanStereo = true;
    swp.allBondsExplicit = false;
    swp.allHsExplicit = false;
    swp.doRandom = false;
    swp.rootedAtAtom = -1;
    swp.includeDativeBonds = true;
    swp.ignoreAtomMapNumbers = false;

    try {
      outMolSmi =
          MolToCXSmiles(*mol, swp, RDKit::SmilesWrite::CXSmilesFields::CX_ALL,
                        RestoreBondDirOption::RestoreBondDirOptionTrue);

    } catch (const RDKit::KekulizeException &) {
      outMolSmi = "";
    }

    CHECK(expectedOutSmi == outMolSmi);
  }
}

TEST_CASE("Github 8328",
          "MolToSmiles with rootedAtAtom for multiple fragments") {
  SECTION("basics") {
    auto mol = "[C:1][C:2].[N:3]([C:4])=[O:5]"_smiles;
    auto ps = SmilesWriteParams();
    ps.rootedAtAtom = 0;
    CHECK(MolToSmiles(*mol, ps) == "[C:1][C:2].[N:3]([C:4])=[O:5]");
    ps.rootedAtAtom = 1;
    CHECK(MolToSmiles(*mol, ps) == "[C:2][C:1].[N:3]([C:4])=[O:5]");
    ps.rootedAtAtom = 2;
    CHECK(MolToSmiles(*mol, ps) == "[C:1][C:2].[N:3]([C:4])=[O:5]");
    ps.rootedAtAtom = 3;
    CHECK(MolToSmiles(*mol, ps) == "[C:1][C:2].[C:4][N:3]=[O:5]");
    ps.rootedAtAtom = 4;
    CHECK(MolToSmiles(*mol, ps) == "[C:1][C:2].[O:5]=[N:3][C:4]");
    ps.rootedAtAtom = 5;
    CHECK_THROWS_AS(MolToSmiles(*mol, ps), Invar::Invariant);
  }
  SECTION("Compare with and without canonicalization") {
    auto mol =
        "[Al+3].[Na+2].[O-]S(=O)(=O)[O-].CC(=O)OCC.NC(=O)Cc1ccccc1"_smiles;
    auto ps = SmilesWriteParams();

    ps.canonical = true;
    ps.rootedAtAtom = -1;
    CHECK(MolToSmiles(*mol, ps) ==
          "CCOC(C)=O.NC(=O)Cc1ccccc1.O=S(=O)([O-])[O-].[Al+3].[Na+2]");

    ps.canonical = true;
    ps.rootedAtAtom = 10;
    CHECK(MolToSmiles(*mol, ps) ==
          "NC(=O)Cc1ccccc1.O(CC)C(C)=O.O=S(=O)([O-])[O-].[Al+3].[Na+2]");

    ps.canonical = true;
    ps.rootedAtAtom = 21;
    CHECK(MolToSmiles(*mol, ps) ==
          "CCOC(C)=O.O=S(=O)([O-])[O-].[Al+3].[Na+2].c1cccc(CC(N)=O)c1");

    ps.canonical = false;
    ps.rootedAtAtom = 21;
    CHECK(MolToSmiles(*mol, ps) ==
          "[Al+3].[Na+2].[O-]S(=O)(=O)[O-].CC(=O)OCC.c1cccc(CC(N)=O)c1");
  }
}

TEST_CASE("atoms bound to metals should always have Hs specified") {
  SECTION("basics") {
    std::vector<std::pair<std::string, std::string>> smileses = {
        {"Cl[Pt](F)([NH2])[OH]", "[NH2][Pt]([OH])([F])[Cl]"},
        {"Cl[Pt](F)(<-[NH3])[OH]", "[NH3]->[Pt]([OH])([F])[Cl]"},
    };
    for (const auto &[smi, expected] : smileses) {
      auto m = v2::SmilesParse::MolFromSmiles(smi);
      REQUIRE(m);
      auto osmi = MolToSmiles(*m);
      INFO(smi);
      CHECK(osmi == expected);
    }
  }
}

TEST_CASE("ZOB cx smiles extension", "[smiles][cxsmiles]") {
  SECTION("basics") {
    auto m = "CC"_smiles;
    REQUIRE(m);

    auto b = m->getBondWithIdx(0);
    b->setBondType(Bond::ZERO);

    auto smi = MolToCXSmiles(*m);
    REQUIRE(smi == "C~C |Z:0|");

    auto m2 = RDKit::v2::SmilesParse::MolFromSmiles(smi);
    REQUIRE(m2);

    CHECK(m2->getBondWithIdx(0)->getBondType() == Bond::ZERO);
  }
  SECTION("Reverse") {
    constexpr const char *smi = "FB1(F)N2CCCC/C2=N/C2=[NH+]~1CCC=C2 |Z:12|";

    auto p = v2::SmilesParse::SmilesParserParams();
    p.sanitize = false;
    auto m = v2::SmilesParse::MolFromSmiles(smi, p);
    REQUIRE(m);

    auto b = m->getBondWithIdx(15);
    CHECK(b->getBondType() == Bond::BondType::ZERO);
    CHECK(b->getBeginAtom()->getAtomicNum() == 7);
    CHECK(b->getEndAtom()->getAtomicNum() == 5);

    REQUIRE(MolToCXSmiles(*m) == smi);
  }
}

TEST_CASE("github #8471: fail on bad characters in SMILES") {
  SECTION("as reported") {
    v2::SmilesParse::SmilesParserParams sp;
    std::vector<std::string> badSmiles = {
        "CCl₂O",
        "CCl₂OZr"
        "₂CClO",
    };
    for (const auto &smi : badSmiles) {
      auto m = v2::SmilesParse::MolFromSmiles(smi, sp);
      REQUIRE(!m);
      m = v2::SmilesParse::MolFromSmarts(smi);
      REQUIRE(!m);
    }
  }
}

TEST_CASE(
    "github #8654: stereogroups used in canonicalization even if not in output") {
  SECTION("as reported") {
    auto m =
        "O=C(N[C@H]1C[C@@H](C(=O)O)[C@@H]2C[C@H]12)C1CC(=O)N(Cc2ccccn2)C1 |&1:3,5,9,11|"_smiles;
    REQUIRE(m);
    SmilesWriteParams ps;
    unsigned int cxFlags = 0;
    auto smi = MolToCXSmiles(*m, ps, cxFlags);
    CHECK(smi ==
          "O=C(N[C@H]1C[C@@H](C(=O)O)[C@@H]2C[C@@H]21)C1CC(=O)N(Cc2ccccn2)C1");
  }
}

TEST_CASE("bond labels not being used in fragment canonicalization") {
  SECTION("as reported") {
    auto m1 = "CCC"_smiles;
    REQUIRE(m1);
    std::vector<std::string> bondLabels = {"A", "B"};
    std::vector<std::string> bondLabels2 = {"B", "A"};
    std::vector<int> atoms = {0, 1, 2};
    std::vector<int> bonds = {0, 1};
    SmilesWriteParams ps;
    std::vector<std::string> *atomLabels = nullptr;
    auto smi1 =
        MolFragmentToSmiles(*m1, ps, atoms, &bonds, atomLabels, &bondLabels);
    auto smi2 =
        MolFragmentToSmiles(*m1, ps, atoms, &bonds, atomLabels, &bondLabels2);
    CHECK(smi1 == smi2);
  }
}