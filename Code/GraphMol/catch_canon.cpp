//
//  Copyright (C) 2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/StereoGroup.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Canon.h>
#include <GraphMol/test_fixtures.h>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/CanonicalizeStereoGroups.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>

#include <random>
#include <tuple>

using namespace RDKit;

TEST_CASE("chirality and canonicalization") {
  SECTION("basics") {
    auto mol = "F[C@](O)(Cl)C[C@](F)(O)Cl"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool force = true;
    MolOps::assignStereochemistry(*mol, cleanIt, force);
    CIPLabeler::assignCIPLabels(*mol);
    std::string cip;
    CHECK(mol->getAtomWithIdx(1)->getPropIfPresent(common_properties::_CIPCode,
                                                   cip));
    CHECK(cip == "S");
    CHECK(mol->getAtomWithIdx(5)->getPropIfPresent(common_properties::_CIPCode,
                                                   cip));
    CHECK(cip == "R");
    std::vector<unsigned int> ranks;
    bool breakTies = false;
    Canon::rankMolAtoms(*mol, ranks, breakTies);
    CHECK(ranks[1] > ranks[5]);
    mol->getAtomWithIdx(1)->clearProp(common_properties::_CIPCode);
    mol->getAtomWithIdx(5)->clearProp(common_properties::_CIPCode);
    Canon::rankMolAtoms(*mol, ranks, breakTies);
    CHECK(ranks[1] > ranks[5]);
  }
  SECTION("same") {
    auto mol = "F[C@](O)(Cl)C[C@](O)(F)Cl"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool force = true;
    MolOps::assignStereochemistry(*mol, cleanIt, force);
    CIPLabeler::assignCIPLabels(*mol);

    std::string cip;
    CHECK(mol->getAtomWithIdx(1)->getPropIfPresent(common_properties::_CIPCode,
                                                   cip));
    CHECK(cip == "S");
    CHECK(mol->getAtomWithIdx(5)->getPropIfPresent(common_properties::_CIPCode,
                                                   cip));
    CHECK(cip == "S");
    std::vector<unsigned int> ranks;
    bool breakTies = false;
    Canon::rankMolAtoms(*mol, ranks, breakTies);
    CHECK(ranks[1] == ranks[5]);
    mol->getAtomWithIdx(1)->clearProp(common_properties::_CIPCode);
    mol->getAtomWithIdx(5)->clearProp(common_properties::_CIPCode);
    Canon::rankMolAtoms(*mol, ranks, breakTies);
    CHECK(ranks[1] == ranks[5]);
  }
  SECTION("dependent") {
    auto mol = "F[C@](O)(Cl)[C@](F)(O)[C@](F)(O)Cl"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool force = true;
    MolOps::assignStereochemistry(*mol, cleanIt, force);
    CIPLabeler::assignCIPLabels(*mol);

    std::string cip;
    CHECK(mol->getAtomWithIdx(1)->getPropIfPresent(common_properties::_CIPCode,
                                                   cip));
    CHECK(cip == "S");
    CHECK(mol->getAtomWithIdx(7)->getPropIfPresent(common_properties::_CIPCode,
                                                   cip));
    CHECK(cip == "R");
    CHECK(mol->getAtomWithIdx(4)->getPropIfPresent(common_properties::_CIPCode,
                                                   cip));
    CHECK(cip == "r");
    std::vector<unsigned int> ranks;
    bool breakTies = false;
    Canon::rankMolAtoms(*mol, ranks, breakTies);
    CHECK(ranks[1] > ranks[7]);
    mol->getAtomWithIdx(1)->clearProp(common_properties::_CIPCode);
    mol->getAtomWithIdx(4)->clearProp(common_properties::_CIPCode);
    mol->getAtomWithIdx(7)->clearProp(common_properties::_CIPCode);
    Canon::rankMolAtoms(*mol, ranks, breakTies);
    CHECK(ranks[1] > ranks[7]);
  }
  SECTION("dependent non-chiral") {
    auto mol = "F[C@](O)(Cl)[C@](F)(O)[C@](O)(F)Cl"_smiles;
    REQUIRE(mol);
    bool cleanIt = false;
    bool force = true;
    mol->getAtomWithIdx(4)->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    MolOps::assignStereochemistry(*mol, cleanIt, force);
    CIPLabeler::assignCIPLabels(*mol);

    std::string cip;
    CHECK(mol->getAtomWithIdx(1)->getPropIfPresent(common_properties::_CIPCode,
                                                   cip));
    CHECK(cip == "S");
    CHECK(mol->getAtomWithIdx(7)->getPropIfPresent(common_properties::_CIPCode,
                                                   cip));
    CHECK(cip == "S");
    std::vector<unsigned int> ranks;
    bool breakTies = false;
    Canon::rankMolAtoms(*mol, ranks, breakTies);
    CHECK(ranks[1] == ranks[7]);
    mol->getAtomWithIdx(1)->clearProp(common_properties::_CIPCode);
    mol->getAtomWithIdx(7)->clearProp(common_properties::_CIPCode);
    Canon::rankMolAtoms(*mol, ranks, breakTies);
    CHECK(ranks[1] == ranks[7]);
  }

  SECTION("swap parity") {
    auto mol = "F[C@](O)(Cl)C[C@@](O)(Cl)F"_smiles;
    REQUIRE(mol);
    bool cleanIt = false;
    bool force = true;
    MolOps::assignStereochemistry(*mol, cleanIt, force);
    CIPLabeler::assignCIPLabels(*mol);

    std::string cip;
    CHECK(mol->getAtomWithIdx(1)->getPropIfPresent(common_properties::_CIPCode,
                                                   cip));
    CHECK(cip == "S");
    CHECK(mol->getAtomWithIdx(5)->getPropIfPresent(common_properties::_CIPCode,
                                                   cip));
    CHECK(cip == "S");
    std::vector<unsigned int> ranks;
    bool breakTies = false;
    Canon::rankMolAtoms(*mol, ranks, breakTies);
    CHECK(ranks[1] == ranks[5]);
    mol->getAtomWithIdx(1)->clearProp(common_properties::_CIPCode);
    mol->getAtomWithIdx(5)->clearProp(common_properties::_CIPCode);
    Canon::rankMolAtoms(*mol, ranks, breakTies);
    CHECK(ranks[1] == ranks[5]);
  }
}
TEST_CASE("double bond stereo and canonicalization") {
  SECTION("basics") {
    auto mol = "CC=C(F)C(B)C(F)=CC"_smiles;
    REQUIRE(mol);
    mol->getBondWithIdx(1)->setStereoAtoms(0, 4);
    mol->getBondWithIdx(1)->setStereo(Bond::BondStereo::STEREOTRANS);
    mol->getBondWithIdx(7)->setStereoAtoms(4, 9);
    mol->getBondWithIdx(7)->setStereo(Bond::BondStereo::STEREOCIS);
    bool breakTies = false;
    std::vector<unsigned int> ranks;
    Canon::rankMolAtoms(*mol, ranks, breakTies);
    mol->getBondWithIdx(1)->setStereo(Bond::BondStereo::STEREOCIS);
    mol->getBondWithIdx(7)->setStereo(Bond::BondStereo::STEREOTRANS);
    std::vector<unsigned int> ranks2;
    Canon::rankMolAtoms(*mol, ranks2, breakTies);
    CHECK(ranks[0] == ranks2[9]);
    CHECK(ranks[1] == ranks2[8]);
    CHECK(ranks[2] == ranks2[6]);
    CHECK(ranks[3] == ranks2[7]);
    CHECK(ranks[4] == ranks2[4]);
    CHECK(ranks[5] == ranks2[5]);

    // same as previous example, different controlling atoms
    mol->getBondWithIdx(7)->setStereoAtoms(7, 9);
    mol->getBondWithIdx(7)->setStereo(Bond::BondStereo::STEREOCIS);
    std::vector<unsigned int> ranks3;
    Canon::rankMolAtoms(*mol, ranks3, breakTies);
    CHECK(ranks2 == ranks3);
  }
  SECTION("STEREOANY is higher priority than STEREONONE") {
    auto mol = "CC=C(F)C(B)C(F)=CC"_smiles;
    REQUIRE(mol);
    mol->getBondWithIdx(7)->setStereoAtoms(4, 9);
    mol->getBondWithIdx(7)->setStereo(Bond::BondStereo::STEREOANY);
    bool breakTies = false;
    std::vector<unsigned int> ranks;
    Canon::rankMolAtoms(*mol, ranks, breakTies);
    CHECK(ranks[0] < ranks[9]);
  }
}

TEST_CASE("enhanced stereo canonicalization") {
  SECTION("simple chiral tags") {
    std::vector<std::pair<std::string, std::string>> tests = {
        {"C[C@H](F)Cl |&1:1|", "C[C@@H](F)Cl |&1:1|"},
        {"C[C@H](F)Cl |o1:1|", "C[C@@H](F)Cl |o1:1|"},
    };
    for (const auto &[smi1, smi2] : tests) {
      INFO(smi1 + " : " + smi2);
      std::unique_ptr<RWMol> mol1{SmilesToMol(smi1)};
      REQUIRE(mol1);
      std::unique_ptr<RWMol> mol2{SmilesToMol(smi2)};
      REQUIRE(mol2);

      Canon::canonicalizeEnhancedStereo(*mol1);
      Canon::canonicalizeEnhancedStereo(*mol2);

      CHECK(mol1->getAtomWithIdx(1)->getChiralTag() ==
            mol2->getAtomWithIdx(1)->getChiralTag());
    }
  }
  SECTION("abs groups are not modified") {
    std::vector<std::pair<std::string, std::string>> tests = {
        {"C[C@H](F)Cl |a:1|", "C[C@@H](F)Cl |a:1|"},
    };
    for (const auto &[smi1, smi2] : tests) {
      INFO(smi1 + " : " + smi2);
      std::unique_ptr<RWMol> mol1{SmilesToMol(smi1)};
      REQUIRE(mol1);
      std::unique_ptr<RWMol> mol2{SmilesToMol(smi2)};
      REQUIRE(mol2);

      Canon::canonicalizeEnhancedStereo(*mol1);
      Canon::canonicalizeEnhancedStereo(*mol2);

      CHECK(mol1->getAtomWithIdx(1)->getChiralTag() !=
            mol2->getAtomWithIdx(1)->getChiralTag());
    }
  }
  SECTION("relative chiral tags") {
    std::vector<std::pair<std::string, std::string>> tests = {
        {"C[C@H](F)[C@H](Br)O |&1:1,3|", "C[C@@H](F)[C@@H](Br)O |&1:1,3|"},
        {"C[C@H](F)[C@@H](Br)O |&1:1,3|", "C[C@@H](F)[C@H](Br)O |&1:1,3|"},
        {"O[C@H](Br)[C@H](F)C |&1:1,3|", "O[C@@H](Br)[C@@H](F)C |&1:1,3|"},
        {"O[C@H](Br)[C@@H](F)C |&1:1,3|", "O[C@@H](Br)[C@H](F)C |&1:1,3|"},
    };
    for (const auto &[smi1, smi2] : tests) {
      INFO(smi1 + " : " + smi2);
      std::unique_ptr<RWMol> mol1{SmilesToMol(smi1)};
      REQUIRE(mol1);
      std::unique_ptr<RWMol> mol2{SmilesToMol(smi2)};
      REQUIRE(mol2);

      Canon::canonicalizeEnhancedStereo(*mol1);
      Canon::canonicalizeEnhancedStereo(*mol2);

      CHECK(mol1->getAtomWithIdx(1)->getChiralTag() ==
            mol2->getAtomWithIdx(1)->getChiralTag());
      CHECK(mol1->getAtomWithIdx(3)->getChiralTag() ==
            mol2->getAtomWithIdx(3)->getChiralTag());
    }
  }
  SECTION("multiple groups") {
    std::vector<std::pair<std::string, std::string>> tests = {
        {"C[C@H](F)[C@H](Br)O |&1:1,o2:3|",
         "C[C@@H](F)[C@@H](Br)O |&1:1,o2:3|"},
        {"C[C@H](F)[C@H](Br)O |&1:1,o2:3|",
         "C[C@@H](F)[C@@H](Br)O |o1:3,&1:1|"},
    };
    for (const auto &[smi1, smi2] : tests) {
      INFO(smi1 + " : " + smi2);
      std::unique_ptr<RWMol> mol1{SmilesToMol(smi1)};
      REQUIRE(mol1);
      std::unique_ptr<RWMol> mol2{SmilesToMol(smi2)};
      REQUIRE(mol2);

      Canon::canonicalizeEnhancedStereo(*mol1);
      Canon::canonicalizeEnhancedStereo(*mol2);

      CHECK(mol1->getAtomWithIdx(1)->getChiralTag() ==
            mol2->getAtomWithIdx(1)->getChiralTag());
      CHECK(mol1->getAtomWithIdx(3)->getChiralTag() ==
            mol2->getAtomWithIdx(3)->getChiralTag());
    }
  }
}

TEST_CASE("chiralRandomTest") {
  SECTION("chiralRandomTest") {
    std::vector<std::string> tests = {
        "CNC(=O)[C@@H](NC(=O)[C@H](OCc1cc(F)ccc1F)[C@@H](O)[C@@H](O)[C@H](OCc1cc(F)ccc1F)C(=O)N[C@H](C(=O)NC)C(C)C)C(C)C |o1:19,23,&1:8,37|",
        "CNC(=O)C(NC(=O)[C@H](OCc1ccccc1)[C@@H](O)[C@@H](O)[C@H](OCc1ccccc1)C(=O)NC(C(=O)NC)[C@@H](C)O)[C@H](C)O o1:38,41,&1:17,21|",
        "CCC(=O)NC[C@H]1CN(c2cc(F)c(N3C[C@H]4[C@H](N)[C@H]4C3)c(F)c2)C(=O)O1 |o1:16,17,&1:6,19|",
        "Cc1cc(O)c([C@@]2(C)CC[C@@H]3C[C@@]32C)cc1-c1cc([C@]2(C)CC[C@H]3C[C@@]32C)c(O)cc1C |o1:6,23,&1:19,25|",
        "C[N+]12CC[C@@]34c5ccccc5N5[C@@H]6OCC=C7C[N+]8(C)CC[C@]9%10c%11ccccc%11N([C@@H]%11OCC=C(C1)[C@H](C[C@@H]32)[C@@H]%11[C@H]54)[C@H]9[C@H]6[C@H]7C[C@@H]%108.[Cl-].[Cl-] |o1:42,45,&1:22,41|",
        "Cl.Oc1ccc2c3c1O[C@H]1c4[nH]c5c(c4C[C@@]4(O)[C@@H](C2)N(CC2CC2)CC[C@]314)C[C@]1(O)[C@H]2Cc3ccc(O)c4c3[C@]1(CCN2CC1CC1)[C@H]5O4 |o1:29,40,&1:16,48|",
        "O=C1NC(=O)c2c1c1c3ccccc3n3c1c1c2c2ccccc2n1[C@@H]1[C@H](O)[C@H](O)[C@H](O)[C@H]3N1Cc1ccccc1 |o1:26,32,&1:25,30|",

        "CC1=C(N2C=CC=C2[C@H](C)Cl)C(C)CCC1 |(2.679,0.4142,;1.3509,1.181,;0.0229,0.4141,;0.0229,-1.1195,;1.2645,-2.0302,;0.7901,-3.4813,;-0.7446,-3.4813,;-1.219,-2.0302,;-2.679,-1.5609,;-3.0039,-0.0556,;-3.8202,-2.595,;-1.3054,1.1809,;-2.6335,0.4141,;-1.3054,2.7145,;0.0229,3.4813,;1.3509,2.7146,),wD:2.11,wU:8.10,&1:8|",
        "O=C(CCCc1ccccc1)OC[C@H]1C[C@@H]2O[C@H]1[C@@H]1[C@H]2C(=O)OC1=O |a:18,19,o1:13,15,17|",
        "CC1=CC[C@H]2[C@H](C1)c1c(O)cc(C=C3C4CC5CC(C4)CC3C5)cc1OC2(C)C |&1:4,5|",
        "O=C(O)C(F)(F)F.O=C(O)[C@@]1(N(CCn2cnc3ncncc32)S(=O)(=O)c2ccc(-c3ccc(Cl)cc3)cc2)C[C@H]1c1ccccc1 |o1:10,40|",
        "CCCCNC1N[C@H](CO)[C@@H](O)[C@@H](O)[C@H](CO)N1 |a:12,&2:10,&1:7,14|",
        "CCCCNC1N[C@H](CO)[C@@H](O)[C@@H](O)[C@H](CO)N1 |a:10,12,&1:7,14|",
        "CCCCNC1N[C@H](CO)C(O)[C@@H](O)[C@H](CO)N1 |a:12,&1:7,14|",

        "CCNC(=O)c1ccc(/C(=C2/C[C@H]3CC[C@@H](C2)N3CCc2ccccc2)c2ccccc2)cc1 |&1:12,15|",
        "C[C@H](O)[C@@H](C)[C@H](C)[C@H](C)O |&1:1,&2:3,5,&3:7|",

        "C[C@@H](Cl)C[C@H](C)Cl |a:1,4,|",
        "C[C@H](Cl)C[C@@H](C)Cl |a:1,4,|",

        "N[C@H]1CC[C@@H](O)CC1 |a:1,4|",
        "N[C@H]1CC[C@@H](O)CC1 |a:1,4|",
        "N[C@H]1CC[C@@H](O)CC1 |a:1,4|",
        "N[C@H]1CC[C@@H](O)CC1 |a:1,4|",
        "N[C@H]1CC[C@@H](O)CC1 |a:1,4|",

        "C[C@H]1CC[C@H](CC1)C1CC(CC(C1)[C@@H]1CC[C@H](C)CC1)[C@@H]1CC[C@H](C)CC1 |o1:23,o2:20,o3:13,o4:16,o5: ,&6:1|",

        "C[C@H]1CC[C@H](CC1)C1CC(CC(C1)[C@@H]1CC[C@H](C)CC1)[C@@H]1CC[C@H](C)CC1 |o1:23,o2:20,o3:13,o4:16,o5:4,o6:1|",

        "C[C@H]1CC[C@@H](C[C@@H]2CC[C@H](C)CC2)CC1 |a:9,o1:6,&1:4,&2:1|",

        "C[C@H]1CC[C@@H](C[C@@H]2CC[C@H](C)CC2)CC1 |o1:6,o2:1,9,o3:4|",

        "C[C@H]1CC[C@@H](C[C@@H]2CC[C@H](C)CC2)CC1 |o1:6,o2:1,9,o3:4|",

        "N[C@H]1CC[C@@H](O)CC1",
        "N[C@H]1CC[C@@H](O)CC1 |o2:1,4|",
        "N[C@H]1CC[C@@H](O)CC1 |&2:1,4|",
        "N[C@H]1CC[C@@H](O)CC1 |a:1,4|",

        // no enhanced stereo
        "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1",

        // enhance stereo abs,abs,abs
        "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |a:3,11,14|",

        // abs, abs, or and abs, or abs

        "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |a:3,11,o1:14|",
        "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |a:3,11,o1:14|",
        "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |a:3,11,o1:14|",
        "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |a:3,11,&1:14|",
        "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |a:3,11,&1:14|",
        "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |a:3,11,&1:14|",

    };

    UseLegacyStereoPerceptionFixture reset_stereo_perception{false};

    for (const auto &smi1 : tests) {
      INFO(smi1);

      SmilesParserParams rps;
      rps.sanitize = true;
      rps.removeHs = false;

      std::unique_ptr<ROMol> mol1{SmilesToMol(smi1, rps)};
      REQUIRE(mol1);

      std::vector<unsigned int> idxV(mol1->getNumAtoms());
      for (unsigned int i = 0; i < mol1->getNumAtoms(); ++i) {
        idxV[i] = i;
      }
      std::string smiExpected = "";
      auto randomGen = std::mt19937(0xf00d);
      for (auto i = 0; i < 100; ++i) {
        std::vector<unsigned int> nVect(idxV);
        if (i > 0) {
          std::shuffle(nVect.begin(), nVect.end(), randomGen);
        }
        std::unique_ptr<ROMol> nmol{MolOps::renumberAtoms(*mol1, nVect)};
        RDKit::canonicalizeStereoGroups(nmol);

        auto outSmi1 = MolToCXSmiles(*nmol);

        if (smiExpected == "") {
          smiExpected = outSmi1;
        } else {
          CHECK(outSmi1 == smiExpected);
        }
      }
    }
  }
}

TEST_CASE("pseudoTestCanonFailure") {
  SECTION("canonFailure") {
    UseLegacyStereoPerceptionFixture reset_stereo_perception{false};

    std::string smi1 = "C[C@H]1CC[C@@H](C[C@@H]2CC[C@H](C)CC2)CC1 |&1:4,&2:1|";

    SmilesParserParams rps;
    rps.sanitize = true;
    rps.removeHs = false;

    std::unique_ptr<ROMol> mol1{SmilesToMol(smi1, rps)};
    REQUIRE(mol1);

    std::vector<unsigned int> idxV(mol1->getNumAtoms());
    for (unsigned int i = 0; i < mol1->getNumAtoms(); ++i) {
      idxV[i] = i;
    }

    auto randomGen = std::mt19937(0xf00d);
    std::vector<unsigned int> nVect(idxV);
    std::shuffle(nVect.begin(), nVect.end(), randomGen);

    std::unique_ptr<ROMol> nmol{MolOps::renumberAtoms(*mol1, nVect)};

    RDKit::canonicalizeStereoGroups(mol1);
    RDKit::canonicalizeStereoGroups(nmol);

    SmilesWriteParams wp;
    wp.canonical = false;
    auto outSmi1Check = MolToCXSmiles(*mol1, wp);
    auto outSmi1Check2 = MolToCXSmiles(*nmol, wp);
    CHECK(outSmi1Check == outSmi1Check2);

    wp.canonical = true;

    outSmi1Check = MolToCXSmiles(*mol1, wp);
    outSmi1Check2 = MolToCXSmiles(*nmol, wp);
    CHECK(outSmi1Check == outSmi1Check2);
  }
}
TEST_CASE("pseudoTest1") {
  SECTION("pseudoTest1") {
    std::vector<std::tuple<std::string, std::string, std::string>> tests = {

        {"CC1=C(N2C=CC=C2[C@H](C)Cl)C(C)CCC1 |(2.679,0.4142,;1.3509,1.181,;0.0229,0.4141,;0.0229,-1.1195,;1.2645,-2.0302,;0.7901,-3.4813,;-0.7446,-3.4813,;-1.219,-2.0302,;-2.679,-1.5609,;-3.0039,-0.0556,;-3.8202,-2.595,;-1.3054,1.1809,;-2.6335,0.4141,;-1.3054,2.7145,;0.0229,3.4813,;1.3509,2.7146,),wD:2.11,wU:8.10,&1:8|",
         "CC1=C(N2C=CC=C2[C@H](C)Cl)C(C)CCC1 |(2.679,0.4142,;1.3509,1.181,;0.0229,0.4141,;0.0229,-1.1195,;1.2645,-2.0302,;0.7901,-3.4813,;-0.7446,-3.4813,;-1.219,-2.0302,;-2.679,-1.5609,;-3.0039,-0.0556,;-3.8202,-2.595,;-1.3054,1.1809,;-2.6335,0.4141,;-1.3054,2.7145,;0.0229,3.4813,;1.3509,2.7146,),wD:2.11,wU:8.10,a:2,&1:8|",
         "CC1=C(n2cccc2[C@H](C)Cl)C(C)CCC1 |(2.679,0.4142,;1.3509,1.181,;0.0229,0.4141,;0.0229,-1.1195,;1.2645,-2.0302,;0.7901,-3.4813,;-0.7446,-3.4813,;-1.219,-2.0302,;-2.679,-1.5609,;-3.0039,-0.0556,;-3.8202,-2.595,;-1.3054,1.1809,;-2.6335,0.4141,;-1.3054,2.7145,;0.0229,3.4813,;1.3509,2.7146,),wD:8.9,2.11,a:2,&1:8|"},
        {"O=C(CCCc1ccccc1)OC[C@H]1C[C@@H]2O[C@H]1[C@@H]1[C@H]2C(=O)OC1=O |a:18,19,o1:13,15,17|",
         "O=C(CCCc1ccccc1)OC[C@@H]1C[C@H]2O[C@@H]1[C@@H]1[C@H]2C(=O)OC1=O |a:18,19,o1:13,15,17|",
         "O=C(CCCc1ccccc1)OC[C@@H]1C[C@H]2O[C@@H]1[C@H]1C(=O)OC(=O)[C@H]12 |a:18,24,o1:13,15,17|"},
        {"CC1=CC[C@H]2[C@H](C1)c1c(O)cc(C=C3C4CC5CC(C4)CC3C5)cc1OC2(C)C |&1:4,5|",
         "CC1=CC[C@@H]2[C@@H](C1)c1c(O)cc(C=C3C4CC5CC(C4)CC3C5)cc1OC2(C)C |&1:4,5|",
         "CC1=CC[C@H]2[C@H](C1)c1c(O)cc(C=C3C4CC5CC(C4)CC3C5)cc1OC2(C)C |&1:4,5|"},
        {"O=C(O)C(F)(F)F.O=C(O)[C@@]1(N(CCn2cnc3ncncc32)S(=O)(=O)c2ccc(-c3ccc(Cl)cc3)cc2)C[C@H]1c1ccccc1 |o1:10,40|",
         "O=C(O)[C@@]1(N(CCn2cnc3ncncc32)S(=O)(=O)c2ccc(-c3ccc(Cl)cc3)cc2)C[C@H]1c1ccccc1.O=C(O)C(F)(F)F |o1:3,33|",
         "O=C(O)C(F)(F)F.O=C(O)[C@@]1(N(CCn2cnc3ncncc32)S(=O)(=O)c2ccc(-c3ccc(Cl)cc3)cc2)C[C@H]1c1ccccc1 |o1:10,40|"},
        {"CCCCNC1N[C@H](CO)[C@@H](O)[C@@H](O)[C@H](CO)N1 |a:12,&2:10,&1:7,14|",
         "CCCCNC1N[C@@H](CO)[C@H](O)[C@H](O)[C@@H](CO)N1 |a:7,14,&1:10,&2:12|",
         "CCCCNC1N[C@@H](CO)[C@H](O)[C@H](O)[C@@H](CO)N1 |a:7,14,&1:10,&2:12|"},
        {"CCCCNC1N[C@H](CO)[C@@H](O)[C@@H](O)[C@H](CO)N1 |a:10,12,&1:7,14|",
         "CCCCNC1N[C@@H](CO)[C@@H](O)[C@@H](O)[C@@H](CO)N1 |a:10,12,&1:7,14|",
         "CCCCNC1N[C@@H](CO)[C@H](O)[C@H](O)[C@@H](CO)N1 |a:7,14,&1:10,12|"},
        {"CCCCNC1N[C@H](CO)C(O)[C@@H](O)[C@H](CO)N1 |a:12,&1:7,14|",
         "CCCCNC1N[C@@H](CO)C(O)[C@@H](O)[C@@H](CO)N1 |a:12,&1:7,14|",
         "CCCCNC1N[C@H](CO)C(O)[C@@H](O)[C@H](CO)N1 |a:12,&1:7,14|"},

        {"CCNC(=O)c1ccc(/C(=C2/C[C@H]3CC[C@@H](C2)N3CCc2ccccc2)c2ccccc2)cc1 |&1:12,15|",
         "CCNC(=O)c1ccc(\\C(=C2/C[C@@H]3CC[C@H](C2)N3CCc2ccccc2)c2ccccc2)cc1 |&1:12,15|",
         "CCNC(=O)c1ccc(/C(=C2/C[C@H]3CC[C@@H](C2)N3CCc2ccccc2)c2ccccc2)cc1"},
        {"C[C@H](O)[C@@H](C)[C@H](C)[C@H](C)O |&1:1,&2:3,5,&3:7|",
         "C[C@H](O)[C@H](C)[C@@H](C)[C@H](C)O |&1:1,&2:3,5,&3:7|",
         "C[C@H](O)[C@H](C)[C@@H](C)[C@H](C)O |&1:1,&2:3,5,&3:7|"},

        {"C[C@@H](Cl)C[C@H](C)Cl |a:1,4,|", "C[C@@H](Cl)C[C@H](C)Cl |o1:1,4,|",
         "C[C@H](Cl)C[C@@H](C)Cl"},
        {"C[C@H](Cl)C[C@@H](C)Cl |a:1,4,|", "C[C@@H](Cl)C[C@H](C)Cl |&1:1,4,|",
         "C[C@H](Cl)C[C@@H](C)Cl"},

        {"N[C@H]1CC[C@@H](O)CC1 |a:1,4|", "N[C@H]1CC[C@@H](O)CC1 |o1:1,4|",
         "N[C@H]1CC[C@@H](O)CC1"},
        {"N[C@H]1CC[C@@H](O)CC1 |a:1,4|", "N[C@H]1CC[C@@H](O)CC1 |&1:1,4|",
         "N[C@H]1CC[C@@H](O)CC1"},
        {"N[C@H]1CC[C@@H](O)CC1 |a:1,4|", "N[C@@H]1CC[C@H](O)CC1 |a:1,4|",
         "N[C@H]1CC[C@@H](O)CC1"},
        {"N[C@H]1CC[C@@H](O)CC1 |a:1,4|", "N[C@@H]1CC[C@H](O)CC1 |o1:1,4|",
         "N[C@H]1CC[C@@H](O)CC1"},
        {"N[C@H]1CC[C@@H](O)CC1 |a:1,4|", "N[C@@H]1CC[C@H](O)CC1 |&1:1,4|",
         "N[C@H]1CC[C@@H](O)CC1"},

        // {"C[C@H]1C[C@@H](C)C[C@@H](C)C1 |o1:1,o2:6,o3:3|",
        //  "C[C@@H]1C[C@H](C)C[C@@H](C)C1 |a:3,o1:6,o3:1|",
        //  "CC1CC(C)CC(C)C1"},  // this excercises a bug in the new canon -
        //  these
        //                       // are ring stereo atoms but the chiral
        //                       markers
        //                       // are removed!

        {"C[C@H]1CC[C@H](CC1)C1CC(CC(C1)[C@@H]1CC[C@H](C)CC1)[C@@H]1CC[C@H](C)CC1 |o1:23,o2:20,o3:13,o4:16,o5:4,&6:1|",
         "C[C@H]1CC[C@@H](CC1)C1CC(CC(C1)[C@H]1CC[C@H](C)CC1)[C@H]1CC[C@H](C)CC1 |a:1,13,20,o2:4,o6:16,&4:23|",
         "C[C@H]1CC[C@@H](C2CC([C@H]3CC[C@@H](C)CC3)CC([C@H]3CC[C@@H](C)CC3)C2)CC1 |a:4,8,17,o1:1,o2:11,&1:20|"},

        {"C[C@H]1CC[C@H](CC1)C1CC(CC(C1)[C@@H]1CC[C@H](C)CC1)[C@@H]1CC[C@H](C)CC1 |o1:23,o2:20,o3:13,o4:16,o5:4,o6:1|",
         "C[C@H]1CC[C@@H](CC1)C1CC(CC(C1)[C@H]1CC[C@H](C)CC1)[C@H]1CC[C@H](C)CC1 |a:4,13,23,o2:20,o4:16,o6:1|",
         "C[C@H]1CC[C@@H](C2CC([C@H]3CC[C@@H](C)CC3)CC([C@H]3CC[C@@H](C)CC3)C2)CC1 |a:4,8,17,o1:1,o2:11,o3:20|"},

        {"C[C@H]1CC[C@@H](C[C@@H]2CC[C@H](C)CC2)CC1 |a:9,o1:6,&1:4,&2:1|",
         "C[C@H]1CC[C@H](C[C@H]2CC[C@H](C)CC2)CC1 |a:4,9,&1:6,o2:1|",
         "C[C@H]1CC[C@@H](C[C@H]2CC[C@@H](C)CC2)CC1 |a:4,6,o1:1,&1:9|"},

        {"C[C@H]1CC[C@@H](C[C@@H]2CC[C@H](C)CC2)CC1 |o1:6,o2:1,9,o3:4|",
         "C[C@H]1CC[C@H](C[C@H]2CC[C@H](C)CC2)CC1 |o1:6,o2:1,9,o3:4|",
         "C[C@H]1CC[C@@H](C[C@H]2CC[C@@H](C)CC2)CC1 |a:4,6,o1:1,o2:9|"},

        {"C[C@H]1CC[C@@H](C[C@@H]2CC[C@H](C)CC2)CC1 |o1:6,o2:1,9,o3:4|",
         "C[C@H]1CC[C@H](C[C@H]2CC[C@H](C)CC2)CC1 |a:4,9,o1:6,o2:1|",
         "C[C@H]1CC[C@@H](C[C@H]2CC[C@@H](C)CC2)CC1 |a:4,6,o1:1,o2:9|"},

        {"N[C@H]1CC[C@@H](O)CC1", "N[C@@H]1CC[C@H](O)CC1",
         "N[C@H]1CC[C@@H](O)CC1"},
        {"N[C@H]1CC[C@@H](O)CC1 |o2:1,4|", "N[C@@H]1CC[C@H](O)CC1 |o2:1,4|",
         "N[C@H]1CC[C@@H](O)CC1"},
        {"N[C@H]1CC[C@@H](O)CC1 |&2:1,4|", "N[C@@H]1CC[C@H](O)CC1 |&2:1,4|",
         "N[C@H]1CC[C@@H](O)CC1"},
        {"N[C@H]1CC[C@@H](O)CC1 |a:1,4|", "N[C@@H]1CC[C@H](O)CC1 |a:1,4|",
         "N[C@H]1CC[C@@H](O)CC1"},

        // no enhanced stereo
        {"CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1",
         "CC(C)[C@H]1CCCCN1C(=O)[C@@H]1CC[C@H](C)CC1",
         "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1"},

        // enhance stereo abs,abs,abs
        {"CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |a:3,11,14|",
         "CC(C)[C@H]1CCCCN1C(=O)[C@@H]1CC[C@H](C)CC1 |a:3,11,14|",
         "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1"},

        // abs, abs, or and abs, or abs

        {"CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |a:3,11,o1:14|",
         "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@H](C)CC1 |a:3,11,o1:14|",
         "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@H](C)CC1 |a:3,11,o1:14|"},
        {"CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |a:3,11,o1:14|",
         "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@H](C)CC1 |a:3,14,o1:11|",
         "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@H](C)CC1 |a:3,11,o1:14|"},
        {"CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |a:3,11,o1:14|",
         "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |a:3,o1:11,o2:14|",
         "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@H](C)CC1 |a:3,11,o1:14|"},
        {"CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |a:3,11,&1:14|",
         "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@H](C)CC1 |a:3,11,&1:14|",
         "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@H](C)CC1 |a:3,11,&1:14|"},
        {"CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |a:3,11,&1:14|",
         "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@H](C)CC1 |a:3,14,&1:11|",
         "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@H](C)CC1 |a:3,11,&1:14|"},
        {"CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |a:3,11,&1:14|",
         "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |a:3,&1:11,&2:14|",
         "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@H](C)CC1 |a:3,11,&1:14|"},

    };

    UseLegacyStereoPerceptionFixture reset_stereo_perception{false};

    for (const auto &[smi1, smi2, smiExpected] : tests) {
      INFO(smi1 + " : " + smi2 + " : " + smiExpected);

      SmilesParserParams rps;
      rps.sanitize = true;
      rps.removeHs = false;

      std::unique_ptr<ROMol> mol1{SmilesToMol(smi1, rps)};
      REQUIRE(mol1);
      std::unique_ptr<ROMol> mol2{SmilesToMol(smi2, rps)};
      REQUIRE(mol2);

      std::vector<unsigned int> idxV(mol1->getNumAtoms());
      for (unsigned int i = 0; i < mol1->getNumAtoms(); ++i) {
        idxV[i] = i;
      }

      RDKit::canonicalizeStereoGroups(mol2);
      auto outSmi2 = MolToCXSmiles(*mol2);
      REQUIRE(outSmi2 == smiExpected);
      auto randomGen = std::mt19937(0xf00d);
      for (auto i = 0; i < 100; ++i) {
        INFO("i: " << i);
        std::vector<unsigned int> nVect(idxV);
        std::shuffle(nVect.begin(), nVect.end(), randomGen);

        std::unique_ptr<ROMol> nmol{MolOps::renumberAtoms(*mol1, nVect)};

        RDKit::canonicalizeStereoGroups(nmol);

        auto outSmi1 = MolToCXSmiles(*nmol);

        CHECK(outSmi1 == outSmi2);
        CHECK(outSmi1 == smiExpected);
      }
    }
  }

  SECTION("pseudoTestFails") {
    std::vector<std::string> tests = {
        "CNC(=O)[C@@H](NC(=O)[C@H](OCc1cc(F)ccc1F)[C@@H](O)[C@@H](O)[C@H](OCc1cc(F)ccc1F)C(=O)N[C@H](C(=O)NC)C(C)C)C(C)C |o1:19,23,&1:8,37|",
        "CNC(=O)C(NC(=O)[C@H](OCc1ccccc1)[C@@H](O)[C@@H](O)[C@H](OCc1ccccc1)C(=O)NC(C(=O)NC)[C@@H](C)O)[C@H](C)O |o1:38,41,&1:17,21|",
        "CCC(=O)NC[C@H]1CN(c2cc(F)c(N3C[C@H]4[C@H](N)[C@H]4C3)c(F)c2)C(=O)O1 |o1:16,17,&1:6,19|",
        "Cc1cc(O)c([C@@]2(C)CC[C@@H]3C[C@@]32C)cc1-c1cc([C@]2(C)CC[C@H]3C[C@@]32C)c(O)cc1C |o1:6,23,&1:19,25|",
        "C[N+]12CC[C@@]34c5ccccc5N5[C@@H]6OCC=C7C[N+]8(C)CC[C@]9%10c%11ccccc%11N([C@@H]%11OCC=C(C1)[C@H](C[C@@H]32)[C@@H]%11[C@H]54)[C@H]9[C@H]6[C@H]7C[C@@H]%108.[Cl-].[Cl-] |o1:42,45,&1:22,41|",
        "Cl.Oc1ccc2c3c1O[C@H]1c4[nH]c5c(c4C[C@@]4(O)[C@@H](C2)N(CC2CC2)CC[C@]314)C[C@]1(O)[C@H]2Cc3ccc(O)c4c3[C@]1(CCN2CC1CC1)[C@H]5O4 |o1:29,40,&1:16,48|",
        "O=C1NC(=O)c2c1c1c3ccccc3n3c1c1c2c2ccccc2n1[C@@H]1[C@H](O)[C@H](O)[C@H](O)[C@H]3N1Cc1ccccc1 |o1:26,32,&1:25,30|",
    };

    UseLegacyStereoPerceptionFixture reset_stereo_perception{false};

    for (const auto &smi1 : tests) {
      INFO(smi1);

      SmilesParserParams rps;
      rps.sanitize = true;
      rps.removeHs = false;

      std::unique_ptr<ROMol> mol1{SmilesToMol(smi1, rps)};
      REQUIRE(mol1);

      std::vector<unsigned int> idxV(mol1->getNumAtoms());
      for (unsigned int i = 0; i < mol1->getNumAtoms(); ++i) {
        idxV[i] = i;
      }

      for (auto i = 0; i < 100; ++i) {
        std::vector<unsigned int> nVect(idxV);
        std::shuffle(nVect.begin(), nVect.end(), std::mt19937(0xf00d));
        std::unique_ptr<ROMol> nmol{MolOps::renumberAtoms(*mol1, nVect)};

        RDKit::canonicalizeStereoGroups(nmol);

        auto outSmi1 = MolToCXSmiles(*nmol);

        CHECK(outSmi1 == smi1);
      }
    }
  }
}
TEST_CASE("using enhanced stereo in rankMolAtoms") {
  SECTION("basics: different ranks") {
    std::vector<std::string> smis{
        "C[C@H](F)[C@@H](F)C |a:1|",
        "C[C@H](F)[C@@H](F)C |&1:1|"
        "C[C@H](F)[C@@H](F)C |o1:1|",
        "C[C@H](F)[C@@H](F)C |o1:1,a:3|",
        "C[C@H](F)[C@@H](F)C |o1:1,&:3|",
        "C[C@H](F)[C@@H](F)C |a:1,&2:3|",
    };
    for (auto &smi : smis) {
      INFO(smi);
      std::unique_ptr<RWMol> mol{SmilesToMol(smi)};
      REQUIRE(mol);
      bool breakTies = false;
      std::vector<unsigned int> atomRanks;
      Canon::rankMolAtoms(*mol, atomRanks, breakTies);
      CHECK(atomRanks[1] != atomRanks[3]);
    }
  }
  SECTION("basics: same ranks") {
    std::vector<std::string> smis{
        "C[C@H](F)[C@@H](F)C |o1:1,o2:3|",
        "C[C@H](F)[C@@H](F)C |&1:1,&2:3|",
        "C[C@H](F)[C@@H](F)C |a:1,a:3|",
    };
    for (auto &smi : smis) {
      INFO(smi);
      std::unique_ptr<RWMol> mol{SmilesToMol(smi)};
      REQUIRE(mol);
      bool breakTies = false;
      std::vector<unsigned int> atomRanks;
      Canon::rankMolAtoms(*mol, atomRanks, breakTies);
      CHECK(atomRanks[1] == atomRanks[3]);
    }
  }
  SECTION("more complex, include group membership") {
    auto m1 =
        "C[C@@H]1CC(C[C@@H](C)[C@@H](C)O)C[C@H](C)C1 |o1:1,7,o2:5,11|"_smiles;
    REQUIRE(m1);
    std::vector<unsigned int> atomRanks;
    bool breakTies = false;
    Canon::rankMolAtoms(*m1, atomRanks, breakTies);
    CHECK(atomRanks[11] > atomRanks[1]);

    auto m2 =
        "C[C@H](O)[C@H](C)CC1C[C@H](C)C[C@@H](C)C1 |o1:1,8;o2:3,11|"_smiles;
    REQUIRE(m2);
    Canon::rankMolAtoms(*m2, atomRanks, breakTies);
    CHECK(atomRanks[11] > atomRanks[8]);
  }
}

TEST_CASE("more enhanced stereo canonicalization") {
  // FIX: add tests for ring stereo in an s group
  SECTION("case 1") {
    auto m1 =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)F |a:3,&1:1,7,&2:5,r|"_smiles;
    REQUIRE(m1);
    auto m2 = "C[C@H](O)[C@H](C)[C@@H](C)[C@H](C)F |a:3,&1:1,7,&2:5,r|"_smiles;
    REQUIRE(m2);
    Canon::canonicalizeEnhancedStereo(*m1);
    Canon::canonicalizeEnhancedStereo(*m2);
    CHECK(MolToCXSmiles(*m1) == MolToCXSmiles(*m2));
  }
  SECTION("case 2") {
    auto m1 =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&3:3,5,o1:7,&2:1,r|"_smiles;
    REQUIRE(m1);
    auto m2 =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&3:3,5,&2:7,o1:1,r|"_smiles;
    REQUIRE(m2);
    Canon::canonicalizeEnhancedStereo(*m1);
    Canon::canonicalizeEnhancedStereo(*m2);
    CHECK(MolToCXSmiles(*m1) == MolToCXSmiles(*m2));
  }
  SECTION("case 3") {
    auto m1 =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&8:3,5,o1:7,&7:1,r|"_smiles;
    REQUIRE(m1);
    auto m2 =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&3:3,5,&2:7,o1:1,r|"_smiles;
    REQUIRE(m2);
    Canon::canonicalizeEnhancedStereo(*m1);
    Canon::canonicalizeEnhancedStereo(*m2);
    CHECK(MolToCXSmiles(*m1) == MolToCXSmiles(*m2));
  }
  SECTION("case 4") {
    auto m1 =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&8:3,5,o1:7,&7:1,r|"_smiles;
    REQUIRE(m1);
    auto m2 =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&3:3,5,&2:7,o1:1,r|"_smiles;
    REQUIRE(m2);

    CHECK(MolToCXSmiles(*m1) == MolToCXSmiles(*m2));
  }
  SECTION("case 5") {
    auto m1 =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&8:3,5,o1:7,&7:1,r|"_smiles;
    REQUIRE(m1);
    auto m2 =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&3:3,5,&2:7,o1:1,r|"_smiles;
    REQUIRE(m2);

    forwardStereoGroupIds(*m1);
    forwardStereoGroupIds(*m2);

    auto cx1 = MolToCXSmiles(*m1);
    auto cx2 = MolToCXSmiles(*m2);
    CHECK(cx1 == "C[C@H](O)[C@H](C)[C@@H](C)[C@H](C)O |o1:1,&7:7,&8:3,5|");
    CHECK(cx2 == "C[C@H](O)[C@H](C)[C@@H](C)[C@H](C)O |o1:1,&2:7,&3:3,5|");
  }

  SECTION("case 5a") {
    std::unique_ptr<ROMol> m1 =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&8:3,5,o1:7,&7:1,r|"_smiles;
    REQUIRE(m1);
    std::unique_ptr<ROMol> m2 =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&3:3,5,&2:7,o1:1,r|"_smiles;
    REQUIRE(m2);

    forwardStereoGroupIds(*m1);
    forwardStereoGroupIds(*m2);

    RDKit::canonicalizeStereoGroups(m1);

    RDKit::canonicalizeStereoGroups(m2);

    auto cx1 = MolToCXSmiles(*m1);
    auto cx2 = MolToCXSmiles(*m2);
    CHECK(cx1 == "C[C@H](O)[C@H](C)[C@@H](C)[C@H](C)O |o1:1,&1:3,5,&2:7|");
    CHECK(cx2 == "C[C@H](O)[C@H](C)[C@@H](C)[C@H](C)O |o1:1,&1:3,5,&2:7|");
  }
  SECTION("case 6") {
    auto m1 =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&8:3,5,o1:7,&7:1,r|"_smiles;
    REQUIRE(m1);
    auto m2 =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&3:3,5,&2:7,o1:1,r|"_smiles;
    REQUIRE(m2);

    forwardStereoGroupIds(*m1);
    forwardStereoGroupIds(*m2);

    // Canonicalization resets the Stereo Group IDs
    Canon::canonicalizeEnhancedStereo(*m1);
    Canon::canonicalizeEnhancedStereo(*m2);
    CHECK(MolToCXSmiles(*m1) == MolToCXSmiles(*m2));

    // "read" ids are also reset!
    forwardStereoGroupIds(*m1);
    forwardStereoGroupIds(*m2);

    CHECK(MolToCXSmiles(*m1) == MolToCXSmiles(*m2));
  }
}

TEST_CASE("ensure unused features are not used") {
  SECTION("isotopes") {
    auto mol = "[13CH3]OC"_smiles;
    REQUIRE(mol);
    std::vector<unsigned int> ranks;
    bool breakTies = false;
    bool includeChirality = true;
    bool includeIsotopes = true;
    Canon::rankMolAtoms(*mol, ranks, breakTies, includeChirality,
                        includeIsotopes);
    CHECK(ranks[0] != ranks[2]);

    includeIsotopes = false;
    Canon::rankMolAtoms(*mol, ranks, breakTies, includeChirality,
                        includeIsotopes);
    CHECK(ranks[0] == ranks[2]);
  }
  SECTION("chirality") {
    auto mol = "F[C@H](Cl)OC(F)Cl"_smiles;
    REQUIRE(mol);
    std::vector<unsigned int> ranks;
    bool breakTies = false;
    bool includeChirality = true;
    bool includeIsotopes = true;
    Canon::rankMolAtoms(*mol, ranks, breakTies, includeChirality,
                        includeIsotopes);
    CHECK(ranks[1] != ranks[4]);

    includeChirality = false;
    Canon::rankMolAtoms(*mol, ranks, breakTies, includeChirality,
                        includeIsotopes);
    CHECK(ranks[1] == ranks[4]);
  }
  SECTION("chirality and stereogroups") {
    auto mol = "F[C@H](Cl)O[C@H](F)Cl |o1:1|"_smiles;
    REQUIRE(mol);
    std::vector<unsigned int> ranks;
    bool breakTies = false;
    bool includeChirality = true;
    bool includeIsotopes = true;
    Canon::rankMolAtoms(*mol, ranks, breakTies, includeChirality,
                        includeIsotopes);
    CHECK(ranks[1] != ranks[4]);

    includeChirality = false;
    Canon::rankMolAtoms(*mol, ranks, breakTies, includeChirality,
                        includeIsotopes);
    CHECK(ranks[1] == ranks[4]);
  }
}

TEST_CASE(
    "GitHub Issue #6633: Pre-condition violation in canonicalization of dative bond adjacent to double bond",
    "[bug][canonicalization]") {
  auto mb = R"CTAB(
                    3D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 16 16 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C  -2.0033 -1.4133 -0.0473 0
M  V30 2 C  -2.9101 -0.3985 -0.2677 0
M  V30 3 O  -2.7092 0.8645 -0.2504 0
M  V30 4 Ir -0.9429 1.8106 0.2184 0
M  V30 5 N  0.0151 -0.0816 0.3618 0
M  V30 6 C  1.4929 -0.0477 0.5631 0
M  V30 7 C  -0.6236 -1.2309 0.2291 0
M  V30 8 C  -4.3730 -0.7437 -0.5877 0
M  V30 9 H  -2.3752 -2.4232 -0.1048 0
M  V30 10 H  1.8628 -0.9806 0.9803 0
M  V30 11 H  1.6928 0.7152 1.3165 0
M  V30 12 H  2.0044 0.1878 -0.3701 0
M  V30 13 H  -4.9409 0.1756 -0.7308 0
M  V30 14 H  -4.4149 -1.3416 -1.4982 0
M  V30 15 H  -4.8022 -1.3104 0.2386 0
M  V30 16 H  0.0202 -2.0891 0.3538 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 7
M  V30 2 2 1 2
M  V30 3 1 1 9
M  V30 4 1 2 3
M  V30 5 1 2 8
M  V30 6 1 3 4
M  V30 7 9 5 4
M  V30 8 1 5 6
M  V30 9 2 5 7
M  V30 10 1 6 10
M  V30 11 1 6 11
M  V30 12 1 6 12
M  V30 13 1 7 16
M  V30 14 1 8 13
M  V30 15 1 8 14
M  V30 16 1 8 15
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB";

  auto countStereoBonds = [](const auto &mol) {
    unsigned num_stereo_bonds = 0;
    for (const auto bond : mol.bonds()) {
      if (bond->getBondType() == Bond::BondType::DOUBLE &&
          bond->getStereo() != Bond::BondStereo::STEREONONE) {
        ++num_stereo_bonds;
      }
    }
    return num_stereo_bonds;
  };

  auto sanitize = true;
  auto removeHs = false;
  std::unique_ptr<ROMol> mol(MolBlockToMol(mb, sanitize, removeHs));

  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 16);
  REQUIRE(countStereoBonds(*mol) == 2);

  CHECK_NOTHROW(MolToSmiles(*mol));

  CHECK(countStereoBonds(*mol) == 2);
}

TEST_CASE("atom mapping in canonicalization") {
  SECTION("basics") {
    auto m = "[F:1]C([F:2])O"_smiles;
    REQUIRE(m);
    std::vector<unsigned int> ranks;
    bool breakTies = false;
    bool includeChirality = true;
    bool includeIsotopes = true;
    bool includeAtomMaps = true;
    Canon::rankMolAtoms(*m, ranks, breakTies, includeChirality, includeIsotopes,
                        includeAtomMaps);
    CHECK(ranks[0] != ranks[2]);
    includeAtomMaps = false;
    Canon::rankMolAtoms(*m, ranks, breakTies, includeChirality, includeIsotopes,
                        includeAtomMaps);
    CHECK(ranks[0] == ranks[2]);
  }
}

TEST_CASE(
    "GitHub Issue #7023: \"Inconsistent state\" when manually sanitizing and assigning stereo when using the new stereo algorithm",
    "[bug]") {
  UseLegacyStereoPerceptionFixture reset_stereo_perception{false};

  const auto molb = R"CTAB("
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 40 47 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C -1.412000 -2.520800 0.000000 0
M  V30 2 N -2.236600 -2.495700 0.000000 0
M  V30 3 N -3.801500 -2.612400 0.000000 0
M  V30 4 C -4.626200 -2.587400 0.000000 0
M  V30 5 C -4.904800 -3.363900 0.000000 0
M  V30 6 C -5.696800 -3.595000 0.000000 0
M  V30 7 C -4.966800 -0.807100 0.000000 0
M  V30 8 N -4.155000 -0.660100 0.000000 0
M  V30 9 C -4.044000 0.157400 0.000000 0
M  V30 10 C -3.039200 0.506000 0.000000 0
M  V30 11 C -2.237800 0.310100 0.000000 0
M  V30 12 N -2.176500 -0.512600 0.000000 0
M  V30 13 C -1.375100 -0.708600 0.000000 0
M  V30 14 C -0.941100 -0.006900 0.000000 0
M  V30 15 C -0.118400 0.054400 0.000000 0
M  V30 16 C -1.474300 0.622600 0.000000 0
M  V30 17 C -0.910300 1.224700 0.000000 0
M  V30 18 C -1.697100 1.417000 0.000000 0
M  V30 19 C -5.077900 -1.624600 0.000000 0
M  V30 20 C -5.893400 -1.749400 0.000000 0
M  V30 21 C -1.133400 -1.744200 0.000000 0
M  V30 22 C -0.309800 -1.791700 0.000000 0
M  V30 23 C -2.515200 -3.272300 0.000000 0
M  V30 24 C -3.570500 -3.404400 0.000000 0
M  V30 25 C -3.278900 -4.176100 0.000000 0
M  V30 26 C 0.225100 -3.058100 0.000000 0
M  V30 27 C -0.404400 -3.591300 0.000000 0
M  V30 28 C -1.181000 -3.312700 0.000000 0
M  V30 29 C -1.076700 -4.131100 0.000000 0
M  V30 30 Co -0.633100 -1.614100 0.000000 0 CHG=1 VAL=6
M  V30 31 C -4.787100 0.515600 0.000000 0
M  V30 32 C -4.934100 1.327400 0.000000 0
M  V30 33 C -1.862800 -3.777200 0.000000 0
M  V30 34 C -1.887800 -4.601800 0.000000 0
M  V30 35 C -4.252300 -3.868900 0.000000 0
M  V30 36 C -4.678900 -4.575000 0.000000 0
M  V30 37 C -3.869300 -4.599600 0.000000 0
M  V30 38 C -5.357500 -0.080500 0.000000 0
M  V30 39 C -6.124200 -0.385000 0.000000 0
M  V30 40 C -6.015200 0.417600 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 21 1
M  V30 2 1 1 2
M  V30 3 1 28 1
M  V30 4 1 23 2 CFG=3
M  V30 5 1 2 30
M  V30 6 1 23 33
M  V30 7 1 23 24
M  V30 8 1 33 28
M  V30 9 1 24 3
M  V30 10 1 24 35
M  V30 11 2 3 4
M  V30 12 1 5 4
M  V30 13 1 4 19
M  V30 14 1 5 6 CFG=3
M  V30 15 1 35 5
M  V30 16 2 19 7
M  V30 17 1 7 8
M  V30 18 1 38 7
M  V30 19 2 8 9
M  V30 20 1 9 10
M  V30 21 1 31 9
M  V30 22 2 10 11
M  V30 23 1 11 12
M  V30 24 1 11 16
M  V30 25 2 12 13
M  V30 26 1 14 13
M  V30 27 1 14 15 CFG=3
M  V30 28 1 14 16
M  V30 29 1 16 17
M  V30 30 1 16 18
M  V30 31 1 31 38
M  V30 32 9 8 30
M  V30 33 9 3 30
M  V30 34 9 12 30
M  V30 35 1 19 20
M  V30 36 1 21 22
M  V30 37 1 24 25 CFG=3
M  V30 38 1 27 26
M  V30 39 1 28 27
M  V30 40 1 28 29 CFG=1
M  V30 41 1 21 13
M  V30 42 1 31 32 CFG=3
M  V30 43 1 33 34 CFG=1
M  V30 44 1 35 36
M  V30 45 1 35 37
M  V30 46 1 38 39
M  V30 47 1 38 40
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB";

  bool sanitize = false;
  bool removeHs = false;
  std::unique_ptr<RWMol> m(MolBlockToMol(molb, sanitize, removeHs));

  MolOps::sanitizeMol(*m);
  MolOps::assignStereochemistry(*m);

  // This should not throw an invariant violation
  auto smiles = MolToSmiles(*m);
  CHECK(
      smiles ==
      R"SMI(CC[C@@]1(C)/C2=C(C)/C3=[N]4->[CoH2+]56[N]2[C@H]([C@@H]1C)[C@]1(C)[N]->5=C(/C(C)=C2[N]->6=C(/C=C4/C(C)(C)[C@@H]3C)[C@@H](C)C\2(C)C)[C@@H](C)C1(C)C)SMI");
}

TEST_CASE("chiral presence and ranking") {
  SECTION("basics") {
    auto mol = "OC(F)C([C@H](F)O)[C@@H](F)O"_smiles;
    REQUIRE(mol);
    std::vector<unsigned int> ranks;
    bool breakTies = false;
    // default: all three centers are different
    Canon::rankMolAtoms(*mol, ranks, breakTies);
    CHECK(ranks[1] != ranks[4]);
    CHECK(ranks[1] != ranks[7]);
    CHECK(ranks[4] != ranks[7]);

    // if we don't include chirality, the ranks should be the same
    bool includeChirality = false;
    Canon::rankMolAtoms(*mol, ranks, breakTies, includeChirality);
    CHECK(ranks[1] == ranks[4]);
    CHECK(ranks[1] == ranks[7]);

    // if we include chiral presence, 4 and 7 are the same, but 1 is
    // different
    includeChirality = false;
    bool includeChiralPresence = true;
    bool includeIsotopes = true;
    bool includeAtomMaps = true;
    const bool includeStereoGroups = includeChirality;
    const bool useNonStereoRanks = false;

    Canon::rankMolAtoms(*mol, ranks, breakTies, includeChirality,
                        includeIsotopes, includeAtomMaps, includeChiralPresence,
                        includeStereoGroups, useNonStereoRanks);
    CHECK(ranks[1] != ranks[4]);
    CHECK(ranks[1] != ranks[7]);
    CHECK(ranks[4] == ranks[7]);
  }
}
TEST_CASE("meso impact on atom ranking") {
  SECTION("basics") {
    {
      auto m = "C1[C@H](O)C[C@@H]1[C@H](F)[C@@H]1C[C@H](O)C1"_smiles;
      REQUIRE(m);
      std::vector<unsigned int> ranks;
      bool breakTies = false;
      Canon::rankMolAtoms(*m, ranks, breakTies);
      CHECK(ranks[4] == ranks[7]);
    }
    {
      auto m = "C1[C@H](O)C[C@@H]1[C@H](F)[C@@H]1C[C@@H](O)C1"_smiles;
      REQUIRE(m);
      std::vector<unsigned int> ranks;
      bool breakTies = false;
      Canon::rankMolAtoms(*m, ranks, breakTies);
      CHECK(ranks[4] != ranks[7]);
    }
    {
      auto m =
          "C[C@@H](Cl)C1[C@H](C)Cl.C[C@@H](Cl)C2[C@H](C)Cl.[C@H]12F"_smiles;
      REQUIRE(m);
      std::vector<unsigned int> ranks;
      bool breakTies = false;
      Canon::rankMolAtoms(*m, ranks, breakTies);
      CHECK(ranks[3] == ranks[10]);
    }
    {
      auto m =
          "C[C@@H](Cl)C1[C@@H](C)Cl.C[C@@H](Cl)C2[C@H](C)Cl.[C@H]12F"_smiles;
      REQUIRE(m);
      std::vector<unsigned int> ranks;
      bool breakTies = false;
      Canon::rankMolAtoms(*m, ranks, breakTies);
      CHECK(ranks[3] != ranks[10]);
    }
  }
}

TEST_CASE("allow disabling ring stereo in ranking") {
  UseLegacyStereoPerceptionFixture reset_stereo_perception(false);

  std::string smi = "C[C@@H]1CC[C@@H](C)CC1";
  auto m = v2::SmilesParse::MolFromSmiles(smi);
  REQUIRE(m);

  bool breakTies = false;
  bool includeChirality = true;
  bool includeIsotopes = true;
  bool includeAtomMaps = true;
  bool includeChiralPresence = false;
  bool includeStereoGroups = true;
  bool useNonStereoRanks = false;
  bool includeRingStereo = true;
  std::vector<unsigned int> res1;
  Canon::rankMolAtoms(*m, res1, breakTies, includeChirality, includeIsotopes,
                      includeAtomMaps, includeChiralPresence,
                      includeStereoGroups, useNonStereoRanks,
                      includeRingStereo);
  CHECK(res1.size() == m->getNumAtoms());
  CHECK(res1[2] != res1[7]);
  CHECK(res1[2] == res1[3]);
  CHECK(res1[3] != res1[6]);
  CHECK(res1[6] == res1[7]);

  includeRingStereo = false;
  Canon::rankMolAtoms(*m, res1, breakTies, includeChirality, includeIsotopes,
                      includeAtomMaps, includeChiralPresence,
                      includeStereoGroups, useNonStereoRanks,
                      includeRingStereo);
  CHECK(res1.size() == m->getNumAtoms());
  CHECK(res1[2] == res1[7]);
  CHECK(res1[2] == res1[3]);
  CHECK(res1[3] == res1[6]);
  CHECK(res1[6] == res1[7]);
}

TEST_CASE("Canonicalization issues watch (see GitHub Issue #8775)") {
  // This is a check about the state of things with canonicalization.
  // The "samples" below initially come from the list compiled in GitHub
  // Issue #8775 by bp-kelly. Please update this list when things get
  // fixed (if they break, please reconsider the changes you are making!).
  // The comment on each sample is the GitHub issue where the issue was
  // first reported.

  const static std::initializer_list<std::tuple<std::string, bool, bool>> samples = {
      {R"smi(C/C=C\C=C(/C=C\C)C(/C=C\C)=C/C)smi", false, false},        // #8759
      {R"smi(C1=C\CCCCCC/C=C/C=C/1)smi", true, true},                   // #8759
      {R"smi(O=C=NC1=CC2C3=C(C=C1)C2=C(N=C=O)C=C3)smi", false, false},  // #8721
      {R"smi(O=C(c1ccccc1C(=O)N1C(=O)c2ccccc2C1=O)N1C(=O)c2ccccc2C1=O)smi",
       false, false},                                                   // #8721
      {R"smi(O=C=NC1=CC2C3=C(C=C1)C2=C(N=C=O)C=C3)smi", false, false},  // #8721
      {R"smi(O=[N+]([O-])c1cc/c2c(c1)=C(c1ccccc1)/N=c1\\ccc([N+](=O)[O-])cc1=C(c1ccccc1)/N=2)smi",
       false, true},  // #8721
      {R"smi(C=Cc1c(C)/c2[n-]c1=C=c1[n-]/c(c(CC)c1C)=C\\c1[n-]c3c(c1C)C(=O)[C@H](C(=O)OC)/C3=C1/[NH+]=C(/C=2)[C@@H](C)[C@@H]1CCC(=O)OC/C=C(\\C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C.[Mg+2])smi",
       false, true},  // #8721
      {R"smi(CC1=C(/C=C2\C(C)=C3/C(C)=C(/C=C4\C(C)=C5/C(CCC(=O)O)=C(C(C)=C5N4)C=C6\C(CCC(=O)O)=C(C)C(=C6N2)C=C1)N3)C=C(\C=C)C)smi",
       false, false},  // #8089
      {R"smi(COC1=N\C2=CC(=O)c3c(c(O)c(C)c4c3C(=O)C(C)(O/C=C/C(OC)C(C)C(OC(C)=O)C(C)C(O)C(C)C(O)C(C)/C=C\C=C/1C)O4)C2=O)smi",
       false, false},  // #8089
      {R"smi(CC1C2c3cc4/c5c6c7c8c9c%10c6c4c4c3C3=C6c%11c%12c%13c%14c%15c(c9c9c%16c8c(c8/c(c%17c%18c%19c(c(c%11c%11c%19c%19c%17c8c%16c(c%149)c%19c%13%11)C62)C1C%181C[N+](C)(C)C1)=C\C\C=5)C7)C%10C4C31C[N+](C)(C)CC%12%151)smi",
       false, true},  // #8089
      {R"smi(CC1=C\[C@H](C)C[C@@]2(C)CC[C@@H](O2)[C@@]23CC[C@@](C)(C[C@@H](O2)[C@H]2O[C@](C)(CC2=O)[C@@H](O)[C@@H]2CC[C@@]4(CCC[C@H](O4)[C@@H](C)C(=O)O[C@@H]4C[C@@H]([C@@]5(O)OCC[C@@H](C)[C@H]5O)O[C@@H]4/C=C/1)O2)O3)smi",
       false, false},  // #8089
      {R"smi(CC1=C/[C@H]2O[C@@H](C/C=C/C=C/C(=O)O[C@@H]3C[C@@H](/C=C/C/C=C/1)O[C@@H](C/C=C\CCO)[C@]3(C)CO)C[C@H](O)[C@H]2C)smi",
       false, false},  // #8089
      {R"smi(CC(=O)OCC1=C\CC/C(C)=C/CC[C@@]2(C)CC[C@@](C(C)C)(/C=C/1)O2)smi",
       false, false},  // #8089
      {R"smi(CC1=C\C/C=C(\C)CC[C@H]2C(C)(C)[C@@H](\C=C/1)CC[C@]2(C)O)smi",
       false, true},  // #8089
      {R"smi(CC(=O)OCC1=C/[C@@H]2OC(=O)[C@H](C)[C@@]2(O)[C@@H](OC(C)=O)[C@H]2[C@]3(CC[C@H](OC(C)=O)[C@]2(C)[C@@H](OC(=O)COC(=O)CC(C)C)\C=C/1)CO3)smi",
       false, false},  // #8089
      {R"smi(CC(=O)OCC1=C/[C@@H]2OC(=O)[C@H](C)[C@@]2(O)[C@@H](OC(C)=O)[C@H]2[C@]3(CC[C@H](OC(C)=O)[C@]2(C)[C@@H](OC(C)=O)\C=C/1)CO3)smi",
       false, false},  // #8089
      {R"smi(C=Cc1c(C)/c2[nH]/c1=C\C1=N/C(=C\C3=C(C)C4C(=O)N(Cc5cccc(C#Cc6cccc(Nc7ncnc8cc(OCCOC)c(OCCOC)cc78)c6)c5)C(=O)/C(=C5/N=C(/C=2)[C@@H](C)[C@@H]5CCC(=O)OC)C4N3)C(CC)=C1C)smi",
       false, true},  // #8089
      {R"smi(CC1=C\C[C@H](O)/C=C/C(C)=C/[C@@H](NC(=O)[C@H](C)O)[C@]2(C)C(=O)O[C@H](C[C@H](O)/C=C/1)[C@@H](C)C2=O)smi",
       false, false},  // #8089
      {R"smi(CC1=C/[C@H]2O[C@@H](C/C=C/C=C/C(=O)O[C@@H]3C[C@@H](/C=C/C/C=C/1)O[C@@H](C/C=C\C[C@@H](O)C(=O)O)[C@]3(C)CO)C[C@H](O)[C@H]2C)smi",
       false, false},  // #8089
      {R"smi(c1ccc2/c3[nH]/c(c2c1)=N\c1ccc(cc1)-c1nc2cc(ccc2o1)/N=c1/[nH]/c(c2ccccc12)=N/c1ccc2nc(oc2c1)-c1ccc(cc1)/N=3)smi",
       false, false},  // #8089
      {R"smi(c1ccc2/c3[nH]/c(c2c1)=N\c1ccc(cc1)-c1nc2ccc(cc2o1)/N=c1/[nH]/c(c2ccccc12)=N/c1ccc2nc(oc2c1)-c1ccc(cc1)/N=3)smi",
       false, false},  // #8089
      {R"smi(CC1=C/[C@@H](C)C[C@]2(C)CC[C@H](O2)[C@]23CC[C@](C)(C[C@H](O2)[C@@H]2O[C@@](C)(CC2=O)[C@@H](O)[C@H]2CC[C@@]4(CCC[C@@H](O4)[C@H](C)C(=O)O[C@H]4C[C@H]([C@]5(O)OCC[C@H](C)[C@@H]5O)O[C@H]4\C=C/1)O2)O3)smi",
       false, false},  // #8089
      {R"smi(CC1=C/[C@@H](C)C[C@]2(C)CC[C@H](O2)[C@]23CC[C@](C(=O)O)(C[C@H](O2)[C@@H]2O[C@@](C)(CC2=O)[C@@H](O)[C@H]2CC[C@@]4(CCC[C@@H](O4)[C@H](C)C(=O)O[C@H]4C[C@H]([C@]5(O)OCC[C@H](C)[C@@H]5O)O[C@H]4\C=C/1)O2)O3)smi",
       false, false},  // #8089
      {R"smi(CC1=C/[C@@H](C)C[C@]2(C)CC[C@H](O2)[C@]23CC[C@](CO)(C[C@H](O2)[C@@H]2O[C@@](C)(CC2=O)[C@@H](O)[C@H]2CC[C@@]4(CCC[C@@H](O4)[C@H](C)C(=O)O[C@H]4C[C@H]([C@]5(O)OCC[C@H](C)[C@@H]5O)O[C@H]4\C=C/1)O2)O3)smi",
       false, false},  // #8089
      {R"smi(COC(=O)C1=C/c2cc3c(cc2-c2c(cc(OC)c(OC)c2OC)\C=C/1C(=O)OC)OCO3)smi",
       false, false},                                           // #8089
      {R"smi(N#C[P@@H]/C=C(/P=O)[P@@H]C#N)smi", false, false},  // #8089
      {R"smi(C1=C\C/C=C(\C)CC[C@H]2C(C)(C)[C@@H](\C=C/1)CC[C@]2(C)O)smi", true,
       true},                                                     // #8089
      {R"smi([H]/N=C(C=C)\C(/N=C\[O-])=N\[H])smi", false, true},  // #7759
      {R"smi([H]N=C(/N=C\[O-])/C(C=C)=N\[H])smi", true, true},    // #7759
      {R"smi([H]/N=C(/C=C)C(=N)/N=C\[O-])smi", true, true},       // #7759
      {R"smi([H]/N=C(/C=C)C(=N)/N=C\[O-])smi", true, true},       // #7759
      {R"smi([C@H]12[C@H]3[C@@H]4[C@H]5[C@@H]([C@H]1N24)N53)smi", false,
       true},  // #7759
      {R"smi([C@H]12[C@H]3[C@H]4[C@@H]5[C@@H]([C@H]1N25)N34)smi", false,
       true},  // #7759
      {R"smi([C@H]12[C@H]3[C@@H]4[C@H]5[C@@H]([C@H]1N24)N53)smi", false,
       true},  // #7759
      {R"smi([C@H]12[C@H]3[C@H]4[C@@H]5[C@@H]([C@H]1N25)N34)smi", false,
       true},                                                      // #7759
      {R"smi([CH]1C[C@H]2C[C@H](C2)[C@]12CC2)smi", false, true},   // #7759
      {R"smi([CH]1O[C@H]2C[C@H](C2)[C@]12CC2)smi", false, true},   // #7759
      {R"smi([CH]1O[C@H]2CN(C2)[C@]12CC2)smi", false, true},       // #7759
      {R"smi([CH]1[C@H]2C[C@H](C2)[C@]12CCC2)smi", false, true},   // #7759
      {R"smi([CH]1[C@H]2C[C@H](C2)[C@@]12COC2)smi", false, true},  // #7759
      {R"smi(O=C1[C]2C[C@H](C2)[C@@]12CC2)smi", false, false},     // #7759
      {R"smi(O=C1[C@H]2C[C](C2)[C@@]12CC2)smi", false, true},      // #7759
      {R"smi(C1[C@H]2[C@@H]3C[C@H]4[C@@H]3[C@@H]1[C@@H]24)smi", false,
       true},  // #7759
      {R"smi([C@H]12[C@H]3[C@@H]4[C@H]5[C@@H]([C@H]1N24)N53)smi", false,
       true},  // #7759
      {R"smi(C[C@]12[C@@H]3[C@@H]1C[C@H]1[C@H]2[C@]31C)smi", false,
       true},                                                      // #7759
      {R"smi([O][C@]12C[C@H](C1)[C@@]1(CC1)C2)smi", false, true},  // #7759
      {R"smi(O[C@]12[C@@H]3[C@@H]4[CH][C@H]1[C@H]2[C@]34O)smi", false,
       true},                                                        // #7759
      {R"smi([CH2][C@]12C[C@H](C1)[C@]1(CC1)C2)smi", false, true},   // #7759
      {R"smi([CH2][C@]12C[C@H](C1)[C@@]1(CC1)O2)smi", false, true},  // #7759
      {R"smi([CH2][C@]1(C)[C@H]2[C@@H]3[C@@H](C)[C@H]2[C@@H]31)smi", false,
       true},  // #7759
      {R"smi([CH2][C@]1(C)[C@H]2[C@@H]3[C@H](O)[C@H]2[C@@H]31)smi", false,
       true},  // #7759
      {R"smi(C[C@H]1[C@@H]2[C@@H]3[C@H]1[C@H]2[C@@]31[CH]C1)smi", false,
       true},  // #7759
      {R"smi(O[C@H]1[C@@H]2[C@@H]3[C@H]1[C@H]2[C@@]31[CH]C1)smi", false,
       true},  // #7759
      {R"smi([O][C@H]1[C@@H]2[C@@H]3[C@H]1[C@H]2[C@]31CN1)smi", false,
       true},  // #7759
      {R"smi(O[C@H]1[C@@H]2[C@@H]3[C@H]1[C@H]2[C@]31[CH]N1)smi", false,
       true},  // #7759
      {R"smi(O[C@H]1[C@@H]2[C@@H]3[C@H]1[C@H]2[C@]31C[N]1)smi", false,
       true},  // #7759
      {R"smi([CH2][C@]1(O)[C@H]2[C@@H]3[C@@H](O)[C@H]2[C@@H]31)smi", false,
       true},  // #7759
      {R"smi(C[C@]1([O])[C@H]2[C@@H]3[C@@H](O)[C@H]2[C@@H]31)smi", false,
       true},  // #7759
      {R"smi([CH2][C@H]1[C@@H]2[C@@H]3[C@H]1[C@H]2[C@]31CO1)smi", false,
       true},  // #7759
      {R"smi(C[C@]1(O)[C@@H]2[C@H]3[C@@H]([O])[C@@H]2[C@H]31)smi", false,
       true},  // #7759
      {R"smi(C[C@H]1[C@H]2[C@H]3[C@@H]1[C@@H]2[C@@]31[CH]O1)smi", false,
       true},  // #7759
      {R"smi([O][C@H]1[C@H]2[C@H]3[C@@H]1[C@@H]2[C@@]31CO1)smi", false,
       true},  // #7759
      {R"smi(O[C@H]1[C@@H]2[C@@H]3[C@H]1[C@H]2[C@]31[CH]O1)smi", false,
       true},  // #7759
      {R"smi([CH2][C@]1(O)[C@@H]2[C@H]3[C@H](C)[C@@H]2[C@H]31)smi", false,
       true},  // #7759
      {R"smi(C[C@H]1[C@H]2[C@H]3[C@@H]1[C@@H]2[C@]3(C)[O])smi", false,
       true},  // #7759
      {R"smi([CH2][C@H]1[C@@H]2[C@@H]3[C@H]1[C@H]2[C@@]3(C)O)smi", false,
       true},                                                           // #7759
      {R"smi([O][C@H]1[C@H]2C[C@H](C2)[C@@]12CC2)smi", false, true},    // #7759
      {R"smi(O[C]1[C@H]2C[C@H](C2)[C@@]12CC2)smi", false, true},        // #7759
      {R"smi([CH2][C@H]1[C@H]2C[C@H](C2)[C@@]12CC2)smi", false, true},  // #7759
      {R"smi(O[C@H]1[C]2C[C@H](C2)[C@]12CC2)smi", false, true},         // #7759
      {R"smi(C[C]1[C@H]2C[C@H](C2)[C@@]12CC2)smi", false, true},        // #7759
      {R"smi(C[C@H]1[C]2C[C@H](C2)[C@]12CC2)smi", false, true},         // #7759
      {R"smi(O[C@H]1[C@H]2C[C](C2)[C@@]12CC2)smi", false, true},        // #7759
      {R"smi(C[C@H]1[C@H]2C[C](C2)[C@@]12CC2)smi", false, true},        // #7759
      {R"smi(C[C@]12[C@@H]3[C@@H]1[CH-][C@H]1[C@H]2[C@]31C)smi", false,
       true},                                                           // #7759
      {R"smi(C[C@]12C[C@H](C1)[N@+]1(CC1)C2)smi", false, true},         // #7759
      {R"smi(C[C@H]1[C@H]2C[C-](C2)[C@@]12CC2)smi", false, true},       // #7759
      {R"smi([CH]1[C@H]2C[C@H](C2)[C@]12CC2)smi", false, true},         // #7759
      {R"smi([CH]1N2C[C@H](C2)O[C@]12CC2)smi", false, true},            // #7759
      {R"smi([C]1=C[C@H]2C[C@H](C2)[C@@]12CC2)smi", false, true},       // #7759
      {R"smi([CH-]1O[C@H]2CN(C2)[C@]12CC2)smi", false, true},           // #7759
      {R"smi([CH-]1C[C@H]2C[C@H](C2)[C@]12CC2)smi", false, true},       // #7759
      {R"smi([CH-]1O[C@H]2C[C@H](C2)[C@]12CC2)smi", false, true},       // #7759
      {R"smi([CH-]1[C@H]2C[C@H](C2)[C@]12CCC2)smi", false, true},       // #7759
      {R"smi([CH-]1[C@H]2C[C@H](C2)[C@@]12COC2)smi", false, true},      // #7759
      {R"smi([O-]C1=C2C[C@H](C2)[C@@]12CC2)smi", false, false},         // #7759
      {R"smi(O=C1[C@H]2C[C-](C2)[C@@]12CC2)smi", false, true},          // #7759
      {R"smi([O-][C@]12C[C@H](C1)[C@@]1(CC1)C2)smi", false, true},      // #7759
      {R"smi([CH2-][C@]12C[C@H](C1)[C@]1(CC1)C2)smi", false, true},     // #7759
      {R"smi(O[C@H]1[C-]2C[C@H](C2)[C@@]12CC2)smi", false, true},       // #7759
      {R"smi([CH2-][C@H]1[C@H]2C[C@H](C2)[C@]12CC2)smi", false, true},  // #7759
      {R"smi([O-][C@H]1[C@H]2C[C@H](C2)[C@@]12CC2)smi", false, true},   // #7759
      {R"smi(O[C-]1[C@H]2C[C@H](C2)[C@]12CC2)smi", false, true},        // #7759
      {R"smi(O[C@H]1[C@H]2C[C-](C2)[C@@]12CC2)smi", false, true},       // #7759
      {R"smi(C[C-]1[C@H]2C[C@H](C2)[C@]12CC2)smi", false, true},        // #7759
      {R"smi(C[C@H]1[C-]2C[C@H](C2)[C@@]12CC2)smi", false, true},       // #7759
      {R"smi([CH-]1[C@H]2C[C@H](C2)[C@]12CC2)smi", false, true},        // #7759
      {R"smi([C-]1=C[C@H]2C[C@H](C2)[C@@]12CC2)smi", false, true},      // #7759
      {R"smi([CH-]1N2C[C@H](C2)O[C@]12CC2)smi", false, true},           // #7759
      {R"smi(C1=[O+][C@H]2C[C@H](C2)[C@]12CC2)smi", false, true},       // #7759
      {R"smi(C1=[O+][C@H]2CN(C2)[C@]12CC2)smi", false, true},           // #7759
      {R"smi([CH+]1[C@H]2C[C@H](C2)[C@]12CCC2)smi", false, true},       // #7759
      {R"smi([CH+]1[C@H]2C[C@H](C2)[C@@]12COC2)smi", false, true},      // #7759
      {R"smi(C[C@]12[C@@H]3[C@@H]1[CH][C@H]1[C@H]2[C@]31C)smi", false,
       true},                                                           // #7759
      {R"smi(C[C+]1[C@H]2C[C@H](C2)[C@]12CC2)smi", false, true},        // #7759
      {R"smi([H]/[O+]=C1/[C@H]2C[C@H](C2)[C@]12CC2)smi", false, true},  // #7759
      {R"smi([CH+]1[C@H]2C[C@H](C2)[C@]12CC2)smi", false, true},        // #7759
      {R"smi(C1[C@H]2O[C@@H]3C[C@H]4[C@@H]3[C@@H]1[C@@H]24)smi", false,
       true},                                                          // #7759
      {R"smi(O=C1N2[C@H]3[C@@H]2[C@@H]2[C@H]3N12)smi", false, false},  // #7759
      {R"smi(O1[C@H]2[C@@H]3O[C@H]4[C@@H]3[C@@H]1[C@@H]24)smi", false,
       true},                                                    // #7759
      {R"smi(O=C1[C@H]2C[C@H](C2)[C@@]12CC2)smi", false, true},  // #7759
      {R"smi(O=C1[C@H]2CN(C2)[C@@]12CC2)smi", false, true},      // #7759
      {R"smi(O[C@]12C[C@H](C1)[C@]1(CC1)C2)smi", false, true},   // #7759
      {R"smi(O[C@]12[C@@H]3[C@@H]4O[C@H]1[C@H]2[C@]34O)smi", false,
       true},                                                    // #7759
      {R"smi(C[C@]12C[C@H](C1)[C@@]1(CC1)O2)smi", false, true},  // #7759
      {R"smi(C[C@]1(O)[C@H]2[C@@H]3[C@@H](O)[C@H]2[C@@H]31)smi", false,
       true},  // #7759
      {R"smi(C[C@H]1[C@@H]2[C@@H]3[C@H]1[C@H]2[C@@]3(C)O)smi", false,
       true},  // #7759
      {R"smi(O[C@H]1[C@@H]2[C@@H]3[C@H]1[C@H]2[C@]31CN1)smi", false,
       true},                                                      // #7759
      {R"smi(C[C@H]1[C@H]2C[C@H](C2)[C@]12CC2)smi", false, true},  // #7759
      {R"smi(C[C@]12[C@@H]3[C@@H]4O[C@H]1[C@H]2[C@]34C)smi", false,
       true},                                                      // #7759
      {R"smi(O[C@H]1[C@H]2C[C@H](C2)[C@]12CC2)smi", false, true},  // #7759
  };

  auto roundtrip_smiles = [](const std::string &smi,
                             bool use_legacy_stereo) -> std::string {
    UseLegacyStereoPerceptionFixture use_legacy(use_legacy_stereo);
    auto m = v2::SmilesParse::MolFromSmiles(smi);
    REQUIRE(m);
    return MolToSmiles(*m);
  };

  const auto &[smiles, legacy_state, modern_state] =
      GENERATE_REF(values(samples));
  auto check_legacy_stereo = GENERATE(false, true);
  CAPTURE(smiles, check_legacy_stereo);

  // pre-canonicalize SMILES: the inputs get outdated when
  // we make changes to the canonicalization algorithm
  const auto first_roundtrip = roundtrip_smiles(smiles, check_legacy_stereo);

  auto should_match = check_legacy_stereo ? legacy_state : modern_state;
  if (should_match) {
    CHECK(first_roundtrip ==
          roundtrip_smiles(first_roundtrip, check_legacy_stereo));
  } else {
    CHECK(first_roundtrip !=
          roundtrip_smiles(first_roundtrip, check_legacy_stereo));
  }
}