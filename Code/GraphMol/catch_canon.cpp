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

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("chirality and canonicalization") {
  SECTION("basics") {
    auto mol = "F[C@](O)(Cl)C[C@](F)(O)Cl"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool force = true;
    MolOps::assignStereochemistry(*mol, cleanIt, force);
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
    std::string cip;
    CHECK(mol->getAtomWithIdx(1)->getPropIfPresent(common_properties::_CIPCode,
                                                   cip));
    CHECK(cip == "S");
    CHECK(mol->getAtomWithIdx(7)->getPropIfPresent(common_properties::_CIPCode,
                                                   cip));
    CHECK(cip == "R");
    CHECK(mol->getAtomWithIdx(4)->getPropIfPresent(common_properties::_CIPCode,
                                                   cip));
    CHECK(cip == "R");
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
    for (const auto& [smi1, smi2] : tests) {
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
    for (const auto& [smi1, smi2] : tests) {
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
    for (const auto& [smi1, smi2] : tests) {
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
    for (const auto& [smi1, smi2] : tests) {
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

  SECTION("psuedoTest1") {
    std::vector<std::pair<std::string, std::string>> tests = {
        // {"N[C@H]1CC[C@@H](O)CC1", "N[C@@H]1CC[C@H](O)CC1"},
        // {"N[C@H]1CC[C@@H](O)CC1 |o2:1,4|", "N[C@@H]1CC[C@H](O)CC1 |o2:1,4|"},
        // {"N[C@H]1CC[C@@H](O)CC1 |&2:1,4|", "N[C@@H]1CC[C@H](O)CC1 |&2:1,4|"},
        // {"N[C@H]1CC[C@@H](O)CC1 |a:1,4|", "N[C@@H]1CC[C@H](O)CC1 |a:1,4|"},
        {"CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |(25.6031,-3.1625,;26.9021,-3.9125,;28.201,-3.1625,;26.9021,-5.4125,;28.201,-6.1625,;28.201,-7.6625,;26.9021,-8.4125,;25.6031,-7.6625,;25.6031,-6.1625,;24.3042,-5.4125,;24.3042,-3.9125,;23.0052,-6.1625,;23.0052,-7.6625,;21.7062,-8.4125,;20.4073,-7.6625,;19.1083,-8.4125,;20.4073,-6.1625,;21.7062,-5.4125,),wD:11.11,14.15,wU:3.2,a:3,11,14|",
         "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |(12.3219,-4.25625,;13.6208,-5.00625,;14.9198,-4.25625,;13.6208,-6.50625,;14.9198,-7.25625,;14.9198,-8.75625,;13.6208,-9.50625,;12.3219,-8.75625,;12.3219,-7.25625,;11.0229,-6.50625,;11.0229,-5.00625,;9.72396,-7.25625,;8.425,-6.50625,;7.12604,-7.25625,;7.12604,-8.75625,;5.82708,-9.50625,;8.425,-9.50625,;9.72396,-8.75625,),wU:3.2,11.11,14.15,a:3,11,14|"},
        {"CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |(25.6031,-3.1625,;26.9021,-3.9125,;28.201,-3.1625,;26.9021,-5.4125,;28.201,-6.1625,;28.201,-7.6625,;26.9021,-8.4125,;25.6031,-7.6625,;25.6031,-6.1625,;24.3042,-5.4125,;24.3042,-3.9125,;23.0052,-6.1625,;23.0052,-7.6625,;21.7062,-8.4125,;20.4073,-7.6625,;19.1083,-8.4125,;20.4073,-6.1625,;21.7062,-5.4125,),wD:11.11,14.15,wU:3.2|",
         "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |(12.3219,-4.25625,;13.6208,-5.00625,;14.9198,-4.25625,;13.6208,-6.50625,;14.9198,-7.25625,;14.9198,-8.75625,;13.6208,-9.50625,;12.3219,-8.75625,;12.3219,-7.25625,;11.0229,-6.50625,;11.0229,-5.00625,;9.72396,-7.25625,;8.425,-6.50625,;7.12604,-7.25625,;7.12604,-8.75625,;5.82708,-9.50625,;8.425,-9.50625,;9.72396,-8.75625,),wU:3.2,11.11,14.15|"},
        {"CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@@H](C)CC1 |(11.2281,-31.6,;12.5271,-32.35,;13.826,-31.6,;12.5271,-33.85,;13.826,-34.6,;13.826,-36.1,;12.5271,-36.85,;11.2281,-36.1,;11.2281,-34.6,;9.92917,-33.85,;9.92917,-32.35,;8.63021,-34.6,;7.33125,-33.85,;6.03229,-34.6,;6.03229,-36.1,;4.73333,-36.85,;7.33125,-36.85,;8.63021,-36.1,),wU:3.2,11.11,14.15,&1:14|", "CC(C)[C@H]1CCCCN1C(=O)[C@H]1CC[C@H](C)CC1 |(24.1187,-31.2875,;25.4177,-32.0375,;26.7167,-31.2875,;25.4177,-33.5375,;26.7167,-34.2875,;26.7167,-35.7875,;25.4177,-36.5375,;24.1187,-35.7875,;24.1187,-34.2875,;22.8198,-33.5375,;22.8198,-32.0375,;21.5208,-34.2875,;21.5208,-35.7875,;20.2219,-36.5375,;18.9229,-35.7875,;17.624,-36.5375,;18.9229,-34.2875,;20.2219,-33.5375,),wD:14.15,wU:3.2,11.11,&1:14|"}};
    for (const auto& [smi1, smi2] : tests) {
      INFO(smi1 + " : " + smi2);
      std::unique_ptr<RWMol> mol1{SmilesToMol(smi1)};
      REQUIRE(mol1);
      std::unique_ptr<RWMol> mol2{SmilesToMol(smi2)};
      REQUIRE(mol2);

      Canon::canonicalizeEnhancedStereo(*mol1);
      Canon::canonicalizeEnhancedStereo(*mol2);
      auto outSmi1 = MolToCXSmiles(*mol1);
      auto outSmi2 = MolToCXSmiles(*mol2);

      CHECK(outSmi1 == outSmi2);
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
    for (auto& smi : smis) {
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
    for (auto& smi : smis) {
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
    CHECK(cx1 != cx2);
    CHECK(cx1.find("&7:") != std::string::npos);
    CHECK(cx1.find("&8:") != std::string::npos);
    CHECK(cx2.find("&2:") != std::string::npos);
    CHECK(cx2.find("&3:") != std::string::npos);
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

    auto cx1 = MolToCXSmiles(*m1);
    auto cx2 = MolToCXSmiles(*m2);
    CHECK(MolToCXSmiles(*m1) == MolToCXSmiles(*m2));

    // "read" ids are also reset!
    forwardStereoGroupIds(*m1);
    forwardStereoGroupIds(*m2);

    cx1 = MolToCXSmiles(*m1);
    cx2 = MolToCXSmiles(*m2);
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

  auto countStereoBonds = [](const auto& mol) {
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