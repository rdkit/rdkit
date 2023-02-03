//
//  Copyright (C) 2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"

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
}
