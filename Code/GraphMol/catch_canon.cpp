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
      R"SMI(CC[C@@]1(C)/C2=C(C)/C3=N4->[CoH2+]56N2[C@H]([C@@H]1C)[C@]1(C)N->5=C(/C(C)=C2N->6=C(/C=C4/C(C)(C)[C@@H]3C)[C@@H](C)C\2(C)C)[C@@H](C)C1(C)C)SMI");
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

    // if we include chiral presence, 4 and 7 are the same, but 1 is different
    includeChirality = false;
    bool includeChiralPresence = true;
    bool includeIsotopes = true;
    bool includeAtomMaps = true;
    Canon::rankMolAtoms(*mol, ranks, breakTies, includeChirality,
                        includeIsotopes, includeAtomMaps,
                        includeChiralPresence);
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