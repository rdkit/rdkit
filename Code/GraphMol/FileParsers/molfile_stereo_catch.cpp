//
//  Copyright (C) 2022 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>

using namespace RDKit;

TEST_CASE("Github #5863: failure in WedgeMolBonds") {
  SECTION("as reported") {
    auto mol =
        "C[C@]12C=CC(=O)C=C1CC[C@@H]1[C@@H]2[C@@H](O)C[C@@]2(C)[C@H]1CC[C@]2(O)C(=O)CO |(1.94354,1.43772,;2.70098,0.14301,;3.44351,1.44633,;4.94349,1.45494,;5.70093,0.160228,;7.20091,0.168838,;4.9584,-1.14309,;3.45842,-1.1517,;2.71589,-2.45502,;1.21592,-2.46363,;0.458474,-1.16892,;1.20101,0.1344,;0.443562,1.42911,;1.18609,2.73243,;-1.05641,1.4205,;-1.79895,0.117181,;-2.53364,1.42493,;-1.0415,-1.17753,;-2.03878,-2.29799,;-3.41258,-1.69576,;-3.26435,-0.203103,;-3.6102,1.25648,;-4.76433,-0.211712,;-5.50686,-1.51503,;-5.52177,1.083,;-7.02175,1.07439,)|"_smiles;
    REQUIRE(mol);
    unsigned int radius = 2;
    unsigned int atomId = 2;
    auto env = findAtomEnvironmentOfRadiusN(*mol, radius + 1, atomId);
    std::unique_ptr<ROMol> frag(Subgraphs::pathToSubmol(*mol, env));
    REQUIRE(frag);
    Chirality::wedgeMolBonds(*frag, &frag->getConformer());
    INFO(MolToV3KMolBlock(*frag));
    CHECK(frag->getBondBetweenAtoms(9, 10)->getBondDir() !=
          Bond::BondDir::NONE);
  }
}

TEST_CASE("translating the chiral flag to stereo groups") {
  SECTION("basics") {
    auto withFlag = R"CTAB(
  Mrv2211 03302308372D          

  5  4  0  0  1  0            999 V2000
   -6.5625    3.9286    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8480    4.3411    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   -5.1336    3.9286    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1775    5.0826    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4384    5.0893    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  5  1  1  0  0  0
  2  4  1  6  0  0  0
M  END
)CTAB"_ctab;
    REQUIRE(withFlag);
    int flag = 0;
    CHECK(withFlag->getPropIfPresent(common_properties::_MolFileChiralFlag,
                                     flag));
    CHECK(flag == 1);
    RWMol zeroFlag(*withFlag);
    zeroFlag.setProp(common_properties::_MolFileChiralFlag, 0);
    RWMol noFlag(*withFlag);
    noFlag.clearProp(common_properties::_MolFileChiralFlag);
    CHECK(withFlag->getStereoGroups().empty());

    translateChiralFlagToStereoGroups(*withFlag);
    CHECK(!withFlag->hasProp(common_properties::_MolFileChiralFlag));
    auto sgs = withFlag->getStereoGroups();
    REQUIRE(sgs.size() == 1);
    CHECK(sgs[0].getGroupType() == StereoGroupType::STEREO_ABSOLUTE);
    CHECK(sgs[0].getAtoms().size() == 1);
    CHECK(sgs[0].getAtoms()[0]->getIdx() == 1);

    {
      RWMol cp(zeroFlag);
      translateChiralFlagToStereoGroups(cp);
      CHECK(!cp.hasProp(common_properties::_MolFileChiralFlag));
      auto sgs = cp.getStereoGroups();
      REQUIRE(sgs.size() == 1);
      CHECK(sgs[0].getGroupType() == StereoGroupType::STEREO_AND);
      CHECK(sgs[0].getAtoms().size() == 1);
      CHECK(sgs[0].getAtoms()[0]->getIdx() == 1);
    }
    {
      RWMol cp(zeroFlag);
      translateChiralFlagToStereoGroups(cp, StereoGroupType::STEREO_OR);
      CHECK(!cp.hasProp(common_properties::_MolFileChiralFlag));
      auto sgs = cp.getStereoGroups();
      REQUIRE(sgs.size() == 1);
      CHECK(sgs[0].getGroupType() == StereoGroupType::STEREO_OR);
      CHECK(sgs[0].getAtoms().size() == 1);
      CHECK(sgs[0].getAtoms()[0]->getIdx() == 1);
    }

    translateChiralFlagToStereoGroups(noFlag);
    CHECK(!noFlag.hasProp(common_properties::_MolFileChiralFlag));
    CHECK(noFlag.getStereoGroups().empty());
  }

  SECTION("explicit zero chiral flag") {
    auto zeroFlag = R"CTAB(
  Mrv2211 03302308372D          

  5  4  0  0  o  0            999 V2000
   -6.5625    3.9286    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8480    4.3411    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   -5.1336    3.9286    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1775    5.0826    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4384    5.0893    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  5  1  1  0  0  0
  2  4  1  6  0  0  0
M  END
)CTAB"_ctab;
    REQUIRE(zeroFlag);
    int flag = 0;
    CHECK(zeroFlag->getPropIfPresent(common_properties::_MolFileChiralFlag,
                                     flag));
    CHECK(flag == 0);
    {
      RWMol cp(*zeroFlag);
      translateChiralFlagToStereoGroups(cp);
      CHECK(!cp.hasProp(common_properties::_MolFileChiralFlag));
      auto sgs = cp.getStereoGroups();
      REQUIRE(sgs.size() == 1);
      CHECK(sgs[0].getGroupType() == StereoGroupType::STEREO_AND);
      CHECK(sgs[0].getAtoms().size() == 1);
      CHECK(sgs[0].getAtoms()[0]->getIdx() == 1);
    }
  }

  SECTION("multiple stereocenters") {
    auto noFlag = "C[C@@](N)(F)C[C@](C)(O)F"_smiles;
    REQUIRE(noFlag);
    RWMol withFlag(*noFlag);
    withFlag.setProp(common_properties::_MolFileChiralFlag, 1);

    CHECK(withFlag.getStereoGroups().empty());
    translateChiralFlagToStereoGroups(withFlag);
    CHECK(!withFlag.hasProp(common_properties::_MolFileChiralFlag));
    auto sgs = withFlag.getStereoGroups();
    REQUIRE(sgs.size() == 1);
    CHECK(sgs[0].getGroupType() == StereoGroupType::STEREO_ABSOLUTE);
    CHECK(sgs[0].getAtoms().size() == 2);
    CHECK(sgs[0].getAtoms()[0]->getIdx() == 1);
    CHECK(sgs[0].getAtoms()[1]->getIdx() == 5);

    {
      RWMol zeroFlag(*noFlag);
      zeroFlag.setProp(common_properties::_MolFileChiralFlag, 0);
      translateChiralFlagToStereoGroups(zeroFlag);
      CHECK(!zeroFlag.hasProp(common_properties::_MolFileChiralFlag));
      auto sgs = zeroFlag.getStereoGroups();
      REQUIRE(sgs.size() == 1);
      CHECK(sgs[0].getGroupType() == StereoGroupType::STEREO_AND);
      CHECK(sgs[0].getAtoms().size() == 2);
      CHECK(sgs[0].getAtoms()[0]->getIdx() == 1);
      CHECK(sgs[0].getAtoms()[1]->getIdx() == 5);
    }
  }

  SECTION("pre-existing ABS stereogroup") {
    auto withFlag = "C[C@@](N)(F)C[C@](C)(O)F |a:1|"_smiles;
    REQUIRE(withFlag);
    withFlag->setProp(common_properties::_MolFileChiralFlag, 1);

    CHECK(withFlag->getStereoGroups().size() == 1);
    translateChiralFlagToStereoGroups(*withFlag);
    CHECK(!withFlag->hasProp(common_properties::_MolFileChiralFlag));
    auto sgs = withFlag->getStereoGroups();
    REQUIRE(sgs.size() == 1);
    CHECK(sgs[0].getGroupType() == StereoGroupType::STEREO_ABSOLUTE);
    CHECK(sgs[0].getAtoms().size() == 2);
    CHECK(sgs[0].getAtoms()[0]->getIdx() == 1);
    CHECK(sgs[0].getAtoms()[1]->getIdx() == 5);
  }
}

TEST_CASE("github #6310: crossed ring bonds not written to mol blocks") {
  SECTION("as reported") {
    std::vector<std::string> smileses = {"C1C=CCCCCCCCCC1", "C1C=CCCCCC1"};

    for (const auto &smiles : smileses) {
      auto m = std::unique_ptr<RWMol>(SmilesToMol(smiles));
      REQUIRE(m);
      CHECK(m->getBondWithIdx(1)->getBondType() == Bond::BondType::DOUBLE);
      CHECK(m->getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREONONE);
      auto mb = MolToV3KMolBlock(*m);
      CHECK(mb.find("CFG=2") != std::string::npos);

      mb = MolToMolBlock(*m);
      CHECK(mb.find("2  3  2  3") != std::string::npos);

      // make sure we also write the crossed bond when it's explicitly crossed
      m->getBondWithIdx(1)->setStereo(Bond::BondStereo::STEREOANY);
      mb = MolToV3KMolBlock(*m);
      CHECK(mb.find("CFG=2") != std::string::npos);

      mb = MolToMolBlock(*m);
      CHECK(mb.find("2  3  2  3") != std::string::npos);
    }
  }
  SECTION("don't cross bonds in small rings") {
    auto m = "C1C=CCCCC1"_smiles;
    REQUIRE(m);
    CHECK(m->getBondWithIdx(1)->getBondType() == Bond::BondType::DOUBLE);
    CHECK(m->getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREONONE);
    auto mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("CFG=2") == std::string::npos);

    mb = MolToMolBlock(*m);
    CHECK(mb.find("2  3  2  3") == std::string::npos);

    m->getBondWithIdx(1)->setStereo(Bond::BondStereo::STEREOANY);
    mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("CFG=2") != std::string::npos);

    mb = MolToMolBlock(*m);
    CHECK(mb.find("2  3  2  3") != std::string::npos);
  }
}

void testStereoExample(const std::string &mb, unsigned int aidx,
                       Atom::ChiralType expected,
                       const std::string &expectedCIP) {
  INFO(mb);
  std::unique_ptr<RWMol> m(MolBlockToMol(mb));
  REQUIRE(m);
  CHECK(m->getAtomWithIdx(aidx)->getChiralTag() == expected);
  CIPLabeler::assignCIPLabels(*m);
  std::string CIP = "";
  CHECK(m->getAtomWithIdx(aidx)->getPropIfPresent(common_properties::_CIPCode,
                                                  CIP));
  CHECK(CIP == expectedCIP);
}

TEST_CASE("IUPAC recommendations") {
#if 1
  SECTION("simple examples") {
    std::vector<std::string> mbs = {
        R"CTAB(
  Mrv2211 06082308462D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.875 5.3333 0 0
M  V30 2 C -6.5413 6.1033 0 0
M  V30 3 O -5.2076 5.3333 0 0
M  V30 4 F -7.3292 7.3957 0 0
M  V30 5 Cl -5.4654 7.2052 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4 CFG=1
M  V30 4 1 2 5 CFG=3
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB",
        R"CTAB(
  Mrv2211 06082309052D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.875 5.3333 0 0
M  V30 2 C -6.5413 6.1033 0 0 CFG=2
M  V30 3 O -5.2076 5.3333 0 0
M  V30 4 F -6.5413 8 0 0
M  V30 5 Cl -5.4654 7.2052 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4 CFG=1
M  V30 4 1 2 5 CFG=3
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB",
        R"CTAB(
  Mrv2211 06082309052D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.875 5.3333 0 0
M  V30 2 C -6.5413 6.1033 0 0 CFG=2
M  V30 3 O -5.2076 5.3333 0 0
M  V30 4 F -4.9369 6.9639 0 0
M  V30 5 Cl -6.5413 7.9602 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4 CFG=1
M  V30 4 1 2 5 CFG=3
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB",
        R"CTAB(IUPAC does not like this one
  Mrv2211 06082309142D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.875 5.3333 0 0
M  V30 2 C -6.5413 6.1033 0 0
M  V30 3 O -5.2076 5.3333 0 0
M  V30 4 F -6.5413 4.2897 0 0
M  V30 5 Cl -6.5413 7.9602 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4 CFG=1
M  V30 4 1 2 5
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB",
        R"CTAB(IUPAC does not like this one2
  Mrv2211 06082309142D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.875 5.3333 0 0
M  V30 2 C -6.5413 6.1033 0 0
M  V30 3 O -5.2076 5.3333 0 0
M  V30 4 F -6.5413 4.2897 0 0
M  V30 5 Cl -6.5413 7.9602 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4 
M  V30 4 1 2 5 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB",
    };
    for (const auto &mb : mbs) {
      testStereoExample(mb, 1, Atom::ChiralType::CHI_TETRAHEDRAL_CCW, "S");
      if (mb.find("CFG=3") != std::string::npos &&
          mb.find("CFG=1") != std::string::npos) {
        std::string cp = mb;
        testStereoExample(cp.replace(mb.find("CFG=1"), 5, "     "), 1,
                          Atom::ChiralType::CHI_TETRAHEDRAL_CCW, "S");
        cp = mb;
        testStereoExample(cp.replace(mb.find("CFG=3"), 5, "     "), 1,
                          Atom::ChiralType::CHI_TETRAHEDRAL_CCW, "S");
      }
    }
  }
  SECTION("three coordinate") {
    auto m = R"CTAB(
  Mrv2108 01192209042D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 2.31 4.001 0 0
M  V30 2 C 1.54 2.6674 0 0 CFG=2
M  V30 3 O -0 2.6674 0 0
M  V30 4 C 2.31 1.3337 0 0 CFG=1
M  V30 5 F 3.85 1.3337 0 0
M  V30 6 C 1.54 0 0 0 CFG=2
M  V30 7 C 2.31 -1.3337 0 0
M  V30 8 O 0 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 2 3 CFG=3
M  V30 3 1 4 2
M  V30 4 1 4 5 CFG=1
M  V30 5 1 4 6
M  V30 6 1 6 7
M  V30 7 1 6 8 CFG=3
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
    CHECK(m->getAtomWithIdx(3)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
    CHECK(m->getAtomWithIdx(5)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
  }
#endif
  SECTION("this came up") {
    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase + "/Code/GraphMol/test_data/github87.mol";
    std::unique_ptr<RWMol> m{MolFileToMol(fName)};
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(0)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
  }

  SECTION("narrow angle") {
    auto m = R"CTAB(
  Mrv2211 06092305312D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.4875 5.1589 0 0
M  V30 2 C -1.5412 6.698 0 0
M  V30 3 F -2.2685 4.0178 0 0
M  V30 4 Br -0.7027 4.0359 0 0
M  V30 5 Cl -1.4337 3.6198 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 CFG=1
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(0)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
  }
  SECTION("track bond starting points") {
    auto m = R"CTAB(blah
     RDKit          2D

  7  7  0  0  0  0  0  0  0  0999 V2000
   -3.1960  -20.3395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1960  -19.5114    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1921  -21.7775    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7763  -21.0619    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0577  -21.4777    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -4.0179  -20.3565    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -2.3411  -20.6109    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  4  1  6
  1  6  1  0
  1  7  1  0
  4  3  1  1
  4  5  1  0
  4  7  1  0
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(0)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(3)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
  SECTION("linear arrangements") {
    {
      auto m = R"CTAB(
  Mrv2211 06102314502D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 0 0 1
M  V30 BEGIN ATOM
M  V30 1 O 9.0665 0.9156 0 0
M  V30 2 N 9.6304 4.5593 0 0
M  V30 3 C 8.2965 2.2493 0 0 MASS=14
M  V30 4 C 6.9628 1.4792 0 0
M  V30 5 C 9.6304 3.0193 0 0
M  V30 6 H 7.8191 3.0761 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 3 1 CFG=1
M  V30 2 1 2 5
M  V30 3 1 3 4
M  V30 4 1 3 5
M  V30 5 1 3 6
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(2)->getChiralTag() ==
            Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    }
    {
      auto m = R"CTAB(opposing stereo
  Mrv2211 06102314502D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 0 0 1
M  V30 BEGIN ATOM
M  V30 1 O 9.0665 0.9156 0 0
M  V30 2 N 9.6304 4.5593 0 0
M  V30 3 C 8.2965 2.2493 0 0 MASS=14
M  V30 4 C 6.9628 1.4792 0 0
M  V30 5 C 9.6304 3.0193 0 0
M  V30 6 H 7.8191 3.0761 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 3 1 CFG=1
M  V30 2 1 2 5
M  V30 3 1 3 4
M  V30 4 1 3 5
M  V30 5 1 3 6 CFG=3
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(2)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
    {
      // std::cerr<<"11111111111111"<<std::endl;
      auto m = R"CTAB(opposing stereo, order change
  Mrv2211 06102314502D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 0 0 1
M  V30 BEGIN ATOM
M  V30 1 O 9.0665 0.9156 0 0
M  V30 2 N 9.6304 4.5593 0 0
M  V30 3 C 8.2965 2.2493 0 0 MASS=14
M  V30 4 C 6.9628 1.4792 0 0
M  V30 5 C 9.6304 3.0193 0 0
M  V30 6 H 7.8191 3.0761 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 3 1 CFG=1
M  V30 2 1 2 5
M  V30 3 1 3 4
M  V30 4 1 3 6 CFG=3
M  V30 5 1 3 5
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(2)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
    {
      // IUPAC (ST-1.2.12) says this one is wrong. It definitely requires making
      // an assumption about where the H is.
      auto m = R"CTAB(three-coordinate, T shaped, wedge in the middle
  Mrv2211 06102314502D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 1
M  V30 BEGIN ATOM
M  V30 1 O 9.0665 0.9156 0 0
M  V30 2 N 9.6304 4.5593 0 0
M  V30 3 C 8.2965 2.2493 0 0 MASS=14
M  V30 4 C 6.9628 1.4792 0 0
M  V30 5 C 9.6304 3.0193 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 3 1 CFG=1
M  V30 2 1 2 5
M  V30 3 1 3 4
M  V30 4 1 3 5
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;

      REQUIRE(m);
      CHECK(m->getAtomWithIdx(2)->getChiralTag() ==
            Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    }
  }
}

TEST_CASE(
    "GitHub Issue #6502: MolToMolBlock writes \"either\" stereo for double bonds"
    "which shouldn't be stereo.",
    "[bug][molblock][stereo]") {
  auto m = "CP1(O)=NP(C)(O)=NP(C)(O)=NP(C)(O)=N1"_smiles;
  REQUIRE(m);

  auto mb = MolToV3KMolBlock(*m);
  CHECK(mb.find("CFG=2") == std::string::npos);
}

TEST_CASE("stereo in ring", "[molblock][stereo]") {
  SECTION("test 1") {
    auto molblock = R"CTAB(
  Mrv2311 10242314442D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 10 11 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C -2.6673 -0.77 0 0
M  V30 2 C -2.6673 0.77 0 0
M  V30 3 C -1.3337 1.54 0 0
M  V30 4 C 0 0.77 0 0
M  V30 5 C 1.3336 1.54 0 0
M  V30 6 C 2.6673 0.77 0 0
M  V30 7 C 2.6673 -0.77 0 0
M  V30 8 C 1.3336 -1.54 0 0
M  V30 9 C 0 -0.77 0 0
M  V30 10 C -1.3337 -1.54 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 7 1 7 8
M  V30 8 1 8 9
M  V30 9 1 4 9
M  V30 10 1 9 10
M  V30 11 1 1 10
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB";

    auto m = MolBlockToMol(molblock, true, false, false);

    REQUIRE(m);
    CHECK(m->getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREONONE);
  }
}
