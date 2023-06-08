//
//  Copyright (C) 2022 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RDGeneral/test.h"
#include "catch.hpp"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolWriters.h>

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
    WedgeMolBonds(*frag, &frag->getConformer());
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
    int flag;
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
    int flag;
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