//
//  Copyright (C) 2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>
#include "RDDepictor.h"
#include "DepictUtils.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>

using namespace RDKit;

TEST_CASE(
    "github #4504: overlapping coordinates with 1,1-disubstituted "
    "cyclobutanes") {
  SECTION("basics") {
    auto m = "CCC1(CCC1)CC1CCCCC1"_smiles;
    REQUIRE(m);
    CHECK(RDDepict::compute2DCoords(*m) == 0);
    auto &conf = m->getConformer();
    auto v = conf.getAtomPos(1) - conf.getAtomPos(3);
    CHECK(v.length() > 0.1);
    v = conf.getAtomPos(1) - conf.getAtomPos(5);
    CHECK(v.length() > 0.1);
  }
  SECTION("this one was ok") {
    auto m = "CCC1(CCC1)C1CCCCC1"_smiles;
    REQUIRE(m);
    CHECK(RDDepict::compute2DCoords(*m) == 0);
    auto &conf = m->getConformer();
    auto v = conf.getAtomPos(1) - conf.getAtomPos(3);
    CHECK(v.length() > 0.1);
    v = conf.getAtomPos(1) - conf.getAtomPos(5);
    CHECK(v.length() > 0.1);
  }
}

TEST_CASE("square planar", "[nontetrahedral]") {
  SECTION("cis-platin") {
    auto m = "Cl[Pt@SP1](Cl)(<-[NH3])<-[NH3]"_smiles;
    REQUIRE(m);
    CHECK(RDDepict::compute2DCoords(*m) == 0);
    auto &conf = m->getConformer();
    auto v1 = conf.getAtomPos(0) - conf.getAtomPos(2);
    auto v2 = conf.getAtomPos(0) - conf.getAtomPos(3);
    CHECK(v1.length() < v2.length());
  }
  SECTION("trans-platin") {
    auto m = "Cl[Pt@SP2](Cl)(<-[NH3])<-[NH3]"_smiles;
    REQUIRE(m);
    CHECK(RDDepict::compute2DCoords(*m) == 0);
    auto &conf = m->getConformer();
    auto v1 = conf.getAtomPos(0) - conf.getAtomPos(2);
    auto v2 = conf.getAtomPos(0) - conf.getAtomPos(3);
    CHECK(v1.length() > v2.length());
  }
  SECTION("trans-metal in a ring") {
    auto m = "C1[Pt@SP2](CCC1)(<-[NH3])<-[NH3]"_smiles;
    REQUIRE(m);
    CHECK(RDDepict::compute2DCoords(*m) == 0);
    // std::cerr << MolToV3KMolBlock(*m) << std::endl;
    auto &conf = m->getConformer();
    auto v1 = conf.getAtomPos(0) - conf.getAtomPos(2);
    auto v2 = conf.getAtomPos(0) - conf.getAtomPos(5);
    CHECK(v1.length() > v2.length());
  }
  SECTION("cis-metal in a ring") {
    auto m = "C1[Pt@SP1](CCC1)(<-[NH3])<-[NH3]"_smiles;
    REQUIRE(m);
    CHECK(RDDepict::compute2DCoords(*m) == 0);
    // std::cerr << MolToV3KMolBlock(*m) << std::endl;
    auto &conf = m->getConformer();
    auto v1 = conf.getAtomPos(0) - conf.getAtomPos(2);
    auto v2 = conf.getAtomPos(0) - conf.getAtomPos(5);
    CHECK(v1.length() < v2.length());
  }
}

TEST_CASE("trigonal bipyramidal", "[nontetrahedral]") {
  SECTION("TB1") {
    auto m = "S[As@TB1](F)(Cl)(Br)N"_smiles;
    REQUIRE(m);
    CHECK(RDDepict::compute2DCoords(*m) == 0);
    auto &conf = m->getConformer();
    auto v1 = conf.getAtomPos(0) - conf.getAtomPos(5);  // ax - ax
    auto v2 = conf.getAtomPos(0) - conf.getAtomPos(4);  // ax - eq
    auto v3 = conf.getAtomPos(2) - conf.getAtomPos(4);  // eq - eq long
    auto v4 = conf.getAtomPos(2) - conf.getAtomPos(3);  // eq - eq short
    CHECK(v1.length() > v2.length());
    CHECK(v1.length() > v3.length());
    CHECK(v3.length() > v2.length());
    CHECK(v3.length() > v4.length());
    CHECK(v2.length() > v4.length());
  }
  SECTION("TB3") {
    auto m = "S[As@TB3](F)(Cl)(N)Br"_smiles;
    REQUIRE(m);
    CHECK(RDDepict::compute2DCoords(*m) == 0);
    auto &conf = m->getConformer();
    auto v1 = conf.getAtomPos(0) - conf.getAtomPos(4);  // ax - ax
    auto v2 = conf.getAtomPos(0) - conf.getAtomPos(5);  // ax - eq
    auto v3 = conf.getAtomPos(2) - conf.getAtomPos(5);  // eq - eq long
    auto v4 = conf.getAtomPos(2) - conf.getAtomPos(3);  // eq - eq short
    CHECK(v1.length() > v2.length());
    CHECK(v1.length() > v3.length());
    CHECK(v3.length() > v2.length());
    CHECK(v3.length() > v4.length());
    CHECK(v2.length() > v4.length());
  }
  SECTION("TB1 missing ax") {
    auto m = "S[As@TB1](F)(Cl)Br"_smiles;
    REQUIRE(m);

    CHECK_THAT(
        Chirality::getIdealAngleBetweenLigands(
            m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(2)),
        Catch::Matchers::WithinAbs(90, 0.001));
    CHECK_THAT(
        Chirality::getIdealAngleBetweenLigands(
            m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(3)),
        Catch::Matchers::WithinAbs(90, 0.001));
    CHECK_THAT(
        Chirality::getIdealAngleBetweenLigands(
            m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(4)),
        Catch::Matchers::WithinAbs(90, 0.001));
    CHECK_THAT(
        Chirality::getIdealAngleBetweenLigands(
            m->getAtomWithIdx(1), m->getAtomWithIdx(2), m->getAtomWithIdx(3)),
        Catch::Matchers::WithinAbs(120, 0.001));

    CHECK(RDDepict::compute2DCoords(*m) == 0);
    auto &conf = m->getConformer();
    auto v2 = conf.getAtomPos(0) - conf.getAtomPos(4);  // ax - eq
    auto v3 = conf.getAtomPos(2) - conf.getAtomPos(4);  // eq - eq long
    auto v4 = conf.getAtomPos(2) - conf.getAtomPos(3);  // eq - eq short
    CHECK(v3.length() > v2.length());
    CHECK(v3.length() > v4.length());
    CHECK(v2.length() > v4.length());
  }
}

TEST_CASE("octahedral", "[nontetrahedral]") {
  SECTION("OH1") {
    auto m = "O[Co@OH1](Cl)(C)(N)(F)P"_smiles;
    REQUIRE(m);
    CHECK(RDDepict::compute2DCoords(*m) == 0);
    // std::cerr << MolToV3KMolBlock(*m) << std::endl;
    auto &conf = m->getConformer();
    auto v1 = conf.getAtomPos(3) - conf.getAtomPos(5);  // ax - ax
    auto v2 = conf.getAtomPos(3) - conf.getAtomPos(4);  // ax - eq
    auto v3 = conf.getAtomPos(3) - conf.getAtomPos(0);  // ax - eq
    auto v4 = conf.getAtomPos(0) - conf.getAtomPos(4);  // eq - eq nbr
    auto v5 = conf.getAtomPos(0) - conf.getAtomPos(6);  // eq - eq cross
    auto v6 = conf.getAtomPos(0) - conf.getAtomPos(2);  // eq - eq longnbr

    CHECK(v1.length() > v2.length());
    CHECK(v1.length() > v3.length());
    CHECK(v1.length() > v4.length());
    CHECK(v1.length() > v6.length());
    CHECK(v5.length() > v4.length());
    CHECK(v5.length() > v6.length());
  }
  SECTION("OH3") {
    auto m = "O[Co@OH3](Cl)(C)(N)(P)F"_smiles;
    REQUIRE(m);
    CHECK(RDDepict::compute2DCoords(*m) == 0);
    // std::cerr << MolToV3KMolBlock(*m) << std::endl;
    auto &conf = m->getConformer();
    auto v1 = conf.getAtomPos(3) - conf.getAtomPos(6);  // ax - ax
    auto v2 = conf.getAtomPos(3) - conf.getAtomPos(4);  // ax - eq
    auto v3 = conf.getAtomPos(3) - conf.getAtomPos(0);  // ax - eq
    auto v4 = conf.getAtomPos(0) - conf.getAtomPos(4);  // eq - eq nbr
    auto v5 = conf.getAtomPos(0) - conf.getAtomPos(5);  // eq - eq cross
    auto v6 = conf.getAtomPos(0) - conf.getAtomPos(2);  // eq - eq longnbr

    CHECK(v1.length() > v2.length());
    CHECK(v1.length() > v3.length());
    CHECK(v1.length() > v4.length());
    CHECK(v1.length() > v6.length());
    CHECK(v5.length() > v4.length());
    CHECK(v5.length() > v6.length());
  }
  SECTION("OH1 missing one ligand") {
    auto m = "O[Co@OH1](Cl)(C)(N)F"_smiles;
    REQUIRE(m);
    CHECK(RDDepict::compute2DCoords(*m) == 0);
    // std::cerr << MolToV3KMolBlock(*m) << std::endl;
    auto &conf = m->getConformer();
    auto v1 = conf.getAtomPos(3) - conf.getAtomPos(5);  // ax - ax
    auto v2 = conf.getAtomPos(3) - conf.getAtomPos(4);  // ax - eq
    auto v3 = conf.getAtomPos(3) - conf.getAtomPos(0);  // ax - eq
    auto v4 = conf.getAtomPos(0) - conf.getAtomPos(4);  // eq - eq nbr
    auto v6 = conf.getAtomPos(0) - conf.getAtomPos(2);  // eq - eq longnbr

    CHECK(v1.length() > v2.length());
    CHECK(v1.length() > v3.length());
    CHECK(v1.length() > v4.length());
    CHECK(v1.length() > v6.length());
  }
}

TEST_CASE("use ring system templates") {
  auto mol = "C1CCC2C(C1)C1CCN2NN1"_smiles;
  RDDepict::Compute2DCoordParameters params;
  RDDepict::compute2DCoords(*mol, params);
  auto diff =
      mol->getConformer().getAtomPos(10) - mol->getConformer().getAtomPos(11);
  // when templates are not used, bond from 10-11 is very short
  TEST_ASSERT(RDKit::feq(diff.length(), 0.116, .1));

  params.useRingTemplates = true;
  RDDepict::compute2DCoords(*mol, params);
  diff =
      mol->getConformer().getAtomPos(10) - mol->getConformer().getAtomPos(11);
  TEST_ASSERT(RDKit::feq(diff.length(), 1.0, .1))
}

TEST_CASE("dative bonds and rings") {
  auto mol = "O->[Pt]1(<-O)<-NC2CCC2N->1"_smiles;
  REQUIRE(mol);
  auto rings = mol->getRingInfo();
  CHECK(rings->numRings() == 1);  // the dative bonds are ignored
  RDDepict::compute2DCoords(*mol);
  CHECK(rings->numRings() == 1);  // ensure the ring count hasn't changed
  auto conf = mol->getConformer();
  auto v1 = conf.getAtomPos(1) - conf.getAtomPos(3);
  auto v2 = conf.getAtomPos(1) - conf.getAtomPos(8);
  CHECK_THAT(v1.length(), Catch::Matchers::WithinAbs(v2.length(), 0.01));
}

TEST_CASE("vicinal R groups can match an aromatic ring") {
  auto benzene = R"CTAB(
  MJ201100                      

  8  8  0  0  0  0  0  0  0  0999 V2000
   -1.0263   -0.3133    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -2.4553    0.5116    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -1.7408   -0.7258    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7408   -1.5509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4553   -1.9633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1698   -1.5509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1698   -0.7258    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4553   -0.3133    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  3  1  1  0  0  0  0
  8  2  1  0  0  0  0
  4  3  2  0  0  0  0
  5  4  1  0  0  0  0
  6  5  2  0  0  0  0
  7  6  1  0  0  0  0
  8  3  1  0  0  0  0
  8  7  2  0  0  0  0
M  RGP  2   1   2   2   1
M  END)CTAB"_ctab;
  REQUIRE(benzene);
  auto biphenyl = R"CTAB(
  MJ201100                      

 14 15  0  0  0  0  0  0  0  0999 V2000
   -0.6027    2.4098    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3171    1.9973    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3171    1.1722    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6027    0.7597    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1117    1.1722    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1117    1.9973    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6027   -0.0652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3171   -0.4777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3171   -1.3028    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6028   -1.7153    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1117   -1.3028    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1117   -0.4777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0316    0.7597    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -2.0316   -0.0652    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  4  7  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  2  0  0  0  0
 10 11  1  0  0  0  0
 11 12  2  0  0  0  0
  7  8  2  0  0  0  0
 12  7  1  0  0  0  0
  3 13  1  0  0  0  0
  8 14  1  0  0  0  0
M  RGP  2  13   1  14   2
M  END)CTAB"_ctab;
  REQUIRE(biphenyl);
  SECTION("R groups on benzene match quinoxaline") {
    auto quinoxaline = "c1ccc2nccnc2c1"_smiles;
    REQUIRE(quinoxaline);
    auto match = RDDepict::generateDepictionMatching2DStructure(
      *quinoxaline, *benzene, -1, nullptr, false, false, true);
    CHECK(match.size() == 8);
  }
  SECTION("R groups on benzene match tetralin") {
    auto tetralin = "c1cccc2CCCCc12"_smiles;
    REQUIRE(tetralin);
    auto match = RDDepict::generateDepictionMatching2DStructure(
      *tetralin, *benzene, -1, nullptr, false, false, true);
    CHECK(match.size() == 8);
  }
  SECTION("R groups on biphenyl match phenantridine") {
    auto phenantridine = "c1cccc2ncc3ccccc3c12"_smiles;
    REQUIRE(phenantridine);
    auto match = RDDepict::generateDepictionMatching2DStructure(
      *phenantridine, *biphenyl, -1, nullptr, false, false, true);
    CHECK(match.size() == 14);
  }
}
