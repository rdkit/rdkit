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
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

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

TEST_CASE("trans bonds in large rings") {
  // In large rings, we need to retain a trans geometry for double bonds.
  // This simulates the case where we write to SDF and read again.
  auto mol = "C1=C/CCCCCCCCCCCCC/1"_smiles;
  RDDepict::compute2DCoords(*mol);
  // simulate writing to SDF and reading again:
  RDKit::MolOps::removeStereochemistry(*mol);
  mol->getConformer().set3D(true);
  RDKit::MolOps::assignStereochemistryFrom3D(*mol);
  CHECK(RDKit::MolToSmiles(*mol) == "C1=C/CCCCCCCCCCCCC/1");
}

TEST_CASE("generate aligned coords accept failure") {
  auto template_ref_molblock = R"CTAB(
     RDKit          2D

  9  9  0  0  0  0  0  0  0  0999 V2000
   -0.8929    1.0942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8929   -1.9059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910    1.0942    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0
    1.7051    1.0942    0.0000 R2  0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910   -1.9059    0.0000 R3  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  6  8  1  0
  3  9  1  0
  2  7  1  0
M  RGP  3   7   1   8   2   9   3
M  END
)CTAB";
  std::unique_ptr<RWMol> template_ref(MolBlockToMol(template_ref_molblock));
  REQUIRE(template_ref);
  auto mol_molblock = R"CTAB(
     RDKit          2D

  9  9  0  0  0  0  0  0  0  0999 V2000
   -0.8929    1.0942    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8929   -1.9059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910    1.0942    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    1.7051    1.0942    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910   -1.9059    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  6  8  1  0
  3  9  1  0
  2  7  1  0
M  END
)CTAB";
  std::unique_ptr<RWMol> mol(MolBlockToMol(mol_molblock));
  REQUIRE(mol);
  SECTION("acceptFailure false, existing coords preserved") {
    REQUIRE_THROWS_AS(RDDepict::generateDepictionMatching2DStructure(
                          *mol, *template_ref, -1, nullptr, false, false, true),
                      RDDepict::DepictException);
    CHECK(MolToMolBlock(*mol) == mol_molblock);
  }
  SECTION("acceptFailure true, existing coords overwritten") {
    CHECK(RDDepict::generateDepictionMatching2DStructure(
              *mol, *template_ref, -1, nullptr, true, false, true)
              .empty());
    CHECK(MolToMolBlock(*mol) != mol_molblock);
  }
  SECTION("acceptFailure false, no existing coords") {
    mol->removeConformer(0);
    RDDepict::ConstrainedDepictionParams p;
    p.allowRGroups = true;
    p.acceptFailure = false;
    CHECK(mol->getNumConformers() == 0);
    REQUIRE_THROWS_AS(RDDepict::generateDepictionMatching2DStructure(
                          *mol, *template_ref, -1, nullptr, p),
                      RDDepict::DepictException);
    CHECK(mol->getNumConformers() == 0);
  }
  SECTION("acceptFailure true, no existing coords") {
    mol->removeConformer(0);
    RDDepict::ConstrainedDepictionParams p;
    p.allowRGroups = true;
    p.acceptFailure = true;
    CHECK(mol->getNumConformers() == 0);
    CHECK(RDDepict::generateDepictionMatching2DStructure(*mol, *template_ref,
                                                         -1, nullptr, p)
              .empty());
    CHECK(mol->getNumConformers() == 1);
  }
}

TEST_CASE("generate aligned coords alignOnly") {
  auto template_ref_molblock = R"CTAB(
     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
  -13.7477    6.3036    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  -13.7477    4.7567    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -12.6540    3.6628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.7477    2.5691    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -14.8414    3.6628    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  -11.1071    3.6628    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  2  5  1  0
  3  6  1  0
M  RGP  2   1   1   6   2
M  END
)CTAB";
  std::unique_ptr<RWMol> template_ref(MolBlockToMol(template_ref_molblock));
  REQUIRE(template_ref);
  auto mol_molblock = R"CTAB(
     RDKit          2D

 18 22  0  0  0  0  0  0  0  0999 V2000
    4.3922   -1.5699    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9211   -2.0479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5995   -0.5349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3731    0.8046    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8441    1.2825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0704   -0.0568    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8666    0.7748    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7736   -0.3197    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7749   -1.8666    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7718   -1.8679    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7731   -0.3208    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8679    0.7718    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0718   -0.0598    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3933   -1.5729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9222   -2.0509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6008   -0.5379    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3744    0.8016    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8454    1.2795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  9 10  1  0
 11 10  1  0
 11  8  1  0
  8  9  1  0
  4  5  1  0
  6  5  1  0
  7  6  1  0
  3  4  1  0
  3  7  1  0
  1  6  1  0
  2  3  1  0
  2  1  1  0
 17 18  1  0
 13 18  1  0
 12 13  1  0
 16 17  1  0
 16 12  1  0
 14 13  1  0
 15 16  1  0
 15 14  1  0
 12 11  1  0
  8  7  1  0
M  END
)CTAB";
  std::unique_ptr<RWMol> mol(MolBlockToMol(mol_molblock));
  REQUIRE(mol);
  auto bondLength11_12 =
      MolTransforms::getBondLength(mol->getConformer(), 11, 12);
  auto bondLength5_6 = MolTransforms::getBondLength(mol->getConformer(), 5, 6);
  REQUIRE(fabs(bondLength11_12 - bondLength5_6) < 1.e-4);
  REQUIRE(bondLength11_12 > 2.3);
  SECTION("alignOnly false/true") {
    for (auto alignOnly : {false, true}) {
      mol.reset(MolBlockToMol(mol_molblock));
      REQUIRE(mol);
      RDDepict::ConstrainedDepictionParams p;
      p.allowRGroups = true;
      p.alignOnly = alignOnly;
      auto res = RDDepict::generateDepictionMatching2DStructure(
          *mol, *template_ref, -1, nullptr, p);
      std::vector<int> expectedMolIndices{11, 10, 7, 8, 9, 6};
      auto sameIndices = std::all_of(
          res.begin(), res.end(), [&expectedMolIndices](const auto &pair) {
            return pair.second == expectedMolIndices.at(pair.first);
          });
      CHECK(sameIndices);
      CHECK(MolToSmiles(*mol) == "C1CC2CCC1N2C1CNC1N1C2CCC1CC2");
      auto samePositions = std::all_of(
          res.begin(), res.end(), [&mol, &template_ref](const auto &pair) {
            return (mol->getConformer().getAtomPos(pair.second) -
                    template_ref->getConformer().getAtomPos(pair.first))
                       .lengthSq() < 1.e-4;
          });
      CHECK(samePositions);
      auto bondLengthAli11_12 =
          MolTransforms::getBondLength(mol->getConformer(), 11, 12);
      auto bondLengthAli5_6 =
          MolTransforms::getBondLength(mol->getConformer(), 5, 6);
      CHECK(fabs(bondLengthAli11_12 - bondLengthAli5_6) < 1.e-4);
      if (alignOnly) {
        CHECK(bondLengthAli11_12 > 2.3);
      } else {
        CHECK(bondLengthAli11_12 < 1.6);
      }
    }
  }
}

TEST_CASE("generate aligned coords and wedging") {
  auto wedgedMol = R"CTAB(
     RDKit          2D

 29 34  0  0  1  0  0  0  0  0999 V2000
    1.3719    5.1304    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5985    3.7907    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9482    3.7907    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7216    5.1304    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2685    5.1304    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8994    3.5835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5597    4.3569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5597    5.9038    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8994    6.6771    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2389    5.9038    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5784    6.6771    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2389    4.3569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3719    2.4510    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5985    1.1115    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3719   -0.2276    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9188   -0.2276    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6921    1.1115    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9188    2.4510    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2389    1.1115    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0124   -0.2276    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2389   -1.5673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6921   -1.5673    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.8996   -5.0201    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2391   -4.2467    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5777   -6.5331    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9909   -5.9040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0124   -2.9070    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3306   -6.6772    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5784   -5.0201    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1
  2  3  1  0
  3  4  1  0
  5  4  1  6
  5  6  1  0
  6  7  1  0
  7  8  1  0
  9  8  1  1
  5  9  1  0
  9 10  1  0
 10 11  1  1
 10 12  1  0
  6 12  1  1
  2 13  1  0
 13 14  2  0
 14 15  1  0
 15 16  2  0
 16 17  1  0
 17 18  2  0
 13 18  1  0
 17 19  1  0
 19 20  1  0
 20 21  1  0
 21 22  1  0
 16 22  1  0
 23 24  1  0
 23 25  1  0
 25 26  1  0
 24 27  1  0
 27 26  1  0
 26 28  1  0
 24 29  1  0
 28 29  1  0
 21 27  1  0
M  END
)CTAB"_ctab;
  REQUIRE(wedgedMol);
  auto originalWedges = R"CTAB(  2  1  1  1
  2  3  1  0
  3  4  1  0
  5  4  1  6
  5  6  1  0
  6  7  1  0
  7  8  1  0
  9  8  1  1
  5  9  1  0
  9 10  1  0
 10 11  1  1
 10 12  1  0
  6 12  1  1
  2 13  1  0
 13 14  2  0
 14 15  1  0
 15 16  2  0
 16 17  1  0
 17 18  2  0
 13 18  1  0
 17 19  1  0
 19 20  1  0
 20 21  1  0
 21 22  1  0
 16 22  1  0
 23 24  1  0
 23 25  1  0
 25 26  1  0
 24 27  1  0
 27 26  1  0
 26 28  1  0
 24 29  1  0
 28 29  1  0
 21 27  1  0
M  END
)CTAB";
  auto invertedWedges = R"CTAB(  2  1  1  6
  2  3  1  0
  3  4  1  0
  5  4  1  1
  5  6  1  0
  6  7  1  0
  7  8  1  0
  9  8  1  6
  5  9  1  0
  9 10  1  0
 10 11  1  6
 10 12  1  0
  6 12  1  6
  2 13  1  0
 13 14  2  0
 14 15  1  0
 15 16  2  0
 16 17  1  0
 17 18  2  0
 13 18  1  0
 17 19  1  0
 19 20  1  0
 20 21  1  0
 21 22  1  0
 16 22  1  0
 23 24  1  0
 23 25  1  0
 25 26  1  0
 24 27  1  0
 27 26  1  0
 26 28  1  0
 24 29  1  0
 28 29  1  0
 21 27  1  0
)CTAB";
  SECTION("wedging all within scaffold") {
    auto scaffold = R"CTAB(
     RDKit          2D

 13 14  0  0  1  0  0  0  0  0999 V2000
   -1.6549    2.5755    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8814    1.2358    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6653    1.2358    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4385    2.5755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9854    2.5755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6161    1.0286    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2766    1.8019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2766    3.3487    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6161    4.1222    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9558    3.3487    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.2953    4.1222    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    4.9558    1.8019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6549   -0.1037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0
  2  3  1  0
  3  4  1  0
  5  4  1  1
  5  6  1  0
  6  7  1  6
  7  8  1  0
  9  8  1  6
  5  9  1  0
  9 10  1  0
 10 11  1  6
 10 12  1  0
  6 12  1  0
  2 13  1  6
M  END
)CTAB"_ctab;
    REQUIRE(scaffold);
    // the "alignOnly" alignment should succeed and preserve molblock wedging
    // (inverted with respect to the original molecule)
    // it should feature a narrow angle between the bridge bonds
    // as the original geometry of the bridge is preserved
    {
      ROMol wedgedMolCopy(*wedgedMol);
      RDDepict::ConstrainedDepictionParams p;
      p.alignOnly = true;
      REQUIRE(!RDDepict::generateDepictionMatching2DStructure(
                   wedgedMolCopy, *scaffold, -1, nullptr, p)
                   .empty());
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto angle =
          MolTransforms::getAngleDeg(wedgedMolCopy.getConformer(), 23, 26, 25);
      CHECK((angle > 10. && angle < 15.));
      auto mbAlignOnly = MolToMolBlock(wedgedMolCopy);
      CHECK(mbAlignOnly.find(invertedWedges) != std::string::npos);
    }
    // the "alignOnly" alignment should succeed and preserve molblock wedging
    // (same as original molecule)
    // it should feature a narrow angle between the bridge bonds
    // as the original geometry of the bridge is preserved
    {
      ROMol wedgedMolCopy(*wedgedMol);
      RDDepict::ConstrainedDepictionParams p;
      p.alignOnly = true;
      p.adjustMolBlockWedging = false;
      REQUIRE(!RDDepict::generateDepictionMatching2DStructure(
                   wedgedMolCopy, *scaffold, -1, nullptr, p)
                   .empty());
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto angle =
          MolTransforms::getAngleDeg(wedgedMolCopy.getConformer(), 23, 26, 25);
      CHECK((angle > 10. && angle < 15.));
      auto mbAlignOnly = MolToMolBlock(wedgedMolCopy);
      CHECK(mbAlignOnly.find(originalWedges) != std::string::npos);
    }
    // the "rebuild" alignment should succeed and preserve molblock wedging
    // (inverted with respect to the original molecule)
    // it should feature a much wider angle between the bridge bonds as the
    // bridged system is entirely rebuilt since it is not part of the scaffold
    {
      ROMol wedgedMolCopy(*wedgedMol);
      REQUIRE(!RDDepict::generateDepictionMatching2DStructure(wedgedMolCopy,
                                                              *scaffold)
                   .empty());
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto angle =
          MolTransforms::getAngleDeg(wedgedMolCopy.getConformer(), 23, 26, 25);
      CHECK((angle > 105. && angle < 110.));
      auto mbAlignOnly = MolToMolBlock(wedgedMolCopy);
      CHECK(mbAlignOnly.find(invertedWedges) != std::string::npos);
    }
    // the "rebuild" alignment should succeed and preserve molblock wedging
    // (same as the original molecule)
    // it should feature a much wider angle between the bridge bonds as the
    // bridged system is entirely rebuilt since it is not part of the scaffold
    {
      ROMol wedgedMolCopy(*wedgedMol);
      RDDepict::ConstrainedDepictionParams p;
      p.adjustMolBlockWedging = false;
      REQUIRE(!RDDepict::generateDepictionMatching2DStructure(
                   wedgedMolCopy, *scaffold, -1, nullptr, p)
                   .empty());
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto angle =
          MolTransforms::getAngleDeg(wedgedMolCopy.getConformer(), 23, 26, 25);
      CHECK((angle > 105. && angle < 110.));
      auto mbAlignOnly = MolToMolBlock(wedgedMolCopy);
      CHECK(mbAlignOnly.find(originalWedges) != std::string::npos);
    }
#ifdef RDK_BUILD_COORDGEN_SUPPORT
    // the "rebuildCoordGen" alignment should succeed and clear original wedging
    // it should feature an even wider angle between the bridge bonds as
    // CoordGen has a template for the bridged system. Additionally, CoordGen
    // also rebuilds the scaffold, therefore original wedging should be cleared
    RDDepict::preferCoordGen = true;
    {
      ROMol wedgedMolCopy(*wedgedMol);
      REQUIRE(!RDDepict::generateDepictionMatching2DStructure(wedgedMolCopy,
                                                              *scaffold)
                   .empty());
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto angle =
          MolTransforms::getAngleDeg(wedgedMolCopy.getConformer(), 23, 26, 25);
      CHECK((angle > 145. && angle < 150.));
      auto mbAlignOnly = MolToMolBlock(wedgedMolCopy);
      CHECK(mbAlignOnly.find(invertedWedges) == std::string::npos);
      CHECK(mbAlignOnly.find(originalWedges) == std::string::npos);
    }
    // the "rebuildCoordGen" alignment should succeed and keep original wedging
    // unaltered.
    // it should feature an even wider angle between the bridge bonds as
    // CoordGen has a template for the bridged system.
    {
      ROMol wedgedMolCopy(*wedgedMol);
      RDDepict::ConstrainedDepictionParams p;
      p.adjustMolBlockWedging = false;
      REQUIRE(!RDDepict::generateDepictionMatching2DStructure(
                   wedgedMolCopy, *scaffold, -1, nullptr, p)
                   .empty());
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto angle =
          MolTransforms::getAngleDeg(wedgedMolCopy.getConformer(), 23, 26, 25);
      CHECK((angle > 145. && angle < 150.));
      auto mbAlignOnly = MolToMolBlock(wedgedMolCopy);
      CHECK(mbAlignOnly.find(originalWedges) != std::string::npos);
    }
    RDDepict::preferCoordGen = false;
#endif
  }
  SECTION("wedging outside scaffold") {
    auto scaffold = R"CTAB(
     RDKit          2D

  9 10  0  0  1  0  0  0  0  0999 V2000
   -0.8816    0.5663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6651    0.5663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2958   -0.9804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0435   -0.2072    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0435    1.3395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2958    2.1129    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6355    1.3395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9750    2.1129    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    2.6355   -0.2072    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1
  2  3  1  0
  3  4  1  6
  4  5  1  0
  6  5  1  6
  2  6  1  0
  6  7  1  0
  7  8  1  6
  7  9  1  0
  3  9  1  0
M  END
)CTAB"_ctab;
    REQUIRE(scaffold);
    // the "alignOnly" alignment should succeed and preserve molblock wedging
    // (inverted with respect to the original molecule)
    // it should feature a narrow angle between the bridge bonds
    // as the original geometry of the bridge is preserved
    {
      ROMol wedgedMolCopy(*wedgedMol);
      RDDepict::ConstrainedDepictionParams p;
      p.alignOnly = true;
      REQUIRE(!RDDepict::generateDepictionMatching2DStructure(
                   wedgedMolCopy, *scaffold, -1, nullptr, p)
                   .empty());
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto angle =
          MolTransforms::getAngleDeg(wedgedMolCopy.getConformer(), 23, 26, 25);
      CHECK((angle > 10. && angle < 15.));
      auto mbAlignOnly = MolToMolBlock(wedgedMolCopy);
      CHECK(mbAlignOnly.find(invertedWedges) != std::string::npos);
    }
    // the "alignOnly" alignment should succeed and preserve molblock wedging
    // (same as original molecule)
    // it should feature a narrow angle between the bridge bonds
    // as the original geometry of the bridge is preserved
    {
      ROMol wedgedMolCopy(*wedgedMol);
      RDDepict::ConstrainedDepictionParams p;
      p.alignOnly = true;
      p.adjustMolBlockWedging = false;
      REQUIRE(!RDDepict::generateDepictionMatching2DStructure(
                   wedgedMolCopy, *scaffold, -1, nullptr, p)
                   .empty());
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto angle =
          MolTransforms::getAngleDeg(wedgedMolCopy.getConformer(), 23, 26, 25);
      CHECK((angle > 10. && angle < 15.));
      auto mbAlignOnly = MolToMolBlock(wedgedMolCopy);
      CHECK(mbAlignOnly.find(originalWedges) != std::string::npos);
    }
    // the "rebuild" alignment should succeed and clear original wedging
    // it should feature a much wider angle between the bridge bonds as the
    // bridged system is entirely rebuilt since it is not part of the scaffold
    {
      ROMol wedgedMolCopy(*wedgedMol);
      REQUIRE(!RDDepict::generateDepictionMatching2DStructure(wedgedMolCopy,
                                                              *scaffold)
                   .empty());
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto angle =
          MolTransforms::getAngleDeg(wedgedMolCopy.getConformer(), 23, 26, 25);
      CHECK((angle > 105. && angle < 110.));
      auto mbAlignOnly = MolToMolBlock(wedgedMolCopy);
      CHECK(mbAlignOnly.find(originalWedges) == std::string::npos);
      CHECK(mbAlignOnly.find(invertedWedges) == std::string::npos);
    }
    // the "rebuild" alignment should succeed and preserve molblock wedging
    // (same as the original molecule)
    // it should feature a much wider angle between the bridge bonds as the
    // bridged system is entirely rebuilt since it is not part of the scaffold
    {
      ROMol wedgedMolCopy(*wedgedMol);
      RDDepict::ConstrainedDepictionParams p;
      p.adjustMolBlockWedging = false;
      REQUIRE(!RDDepict::generateDepictionMatching2DStructure(
                   wedgedMolCopy, *scaffold, -1, nullptr, p)
                   .empty());
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto angle =
          MolTransforms::getAngleDeg(wedgedMolCopy.getConformer(), 23, 26, 25);
      CHECK((angle > 105. && angle < 110.));
      auto mbAlignOnly = MolToMolBlock(wedgedMolCopy);
      CHECK(mbAlignOnly.find(originalWedges) != std::string::npos);
    }
#ifdef RDK_BUILD_COORDGEN_SUPPORT
    // the "rebuildCoordGen" alignment should succeed and clear original wedging
    // it should feature an even wider angle between the bridge bonds as
    // CoordGen has a template for the bridged system. Additionally, CoordGen
    // also rebuilds the scaffold, therefore original wedging should be cleared
    RDDepict::preferCoordGen = true;
    {
      ROMol wedgedMolCopy(*wedgedMol);
      REQUIRE(!RDDepict::generateDepictionMatching2DStructure(wedgedMolCopy,
                                                              *scaffold)
                   .empty());
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto angle =
          MolTransforms::getAngleDeg(wedgedMolCopy.getConformer(), 23, 26, 25);
      CHECK((angle > 145. && angle < 150.));
      auto mbAlignOnly = MolToMolBlock(wedgedMolCopy);
      CHECK(mbAlignOnly.find(invertedWedges) == std::string::npos);
      CHECK(mbAlignOnly.find(originalWedges) == std::string::npos);
    }
    // the "rebuildCoordGen" alignment should succeed and keep original wedging
    // unaltered.
    // it should feature an even wider angle between the bridge bonds as
    // CoordGen has a template for the bridged system.
    {
      ROMol wedgedMolCopy(*wedgedMol);
      RDDepict::ConstrainedDepictionParams p;
      p.adjustMolBlockWedging = false;
      REQUIRE(!RDDepict::generateDepictionMatching2DStructure(
                   wedgedMolCopy, *scaffold, -1, nullptr, p)
                   .empty());
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto angle =
          MolTransforms::getAngleDeg(wedgedMolCopy.getConformer(), 23, 26, 25);
      CHECK((angle > 145. && angle < 150.));
      auto mbAlignOnly = MolToMolBlock(wedgedMolCopy);
      CHECK(mbAlignOnly.find(originalWedges) != std::string::npos);
    }
    RDDepict::preferCoordGen = false;
#endif
  }
  SECTION("wedging no match") {
    auto scaffoldNoMatch = R"CTAB(
     RDKit          2D

 13 14  0  0  1  0  0  0  0  0999 V2000
   -1.6549    2.5755    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8814    1.2358    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6653    1.2358    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4385    2.5755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9854    2.5755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6161    1.0286    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2766    1.8019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2766    3.3487    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6161    4.1222    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9558    3.3487    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.2953    4.1222    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    4.9558    1.8019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6549   -0.1037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0
  2  3  1  0
  3  4  1  0
  5  4  1  1
  5  6  1  0
  6  7  1  6
  7  8  1  0
  9  8  1  6
  5  9  1  0
  9 10  1  0
 10 11  1  6
 10 12  1  0
  6 12  1  0
  2 13  1  6
M  END
)CTAB"_ctab;
    REQUIRE(scaffoldNoMatch);
    std::string origMolBlock;
    {
      ROMol wedgedMolCopy(*wedgedMol);
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      origMolBlock = MolToMolBlock(wedgedMolCopy);
    }
    REQUIRE(!origMolBlock.empty());
    // the "alignOnly" alignment should throw if acceptFailure is false
    // and preserve the original coordinates
    {
      ROMol wedgedMolCopy(*wedgedMol);
      RDDepict::ConstrainedDepictionParams p;
      p.alignOnly = true;
      REQUIRE_THROWS_AS(RDDepict::generateDepictionMatching2DStructure(
                            wedgedMolCopy, *scaffoldNoMatch, -1, nullptr, p),
                        RDDepict::DepictException);
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto currMolBlock = MolToMolBlock(wedgedMolCopy);
      CHECK(currMolBlock == origMolBlock);
      CHECK(currMolBlock.find(invertedWedges) == std::string::npos);
    }
    // the "alignOnly" alignment should return an empty MatchVect if
    // acceptFailure is true and generate new coordinates, hence wedging should
    // be cleared
    {
      ROMol wedgedMolCopy(*wedgedMol);
      RDDepict::ConstrainedDepictionParams p;
      p.alignOnly = true;
      p.acceptFailure = true;
      REQUIRE(RDDepict::generateDepictionMatching2DStructure(
                  wedgedMolCopy, *scaffoldNoMatch, -1, nullptr, p)
                  .empty());
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto currMolBlock = MolToMolBlock(wedgedMolCopy);
      CHECK(currMolBlock != origMolBlock);
      CHECK(currMolBlock.find(invertedWedges) == std::string::npos);
    }
    // the "rebuild" alignment should throw if acceptFailure is false
    // and preserve the original coordinates
    {
      ROMol wedgedMolCopy(*wedgedMol);
      REQUIRE_THROWS_AS(RDDepict::generateDepictionMatching2DStructure(
                            wedgedMolCopy, *scaffoldNoMatch),
                        RDDepict::DepictException);
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto currMolBlock = MolToMolBlock(wedgedMolCopy);
      CHECK(currMolBlock == origMolBlock);
      CHECK(currMolBlock.find(invertedWedges) == std::string::npos);
    }
    // the "rebuild" alignment should return an empty MatchVect if acceptFailure
    // is true and generate new coordinates, hence wedging should be cleared
    {
      ROMol wedgedMolCopy(*wedgedMol);
      RDDepict::ConstrainedDepictionParams p;
      p.acceptFailure = true;
      REQUIRE(RDDepict::generateDepictionMatching2DStructure(
                  wedgedMolCopy, *scaffoldNoMatch, -1, nullptr, p)
                  .empty());
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto currMolBlock = MolToMolBlock(wedgedMolCopy);
      CHECK(currMolBlock != origMolBlock);
      CHECK(currMolBlock.find(invertedWedges) == std::string::npos);
    }
#ifdef RDK_BUILD_COORDGEN_SUPPORT
    // the "rebuildCoordGen" alignment should throw if acceptFailure is false
    // and preserve the original coordinates
    RDDepict::preferCoordGen = true;
    {
      ROMol wedgedMolCopy(*wedgedMol);
      REQUIRE_THROWS_AS(RDDepict::generateDepictionMatching2DStructure(
                            wedgedMolCopy, *scaffoldNoMatch),
                        RDDepict::DepictException);
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto currMolBlock = MolToMolBlock(wedgedMolCopy);
      CHECK(currMolBlock == origMolBlock);
      CHECK(currMolBlock.find(invertedWedges) == std::string::npos);
    }
    // the "rebuildCoordGen" alignment should return an empty MatchVect if
    // acceptFailure is true and generate new coordinates, hence wedging should
    // be cleared
    {
      ROMol wedgedMolCopy(*wedgedMol);
      RDDepict::ConstrainedDepictionParams p;
      p.acceptFailure = true;
      REQUIRE(RDDepict::generateDepictionMatching2DStructure(
                  wedgedMolCopy, *scaffoldNoMatch, -1, nullptr, p)
                  .empty());
      Chirality::reapplyMolBlockWedging(wedgedMolCopy);
      auto currMolBlock = MolToMolBlock(wedgedMolCopy);
      CHECK(currMolBlock != origMolBlock);
      CHECK(currMolBlock.find(invertedWedges) == std::string::npos);
    }
    RDDepict::preferCoordGen = false;
#endif
  }
}

TEST_CASE("generate aligned coords R group match") {
  auto templateRef = R"CTAB(
  MJ201100                      

  7  7  0  0  0  0  0  0  0  0999 V2000
   -0.5804    1.2045    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2948    0.7920    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2948   -0.0330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5804   -0.4455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1340   -0.0330    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
    0.1340    0.7920    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8485   -0.4455    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  6  1  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  5  7  1  0  0  0  0
M  RGP  1   7   1
M  END
)CTAB"_ctab;
  REQUIRE(templateRef);
  RDDepict::ConstrainedDepictionParams p;
  p.allowRGroups = true;
  SECTION("heavy") {
    for (auto alignOnly : {true, false}) {
      p.alignOnly = alignOnly;
      auto mol = "Cc1ccccc1"_smiles;
      REQUIRE(mol);
      REQUIRE(mol->getNumAtoms() == 7);
      CHECK(RDDepict::generateDepictionMatching2DStructure(*mol, *templateRef,
                                                           -1, nullptr, p)
                .size() == 7);
    }
  }
  SECTION("implicit hydrogen") {
    for (auto alignOnly : {true, false}) {
      p.alignOnly = alignOnly;
      auto mol = "c1ccccc1"_smiles;
      REQUIRE(mol);
      REQUIRE(mol->getNumAtoms() == 6);
      CHECK(RDDepict::generateDepictionMatching2DStructure(*mol, *templateRef,
                                                           -1, nullptr, p)
                .size() == 6);
    }
  }
  SECTION("explicit hydrogen") {
    auto smi = "[H]c1ccccc1";
    SmilesParserParams smilesParams;
    smilesParams.removeHs = false;
    for (auto alignOnly : {true, false}) {
      p.alignOnly = alignOnly;
      std::unique_ptr<RWMol> mol(SmilesToMol(smi, smilesParams));
      REQUIRE(mol);
      REQUIRE(mol->getNumAtoms() == 7);
      CHECK(RDDepict::generateDepictionMatching2DStructure(*mol, *templateRef,
                                                           -1, nullptr, p)
                .size() == 7);
    }
  }
  SECTION("no atom") {
    for (auto alignOnly : {true, false}) {
      p.alignOnly = alignOnly;
      auto mol = "n1ccccc1"_smiles;
      REQUIRE(mol);
      REQUIRE(mol->getNumAtoms() == 6);
      CHECK(RDDepict::generateDepictionMatching2DStructure(*mol, *templateRef,
                                                           -1, nullptr, p)
                .size() == 6);
    }
  }
  SECTION("charged") {
    for (auto alignOnly : {true, false}) {
      p.alignOnly = alignOnly;
      auto mol = "C[n+]1ccccc1"_smiles;
      REQUIRE(mol);
      REQUIRE(mol->getNumAtoms() == 7);
      CHECK(RDDepict::generateDepictionMatching2DStructure(*mol, *templateRef,
                                                           -1, nullptr, p)
                .size() == 7);
    }
  }
}
