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
