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
}
