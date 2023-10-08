//
//  Copyright (C) 2023 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <string_view>
#include <streambuf>

#include "RDGeneral/test.h"
#include "catch.hpp"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/ProximityBonds.h>

using namespace RDKit;

TEST_CASE(
    "Github #6756: ConnectTheDots can segfault if atoms do not have residue info") {
  SECTION("small molecule, was not a problem") {
    std::string smiles =
        "[H]C([H])=C([H])C([H])([H])[H] |(1.35283,-0.241953,-1.28693;1.35436,-0.216135,-0.201646;2.30673,-0.352635,0.30154;0.23076,-0.0239947,0.500064;0.281106,-0.00576639,1.58653;-1.11716,0.171588,-0.112486;-1.08265,0.139951,-1.20609;-1.52627,1.14083,0.188471;-1.79971,-0.611884,0.230548)|";
    SmilesParserParams ps;
    ps.removeHs = false;
    std::unique_ptr<RWMol> m(SmilesToMol(smiles, ps));
    REQUIRE(m);
    CHECK(m->getNumBonds() == 8);
    // start without bonds
    m->beginBatchEdit();
    for (auto bond : m->bonds()) {
      m->removeBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
    }
    m->commitBatchEdit();
    CHECK(m->getNumBonds() == 0);
    ConnectTheDots(m.get());
    CHECK(m->getNumBonds() == 8);
  }
  SECTION("small molecule with bivalent H, this did crash") {
    std::string smiles = "F[H]F |(1.1,0,0;0,0,0;-1.1,0,0)|";
    SmilesParserParams ps;
    ps.removeHs = false;
    ps.sanitize = false;
    std::unique_ptr<RWMol> m(SmilesToMol(smiles, ps));
    REQUIRE(m);
    m->updatePropertyCache(false);
    CHECK(m->getNumBonds() == 2);
    // start without bonds
    m->beginBatchEdit();
    for (auto bond : m->bonds()) {
      m->removeBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
    }
    m->commitBatchEdit();
    CHECK(m->getNumBonds() == 0);
    ConnectTheDots(m.get());
    CHECK(m->getNumBonds() == 1);
  }
}