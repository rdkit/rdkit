//
//  Copyright (C) 2023 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>
#ifdef RDK_BUILD_THREADSAFE_SSS
#include <future>
#include <thread>
#endif

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

using namespace RDKit::v2;

TEST_CASE("v2 basics") {
  {
    auto mol = SmilesParse::SmilesToMol("CCO[H]");
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 3);
  }
  {
    SmilesParse::SmilesParserParams ps;
    ps.removeHs = false;
    auto mol = SmilesParse::SmilesToMol("CCO[H]", ps);
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 4);
  }
  {
    auto mol = SmilesParse::SmartsToMol("[H]CC[R]");
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 4);
  }
  {
    SmilesParse::SmartsParserParams ps;
    ps.mergeHs = true;
    auto mol = SmilesParse::SmartsToMol("[H]CC[R]", ps);
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 3);
  }
  {
    auto atm = SmilesParse::SmilesToAtom("C");
    REQUIRE(atm);
  }
  {
    auto bnd = SmilesParse::SmilesToBond("-");
    REQUIRE(bnd);
  }
  {
    auto atm = SmilesParse::SmartsToAtom("[R]");
    REQUIRE(atm);
  }
  {
    auto bnd = SmilesParse::SmartsToBond("@");
    REQUIRE(bnd);
  }
}