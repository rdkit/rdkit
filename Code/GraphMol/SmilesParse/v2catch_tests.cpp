//
//  Copyright (C) 2023 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"
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

TEST_CASE("v2 basics") {
  {
    auto mol = RDKitv2::SmilesParse::SmilesToMol("CCC");
    REQUIRE(mol);
  }
  {
    auto mol = RDKitv2::SmilesParse::SmartsToMol("CC[R]");
    REQUIRE(mol);
  }
  {
    auto atm = RDKitv2::SmilesParse::SmilesToAtom("C");
    REQUIRE(atm);
  }
  {
    auto bnd = RDKitv2::SmilesParse::SmilesToBond("-");
    REQUIRE(bnd);
  }
  {
    auto atm = RDKitv2::SmilesParse::SmartsToAtom("[R]");
    REQUIRE(atm);
  }
  {
    auto bnd = RDKitv2::SmilesParse::SmartsToBond("@");
    REQUIRE(bnd);
  }
}