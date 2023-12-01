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

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Chirality.h>
#include <algorithm>

using namespace RDKit;

#ifdef RDK_USE_BOOST_SERIALIZATION

TEST_CASE("MolBundle serialization") {
  SECTION("basics") {
    MolBundle bundle;
    bundle.addMol(ROMOL_SPTR(SmilesToMol("CCC")));
    bundle.addMol(ROMOL_SPTR(SmilesToMol("CCN")));
    CHECK(!bundle.empty());
    auto pkl = bundle.serialize();
    MolBundle nbundle(pkl);
    REQUIRE(bundle.size() == nbundle.size());
    for (auto i = 0u; i < bundle.size(); ++i) {
      CHECK(MolToSmiles(*bundle[i]) == MolToSmiles(*nbundle[i]));
    }
  }
  SECTION("empty") {
    MolBundle bundle;
    CHECK(bundle.empty());
    auto pkl = bundle.serialize();
    MolBundle nbundle(pkl);
    REQUIRE(bundle.size() == nbundle.size());
  }
}
#else
TEST_CASE("MolBundle serialization") {}
#endif