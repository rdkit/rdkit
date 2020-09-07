//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"
#include "RDGeneral/test.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Abbreviations/Abbreviations.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

TEST_CASE("parsing", "[abbreviations]") {
  SECTION("abbreviations") {
    auto abbrevs = Abbreviations::Utils::getDefaultAbbreviations();
    CHECK(abbrevs.size() == 29);
    CHECK(abbrevs[0].llabel == "CO2Et");
    CHECK(abbrevs[0].rlabel == "EtO2C");
    CHECK(abbrevs[0].smarts == "C(=O)OCC");
    REQUIRE(abbrevs[0].mol);
    CHECK(abbrevs[0].mol->getNumAtoms() == 6);
    unsigned int nDummies = 0;
    CHECK(abbrevs[0].mol->getPropIfPresent(
        Abbreviations::common_properties::numDummies, nDummies));
    CHECK(nDummies == 1);
  }
  SECTION("linkers") {
    auto abbrevs = Abbreviations::Utils::getDefaultLinkers();
    CHECK(abbrevs.size() == 24);
    CHECK(abbrevs[0].llabel == "PEG4");
    CHECK(abbrevs[0].rlabel == "PEG4");
    CHECK(abbrevs[0].smarts == "*OCCOCCOCCOCC*");
    REQUIRE(abbrevs[0].mol);
    CHECK(abbrevs[0].mol->getNumAtoms() == 13);
    unsigned int nDummies = 0;
    CHECK(abbrevs[0].mol->getPropIfPresent(
        Abbreviations::common_properties::numDummies, nDummies));
    CHECK(nDummies == 1);
  }
}