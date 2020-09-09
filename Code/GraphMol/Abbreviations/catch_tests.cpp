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

TEST_CASE("parsing") {
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
  SECTION("bad SMILES in defintions") {
    const std::string defns = R"ABBREVS(CO2Et    EtO2C    C(=O)OCC
COOEt    EtOOC    fail
OiBu     iBuO     OCC(C)C)ABBREVS";
    auto abbrevs = Abbreviations::Utils::parseAbbreviations(defns);
    REQUIRE(abbrevs.size() == 2);
    CHECK(abbrevs[0].llabel == "CO2Et");
    CHECK(abbrevs[1].llabel == "OiBu");
  }
}

TEST_CASE("findApplicableMatches") {
  auto abbrevs = Abbreviations::Utils::getDefaultAbbreviations();
  SECTION("basics") {
    auto m = "NCCC(F)(F)F"_smiles;
    REQUIRE(m);
    {
      double maxCoverage = 0.4;
      auto matches = Abbreviations::findApplicableAbbreviationMatches(
          *m, abbrevs, maxCoverage);
      CHECK(matches.empty());
    }
    {
      double maxCoverage = 1.0;
      auto matches = Abbreviations::findApplicableAbbreviationMatches(
          *m, abbrevs, maxCoverage);
      CHECK(matches.size() == 1);
      CHECK(matches[0].abbrev.llabel == "CF3");
      CHECK(matches[0].match[0].second == 2);
      CHECK(matches[0].match[1].second == 3);
    }
  }
  SECTION("multiple abbreviations") {
    {
      auto m = "FC(F)(F)CC(=O)O"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      auto matches = Abbreviations::findApplicableAbbreviationMatches(
          *m, abbrevs, maxCoverage);
      CHECK(matches.size() == 2);
      CHECK(matches[0].abbrev.llabel == "CF3");
      CHECK(matches[1].abbrev.llabel == "CO2H");
    }
    {  // overlapping
      auto m = "FC(F)(F)C(=O)O"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      auto matches = Abbreviations::findApplicableAbbreviationMatches(
          *m, abbrevs, maxCoverage);
      CHECK(matches.empty());
    }
    {  // overlapping
      auto m = "FC(F)(F)C(F)(F)F"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      auto matches = Abbreviations::findApplicableAbbreviationMatches(
          *m, abbrevs, maxCoverage);
      CHECK(matches.empty());
    }
    {  // overlapping, one is too big, so there is an abbreviation for the other
      auto m = "CCC(F)(F)F"_smiles;
      REQUIRE(m);
      double maxCoverage = 0.4;
      auto matches = Abbreviations::findApplicableAbbreviationMatches(
          *m, abbrevs, maxCoverage);
      CHECK(matches.size() == 1);
      CHECK(matches[0].abbrev.llabel == "Et");
      // remove the size constraint and there's no abbreviation:
      maxCoverage = 1.0;
      matches = Abbreviations::findApplicableAbbreviationMatches(*m, abbrevs,
                                                                 maxCoverage);
      CHECK(matches.empty());
    }
  }
}
TEST_CASE("findApplicableMatches linkers") {
  auto linkers = Abbreviations::Utils::getDefaultLinkers();
  SECTION("basics") {
    {
      auto m = "FCOCCOCCOCCNCCCCCCl"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      auto matches = Abbreviations::findApplicableAbbreviationMatches(
          *m, linkers, maxCoverage);
      CHECK(matches.size() == 2);
      CHECK(matches[0].abbrev.llabel == "PEG3");
      CHECK(matches[1].abbrev.llabel == "pentyl");
    }
    {  // directly connected
      auto m = "FCOCCOCCOCCCCCCCCl"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      auto matches = Abbreviations::findApplicableAbbreviationMatches(
          *m, linkers, maxCoverage);
      CHECK(matches.size() == 2);
      CHECK(matches[0].abbrev.llabel == "PEG3");
      CHECK(matches[1].abbrev.llabel == "pentyl");
      CHECK(matches[0].match[9].second == 10);
      CHECK(matches[1].match[0].second == 10);
    }
  }
}