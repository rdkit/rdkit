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
#include <GraphMol/FileParsers/SequenceParsers.h>

using namespace RDKit;

TEST_CASE("parsing") {
  SECTION("abbreviations") {
    auto abbrevs = Abbreviations::Utils::getDefaultAbbreviations();
    CHECK(abbrevs.size() == 29);
    CHECK(abbrevs[0].label == "CO2Et");
    CHECK(abbrevs[0].displayLabel == "CO<sub>2</sub>Et");
    CHECK(abbrevs[0].displayLabelW == "EtO<sub>2</sub>C");
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
    CHECK(abbrevs[0].label == "PEG4");
    CHECK(abbrevs[0].displayLabel == "PEG4");
    CHECK(abbrevs[0].displayLabelW.empty());
    CHECK(abbrevs[0].smarts == "*OCCOCCOCCOCC*");
    REQUIRE(abbrevs[0].mol);
    CHECK(abbrevs[0].mol->getNumAtoms() == 13);
    unsigned int nDummies = 0;
    CHECK(abbrevs[0].mol->getPropIfPresent(
        Abbreviations::common_properties::numDummies, nDummies));
    CHECK(nDummies == 1);
  }
  SECTION("bad SMILES in defintions") {
    const std::string defns = R"ABBREVS(CO2Et    C(=O)OCC
COOEt    fail
OiBu     OCC(C)C)ABBREVS";
    auto abbrevs = Abbreviations::Utils::parseAbbreviations(defns);
    REQUIRE(abbrevs.size() == 2);
    CHECK(abbrevs[0].label == "CO2Et");
    CHECK(abbrevs[1].label == "OiBu");
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
      CHECK(matches[0].abbrev.label == "CF3");
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
      CHECK(matches[0].abbrev.label == "CF3");
      CHECK(matches[1].abbrev.label == "CO2H");
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
      CHECK(matches[0].abbrev.label == "Et");
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
      CHECK(matches[0].abbrev.label == "PEG3");
      CHECK(matches[1].abbrev.label == "pentyl");
    }
    {  // directly connected
      auto m = "FCOCCOCCOCCCCCCCCl"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      auto matches = Abbreviations::findApplicableAbbreviationMatches(
          *m, linkers, maxCoverage);
      CHECK(matches.size() == 2);
      CHECK(matches[0].abbrev.label == "PEG3");
      CHECK(matches[1].abbrev.label == "pentyl");
      CHECK(matches[0].match[9].second == 10);
      CHECK(matches[1].match[0].second == 10);
    }
  }
}

TEST_CASE("applyMatches") {
  auto abbrevs = Abbreviations::Utils::getDefaultAbbreviations();
  SECTION("basics") {
    {
      auto m = "FC(F)(F)CC(=O)O"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      auto matches = Abbreviations::findApplicableAbbreviationMatches(
          *m, abbrevs, maxCoverage);
      CHECK(matches.size() == 2);
      Abbreviations::applyMatches(*m, matches);
      CHECK(m->getNumAtoms() == 3);
      CHECK(MolToCXSmiles(*m) == "*C* |$CF3;;CO2H$|");
    }
  }
}

TEST_CASE("applyMatches linkers") {
  auto linkers = Abbreviations::Utils::getDefaultLinkers();
  SECTION("basics") {
    {
      auto m = "FCOCCOCCOCCCCCCCCl"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      auto matches = Abbreviations::findApplicableAbbreviationMatches(
          *m, linkers, maxCoverage);
      CHECK(matches.size() == 2);
      Abbreviations::applyMatches(*m, matches);
      CHECK(m->getNumAtoms() == 5);
      CHECK(MolToCXSmiles(*m) == "FC**Cl |$;;PEG3;pentyl;$|");
    }
    {
      auto m = "COC1CCC(C)CC1"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      auto matches = Abbreviations::findApplicableAbbreviationMatches(
          *m, linkers, maxCoverage);
      CHECK(matches.size() == 1);
      Abbreviations::applyMatches(*m, matches);
      CHECK(m->getNumAtoms() == 4);
      CHECK(MolToCXSmiles(*m) == "C*OC |$;cyhex;;$|");
    }
  }
}

TEST_CASE("condense abbreviations") {
  auto abbrevs = Abbreviations::Utils::getDefaultAbbreviations();
  SECTION("basics") {
    {
      auto m = "FC(F)(F)CC(=O)O"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      Abbreviations::condenseMolAbbreviations(*m, abbrevs, maxCoverage);
      CHECK(MolToCXSmiles(*m) == "*C* |$CF3;;CO2H$|");
    }
  }
}

TEST_CASE("condense abbreviations linkers") {
  auto linkers = Abbreviations::Utils::getDefaultLinkers();
  SECTION("basics") {
    {
      auto m = "FCOCCOCCOCCCCCCCCl"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      Abbreviations::condenseMolAbbreviations(*m, linkers, maxCoverage);
      CHECK(m->getNumAtoms() == 5);
      CHECK(MolToCXSmiles(*m) == "FC**Cl |$;;PEG3;pentyl;$|");
    }
    {
      auto m = "COC1CCC(C)CC1"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      Abbreviations::condenseMolAbbreviations(*m, linkers, maxCoverage);
      CHECK(m->getNumAtoms() == 4);
      CHECK(MolToCXSmiles(*m) == "C*OC |$;cyhex;;$|");
    }
  }
  SECTION("peptides") {
    std::unique_ptr<RWMol> m(SequenceToMol("GYTKC"));
    REQUIRE(m);
    double maxCoverage = 1.0;
    Abbreviations::condenseMolAbbreviations(*m, linkers, maxCoverage);
    CHECK(MolToCXSmiles(*m) == "NCC(=O)****O |$;;;;tyr;thr;lys;cys;$|");
  }
}

TEST_CASE("abbreviations and linkers") {
  auto abbrevs = Abbreviations::Utils::getDefaultAbbreviations();
  auto linkers = Abbreviations::Utils::getDefaultLinkers();
  SECTION("basics") {
    {  // this isn't the order we'd normally do this in:
      auto m = "COC1CCC(C)CC1"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      Abbreviations::condenseMolAbbreviations(*m, abbrevs, maxCoverage);
      CHECK(m->getNumAtoms() == 8);
      CHECK(MolToCXSmiles(*m) == "*C1CCC(C)CC1 |$OMe;;;;;;;$|");
      Abbreviations::condenseMolAbbreviations(*m, linkers, maxCoverage);
      CHECK(m->getNumAtoms() == 3);
      CHECK(MolToCXSmiles(*m) == "**C |$OMe;cyhex;$|");
    }
    {  // a more sensible order
      auto m = "COC1CCC(C)CC1"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      Abbreviations::condenseMolAbbreviations(*m, linkers, maxCoverage);
      CHECK(m->getNumAtoms() == 4);
      CHECK(MolToCXSmiles(*m) == "C*OC |$;cyhex;;$|");
      Abbreviations::condenseMolAbbreviations(*m, abbrevs, maxCoverage);
      CHECK(m->getNumAtoms() == 4);
      CHECK(MolToCXSmiles(*m) == "C*OC |$;cyhex;;$|");
    }
  }
}

TEST_CASE("labelMatches") {
  auto abbrevs = Abbreviations::Utils::getDefaultAbbreviations();
  SECTION("basics") {
    {
      auto m = "CC(C)CC(F)(F)F"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      auto matches = Abbreviations::findApplicableAbbreviationMatches(
          *m, abbrevs, maxCoverage);
      CHECK(matches.size() == 2);
      Abbreviations::labelMatches(*m, matches);
      CHECK(m->getNumAtoms() == 8);
      const auto &sgs = getSubstanceGroups(*m);
      REQUIRE(sgs.size() == 2);
      CHECK(sgs[0].getProp<std::string>("TYPE") == "SUP");
      CHECK(sgs[0].getProp<std::string>("LABEL") == "iPr");
      CHECK(sgs[0].getBonds() == std::vector<unsigned int>({2}));
      CHECK(sgs[0].getAtoms() == std::vector<unsigned int>({1, 0, 2}));
      CHECK(sgs[0].getAttachPoints().size() == 1);
      CHECK(sgs[0].getAttachPoints()[0].aIdx == 1);
      CHECK(sgs[0].getAttachPoints()[0].lvIdx == 3);

      CHECK(sgs[1].getProp<std::string>("TYPE") == "SUP");
      CHECK(sgs[1].getProp<std::string>("LABEL") == "CF3");
      CHECK(sgs[1].getBonds() == std::vector<unsigned int>({3}));
      CHECK(sgs[1].getAtoms() == std::vector<unsigned int>({4, 5, 6, 7}));
      CHECK(sgs[1].getAttachPoints().size() == 1);
      CHECK(sgs[1].getAttachPoints()[0].aIdx == 4);
      CHECK(sgs[1].getAttachPoints()[0].lvIdx == 3);
    }
  }
}

TEST_CASE("labelMolAbbreviations") {
  auto abbrevs = Abbreviations::Utils::getDefaultAbbreviations();
  SECTION("basics") {
    {
      auto m = "CC(C)CC(F)(F)F"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      Abbreviations::labelMolAbbreviations(*m, abbrevs, maxCoverage);
      CHECK(m->getNumAtoms() == 8);
      const auto &sgs = getSubstanceGroups(*m);
      REQUIRE(sgs.size() == 2);
      CHECK(sgs[0].getProp<std::string>("TYPE") == "SUP");
      CHECK(sgs[0].getProp<std::string>("LABEL") == "iPr");
      CHECK(sgs[0].getBonds() == std::vector<unsigned int>({2}));
      CHECK(sgs[0].getAtoms() == std::vector<unsigned int>({1, 0, 2}));
      CHECK(sgs[0].getAttachPoints().size() == 1);
      CHECK(sgs[0].getAttachPoints()[0].aIdx == 1);
      CHECK(sgs[0].getAttachPoints()[0].lvIdx == 3);

      CHECK(sgs[1].getProp<std::string>("TYPE") == "SUP");
      CHECK(sgs[1].getProp<std::string>("LABEL") == "CF3");
      CHECK(sgs[1].getBonds() == std::vector<unsigned int>({3}));
      CHECK(sgs[1].getAtoms() == std::vector<unsigned int>({4, 5, 6, 7}));
      CHECK(sgs[1].getAttachPoints().size() == 1);
      CHECK(sgs[1].getAttachPoints()[0].aIdx == 4);
      CHECK(sgs[1].getAttachPoints()[0].lvIdx == 3);
    }
  }
}
