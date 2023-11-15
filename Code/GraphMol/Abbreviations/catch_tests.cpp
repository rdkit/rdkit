//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>
#include "RDGeneral/test.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Abbreviations/Abbreviations.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/SequenceParsers.h>
#include <GraphMol/FileParsers/FileParsers.h>

using namespace RDKit;

TEST_CASE("parsing") {
  SECTION("abbreviations") {
    auto abbrevs = Abbreviations::Utils::getDefaultAbbreviations();
    CHECK(abbrevs.size() == 37);
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
    CHECK(abbrevs.size() == 8);
    CHECK(abbrevs[0].label == "PEG6");
    CHECK(abbrevs[0].displayLabel == "PEG6");
    CHECK(abbrevs[0].displayLabelW.empty());
    CHECK(abbrevs[0].smarts == "*OCCOCCOCCOCCOCCOCC*");
    REQUIRE(abbrevs[0].mol);
    CHECK(abbrevs[0].mol->getNumAtoms() == 19);
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
      auto m = "FCOCCOCCOCCNCCCCCCCCl"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      auto matches = Abbreviations::findApplicableAbbreviationMatches(
          *m, linkers, maxCoverage);
      CHECK(matches.size() == 2);
      CHECK(matches[0].abbrev.label == "PEG3");
      CHECK(matches[1].abbrev.label == "Hept");
    }
    {  // directly connected
      auto m = "FCOCCOCCOCCCCCCCCCCl"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      auto matches = Abbreviations::findApplicableAbbreviationMatches(
          *m, linkers, maxCoverage);
      CHECK(matches.size() == 2);
      CHECK(matches[0].abbrev.label == "PEG3");
      CHECK(matches[1].abbrev.label == "Hept");
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
      std::vector<unsigned int> atomMapping;
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origAtomMapping, atomMapping));
      CHECK(atomMapping == std::vector<unsigned int>{1, 4, 5});
      std::vector<unsigned int> bondMapping;
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origBondMapping, bondMapping));
      CHECK(bondMapping == std::vector<unsigned int>{3, 4});
    }
  }
}

TEST_CASE("applyMatches linkers") {
  auto linkers =
      Abbreviations::Utils::parseLinkers(R"ABBREV(PEG3  *OCCOCCOCC* PEG3
Pent  *CCCCC*
Cy   *C1CCC(*)CC1  Cy)ABBREV");
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
      CHECK(MolToCXSmiles(*m) == "FC**Cl |$;;PEG3;Pent;$|");
      std::vector<unsigned int> atomMapping;
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origAtomMapping, atomMapping));
      CHECK(atomMapping == std::vector<unsigned int>{0, 1, 2, 11, 16});
      std::vector<unsigned int> bondMapping;
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origBondMapping, bondMapping));
      CHECK(bondMapping == std::vector<unsigned int>{0, 1, 10, 15});
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
      CHECK(MolToCXSmiles(*m) == "C*OC |$;Cy;;$|");
      std::vector<unsigned int> atomMapping;
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origAtomMapping, atomMapping));
      CHECK(atomMapping == std::vector<unsigned int>{0, 1, 2, 6});
      std::vector<unsigned int> bondMapping;
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origBondMapping, bondMapping));
      CHECK(bondMapping == std::vector<unsigned int>{0, 1, 5});
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
      std::vector<unsigned int> atomMapping;
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origAtomMapping, atomMapping));
      CHECK(atomMapping == std::vector<unsigned int>{1, 4, 5});
      std::vector<unsigned int> bondMapping;
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origBondMapping, bondMapping));
      CHECK(bondMapping == std::vector<unsigned int>{3, 4});
    }
  }
}

TEST_CASE("condense abbreviations linkers") {
  auto linkers = Abbreviations::Utils::getDefaultLinkers();
  auto customLinkers =
      Abbreviations::Utils::parseLinkers(R"ABBREV(PEG3  *OCCOCCOCC* PEG3
Pent  *CCCCC*
Cy   *C1CCC(*)CC1  Cy
ala *N[C@@H](C)C(=O)* ala
arg *N[C@@H](CCCNC(N)=[NH])C(=O)* arg
asn *N[C@@H](CC(N)=O)C(=O)* asn
asp *N[C@@H](CC(O)=O)C(=O)* asp
cys *N[C@@H](CS)C(=O)* cys
gln *N[C@@H](CCC(N)=O)C(=O)* gln
glu *N[C@@H](CCC(O)=O)C(=O)* glu
gly *NCC(=O)* gly
his *N[C@@H](Cc1c[nH]cn1)C(=O)* his
ile *N[C@@H](C(C)CC)C(=O)* ile
leu *N[C@@H](CC(C)C)C(=O)* leu
lys *N[C@@H](CCCCN)C(=O)* lys
met *N[C@@H](CCSC)C(=O)* met
phe *N[C@@H](Cc1ccccc1)C(=O)* phe
pro *N1[C@@H](CCC1)C(=O)* pro
ser *N[C@@H](CO)C(=O)* ser
thr *N[C@@H](C(O)C)C(=O)* thr
trp *N[C@@H](Cc1c[nH]c2ccccc21)C(=O)* trp
tyr *N[C@@H](Cc1ccc(O)cc1)C(=O)* tyr
val *N[C@@H](C(C)C)C(=O)* val)ABBREV");
  SECTION("basics") {
    {
      auto m = "FCOCCOCCOCCCCCCCCCCl"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      Abbreviations::condenseMolAbbreviations(*m, linkers, maxCoverage);
      CHECK(m->getNumAtoms() == 5);
      CHECK(MolToCXSmiles(*m) == "FC**Cl |$;;PEG3;Hept;$|");
      std::vector<unsigned int> atomMapping;
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origAtomMapping, atomMapping));
      CHECK(atomMapping == std::vector<unsigned int>{0, 1, 2, 11, 18});
      std::vector<unsigned int> bondMapping;
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origBondMapping, bondMapping));
      CHECK(bondMapping == std::vector<unsigned int>{0, 1, 10, 17});
    }
    {
      auto m = "COC1CCC(C)CC1"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      Abbreviations::condenseMolAbbreviations(*m, customLinkers, maxCoverage);
      CHECK(m->getNumAtoms() == 4);
      CHECK(MolToCXSmiles(*m) == "C*OC |$;Cy;;$|");
      std::vector<unsigned int> atomMapping;
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origAtomMapping, atomMapping));
      CHECK(atomMapping == std::vector<unsigned int>{0, 1, 2, 6});
      std::vector<unsigned int> bondMapping;
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origBondMapping, bondMapping));
      CHECK(bondMapping == std::vector<unsigned int>{0, 1, 5});
    }
  }
  SECTION("peptides") {
    std::unique_ptr<RWMol> m(SequenceToMol("GYTKC"));
    REQUIRE(m);
    double maxCoverage = 1.0;
    Abbreviations::condenseMolAbbreviations(*m, customLinkers, maxCoverage);
    CHECK(MolToCXSmiles(*m) == "NCC(=O)****O |$;;;;tyr;thr;lys;cys;$|");
    std::vector<unsigned int> atomMapping;
    CHECK(m->getPropIfPresent(Abbreviations::common_properties::origAtomMapping,
                              atomMapping));
    CHECK(atomMapping ==
          std::vector<unsigned int>{0, 1, 2, 3, 4, 16, 23, 32, 38});
    std::vector<unsigned int> bondMapping;
    CHECK(m->getPropIfPresent(Abbreviations::common_properties::origBondMapping,
                              bondMapping));
    CHECK(bondMapping ==
          std::vector<unsigned int>{0, 1, 2, 15, 38, 37, 31, 22});
  }
}

TEST_CASE("abbreviations and linkers") {
  auto abbrevs = Abbreviations::Utils::getDefaultAbbreviations();
  auto linkers = Abbreviations::Utils::parseLinkers(
      R"ABBREV(Cy   *C1CCC(*)CC1  Cy)ABBREV");
  SECTION("basics") {
    {  // this isn't the order we'd normally do this in:
      auto m = "COC1CCC(C)CC1"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      Abbreviations::condenseMolAbbreviations(*m, abbrevs, maxCoverage);
      CHECK(m->getNumAtoms() == 8);
      CHECK(MolToCXSmiles(*m) == "*C1CCC(C)CC1 |$OMe;;;;;;;$|");
      std::vector<unsigned int> atomMapping;
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origAtomMapping, atomMapping));
      CHECK(atomMapping == std::vector<unsigned int>{1, 2, 3, 4, 5, 6, 7, 8});
      std::vector<unsigned int> bondMapping;
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origBondMapping, bondMapping));
      CHECK(bondMapping == std::vector<unsigned int>{1, 2, 3, 4, 5, 6, 7, 8});
      Abbreviations::condenseMolAbbreviations(*m, linkers, maxCoverage);
      CHECK(m->getNumAtoms() == 3);
      CHECK(MolToCXSmiles(*m) == "**C |$OMe;Cy;$|");
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origAtomMapping, atomMapping));
      CHECK(atomMapping == std::vector<unsigned int>{1, 2, 6});
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origBondMapping, bondMapping));
      CHECK(bondMapping == std::vector<unsigned int>{1, 5});
    }
    {  // a more sensible order
      auto m = "COC1CCC(C)CC1"_smiles;
      REQUIRE(m);
      double maxCoverage = 1.0;
      Abbreviations::condenseMolAbbreviations(*m, linkers, maxCoverage);
      CHECK(m->getNumAtoms() == 4);
      CHECK(MolToCXSmiles(*m) == "C*OC |$;Cy;;$|");
      std::vector<unsigned int> atomMapping;
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origAtomMapping, atomMapping));
      CHECK(atomMapping == std::vector<unsigned int>{0, 1, 2, 6});
      std::vector<unsigned int> bondMapping;
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origBondMapping, bondMapping));
      CHECK(bondMapping == std::vector<unsigned int>{0, 1, 5});
      Abbreviations::condenseMolAbbreviations(*m, abbrevs, maxCoverage);
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origAtomMapping, atomMapping));
      CHECK(m->getNumAtoms() == 4);
      CHECK(MolToCXSmiles(*m) == "C*OC |$;Cy;;$|");
      CHECK(atomMapping == std::vector<unsigned int>{0, 1, 2, 6});
      CHECK(m->getPropIfPresent(
          Abbreviations::common_properties::origBondMapping, bondMapping));
      CHECK(bondMapping == std::vector<unsigned int>{0, 1, 5});
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

TEST_CASE("condenseAbbreviationSubstanceGroups") {
  SECTION("abbreviations") {
    auto m = R"CTAB(
  ACCLDraw09152005292D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 10 10 2 0 0
M  V30 BEGIN ATOM
M  V30 1 C 12.8333 -9.32 0 0 CFG=3 
M  V30 2 C 13.8565 -8.7293 0 0 
M  V30 3 O 14.8802 -9.3201 0 0 
M  V30 4 O 13.8565 -7.5471 0 0 
M  V30 5 C 11.6489 -9.32 0 0 
M  V30 6 C 12.241 -10.3432 0 0 CFG=3 
M  V30 7 C 12.241 -11.5253 0 0 CFG=3 
M  V30 8 F 12.241 -12.5874 0 0 
M  V30 9 F 11.0366 -11.5253 0 0 
M  V30 10 F 13.4231 -11.5253 0 0 
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 2 4 
M  V30 2 1 2 3 
M  V30 3 1 1 2 
M  V30 4 1 5 6 
M  V30 5 1 5 1 
M  V30 6 1 1 6 
M  V30 7 1 7 10 
M  V30 8 1 7 9 
M  V30 9 1 7 8 
M  V30 10 1 6 7 
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 1 ATOMS=(3 2 3 4) XBONDS=(1 3) CSTATE=(4 3 -1.02 -0.59 0) LABEL=-
M  V30 CO2H 
M  V30 2 SUP 2 ATOMS=(4 7 8 9 10) XBONDS=(1 10) CSTATE=(4 10 0 1.18 0) LABEL=-
M  V30 CF3 
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 10);
    Abbreviations::condenseAbbreviationSubstanceGroups(*m);
    CHECK(m->getNumAtoms() == 5);
    std::vector<unsigned int> atomMapping;
    CHECK(m->getPropIfPresent(Abbreviations::common_properties::origAtomMapping,
                              atomMapping));
    CHECK(atomMapping == std::vector<unsigned int>{0, 1, 4, 5, 6});
    std::vector<unsigned int> bondMapping;
    CHECK(m->getPropIfPresent(Abbreviations::common_properties::origBondMapping,
                              bondMapping));
    CHECK(bondMapping == std::vector<unsigned int>{2, 3, 4, 5, 9});
    // remove the conformer before generating CXSMILES
    m->clearConformers();
    CHECK(MolToCXSmiles(*m) == "*C1CC1* |$CO2H;;;;CF3$|");
  }
  SECTION("abbreviations MRV") {
    auto m = R"CTAB(
  Mrv2014 09152006492D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 7 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 5.25 -5.9858 0 0
M  V30 2 C 4.48 -7.3196 0 0
M  V30 3 C 6.02 -7.3196 0 0
M  V30 4 F 8.6873 -8.8596 0 0
M  V30 5 C 7.3537 -8.0896 0 0
M  V30 6 F 6.02 -8.8596 0 0
M  V30 7 F 7.3537 -6.5496 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 1
M  V30 3 1 2 3
M  V30 4 1 3 5
M  V30 5 1 4 5
M  V30 6 1 5 6
M  V30 7 1 5 7
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 ATOMS=(4 4 5 6 7) SAP=(3 5 3 1) XBONDS=(1 4) LABEL=CF3
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 7);
    Abbreviations::condenseAbbreviationSubstanceGroups(*m);
    CHECK(m->getNumAtoms() == 4);
    // remove the conformer before generating CXSMILES
    Abbreviations::condenseAbbreviationSubstanceGroups(*m);
    std::vector<unsigned int> atomMapping;
    CHECK(m->getPropIfPresent(Abbreviations::common_properties::origAtomMapping,
                              atomMapping));
    CHECK(atomMapping == std::vector<unsigned int>{0, 1, 2, 4});
    std::vector<unsigned int> bondMapping;
    CHECK(m->getPropIfPresent(Abbreviations::common_properties::origBondMapping,
                              bondMapping));
    CHECK(bondMapping == std::vector<unsigned int>{0, 1, 2, 3});
    m->clearConformers();
    CHECK(MolToCXSmiles(*m) == "*C1CC1 |$CF3;;;$|");
  }

  SECTION("linker") {
    auto m = R"CTAB(
  ACCLDraw09152006102D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 7.2482 -5.1911 0 0 
M  V30 2 O 5.8143 -6.2327 0 0 
M  V30 3 C 6.77 -5.5382 0 0 
M  V30 4 C 7.8494 -6.0186 0 0 
M  V30 5 O 8.8052 -5.3241 0 0 
M  V30 6 C 9.8845 -5.8046 0 0 
M  V30 7 C 10.8403 -5.1101 0 0 
M  V30 8 C 9.4066 -6.1518 0 0 
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 
M  V30 2 1 2 3 
M  V30 3 1 3 4 
M  V30 4 1 4 5 
M  V30 5 1 5 6 
M  V30 6 1 6 7 
M  V30 7 1 7 8 
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 1 ATOMS=(6 2 3 4 5 6 7) XBONDS=(2 1 7) CSTATE=(4 1 -1.08 0.48 0) -
M  V30 CSTATE=(4 7 1.08 -0.48 0) LABEL=PEG2 
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 8);
    Abbreviations::condenseAbbreviationSubstanceGroups(*m);
    std::vector<unsigned int> atomMapping;
    CHECK(m->getPropIfPresent(Abbreviations::common_properties::origAtomMapping,
                              atomMapping));
    CHECK(atomMapping == std::vector<unsigned int>{0, 1, 7});
    std::vector<unsigned int> bondMapping;
    CHECK(m->getPropIfPresent(Abbreviations::common_properties::origBondMapping,
                              bondMapping));
    CHECK(bondMapping == std::vector<unsigned int>{0, 6});
    CHECK(m->getNumAtoms() == 3);
    // remove the conformer before generating CXSMILES
    m->clearConformers();
    CHECK(MolToCXSmiles(*m) == "C*C |$;PEG2;$|");
  }
  SECTION("linker MRV") {
    auto m = R"CTAB(
  Mrv2014 09152006522D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1.625 -8.9167 0 0
M  V30 2 O 2.9587 -8.1467 0 0
M  V30 3 C 4.2924 -8.9167 0 0
M  V30 4 C 5.626 -8.1467 0 0
M  V30 5 O 6.9597 -8.9167 0 0
M  V30 6 C 8.2934 -8.1467 0 0
M  V30 7 C 9.6271 -8.9167 0 0
M  V30 8 C 10.9608 -8.1467 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 7 1 7 8
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 ATOMS=(6 2 3 4 5 6 7) SAP=(3 2 1 1) SAP=(3 7 8 2) XBONDS=(2 1 -
M  V30 7) LABEL=PEG2 ESTATE=E
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 8);
    Abbreviations::condenseAbbreviationSubstanceGroups(*m);
    std::vector<unsigned int> atomMapping;
    CHECK(m->getPropIfPresent(Abbreviations::common_properties::origAtomMapping,
                              atomMapping));
    CHECK(atomMapping == std::vector<unsigned int>{0, 1, 7});
    std::vector<unsigned int> bondMapping;
    CHECK(m->getPropIfPresent(Abbreviations::common_properties::origBondMapping,
                              bondMapping));
    CHECK(bondMapping == std::vector<unsigned int>{0, 6});
    CHECK(m->getNumAtoms() == 3);
    // remove the conformer before generating CXSMILES
    m->clearConformers();
    CHECK(MolToCXSmiles(*m) == "C*C |$;PEG2;$|");
  }
}

TEST_CASE("comparison") {
  auto abbrevs = Abbreviations::Utils::getDefaultAbbreviations();
  Abbreviations::AbbreviationDefinition cp = abbrevs[0];
  CHECK(cp == abbrevs[0]);
  CHECK(cp != abbrevs[1]);
  CHECK(abbrevs[1] == abbrevs[1]);
}