//
//  Copyright (c) 2022 Brian P Kelley
//  All rights reserved.
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <RDGeneral/FileParseException.h>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/SmilesParse/CanonicalizeStereoGroups.h>
#include <filesystem>

using namespace RDKit;
using namespace RDKit::v2::CDXMLParser;

std::string canon(const std::string &smi) {
  auto *m = SmilesToMol(smi);
  auto res = MolToSmiles(*m);
  delete m;
  return res;
}

void check_smiles_and_roundtrip(const RWMol &m, const std::string &expected) {
  CHECK(MolToSmiles(m) == expected);
  //  std::cout << "*********" << std::endl;
  //  std::cout << MolToMolBlock(m) << std::endl;
  std::unique_ptr<RWMol> mol(MolBlockToMol(MolToMolBlock(m)));
  CHECK(MolToSmiles(*mol) == expected);
}

TEST_CASE("CDXML") {
  std::string cdxmlbase =
      std::string(getenv("RDBASE")) + "/Code/GraphMol/test_data/CDXML/";
  SECTION("SIMPLE") {
    std::string cdxml1 = R"(<?xml version="1.0" encoding="UTF-8" ?>
        <!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >
        <CDXML
         CreationProgram="ChemDraw JS 2.0.0.9"
         Name="ACS Document 1996"
         BoundingBox="94.75 178.16 154.89 211.17"
         WindowPosition="0 0"
         WindowSize="0 0"
         FractionalWidths="yes"
         InterpretChemically="yes"
         ShowAtomQuery="yes"
         ShowAtomStereo="no"
         ShowAtomEnhancedStereo="yes"
         ShowAtomNumber="no"
         ShowResidueID="no"
         ShowBondQuery="yes"
         ShowBondRxn="yes"
         ShowBondStereo="no"
         ShowTerminalCarbonLabels="no"
         ShowNonTerminalCarbonLabels="no"
         HideImplicitHydrogens="no"
         Magnification="666"
         LabelFont="24"
         LabelSize="10"
         LabelFace="96"
         CaptionFont="24"
         CaptionSize="10"
         HashSpacing="2.50"
         MarginWidth="1.60"
         LineWidth="0.60"
         BoldWidth="2"
         BondLength="14.40"
         BondSpacing="18"
         ChainAngle="120"
         LabelJustification="Auto"
         CaptionJustification="Left"
         AminoAcidTermini="HOH"
         ShowSequenceTermini="yes"
         ShowSequenceBonds="yes"
         ShowSequenceUnlinkedBranches="no"
         ResidueWrapCount="40"
         ResidueBlockCount="10"
         ResidueZigZag="yes"
         NumberResidueBlocks="no"
         PrintMargins="36 36 36 36"
         MacPrintInfo="0003000001200120000000000B6608A0FF84FF880BE309180367052703FC0002000001200120000000000B6608A0000100000064000000010001010100000001270F000100010000000000000000000000000002001901900000000000400000000000000000000100000000000000000000000000000000"
         ChemPropName=""
         ChemPropFormula="Chemical Formula: "
         ChemPropExactMass="Exact Mass: "
         ChemPropMolWt="Molecular Weight: "
         ChemPropMOverZ="m/z: "
         ChemPropAnalysis="Elemental Analysis: "
         ChemPropBoilingPt="Boiling Point: "
         ChemPropMeltingPt="Melting Point: "
         ChemPropCritTemp="Critical Temp: "
         ChemPropCritPres="Critical Pres: "
         ChemPropCritVol="Critical Vol: "
         ChemPropGibbs="Gibbs Energy: "
         ChemPropLogP="Log P: "
         ChemPropMR="MR: "
         ChemPropHenry="Henry&apos;s Law: "
         ChemPropEForm="Heat of Form: "
         ChemProptPSA="tPSA: "
         ChemPropID=""
         ChemPropFragmentLabel=""
         color="0"
         bgcolor="1"
         RxnAutonumberStart="1"
         RxnAutonumberConditions="no"
         RxnAutonumberStyle="Roman"
        ><colortable>
        <color r="1" g="1" b="1"/>
        <color r="0" g="0" b="0"/>
        <color r="1" g="0" b="0"/>
        <color r="1" g="1" b="0"/>
        <color r="0" g="1" b="0"/>
        <color r="0" g="1" b="1"/>
        <color r="0" g="0" b="1"/>
        <color r="1" g="0" b="1"/>
        </colortable><fonttable>
        <font id="24" charset="utf-8" name="Arial"/>
        </fonttable><page
         id="32"
         BoundingBox="0 0 542 354"
         Width="542"
         Height="354"
         HeaderPosition="36"
         FooterPosition="36"
         PageOverlap="0"
         PrintTrimMarks="yes"
         HeightPages="1"
         WidthPages="1"
         DrawingSpace="poster"
        ><fragment
         id="10"
         BoundingBox="94.75 178.16 154.89 211.17"
         Z="4"
        ><n
         id="7"
         p="95.05 187.47"
         Z="1"
         AS="N"
        /><n
         id="9"
         p="95.05 201.87"
         Z="3"
         AS="N"
        /><n
         id="11"
         p="106.31 210.84"
         Z="5"
         AS="N"
        /><n
         id="13"
         p="120.35 207.64"
         Z="7"
         AS="N"
        /><n
         id="15"
         p="126.59 194.67"
         Z="9"
         AS="N"
        /><n
         id="17"
         p="120.35 181.69"
         Z="11"
         AS="N"
        /><n
         id="19"
         p="106.31 178.49"
         Z="13"
         AS="N"
        /><n
         id="28"
         p="140.99 194.67"
         Z="22"
         NodeType="Nickname"
         NeedsClean="yes"
         AS="N"
        ><fragment
         id="33"
        ><n
         id="34"
         p="148.17 207.09"
         Element="8"
         NumHydrogens="0"
        /><n
         id="35"
         p="162.52 207.09"
        /><n
         id="36"
         p="176.87 207.09"
        /><n
         id="37"
         p="169.69 194.67"
        /><n
         id="38"
         p="169.69 219.52"
        /><n
         id="39"
         p="140.99 194.67"
        /><n
         id="40"
         p="148.17 182.24"
         Element="8"
         NumHydrogens="0"
        /><n
         id="41"
         p="126.64 194.67"
         NodeType="ExternalConnectionPoint"
        /><b
         id="42"
         B="39"
         E="40"
         Order="2"
        /><b
         id="43"
         B="35"
         E="38"
        /><b
         id="44"
         B="35"
         E="36"
        /><b
         id="45"
         B="35"
         E="37"
        /><b
         id="46"
         B="34"
         E="35"
        /><b
         id="47"
         B="34"
         E="39"
        /><b
         id="48"
         B="41"
         E="39"
        /></fragment><t
         p="137.66 198.28"
         BoundingBox="137.66 189.64 154.89 198.28"
         LabelJustification="Left"
         LabelAlignment="Left"
        ><s font="24" size="9.95" color="0" face="96">Boc</s></t></n><b
         id="21"
         Z="15"
         B="7"
         E="9"
         BS="N"
        /><b
         id="22"
         Z="16"
         B="9"
         E="11"
         BS="N"
        /><b
         id="23"
         Z="17"
         B="11"
         E="13"
         BS="N"
        /><b
         id="24"
         Z="18"
         B="13"
         E="15"
         BS="N"
        /><b
         id="25"
         Z="19"
         B="15"
         E="17"
         BS="N"
        /><b
         id="26"
         Z="20"
         B="17"
         E="19"
         BS="N"
        /><b
         id="27"
         Z="21"
         B="19"
         E="7"
         BS="N"
        /><b
         id="29"
         Z="23"
         B="15"
         E="28"
         BS="N"
        /></fragment></page></CDXML>)";

    {
      std::stringstream iss(cdxml1);
      auto mols = MolsFromCDXMLDataStream(iss);
      for (auto &mol : mols) {
        CHECK(MolToSmiles(*mol) == "CC(C)(C)OC(=O)C1CCCCCC1");
      }
    }
    {
      std::stringstream iss(cdxml1);
      // v1 api
      auto mols = CDXMLDataStreamToMols(iss);
      for (auto &mol : mols) {
        CHECK(MolToSmiles(*mol) == "CC(C)(C)OC(=O)C1CCCCCC1");
      }
    }
  }
  SECTION("REACTION") {
    std::string fname = cdxmlbase + "reaction-with-boc.cdxml";
    std::vector<std::string> expected = {"CC(C)(C)OC(=O)C1CCCCCC1[*:1]",
                                         "c1ccc([*:1])cc1", "C1CC1", "C1CCC1"};
    {
      auto mols = MolsFromCDXMLFile(fname);
      CHECK(mols.size() == expected.size());
      int i = 0;
      for (auto &mol : mols) {
        CHECK(MolToSmiles(*mol) == expected[i++]);
      }
    }
    {
      // v1 api
      auto mols = CDXMLFileToMols(fname);
      CHECK(mols.size() == expected.size());
      int i = 0;
      for (auto &mol : mols) {
        CHECK(MolToSmiles(*mol) == expected[i++]);
      }
    }
  }
  SECTION("RING CHIRALITY") {
    std::string fname = cdxmlbase + "ring-stereo1.cdxml";
    std::vector<std::string> expected = {"C1CC[C@H]2CCCC[C@H]2C1"};
    auto mols = MolsFromCDXMLFile(fname);
    CHECK(mols.size() == expected.size());
    int i = 0;
    for (auto &mol : mols) {
      CHECK(MolToSmiles(*mol) == expected[i++]);
    }
  }
  SECTION("SIMPLE CHIRAL") {
    std::string fname = cdxmlbase + "chirality1.cdxml";
    std::vector<std::string> expected = {"C[C@H](N)C[C@H](C)N"};
    auto mols = MolsFromCDXMLFile(fname);
    CHECK(mols.size() == expected.size());
    int i = 0;
    for (auto &mol : mols) {
      CHECK(MolToSmiles(*mol) == expected[i++]);
    }
  }
  SECTION("CDXML-CISTRANS") {
    auto fname = cdxmlbase + "cistrans1.cdxml";
    std::vector<std::string> expected = {"F/C(I)=C(\\Cl)Br"};
    auto mols = MolsFromCDXMLFile(fname);
    CHECK(mols.size() == expected.size());
    int i = 0;
    for (auto &mol : mols) {
      CHECK(MolToSmiles(*mol) == expected[i++]);
    }
  }
  SECTION("REACTION") {
    auto fname = cdxmlbase + "rxn1.cdxml";
    std::vector<std::string> expected = {
        "Cl[c:1]1[cH:4][cH:3][cH:2][cH:6][cH:5]1",
        "OC(O)B[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1",
        "[cH:1]1[cH:4][cH:3][cH:2][c:6](-[c:7]2[cH:8][cH:9][cH:10][cH:11][cH:12]2)[cH:5]1"};
    auto mols = MolsFromCDXMLFile(fname);
    int i = 0;
    for (auto &mol : mols) {
      CHECK(mol->getProp<unsigned int>("CDX_SCHEME_ID") == 397);
      CHECK(mol->getProp<unsigned int>("CDX_STEP_ID") == 398);
      if (i == 0) {
        CHECK(mol->getProp<unsigned int>("CDX_REAGENT_ID") == 0);
      } else if (i == 1) {
        CHECK(mol->getProp<unsigned int>("CDX_REAGENT_ID") == 1);
      } else if (i == 2) {
        CHECK(mol->getProp<unsigned int>("CDX_PRODUCT_ID") == 0);
      }
      CHECK(MolToSmiles(*mol) == expected[i++]);
    }
  }
  SECTION("DEUTERIUM") {
    auto fname = cdxmlbase + "deuterium.cdxml";
    {
      std::vector<std::string> expected = {
          "[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]"};
      auto mols = MolsFromCDXMLFile(fname);
      CHECK(mols.size() == expected.size());
      int i = 0;
      for (auto &mol : mols) {
        CHECK(MolToSmiles(*mol) == expected[i++]);
      }
    }
    {
      std::vector<std::string> expected = {
          "[2H]C1=C([2H])C([2H])=C([2H])C([2H])=C1[2H]"};
      CDXMLParserParams params;
      params.sanitize = false;
      params.removeHs = false;
      auto mols = MolsFromCDXMLFile(fname, params);
      CHECK(mols.size() == expected.size());
      int i = 0;
      for (auto &mol : mols) {
        CHECK(MolToSmiles(*mol) == expected[i++]);
      }
    }
    {
      std::vector<std::string> expected = {
          "[2H]C1=C([2H])C([2H])=C([2H])C([2H])=C1[2H]"};
      CDXMLParserParams params;
      params.sanitize = false;
      params.removeHs = true;
      auto mols = MolsFromCDXMLFile(fname, params);
      CHECK(mols.size() == expected.size());
      int i = 0;
      for (auto &mol : mols) {
        CHECK(MolToSmiles(*mol) == expected[i++]);
      }
    }
  }
  SECTION("Queries") {
    {
      auto fname = cdxmlbase + "query-atoms.cdxml";

      std::vector<std::string> expected = {"*c1ccccc1", "*c1ccccc1",
                                           "*c1ccccc1"};
      std::vector<std::string> expected_smarts = {
          "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-*",
          "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[!#1]",
          "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[!#6&!#1]"};
      auto mols = MolsFromCDXMLFile(fname);
      CHECK(mols.size() == expected.size());
      int i = 0;
      for (auto &mol : mols) {
        CHECK(MolToSmarts(*mol) == expected_smarts[i]);
        CHECK(MolToSmiles(*mol) == expected[i++]);
      }
    }
    {
      auto fname = cdxmlbase + "anybond.cdxml";
      auto mols = MolsFromCDXMLFile(fname);
      CHECK(mols.size() == 1);
      CHECK(MolToSmiles(*mols[0]) == "C1CCC~CC1");
      CHECK(MolToSmarts(*mols[0]) == "[#6]1~[#6]-[#6]-[#6]-[#6]-[#6]-1");
    }
  }
  SECTION("ElementList") {
    auto fname = cdxmlbase + "element-list.cdxml";

    std::vector<std::string> expected = {"[C]CC"};
    std::vector<std::string> expected_smarts = {"[#6]-[#6]-[#6,#7,#8,#16]"};
    auto mols = MolsFromCDXMLFile(fname);
    CHECK(mols.size() == expected.size());
    int i = 0;
    for (auto &mol : mols) {
      CHECK(MolToSmarts(*mol) == expected_smarts[i]);
      CHECK(MolToSmiles(*mol) == expected[i++]);
    }
  }
  SECTION("Enhanced Stereo") {
    auto fname = cdxmlbase + "beta-cypermethrin.cdxml";
    std::vector<std::string> expected = {
        "CC1(C)[C@H](C=C(Cl)Cl)[C@H]1C(=O)O[C@@H](C#N)c1cccc(Oc2ccccc2)c1"};
    std::vector<std::string> expected_cx = {
        "CC1(C)[C@H](C=C(Cl)Cl)[C@@H]1C(=O)O[C@H](C#N)c1cccc(Oc2ccccc2)c1 |&1:3,&2:8,12|"};
    auto mols = MolsFromCDXMLFile(fname);
    CHECK(mols.size() == expected.size());
    int i = 0;
    for (auto &mol : mols) {
      mol.get()->clearConformers();
      CHECK(MolToSmiles(*mol) == expected[i]);
      CHECK(MolToCXSmiles(*mol) == expected_cx[i++]);
    }
  }
  SECTION("Enhanced Stereo 2") {
    auto fname = cdxmlbase + "beta-cypermethrin-or-abs.cdxml";
    std::vector<std::string> expected = {
        "CC1(C)[C@H](C=C(Cl)Cl)[C@H]1C(=O)O[C@@H](C#N)c1cccc(Oc2ccccc2)c1"};
    std::vector<std::string> expected_cx = {
        "CC1(C)[C@H](C=C(Cl)Cl)[C@@H]1C(=O)O[C@H](C#N)c1cccc(Oc2ccccc2)c1 |a:3,o1:8,12|"};
    auto mols = MolsFromCDXMLFile(fname);
    CHECK(mols.size() == expected.size());
    int i = 0;
    SmilesWriteParams wp;
    for (auto &mol : mols) {
      auto tomol = std::unique_ptr<ROMol>(mol.release());
      tomol.get()->clearConformers();
      RDKit::canonicalizeStereoGroups(tomol);

      CHECK(MolToSmiles(*tomol) == expected[i]);
      CHECK(MolToCXSmiles(*tomol, wp) == expected_cx[i++]);
    }
  }
  SECTION("Bad CDXML") {
    auto fname = cdxmlbase + "bad-cdxml.cdxml";
    // Only one passes sanitization
    {
      std::vector<std::string> expected = {"*c1ccccc1"};
      std::vector<std::string> expected_smarts = {
          "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-*",
      };
      auto mols = MolsFromCDXMLFile(fname);
      CHECK(mols.size() == expected.size());
      int i = 0;
      for (auto &mol : mols) {
        CHECK(MolToSmarts(*mol) == expected_smarts[i]);
        CHECK(MolToSmiles(*mol) == expected[i++]);
      }
    }
    // setting sanitization to false, we get both
    std::vector<std::string> expected = {"*C1=C([H])C([H])=C([H])C([H])=C1[H]",
                                         "*C1=C([H])N([H])=C([H])C([H])=C1[H]"};
    std::vector<std::string> expected_smarts = {
        "[#6]1(=[#6](-[#6](=[#6](-[#6](=[#6]-1-*)-[H])-[H])-[H])-[H])-[H]",
        "[#6]1(=[#6](-[#6](=[#7](-[#6](=[#6]-1-[!#1])-[H])-[H])-[H])-[H])-[H]",
    };
    CDXMLParserParams params;
    params.sanitize = false;
    auto mols = MolsFromCDXMLFile(fname, params);
    CHECK(mols.size() == expected.size());
    int i = 0;
    for (auto &mol : mols) {
      CHECK(MolToSmarts(*mol) == expected_smarts[i]);
      CHECK(MolToSmiles(*mol) == expected[i++]);
    }
  }

  SECTION("Fusion with chirality") {
    auto fname = cdxmlbase + "fusion-chiral.cdxml";
    auto mols = MolsFromCDXMLFile(fname);
    std::vector<std::string> expected = {"C[C@@H](O)[C@@H](C)O"};
    CHECK(mols.size() == expected.size());
    int i = 0;
    for (auto &mol : mols) {
      CHECK(MolToSmiles(*mol) == expected[i++]);
    }
  }

  SECTION("ChemDraw Template from the synthesis-workshop") {
    // this was hella fun to validate the stereo-chemistry...
    auto fname = cdxmlbase + "chemdraw_template1.cdxml";
    auto mols = MolsFromCDXMLFile(fname);
#ifndef RDK_BUILD_CHEMDRAW_SUPPORT
    std::vector<std::string> expected = {
        "CCC/C=C/C=C/C(=O)O[C@H]1/C(=C/C(=O)OC)C[C@H]2C[C@H]([C@@H](C)O)OC(=O)C[C@H](O)C[C@@H]3C[C@H](OC(C)=O)C(C)(C)[C@](O)(C[C@@H]4C/C(=C/C(=O)OC)C[C@H](/C=C/C(C)(C)[C@]1(O)O2)O4)O3",
        "[B]",
        "*",
        "[C]",
        "Cc1ccc2n1[C@@H]1[C@@H]3O[C@]([C@H](C)O)(C=C2)[C@H]1c1ccc(C)n1[C@@H]3C",  // this is may or may not be correct, but the structure is drawn incorrectly. There's a test below which fixes this
        "Cc1ccc2n1[C@H](C)C(=O)[C@@H]1[C@H]2C(=O)C=Cc2ccc(C)n21",
        "Cc1ccc2ccc(=O)ccn12",
        "Cc1cccn1[C@H](C)C=O",
        "Cc1ccc2ccc([O-])cc[n+]1-2",
        "Cc1ccc2ccc(=O)ccn12",
        "Cc1cccn1[C@H](C)C(C#N)O[Si](C)(C)C",
        "CC1CCC2(O)C3(OC4(O)C[C@]2(C)C2(O)[C@H](OC(=O)c5ccc[nH]5)C(O)(C(C)C)C4(C)C32O)C1O",
        "C=C(C)[C@H]1CC(=O)CC2=C(C1)[C@H]1C(=O)O[C@H]3C[C@@](C)(O)[C@@H](C2=O)[C@@H]13"};
#else
    // the new cdxml parser handles stereo a lot better
    std::vector<std::string> expected = {
        "CCC/C=C/C=C/C(=O)O[C@H]1/C(=C/C(=O)OC)C[C@H]2C[C@H]([C@@H](C)O)OC(=O)C[C@H](O)C[C@@H]3C[C@H](OC(C)=O)C(C)(C)[C@](O)(C[C@@H]4C/C(=C/C(=O)OC)C[C@H](/C=C/C(C)(C)[C@]1(O)O2)O4)O3",
        "[B]",
	"*",
	"[C]",
        "Cc1ccc2n1[C@@H]1[C@@H]3O[C@]([C@H](C)O)(C=C2)[C@H]1c1ccc(C)n1[C@@H]3C",
        // this is may or may not be correct, but the structure is drawn
        // incorrectly.
        // There's a test below which fixes this
        "Cc1ccc2n1[C@H](C)C(=O)[C@@H]1[C@H]2C(=O)C=Cc2ccc(C)n21",
        "Cc1ccc2ccc(=O)ccn12",
	"Cc1cccn1[C@H](C)C=O",
        "Cc1ccc2ccc([O-])cc[n+]1-2",
	"Cc1ccc2ccc(=O)ccn12",
        "Cc1cccn1[C@H](C)C(C#N)O[Si](C)(C)C",
        "CC1CC[C@]2(O)[C@]3(C)C[C@]4(O)O[C@@]2([C@@H]1O)C1(O)C4(C)C(O)(C(C)C)[C@@H](OC(=O)c2ccc[nH]2)[C@]13O",
        "C=C(C)[C@H]1CC(=O)CC2=C(C1)[C@H]1C(=O)O[C@H]3C[C@@](C)(O)[C@@H](C2=O)[C@@H]13"};
    
#endif
    CHECK(mols.size() == expected.size());
    int i = 0;
    for (auto &mol : mols) {
      CHECK(MolToSmiles(*mol) == expected[i++]);
    }
  }

  SECTION("deuterium atom") {
    auto fname = cdxmlbase + "deuterium-atom.cdxml";
    CDXMLParserParams params;
    params.sanitize = false;
    params.removeHs = false;
    auto mols = MolsFromCDXMLFile(fname, params);
    std::vector<std::string> expected = {"[2H]"};
    CHECK(mols.size() == expected.size());
    int i = 0;
    for (auto &mol : mols) {
      CHECK(MolToSmiles(*mol) == expected[i++]);
    }
  }
  SECTION("ChemDraw Template 2 from the synthesis-workshop") {
    // this was another hella fun to validate the stereo-chemistry...
    //   there were so many stereo warnings in chemdraw, I'm just going to
    //   assume
    //    the rdkit is correct here...
    auto fname = cdxmlbase + "chemdraw_template2.cdxml";
    auto mols = MolsFromCDXMLFile(fname);
#ifndef RDK_BUILD_CHEMDRAW_SUPPORT    
    std::vector<std::string> expected = {
        "CCN1CC2(COC)CCC(OC)C34C5CC6C(OC)CC(O)(C(CC23)C14)C5C6O",
        "*",
        "C",
        "[F]",
        "[B]",
        "[C]",
        "[2H]",
        "CCN1CC2(COC)CCC(OC)C34C5CC6C(OC)CC(O)(C(CC23)C14)C5C6O",
        "*",
        "C",
        "[F]",
        "[B]",
        "[C]",
        "[2H]",
        "CCN1CC2(COC)CCC(OC)C34C5CC6C(OC)CC(O)(C(CC23)C14)C5C6O",
        "*",
        "C",
        "[F]",
        "[B]",
        "[C]",
        "[2H]",
        "CCN1CC2(COC)CCC(OC)C34C5CC6C(OC)CC(O)(C(CC23)C14)C5C6O",
        "*",
        "C",
        "[F]",
        "[B]",
        "[C]",
        "[2H]",
        "CCN1CC2(COC)CCC(OC)C34C5CC6C(OC)CC(O)(C(CC23)C14)C5C6O",
        "CCN1CC2(COC)CCC(OC)C34C5CC6C(OC)CC(O)(C(CC23)C14)C5C6O",
        "*",
        "[B]",
        "[C]",
        "[2H]",
        "C",
        "[F]",
        "*",
        "C",
        "[F]",
        "[B]",
        "[C]",
        "[2H]",
        "CCN1CC2(COC)CCC(OC)C34C5CC6C(OC)CC(O)(C(CC23)C14)C5C6O",
        "*",
        "C",
        "[F]",
        "[B]",
        "[C]",
        "[2H]",
        "CCN1CC2(COC)CCC(OC)C34C5CC6C(OC)CC(O)(C(CC23)C14)C5C6O",
        "*",
        "C",
        "[F]",
        "[B]",
        "[C]",
        "[2H]",
        "CC1CCC2(O)C3(OC4(O)C[C@]2(C)C2(O)[C@H](O)C(O)(C(C)C)C4(C)C32O)C1O",
        "CC1=C(C(C)C)[C@@H](O)C2(O)C1(O)C13OC(=O)C[C@@]2(C)C1(O)CCC(C)C3O",
        "CC1=CC23OC(=O)C[C@@](C)(C2(O)CC1)C1(O)[C@H](O)C2(C(C)C)OC2(C)C31O",
        "*",
        "[B]",
        "[C]",
        "CC1CCC2C3(C1)OC1C[C@]2(C)C2CC(C(C)C)C1(C)C23",
        "[2H]",
        "*",
        "[B]",
        "[C]",
        "C",
        "CC1CCC2(O)C3(OC4(O)C[C@]2(C)C2(O)[C@H](O)C(O)(C(C)C)C4(C)C32O)C1O",
        "[2H]"};
#else
    // The new cdxml parser handles stereo a LOT better
    std::string talatisamine = "CCN1C[C@]2(COC)CCC(OC)[C@@]34[C@@H]5C[C@@H]6C(OC)C[C@@](O)([C@H]5[C@H]6O)[C@@H](C[C@H]23)[C@H]14";
    std::vector<std::string> expected = {
      talatisamine, //0 
        "*",
        "C",
        "[F]",
        "[B]",
        "[C]",
        "[2H]",
        talatisamine,
        "*",
        "C",
        "[F]", // 10
        "[B]",
        "[C]",
        "[2H]",
	talatisamine,
        "*",
        "C",
        "[F]",
        "[B]",
        "[C]",
        "[2H]", // 20
        talatisamine,
        "*",
        "C",
        "[F]",
        "[B]",
        "[C]",
        "[2H]",
        talatisamine,
	"CCN1C[C@]2(COC)CC[C@H](OC)[C@]34C1C(C[C@H]23)[C@@]1(O)CC(OC)[C@H]2C[C@@H]4[C@@H]1[C@H]2O",
        "*", // 30
        "[B]",
        "[C]",
        "[2H]",
        "C",
        "[F]",
        "*",
        "C",
        "[F]",
        "[B]",
        "[C]", // 40
        "[2H]",
        talatisamine,
        "*",
        "C",
        "[F]",
        "[B]",
        "[C]",
        "[2H]",
        talatisamine,
        "*", // 50
        "C",
        "[F]",
        "[B]",
        "[C]",
        "[2H]",
        "CC1CC[C@]2(O)[C@]3(C)C[C@]4(O)O[C@@]2([C@@H]1O)C1(O)C4(C)C(O)(C(C)C)[C@@H](O)[C@]13O",
        "CC1=C(C(C)C)[C@@H](O)[C@@]2(O)[C@@]3(C)CC(=O)O[C@@]4([C@H](O)C(C)CC[C@]34O)[C@@]12O",
        "CC1=C[C@@]23OC(=O)C[C@@](C)([C@@]2(O)CC1)[C@]1(O)[C@H](O)C2(C(C)C)OC2(C)[C@@]31O",
        "*",
        "[B]", // 60
        "[C]",
        "CC1CC[C@@H]2[C@]3(C)C[C@@H]4O[C@@]2(C1)C1[C@@H]3CC(C(C)C)C14C",
        "[2H]",
        "*",
        "[B]",
        "[C]",
        "C",
        "CC1CC[C@]2(O)[C@]3(C)C[C@]4(O)O[C@@]2([C@@H]1O)C1(O)C4(C)C(O)(C(C)C)[C@@H](O)[C@]13O",
        "[2H]"};
#endif
    CHECK(mols.size() == expected.size());
    int i = 0;
    for (auto &mol : mols) {
      INFO(i);
      CHECK(MolToSmiles(*mol) == expected[i++]);
    }
  }
  SECTION("ChemDraw Template 3 from the synthesis-workshop") {
    // this was another hella fun to validate the stereo-chemistry...
    //   there were so many stereo warnings in chemdraw, I'm just going to
    //   assume
    //    the rdkit is correct here...
    auto fname = cdxmlbase + "chemdraw_template3.cdxml";
    auto mols = MolsFromCDXMLFile(fname);
    std::vector<std::string> expected = {
        "CCC/C=C/C=C/C(=O)O[C@H]1/C(=C/C(=O)OC)C[C@H]2C[C@H]([C@@H](C)O)OC(=O)C[C@H](O)C[C@@H]3C[C@H](OC(C)=O)C(C)(C)[C@](O)(C[C@@H]4C/C(=C/C(=O)OC)C[C@H](/C=C/C(C)(C)[C@]1(O)O2)O4)O3",
        "[B]",
        "*",
        "[C]",
        "CCC/C=C/C=C/C(=O)O[C@H]1/C(=C/C(=O)OC)C[C@H]2C[C@H]([C@@H](C)O)OC(=O)C[C@H](O)C[C@@H]3C[C@H](OC(C)=O)C(C)(C)[C@](O)(C[C@@H]4C/C(=C/C(=O)OC)C[C@H](/C=C/C(C)(C)[C@]1(O)O2)O4)O3",
        "[B]",
        "*",
        "[C]",
        "CCC/C=C/C=C/C(=O)O[C@H]1/C(=C/C(=O)OC)C[C@@H](C[C@@H](O)[C@@H](C)O)O[C@@]1(O)C(C)(C)/C=C/C=O",
        "*",
        "[C]",
        "C=C(C[C@H]([O])C[C@]1(O)O[C@H](C[C@@H](O)CC(=O)O)C[C@H](OC(C)=O)C1(C)C)C[Si](C)(C)C",
        "*.CC[Si](CC)CC",
        "CC[Si](C)(CC)CC",
        "CC[Si](C)(CC)CC",
        "CC",
        "CC",
        "*",
        "C=C(C[C@H]([O])C[C@]1(O)O[C@H](C[C@@H](O)CC(=O)O)C[C@H](OC(C)=O)C1(C)C)C[Si](C)(C)C",
        "*.CC[Si](CC)CC",
        "CCC/C=C/C=C/C(=O)O[C@H]1/C(=C/C(=O)OC)C[C@@H](C[C@@H](O)[C@@H](C)O)O[C@@]1(O)C(C)(C)/C=C/C=O",
        "[C]"};
    int i = 0;
    for (auto &mol : mols) {
      INFO(i);
      check_smiles_and_roundtrip(*mol, expected[i++]);
    }
  }
  SECTION("protecting group") {
    auto fname = cdxmlbase + "protecting-groups.cdxml";
    auto mols = MolsFromCDXMLFile(fname);
    std::vector<std::string> expected = {"CC[Si](C)(CC)CC", "CC"};
    CHECK(mols.size() == expected.size());
    int i = 0;
    for (auto &mol : mols) {
      check_smiles_and_roundtrip(*mol, expected[i++]);
    }
  }
  SECTION("protecting group 2") {
    auto fname = cdxmlbase + "protecting-groups2.cdxml";
    auto mols = MolsFromCDXMLFile(fname);
    std::vector<std::string> expected = {"CC[Si](C)(CC)CC", "CC"};
    CHECK(mols.size() == expected.size());
    int i = 0;
    for (auto &mol : mols) {
      check_smiles_and_roundtrip(*mol, expected[i++]);
    }
  }

  SECTION("floating protecting group") {
    auto fname = cdxmlbase + "floating-protecting-group.cdxml";
    // This is a weird one, chemdraw simply ignores the error that causes the
    // bond issue, we should probably drop the floating fragment here if
    //  we are sanitizing
    auto mols = MolsFromCDXMLFile(fname);
    std::vector<std::string> expected = {
        "*",
        "C=C(C[C@H]([O])C[C@]1(O)O[C@H](C[C@@H](O)CC(=O)O)C[C@H](OC(C)=O)C1(C)C)C[Si](C)(C)C",
        "*.CC[Si](CC)CC"};
    CHECK(mols.size() == expected.size());
    int i = 0;
    for (auto &mol : mols) {
      INFO(i);
      check_smiles_and_roundtrip(*mol, expected[i++]);
    }
  }

  SECTION("Missing File Name") {
    try {
      auto mols = MolsFromCDXMLFile("missing file");
      CHECK(0);  // Bad file exception not caught
    } catch (RDKit::BadFileException &) {
    }
  }

  SECTION("Aromatic ring (bondorder==4") {
    auto fname = cdxmlbase + "aromatic.cdxml";
    auto mols = MolsFromCDXMLFile(fname);
    std::vector<std::string> expected = {"c1ccccc1"};
    CHECK(mols.size() == expected.size());
    int i = 0;
    for (auto &mol : mols) {
      check_smiles_and_roundtrip(*mol, expected[i++]);
    }
  }
  SECTION("Malformed") {
    auto fname = cdxmlbase + "malformed.cdxml";
    try {
      auto mols = MolsFromCDXMLFile(fname);
      CHECK(0);
    } catch (FileParseException &e) {
#ifndef RDK_BUILD_CHEMDRAW_SUPPORT      
      CHECK(std::string(e.what()) == "expected > at line: 373");
#else
      CHECK(std::string(e.what()).find("Failed parsing XML with error code 5") !=
	    std::string::npos);
#endif
    }
  }
  SECTION("Lots of stereo") {
    {
      auto fname = cdxmlbase + "stereo.cdxml";
      std::vector<std::string> expected = {
          "C[C@@H](Cl)[C@H](N)O.C[C@@H](F)[C@H](N)O.C[C@H](Br)[C@@H](N)O.C[C@H](I)[C@@H](N)O"};
      auto mols = MolsFromCDXMLFile(fname);
      CHECK(mols.size() == expected.size());
      int i = 0;
      for (auto &mol : mols) {
        check_smiles_and_roundtrip(*mol, expected[i++]);
      }
    }

    {  // The above, but broken out for easier testing
      std::vector<std::string> filenames = {"stereo1.cdxml", "stereo2.cdxml",
                                            "stereo3.cdxml", "stereo4.cdxml"};
      std::vector<std::string> expected = {
          "C[C@H](I)[C@@H](N)O", "C[C@@H](I)[C@H](N)O", "C[C@@H](Cl)[C@H](N)O",
          "C[C@H](Br)[C@@H](N)O"};

      for (auto i = 0u; i < filenames.size(); ++i) {
        auto fname = cdxmlbase + filenames[i];
        auto mols = MolsFromCDXMLFile(fname);
        CHECK(mols.size() == 1);
        auto &m = *mols.back();
        check_smiles_and_roundtrip(m, expected[i]);
      }
    }

    {
      auto fname = cdxmlbase + "wavy.cdxml";
      std::vector<std::string> expected = {"Cc1cccc(C)c1NC(=O)N=C1CCCN1C",
                                           "Cc1cccc(C)c1NC(=O)/N=C1\\CCCN1C"};
      auto mols = MolsFromCDXMLFile(fname);
      CHECK(mols.size() == expected.size());
      int i = 0;
      for (auto &mol : mols) {
        if (i == 0) {
          CHECK(mol->getBondWithIdx(11)->getStereo() ==
                Bond::BondStereo::STEREOANY);
        }
        check_smiles_and_roundtrip(*mol, expected[i++]);
      }
    }
    {
      auto fname = cdxmlbase + "wavy-single.cdxml";
      std::vector<std::string> expected = {"CCCC"};
      auto mols = MolsFromCDXMLFile(fname);
      CHECK(mols.size() == expected.size());
      int i = 0;
      for (auto &mol : mols) {
        CHECK(mol->getBondWithIdx(0)->getBondDir() == Bond::BondDir::NONE);
        CHECK(mol->getBondWithIdx(1)->getBondDir() == Bond::BondDir::NONE);
        CHECK(mol->getBondWithIdx(2)->getBondDir() == Bond::BondDir::NONE);
        check_smiles_and_roundtrip(*mol, expected[i++]);
      }
    }
  }
  SECTION("Lots of bad molecules") {
    {
      auto fname = cdxmlbase + "bad-id.cdxml";
      auto mols = MolsFromCDXMLFile(fname);
      CHECK(mols.size() == 0);
    }
    {
      auto fname = cdxmlbase + "bad-coords.cdxml";
      auto mols = MolsFromCDXMLFile(fname);
      CHECK(mols.size() == 0);
    }
    {
      auto fname = cdxmlbase + "bad-bondorder.cdxml";
      auto mols = MolsFromCDXMLFile(fname);
#ifndef RDK_BUILD_CHEMDRAW_SUPPORT      
      CHECK(mols.size() == 0);
#else      
      CHECK(mols.size() == 1); // The original chemdraw reader makes unknowns single bonds
#endif
    }
    {
      auto fname = cdxmlbase + "bad-bondorder2.cdxml";
      auto mols = MolsFromCDXMLFile(fname);
      CHECK(mols.size() == 0);
    }
  }
}

TEST_CASE("atropisomers") {
  std::string cdxmlbase =
      std::string(getenv("RDBASE")) + "/Code/GraphMol/test_data/CDXML/";

  SECTION("atropisomer") {
    {
      std::vector<std::string> filenames = {"atrop1.cdxml"};

#ifndef RDK_BUILD_CHEMDRAW_SUPPORT            
      std::vector<std::string> expected = {
          "C[C]1[C][CH]C(Cl)C(C)=C1c1c(C)ccc(Cl)c1C |(-2.936,-0.12,;-2.936,-1.66,;-1.602,-2.43,;-1.602,-3.97,;-2.936,-4.74,;-2.93,-6.28,;-4.27,-3.97,;-5.603,-4.74,;-4.27,-2.43,;-5.603,-1.66,;-5.603,-0.12,;-4.27,0.64,;-6.937,0.64,;-8.271,-0.12,;-8.271,-1.66,;-9.604,-2.43,;-6.937,-2.43,;-6.937,-3.97,),^1:1,3,^2:2,wU:8.8|"};
#else
      std::vector<std::string> expected = {
	"C[C]1[C][CH]C(Cl)C(C)=C1c1c(C)ccc(Cl)c1C |(-2.936,-0.12,;-2.936,-1.66,;-1.602,-2.43,;-1.602,-3.97,;-2.936,-4.74,;-2.93,-6.28,;-4.27,-3.97,;-5.603,-4.74,;-4.27,-2.43,;-5.603,-1.66,;-5.603,-0.12,;-4.27,0.639999,;-6.937,0.639999,;-8.271,-0.12,;-8.271,-1.66,;-9.604,-2.43,;-6.937,-2.43,;-6.937,-3.97,),^1:1,3,^2:2,wU:8.8|"};
      
#endif
      for (auto i = 0u; i < filenames.size(); ++i) {
        auto fname = cdxmlbase + filenames[i];
        auto mol = MolsFromCDXMLFile(fname);

        SmilesWriteParams ps;
        auto smi = MolToCXSmiles(*(mol[0].get()), ps,
                                 SmilesWrite::CXSmilesFields::CX_ALL);
        CHECK(smi == expected[i]);
      }
    }
  }
}

TEST_CASE("bad stereo in a natural product") {
  std::string cdxmlbase =
      std::string(getenv("RDBASE")) + "/Code/GraphMol/test_data/CDXML/";
  SECTION("case 1") {
    auto fname = cdxmlbase + "stereo5.cdxml";
    auto mols = MolsFromCDXMLFile(fname);
    REQUIRE(mols.size() == 1);
    CHECK(
        MolToSmiles(*mols[0]) ==
        "Cc1ccc2n1[C@@H]1[C@@H]3O[C@]([C@H](C)O)(C=C2)[C@H]1c1ccc(C)n1[C@@H]3C");
  }
}

TEST_CASE("Github #6262: preserve bond wedging") {
  std::string cdxmlbase =
      std::string(getenv("RDBASE")) + "/Code/GraphMol/test_data/CDXML/";
  SECTION("case 1") {
    auto fname = cdxmlbase + "stereo6.cdxml";
    auto mols = MolsFromCDXMLFile(fname);
    REQUIRE(mols.size() == 1);
    {
      REQUIRE(mols[0]->getBondBetweenAtoms(2, 5));
      unsigned int cfg = 0;
      CHECK(mols[0]->getBondBetweenAtoms(2, 5)->getPropIfPresent(
          common_properties::_MolFileBondCfg, cfg));
      CHECK(cfg == 1);
    }
    {
      REQUIRE(mols[0]->getBondBetweenAtoms(1, 4));
      unsigned int cfg = 0;
      CHECK(mols[0]->getBondBetweenAtoms(1, 4)->getPropIfPresent(
          common_properties::_MolFileBondCfg, cfg));
      CHECK(cfg == 3);
    }
    {
      REQUIRE(mols[0]->getBondBetweenAtoms(3, 8));
      unsigned int cfg = 0;
      CHECK(mols[0]->getBondBetweenAtoms(3, 8)->getPropIfPresent(
          common_properties::_MolFileBondCfg, cfg));
      CHECK(cfg == 2);
    }
  }
}

TEST_CASE("Github #6887: and1 or1 in same mol") {
  SECTION("case 1") {
    std::string cdxml1 = R"(<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >
<CDXML
 CreationProgram="ChemDraw 21.0.0.28"
 Name="Untitled A4 Document-2"
 BoundingBox="115.75 357.32 272.13 406.58"
 WindowPosition="0 0"
 WindowSize="0 0"
 FractionalWidths="yes"
 InterpretChemically="yes"
 ShowAtomQuery="yes"
 ShowAtomStereo="no"
 ShowAtomEnhancedStereo="yes"
 ShowAtomNumber="no"
 ShowResidueID="no"
 ShowBondQuery="yes"
 ShowBondRxn="yes"
 ShowBondStereo="no"
 ShowTerminalCarbonLabels="no"
 ShowNonTerminalCarbonLabels="no"
 HideImplicitHydrogens="no"
 LabelFont="21"
 LabelSize="10"
 LabelFace="96"
 CaptionFont="510"
 CaptionSize="12"
 HashSpacing="2.70"
 MarginWidth="2"
 LineWidth="1"
 BoldWidth="4"
 BondLength="30"
 BondSpacing="12"
 ChainAngle="120"
 LabelJustification="Auto"
 CaptionJustification="Left"
 AminoAcidTermini="HOH"
 ShowSequenceTermini="yes"
 ShowSequenceBonds="yes"
 ShowSequenceUnlinkedBranches="no"
 ResidueWrapCount="40"
 ResidueBlockCount="10"
 ResidueZigZag="yes"
 NumberResidueBlocks="no"
 PrintMargins="36 36 36 36"
 color="0"
 bgcolor="1"
><colortable>
<color r="1" g="1" b="1"/>
<color r="0" g="0" b="0"/>
<color r="1" g="0" b="0"/>
<color r="1" g="1" b="0"/>
<color r="0" g="1" b="0"/>
<color r="0" g="1" b="1"/>
<color r="0" g="0" b="1"/>
<color r="1" g="0" b="1"/>
</colortable><fonttable>
<font id="21" charset="x-mac-roman" name="Helvetica"/>
<font id="510" charset="x-mac-roman" name="Times New Roman"/>
</fonttable><page
 id="24"
 BoundingBox="0 0 523 770"
 HeaderPosition="36"
 FooterPosition="36"
 PrintTrimMarks="yes"
 HeightPages="1"
 WidthPages="1"
><fragment
 id="3"
 BoundingBox="115.75 357.32 272.13 406.58"
 Z="2"
><n
 id="2"
 p="116 406"
 Z="1"
 AS="N"
/><n
 id="4"
 p="141.98 391"
 Z="3"
 Geometry="Tetrahedral"
 AS="N"
 BondOrdering="5 7 0 17"
 EnhancedStereoType="And"
 EnhancedStereoGroupNum="1"
><objecttag
 TagType="Unknown"
 Name="enhancedstereo"
><t
 p="137.96 401.37"
 BoundingBox="138.29 396 145.62 401.50"
 CaptionLineHeight="variable"
><s font="21" size="7.5" color="0">&amp;1</s></t></objecttag></n><n
 id="6"
 p="167.96 406"
 Z="5"
 AS="N"
/><n
 id="8"
 p="193.94 391"
 Z="7"
 Geometry="Tetrahedral"
 AS="N"
 BondOrdering="9 11 0 19"
 EnhancedStereoType="Or"
 EnhancedStereoGroupNum="2"
><objecttag
 TagType="Unknown"
 Name="enhancedstereo"
><t
 p="188.54 401.26"
 BoundingBox="188.76 396 199.07 401.41"
 CaptionLineHeight="variable"
><s font="21" size="7.5" color="0">or2</s></t></objecttag></n><n
 id="10"
 p="219.92 406"
 Z="9"
 AS="N"
/><n
 id="12"
 p="245.90 391"
 Z="11"
 Geometry="Tetrahedral"
 AS="N"
 BondOrdering="13 15 0 21"
 EnhancedStereoType="Or"
 EnhancedStereoGroupNum="1"
><objecttag
 TagType="Unknown"
 Name="enhancedstereo"
><t
 p="241.11 401.22"
 BoundingBox="241.32 396 250.43 401.36"
 CaptionLineHeight="variable"
><s font="21" size="7.5" color="0">or1</s></t></objecttag></n><n
 id="14"
 p="271.88 406"
 Z="13"
 AS="N"
/><n
 id="16"
 p="141.98 361"
 Z="15"
 NodeType="Fragment"
 NeedsClean="yes"
 AS="N"
><fragment
 id="25"
><n
 id="26"
 p="141.98 361"
 Element="8"
 NumHydrogens="0"
/><n
 id="27"
 p="116 346"
 NumHydrogens="3"
/><n
 id="28"
 p="141.98 391"
 NodeType="ExternalConnectionPoint"
/><b
 id="29"
 B="26"
 E="27"
/><b
 id="30"
 B="28"
 E="26"
 Display="WedgeBegin"
/></fragment><t
 p="138.09 364.68"
 BoundingBox="138.48 357.32 159.33 364.89"
 LabelJustification="Left"
 LabelAlignment="Left"
><s font="21" size="10" color="0" face="96">OMe</s></t></n><n
 id="18"
 p="193.94 361"
 Z="17"
 Element="17"
 NumHydrogens="0"
 NeedsClean="yes"
 AS="N"
><t
 p="190.33 364.68"
 BoundingBox="190.77 357.32 199.10 364.87"
 LabelJustification="Left"
 LabelAlignment="Left"
><s font="21" size="10" color="0" face="96">Cl</s></t></n><n
 id="20"
 p="245.90 361"
 Z="19"
 Element="35"
 NumHydrogens="0"
 NeedsClean="yes"
 AS="N"
><t
 p="242.57 364.59"
 BoundingBox="243.31 357.41 252.45 364.59"
 LabelJustification="Left"
 LabelAlignment="Left"
><s font="21" size="10" color="0" face="96">Br</s></t></n><b
 id="5"
 Z="4"
 B="2"
 E="4"
 BS="N"
/><b
 id="7"
 Z="6"
 B="4"
 E="6"
 BS="N"
/><b
 id="9"
 Z="8"
 B="6"
 E="8"
 BS="N"
/><b
 id="11"
 Z="10"
 B="8"
 E="10"
 BS="N"
/><b
 id="13"
 Z="12"
 B="10"
 E="12"
 BS="N"
/><b
 id="15"
 Z="14"
 B="12"
 E="14"
 BS="N"
/><b
 id="17"
 Z="16"
 B="4"
 E="16"
 Display="WedgeBegin"
 BS="N"
/><b
 id="19"
 Z="18"
 B="8"
 E="18"
 Display="WedgeBegin"
 BS="N"
/><b
 id="21"
 Z="20"
 B="12"
 E="20"
 Display="WedgeBegin"
 BS="N"
/></fragment></page></CDXML>
)";
    auto mols = MolsFromCDXML(cdxml1);
    mols[0]->clearConformers();
    CHECK(MolToCXSmiles(*mols[0]) ==
          "CO[C@H](C)C[C@H](Cl)C[C@H](C)Br |o1:5,o2:8,&1:2|");
  }
}

TEST_CASE("Github #7528 - read fragments in groups") {
  std::string cdxmlbase =
      std::string(getenv("RDBASE")) + "/Code/GraphMol/test_data/CDXML/";
  SECTION("case 1") {
    auto fname = cdxmlbase + "github7467-grouped-fragments.cdxml";
    CDXMLParserParams params;
    params.sanitize = false;
    auto mols = MolsFromCDXMLFile(fname, params);
    REQUIRE(mols.size() == 2);
  }
}

TEST_CASE("Github #7501 - dative bonds") {
  std::string cdxmlbase =
      std::string(getenv("RDBASE")) + "/Code/GraphMol/test_data/CDXML/";
  SECTION("case 1") {
    auto fname = cdxmlbase + "github7501-dative.cdxml";
    CDXMLParserParams params;
    auto mols = MolsFromCDXMLFile(fname, params);
    CHECK(MolToSmiles(*mols[0]) ==
          "C[CH2](C)->[Os]12<-[CH3]CC[NH]->1CC=[NH]->2");  // All datives to the
                                                           // Osmium
  }
}

#ifdef RDK_BUILD_CHEMDRAW_SUPPORT
struct format_check {
  std::string filename;
  bool stream, iscdx, cdxres, cdxmlres, autores;
};

TEST_CASE("CDX and Formats") {
  std::string cdxmlbase =
      std::string(getenv("RDBASE")) + "/Code/GraphMol/test_data/CDXML/";
  
  std::string cdxbase =
      std::string(getenv("RDBASE")) + "/Code/GraphMol/test_data/CDX/";

  SECTION("READ CDX") {
    auto cdxfname = cdxmlbase + "ring-stereo1.cdx";
    auto cdxmlfname = cdxmlbase + "ring-stereo1.cdxml";
    // should default to CDXFormat Auto
    auto mols1 = MolsFromCDXMLFile(cdxfname);
    auto mols2 = MolsFromCDXMLFile(cdxmlfname);
    CHECK(MolToSmiles(*mols1[0]) == MolToSmiles(*mols2[0]));
  }
  
  SECTION("Read CDX Files/Streams") {
    const std::vector<std::pair<std::string, std::string>> tests {
      {"structure_1.cdx", "C1CCOC1"},     
      {"structure_2.cdx", "C1=CCN=C1"},   
      {"structure_3.cdx", "CS(C)=O"},     
      {"structure_4.cdx", "CCO"},         
      {"structure_5.cdx", "c1cc[nH]c1"},  
      {"structure_6.cdx", "c1ccoc1"}
    };
    
    for(const auto &test :  tests) {
      auto fname = cdxbase + test.first;
      // Read the file
      auto m = MolsFromCDXMLFile(fname);
      CHECK(MolToSmiles(*m[0]) == test.second);

      // Read the CDX stream
      auto size = std::filesystem::file_size(fname);
      std::string content(size, '\0');
      std::ifstream in(fname, std::ios::binary);
      in.read(&content[0], size);
      
      auto m2 = MolsFromCDXML(content);
      CHECK(MolToSmiles(*m2[0]) == test.second);
    }
  }

  SECTION("READ CDX/CDXML Blocks") {
    auto cdxfname = cdxmlbase + "ring-stereo1.cdx";
    auto cdxmlfname = cdxmlbase + "ring-stereo1.cdxml";
    
    auto size = std::filesystem::file_size(cdxfname);
    std::string content(size, '\0');
    std::ifstream in(cdxfname, std::ios::binary);
    in.read(&content[0], size);
    auto mols1 = MolsFromCDXML(content);

    auto size2 = std::filesystem::file_size(cdxmlfname);
    std::string content2(size2, ' ');
    std::ifstream in2(cdxmlfname, std::ios::binary);
    in2.read(&content2[0], size2);
    auto mols2 = MolsFromCDXML(content2);
    CHECK(MolToSmiles(*mols1[0]) == MolToSmiles(*mols2[0]));
  }
  
  SECTION("Check Formats") {
    auto cdxfilename = cdxmlbase + "ring-stereo1.cdx";
    auto cdxmlfilename = cdxmlbase + "ring-stereo1.cdxml";
    std::vector<format_check> checks { {cdxfilename, true, true, true, false, false},
				       {cdxfilename, false, true, true, false, true},
				       {cdxmlfilename, true, false, false, true, true},
				       {cdxmlfilename, false, false, false, true, true},
    };
    std::vector<CDXMLFormat> formats = { CDXMLFormat::CDX, CDXMLFormat::CDXML, CDXMLFormat::Auto };

    for(auto &check : checks) {
      if(check.stream) {
      } else {
	for(auto format : formats) {
	  bool hasmols = false;
	  bool exception = false;
	  try {
	    auto mols = MolsFromCDXMLFile(check.filename, CDXMLParserParams(true, true, format));	    
	    hasmols = mols.size() > 0;
	    //	    std::cerr << check.filename << " not stream " << (unsigned)format << " hasmols: " << hasmols << std::endl;

	  } catch (...) {
	    exception = true;
	    //	    std::cerr << check.filename << " not stream " << (unsigned)format << " exception" << std::endl;
	  }

	  bool expected = false;
	  if(format == CDXMLFormat::CDX)
	    expected = check.cdxres;
	  else if(format == CDXMLFormat::CDXML)
	    expected = check.cdxmlres;
	  else
	    expected = check.autores;

	  if (exception)
	    CHECK(expected == false);
	  else
	    CHECK(expected == hasmols);
	}
      }
    }
  }
}
#endif
	 
