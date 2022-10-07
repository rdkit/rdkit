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
#include "catch.hpp"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <RDGeneral/FileParseException.h>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BadFileException.h>

using namespace RDKit;

TEST_CASE("CDXML") {
    std::string cdxmlbase = std::string(getenv("RDBASE")) + "/Code/GraphMol/test_data/CDXML/";
    
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
        std::stringstream iss(cdxml1);
        auto mols = CDXMLDataStreamToMols(iss);
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == "CC(C)(C)OC(=O)C1CCCCCC1");
        }
    }
    SECTION("REACTION") {
        std::string fname =  cdxmlbase + "reaction-with-boc.cdxml";
        std::vector<std::string> expected = {"CC(C)(C)OC(=O)C1CCCCCC1[*:1]",
            "c1ccc([*:1])cc1",
            "C1CC1",
            "C1CCC1"};
        auto mols = CDXMLFileToMols(fname);
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    SECTION("RING CHIRALITY") {
        std::string fname =  cdxmlbase + "ring-stereo1.cdxml";
        std::vector<std::string> expected = {"C1CC[C@H]2CCCC[C@H]2C1"};
        auto mols = CDXMLFileToMols(fname);
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    SECTION("SIMPLE CHIRAL") {
        std::string fname = cdxmlbase + "chirality1.cdxml";
        std::vector<std::string> expected = {"C[C@H](N)C[C@H](C)N"};
        auto mols = CDXMLFileToMols(fname);
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    SECTION("CDXML-CISTRANS") {
        auto fname = cdxmlbase + "cistrans1.cdxml";
        std::vector<std::string> expected = {"F/C(I)=C(\\Cl)Br"};
        auto mols = CDXMLFileToMols(fname);
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    SECTION("REACTION") {
        auto fname = cdxmlbase + "rxn1.cdxml";
        std::vector<std::string> expected = {"Cl[c:1]1[cH:4][cH:3][cH:2][cH:6][cH:5]1",
            "OC(O)B[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1",
            "[cH:1]1[cH:4][cH:3][cH:2][c:6](-[c:7]2[cH:8][cH:9][cH:10][cH:11][cH:12]2)[cH:5]1"};
        auto mols = CDXMLFileToMols(fname);
        int i=0;
        for(auto &mol : mols) {
            CHECK(mol->getProp<unsigned int>("CDX_SCHEME_ID") == 397);
            CHECK(mol->getProp<unsigned int>("CDX_STEP_ID") == 398);
            if (i == 0) {
                CHECK(mol->getProp<unsigned int>("CDX_REAGENT_ID") == 0);
            }
            else if (i == 1) {
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
            std::vector<std::string> expected = {"[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]"};
            auto mols = CDXMLFileToMols(fname);
            CHECK(mols.size()==expected.size());
            int i=0;
            for(auto &mol : mols) {
                CHECK(MolToSmiles(*mol) == expected[i++]);
            }
        }
        {
            std::vector<std::string> expected = {"[2H]C1=C([2H])C([2H])=C([2H])C([2H])=C1[2H]"};
            auto mols = CDXMLFileToMols(fname, false, false);
            CHECK(mols.size()==expected.size());
            int i=0;
            for(auto &mol : mols) {
                CHECK(MolToSmiles(*mol) == expected[i++]);
            }
        }
        {
            std::vector<std::string> expected = {"[2H]C1=C([2H])C([2H])=C([2H])C([2H])=C1[2H]"};
            auto mols = CDXMLFileToMols(fname, false, true);
            CHECK(mols.size()==expected.size());
            int i=0;
            for(auto &mol : mols) {
                CHECK(MolToSmiles(*mol) == expected[i++]);
            }
        }
    }
    SECTION("Queries") {
        auto fname = cdxmlbase + "query-atoms.cdxml";
        
        std::vector<std::string> expected = {"*c1ccccc1", "*c1ccccc1", "*c1ccccc1"};
        std::vector<std::string> expected_smarts = {
            "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-*",
            "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[!#1]",
            "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[!#6&!#1]"};
        auto mols = CDXMLFileToMols(fname);
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmarts(*mol) == expected_smarts[i]);
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    SECTION("ElementList") {
        auto fname = cdxmlbase + "element-list.cdxml";
        
        std::vector<std::string> expected = {"[C]CC"};
        std::vector<std::string> expected_smarts = { "[#6]-[#6]-[#6,#7,#8,#16]" };
        auto mols = CDXMLFileToMols(fname);
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmarts(*mol) == expected_smarts[i]);
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    SECTION("Enhanced Stereo") {
        auto fname = cdxmlbase + "beta-cypermethrin.cdxml";
        std::vector<std::string> expected = {"CC1(C)[C@H](C=C(Cl)Cl)[C@H]1C(=O)O[C@@H](C#N)c1cccc(Oc2ccccc2)c1"};
        std::vector<std::string> expected_cx = {"CC1(C)[C@H](C=C(Cl)Cl)[C@H]1C(=O)O[C@@H](C#N)c1cccc(Oc2ccccc2)c1 |&1:3,&2:8,12|"};
        auto mols = CDXMLFileToMols(fname);
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            mol.get()->clearConformers();
            CHECK(MolToSmiles(*mol) == expected[i]);
            CHECK(MolToCXSmiles(*mol) == expected_cx[i++]);
        }
    }
    SECTION("Enhanced Stereo 2") {
        auto fname = cdxmlbase + "beta-cypermethrin-or-abs.cdxml";
        std::vector<std::string> expected = {"CC1(C)[C@H](C=C(Cl)Cl)[C@H]1C(=O)O[C@@H](C#N)c1cccc(Oc2ccccc2)c1"};
        std::vector<std::string> expected_cx = {"CC1(C)[C@H](C=C(Cl)Cl)[C@H]1C(=O)O[C@@H](C#N)c1cccc(Oc2ccccc2)c1 |o1:8,12|"};
        auto mols = CDXMLFileToMols(fname);
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            mol.get()->clearConformers();
            CHECK(MolToSmiles(*mol) == expected[i]);
            CHECK(MolToCXSmiles(*mol) == expected_cx[i++]);
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
          auto mols = CDXMLFileToMols(fname);
          CHECK(mols.size()==expected.size());
          int i=0;
          for(auto &mol : mols) {
              CHECK(MolToSmarts(*mol) == expected_smarts[i]);
              CHECK(MolToSmiles(*mol) == expected[i++]);
          }
        }
        // setting sanitization to false, we get both
        std::vector<std::string> expected = {
            "*C1=C([H])C([H])=C([H])C([H])=C1[H]",
            "*C1=C([H])N([H])=C([H])C([H])=C1[H]"};
        std::vector<std::string> expected_smarts = {
            "[#6]1(=[#6](-[#6](=[#6](-[#6](=[#6]-1-*)-[H])-[H])-[H])-[H])-[H]",
            "[#6]1(=[#6](-[#6](=[#7](-[#6](=[#6]-1-[!#1])-[H])-[H])-[H])-[H])-[H]",
            };
        auto mols = CDXMLFileToMols(fname, false);
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmarts(*mol) == expected_smarts[i]);
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
     
    SECTION("Fusion with chirality") {
        auto fname = cdxmlbase + "fusion-chiral.cdxml";
        auto mols = CDXMLFileToMols(fname);
        std::vector<std::string> expected = {"C[C@@H](O)[C@@H](C)O"};
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    
    SECTION("ChemDraw Template from the synthesis-workshop") {
        // this was hella fun to validate the stereo-chemistry...
        auto fname = cdxmlbase + "chemdraw_template1.cdxml";
        auto mols = CDXMLFileToMols(fname);
        std::vector<std::string> expected = {
            "CCC/C=C/C=C/C(=O)O[C@H]1/C(=C/C(=O)OC)C[C@H]2C[C@H]([C@@H](C)O)OC(=O)C[C@H](O)C[C@@H]3C[C@H](OC(C)=O)C(C)(C)[C@](O)(C[C@@H]4C/C(=C/C(=O)OC)C[C@H](/C=C/C(C)(C)[C@]1(O)O2)O4)O3",
            "[B]",
            "*",
            "[C]",
            "Cc1ccc2n1[C@@H]1[C@@H]3O[C@]([C@H](C)O)(C=C2)[C@H]1c1ccc(C)n1[C@@H]3C",
            "Cc1ccc2n1[C@H](C)C(=O)[C@@H]1[C@H]2C(=O)C=Cc2ccc(C)n21",
            "Cc1ccc2ccc(=O)ccn12",
            "Cc1cccn1[C@H](C)C=O",
            "Cc1ccc2ccc([O-])cc[n+]1-2",
            "Cc1ccc2ccc(=O)ccn12",
            "Cc1cccn1[C@H](C)C(C#N)O[Si](C)(C)C",
            "CC1CCC2(O)C3(OC4(O)C[C@]2(C)C2(O)[C@H](OC(=O)c5ccc[nH]5)C(O)(C(C)C)C4(C)C32O)C1O",
            "C=C(C)[C@H]1CC(=O)CC2=C(C1)[C@H]1C(=O)O[C@H]3C[C@@](C)(O)[C@@H](C2=O)[C@@H]13"};
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
     
    SECTION("deuterium atom") {
        auto fname = cdxmlbase + "deuterium-atom.cdxml";
        auto mols = CDXMLFileToMols(fname, false, false);
        std::vector<std::string> expected = {
            "[2H]"};
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    SECTION("ChemDraw Template 2 from the synthesis-workshop") {
        // this was another hella fun to validate the stereo-chemistry...
        //   there were so many stereo warnings in chemdraw, I'm just going to assume
        //    the rdkit is correct here...
        auto fname = cdxmlbase + "chemdraw_template2.cdxml";
        auto mols = CDXMLFileToMols(fname);
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
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    SECTION("ChemDraw Template 3 from the synthesis-workshop") {
        // this was another hella fun to validate the stereo-chemistry...
        //   there were so many stereo warnings in chemdraw, I'm just going to assume
        //    the rdkit is correct here...
        auto fname = cdxmlbase + "chemdraw_template3.cdxml";
        auto mols = CDXMLFileToMols(fname);
        std::vector<std::string> expected = {"CCC/C=C/C=C/C(=O)O[C@H]1/C(=C/C(=O)OC)C[C@H]2C[C@H]([C@@H](C)O)OC(=O)C[C@H](O)C[C@@H]3C[C@H](OC(C)=O)C(C)(C)[C@](O)(C[C@@H]4C/C(=C/C(=O)OC)C[C@H](/C=C/C(C)(C)[C@]1(O)O2)O4)O3",
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
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    SECTION("protecting group") {
        auto fname = cdxmlbase + "protecting-groups.cdxml";
        auto mols = CDXMLFileToMols(fname);
        std::vector<std::string> expected = {"CC[Si](C)(CC)CC", "CC"};
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    SECTION("protecting group 2") {
        auto fname = cdxmlbase + "protecting-groups2.cdxml";
        auto mols = CDXMLFileToMols(fname);
        std::vector<std::string> expected = {"CC[Si](C)(CC)CC", "CC"};
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    
    SECTION("floating protecting group") {
        auto fname = cdxmlbase + "floating-protecting-group.cdxml";
        // This is a weird one, chemdraw simply ignores the error that causes the bond issue, we should probably drop the floating fragment here if
        //  we are sanitizing
        auto mols = CDXMLFileToMols(fname);
        std::vector<std::string> expected = { "*", "C=C(C[C@H]([O])C[C@]1(O)O[C@H](C[C@@H](O)CC(=O)O)C[C@H](OC(C)=O)C1(C)C)C[Si](C)(C)C", "*.CC[Si](CC)CC"};
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    
    SECTION("Missing File Name") {
        try {
            auto mols = CDXMLFileToMols("missing file");
            CHECK(0); // Bad file exception not caught
        }
        catch (RDKit::BadFileException &) {
        
        }
    }
    
    SECTION("Aromatic ring (bondorder==4") {
        auto fname = cdxmlbase + "aromatic.cdxml";
        auto mols = CDXMLFileToMols(fname);
        std::vector<std::string> expected = { "c1ccccc1" };
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    SECTION("Malformed") {
        auto fname = cdxmlbase + "malformed.cdxml";
        try {
            auto mols = CDXMLFileToMols(fname);
            CHECK(0);
        } catch (FileParseException &e) {
            CHECK(std::string(e.what()) == "expected > at line: 373");
        }
    }
    SECTION("Lots of stereo") {
        auto fname = cdxmlbase + "stereo.cdxml";
        std::vector<std::string> expected = { "C[C@@H](Cl)[C@H](N)O.C[C@@H](F)[C@H](N)O.C[C@H](Br)[C@@H](N)O.C[C@H](I)[C@@H](N)O"};
        auto mols = CDXMLFileToMols(fname);
        CHECK(mols.size()==expected.size());
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    SECTION("Lots of stereo") {
        {
            auto fname = cdxmlbase + "bad-id.cdxml";
            auto mols = CDXMLFileToMols(fname);
            CHECK(mols.size() == 0);
        }
        {
            auto fname = cdxmlbase + "bad-coords.cdxml";
            auto mols = CDXMLFileToMols(fname);
            CHECK(mols.size() == 0);
        }
        {
            auto fname = cdxmlbase + "bad-bondorder.cdxml";
            auto mols = CDXMLFileToMols(fname);
            CHECK(mols.size() == 0);
        }
        {
            auto fname = cdxmlbase + "bad-bondorder2.cdxml";
            auto mols = CDXMLFileToMols(fname);
            CHECK(mols.size() == 0);
        }
    }
}
