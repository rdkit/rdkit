//
//  Copyright (C) 2020-2021 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
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
#include <RDGeneral/FileParseException.h>
#include <boost/algorithm/string.hpp>

using namespace RDKit;

TEST_CASE("CDXML") {
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
        auto mols = CDXMLToMols(iss);
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == "CC(C)(C)OC(=O)C1CCCCCC1");
        }
    }
    SECTION("REACTION") {
        std::string cdxml1 = R"(<?xml version="1.0" encoding="UTF-8" ?>
 <!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >
 <CDXML
  CreationProgram="ChemDraw JS 2.0.0.9"
  Name="ACS Document 1996"
  BoundingBox="94.75 160.95 310.72 223.50"
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
  id="76"
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
  BoundingBox="94.75 165.41 154.89 211.17"
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
  id="77"
 ><n
  id="78"
  p="148.17 207.09"
  Element="8"
  NumHydrogens="0"
 /><n
  id="79"
  p="162.52 207.09"
 /><n
  id="80"
  p="176.87 207.09"
 /><n
  id="81"
  p="169.69 194.67"
 /><n
  id="82"
  p="169.69 219.52"
 /><n
  id="83"
  p="140.99 194.67"
 /><n
  id="84"
  p="148.17 182.24"
  Element="8"
  NumHydrogens="0"
 /><n
  id="85"
  p="126.64 194.67"
  NodeType="ExternalConnectionPoint"
 /><b
  id="86"
  B="83"
  E="84"
  Order="2"
 /><b
  id="87"
  B="79"
  E="82"
 /><b
  id="88"
  B="79"
  E="80"
 /><b
  id="89"
  B="79"
  E="81"
 /><b
  id="90"
  B="78"
  E="79"
 /><b
  id="91"
  B="78"
  E="83"
 /><b
  id="92"
  B="85"
  E="83"
 /></fragment><t
  p="137.66 198.28"
  BoundingBox="137.66 189.64 154.89 198.28"
  LabelJustification="Left"
  LabelAlignment="Left"
 ><s font="24" size="9.95" color="0" face="96">Boc</s></t></n><n
  id="50"
  p="129.32 170.43"
  Z="44"
  NodeType="GenericNickname"
  NumHydrogens="0"
  NeedsClean="yes"
  AS="N"
  GenericNickname="R"
 ><t
  p="125.71 174.05"
  BoundingBox="125.71 165.41 137.11 176.10"
  LabelJustification="Left"
  LabelAlignment="Left"
 ><s font="24" size="9.95" color="0" face="96">R1</s></t></n><b
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
 /><b
  id="51"
  Z="45"
  B="17"
  E="50"
  BS="N"
 /></fragment><fragment
  id="35"
  BoundingBox="265.23 175.24 310.72 209.41"
  Z="29"
 ><n
  id="32"
  p="265.53 187.47"
  Z="26"
  AS="N"
 /><n
  id="34"
  p="265.53 201.87"
  Z="28"
  AS="N"
 /><n
  id="36"
  p="278 209.07"
  Z="30"
  AS="N"
 /><n
  id="38"
  p="290.47 201.87"
  Z="32"
  AS="N"
 /><n
  id="40"
  p="290.47 187.47"
  Z="34"
  AS="N"
 /><n
  id="42"
  p="278 180.27"
  Z="36"
  AS="N"
 /><n
  id="52"
  p="302.94 180.27"
  Z="46"
  NodeType="GenericNickname"
  NumHydrogens="0"
  NeedsClean="yes"
  AS="N"
  GenericNickname="R"
 ><t
  p="299.33 183.88"
  BoundingBox="299.33 175.24 310.72 185.93"
  LabelJustification="Left"
  LabelAlignment="Left"
 ><s font="24" size="9.95" color="0" face="96">R1</s></t></n><b
  id="44"
  Z="38"
  B="32"
  E="34"
  BS="N"
 /><b
  id="45"
  Z="39"
  B="34"
  E="36"
  Order="2"
  BS="N"
  BondCircularOrdering="44 0 0 46"
 /><b
  id="46"
  Z="40"
  B="36"
  E="38"
  BS="N"
 /><b
  id="47"
  Z="41"
  B="38"
  E="40"
  Order="2"
  BS="N"
  BondCircularOrdering="46 0 53 48"
 /><b
  id="48"
  Z="42"
  B="40"
  E="42"
  BS="N"
 /><b
  id="49"
  Z="43"
  B="42"
  E="32"
  Order="2"
  BS="N"
  BondCircularOrdering="48 0 0 44"
 /><b
  id="53"
  Z="47"
  B="40"
  E="52"
  BS="N"
 /></fragment><fragment
  id="57"
  BoundingBox="206.21 160.95 219.58 176.39"
  Z="51"
 ><n
  id="54"
  p="206.51 161.47"
  Z="48"
  AS="N"
 /><n
  id="56"
  p="206.51 175.87"
  Z="50"
  AS="N"
 /><n
  id="58"
  p="218.98 168.67"
  Z="52"
  AS="N"
 /><b
  id="60"
  Z="54"
  B="54"
  E="56"
  BS="N"
 /><b
  id="61"
  Z="55"
  B="56"
  E="58"
  BS="N"
 /><b
  id="62"
  Z="56"
  B="58"
  E="54"
  BS="N"
 /></fragment><fragment
  id="66"
  BoundingBox="205.83 208.50 220.83 223.50"
  Z="60"
 ><n
  id="63"
  p="206.13 208.80"
  Z="57"
  AS="N"
 /><n
  id="65"
  p="206.13 223.20"
  Z="59"
  AS="N"
 /><n
  id="67"
  p="220.53 223.20"
  Z="61"
  AS="N"
 /><n
  id="69"
  p="220.53 208.80"
  Z="63"
  AS="N"
 /><b
  id="71"
  Z="65"
  B="63"
  E="65"
  BS="N"
 /><b
  id="72"
  Z="66"
  B="65"
  E="67"
  BS="N"
 /><b
  id="73"
  Z="67"
  B="67"
  E="69"
  BS="N"
 /><b
  id="74"
  Z="68"
  B="69"
  E="63"
  BS="N"
 /></fragment><graphic
  id="31"
  SupersededBy="94"
  BoundingBox="242 198 189.33 198"
  Z="25"
  GraphicType="Line"
  ArrowType="FullHead"
  HeadSize="2250"
 /><scheme
  id="95"
 ><step
  id="96"
  ReactionStepReactants="10"
  ReactionStepProducts="35"
  ReactionStepArrows="31"
  ReactionStepObjectsAboveArrow="57"
  ReactionStepObjectsBelowArrow="66"
 /></scheme><arrow
  id="94"
  BoundingBox="189.33 193.95 242 201.38"
  Z="25"
  FillType="None"
  ArrowheadHead="Full"
  ArrowheadType="Solid"
  HeadSize="2250"
  ArrowheadCenterSize="1969"
  ArrowheadWidth="563"
  Head3D="242 198 0"
  Tail3D="189.33 198 0"
  Center3D="336.67 297 0"
  MajorAxisEnd3D="389.33 297 0"
  MinorAxisEnd3D="336.67 349.67 0"
 /></page></CDXML>)";
        std::vector<std::string> expected = {"CC(C)(C)OC(=O)C1CCCCCC1[*:1]",
            "c1ccc([*:1])cc1",
            "C1CC1",
            "C1CCC1"};
        std::stringstream iss(cdxml1);
        auto mols = CDXMLToMols(iss);
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    SECTION("RINC CHIRALITY") {
        std::string cdx = R"(<?xml version="1.0" encoding="UTF-8" ?>
        <!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >
        <CDXML
         CreationProgram="ChemDraw 21.0.0.28"
         Name="test1.cdxml"
         BoundingBox="217.54 328.75 322.46 391.25"
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
         CaptionFont="21"
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
         MacPrintInfo="00030000004800480000000002DE0240FFEEFFEE030602520367052803FC00020000004800480000000002DE0240000100000064000000010001010100000001270F000100010000000000000000000000000002001901900000000000400000000000000000000100000000000000000000000000000000"
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
        <font id="21" charset="x-mac-roman" name="Helvetica"/>
        </fonttable><page
         id="108"
         BoundingBox="0 0 540 720"
         HeaderPosition="36"
         FooterPosition="36"
         PrintTrimMarks="yes"
         HeightPages="1"
         WidthPages="1"
        ><fragment
         id="24"
         BoundingBox="217.54 328.75 322.46 391.25"
         Z="16"
        ><n
         id="1"
         p="321.96 345"
         Z="17"
         AS="N"
        /><n
         id="2"
         p="295.98 330"
         Z="18"
         AS="N"
        /><n
         id="3"
         p="270 345"
         Z="19"
         Geometry="Tetrahedral"
         AS="s"
         BondOrdering="13 14 19 0"
        /><n
         id="4"
         p="244.02 330"
         Z="20"
         AS="N"
        /><n
         id="5"
         p="218.04 345"
         Z="21"
         AS="N"
        /><n
         id="6"
         p="218.04 375"
         Z="22"
         AS="N"
        /><n
         id="7"
         p="244.02 390"
         Z="23"
         AS="N"
        /><n
         id="8"
         p="270 375"
         Z="24"
         Geometry="Tetrahedral"
         AS="s"
         BondOrdering="19 18 20 0"
        /><n
         id="9"
         p="295.98 390"
         Z="25"
         AS="N"
        /><n
         id="10"
         p="321.96 375"
         Z="26"
         AS="N"
        /><b
         id="12"
         Z="27"
         B="1"
         E="2"
         BS="N"
        /><b
         id="13"
         Z="28"
         B="2"
         E="3"
         Display="WedgedHashEnd"
         BS="N"
        /><b
         id="14"
         Z="29"
         B="3"
         E="4"
         BS="N"
        /><b
         id="15"
         Z="30"
         B="4"
         E="5"
         BS="N"
        /><b
         id="16"
         Z="31"
         B="5"
         E="6"
         BS="N"
        /><b
         id="17"
         Z="32"
         B="6"
         E="7"
         BS="N"
        /><b
         id="18"
         Z="33"
         B="7"
         E="8"
         Display="WedgedHashEnd"
         BS="N"
        /><b
         id="19"
         Z="34"
         B="3"
         E="8"
         BS="N"
        /><b
         id="20"
         Z="35"
         B="8"
         E="9"
         BS="N"
        /><b
         id="21"
         Z="36"
         B="9"
         E="10"
         BS="N"
        /><b
         id="22"
         Z="37"
         B="1"
         E="10"
         BS="N"
        /></fragment></page></CDXML>)";
        std::vector<std::string> expected = {"C1CC[C@H]2CCCC[C@H]2C1"};
        std::stringstream iss(cdx);
        auto mols = CDXMLToMols(iss);
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    SECTION("SIMPLE CHIRAL") {
        std::string cdx = R"(<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >
<CDXML
 CreationProgram="ChemDraw 21.0.0.28"
 Name="test1.cdxml"
 BoundingBox="217.79 329.73 322.21 393.08"
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
 CaptionFont="21"
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
 MacPrintInfo="00030000004800480000000002DE0240FFEEFFEE030602520367052803FC00020000004800480000000002DE0240000100000064000000010001010100000001270F000100010000000000000000000000000002001901900000000000400000000000000000000100000000000000000000000000000000"
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
<font id="21" charset="x-mac-roman" name="Helvetica"/>
</fonttable><page
 id="280"
 BoundingBox="0 0 540 720"
 HeaderPosition="36"
 FooterPosition="36"
 PrintTrimMarks="yes"
 HeightPages="1"
 WidthPages="1"
><fragment
 id="255"
 BoundingBox="217.79 329.73 322.21 393.08"
 Z="193"
><n
 id="256"
 p="295.98 333.32"
 Z="194"
 Element="7"
 NumHydrogens="2"
 NeedsClean="yes"
 AS="N"
><t
 p="292.37 336.90"
 BoundingBox="293.13 329.73 310.67 339.20"
 LabelJustification="Left"
 LabelAlignment="Left"
><s font="21" size="10" color="0" face="96">NH2</s></t></n><n
 id="257"
 p="295.98 363.32"
 Z="195"
 Geometry="Tetrahedral"
 AS="S"
 BondOrdering="263 264 265 0"
/><n
 id="258"
 p="321.96 378.32"
 Z="196"
 AS="N"
/><n
 id="259"
 p="270 378.32"
 Z="197"
 NumHydrogens="2"
 NeedsClean="yes"
 AS="N"
><t
 p="266.39 382"
 BoundingBox="266.83 374.63 277.47 393.08"
 LabelJustification="Left"
 LabelAlignment="Below"
 LineStarts="2 4"
><s font="21" size="10" color="0" face="96">CH2</s></t></n><n
 id="260"
 p="244.02 363.32"
 Z="198"
 NumHydrogens="1"
 NeedsClean="yes"
 Geometry="Tetrahedral"
 AS="S"
 BondOrdering="266 267 0 268"
><t
 p="240.41 367"
 BoundingBox="240.85 359.63 247.21 375.78"
 LabelJustification="Left"
 LabelAlignment="Below"
 LineStarts="2 3"
><s font="21" size="10" color="0" face="96">CH</s></t></n><n
 id="261"
 p="218.04 378.32"
 Z="199"
 AS="N"
/><n
 id="262"
 p="244.02 333.32"
 Z="200"
 Element="7"
 NumHydrogens="2"
 NeedsClean="yes"
 AS="N"
><t
 p="240.41 336.90"
 BoundingBox="241.17 329.73 258.71 339.20"
 LabelJustification="Left"
 LabelAlignment="Left"
><s font="21" size="10" color="0" face="96">NH2</s></t></n><b
 id="263"
 Z="201"
 B="256"
 E="257"
 Display="WedgeEnd"
 BS="N"
/><b
 id="264"
 Z="202"
 B="257"
 E="258"
 BS="N"
/><b
 id="265"
 Z="203"
 B="257"
 E="259"
 BS="N"
/><b
 id="266"
 Z="204"
 B="259"
 E="260"
 BS="N"
/><b
 id="267"
 Z="205"
 B="260"
 E="261"
 BS="N"
/><b
 id="268"
 Z="206"
 B="260"
 E="262"
 Display="WedgedHashBegin"
 BS="N"
/></fragment></page></CDXML>)";
        std::vector<std::string> expected = {"C[C@H](N)C[C@H](C)N"};
        std::stringstream iss(cdx);
        auto mols = CDXMLToMols(iss);
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
    SECTION("CDXML-CISTRANS") {
        auto str = R"(<?xml version="1.0" encoding="UTF-8" ?>
 <!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >
 <CDXML
  CreationProgram="ChemDraw 21.0.0.28"
  Name="Untitled Document-1"
  BoundingBox="228.25 318.91 311.75 401.09"
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
  CaptionFont="21"
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
  MacPrintInfo="0003000000480048000000000300024CFFF4FFF4030C02580367052803FC0002000000480048000000000300024C000100000064000000010001010100000001270F000100010000000000000000000000000002001901900000000000400000000000000000000100000000000000000000000000000000"
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
 <font id="21" charset="x-mac-roman" name="Helvetica"/>
 </fonttable><page
  id="296"
  BoundingBox="0 0 540 720"
  HeaderPosition="36"
  FooterPosition="36"
  PrintTrimMarks="yes"
  HeightPages="1"
  WidthPages="1"
 ><fragment
  id="14"
  BoundingBox="228.25 318.91 311.75 401.09"
  Z="41"
 ><n
  id="1"
  p="254.63 397.50"
  Z="42"
  Element="9"
  NumHydrogens="0"
  NeedsClean="yes"
  AS="N"
 ><t
  p="251.58 401.09"
  BoundingBox="252.43 393.91 257.41 401.09"
  LabelJustification="Left"
 ><s font="21" size="10" color="0" face="96">F</s></t></n><n
  id="2"
  p="254.63 367.50"
  Z="43"
  AS="N"
 /><n
  id="3"
  p="228.65 352.50"
  Z="44"
  Element="53"
  NumHydrogens="0"
  NeedsClean="yes"
  AS="N"
 ><t
  p="227.26 356.09"
  BoundingBox="228.25 348.91 229.23 356.09"
  LabelJustification="Left"
 ><s font="21" size="10" color="0" face="96">I</s></t></n><n
  id="4"
  p="280.61 352.50"
  Z="45"
  AS="N"
 /><n
  id="5"
  p="306.60 367.50"
  Z="46"
  Element="17"
  NumHydrogens="0"
  NeedsClean="yes"
  AS="N"
 ><t
  p="302.98 371.18"
  BoundingBox="303.42 363.82 311.75 371.37"
  LabelJustification="Left"
  LabelAlignment="Left"
 ><s font="21" size="10" color="0" face="96">Cl</s></t></n><n
  id="6"
  p="280.61 322.50"
  Z="47"
  Element="35"
  NumHydrogens="0"
  NeedsClean="yes"
  AS="N"
 ><t
  p="277.28 326.09"
  BoundingBox="278.02 318.91 287.16 326.09"
  LabelJustification="Left"
  LabelAlignment="Left"
 ><s font="21" size="10" color="0" face="96">Br</s></t></n><b
  id="8"
  Z="48"
  B="1"
  E="2"
  BS="N"
 /><b
  id="9"
  Z="49"
  B="2"
  E="3"
  BS="N"
 /><b
  id="10"
  Z="50"
  B="2"
  E="4"
  Order="2"
  BS="Z"
  BondCircularOrdering="9 8 11 12"
 /><b
  id="11"
  Z="51"
  B="4"
  E="5"
  BS="N"
 /><b
  id="12"
  Z="52"
  B="4"
  E="6"
  BS="N"
 /></fragment></page></CDXML>)";
        std::vector<std::string> expected = {"F/C(I)=C(\\Cl)Br"};
        std::stringstream iss(str);
        auto mols = CDXMLToMols(iss);
        int i=0;
        for(auto &mol : mols) {
            CHECK(MolToSmiles(*mol) == expected[i++]);
        }
    }
}
