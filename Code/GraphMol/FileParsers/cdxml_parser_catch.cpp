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
    SECTION("RING CHIRALITY") {
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
    SECTION("REACTION") {
        auto cdx = R"(<?xml version="1.0" encoding="UTF-8" ?>
 <!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >
 <CDXML
  CreationProgram="ChemDraw 21.0.0.28"
  Name="Untitled Document-1"
  BoundingBox="62.52 167.32 472.83 352.47"
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
  id="395"
  BoundingBox="0 0 540 720"
  HeaderPosition="36"
  FooterPosition="36"
  PrintTrimMarks="yes"
  HeightPages="1"
  WidthPages="1"
 ><fragment
  id="303"
  BoundingBox="62.52 210.76 146.12 275.02"
  Z="61"
 ><n
  id="300"
  p="63.02 259.45"
  Z="58"
  AS="N"
 /><n
  id="302"
  p="89 274.45"
  Z="60"
  AS="N"
 /><n
  id="304"
  p="114.98 259.45"
  Z="62"
  AS="N"
 /><n
  id="306"
  p="114.98 229.45"
  Z="64"
  AS="N"
 /><n
  id="308"
  p="89 214.45"
  Z="66"
  AS="N"
 /><n
  id="310"
  p="63.02 229.45"
  Z="68"
  AS="N"
 /><n
  id="328"
  p="140.96 214.45"
  Z="86"
  Element="17"
  NumHydrogens="0"
  NeedsClean="yes"
  AS="N"
 ><t
  p="137.35 218.13"
  BoundingBox="137.79 210.76 146.12 218.32"
  LabelJustification="Left"
  LabelAlignment="Left"
 ><s font="21" size="10" color="0" face="96">Cl</s></t></n><b
  id="312"
  Z="70"
  B="300"
  E="302"
  Order="2"
  BS="N"
  BondCircularOrdering="317 0 0 313"
 /><b
  id="313"
  Z="71"
  B="302"
  E="304"
  BS="N"
 /><b
  id="314"
  Z="72"
  B="304"
  E="306"
  Order="2"
  BS="N"
  BondCircularOrdering="313 0 329 315"
 /><b
  id="315"
  Z="73"
  B="306"
  E="308"
  BS="N"
 /><b
  id="316"
  Z="74"
  B="308"
  E="310"
  Order="2"
  BS="N"
  BondCircularOrdering="315 0 0 317"
 /><b
  id="317"
  Z="75"
  B="310"
  E="300"
  BS="N"
 /><b
  id="329"
  Z="87"
  B="306"
  E="328"
  RxnParticipation="MakeOrBreak"
  BS="N"
 ><objecttag
  TagType="Unknown"
  Name="query"
 ><t
  p="112.90 217.21"
  BoundingBox="113.56 211.83 125.74 217.21"
  CaptionLineHeight="variable"
 ><s font="21" size="7.5" color="0">Rxn</s></t></objecttag></b><annotation
  Keyword="ComponentID"
  Content="1"
 /></fragment><fragment
  id="320"
  BoundingBox="210.65 173.26 299.46 312.52"
  Z="78"
 ><n
  id="319"
  p="263.12 206.95"
  Z="77"
  AS="N"
 /><n
  id="321"
  p="237.13 221.95"
  Z="79"
  Element="5"
  NumHydrogens="1"
  NeedsClean="yes"
  AS="N"
 ><t
  p="240.47 225.53"
  BoundingBox="227.36 218.36 240.07 225.53"
  LabelJustification="Right"
  Justification="Right"
  LabelAlignment="Right"
 ><s font="21" size="10" color="0" face="96">BH</s></t></n><n
  id="323"
  p="263.12 176.95"
  Z="81"
  Element="8"
  NumHydrogens="1"
  NeedsClean="yes"
  AS="N"
 ><t
  p="259.23 180.63"
  BoundingBox="259.62 173.26 273.48 180.84"
  LabelJustification="Left"
  LabelAlignment="Left"
 ><s font="21" size="10" color="0" face="96">OH</s></t></n><n
  id="325"
  p="289.10 221.95"
  Z="83"
  Element="8"
  NumHydrogens="1"
  NeedsClean="yes"
  AS="N"
 ><t
  p="285.21 225.63"
  BoundingBox="285.60 218.26 299.46 225.84"
  LabelJustification="Left"
  LabelAlignment="Left"
 ><s font="21" size="10" color="0" face="96">OH</s></t></n><n
  id="330"
  p="237.13 251.95"
  Z="88"
  AS="N"
 /><n
  id="332"
  p="211.15 266.95"
  Z="90"
  AS="N"
 /><n
  id="334"
  p="211.15 296.95"
  Z="92"
  AS="N"
 /><n
  id="336"
  p="237.13 311.95"
  Z="94"
  AS="N"
 /><n
  id="338"
  p="263.12 296.95"
  Z="96"
  AS="N"
 /><n
  id="340"
  p="263.12 266.95"
  Z="98"
  AS="N"
 /><b
  id="322"
  Z="80"
  B="319"
  E="321"
  BS="N"
 /><b
  id="324"
  Z="82"
  B="319"
  E="323"
  BS="N"
 /><b
  id="326"
  Z="84"
  B="319"
  E="325"
  BS="N"
 /><b
  id="331"
  Z="89"
  B="321"
  E="330"
  RxnParticipation="MakeOrBreak"
  BS="N"
 ><objecttag
  TagType="Unknown"
  Name="query"
 ><t
  p="218.79 241.34"
  BoundingBox="219.45 235.96 231.63 241.34"
  CaptionLineHeight="variable"
 ><s font="21" size="7.5" color="0">Rxn</s></t></objecttag></b><b
  id="342"
  Z="100"
  B="330"
  E="332"
  Order="2"
  BS="N"
  BondCircularOrdering="347 331 0 343"
 /><b
  id="343"
  Z="101"
  B="332"
  E="334"
  BS="N"
 /><b
  id="344"
  Z="102"
  B="334"
  E="336"
  Order="2"
  BS="N"
  BondCircularOrdering="343 0 0 345"
 /><b
  id="345"
  Z="103"
  B="336"
  E="338"
  BS="N"
 /><b
  id="346"
  Z="104"
  B="338"
  E="340"
  Order="2"
  BS="N"
  BondCircularOrdering="345 0 0 347"
 /><b
  id="347"
  Z="105"
  B="340"
  E="330"
  BS="N"
 /><annotation
  Keyword="ComponentID"
  Content="2"
 /></fragment><fragment
  id="369"
  BoundingBox="419.86 167.32 472.83 318.47"
  Z="127"
 ><n
  id="348"
  p="420.36 182.89"
  Z="106"
  AS="N"
 /><n
  id="350"
  p="420.36 212.89"
  Z="108"
  AS="N"
 /><n
  id="352"
  p="446.35 227.89"
  Z="110"
  AS="N"
 /><n
  id="354"
  p="472.33 212.89"
  Z="112"
  AS="N"
 /><n
  id="356"
  p="472.33 182.89"
  Z="114"
  AS="N"
 /><n
  id="358"
  p="446.35 167.89"
  Z="116"
  AS="N"
 /><n
  id="366"
  p="420.36 272.89"
  Z="124"
  AS="N"
 /><n
  id="368"
  p="420.36 302.89"
  Z="126"
  AS="N"
 /><n
  id="370"
  p="446.35 317.89"
  Z="128"
  AS="N"
 /><n
  id="372"
  p="472.33 302.89"
  Z="130"
  AS="N"
 /><n
  id="374"
  p="472.33 272.89"
  Z="132"
  AS="N"
 /><n
  id="376"
  p="446.35 257.89"
  Z="134"
  AS="N"
 /><b
  id="360"
  Z="118"
  B="348"
  E="350"
  Order="2"
  BS="N"
  BondCircularOrdering="365 0 0 361"
 /><b
  id="361"
  Z="119"
  B="350"
  E="352"
  BS="N"
 /><b
  id="362"
  Z="120"
  B="352"
  E="354"
  Order="2"
  BS="N"
  BondCircularOrdering="361 387 0 363"
 /><b
  id="363"
  Z="121"
  B="354"
  E="356"
  BS="N"
 /><b
  id="364"
  Z="122"
  B="356"
  E="358"
  Order="2"
  BS="N"
  BondCircularOrdering="363 0 0 365"
 /><b
  id="365"
  Z="123"
  B="358"
  E="348"
  BS="N"
 /><b
  id="378"
  Z="136"
  B="366"
  E="368"
  Order="2"
  BS="N"
  BondCircularOrdering="383 0 0 379"
 /><b
  id="379"
  Z="137"
  B="368"
  E="370"
  BS="N"
 /><b
  id="380"
  Z="138"
  B="370"
  E="372"
  Order="2"
  BS="N"
  BondCircularOrdering="379 0 0 381"
 /><b
  id="381"
  Z="139"
  B="372"
  E="374"
  BS="N"
 /><b
  id="382"
  Z="140"
  B="374"
  E="376"
  Order="2"
  BS="N"
  BondCircularOrdering="381 0 387 383"
 /><b
  id="383"
  Z="141"
  B="376"
  E="366"
  BS="N"
 /><b
  id="387"
  Z="145"
  B="352"
  E="376"
  RxnParticipation="MakeOrBreak"
  BS="N"
 ><objecttag
  TagType="Unknown"
  Name="query"
 ><t
  p="451.19 244.58"
  BoundingBox="451.85 239.20 464.03 244.58"
  CaptionLineHeight="variable"
 ><s font="21" size="7.5" color="0">Rxn</s></t></objecttag></b><annotation
  Keyword="ComponentID"
  Content="3"
 /></fragment><t
  id="388"
  p="99.48 306.57"
  BoundingBox="94.64 297.83 104.28 309.02"
  Z="146"
  InterpretChemically="no"
  CaptionJustification="Center"
  Justification="Center"
  LineHeight="auto"
 ><s font="21" size="12" color="0">(I)</s></t><t
  id="390"
  p="250.21 344.07"
  BoundingBox="243.72 335.33 256.68 346.52"
  Z="147"
  InterpretChemically="no"
  CaptionJustification="Center"
  Justification="Center"
  LineHeight="auto"
 ><s font="21" size="12" color="0">(II)</s></t><t
  id="392"
  p="448.45 350.02"
  BoundingBox="440.28 341.27 456.58 352.47"
  Z="148"
  InterpretChemically="no"
  CaptionJustification="Center"
  Justification="Center"
  LineHeight="auto"
 ><s font="21" size="12" color="0">(III)</s></t><graphic
  id="318"
  BoundingBox="178.39 242.89 178.39 250.39"
  Z="76"
  GraphicType="Symbol"
  SymbolType="Plus"
 /><graphic
  id="327"
  SupersededBy="396"
  BoundingBox="389.76 242.89 329.56 242.89"
  Z="85"
  GraphicType="Line"
  ArrowType="FullHead"
  HeadSize="2250"
 /><scheme
  id="397"
 ><step
  id="398"
  ReactionStepReactants="303 320"
  ReactionStepProducts="369"
  ReactionStepArrows="327"
  ReactionStepAtomMap="306 348 300 354 302 356 304 358 308 350 310 352 330 376 332 366 334 368 336 370 338 372 340 374"
  ReactionStepAtomMapManual="306 348"
  ReactionStepAtomMapAuto="300 354 302 356 304 358 308 350 310 352 330 376 332 366 334 368 336 370 338 372 340 374"
 /></scheme><step
  id="399"
 /><chemicalproperty
  id="389"
  ChemicalPropertyDisplayID="388"
  ChemicalPropertyIsActive="yes"
  ChemicalPropertyType="22"
  BasisObjects="300 302 303 304 306 308 310 312 313 314 315 316 317 328 329"
 /><chemicalproperty
  id="391"
  ChemicalPropertyDisplayID="390"
  ChemicalPropertyIsActive="yes"
  ChemicalPropertyType="22"
  BasisObjects="319 320 321 322 323 324 325 326 330 331 332 334 336 338 340 342 343 344 345 346 347"
 /><chemicalproperty
  id="393"
  ChemicalPropertyDisplayID="392"
  ChemicalPropertyIsActive="yes"
  ChemicalPropertyType="22"
  BasisObjects="348 350 352 354 356 358 360 361 362 363 364 365 366 368 369 370 372 374 376 378 379 380 381 382 383 387"
 /><arrow
  id="396"
  BoundingBox="329.56 236.27 389.76 248.52"
  Z="85"
  FillType="None"
  ArrowheadHead="Full"
  ArrowheadType="Solid"
  HeadSize="2250"
  ArrowheadCenterSize="1969"
  ArrowheadWidth="563"
  Head3D="389.76 242.89 0"
  Tail3D="329.56 242.89 0"
  Center3D="536.80 393.89 0"
  MajorAxisEnd3D="597 393.89 0"
  MinorAxisEnd3D="536.80 500.82 0"
 /></page></CDXML>)";
        std::vector<std::string> expected = {"Cl[c:1]1[cH:4][cH:3][cH:2][cH:6][cH:5]1",
            "OC(O)B[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1",
            "[cH:1]1[cH:4][cH:3][cH:2][c:6](-[c:7]2[cH:8][cH:9][cH:10][cH:11][cH:12]2)[cH:5]1"};
        std::stringstream iss(cdx);
        auto mols = CDXMLToMols(iss);
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
}
