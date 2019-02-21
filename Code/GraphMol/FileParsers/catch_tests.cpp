//
//  Copyright (C) 2018 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "RDGeneral/test.h"
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

using namespace RDKit;

TEST_CASE("Basic SVG Parsing", "[SVG,parser]") {
  SECTION("basics") {
    std::string svg = R"SVG(<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='200px' height='200px' >
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200' height='200' x='0' y='0'> </rect>
<path d='M 9.09091,89.4974 24.2916,84.7462' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 24.2916,84.7462 39.4923,79.9949' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 86.2908,106.814 75.1709,93.4683 72.0765,96.8285 86.2908,106.814' style='fill:#000000;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 75.1709,93.4683 57.8622,86.8431 64.051,80.1229 75.1709,93.4683' style='fill:#0000FF;fill-rule:evenodd;stroke:#0000FF;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 75.1709,93.4683 72.0765,96.8285 57.8622,86.8431 75.1709,93.4683' style='fill:#0000FF;fill-rule:evenodd;stroke:#0000FF;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 86.2908,106.814 82.1459,125.293' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 82.1459,125.293 78.0009,143.772' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 86.2908,106.814 129.89,93.1862' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 134.347,94.186 138.492,75.7069' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 138.492,75.7069 142.637,57.2277' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 125.432,92.1865 129.577,73.7074' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 129.577,73.7074 133.722,55.2282' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 129.89,93.1862 142.557,104.852' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 142.557,104.852 155.224,116.517' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<text x='39.4923' y='83.483' style='font-size:15px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#0000FF' ><tspan>NH</tspan></text>
<text x='67.6656' y='158.998' style='font-size:15px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#00CC00' ><tspan>Cl</tspan></text>
<text x='132.777' y='56.228' style='font-size:15px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#FF0000' ><tspan>O</tspan></text>
<text x='149.782' y='131.743' style='font-size:15px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#FF0000' ><tspan>OH</tspan></text>
<text x='89.9952' y='194' style='font-size:12px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' ><tspan>m1</tspan></text>
<metadata>
<rdkit:mol xmlns:rdkit = "http://www.rdkit.org/xml" version="0.9">
<rdkit:atom idx="1" atom-smiles="[CH3]" drawing-x="9.09091" drawing-y="89.4974" x="-2.78651" y="0.295614" z="0" />
<rdkit:atom idx="2" atom-smiles="[NH]" drawing-x="52.6897" drawing-y="75.8699" x="-1.35482" y="0.743114" z="0" />
<rdkit:atom idx="3" atom-smiles="[C@H]" drawing-x="86.2908" drawing-y="106.814" x="-0.251428" y="-0.273019" z="0" />
<rdkit:atom idx="4" atom-smiles="[Cl]" drawing-x="76.2932" drawing-y="151.385" x="-0.579728" y="-1.73665" z="0" />
<rdkit:atom idx="5" atom-smiles="[C]" drawing-x="129.89" drawing-y="93.1862" x="1.18027" y="0.174481" z="0" />
<rdkit:atom idx="6" atom-smiles="[O]" drawing-x="139.887" drawing-y="48.6148" x="1.50857" y="1.63811" z="0" />
<rdkit:atom idx="7" atom-smiles="[OH]" drawing-x="163.491" drawing-y="124.13" x="2.28366" y="-0.841652" z="0" />
<rdkit:bond idx="1" begin-atom-idx="1" end-atom-idx="2" bond-smiles="-" />
<rdkit:bond idx="2" begin-atom-idx="2" end-atom-idx="3" bond-smiles="-" />
<rdkit:bond idx="3" begin-atom-idx="3" end-atom-idx="4" bond-smiles="-" />
<rdkit:bond idx="4" begin-atom-idx="3" end-atom-idx="5" bond-smiles="-" />
<rdkit:bond idx="5" begin-atom-idx="5" end-atom-idx="6" bond-smiles="=" />
<rdkit:bond idx="6" begin-atom-idx="5" end-atom-idx="7" bond-smiles="-" />
</rdkit:mol></metadata>
</svg>)SVG";

    std::unique_ptr<RWMol> mol(RDKitSVGToMol(svg));
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 7);
    CHECK(mol->getNumConformers() == 1);
    CHECK_FALSE(mol->getConformer().is3D());
    auto smiles = MolToSmiles(*mol);
    CHECK(smiles == "CN[C@H](Cl)C(=O)O");
  }
}

TEST_CASE(
    "Github #2040: Failure to parse V3K mol file with bonds to multi-center "
    "linkage points",
    "[bug,parser]") {
  std::string rdbase = getenv("RDBASE");
  SECTION("basics") {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/github2040_1.mol";
    std::unique_ptr<RWMol> mol(
        MolFileToMol(fName, false));  // don't sanitize yet
    REQUIRE(mol);
    CHECK(mol->getBondWithIdx(0)->getBondType() == Bond::SINGLE);
    CHECK(
        mol->getBondWithIdx(0)->hasProp(common_properties::_MolFileBondEndPts));
    CHECK(mol->getBondWithIdx(0)->getProp<std::string>(
              common_properties::_MolFileBondEndPts) == "(3 5 4 3)");
    CHECK(
        mol->getBondWithIdx(0)->hasProp(common_properties::_MolFileBondAttach));
    CHECK(mol->getBondWithIdx(0)->getProp<std::string>(
              common_properties::_MolFileBondAttach) == "ANY");
    CHECK(mol->getBondWithIdx(1)->getBondType() == Bond::AROMATIC);
  }
}

TEST_CASE("Github #2225: failure round-tripping mol block with Q atoms",
          "[bug,writer]") {
  std::string rdbase = getenv("RDBASE");
  SECTION("basics") {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/github2225_1.mol";
    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 7);
    REQUIRE(!mol->getAtomWithIdx(0)->hasQuery());
    REQUIRE(mol->getAtomWithIdx(6)->hasQuery());
    auto outBlock = MolToMolBlock(*mol);
    REQUIRE(outBlock.find(" Q ") != std::string::npos);
    REQUIRE(outBlock.find(" ALS ") == std::string::npos);
    std::unique_ptr<RWMol> mol2(MolBlockToMol(outBlock));
    REQUIRE(mol2);
    REQUIRE(mol2->getNumAtoms() == 7);
    REQUIRE(!mol2->getAtomWithIdx(0)->hasQuery());
    REQUIRE(mol2->getAtomWithIdx(6)->hasQuery());
    auto outBlock2 = MolToMolBlock(*mol2);
    REQUIRE(outBlock2.find(" Q ") != std::string::npos);
    REQUIRE(outBlock2.find(" ALS ") == std::string::npos);
  }
  SECTION("check that SMARTS still works") {
    std::unique_ptr<RWMol> mol(SmartsToMol("C[#8,#7]"));
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 2);
    auto outBlock = MolToMolBlock(*mol);
    REQUIRE(outBlock.find(" Q ") == std::string::npos);
    REQUIRE(outBlock.find(" ALS ") != std::string::npos);
    std::unique_ptr<RWMol> mol2(MolBlockToMol(outBlock));
    REQUIRE(mol2);
    auto smarts = MolToSmarts(*mol2);
    REQUIRE(smarts == "[#6][#8,#7]");
  }
  SECTION("basics with v3K") {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/github2225_2.mol";
    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 7);
    REQUIRE(!mol->getAtomWithIdx(0)->hasQuery());
    REQUIRE(mol->getAtomWithIdx(6)->hasQuery());

    bool includeStereo = true;
    int confId = -1;
    bool kekulize = true;
    bool forceV3000 = true;
    auto outBlock =
        MolToMolBlock(*mol, includeStereo, confId, kekulize, forceV3000);
    REQUIRE(outBlock.find(" Q ") != std::string::npos);
    REQUIRE(outBlock.find(" ALS ") == std::string::npos);
    std::unique_ptr<RWMol> mol2(MolBlockToMol(outBlock));
    REQUIRE(mol2);
    REQUIRE(mol2->getNumAtoms() == 7);
    REQUIRE(!mol2->getAtomWithIdx(0)->hasQuery());
    REQUIRE(mol2->getAtomWithIdx(6)->hasQuery());
    auto outBlock2 =
        MolToMolBlock(*mol2, includeStereo, confId, kekulize, forceV3000);
    REQUIRE(outBlock2.find(" Q ") != std::string::npos);
    REQUIRE(outBlock2.find(" ALS ") == std::string::npos);
  }
  SECTION("check that SMARTS still works with v3K output") {
    std::unique_ptr<RWMol> mol(SmartsToMol("C[#8,#7]"));
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 2);
    bool includeStereo = true;
    int confId = -1;
    bool kekulize = true;
    bool forceV3000 = true;
    auto outBlock =
        MolToMolBlock(*mol, includeStereo, confId, kekulize, forceV3000);
    REQUIRE(outBlock.find(" Q ") == std::string::npos);
    REQUIRE(outBlock.find(" [O,N] ") != std::string::npos);
    std::unique_ptr<RWMol> mol2(MolBlockToMol(outBlock));
    REQUIRE(mol2);
    auto smarts = MolToSmarts(*mol2);
    REQUIRE(smarts == "[#6][#8,#7]");
  }
}
TEST_CASE(
    "Github #2229: problem round-tripping mol files with bond topology info",
    "[bug,writer]") {
  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/github2229_1.mol";
  std::unique_ptr<RWMol> mol(MolFileToMol(fName));
  REQUIRE(mol);
  REQUIRE(mol->getNumBonds() == 9);
  REQUIRE(!mol->getBondWithIdx(0)->hasQuery());
  REQUIRE(mol->getBondWithIdx(7)->hasQuery());
  SECTION("basics") {
    auto outBlock = MolToMolBlock(*mol);
    REQUIRE(outBlock.find(" 7  8  1  0  0  2") != std::string::npos);
    std::unique_ptr<RWMol> mol2(MolBlockToMol(outBlock));
    REQUIRE(mol2);
    REQUIRE(mol2->getNumBonds() == 9);
    REQUIRE(!mol2->getBondWithIdx(0)->hasQuery());
    REQUIRE(mol2->getBondWithIdx(7)->hasQuery());
  }
  SECTION("basics with v3k") {
    bool includeStereo = true;
    int confId = -1;
    bool kekulize = true;
    bool forceV3000 = true;
    auto outBlock =
        MolToMolBlock(*mol, includeStereo, confId, kekulize, forceV3000);
    REQUIRE(outBlock.find("1 7 8 TOPO=2") != std::string::npos);
    std::unique_ptr<RWMol> mol2(MolBlockToMol(outBlock));
    REQUIRE(mol2);
    REQUIRE(mol2->getNumBonds() == 9);
    REQUIRE(!mol2->getBondWithIdx(0)->hasQuery());
    REQUIRE(mol2->getBondWithIdx(7)->hasQuery());
  }
}
TEST_CASE("preserve mol file properties on bonds", "[parser,ctab]") {
  SECTION("basics") {
    std::string molblock = R"CTAB(
  Mrv1810 02111915042D          

  4  3  0  0  0  0            999 V2000
   -1.5625    1.6071    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8480    2.0196    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2770    2.0196    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5625    0.7821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0  0  0  0
  1  2  6  0  0  0  0
  1  4  1  1  0  0  0
M  END
      )CTAB";
    std::unique_ptr<ROMol> mol(MolBlockToMol(molblock));
    REQUIRE(mol);
    CHECK(mol->getBondWithIdx(1)->getProp<unsigned int>(
              common_properties::_MolFileBondType) == 6);
    CHECK(mol->getBondWithIdx(2)->getProp<unsigned int>(
              common_properties::_MolFileBondType) == 1);
    CHECK(mol->getBondWithIdx(2)->getProp<unsigned int>(
              common_properties::_MolFileBondStereo) == 1);
  }
  SECTION("basics-v3k") {
    std::string molblock = R"CTAB(
  Mrv1810 02111915102D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -2.9167 3 0 0
M  V30 2 C -1.583 3.77 0 0
M  V30 3 C -4.2503 3.77 0 0
M  V30 4 C -2.9167 1.46 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 3
M  V30 2 6 1 2
M  V30 3 1 1 4 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB";
    std::unique_ptr<ROMol> mol(MolBlockToMol(molblock));
    REQUIRE(mol);
    CHECK(mol->getBondWithIdx(1)->getProp<unsigned int>(
              common_properties::_MolFileBondType) == 6);
    CHECK(mol->getBondWithIdx(2)->getProp<unsigned int>(
              common_properties::_MolFileBondType) == 1);
    CHECK(mol->getBondWithIdx(2)->getProp<unsigned int>(
              common_properties::_MolFileBondCfg) == 1);
  }
}

TEST_CASE("github #2277 : Failure when parsing mol block with M PXA",
          "[parser,ctab]") {
  std::string molblock = R"CTAB(
  Mrv1810 02151911552D          

 13 12  0  0  1  0            999 V2000
   -3.6588  -26.0592    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -2.9453  -27.2971    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9453  -26.4713    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6588  -25.2360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9467  -24.8200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9467  -23.9968    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2304  -25.2358    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5102  -26.4716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7989  -25.2306    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7989  -26.0582    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3730  -26.4732    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2277  -26.0635    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0839  -26.4708    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0  0  0  0
  3  2  2  0  0  0  0
  1  4  1  1  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  5  7  1  0  0  0  0
  8 10  1  0  0  0  0
 10  9  2  0  0  0  0
  3 12  1  0  0  0  0
 12  8  1  0  0  0  0
 11  1  1  0  0  0  0
 10 13  1  0  0  0  0
M  PXA  11   -5.0817  -26.0408    0.0000 H
M  END
)CTAB";
  std::unique_ptr<ROMol> mol(MolBlockToMol(molblock));
  SECTION("basics, make sure we can parse the original data") {
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(10)->hasProp("_MolFile_PXA"));
    CHECK(!mol->getAtomWithIdx(11)->hasProp("_MolFile_PXA"));
  }
  SECTION("basics, can we write it?") {
    REQUIRE(mol);
    std::string outmb = MolToMolBlock(*mol);
    CHECK(outmb.find("M  PXA  11") != std::string::npos);
  }
}

TEST_CASE(
    "github #2266: missing stereo in adamantyl-like cages with "
    "exocyclic bonds",
    "[bug]") {
  SECTION("basics") {
    std::string molblock = R"CTAB(
        SciTegic12231509382D

 14 16  0  0  0  0            999 V2000
    1.5584   -5.7422    0.0000 C   0  0
    2.2043   -5.0535    0.0000 C   0  0  2  0  0  0
    2.3688   -5.5155    0.0000 C   0  0  1  0  0  0
    2.9210   -5.3181    0.0000 C   0  0
    3.1270   -5.8206    0.0000 C   0  0
    3.6744   -5.1312    0.0000 C   0  0  2  0  0  0
    2.3619   -4.6609    0.0000 C   0  0
    2.9268   -3.9939    0.0000 C   0  0  2  0  0  0
    2.1999   -4.2522    0.0000 C   0  0
    3.6803   -4.3062    0.0000 C   0  0
    2.9436   -3.1692    0.0000 N   0  0
    4.4569   -5.4095    0.0000 H   0  0
    2.3246   -6.3425    0.0000 H   0  0
    1.4365   -4.7500    0.0000 H   0  0
  1  2  1  0
  1  3  1  0
  2  4  1  0
  3  5  1  0
  4  6  1  0
  5  6  1  0
  7  8  1  0
  3  7  1  0
  2  9  1  0
  6 10  1  0
 10  8  1  0
  8  9  1  0
  8 11  1  1
  6 12  1  6
  3 13  1  1
  2 14  1  6
M  END)CTAB";
    {
      std::unique_ptr<ROMol> mol(MolBlockToMol(molblock));
      REQUIRE(mol);
      CHECK(mol->getNumAtoms() == 11);
      CHECK(mol->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
      CHECK(mol->getAtomWithIdx(2)->getChiralTag() != Atom::CHI_UNSPECIFIED);
      CHECK(mol->getAtomWithIdx(5)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    }
    {
      bool sanitize = true;
      bool removeHs = false;
      std::unique_ptr<ROMol> mol(MolBlockToMol(molblock, sanitize, removeHs));
      REQUIRE(mol);
      CHECK(mol->getNumAtoms() == 14);

      CHECK(mol->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
      CHECK(mol->getAtomWithIdx(2)->getChiralTag() != Atom::CHI_UNSPECIFIED);
      CHECK(mol->getAtomWithIdx(5)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    }
  }
  SECTION("with F") {
    std::string molblock = R"CTAB(
        SciTegic12231509382D

 14 16  0  0  0  0            999 V2000
    1.5584   -5.7422    0.0000 C   0  0
    2.2043   -5.0535    0.0000 C   0  0  2  0  0  0
    2.3688   -5.5155    0.0000 C   0  0  1  0  0  0
    2.9210   -5.3181    0.0000 C   0  0
    3.1270   -5.8206    0.0000 C   0  0
    3.6744   -5.1312    0.0000 C   0  0  2  0  0  0
    2.3619   -4.6609    0.0000 C   0  0
    2.9268   -3.9939    0.0000 C   0  0  2  0  0  0
    2.1999   -4.2522    0.0000 C   0  0
    3.6803   -4.3062    0.0000 C   0  0
    2.9436   -3.1692    0.0000 N   0  0
    4.4569   -5.4095    0.0000 F   0  0
    2.3246   -6.3425    0.0000 F   0  0
    1.4365   -4.7500    0.0000 F   0  0
  1  2  1  0
  1  3  1  0
  2  4  1  0
  3  5  1  0
  4  6  1  0
  5  6  1  0
  7  8  1  0
  3  7  1  0
  2  9  1  0
  6 10  1  0
 10  8  1  0
  8  9  1  0
  8 11  1  1
  6 12  1  6
  3 13  1  1
  2 14  1  6
M  END)CTAB";
    {
      std::unique_ptr<ROMol> mol(MolBlockToMol(molblock));
      REQUIRE(mol);
      CHECK(mol->getNumAtoms() == 14);
      CHECK(mol->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
      CHECK(mol->getAtomWithIdx(2)->getChiralTag() != Atom::CHI_UNSPECIFIED);
      CHECK(mol->getAtomWithIdx(5)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    }
  }
}

TEST_CASE("parsing of SCN lines", "[bug, sgroups]") {
  SECTION("basics") {
    std::string molblock = R"CTAB(
  MJ171200

 76 80  0  0  0  0  0  0  0  0999 V2000
   -6.4802    2.6494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8927    3.3638    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.7177    3.3638    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1302    2.6494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.7177    1.9349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8927    1.9349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2426    1.9349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6552    2.6494    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -6.4802    1.2203    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6552    1.2203    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0051    1.2203    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4176    1.9349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2426    0.5060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4176    0.5060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8927    0.5060    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4802   -0.2085    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6552   -0.2085    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1302    1.2203    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -6.8927   -0.9230    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4802   -1.6374    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6552   -1.6374    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0051   -1.6374    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2426   -2.3519    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4176   -2.3519    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0051   -3.0663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1801   -3.0663    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7676   -3.7808    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4176   -3.7808    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0051   -4.4953    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1801   -4.4953    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3243   -3.7791    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7368   -3.0648    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5619   -3.0648    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9744   -3.7791    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5619   -4.4936    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7368   -4.4936    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3243   -5.2082    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9744   -5.2082    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5619   -5.9226    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7368   -5.9226    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9744   -2.3503    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7994   -2.3503    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7994   -3.7791    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.9487   -5.1497    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3243   -6.6371    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3705   -2.2987    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9580   -1.5842    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1329   -1.5842    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2796   -2.2987    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1329   -3.0132    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9580   -3.0132    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1329   -0.1553    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2796   -0.8698    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3705   -3.7276    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2796   -3.7276    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1329   -4.4420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9580   -4.4420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3705   -5.1566    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8575   -2.2792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4450   -1.5648    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6200   -1.5648    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2073   -2.2792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6200   -2.9937    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4450   -2.9937    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3789   -3.8003    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3797   -5.2304    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8445   -5.2489    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3318   -2.2987    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2075   -3.7082    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8195   -3.7288    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.3448   -2.2597    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1706   -0.8729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3860   -0.6179    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
    3.2997   -0.0580    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1046   -2.2987    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9675   -0.6233    0.0000 Na  0  3  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  1  6  1  0  0  0  0
  9 10  1  0  0  0  0
  7 10  1  0  0  0  0
  8  1  1  0  0  0  0
  6  9  1  0  0  0  0
 11 12  1  0  0  0  0
 13 14  1  0  0  0  0
 11 14  2  0  0  0  0
 10 13  2  0  0  0  0
 12  7  2  0  0  0  0
 15 16  2  0  0  0  0
 16 17  1  0  0  0  0
 13 17  1  0  0  0  0
  5 18  1  0  0  0  0
 19 20  1  0  0  0  0
 20 21  1  0  0  0  0
 16 19  1  0  0  0  0
 23 24  1  0  0  0  0
 22 24  2  0  0  0  0
 21 23  1  0  0  0  0
 25 26  1  0  0  0  0
 25 24  1  1  0  0  0
 28 29  1  0  0  0  0
 29 30  1  0  0  0  0
 27 30  1  0  0  0  0
 25 28  1  0  0  0  0
 27 26  1  0  0  0  0
 31 32  1  0  0  0  0
 32 33  1  0  0  0  0
 33 34  1  0  0  0  0
 34 35  1  0  0  0  0
 35 36  1  0  0  0  0
 31 36  1  0  0  0  0
 39 40  2  0  0  0  0
 37 40  1  0  0  0  0
 35 38  1  1  0  0  0
 36 37  1  6  0  0  0
 41 42  1  0  0  0  0
 34 43  1  6  0  0  0
 33 41  1  1  0  0  0
 44 38  1  1  0  0  0
 40 45  1  0  0  0  0
 46 47  1  0  0  0  0
 47 48  1  0  0  0  0
 48 49  1  0  0  0  0
 49 50  1  0  0  0  0
 50 51  1  0  0  0  0
 46 51  1  0  0  0  0
 52 53  1  0  0  0  0
 48 53  1  1  0  0  0
 56 57  2  0  0  0  0
 54 57  1  0  0  0  0
 51 54  1  6  0  0  0
 50 55  1  1  0  0  0
 57 58  1  0  0  0  0
 59 60  1  0  0  0  0
 60 61  1  0  0  0  0
 61 62  1  0  0  0  0
 62 63  1  0  0  0  0
 63 64  1  0  0  0  0
 59 64  1  0  0  0  0
 28 65  1  6  0  0  0
 31 65  1  1  0  0  0
 29 66  1  1  0  0  0
 30 67  1  6  0  0  0
 27 55  1  1  0  0  0
 46 68  1  1  0  0  0
 62 68  1  6  0  0  0
 63 69  1  1  0  0  0
 64 70  1  6  0  0  0
 59 71  1  1  0  0  0
 61 72  1  1  0  0  0
 72 73  1  0  0  0  0
 72 74  2  0  0  0  0
 49 75  1  6  0  0  0
M  STY  3   1 SRU   2 SRU   3 SRU
M  SCN  3   1 HT    2 HT    3 HT
M  SAL   1 15  55  50  51  54  57  58  56  46  68  62  63  69  64  70  61
M  SAL   1 11  72  74  73  60  47  48  53  52  49  75  59
M  SMT   1 b
M  SBL   1  2  71  76
M  SAL   2 15  27  26  25  28  65  31  36  37  40  45  39  35  34  43  33
M  SAL   2 15  41  42  32  29  66  24  22  23  21  20  19  16  17  13  10
M  SAL   2 15   7  12  11   9   6   1   8   2   3   4   5  18  14  15  30
M  SAL   2  2  67  38
M  SMT   2 a
M  SBL   2  2  71  46
M  SAL   3 15  38  35  36  37  40  45  39  31  65  28  25  24  22  23  21
M  SAL   3 15  20  19  16  17  13  10   7  12  11   9   6   1   8   2   3
M  SAL   3 15   4   5  18  14  15  26  27  55  50  51  54  57  58  56  46
M  SAL   3 15  68  62  63  69  64  70  61  72  74  73  60  47  48  53  52
M  SAL   3 13  49  75  30  67  29  66  32  33  41  42  34  43  59
M  SMT   3 n
M  SBL   3  2  46  76
M  END
)CTAB";
    std::unique_ptr<ROMol> mol(MolBlockToMol(molblock));
    REQUIRE(mol);
  }
}

TEST_CASE("A couple more S group problems", "[bug, sgroups]") {
  std::string molblock = R"CTAB(CHEMBL3666739
      SciTegic05171617282D

 35 40  0  0  0  0            999 V2000
   -3.6559    5.8551    0.0000 O   0  0
   -2.6152    5.2576    0.0000 C   0  0
   -2.6120    3.7568    0.0000 N   0  0
   -3.9097    3.0028    0.0000 C   0  0
   -5.2093    3.7519    0.0000 C   0  0
   -6.5078    3.0010    0.0000 C   0  0
   -6.5067    1.5010    0.0000 C   0  0
   -5.2071    0.7519    0.0000 C   0  0
   -3.9086    1.5029    0.0000 C   0  0
   -2.6111    0.7486    0.0000 C   0  0
   -2.6111   -0.7486    0.0000 N   0  0
   -1.2964   -1.4973    0.0000 C   0  0
   -1.2907   -2.9981    0.0000 N   0  0
   -2.5870   -3.7544    0.0000 C   0  0
   -2.5748   -5.2506    0.0000 C   0  0
   -3.8815   -6.0264    0.0000 C   0  0
   -5.1819   -5.2707    0.0000 C   0  0
   -6.6004   -5.7374    0.0000 N   0  0
   -7.4849   -4.5227    0.0000 N   0  0
   -6.6189   -3.3309    0.0000 C   0  0
   -5.1934   -3.7757    0.0000 C   0  0
   -3.9049   -3.0000    0.0000 C   0  0
    0.0000   -0.7486    0.0000 C   0  0
    1.2964   -1.4973    0.0000 C   0  0
    2.5929   -0.7486    0.0000 C   0  0
    2.5929    0.7486    0.0000 C   0  0
    1.2964    1.4973    0.0000 C   0  0
    0.0000    0.7486    0.0000 C   0  0
   -1.2964    1.4973    0.0000 N   0  0
   -1.3175    6.0116    0.0000 C   0  0
   -1.3185    7.5117    0.0000 C   0  0
   -0.0200    8.2626    0.0000 C   0  0
    1.2795    7.5135    0.0000 N   0  0
    1.2806    6.0135    0.0000 C   0  0
   -0.0178    5.2626    0.0000 C   0  0
  1  2  2  0
  2  3  1  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  7  8  1  0
  8  9  2  0
  4  9  1  0
  9 10  1  0
 10 11  2  0
 11 12  1  0
 12 13  1  0
 13 14  1  0
 14 15  2  0
 15 16  1  0
 16 17  2  0
 17 18  1  0
 18 19  1  0
 19 20  2  0
 20 21  1  0
 17 21  1  0
 21 22  2  0
 14 22  1  0
 12 23  2  0
 23 24  1  0
 24 25  2  0
 25 26  1  0
 26 27  2  0
 27 28  1  0
 23 28  1  0
 28 29  2  0
 10 29  1  0
  2 30  1  0
 30 31  2  0
 31 32  1  0
 32 33  2  0
 33 34  1  0
 34 35  2  0
 30 35  1  0
M  STY  1   1 DAT
M  SLB  1   1   1
M  SAL   1  1  33
M  SDT   1 MRV_IMPLICIT_H
M  SDD   1     0.5304   -0.4125    DR    ALL  0       0
M  SED   1 IMPL_H1
M  END
)CTAB";
  SECTION("spaces in count lines") {
    std::unique_ptr<ROMol> mol(MolBlockToMol(molblock));
    REQUIRE(mol);
  }
  SECTION("short SDT lines") {
    std::unique_ptr<ROMol> mol(MolBlockToMol(molblock));
    REQUIRE(mol);
    const auto &sgroups = getSGroups(*mol);
    CHECK(sgroups.size() == 1);
    CHECK(sgroups[0].hasProp("TYPE"));
    CHECK(sgroups[0].getProp<std::string>("TYPE") == "DAT");
    CHECK(sgroups[0].hasProp("FIELDNAME"));
    CHECK(sgroups[0].getProp<std::string>("FIELDNAME") == "MRV_IMPLICIT_H");
  }
}
