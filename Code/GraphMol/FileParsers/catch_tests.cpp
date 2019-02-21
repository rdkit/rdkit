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
