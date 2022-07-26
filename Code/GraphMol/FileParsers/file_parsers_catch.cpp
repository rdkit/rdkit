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
#include <GraphMol/QueryAtom.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/FileParsers/SequenceParsers.h>
#include <GraphMol/FileParsers/SequenceWriters.h>
#include <GraphMol/FileParsers/PNGParser.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <RDGeneral/FileParseException.h>
#include <boost/algorithm/string.hpp>

using namespace RDKit;

TEST_CASE("Basic SVG Parsing", "[SVG][reader]") {
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
    "[bug][reader]") {
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
          "[bug][writer]") {
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
    "[bug][writer]") {
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
TEST_CASE("preserve mol file properties on bonds", "[reader][ctab]") {
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
          "[reader][ctab]") {
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

TEST_CASE("workaround for broken MJ2009-MJ2011 molblocks",
          "[feature][sgroups]") {
  SECTION("molblock1 strictParsing true/false") {
    std::string molblock1 = R"CTAB(
  MJ201100                      

 10 10  0  0  0  0  0  0  0  0999 V2000
   -1.2946    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2946   -1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467    1.2493    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.1342    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467   -0.1796    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6907    0.5348    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  6  8  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  1  0  0  0  0
  8 10  1  0  0  0  0
M  STY  1   1 SUP
M  SAL   1  4   7   8   9  10
M  SMT   1 CF3
M  SBL   1  1   7
M  SAP   1  1   8
M  END
)CTAB";
    std::string expectedMolblock1 = R"CTAB(
     RDKit          2D

 10 10  0  1  0  0  0  0  0  0999 V2000
   -1.2946    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2946   -1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467    1.2493    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.1342    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467   -0.1796    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6907    0.5348    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  6  8  1  0
  7  8  1  0
  8  9  1  0
  8 10  1  0
M  STY  1   1 SUP
M  SAL   1  4   7   8   9  10
M  SBL   1  1   7
M  SMT   1 CF3
M  SAP   1  1   8   6   
M  END
)CTAB";
    std::unique_ptr<ROMol> mol;
    REQUIRE_THROWS_AS(mol.reset(MolBlockToMol(molblock1)), FileParseException);
    REQUIRE_NOTHROW(mol.reset(MolBlockToMol(molblock1, true, true, false)));
    REQUIRE(mol);
    CHECK(MolToMolBlock(*mol) == expectedMolblock1);
  }
  SECTION("molblock1 strictParsing true/false no/bad SBL group") {
    std::string molblock1NoSBL = R"CTAB(
  MJ201100                      

 10 10  0  0  0  0  0  0  0  0999 V2000
   -1.2946    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2946   -1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467    1.2493    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.1342    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467   -0.1796    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6907    0.5348    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  6  8  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  1  0  0  0  0
  8 10  1  0  0  0  0
M  STY  1   1 SUP
M  SAL   1  4   7   8   9  10
M  SMT   1 CF3
M  SAP   1  1   8
M  END
)CTAB";
    std::string molblock1BadSBL = R"CTAB(
  MJ201100                      

 10 10  0  0  0  0  0  0  0  0999 V2000
   -1.2946    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2946   -1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467    1.2493    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.1342    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467   -0.1796    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6907    0.5348    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  6  8  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  1  0  0  0  0
  8 10  1  0  0  0  0
M  STY  1   1 SUP
M  SAL   1  4   7   8   9  10
M  SMT   1 CF3
M  SAP   1  1   8
M  SBL   1  2   7   8
M  END
)CTAB";
    std::string expectedMolblock1NoSGroups = R"CTAB(
     RDKit          2D

 10 10  0  0  0  0  0  0  0  0999 V2000
   -1.2946    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2946   -1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467    1.2493    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.1342    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467   -0.1796    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6907    0.5348    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  6  8  1  0
  7  8  1  0
  8  9  1  0
  8 10  1  0
M  END
)CTAB";
    std::unique_ptr<ROMol> mol;
    REQUIRE_THROWS_AS(mol.reset(MolBlockToMol(molblock1NoSBL)),
                      FileParseException);
    REQUIRE_NOTHROW(
        mol.reset(MolBlockToMol(molblock1NoSBL, true, true, false)));
    REQUIRE(mol);
    CHECK(MolToMolBlock(*mol) == expectedMolblock1NoSGroups);
    REQUIRE_THROWS_AS(mol.reset(MolBlockToMol(molblock1BadSBL)),
                      FileParseException);
    REQUIRE_NOTHROW(
        mol.reset(MolBlockToMol(molblock1BadSBL, true, true, false)));
    REQUIRE(mol);
    CHECK(MolToMolBlock(*mol) == expectedMolblock1NoSGroups);
  }
  SECTION("molblock2 strictParsing true/false") {
    std::string molblock2 = R"CTAB(
  MJ201100                      

 13 13  0  0  0  0  0  0  0  0999 V2000
   -1.2946    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2946   -1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467    1.2493    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.1342    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467   -0.1796    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6907    0.5348    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4380   -0.7027    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7235   -1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0091   -0.7027    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  6  8  1  0  0  0  0
  3 12  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  1  0  0  0  0
  8 10  1  0  0  0  0
 11 12  1  0  0  0  0
 12 13  2  0  0  0  0
M  STY  2   1 SUP   2 SUP
M  SAL   1  4   7   8   9  10
M  SMT   1 CF3
M  SBL   1  1   7
M  SAP   1  1   8
M  SAL   2  3  11  12  13
M  SMT   2 COOH
M  SBL   2  1   8
M  SAP   2  1  12
M  END
)CTAB";
    std::string expectedMolblock2 = R"CTAB(
     RDKit          2D

 13 13  0  2  0  0  0  0  0  0999 V2000
   -1.2946    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2946   -1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467    1.2493    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.1342    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467   -0.1796    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6907    0.5348    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4380   -0.7027    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7235   -1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0091   -0.7027    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  2  3  1  0
  3  4  2  0
  6  8  1  0
  3 12  1  0
  7  8  1  0
  8  9  1  0
  8 10  1  0
 11 12  1  0
 12 13  2  0
M  STY  2   1 SUP   2 SUP
M  SAL   1  4   7   8   9  10
M  SBL   1  1   7
M  SMT   1 CF3
M  SAP   1  1   8   6   
M  SAL   2  3  11  12  13
M  SBL   2  1   8
M  SMT   2 COOH
M  SAP   2  1  12   3   
M  END
)CTAB";
    std::unique_ptr<ROMol> mol;
    REQUIRE_THROWS_AS(mol.reset(MolBlockToMol(molblock2)), FileParseException);
    REQUIRE_NOTHROW(mol.reset(MolBlockToMol(molblock2, true, true, false)));
    REQUIRE(mol);
    CHECK(MolToMolBlock(*mol) == expectedMolblock2);
  }
  SECTION("molblock2 strictParsing true/false no/bad SBL group1") {
    std::string molblock2NoSBL = R"CTAB(
  MJ201100                      

 13 13  0  0  0  0  0  0  0  0999 V2000
   -1.2946    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2946   -1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467    1.2493    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.1342    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467   -0.1796    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6907    0.5348    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4380   -0.7027    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7235   -1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0091   -0.7027    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  6  8  1  0  0  0  0
  3 12  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  1  0  0  0  0
  8 10  1  0  0  0  0
 11 12  1  0  0  0  0
 12 13  2  0  0  0  0
M  STY  2   1 SUP   2 SUP
M  SAL   1  4   7   8   9  10
M  SMT   1 CF3
M  SAP   1  1   8
M  SAL   2  3  11  12  13
M  SMT   2 COOH
M  SBL   2  1   8
M  SAP   2  1  12
M  END
)CTAB";
    std::string molblock2BadSBL = R"CTAB(
  MJ201100                      

 13 13  0  0  0  0  0  0  0  0999 V2000
   -1.2946    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2946   -1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467    1.2493    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.1342    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467   -0.1796    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6907    0.5348    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4380   -0.7027    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7235   -1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0091   -0.7027    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  6  8  1  0  0  0  0
  3 12  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  1  0  0  0  0
  8 10  1  0  0  0  0
 11 12  1  0  0  0  0
 12 13  2  0  0  0  0
M  STY  2   1 SUP   2 SUP
M  SAL   1  4   7   8   9  10
M  SMT   1 CF3
M  SAP   1  1   8
M  SBL   1  2   7   8
M  SAL   2  3  11  12  13
M  SMT   2 COOH
M  SBL   2  1   8
M  SAP   2  1  12
M  END
)CTAB";
    std::string expectedMolblock2NoSGroup1 = R"CTAB(
     RDKit          2D

 13 13  0  1  0  0  0  0  0  0999 V2000
   -1.2946    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0090   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2946   -1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5801    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467    1.2493    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.1342    0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5467   -0.1796    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6907    0.5348    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4380   -0.7027    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7235   -1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0091   -0.7027    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  2  3  1  0
  3  4  2  0
  6  8  1  0
  3 12  1  0
  7  8  1  0
  8  9  1  0
  8 10  1  0
 11 12  1  0
 12 13  2  0
M  STY  1   1 SUP
M  SAL   1  3  11  12  13
M  SBL   1  1   8
M  SMT   1 COOH
M  SAP   1  1  12   3   
M  END
)CTAB";
    std::unique_ptr<ROMol> mol;
    REQUIRE_THROWS_AS(mol.reset(MolBlockToMol(molblock2NoSBL)),
                      FileParseException);
    REQUIRE_NOTHROW(
        mol.reset(MolBlockToMol(molblock2NoSBL, true, true, false)));
    REQUIRE(mol);
    CHECK(MolToMolBlock(*mol) == expectedMolblock2NoSGroup1);
    REQUIRE_THROWS_AS(mol.reset(MolBlockToMol(molblock2BadSBL)),
                      FileParseException);
    REQUIRE_NOTHROW(
        mol.reset(MolBlockToMol(molblock2BadSBL, true, true, false)));
    REQUIRE(mol);
    CHECK(MolToMolBlock(*mol) == expectedMolblock2NoSGroup1);
  }
}

TEST_CASE(
    "do not throw but remove malformed V2000 SGroups when strictParsing is "
    "false",
    "[feature][sgroups]") {
  std::string molblock = R"CTAB(
  ChemDraw01072117362D

 16 16  0  0  0  0  0  0  0  0999 V2000
   -1.7862    0.7219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7862   -0.1031    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0717   -0.5156    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3572   -0.1031    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3572    0.7219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0717    1.1344    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3572    2.3719    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3572   -0.5156    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0717   -0.1031    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3572   -1.3406    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0717   -1.7531    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0717   -2.5781    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7862   -1.3406    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4289    1.3406    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7145    1.7531    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7145    2.5781    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
 11 12  1  0      
 11 13  1  0      
 10 11  1  0      
  8  9  2  0      
  8 10  1  0      
  4  8  1  0      
 15 16  1  0      
 14 15  1  0      
  6 14  1  0      
  7 16  1  0      
  1  2  2  0      
  2  3  1  0      
  3  4  2  0      
  4  5  1  0      
  5  6  2  0      
  6  1  1  0      
M  STY  1   1 SUP
M  SLB  1   1   1
M  SAL   1  6   8   9  10  11  12  13
M  SBL   1  1   6
M  SMT   1 COOiPr
M  SBV   1   6   -0.7145    0.4125
M  STY  1   2 SUP
M  SLB  1   2   2
M  SAL   2  3  14  15  16
M  SBL   2  2   9  10
M  SMT   2 (CH2)3
M  SBV   2   9    0.3572   -0.2062
M  SBV   2  10    0.3572   -0.2062
M  END
)CTAB";
  SECTION("molblock strictParsing true") {
    std::unique_ptr<ROMol> mol(MolBlockToMol(molblock));
    REQUIRE(mol);
    CHECK(getSubstanceGroups(*mol).size() == 2);
  }
  SECTION("molblock bad sgroup idx") {
    std::string molblockBad = boost::replace_all_copy(
        molblock, "M  SBL   1  1   6", "M  SBL   3  1   6");
    std::unique_ptr<ROMol> mol;
    REQUIRE_THROWS(mol.reset(MolBlockToMol(molblockBad)));
    REQUIRE_NOTHROW(mol.reset(MolBlockToMol(molblockBad, true, true, false)));
    CHECK(getSubstanceGroups(*mol).size() == 1);
  }
  SECTION("molblock sgroup line too short (1)") {
    std::string molblockBad = boost::replace_all_copy(
        molblock, "M  SBV   1   6   -0.7145    0.4125", "M  SBV   1   6");
    std::unique_ptr<ROMol> mol;
    REQUIRE_THROWS(mol.reset(MolBlockToMol(molblockBad)));
    REQUIRE_NOTHROW(mol.reset(MolBlockToMol(molblockBad, true, true, false)));
    CHECK(getSubstanceGroups(*mol).size() == 1);
  }
  SECTION("molblock sgroup line too short (2)") {
    std::string molblockBad = boost::replace_all_copy(
        molblock, "M  SBL   2  2   9  10", "M  SBL   2  3   9  10");
    std::unique_ptr<ROMol> mol;
    REQUIRE_THROWS(mol.reset(MolBlockToMol(molblockBad)));
    REQUIRE_NOTHROW(mol.reset(MolBlockToMol(molblockBad, true, true, false)));
    CHECK(getSubstanceGroups(*mol).size() == 1);
  }
  SECTION("molblock sgroup bad bond idx") {
    std::string molblockBad = boost::replace_all_copy(
        molblock, "M  SBL   2  2   9  10", "M  SBL   2  2   9  99");
    std::unique_ptr<ROMol> mol;
    REQUIRE_THROWS(mol.reset(MolBlockToMol(molblockBad)));
    REQUIRE_NOTHROW(mol.reset(MolBlockToMol(molblockBad, true, true, false)));
    CHECK(getSubstanceGroups(*mol).size() == 1);
  }
  SECTION("molblock sgroup bad atom idx") {
    std::string molblockBad = boost::replace_all_copy(
        molblock, "M  SAL   2  3  14  15  16", "M  SAL   2  3  14  15  99");
    std::unique_ptr<ROMol> mol;
    REQUIRE_THROWS(mol.reset(MolBlockToMol(molblockBad)));
    REQUIRE_NOTHROW(mol.reset(MolBlockToMol(molblockBad, true, true, false)));
    CHECK(getSubstanceGroups(*mol).size() == 1);
  }
}

TEST_CASE(
    "do not throw but remove malformed V3000 SGroups when strictParsing is "
    "false",
    "[feature][sgroups]") {
  std::string molblock = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 16 16 2 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.786200 0.721900 0.000000 0
M  V30 2 C -1.786200 -0.103100 0.000000 0
M  V30 3 C -1.071700 -0.515600 0.000000 0
M  V30 4 C -0.357200 -0.103100 0.000000 0
M  V30 5 C -0.357200 0.721900 0.000000 0
M  V30 6 C -1.071700 1.134400 0.000000 0
M  V30 7 O -0.357200 2.371900 0.000000 0
M  V30 8 C 0.357200 -0.515600 0.000000 0
M  V30 9 O 1.071700 -0.103100 0.000000 0
M  V30 10 O 0.357200 -1.340600 0.000000 0
M  V30 11 C 1.071700 -1.753100 0.000000 0
M  V30 12 C 1.071700 -2.578100 0.000000 0
M  V30 13 C 1.786200 -1.340600 0.000000 0
M  V30 14 C -1.428900 1.340600 0.000000 0
M  V30 15 C -0.714500 1.753100 0.000000 0
M  V30 16 C -0.714500 2.578100 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 11 12
M  V30 2 1 11 13
M  V30 3 1 10 11
M  V30 4 2 8 9
M  V30 5 1 8 10
M  V30 6 1 4 8
M  V30 7 1 15 16
M  V30 8 1 14 15
M  V30 9 1 6 14
M  V30 10 1 7 16
M  V30 11 2 1 2
M  V30 12 1 2 3
M  V30 13 2 3 4
M  V30 14 1 4 5
M  V30 15 2 5 6
M  V30 16 1 6 1
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 1 ATOMS=(6 8 9 10 11 12 13) XBONDS=(1 6) LABEL=COOiPr CSTATE=(4 6 -
M  V30 -0.7145 0.4125 0)
M  V30 2 SUP 2 ATOMS=(3 14 15 16) XBONDS=(2 9 10) LABEL=(CH2)3 CSTATE=(4 9 0.35-
M  V30 72 -0.2062 0) CSTATE=(4 10 0.3572 -0.2062 0)
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB";
  SECTION("molblock strictParsing true") {
    std::unique_ptr<ROMol> mol(MolBlockToMol(molblock));
    REQUIRE(mol);
    CHECK(getSubstanceGroups(*mol).size() == 2);
  }
  SECTION("molblock sgroup line too short (1)") {
    std::string molblockBad =
        boost::replace_all_copy(molblock, "XBONDS=(1 6)", "XBONDS=(2 6)");
    std::unique_ptr<ROMol> mol;
    REQUIRE_THROWS(mol.reset(MolBlockToMol(molblockBad)));
    REQUIRE_NOTHROW(mol.reset(MolBlockToMol(molblockBad, true, true, false)));
    CHECK(getSubstanceGroups(*mol).size() == 1);
  }
  SECTION("molblock sgroup line too short (2)") {
    std::string molblockBad = boost::replace_all_copy(
        molblock, "ATOMS=(6 8 9 10 11 12 13)", "ATOMS=(7 8 9 10 11 12 13)");
    std::unique_ptr<ROMol> mol;
    REQUIRE_THROWS(mol.reset(MolBlockToMol(molblockBad)));
    REQUIRE_NOTHROW(mol.reset(MolBlockToMol(molblockBad, true, true, false)));
    CHECK(getSubstanceGroups(*mol).size() == 1);
  }
  SECTION("molblock sgroup bad bond idx") {
    std::string molblockBad =
        boost::replace_all_copy(molblock, "XBONDS=(2 9 10)", "XBONDS=(2 9 99)");
    std::unique_ptr<ROMol> mol;
    REQUIRE_THROWS(mol.reset(MolBlockToMol(molblockBad)));
    REQUIRE_NOTHROW(mol.reset(MolBlockToMol(molblockBad, true, true, false)));
    CHECK(getSubstanceGroups(*mol).size() == 1);
  }
  SECTION("molblock sgroup bad atom idx") {
    std::string molblockBad = boost::replace_all_copy(
        molblock, "ATOMS=(3 14 15 16)", "ATOMS=(3 14 15 99)");
    std::unique_ptr<ROMol> mol;
    REQUIRE_THROWS(mol.reset(MolBlockToMol(molblockBad)));
    REQUIRE_NOTHROW(mol.reset(MolBlockToMol(molblockBad, true, true, false)));
    CHECK(getSubstanceGroups(*mol).size() == 1);
  }
}

TEST_CASE("parsing of SCN lines", "[bug][sgroups]") {
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

TEST_CASE("A couple more S group problems", "[bug][sgroups]") {
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
M  SDT   1 FAKE_MRV_IMPLICIT_H
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
    const auto &sgroups = getSubstanceGroups(*mol);
    CHECK(sgroups.size() == 1);
    CHECK(sgroups[0].hasProp("TYPE"));
    CHECK(sgroups[0].getProp<std::string>("TYPE") == "DAT");
    CHECK(sgroups[0].hasProp("FIELDNAME"));
    CHECK(sgroups[0].getProp<std::string>("FIELDNAME") ==
          "FAKE_MRV_IMPLICIT_H");
  }
}

TEST_CASE("Github #2527: handling of \"R\" in CTABs", "[rgroups]") {
  std::string molblock = R"CTAB(example
  Mrv1902 07031913362D          

  2  1  0  0  0  0            999 V2000
   -1.1418    0.0687    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
   -1.9668    0.0687    0.0000 R   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
)CTAB";
  SECTION("basics") {
    bool sanitize = false;
    std::unique_ptr<ROMol> mol(MolBlockToMol(molblock, sanitize));
    REQUIRE(mol);
    auto *at = static_cast<QueryAtom *>(mol->getAtomWithIdx(1));
    REQUIRE(at->hasQuery());
    CHECK(at->getQuery()->getDescription() == "AtomNull");
  }
}

TEST_CASE("CML writer", "[CML][writer]") {
  SECTION("basics") {
    std::unique_ptr<RWMol> mol{new RWMol{}};
    mol->setProp(common_properties::_Name, "S-lactic acid");

    for (auto z : {6u, 1u, 1u, 1u, 6u, 1u, 8u, 1u, 6u, 8u, 8u}) {
      auto *a = new Atom{z};
      mol->addAtom(a, false, true);
    }
    mol->getAtomWithIdx(7u)->setIsotope(2u);
    mol->getAtomWithIdx(10u)->setFormalCharge(-1);

    mol->addBond(0u, 1u, Bond::SINGLE);
    mol->addBond(0u, 2u, Bond::SINGLE);
    mol->addBond(0u, 3u, Bond::SINGLE);
    mol->addBond(0u, 4u, Bond::SINGLE);
    mol->addBond(4u, 5u, Bond::SINGLE);
    mol->addBond(4u, 6u, Bond::SINGLE);
    mol->addBond(4u, 8u, Bond::SINGLE);
    mol->addBond(6u, 7u, Bond::SINGLE);
    mol->addBond(8u, 9u, Bond::DOUBLE);
    mol->addBond(8u, 10u, Bond::SINGLE);

    auto *conf = new Conformer{11u};
    conf->setId(0u);

    conf->setAtomPos(0u, RDGeom::Point3D{-0.95330, 0.60416, 1.01609});
    conf->setAtomPos(1u, RDGeom::Point3D{-1.00832, 1.68746, 0.83520});
    conf->setAtomPos(2u, RDGeom::Point3D{-1.96274, 0.16103, 0.94471});
    conf->setAtomPos(3u, RDGeom::Point3D{-0.57701, 0.44737, 2.04167});
    conf->setAtomPos(4u, RDGeom::Point3D{0.00000, 0.00000, 0.00000});
    conf->setAtomPos(5u, RDGeom::Point3D{-0.43038, 0.18596, -1.01377});
    conf->setAtomPos(6u, RDGeom::Point3D{0.22538, -1.36531, 0.19373});
    conf->setAtomPos(7u, RDGeom::Point3D{1.21993, -1.33937, 0.14580});
    conf->setAtomPos(8u, RDGeom::Point3D{1.38490, 0.73003, 0.00000});
    conf->setAtomPos(9u, RDGeom::Point3D{1.38490, 1.96795, 0.00000});
    conf->setAtomPos(10u, RDGeom::Point3D{2.35253, -0.07700, 0.00000});

    mol->addConformer(conf);

    mol->updatePropertyCache();
    MolOps::assignStereochemistryFrom3D(*mol);

    const std::string cmlblock = MolToCMLBlock(*mol);
    const std::string cmlblock_expected =
        R"CML(<?xml version="1.0" encoding="utf-8"?>
<cml xmlns="http://www.xml-cml.org/schema" xmlns:convention="http://www.xml-cml.org/convention/" convention="convention:molecular">
  <molecule id="m-1" formalCharge="-1" spinMultiplicity="1">
    <name>S-lactic acid</name>
    <atomArray>
      <atom id="a0" elementType="C" formalCharge="0" hydrogenCount="3" x3="-0.953300" y3="0.604160" z3="1.016090"/>
      <atom id="a1" elementType="H" formalCharge="0" hydrogenCount="0" x3="-1.008320" y3="1.687460" z3="0.835200"/>
      <atom id="a2" elementType="H" formalCharge="0" hydrogenCount="0" x3="-1.962740" y3="0.161030" z3="0.944710"/>
      <atom id="a3" elementType="H" formalCharge="0" hydrogenCount="0" x3="-0.577010" y3="0.447370" z3="2.041670"/>
      <atom id="a4" elementType="C" formalCharge="0" hydrogenCount="1" x3="0.000000" y3="0.000000" z3="0.000000">
        <atomParity atomRefs4="a0 a5 a6 a8">1</atomParity>
      </atom>
      <atom id="a5" elementType="H" formalCharge="0" hydrogenCount="0" x3="-0.430380" y3="0.185960" z3="-1.013770"/>
      <atom id="a6" elementType="O" formalCharge="0" hydrogenCount="1" x3="0.225380" y3="-1.365310" z3="0.193730"/>
      <atom id="a7" elementType="H" formalCharge="0" hydrogenCount="0" isotopeNumber="2" x3="1.219930" y3="-1.339370" z3="0.145800"/>
      <atom id="a8" elementType="C" formalCharge="0" hydrogenCount="0" x3="1.384900" y3="0.730030" z3="0.000000"/>
      <atom id="a9" elementType="O" formalCharge="0" hydrogenCount="0" x3="1.384900" y3="1.967950" z3="0.000000"/>
      <atom id="a10" elementType="O" formalCharge="-1" hydrogenCount="0" x3="2.352530" y3="-0.077000" z3="0.000000"/>
    </atomArray>
    <bondArray>
      <bond atomRefs2="a0 a1" id="b0" order="S"/>
      <bond atomRefs2="a0 a2" id="b1" order="S"/>
      <bond atomRefs2="a0 a3" id="b2" order="S"/>
      <bond atomRefs2="a0 a4" id="b3" order="S"/>
      <bond atomRefs2="a4 a5" id="b4" order="S" bondStereo="H"/>
      <bond atomRefs2="a4 a6" id="b5" order="S"/>
      <bond atomRefs2="a4 a8" id="b6" order="S"/>
      <bond atomRefs2="a6 a7" id="b7" order="S"/>
      <bond atomRefs2="a8 a9" id="b8" order="D"/>
      <bond atomRefs2="a8 a10" id="b9" order="S"/>
    </bondArray>
  </molecule>
</cml>
)CML";
    CHECK(cmlblock == cmlblock_expected);
  }

  SECTION("chirality1") {
    auto mol = R"CTAB(
  Mrv1921 04232106262D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 1
M  V30 BEGIN ATOM
M  V30 1 I 0.8918 -1.0472 0 0
M  V30 2 C 0.8918 0.4928 0 0
M  V30 3 Br 0.8918 2.0328 0 0
M  V30 4 F 2.4318 0.4928 0 0
M  V30 5 Cl -0.6482 0.4928 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 2 3 CFG=3
M  V30 3 1 2 4 CFG=1
M  V30 4 1 2 5
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    const std::string cmlblock = MolToCMLBlock(*mol);
    const std::string cmlblock_expected =
        R"CML(<?xml version="1.0" encoding="utf-8"?>
<cml xmlns="http://www.xml-cml.org/schema" xmlns:convention="http://www.xml-cml.org/convention/" convention="convention:molecular">
  <molecule id="m-1" formalCharge="0" spinMultiplicity="1">
    <atomArray>
      <atom id="a0" elementType="I" formalCharge="0" hydrogenCount="0" x2="0.891800" y2="-1.047200"/>
      <atom id="a1" elementType="C" formalCharge="0" hydrogenCount="0" x2="0.891800" y2="0.492800">
        <atomParity atomRefs4="a0 a2 a3 a4">1</atomParity>
      </atom>
      <atom id="a2" elementType="Br" formalCharge="0" hydrogenCount="0" x2="0.891800" y2="2.032800"/>
      <atom id="a3" elementType="F" formalCharge="0" hydrogenCount="0" x2="2.431800" y2="0.492800"/>
      <atom id="a4" elementType="Cl" formalCharge="0" hydrogenCount="0" x2="-0.648200" y2="0.492800"/>
    </atomArray>
    <bondArray>
      <bond atomRefs2="a0 a1" id="b0" order="S"/>
      <bond atomRefs2="a1 a2" id="b1" order="S"/>
      <bond atomRefs2="a1 a3" id="b2" order="S" bondStereo="W"/>
      <bond atomRefs2="a1 a4" id="b3" order="S"/>
    </bondArray>
  </molecule>
</cml>
)CML";
    CHECK(cmlblock == cmlblock_expected);
  }
  SECTION("chirality2") {
    auto mol = R"CTAB(
  Mrv1921 04232106262D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 1
M  V30 BEGIN ATOM
M  V30 1 I 0.8918 -1.0472 0 0
M  V30 2 C 0.8918 0.4928 0 0
M  V30 3 Br 0.8918 2.0328 0 0
M  V30 4 F 2.4318 0.4928 0 0
M  V30 5 Cl -0.6482 0.4928 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 2 3 CFG=1
M  V30 3 1 2 4 CFG=3
M  V30 4 1 2 5
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    const std::string cmlblock = MolToCMLBlock(*mol);
    const std::string cmlblock_expected =
        R"CML(<?xml version="1.0" encoding="utf-8"?>
<cml xmlns="http://www.xml-cml.org/schema" xmlns:convention="http://www.xml-cml.org/convention/" convention="convention:molecular">
  <molecule id="m-1" formalCharge="0" spinMultiplicity="1">
    <atomArray>
      <atom id="a0" elementType="I" formalCharge="0" hydrogenCount="0" x2="0.891800" y2="-1.047200"/>
      <atom id="a1" elementType="C" formalCharge="0" hydrogenCount="0" x2="0.891800" y2="0.492800">
        <atomParity atomRefs4="a0 a2 a3 a4">-1</atomParity>
      </atom>
      <atom id="a2" elementType="Br" formalCharge="0" hydrogenCount="0" x2="0.891800" y2="2.032800"/>
      <atom id="a3" elementType="F" formalCharge="0" hydrogenCount="0" x2="2.431800" y2="0.492800"/>
      <atom id="a4" elementType="Cl" formalCharge="0" hydrogenCount="0" x2="-0.648200" y2="0.492800"/>
    </atomArray>
    <bondArray>
      <bond atomRefs2="a0 a1" id="b0" order="S"/>
      <bond atomRefs2="a1 a2" id="b1" order="S"/>
      <bond atomRefs2="a1 a3" id="b2" order="S" bondStereo="H"/>
      <bond atomRefs2="a1 a4" id="b3" order="S"/>
    </bondArray>
  </molecule>
</cml>
)CML";
    CHECK(cmlblock == cmlblock_expected);
  }

  SECTION("no conformer") {
    auto mol = "C[C@](O)(F)Cl"_smiles;
    REQUIRE(mol);
    const std::string cmlblock = MolToCMLBlock(*mol);
    const std::string cmlblock_expected =
        R"CML(<?xml version="1.0" encoding="utf-8"?>
<cml xmlns="http://www.xml-cml.org/schema" xmlns:convention="http://www.xml-cml.org/convention/" convention="convention:molecular">
  <molecule id="m-1" formalCharge="0" spinMultiplicity="1">
    <atomArray>
      <atom id="a0" elementType="C" formalCharge="0" hydrogenCount="3"/>
      <atom id="a1" elementType="C" formalCharge="0" hydrogenCount="0">
        <atomParity atomRefs4="a0 a2 a3 a4">1</atomParity>
      </atom>
      <atom id="a2" elementType="O" formalCharge="0" hydrogenCount="1"/>
      <atom id="a3" elementType="F" formalCharge="0" hydrogenCount="0"/>
      <atom id="a4" elementType="Cl" formalCharge="0" hydrogenCount="0"/>
    </atomArray>
    <bondArray>
      <bond atomRefs2="a0 a1" id="b0" order="S"/>
      <bond atomRefs2="a1 a2" id="b1" order="S"/>
      <bond atomRefs2="a1 a3" id="b2" order="S"/>
      <bond atomRefs2="a1 a4" id="b3" order="S"/>
    </bondArray>
  </molecule>
</cml>
)CML";
    CHECK(cmlblock == cmlblock_expected);
  }
}

TEST_CASE("XYZ", "[XYZ][writer]") {
  SECTION("basics") {
    std::unique_ptr<RWMol> mol{new RWMol{}};
    mol->setProp(common_properties::_Name,
                 "methane\nthis part should not be output");

    for (unsigned z : {6, 1, 1, 1, 1}) {
      auto *a = new Atom{z};
      mol->addAtom(a, false, true);
    }

    auto *conf = new Conformer{5};
    conf->setId(0);
    conf->setAtomPos(0, RDGeom::Point3D{0.000, 0.000, 0.000});
    conf->setAtomPos(1, RDGeom::Point3D{-0.635, -0.635, 0.635});
    conf->setAtomPos(2, RDGeom::Point3D{-0.635, 0.635, -0.635});
    conf->setAtomPos(3, RDGeom::Point3D{0.635, -0.635, -0.635});
    conf->setAtomPos(4, RDGeom::Point3D{0.635, 0.635, 0.635});
    mol->addConformer(conf);

    const std::string xyzblock = MolToXYZBlock(*mol);
    std::string xyzblock_expected = R"XYZ(5
methane
C      0.000000    0.000000    0.000000
H     -0.635000   -0.635000    0.635000
H     -0.635000    0.635000   -0.635000
H      0.635000   -0.635000   -0.635000
H      0.635000    0.635000    0.635000
)XYZ";
    CHECK(xyzblock == xyzblock_expected);
  }
}

TEST_CASE("valence writing 1", "[bug][writer]") {
  SECTION("carbon") {
    std::string molblock = R"CTAB(carbon atom


  1  0  0  0  0  0            999 V2000
   -0.3958   -0.0542    0.0000 C   0  0  0  0  0 15
M  END)CTAB";
    bool sanitize = false;
    bool removeHs = false;
    std::unique_ptr<ROMol> mol(MolBlockToMol(molblock, sanitize, removeHs));
    REQUIRE(mol);
    mol->updatePropertyCache();
    CHECK(mol->getAtomWithIdx(0)->getNoImplicit());
    CHECK(mol->getAtomWithIdx(0)->getExplicitValence() == 0);
    CHECK(mol->getAtomWithIdx(0)->getTotalValence() == 0);
    auto outBlock = MolToMolBlock(*mol);
    REQUIRE(outBlock.find("0  0 15") != std::string::npos);
  }
  SECTION("P valences") {
    std::string molblock = R"CTAB(H2PO2


  3  2  0  0  0  0            999 V2000
    0.2667   -0.4167    0.0000 P   0  0  0  0  0  5
    0.2667    1.1083    0.0000 O   0  0
   -1.0958   -1.0042    0.0000 O   0  0
  2  1  2  0
  3  1  1  0
M  END)CTAB";
    bool sanitize = false;
    bool removeHs = false;
    std::unique_ptr<ROMol> mol(MolBlockToMol(molblock, sanitize, removeHs));
    REQUIRE(mol);
    mol->updatePropertyCache();
    CHECK(mol->getAtomWithIdx(0)->getNoImplicit());
    CHECK(mol->getAtomWithIdx(0)->getExplicitValence() == 5);
    CHECK(mol->getAtomWithIdx(0)->getTotalValence() == 5);
    auto outBlock = MolToMolBlock(*mol);
    REQUIRE(outBlock.find("0  0  5") != std::string::npos);
  }
}

TEST_CASE("Github #2695: Error when a squiggle bond is in an aromatic ring",
          "[bug][reader]") {
  SECTION("reported") {
    auto ctab = R"CTAB(
  -ISIS-  -- StrEd -- 

 19 22  0  0  0  0  0  0  0  0999 V2000
   -3.1355   -0.9331    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6355   -1.7990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6356   -1.7990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1356   -0.9331    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6356   -0.0671    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6355   -0.0671    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1355    0.7991    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6355    1.6651    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6356    1.6651    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1356    0.7991    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   -0.1354    0.7991    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4523    1.6080    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4034    1.2991    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2693    1.7990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1355    1.2991    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1355    0.2991    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.2693   -0.2011    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4034    0.2991    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4523   -0.0099    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  1  6  1  0  0  0  0
  6  7  2  0  0  0  0
  7  8  1  0  0  0  0
  8  9  2  0  0  0  0
 10  9  1  4  0  0  0
  5 10  1  0  0  0  0
 10 11  2  0  0  0  0
 11 12  1  0  0  0  0
 12 13  2  0  0  0  0
 13 14  1  0  0  0  0
 14 15  2  0  0  0  0
 15 16  1  0  0  0  0
 16 17  2  0  0  0  0
 17 18  1  0  0  0  0
 13 18  1  0  0  0  0
 18 19  2  0  0  0  0
 11 19  1  0  0  0  0
M  END)CTAB";
    std::unique_ptr<ROMol> mol(MolBlockToMol(ctab));
    REQUIRE(mol);
  }
}

TEST_CASE("Github #2917: _ctab _mol2 and _pdb support", "[feature][reader]") {
  SECTION("_ctab") {
    auto mol = R"CTAB(
  Mrv1810 01292008292D          

  4  3  0  0  0  0            999 V2000
   -3.7669    1.1053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0524    1.5178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3380    1.1053    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4814    1.5178    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  2  2  0  0  0  0
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 4);
  }
  SECTION("_ctab failure") {
    auto mol = R"CTAB(
  Mrv1810 01292008292D          

  4  3  0  0  0  0            999 V2000
   -3.7669    1.1053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0524    1.5178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3380    1.1053    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4814    1.5178    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  3  1  1  0  0  0  0
M  END
)CTAB"_ctab;
    REQUIRE(!mol);
  }
  SECTION("_pdb") {
    auto mol = R"DATA(HEADER    2VNF_PROTEIN
COMPND    2VNF_PROTEIN
REMARK    GENERATED BY X-TOOL on Wed Nov 21 18:02:19 2012
ATOM      1  N   ALA A 225      10.250 -13.177   9.152  1.00 19.76           N
ATOM      2  H   ALA A 225      10.605 -14.082   8.782  1.00  0.00           H
ATOM      3  CA  ALA A 225      11.136 -12.000   9.236  1.00 21.97           C
ATOM      4  HA  ALA A 225      11.079 -11.589  10.244  1.00  0.00           H
ATOM      5  C   ALA A 225      10.683 -10.934   8.231  1.00 21.61           C
ATOM      6  O   ALA A 225      10.811  -9.723   8.485  1.00 20.83           O
ATOM      7  CB  ALA A 225      12.572 -12.399   8.956  1.00 22.98           C
ATOM      8  HB1 ALA A 225      12.892 -13.138   9.690  1.00  0.00           H
ATOM      9  HB2 ALA A 225      12.641 -12.825   7.955  1.00  0.00           H
ATOM     10  HB3 ALA A 225      13.212 -11.519   9.022  1.00  0.00           H
TER      11      ALA A 225
HETATM   12  O   HOH    43      12.371  -9.746   8.354  0.50 30.13           O
END
)DATA"_pdb;
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 6);
  }
  SECTION("_pdb failure") {
    auto mol_ok = R"DATA(HEADER    TEST
COMPND    ACY
REMARK    invented
HETATM 2779  C   ACY A 404      15.911  -4.912  26.073  1.00 30.04           C  
HETATM 2780  O   ACY A 404      15.855  -4.063  25.124  1.00 24.12           O  
HETATM 2781  OXT ACY A 404      16.514  -4.578  27.173  1.00 34.44           O  
HETATM 2782  CH3 ACY A 404      15.319  -6.258  25.820  1.00 30.60           C  
CONECT 2780 2781
END
)DATA"_pdb;
    REQUIRE(mol_ok);
    CHECK(mol_ok->getNumAtoms() == 4);
    CHECK(mol_ok->getBondBetweenAtoms(1, 2));

    auto mol_fail = R"DATA(HEADER    TEST
COMPND    ACY
REMARK    invented
HETATM 2779  C   ACY A 404      15.911  -4.912  26.073  1.00 30.04           C  
HETATM 2780  O   ACY A 404      15.855  -4.063  25.124  1.00 24.12           O  
HETATM 2781  OXT ACY A 404      16.514  -4.578  27.173  1.00 34.44           O  
HETATM 2782  CH3 ACY A 404      15.319  -6.258  25.820  1.00 30.60           C  
CONECT 2780 2781
CONECT 2780 2779
CONECT 2780 2782
END
)DATA"_pdb;
    REQUIRE(!mol_fail);
  }
  SECTION("_mol2") {
    auto mol = R"DATA(@<TRIPOS>MOLECULE
UNK
6 4 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
   1   Na        0.0000    0.0000    0.0000 Na    1 UNL  1.0000
   2    C        0.0000    0.0000    0.0000 C.3   1 UNL -0.0305
   3    H        0.0000    0.0000    0.0000 H     1 UNL  0.0265
   4    H        0.0000    0.0000    0.0000 H     1 UNL  0.0265
   5    H        0.0000    0.0000    0.0000 H     1 UNL  0.0265
   6    H        0.0000    0.0000    0.0000 H     1 UNL  0.0265
@<TRIPOS>BOND
    1     2     3  1
    2     2     4  1
    3     2     5  1
    4     2     6  1
)DATA"_mol2;
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 2);
  }
  SECTION("_mol2 failure") {
    auto mol = R"DATA(@<TRIPOS>MOLECULE
UNK
6 5 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
   1   Na        0.0000    0.0000    0.0000 Na    1 UNL  1.0000
   2    C        0.0000    0.0000    0.0000 C.3   1 UNL -0.0305
   3    H        0.0000    0.0000    0.0000 H     1 UNL  0.0265
   4    H        0.0000    0.0000    0.0000 H     1 UNL  0.0265
   5    H        0.0000    0.0000    0.0000 H     1 UNL  0.0265
   6    H        0.0000    0.0000    0.0000 H     1 UNL  0.0265
@<TRIPOS>BOND
    1     2     3  1
    2     2     4  1
    3     2     5  1
    4     2     6  1
    5     2     1  1
)DATA"_mol2;
    REQUIRE(!mol);
  }
}

TEST_CASE("handling STBOX properties from v3k ctabs", "[feature][v3k]") {
  SECTION("atoms and bonds") {
    auto mol = R"CTAB(basic test
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0 STBOX=1
M  V30 2 C -5.6979 2.8332 0 0 STBOX=1
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2 STBOX=1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 4);
    int val;
    CHECK(mol->getAtomWithIdx(0)->getPropIfPresent(
        common_properties::molStereoCare, val));
    CHECK(val == 1);
    CHECK(mol->getAtomWithIdx(1)->getPropIfPresent(
        common_properties::molStereoCare, val));
    CHECK(val == 1);
    CHECK(!mol->getAtomWithIdx(2)->hasProp(common_properties::molStereoCare));
    REQUIRE(mol->getBondBetweenAtoms(0, 1));
    CHECK(mol->getBondBetweenAtoms(0, 1)->getPropIfPresent(
        common_properties::molStereoCare, val));
    CHECK(val == 1);
  }
  SECTION("bonds set if the atoms are also set 1") {
    auto mol = R"CTAB(basic test
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0 STBOX=1
M  V30 2 C -5.6979 2.8332 0 0 STBOX=1
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 4);
    int val;
    REQUIRE(mol->getBondBetweenAtoms(0, 1));
    CHECK(mol->getBondBetweenAtoms(0, 1)->getPropIfPresent(
        common_properties::molStereoCare, val));
    CHECK(val == 1);
  }
  SECTION("bonds set if the atoms are also set 2") {
    auto mol = R"CTAB(basic test
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0 STBOX=0
M  V30 2 C -5.6979 2.8332 0 0 STBOX=0
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 4);
    REQUIRE(mol->getBondBetweenAtoms(0, 1));
    CHECK(!mol->getBondBetweenAtoms(0, 1)->hasProp(
        common_properties::molStereoCare));
  }
  SECTION("bonds set if the atoms are also set 2") {
    auto mol = R"CTAB(basic test
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0 STBOX=1
M  V30 2 C -5.6979 2.8332 0 0 STBOX=0
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 4);
    REQUIRE(mol->getBondBetweenAtoms(0, 1));
    CHECK(!mol->getBondBetweenAtoms(0, 1)->hasProp(
        common_properties::molStereoCare));
  }
  SECTION("bonds set if the atoms are also set v2k") {
    auto mol = R"CTAB(basic test
  Mrv1810 01292015042D          

  4  3  0  0  0  0            999 V2000
   -3.7669    1.1053    0.0000 C   0  0  0  0  1  0  0  0  0  0  0  0
   -3.0524    1.5178    0.0000 C   0  0  0  0  1  0  0  0  0  0  0  0
   -2.3380    1.1053    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4814    1.5178    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  2  2  0  0  0  0
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 4);
    int val;
    REQUIRE(mol->getBondBetweenAtoms(0, 1));
    CHECK(mol->getBondBetweenAtoms(0, 1)->getPropIfPresent(
        common_properties::molStereoCare, val));
    CHECK(val == 1);
  }
}

TEST_CASE("github #2829: support MRV_IMPLICIT_H", "[feature][sgroups]") {
  SECTION("basics v2k") {
    auto mol = R"CTAB(
  Mrv1810 01302015262D          
 
  5  5  0  0  0  0            999 V2000
    1.1387   -1.3654    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6538   -0.6979    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1387   -0.0305    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.9233   -0.2854    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9233   -1.1104    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0  0  0  0
  2  3  4  0  0  0  0
  3  4  4  0  0  0  0
  1  5  4  0  0  0  0
  4  5  4  0  0  0  0
M  STY  1   1 DAT
M  SAL   1  1   3
M  SDT   1 MRV_IMPLICIT_H                                        
M  SDD   1     0.0000    0.0000    DR    ALL  0       0  
M  SED   1 IMPL_H1
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(MolToSmiles(*mol) == "c1cc[nH]c1");
  }
  SECTION("basics v3k") {
    auto mol = R"CTAB(
  Mrv1810 01302015452D          
 
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 5 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 2.1256 -2.5487 0 0
M  V30 2 C 1.2204 -1.3027 0 0
M  V30 3 N 2.1256 -0.0569 0 0
M  V30 4 C 3.5902 -0.5327 0 0
M  V30 5 C 3.5902 -2.0727 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 4 1 2
M  V30 2 4 2 3
M  V30 3 4 3 4
M  V30 4 4 1 5
M  V30 5 4 4 5
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 3) FIELDNAME=MRV_IMPLICIT_H -
M  V30 FIELDDISP="    0.0000    0.0000    DR    ALL  0       0" -
M  V30 FIELDDATA=IMPL_H1
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(MolToSmiles(*mol) == "c1cc[nH]c1");
    // we removed all the S groups:
    CHECK(getSubstanceGroups(*mol).empty());
  }
  SECTION("v3k two groups") {
    auto mol = R"CTAB(
  Mrv1810 01302016392D          
 
 12 14  0  0  0  0            999 V2000
    1.1387   -1.3654    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6538   -0.6979    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1387   -0.0305    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.9233   -0.2854    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9233   -1.1104    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6378    0.1271    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6378   -1.5229    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3522   -1.1104    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3522   -0.2854    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1369   -1.3653    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.1369   -0.0305    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6218   -0.6979    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0  0  0  0
  2  3  4  0  0  0  0
  3  4  4  0  0  0  0
  1  5  4  0  0  0  0
  4  5  4  0  0  0  0
  7  8  1  0  0  0  0
  6  9  1  0  0  0  0
  5  7  1  0  0  0  0
  6  4  1  0  0  0  0
 11 12  4  0  0  0  0
 10 12  4  0  0  0  0
  9 11  4  0  0  0  0
 10  8  4  0  0  0  0
  8  9  4  0  0  0  0
M  STY  2   1 DAT   2 DAT
M  SAL   1  1   3
M  SDT   1 MRV_IMPLICIT_H                                        
M  SDD   1     0.0000    0.0000    DR    ALL  0       0  
M  SED   1 IMPL_H1
M  SAL   2  1  10
M  SDT   2 MRV_IMPLICIT_H                                        
M  SDD   2     0.0000    0.0000    DR    ALL  0       0  
M  SED   2 IMPL_H1
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(MolToSmiles(*mol) == "c1cc2c([nH]1)Cc1cc[nH]c1C2");
    // we removed all the S groups:
    CHECK(getSubstanceGroups(*mol).empty());
  }
  SECTION("removal leaves other s groups intact") {
    auto mol = R"CTAB(
  Mrv1810 02022006062D          
 
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 2 0 0
M  V30 BEGIN ATOM
M  V30 1 C 2.1256 -2.5487 0 0
M  V30 2 C 1.2204 -1.3027 0 0
M  V30 3 N 2.1256 -0.0569 0 0
M  V30 4 C 3.5902 -0.5327 0 0
M  V30 5 C 3.5902 -2.0727 0 0
M  V30 6 C 4.8361 -2.9778 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 4 1 2
M  V30 2 4 2 3
M  V30 3 4 3 4
M  V30 4 4 1 5
M  V30 5 4 4 5
M  V30 6 1 5 6
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 6) FIELDNAME=some_data -
M  V30 FIELDDISP="    4.8361   -2.9778    DAU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=foo
M  V30 2 DAT 0 ATOMS=(1 3) FIELDNAME=MRV_IMPLICIT_H -
M  V30 FIELDDISP="    0.0000    0.0000    DR    ALL  0       0" -
M  V30 FIELDDATA=IMPL_H1
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(MolToSmiles(*mol) == "Cc1cc[nH]c1");
    CHECK(getSubstanceGroups(*mol).size() == 1);
  }
}

TEST_CASE("extra v3k mol file properties", "[ctab][v3k]") {
  SECTION("ATTCHPT") {
    auto mol = R"CTAB(
  Mrv2007 03132014352D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -16.625 7.1667 0 0 ATTCHPT=2
M  V30 2 C -15.2913 7.9367 0 0 ATTCHPT=1
M  V30 3 N -13.9576 7.1667 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(0)->getProp<int>(
              common_properties::molAttachPoint) == 2);
    CHECK(mol->getAtomWithIdx(1)->getProp<int>(
              common_properties::molAttachPoint) == 1);
    auto molb = MolToV3KMolBlock(*mol);
    CHECK(molb.find("ATTCHPT=1") != std::string::npos);
    CHECK(molb.find("ATTCHPT=2") != std::string::npos);
  }
  SECTION("others") {
    // this is not reasonable; just there to ensure that the reading/writing is
    // working
    auto mol = R"CTAB(really fake
  Mrv2007 03132015062D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -22.5833 11.0833 0 0 EXACHG=1
M  V30 2 C -21.2497 11.8533 0 0 INVRET=2
M  V30 3 C -23.917 11.8533 0 0 ATTCHORD=3
M  V30 4 C -25.2507 11.0833 0 0 CLASS=foo
M  V30 5 C -26.5844 11.8533 0 0 SEQID=4
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(0)->getProp<int>(
              common_properties::molRxnExactChange) == 1);
    CHECK(mol->getAtomWithIdx(1)->getProp<int>(
              common_properties::molInversionFlag) == 2);
    CHECK(mol->getAtomWithIdx(2)->getProp<int>(
              common_properties::molAttachOrder) == 3);
    CHECK(mol->getAtomWithIdx(3)->getProp<std::string>(
              common_properties::molAtomClass) == "foo");
    CHECK(mol->getAtomWithIdx(4)->getProp<int>(
              common_properties::molAtomSeqId) == 4);
    auto molb = MolToV3KMolBlock(*mol);
    CHECK(molb.find("EXACHG=1") != std::string::npos);
    CHECK(molb.find("INVRET=2") != std::string::npos);
    CHECK(molb.find("ATTCHORD=3") != std::string::npos);
    CHECK(molb.find("CLASS=foo") != std::string::npos);
    CHECK(molb.find("SEQID=4") != std::string::npos);
  }
  SECTION("SUBST") {
    auto mol = R"CTAB(test
  Mrv2007 03132018122D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -16.6248 7.1666 0 0 SUBST=3
M  V30 2 C -15.2911 7.9366 0 0 SUBST=-2
M  V30 3 N -13.9574 7.1666 0 0 SUBST=-1
M  V30 4 C -17.9585 7.9366 0 0 SUBST=6
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 1 4
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);

    auto smarts = MolToSmarts(*mol);
    CHECK(smarts == "[#6&D3](-[#6&D2]-[#7&D0])-[#6&D{6-}]");

    auto molb = MolToV3KMolBlock(*mol);
    CHECK(molb.find("SUBST=3") != std::string::npos);
    CHECK(molb.find("SUBST=-2") != std::string::npos);
    CHECK(molb.find("SUBST=-1") != std::string::npos);
  }
  SECTION("bond props") {
    auto mol = R"CTAB(bogus example
  Mrv2007 03132017102D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -28.125 8.2067 0 0
M  V30 2 C -26.7913 7.4367 0 0
M  V30 3 C -26.7913 5.8967 0 0
M  V30 4 C -28.125 5.1267 0 0
M  V30 5 N -29.4587 5.8967 0 0
M  V30 6 C -29.4587 7.4367 0 0
M  V30 7 * -27.2359 7.18 0 0
M  V30 8 R# -25.2354 8.335 0 0 RGROUPS=(1 1)
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6 
M  V30 6 2 1 6 RXCTR=1
M  V30 7 1 7 8 ENDPTS=(3 1 2 3) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(mol);

    auto molb = MolToV3KMolBlock(*mol);
    CHECK(molb.find("ENDPTS=(3 1 2 3) ATTACH=ANY") != std::string::npos);
    CHECK(molb.find("RXCTR=1") != std::string::npos);
  }
}

TEST_CASE(
    "Problems parsing SGroup abbreviations with multiple attachment points",
    "[bug][reader]") {
  std::string rdbase = getenv("RDBASE");
  SECTION("basics") {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/sgroup_ap_bug.mol";
    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    REQUIRE(mol);

    const auto &sgroups = getSubstanceGroups(*mol);
    CHECK(sgroups.size() == 3);
    CHECK(sgroups[0].hasProp("TYPE"));
    CHECK(sgroups[0].getProp<std::string>("TYPE") == "SUP");
    CHECK(sgroups[0].getAttachPoints().size() == 1);
    CHECK(sgroups[1].hasProp("TYPE"));
    CHECK(sgroups[1].getProp<std::string>("TYPE") == "SUP");
    CHECK(sgroups[1].getAttachPoints().size() == 1);
    CHECK(sgroups[2].hasProp("TYPE"));
    CHECK(sgroups[2].getProp<std::string>("TYPE") == "SUP");
    CHECK(sgroups[2].getAttachPoints().size() == 2);
  }
}

TEST_CASE(
    "github #3207: Attachment point info not being read from V2000 mol blocks",
    "[ctab][bug]") {
  SECTION("ATTCHPT") {
    auto mol = R"CTAB(
  Mrv1824 06092009122D          

  3  2  0  0  0  0            999 V2000
   -8.9061    3.8393    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1917    4.2518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4772    3.8393    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  APO  2   1   2   2   1
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(0)->getProp<int>(
              common_properties::molAttachPoint) == 2);
    CHECK(mol->getAtomWithIdx(1)->getProp<int>(
              common_properties::molAttachPoint) == 1);
    auto molb = MolToV3KMolBlock(*mol);
    CHECK(molb.find("ATTCHPT=1") != std::string::npos);
    CHECK(molb.find("ATTCHPT=2") != std::string::npos);
  }
  SECTION("Val=-1") {
    auto mol = R"CTAB(
  Mrv1824 06092009122D          

  3  2  0  0  0  0            999 V2000
   -8.9061    3.8393    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1917    4.2518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4772    3.8393    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  APO  2   1   3   2   1
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(0)->getProp<int>(
              common_properties::molAttachPoint) == -1);
    CHECK(mol->getAtomWithIdx(1)->getProp<int>(
              common_properties::molAttachPoint) == 1);
    auto molb = MolToV3KMolBlock(*mol);
    CHECK(molb.find("ATTCHPT=1") != std::string::npos);
    CHECK(molb.find("ATTCHPT=-1") != std::string::npos);
  }
}

TEST_CASE("XBHEAD and XBCORR causing parser failures", "[bug][reader]") {
  std::string rdbase = getenv("RDBASE");
  SECTION("basics") {
    std::string fName =
        rdbase +
        "/Code/GraphMol/FileParsers/sgroup_test_data/repeat_groups_query1.mol";
    std::unique_ptr<RWMol> mol(MolFileToMol(fName));
    REQUIRE(mol);
    const auto &sgroups = getSubstanceGroups(*mol);
    CHECK(sgroups.size() == 1);
    CHECK(sgroups[0].hasProp("TYPE"));
    CHECK(sgroups[0].getProp<std::string>("TYPE") == "SRU");
    CHECK(sgroups[0].hasProp("XBHEAD"));
    auto v = sgroups[0].getProp<std::vector<unsigned int>>("XBHEAD");
    CHECK(v.size() == 2);
    CHECK(v[0] == 5);
    CHECK(v[1] == 0);

    CHECK(sgroups[0].hasProp("XBCORR"));
    CHECK(sgroups[0].getProp<std::vector<unsigned int>>("XBCORR").size() == 4);

    auto mb = MolToV3KMolBlock(*mol);
    CHECK(mb.find("XBHEAD=(2 6 1)") != std::string::npos);
    CHECK(mb.find("XBCORR=(4 6 6 1 1)") != std::string::npos);
  }
}
TEST_CASE("LINKNODE information being ignored", "[ctab][bug]") {
  SECTION("v3000") {
    auto mol = R"CTAB(
  Mrv2007 06212005162D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -3.25 12.2683 0 0
M  V30 2 C -4.4959 11.3631 0 0
M  V30 3 C -4.02 9.8986 0 0
M  V30 4 C -2.48 9.8986 0 0
M  V30 5 C -2.0041 11.3631 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 1 5
M  V30 END BOND
M  V30 LINKNODE 1 3 2 1 2 1 5
M  V30 LINKNODE 1 4 2 4 3 4 5
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getProp<std::string>(common_properties::molFileLinkNodes) ==
          "1 3 2 1 2 1 5|1 4 2 4 3 4 5");
    auto molb = MolToV3KMolBlock(*mol);
    CHECK(molb.find("LINKNODE 1 3 2 1 2 1 5") != std::string::npos);
    CHECK(molb.find("LINKNODE 1 4 2 4 3 4 5") != std::string::npos);
  }
  SECTION("v2000") {
    auto mol = R"CTAB(
  Mrv2007 06222015182D          

  5  5  0  0  0  0            999 V2000
   -1.7411    6.5723    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4085    6.0874    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1536    5.3028    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3286    5.3028    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0736    6.0874    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  1  5  1  0  0  0  0
M  LIN  2   1   3   2   5   4   4   3   5
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getProp<std::string>(common_properties::molFileLinkNodes) ==
          "1 3 2 1 2 1 5|1 4 2 4 3 4 5");
    auto molb = MolToV3KMolBlock(*mol);
    CHECK(molb.find("LINKNODE 1 3 2 1 2 1 5") != std::string::npos);
    CHECK(molb.find("LINKNODE 1 4 2 4 3 4 5") != std::string::npos);
  }
}
TEST_CASE("more complex queries in CTAB parsers", "[ctab]") {
  SECTION("v3000") {
    auto mol = R"CTAB(*.*.*.*.*.*.*.* |$;Q_e;M_p;X_p;AH_p;QH_p;MH_p;XH_p$|
  manual  06272007272D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 A -3.2083 5.25 0 0
M  V30 2 Q -0.25 6 0 0
M  V30 3 M 4.5417 6.0417 0 0
M  V30 4 X 1.2917 4.2083 0 0
M  V30 5 AH -4.2083 5.25 0 0
M  V30 6 QH -1.25 6 0 0
M  V30 7 MH 3.5417 6.0417 0 0
M  V30 8 XH 0.2917 4.2083 0 0
M  V30 END ATOM
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    for (const auto atom : mol->atoms()) {
      REQUIRE(atom->hasQuery());
      CHECK(!atom->getQuery()->getTypeLabel().empty());
    }
    std::string pkl;
    MolPickler::pickleMol(*mol, pkl);
    ROMol cp(pkl);
    for (const auto atom : cp.atoms()) {
      REQUIRE(atom->hasQuery());
      CHECK(!atom->getQuery()->getTypeLabel().empty());
      CHECK(atom->getQuery()->getTypeLabel() ==
            mol->getAtomWithIdx(atom->getIdx())->getQuery()->getTypeLabel());
    }
    auto molb = MolToV3KMolBlock(*mol);
    CHECK(molb.find(" A ") != std::string::npos);
    CHECK(molb.find(" AH ") != std::string::npos);
    CHECK(molb.find(" Q ") != std::string::npos);
    CHECK(molb.find(" QH ") != std::string::npos);
    CHECK(molb.find(" M ") != std::string::npos);
    CHECK(molb.find(" MH ") != std::string::npos);
    CHECK(molb.find(" X ") != std::string::npos);
    CHECK(molb.find(" XH ") != std::string::npos);
  }
  SECTION("v2000") {
    auto mol = R"CTAB(*.*.*.*.*.*.*.* |$;Q_e;M_p;X_p;AH_p;QH_p;MH_p;XH_p$|
  manual  06272007272D          

  8  0  0  0  0  0            999 V2000
   -3.2083    5.2500    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2500    6.0000    0.0000 Q   0  0  0  0  0  0  0  0  0  0  0  0
    4.5417    6.0417    0.0000 M   0  0  0  0  0  0  0  0  0  0  0  0
    1.2917    4.2083    0.0000 X   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2083    5.2500    0.0000 AH  0  0  0  0  0  0  0  0  0  0  0  0
   -1.2500    6.0000    0.0000 QH  0  0  0  0  0  0  0  0  0  0  0  0
    3.5417    6.0417    0.0000 MH  0  0  0  0  0  0  0  0  0  0  0  0
    0.2917    4.2083    0.0000 XH  0  0  0  0  0  0  0  0  0  0  0  0
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    for (const auto atom : mol->atoms()) {
      REQUIRE(atom->hasQuery());
      CHECK(!atom->getQuery()->getTypeLabel().empty());
    }
    std::string pkl;
    MolPickler::pickleMol(*mol, pkl);
    ROMol cp(pkl);
    for (const auto atom : cp.atoms()) {
      REQUIRE(atom->hasQuery());
      CHECK(!atom->getQuery()->getTypeLabel().empty());
      CHECK(atom->getQuery()->getTypeLabel() ==
            mol->getAtomWithIdx(atom->getIdx())->getQuery()->getTypeLabel());
    }
    auto molb = MolToMolBlock(*mol);
    CHECK(molb.find(" A ") != std::string::npos);
    CHECK(molb.find(" AH ") != std::string::npos);
    CHECK(molb.find(" Q ") != std::string::npos);
    CHECK(molb.find(" QH ") != std::string::npos);
    CHECK(molb.find(" M ") != std::string::npos);
    CHECK(molb.find(" MH ") != std::string::npos);
    CHECK(molb.find(" X ") != std::string::npos);
    CHECK(molb.find(" XH ") != std::string::npos);
    /// SMARTS-based queries are not written for these:
    CHECK(molb.find("V    ") == std::string::npos);
  }
}

TEST_CASE("read metadata from PNG", "[reader][PNG]") {
  std::string rdbase = getenv("RDBASE");
  SECTION("basics") {
    std::string fname =
        rdbase + "/Code/GraphMol/FileParsers/test_data/colchicine.png";
    auto metadata = PNGFileToMetadata(fname);

    auto iter =
        std::find_if(metadata.begin(), metadata.end(),
                     [](const std::pair<std::string, std::string> &val) {
                       return boost::starts_with(val.first, PNGData::smilesTag);
                     });
    REQUIRE(iter != metadata.end());
    CHECK(
        iter->second ==
        "COc1cc2c(-c3ccc(OC)c(=O)cc3[C@@H](NC(C)=O)CC2)c(OC)c1OC "
        "|(6.46024,1.03002,;5.30621,1.98825,;3.89934,1.46795,;2.74531,2.42618,;"
        "1.33844,1.90588,;1.0856,0.427343,;-0.228013,-0.296833,;0.1857,-1."
        "73865,;-0.683614,-2.96106,;-2.18134,-3.04357,;-2.75685,-4.42878,;-4."
        "24422,-4.62298,;-3.17967,-1.92404,;-4.62149,-2.33775,;-2.92683,-0."
        "445502,;-1.61322,0.278673,;-2.02693,1.72049,;-3.50547,1.97333,;-4."
        "02577,3.3802,;-5.50431,3.63304,;-3.06754,4.53423,;-1.15762,2.9429,;0."
        "340111,3.02541,;2.23963,-0.530891,;1.98679,-2.00943,;3.14082,-2.96766,"
        ";3.6465,-0.0105878,;4.80053,-0.968822,;4.54769,-2.44736,)|");
  }
  SECTION("no metadata") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/colchicine.no_metadata.png";
    auto metadata = PNGFileToMetadata(fname);
    REQUIRE(metadata.empty());
  }
  SECTION("bad PNG") {
    std::string text = "NOT A PNG";
    REQUIRE_THROWS_AS(PNGStringToMetadata(text), FileParseException);
  }

  SECTION("truncated PNG") {
    std::string fname =
        rdbase + "/Code/GraphMol/FileParsers/test_data/colchicine.png";
    auto istr = std::ifstream(fname, std::ios_base::binary);
    istr.seekg(0, istr.end);
    auto sz = istr.tellg();
    istr.seekg(0, istr.beg);
    char *buff = new char[sz];
    istr.read(buff, sz);
    std::string data(buff, sz);
    delete[] buff;
    auto metadata = PNGStringToMetadata(data);
    REQUIRE(!metadata.empty());
    REQUIRE_THROWS_AS(PNGStringToMetadata(data.substr(1000)),
                      FileParseException);
  }
#ifdef RDK_USE_BOOST_IOSTREAMS
  SECTION("compressed metadata") {
    std::string fname =
        rdbase + "/Code/GraphMol/FileParsers/test_data/colchicine.mrv.png";
    auto metadata = PNGFileToMetadata(fname);
    auto iter =
        std::find_if(metadata.begin(), metadata.end(),
                     [](const std::pair<std::string, std::string> &val) {
                       return val.first == "molSource";
                     });
    REQUIRE(iter != metadata.end());
    CHECK(iter->second.find("<MChemicalStruct>") != std::string::npos);
  }
#endif
}

TEST_CASE("write metadata to PNG", "[writer][PNG]") {
  std::string rdbase = getenv("RDBASE");
  SECTION("basics") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/colchicine.no_metadata.png";
    std::vector<std::pair<std::string, std::string>> metadata;
    metadata.push_back(std::make_pair(
        PNGData::smilesTag,
        std::string(
            "COc1cc2c(-c3ccc(OC)c(=O)cc3[C@@H](NC(C)=O)CC2)c(OC)c1OC "
            "|(6.46024,1.03002,;5.30621,1.98825,;3.89934,1.46795,;2.74531,2."
            "42618,;1.33844,1.90588,;1.0856,0.427343,;-0.228013,-0.296833,;0."
            "1857,-1.73865,;-0.683614,-2.96106,;-2.18134,-3.04357,;-2.75685,-4."
            "42878,;-4.24422,-4.62298,;-3.17967,-1.92404,;-4.62149,-2.33775,;-"
            "2.92683,-0.445502,;-1.61322,0.278673,;-2.02693,1.72049,;-3.50547,"
            "1.97333,;-4.02577,3.3802,;-5.50431,3.63304,;-3.06754,4.53423,;-1."
            "15762,2.9429,;0.340111,3.02541,;2.23963,-0.530891,;1.98679,-2."
            "00943,;3.14082,-2.96766,;3.6465,-0.0105878,;4.80053,-0.968822,;4."
            "54769,-2.44736,)|")));
    auto pngData = addMetadataToPNGFile(fname, metadata);
    std::ofstream ofs("write_metadata.png");
    ofs.write(pngData.c_str(), pngData.size());
    ofs.flush();
    auto ometadata = PNGStringToMetadata(pngData);
    REQUIRE(ometadata.size() == metadata.size());
    for (unsigned int i = 0; i < ometadata.size(); ++i) {
      CHECK(ometadata[i].first == metadata[i].first);
      CHECK(ometadata[i].second == metadata[i].second);
    }
  }
}

TEST_CASE("read molecule from PNG", "[reader][PNG]") {
  std::string rdbase = getenv("RDBASE");
  SECTION("smiles") {
    std::string fname =
        rdbase + "/Code/GraphMol/FileParsers/test_data/colchicine.png";
    std::unique_ptr<ROMol> mol(PNGFileToMol(fname));
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 29);
    CHECK(mol->getNumConformers() == 1);
  }
  SECTION("mol") {
    std::string fname =
        rdbase + "/Code/GraphMol/FileParsers/test_data/colchicine.mol.png";
    std::unique_ptr<ROMol> mol(PNGFileToMol(fname));
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 29);
    CHECK(mol->getNumConformers() == 1);
  }
  SECTION("no metadata") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/colchicine.no_metadata.png";
    REQUIRE_THROWS_AS(PNGFileToMol(fname), FileParseException);
  }
}

TEST_CASE("write molecule to PNG", "[writer][PNG]") {
  std::string rdbase = getenv("RDBASE");
  SECTION("basics") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/colchicine.no_metadata.png";
    std::ifstream strm(fname, std::ios::in | std::ios::binary);
    auto colchicine =
        "COc1cc2c(c(OC)c1OC)-c1ccc(OC)c(=O)cc1[C@@H](NC(C)=O)CC2"_smiles;
    REQUIRE(colchicine);
    auto pngString = addMolToPNGStream(*colchicine, strm);
    // read it back out
    std::unique_ptr<ROMol> mol(PNGStringToMol(pngString));
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 29);
    CHECK(mol->getNumConformers() == 0);
  }
  SECTION("use SMILES") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/colchicine.no_metadata.png";
    std::ifstream strm(fname, std::ios::in | std::ios::binary);
    auto colchicine =
        "COc1cc2c(c(OC)c1OC)-c1ccc(OC)c(=O)cc1[C@@H](NC(C)=O)CC2"_smiles;
    REQUIRE(colchicine);
    bool includePkl = false;
    auto pngString = addMolToPNGStream(*colchicine, strm, includePkl);
    // read it back out
    std::unique_ptr<ROMol> mol(PNGStringToMol(pngString));
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 29);
    CHECK(mol->getNumConformers() == 0);
  }
  SECTION("use MOL") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/colchicine.no_metadata.png";
    std::ifstream strm(fname, std::ios::in | std::ios::binary);
    auto colchicine =
        "COc1cc2c(c(OC)c1OC)-c1ccc(OC)c(=O)cc1[C@@H](NC(C)=O)CC2"_smiles;
    REQUIRE(colchicine);
    bool includePkl = false;
    bool includeSmiles = false;
    bool includeMol = true;
    auto pngString = addMolToPNGStream(*colchicine, strm, includePkl,
                                       includeSmiles, includeMol);
    // read it back out
    std::unique_ptr<ROMol> mol(PNGStringToMol(pngString));
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 29);
    CHECK(mol->getNumConformers() == 1);
  }
}
TEST_CASE("multiple molecules in the PNG", "[writer][PNG]") {
  std::string rdbase = getenv("RDBASE");

  std::vector<std::string> smiles = {"c1ccccc1", "CCCOC", "c1ncccc1"};
  std::vector<std::unique_ptr<ROMol>> mols;
  for (const auto &smi : smiles) {
    mols.emplace_back(SmilesToMol(smi));
  }
  SECTION("pickles") {
    std::vector<std::pair<std::string, std::string>> metadata;
    for (const auto &mol : mols) {
      std::string pkl;
      MolPickler::pickleMol(*mol, pkl);
      metadata.push_back(std::make_pair(PNGData::pklTag, pkl));
    }
    // for the purposes of this test we'll add the metadata to an unrelated
    // PNG
    std::string fname =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/colchicine.no_metadata.png";
    std::ifstream strm(fname, std::ios::in | std::ios::binary);
    auto png = addMetadataToPNGStream(strm, metadata);
    std::stringstream pngstrm(png);
    auto molsRead = PNGStreamToMols(pngstrm);
    REQUIRE(molsRead.size() == mols.size());
    for (unsigned i = 0; i < molsRead.size(); ++i) {
      CHECK(MolToSmiles(*molsRead[i]) == MolToSmiles(*mols[i]));
    }
  }
  SECTION("SMILES") {
    std::vector<std::pair<std::string, std::string>> metadata;
    for (const auto &mol : mols) {
      std::string pkl = "BOGUS";
      // add bogus pickle data so we know that's not being read
      metadata.push_back(std::make_pair(PNGData::pklTag, pkl));
      metadata.push_back(std::make_pair(PNGData::smilesTag, MolToSmiles(*mol)));
    }
    // for the purposes of this test we'll add the metadata to an unrelated
    // PNG
    std::string fname =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/colchicine.no_metadata.png";
    std::ifstream strm(fname, std::ios::in | std::ios::binary);
    auto png = addMetadataToPNGStream(strm, metadata);
    std::stringstream pngstrm(png);
    auto molsRead = PNGStreamToMols(pngstrm, PNGData::smilesTag);
    REQUIRE(molsRead.size() == mols.size());
    for (unsigned i = 0; i < molsRead.size(); ++i) {
      CHECK(MolToSmiles(*molsRead[i]) == MolToSmiles(*mols[i]));
    }
  }
}

TEST_CASE("multiple molecules in the PNG, second example", "[writer][PNG]") {
  std::string rdbase = getenv("RDBASE");

  std::vector<std::string> smiles = {"c1ccccc1", "CCO", "CC(=O)O", "c1ccccn1"};
  std::vector<std::unique_ptr<ROMol>> mols;
  for (const auto &smi : smiles) {
    mols.emplace_back(SmilesToMol(smi));
  }
  SECTION("pickles") {
    std::string fname =
        rdbase + "/Code/GraphMol/FileParsers/test_data/multiple_mols.png";
    std::ifstream strm(fname, std::ios::in | std::ios::binary);
    auto molsRead = PNGStreamToMols(strm);
    REQUIRE(molsRead.size() == mols.size());
    for (unsigned i = 0; i < molsRead.size(); ++i) {
      CHECK(MolToSmiles(*molsRead[i]) == MolToSmiles(*mols[i]));
    }
  }
  SECTION("SMILES") {
    std::vector<std::pair<std::string, std::string>> metadata;
    for (const auto &mol : mols) {
      std::string pkl = "BOGUS";
      // add bogus pickle data so we know that's not being read
      metadata.push_back(std::make_pair(PNGData::pklTag, pkl));
      metadata.push_back(std::make_pair(PNGData::smilesTag, MolToSmiles(*mol)));
    }
    // for the purposes of this test we'll add the metadata to an unrelated
    // PNG
    std::string fname =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/colchicine.no_metadata.png";
    std::ifstream strm(fname, std::ios::in | std::ios::binary);
    auto png = addMetadataToPNGStream(strm, metadata);
    std::stringstream pngstrm(png);
    auto molsRead = PNGStreamToMols(pngstrm, PNGData::smilesTag);
    REQUIRE(molsRead.size() == mols.size());
    for (unsigned i = 0; i < molsRead.size(); ++i) {
      CHECK(MolToSmiles(*molsRead[i]) == MolToSmiles(*mols[i]));
    }
  }
}

TEST_CASE("github #3413: V3K mol blocks with no atoms fail to parse", "[bug]") {
  SECTION("basics") {
    auto m = R"CTAB(6065
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 0 0 0 0 0
M  V30 BEGIN ATOM
M  V30 END ATOM
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 0);
    CHECK(m->getNumBonds() == 0);
  }
}

TEST_CASE("github #3415: problem parsing SGroup data containing \" ", "[bug]") {
  SECTION("basics") {
    auto m = R"CTAB(
  Mrv2014 09172018222D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1.3337 2.31 0 0
M  V30 2 C 2.6674 1.54 0 0
M  V30 3 C 2.6674 -0 0 0
M  V30 4 C 1.3337 -0.77 0 0
M  V30 5 C 0 0 0 0
M  V30 6 C 0 1.54 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 1 6
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 1) FIELDNAME=Tempstruct FIELDINFO="""" -
M  V30 FIELDDISP="    2.1037    1.5400    DA    ALL  0       0" QUERYOP="""" -
M  V30 FIELDDATA=Foo1
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 6);
    CHECK(m->getNumBonds() == 6);
    auto sgs = getSubstanceGroups(*m);
    REQUIRE(sgs.size() == 1);
    CHECK(sgs[0].getProp<std::string>("TYPE") == "DAT");
    CHECK(sgs[0].getProp<std::string>("FIELDINFO") == "\"");
    CHECK(sgs[0].getProp<std::string>("QUERYOP") == "\"");
  }
  SECTION("empty string") {
    auto m = R"CTAB(
  Mrv2014 09172018222D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1.3337 2.31 0 0
M  V30 2 C 2.6674 1.54 0 0
M  V30 3 C 2.6674 -0 0 0
M  V30 4 C 1.3337 -0.77 0 0
M  V30 5 C 0 0 0 0
M  V30 6 C 0 1.54 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 1 6
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 1) FIELDNAME=Tempstruct FIELDINFO="" -
M  V30 FIELDDISP="    2.1037    1.5400    DA    ALL  0       0" QUERYOP="""" -
M  V30 FIELDDATA=Foo1
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 6);
    CHECK(m->getNumBonds() == 6);
    auto sgs = getSubstanceGroups(*m);
    REQUIRE(sgs.size() == 1);
    CHECK(sgs[0].getProp<std::string>("TYPE") == "DAT");
    CHECK(sgs[0].getProp<std::string>("FIELDINFO").empty());
    CHECK(sgs[0].getProp<std::string>("QUERYOP") == "\"");
  }
}

TEST_CASE("github #3597: Scientific notation in SDF V3000 files", "[bug]") {
  SECTION("basics") {
    auto m = R"CTAB(
  Mrv2020 11302014062D          

  2  1  0  0  0  0            999 V2000
   -2.8125    1.9196    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0980    2.3321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
)CTAB"_ctab;
    REQUIRE(m);
    m->getConformer().getAtomPos(0).z = 1e-6;
    m->getConformer().getAtomPos(1).z = 1e-4;
    auto mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("1e-06") == std::string::npos);
  }
  SECTION("toosmall") {
    auto m = R"CTAB(
  Mrv2020 11302014062D          

  2  1  0  0  0  0            999 V2000
   -2.8125    1.9196    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0980    2.3321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
)CTAB"_ctab;
    REQUIRE(m);
    m->getConformer().getAtomPos(0).z = 1e-17;
    m->getConformer().getAtomPos(1).z = 1e-4;
    auto mb = MolToV3KMolBlock(*m);
    // std::cerr<<mb<<std::endl;
    CHECK(mb.find("M  V30 1 C -2.812500 1.919600 0.000000 0") !=
          std::string::npos);
  }
}

TEST_CASE("github #3620: V3K mol block parser not saving the chiral flag",
          "[bug]") {
  SECTION("basics") {
    auto m = R"CTAB(
  Mrv2014 12082009582D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C -1.875 6.0417 0 0 CFG=2
M  V30 2 C -0.5413 6.8117 0 0
M  V30 3 F -3.2087 6.8117 0 0
M  V30 4 Cl -1.875 4.5017 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 3
M  V30 2 1 1 4
M  V30 3 1 1 2 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    unsigned int chiralFlag = 0;
    CHECK(
        m->getPropIfPresent(common_properties::_MolFileChiralFlag, chiralFlag));
    CHECK(chiralFlag == 1);
    auto mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("4 3 0 0 1") != std::string::npos);
  }
}

TEST_CASE("test bond flavors when writing PDBs", "[bug]") {
  SECTION("basics") {
    std::unique_ptr<RWMol> m{SequenceToMol("G")};
    REQUIRE(m);
    int confId = -1;
    {
      int flavor = 0;
      auto pdb = MolToPDBBlock(*m, confId, flavor);
      CHECK(pdb.find("CONECT    1    2\n") != std::string::npos);
      CHECK(pdb.find("CONECT    3    4    4    5\n") != std::string::npos);
    }
    {
      int flavor = 2;
      auto pdb = MolToPDBBlock(*m, confId, flavor);
      CHECK(pdb.find("CONECT    1    2\n") == std::string::npos);
      CHECK(pdb.find("CONECT    3    4    4\n") != std::string::npos);
    }
    {
      int flavor = 8;
      auto pdb = MolToPDBBlock(*m, confId, flavor);
      CHECK(pdb.find("CONECT    1    2\n") != std::string::npos);
      CHECK(pdb.find("CONECT    3    4    5\n") != std::string::npos);
    }
    {
      int flavor = 2 | 8;
      auto pdb = MolToPDBBlock(*m, confId, flavor);
      CHECK(pdb.find("CONECT") == std::string::npos);
    }
  }
}

TEST_CASE(
    "github #3599: Add explicit support for remaining CTAB query bond types",
    "[feature]") {
  SECTION("basics V3K") {
    auto m = R"CTAB(
  Mrv2014 11302009242D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 3.7917 -2.96 0 0
M  V30 2 C 2.458 -3.73 0 0
M  V30 3 O 2.458 -5.27 0 0
M  V30 4 C 3.7917 -6.04 0 0
M  V30 5 C 5.1253 -5.27 0 0
M  V30 6 C 5.1253 -3.73 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 4 5
M  V30 3 1 1 6
M  V30 4 5 1 2
M  V30 5 6 5 6
M  V30 6 7 3 4
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(3)->hasQuery());
    CHECK(m->getBondWithIdx(3)->getQuery()->getDescription() ==
          "SingleOrDoubleBond");
    REQUIRE(m->getBondWithIdx(4)->hasQuery());
    CHECK(m->getBondWithIdx(4)->getQuery()->getDescription() ==
          "SingleOrAromaticBond");
    REQUIRE(m->getBondWithIdx(5)->hasQuery());
    CHECK(m->getBondWithIdx(5)->getQuery()->getDescription() ==
          "DoubleOrAromaticBond");
    auto mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("M  V30 4 5 1 2") != std::string::npos);
    CHECK(mb.find("M  V30 5 6 5 6") != std::string::npos);
    CHECK(mb.find("M  V30 6 7 3 4") != std::string::npos);
    std::string pkl;
    MolPickler::pickleMol(*m, pkl);
    m.reset(new RWMol(pkl));
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(3)->hasQuery());
    CHECK(m->getBondWithIdx(3)->getQuery()->getDescription() ==
          "SingleOrDoubleBond");
    REQUIRE(m->getBondWithIdx(4)->hasQuery());
    CHECK(m->getBondWithIdx(4)->getQuery()->getDescription() ==
          "SingleOrAromaticBond");
    REQUIRE(m->getBondWithIdx(5)->hasQuery());
    CHECK(m->getBondWithIdx(5)->getQuery()->getDescription() ==
          "DoubleOrAromaticBond");

    auto smarts = MolToSmarts(*m);
    CHECK(smarts == "[#7]1-,=[#6]-[#8]=,:[#6]-[#6][#6]-1");
  }
  SECTION("basics V2K") {
    auto m = R"CTAB(
  Mrv2014 11302009442D          

  6  6  0  0  0  0            999 V2000
    2.0313   -1.5857    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.3168   -1.9982    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3168   -2.8232    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0313   -3.2357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7457   -2.8232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7457   -1.9982    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  4  5  1  0  0  0  0
  1  6  1  0  0  0  0
  1  2  5  0  0  0  0
  5  6  6  0  0  0  0
  3  4  7  0  0  0  0
M  END
)CTAB"_ctab;
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(3)->hasQuery());
    CHECK(m->getBondWithIdx(3)->getQuery()->getDescription() ==
          "SingleOrDoubleBond");
    REQUIRE(m->getBondWithIdx(4)->hasQuery());
    CHECK(m->getBondWithIdx(4)->getQuery()->getDescription() ==
          "SingleOrAromaticBond");
    REQUIRE(m->getBondWithIdx(5)->hasQuery());
    CHECK(m->getBondWithIdx(5)->getQuery()->getDescription() ==
          "DoubleOrAromaticBond");
    auto mb = MolToMolBlock(*m);
    CHECK(mb.find("  1  2  5") != std::string::npos);
    CHECK(mb.find("  5  6  6") != std::string::npos);
    CHECK(mb.find("  3  4  7") != std::string::npos);
    std::string pkl;
    MolPickler::pickleMol(*m, pkl);
    m.reset(new RWMol(pkl));
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(3)->hasQuery());
    CHECK(m->getBondWithIdx(3)->getQuery()->getDescription() ==
          "SingleOrDoubleBond");
    REQUIRE(m->getBondWithIdx(4)->hasQuery());
    CHECK(m->getBondWithIdx(4)->getQuery()->getDescription() ==
          "SingleOrAromaticBond");
    REQUIRE(m->getBondWithIdx(5)->hasQuery());
    CHECK(m->getBondWithIdx(5)->getQuery()->getDescription() ==
          "DoubleOrAromaticBond");
  }
}

TEST_CASE("supplier close methods") {
  std::string rdbase = getenv("RDBASE");

  SECTION("SDF") {
    std::string fname =
        rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
    {
      SDMolSupplier suppl(fname);
      std::unique_ptr<ROMol> mol{suppl.next()};
      REQUIRE(mol);
      suppl.close();
#if INVARIANT_EXCEPTION_METHOD
      REQUIRE_THROWS_AS(suppl.next(), Invar::Invariant);
#endif
    }
    {
      std::ifstream instr(fname);
      ForwardSDMolSupplier suppl(&instr, false);
      std::unique_ptr<ROMol> mol{suppl.next()};
      REQUIRE(mol);
      suppl.close();
#if INVARIANT_EXCEPTION_METHOD
      REQUIRE_THROWS_AS(suppl.next(), Invar::Invariant);
#endif
    }
  }
  SECTION("SMILES") {
    std::string fname =
        rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.2.csv";
    {
      SmilesMolSupplier suppl(fname, ",", 1, 0, true);
      std::unique_ptr<ROMol> mol{suppl.next()};
      REQUIRE(mol);
      suppl.close();
#if INVARIANT_EXCEPTION_METHOD
      REQUIRE_THROWS_AS(suppl.next(), Invar::Invariant);
#endif
    }
  }
  SECTION("TDT") {
    std::string fname =
        rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
    {
      TDTMolSupplier suppl(fname, "PN");
      std::unique_ptr<ROMol> mol{suppl.next()};
      REQUIRE(mol);
      suppl.close();
#if INVARIANT_EXCEPTION_METHOD
      REQUIRE_THROWS_AS(suppl.next(), Invar::Invariant);
#endif
    }
  }
#ifdef RDK_BUILD_MAEPARSER_SUPPORT
  SECTION("MAE") {
    std::string fname =
        rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.mae";
    {
      MaeMolSupplier suppl(fname);
      std::unique_ptr<ROMol> mol{suppl.next()};
      REQUIRE(mol);
      suppl.close();
#if INVARIANT_EXCEPTION_METHOD
      REQUIRE_THROWS_AS(suppl.next(), Invar::Invariant);
#endif
    }
  }
#endif
}

TEST_CASE(
    "github #3768: SubstanceGroup output doesn't properly quote double "
    "quotes") {
  SECTION("basics") {
    auto m = R"CTAB(
  Mrv2014 01292104542D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.3343 -0.7691 0 0
M  V30 2 C -1.333 0.7709 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 1) FIELDNAME=[DUP]Tempstruct FIELDINFO="""" -
M  V30 FIELDDISP="   -0.1770   -0.5034    DA    ALL  0       0" -
M  V30 QUERYOP="""""" FIELDDATA=Foo1
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    auto sgs = getSubstanceGroups(*m);
    REQUIRE(sgs.size() == 1);
    auto sg = sgs[0];
    CHECK(sg.getProp<std::string>("FIELDINFO") == "\"");
    CHECK(sg.getProp<std::string>("QUERYOP") == "\"\"");
    auto mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("FIELDINFO=\"\"\"\"") != std::string::npos);
    CHECK(mb.find("QUERYOP=\"\"\"\"\"") != std::string::npos);
  }
  SECTION("parens and quote not at beginning") {
    auto m = R"CTAB(
  Mrv2014 01292104542D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.3343 -0.7691 0 0
M  V30 2 C -1.333 0.7709 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 1) FIELDNAME=[DUP]Tempstruct FIELDINFO="foo""" -
M  V30 FIELDDISP="   -0.1770   -0.5034    DA    ALL  0       0" -
M  V30 QUERYOP="(bar)" FIELDDATA=Foo1
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    auto sgs = getSubstanceGroups(*m);
    REQUIRE(sgs.size() == 1);
    auto sg = sgs[0];
    CHECK(sg.getProp<std::string>("FIELDINFO") == "foo\"");
    CHECK(sg.getProp<std::string>("QUERYOP") == "(bar)");
    auto mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("FIELDINFO=\"foo\"\"\"") != std::string::npos);
    CHECK(mb.find("QUERYOP=\"(bar)\"") != std::string::npos);
  }
}

TEST_CASE("github #3216: WedgeMolBonds() should prefer degree-1 atoms") {
  SECTION("basics") {
    auto m = R"CTAB(
  Mrv2007 06082008522D          
 
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -17.7571 16.6703 0 0
M  V30 2 C -16.4234 17.4403 0 0
M  V30 3 C -15.0897 16.6703 0 0 CFG=1
M  V30 4 C -13.7561 17.4403 0 0 CFG=2
M  V30 5 Br -15.0897 15.1303 0 0
M  V30 6 C -12.4225 16.6703 0 0
M  V30 7 Cl -13.7561 18.9803 0 0
M  V30 8 C -11.0888 17.4403 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 5 CFG=1
M  V30 4 1 3 4
M  V30 5 1 4 6
M  V30 6 1 4 7 CFG=1
M  V30 7 1 6 8
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    WedgeMolBonds(*m, &m->getConformer());
    CHECK(m->getBondBetweenAtoms(2, 4)->getBondDir() != Bond::BondDir::NONE);
    CHECK(m->getBondBetweenAtoms(3, 6)->getBondDir() != Bond::BondDir::NONE);
    CHECK(m->getBondBetweenAtoms(2, 1)->getBondDir() == Bond::BondDir::NONE);
    CHECK(m->getBondBetweenAtoms(2, 3)->getBondDir() == Bond::BondDir::NONE);
    CHECK(m->getBondBetweenAtoms(3, 5)->getBondDir() == Bond::BondDir::NONE);
  }
}

TEST_CASE("Hydrogen bonds in CTABs") {
  SECTION("basics") {
    auto m = R"CTAB(
  Mrv2014 03022114422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -5.4583 -0.125 0 0
M  V30 2 C -4.1247 0.645 0 0
M  V30 3 C -2.791 -0.125 0 0
M  V30 4 C -1.4573 0.645 0 0
M  V30 5 O -2.791 -1.665 0 0
M  V30 6 C -6.792 0.645 0 0
M  V30 7 O -5.4583 -1.665 0 0
M  V30 8 H -4.1247 -2.435 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 2 3 5
M  V30 5 1 1 6
M  V30 6 1 1 7
M  V30 7 1 7 8
M  V30 8 10 5 8
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(4, 7));
    CHECK(m->getBondBetweenAtoms(4, 7)->getBondType() ==
          Bond::BondType::HYDROGEN);
    auto mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("V30 8 10 5 8") != std::string::npos);
    CHECK(MolToSmiles(*m) ==
          "CC1=O~[H]OC(C)C1");  // the SMILES writer still doesn't know what to
                                // do with it
  }
}

TEST_CASE("Support empty FIELDNAMES in SDT lines") {
  SECTION("basics") {
    auto m = R"CTAB(
  Mrv2014 03112117322D          

  6  6  0  0  0  0            999 V2000
   -1.8270   -1.5114    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2764   -0.8194    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1002   -0.8626    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4748   -1.5977    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0255   -2.2896    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2016   -2.2465    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  1  6  1  0  0  0  0
M  STY  1   1 DAT
M  SAL   1  6   1   2   3   4   5   6
M  SDT   1                                                       
M  SDD   1    -2.4921   -3.0466    DA    ALL  1       5  
M  SED   1 foo: 1234.6
M  END
)CTAB"_ctab;
    REQUIRE(m);
    auto sgs = getSubstanceGroups(*m);
    REQUIRE(sgs.size() == 1);
    {
      auto outctab = MolToMolBlock(*m);
      CHECK(outctab.find("1234.6") != std::string::npos);
      auto nm = MolBlockToMol(outctab);
      REQUIRE(nm);
      auto sgs = getSubstanceGroups(*nm);
      REQUIRE(sgs.size() == 1);
      delete nm;
    }
    {
      auto outctab = MolToV3KMolBlock(*m);
      CHECK(outctab.find("1234.6") != std::string::npos);
      auto nm = MolBlockToMol(outctab);
      REQUIRE(nm);
      auto sgs = getSubstanceGroups(*nm);
      REQUIRE(sgs.size() == 1);
      delete nm;
    }
  }
}

TEST_CASE("Support reading unambiguous short atom lines") {
  SECTION("basics") {
    std::string mb = R"CTAB(
  Mrv2014 03112117322D          

  2  1  0  0  0  0            999 V2000
   -1.8270   -1.5114    0.0000 C 
   -2.2764   -0.8194    0.0000 C
  1  2  1  0  0  0  0
M  END
)CTAB";
    // we fail when doing strict parsing
    REQUIRE_THROWS_AS(MolBlockToMol(mb), FileParseException);

    bool removeHs = true;
    bool sanitize = true;
    bool strictParsing = false;
    // but can read it with non-strict parsing
    std::unique_ptr<ROMol> m{
        MolBlockToMol(mb, sanitize, removeHs, strictParsing)};
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 2);
    CHECK(m->getAtomWithIdx(0)->getAtomicNum() == 6);
    CHECK(m->getAtomWithIdx(1)->getAtomicNum() == 6);
  }
  SECTION("too short") {
    std::string mb = R"CTAB(
  Mrv2014 03112117322D          

  2  1  0  0  0  0            999 V2000
   -1.8270   -1.5114    0.0000  
   -2.2764   -0.8194    0.0000 C
  1  2  1  0  0  0  0
M  END
)CTAB";
    // we fail when doing strict parsing
    REQUIRE_THROWS_AS(MolBlockToMol(mb), FileParseException);

    bool removeHs = true;
    bool sanitize = true;
    bool strictParsing = false;
    // fail even with non-strict parsing
    REQUIRE_THROWS_AS(MolBlockToMol(mb, sanitize, removeHs, strictParsing),
                      FileParseException);
  }
}

TEST_CASE("Github #4099: HCount field in v2000 mol blocks ignored") {
  SECTION("basics") {
    auto mol = R"CTAB(Test


  3  2  0  0  0  0  0  0  0  0999 V2000
    2.7500   -7.9167   -0.0000 C   0  0  0  3  0  0  0  0  0  0  0  0
    3.6160   -7.4167   -0.0000 N   0  0  0  1  0  0  0  0  0  0  0  0
    4.4821   -7.9167   -0.0000 C   0  0  0  2  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(MolToSmarts(*mol) == "[#6&h{2-}]-[#7&h0]-[#6&h{1-}]");
  }
}

TEST_CASE("Github #4131: HCOUNT from v3000 CTABS incorrectly interpreted") {
  SECTION("basics") {
    auto mol = R"CTAB(
  Mrv2108 05122108272D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 5.458 -2.2567 0 0 HCOUNT=2
M  V30 2 N 6.7916 -1.4867 0 0 HCOUNT=0
M  V30 3 C 8.1254 -2.2567 0 0 HCOUNT=-1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(MolToSmarts(*mol) == "[#6&h{2-}]-[#7]-[#6&h0]");
  }
}

TEST_CASE("sgroups and strict parsing") {
  SECTION("everything ok") {
    std::string ctab = R"CTAB(
  Mrv2108 06052107052D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 2 0 0
M  V30 BEGIN ATOM
M  V30 1 C -6.875 1.0417 0 0
M  V30 2 C -5.5413 1.8117 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 2) FIELDNAME=foo -
M  V30 FIELDDISP="   -5.5413    1.8117    DAU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=bar
M  V30 2 DAT 0 ATOMS=(1 1) FIELDNAME=foo -
M  V30 FIELDDISP="   -6.8750    1.0417    DAU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=baz
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB";
    std::unique_ptr<ROMol> m(MolBlockToMol(ctab));
    CHECK(m);
  }
  SECTION("SGroups totally missing") {
    std::string ctab = R"CTAB(
  Mrv2108 06052107052D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 2 0 0
M  V30 BEGIN ATOM
M  V30 1 C -6.875 1.0417 0 0
M  V30 2 C -5.5413 1.8117 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB";
    std::unique_ptr<ROMol> m;
    CHECK_THROWS_AS(m.reset(MolBlockToMol(ctab)), FileParseException);
    bool sanitize = true;
    bool removeHs = true;
    bool strictParsing = false;
    m.reset(MolBlockToMol(ctab, sanitize, removeHs, strictParsing));
    CHECK(m);
  }
  SECTION("one SGroup missing") {
    std::string ctab = R"CTAB(
  Mrv2108 06052107052D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 2 0 0
M  V30 BEGIN ATOM
M  V30 1 C -6.875 1.0417 0 0
M  V30 2 C -5.5413 1.8117 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 2 DAT 0 ATOMS=(1 1) FIELDNAME=foo -
M  V30 FIELDDISP="   -6.8750    1.0417    DAU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=baz
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB";
    std::unique_ptr<ROMol> m;
    // fails without an exception
    m.reset(MolBlockToMol(ctab));
    CHECK(!m);
    bool sanitize = true;
    bool removeHs = true;
    bool strictParsing = false;
    m.reset(MolBlockToMol(ctab, sanitize, removeHs, strictParsing));
    CHECK(m);
  }

  SECTION("END SGROUPS missing") {
    std::string ctab = R"CTAB(
  Mrv2108 06052107052D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 2 0 0
M  V30 BEGIN ATOM
M  V30 1 C -6.875 1.0417 0 0
M  V30 2 C -5.5413 1.8117 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 2) FIELDNAME=foo -
M  V30 FIELDDISP="   -5.5413    1.8117    DAU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=bar
M  V30 2 DAT 0 ATOMS=(1 1) FIELDNAME=foo -
M  V30 FIELDDISP="   -6.8750    1.0417    DAU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=baz
M  V30 END CTAB
M  END
)CTAB";
    std::unique_ptr<ROMol> m;
    CHECK_THROWS_AS(m.reset(MolBlockToMol(ctab)), FileParseException);
    bool sanitize = true;
    bool removeHs = true;
    bool strictParsing = false;
    m.reset(MolBlockToMol(ctab, sanitize, removeHs, strictParsing));
    CHECK(m);
  }
}

TEST_CASE("double bond stereo should not be set when the coords are all zero") {
  auto m = R"CTAB(
     RDKit          2D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
M  END)CTAB"_ctab;
  REQUIRE(m);
  REQUIRE(m->getBondBetweenAtoms(1, 2));
  CHECK(m->getBondBetweenAtoms(1, 2)->getBondDir() == Bond::EITHERDOUBLE);
}

TEST_CASE("Handle MRV_COORDINATE_BOND_TYPE data Substance Groups") {
  SECTION(
      "Convert SDF V2000 MRV_COORDINATE_BOND_TYPE data Substance Groups "
      "into coordinate bonds") {
    auto m = R"CTAB(
  Mrv2111 06302118332D          

  9  9  0  0  0  0            999 V2000
   -2.9465    0.7804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6609    0.3679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6609   -0.4572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9465   -0.8697    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2320   -0.4572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2320    0.3679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5175    0.7804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9465   -1.6947    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3754    0.7804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  6  7  8  0  0  0  0
  4  8  8  0  0  0  0
  2  9  8  0  0  0  0
M  STY  3   1 DAT   2 DAT   3 DAT
M  SAL   1  2   6   7
M  SDT   1 MRV_COORDINATE_BOND_TYPE                              
M  SDD   1     0.0000    0.0000    DR    ALL  0       0  
M  SED   1 7
M  SAL   2  2   4   8
M  SDT   2 MRV_COORDINATE_BOND_TYPE                              
M  SDD   2     0.0000    0.0000    DR    ALL  0       0  
M  SED   2 8
M  SAL   3  2   2   9
M  SDT   3 MRV_COORDINATE_BOND_TYPE                              
M  SDD   3     0.0000    0.0000    DR    ALL  0       0  
M  SED   3 9
M  END
)CTAB"_ctab;
    REQUIRE(m);

    std::vector<std::pair<unsigned, unsigned>> coordinate_bonds{
        {5, 6}, {3, 7}, {1, 8}};

    for (const auto &bond_atoms : coordinate_bonds) {
      auto bnd = m->getBondBetweenAtoms(bond_atoms.first, bond_atoms.second);
      REQUIRE(bnd);
      CHECK(bnd->getBondType() == Bond::BondType::DATIVE);
      CHECK(typeid(*bnd) == typeid(Bond));
    }

    CHECK(getSubstanceGroups(*m).empty());
  }

  SECTION(
      "GitHub Issue #4473: MRV_COORDINATE_BOND_TYPE SGroup may reference bond "
      "index, instead of atom") {
    // Same input as previous test, just shuffled the bonds and changed
    // the indexes in the SGroups

    auto m1 = R"CTAB(
  Mrv2111 06302118332D          

  9  9  0  0  0  0            999 V2000
   -2.9465    0.7804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6609    0.3679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6609   -0.4572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9465   -0.8697    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2320   -0.4572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2320    0.3679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5175    0.7804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9465   -1.6947    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3754    0.7804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  6  7  8  0  0  0  0
  4  8  8  0  0  0  0
  2  9  8  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
M  STY  3   1 DAT   2 DAT   3 DAT
M  SAL   1  2   6   7
M  SDT   1 MRV_COORDINATE_BOND_TYPE                              
M  SDD   1     0.0000    0.0000    DR    ALL  0       0  
M  SED   1 1
M  SAL   2  2   4   8
M  SDT   2 MRV_COORDINATE_BOND_TYPE                              
M  SDD   2     0.0000    0.0000    DR    ALL  0       0  
M  SED   2 2
M  SAL   3  2   2   9
M  SDT   3 MRV_COORDINATE_BOND_TYPE                              
M  SDD   3     0.0000    0.0000    DR    ALL  0       0  
M  SED   3 3
M  END
)CTAB"_ctab;

    // Same input, but changing the type of 2 of the bonds, and giving
    // a random value to the other SGroup to check that we fail
    auto m2 = R"CTAB(
  Mrv2111 06302118332D          

  9  9  0  0  0  0            999 V2000
   -2.9465    0.7804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6609    0.3679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6609   -0.4572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9465   -0.8697    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2320   -0.4572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2320    0.3679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5175    0.7804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9465   -1.6947    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3754    0.7804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  6  7  1  0  0  0  0
  4  8  1  0  0  0  0
  2  9  1  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
M  STY  3   1 DAT   2 DAT   3 DAT
M  SAL   1  2   6   7
M  SDT   1 MRV_COORDINATE_BOND_TYPE                              
M  SDD   1     0.0000    0.0000    DR    ALL  0       0  
M  SED   1 1
M  SAL   2  2   4   8
M  SDT   2 MRV_COORDINATE_BOND_TYPE                              
M  SDD   2     0.0000    0.0000    DR    ALL  0       0  
M  SED   2 2
M  SAL   3  2   2   9
M  SDT   3 MRV_COORDINATE_BOND_TYPE                              
M  SDD   3     0.0000    0.0000    DR    ALL  0       0  
M  SED   3 100
M  END
)CTAB"_ctab;

    std::vector<std::pair<unsigned, unsigned>> coordinate_bonds{
        {5, 6}, {3, 7}, {1, 8}};

    for (const auto &bond_atoms : coordinate_bonds) {
      auto bnd = m1->getBondBetweenAtoms(bond_atoms.first, bond_atoms.second);
      REQUIRE(bnd);
      CHECK(bnd->getBondType() == Bond::BondType::DATIVE);
      CHECK(typeid(*bnd) == typeid(Bond));
    }
    CHECK(getSubstanceGroups(*m1).empty());

    REQUIRE(m2);
    for (const auto &bond_atoms : coordinate_bonds) {
      auto bnd = m2->getBondBetweenAtoms(bond_atoms.first, bond_atoms.second);
      REQUIRE(bnd);
      CHECK(bnd->getBondType() != Bond::BondType::DATIVE);
    }
    CHECK(getSubstanceGroups(*m2).empty());
  }
}

TEST_CASE(
    "Github #4256: multiple ATTCHPT entries for one atom handled "
    "incorrectly") {
  SECTION("V3000") {
    std::string ctab = R"CTAB(
  Mrv2108 06172117542D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.8333 3.5 0 0
M  V30 2 C -3.4997 4.27 0 0 ATTCHPT=-1 ATTCHPT=3
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB";
    { REQUIRE_THROWS_AS(MolBlockToMol(ctab), FileParseException); }
    {
      bool sanitize = true;
      bool removeHs = true;
      bool strictParsing = false;
      std::unique_ptr<RWMol> m{
          MolBlockToMol(ctab, sanitize, removeHs, strictParsing)};
      REQUIRE(m);
      auto atom = m->getAtomWithIdx(1);
      REQUIRE(atom->hasProp(common_properties::molAttachPoint));
      REQUIRE(atom->getProp<int>(common_properties::molAttachPoint) == -1);
    }
  }
  SECTION("V2000 1") {  // Marvin doesn't actually do this, but might as well
                        // test for it anyway
    std::string ctab = R"CTAB(
  Mrv2108 06212115462D          

  2  1  0  0  0  0            999 V2000
   -2.5894    1.8751    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8749    2.2876    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  APO  2   2   3   2   2
M  END
)CTAB";
  }
  SECTION("V2000 2") {  // Marvin doesn't actually do this, but might as well
                        // test for it anyway
    std::string ctab = R"CTAB(
  Mrv2108 06212115482D          

  2  1  0  0  0  0            999 V2000
   -2.5894    1.8751    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8749    2.2876    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  APO  1   2   3
M  APO  1   2   2
M  END
)CTAB";
  }
}

TEST_CASE("Long lines in V3000 mol blocks") {
  SECTION("basics") {
    auto m = R"CTAB(query
  Mrv2108 07152116102D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 3.5417 -5.875 0 0
M  V30 2 C 4.8753 -5.105 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 2) -
M  V30 FIELDDISP="    4.8753   -5.1050    DA    ALL  0       0" -
M  V30 QUERYTYPE=SMARTSL QUERYOP== FIELDDATA=[#6;R]
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(getSubstanceGroups(*m).size() == 1);
    auto mb = MolToV3KMolBlock(*m);
    std::unique_ptr<RWMol> m2(MolBlockToMol(mb));
    REQUIRE(m2);
    CHECK(getSubstanceGroups(*m2).size() == 1);
  }
  SECTION("long data elements") {
    auto m = R"CTAB(query with bogus sgroups
  Mrv2108 07152116102D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 5 0 0
M  V30 BEGIN ATOM
M  V30 1 C 3.5417 -5.875 0 0
M  V30 2 C 4.8753 -5.105 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 2) -
M  V30 FIELDDISP="    0.0000    0.0000    DR    ALL  0       0" -
M  V30 QUERYTYPE=SMARTSR QUERYOP== -
M  V30 FIELDDATA="quite long piece of text that needs to be broken -
M  V30 across two lines"
M  V30 2 DAT 0 ATOMS=(1 1) -
M  V30 FIELDDISP="    0.0000    0.0000    DR    ALL  0       0" -
M  V30 QUERYTYPE=SMARTSR QUERYOP== -
M  V30 FIELDDATA="quite long piece of text that needs to be broken -
M  V30 across more than two lines because we really want to be sure -
M  V30 that we are doing this right"
M  V30 3 DAT 0 ATOMS=(1 1) -
M  V30 FIELDDISP="    0.0000    0.0000    DR    ALL  0       0" -
M  V30 QUERYTYPE=SMARTSR QUERYOP== -
M  V30 FIELDDATA="quite long piece of text that needs to be broken -
M  V30 across exactly two lines so that we can check the edge case -
M  V30 11111111111111111111"
M  V30 4 DAT 0 ATOMS=(1 1) -
M  V30 FIELDDISP="    0.0000    0.0000    DR    ALL  0       0" -
M  V30 QUERYTYPE=SMARTSR QUERYOP== -
M  V30 FIELDDATA="quite long piece of text that needs to be broken -
M  V30 across more than two lines because we really want to be sure -
M  V30 that we are doing this right" SEQID=1
M  V30 5 DAT 0 ATOMS=(1 1) -
M  V30 FIELDDISP="    0.0000    0.0000    DR    ALL  0       0" -
M  V30 QUERYTYPE=SMARTSR QUERYOP== -
M  V30 FIELDDATA="quite long piece of text that needs to be broken -
M  V30 across exactly two lines so that we can check the edge case -
M  V30 11111111111111111111" SEQID=2
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(getSubstanceGroups(*m).size() == 5);
    auto mb = MolToV3KMolBlock(*m);
    std::unique_ptr<RWMol> m2(MolBlockToMol(mb));
    REQUIRE(m2);
    CHECK(getSubstanceGroups(*m2).size() == 5);
  }

  SECTION(
      "GitHub Issue #4471: SDF SGroups may be missing the final space in the "
      "\"M V30 \" prefix") {
    auto m = R"CTAB(bogus mol with unspaced SGroup field
  Mrv2114 09022123382D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.7917 0 0 0
M  V30 END ATOM
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 1) FIELDNAME=Data -
M  V30 FIELDDISP="    0.0000    0.0000    DRU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 -
M  V30 FIELDDATA=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA-
M  V30 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA-
M  V30 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA-
M  V30 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(getSubstanceGroups(*m).size() == 1);
    auto mb = MolToV3KMolBlock(*m);
    std::unique_ptr<RWMol> m2(MolBlockToMol(mb));
    REQUIRE(m2);
    CHECK(getSubstanceGroups(*m2).size() == 1);
  }

  SECTION(
      "GitHub Issue #4477: Same SDF SGroup lines may be written multiple "
      "times") {
    auto m = R"CTAB(
  Mrv2014 03112117322D

  6  6  0  0  0  0            999 V2000
   -1.8270   -1.5114    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2764   -0.8194    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1002   -0.8626    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4748   -1.5977    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0255   -2.2896    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2016   -2.2465    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  1  6  1  0  0  0  0
M  STY  1   1 DAT
M  SAL   1  6   1   2   3   4   5   6
M  SDT   1
M  SDD   1    -2.4921   -3.0466    DA  123456789012345  ALL  1       5
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(getSubstanceGroups(*m).size() == 1);
    auto mb = MolToV3KMolBlock(*m);

    auto pos = 0u;
    auto count = 0u;
    std::string target{"FIELDDISP"};
    while (pos < mb.size()) {
      pos = mb.find(target, pos);
      if (pos < mb.size()) {
        pos += target.size();
        ++count;
      }
    }

    CHECK(count == 1);
  }
}

TEST_CASE("github #4345: non-stereo bonds written with unspecified parity") {
  SECTION("basics") {
    auto m = "CC=C(F)F"_smiles;
    REQUIRE(m);
    auto mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("CFG=2") == std::string::npos);
    mb = MolToMolBlock(*m);
    CHECK(mb.find("  2  3  2  3") == std::string::npos);
  }
  SECTION("possible chirality") {
    auto m = "CC=C(O)F"_smiles;
    REQUIRE(m);
    auto mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("CFG=2") != std::string::npos);
    mb = MolToMolBlock(*m);
    CHECK(mb.find("  2  3  2  3") != std::string::npos);
  }
  SECTION("terminal") {
    auto m = "CC=C"_smiles;
    REQUIRE(m);
    auto mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("CFG=2") == std::string::npos);
    mb = MolToMolBlock(*m);
    CHECK(mb.find("  2  3  2  3") == std::string::npos);
  }
  SECTION("nitrogen") {
    auto m = "CC(C)=NF"_smiles;
    REQUIRE(m);
    auto mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("CFG=2") == std::string::npos);
    mb = MolToMolBlock(*m);
    CHECK(mb.find("  3  4  2  3") == std::string::npos);
  }
  SECTION("nitrogen with") {
    auto m = "CC=NF"_smiles;
    REQUIRE(m);
    auto mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("CFG=2") != std::string::npos);
    mb = MolToMolBlock(*m);
    CHECK(mb.find("  2  3  2  3") != std::string::npos);
  }
  SECTION("direction explicitly set should be ignored") {
    auto m = "CC=C(F)F"_smiles;
    REQUIRE(m);
    m->getBondWithIdx(0)->setBondDir(Bond::BondDir::ENDUPRIGHT);
    m->getBondWithIdx(2)->setBondDir(Bond::BondDir::ENDUPRIGHT);
    auto mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("CFG=2") == std::string::npos);
    mb = MolToMolBlock(*m);
    CHECK(mb.find("  2  3  2  3") == std::string::npos);
  }
}

TEST_CASE(
    "github #4476: Additional SDT properties not decoded if FIELDNAME is "
    "empty") {
  SECTION("basics") {
    auto m = R"CTAB(query
  Mrv2102 09032106302D          

  2  1  0  0  0  0            999 V2000
   -0.4464    2.4334    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2680    2.8459    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  STY  1   1 DAT
M  SAL   1  1   2
M  SDT   1                                                     PQ=
M  SDD   1     0.0000    0.0000    DR    ALL  0       0  
M  SED   1 [#6;R]
M  END
)CTAB"_ctab;
    REQUIRE(m);
    REQUIRE(getSubstanceGroups(*m).size() == 1);
    const auto sg = getSubstanceGroups(*m)[0];
    CHECK(sg.hasProp("QUERYTYPE"));
    CHECK(sg.getProp<std::string>("QUERYTYPE") == "PQ");
    CHECK(sg.hasProp("QUERYOP"));
    CHECK(sg.getProp<std::string>("QUERYOP") == "=");
  }
}

TEST_CASE("github #4468: decode SMARTS in SGroups") {
  SECTION("parsing v3000") {
    auto m = R"CTAB(query
  Mrv2108 07152116012D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.8333 4.5421 0 0
M  V30 2 C 0.5003 5.3121 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 2) -
M  V30 FIELDDISP="    0.0000    0.0000    DR    ALL  0       0" -
M  V30 QUERYTYPE=SMARTSQ QUERYOP== FIELDDATA=[#6;R]
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->hasQuery());
    CHECK(SmartsWrite::GetAtomSmarts(
              static_cast<QueryAtom *>(m->getAtomWithIdx(1))) == "[#6&R]");
    CHECK(getSubstanceGroups(*m).empty());
  }
  SECTION("ensure bad SMARTS don't break things") {
    auto m = R"CTAB(query
  Mrv2108 07152116012D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.8333 4.5421 0 0
M  V30 2 C 0.5003 5.3121 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 2) -
M  V30 FIELDDISP="    0.0000    0.0000    DR    ALL  0       0" -
M  V30 QUERYTYPE=SMARTSQ QUERYOP== FIELDDATA=[#6;R
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    REQUIRE(!m->getAtomWithIdx(1)->hasQuery());
    CHECK(getSubstanceGroups(*m).empty());
  }
  SECTION("empty SMARTS") {
    auto m = R"CTAB(query
  Mrv2108 07152116012D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.8333 4.5421 0 0
M  V30 2 C 0.5003 5.3121 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 2) -
M  V30 FIELDDISP="    0.0000    0.0000    DR    ALL  0       0" -
M  V30 QUERYTYPE=SMARTSQ QUERYOP== FIELDDATA=
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    REQUIRE(!m->getAtomWithIdx(1)->hasQuery());
    CHECK(getSubstanceGroups(*m).empty());
  }
  SECTION("bad operator") {
    auto m = R"CTAB(query
  Mrv2108 07152116012D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.8333 4.5421 0 0
M  V30 2 C 0.5003 5.3121 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 2) -
M  V30 FIELDDISP="    0.0000    0.0000    DR    ALL  0       0" -
M  V30 QUERYTYPE=SMARTSQ QUERYOP=> FIELDDATA=[#6;R]
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    REQUIRE(!m->getAtomWithIdx(1)->hasQuery());
    CHECK(getSubstanceGroups(*m).empty());
  }
  SECTION("SMARTS with multiple atoms become recursive") {
    auto m = R"CTAB(query
  Mrv2108 07152116012D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.8333 4.5421 0 0
M  V30 2 C 0.5003 5.3121 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 2) -
M  V30 FIELDDISP="    0.0000    0.0000    DR    ALL  0       0" -
M  V30 QUERYTYPE=SMARTSQ QUERYOP== FIELDDATA=[#6;R]-[#8]
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->hasQuery());
    CHECK(SmartsWrite::GetAtomSmarts(static_cast<QueryAtom *>(
              m->getAtomWithIdx(1))) == "[$([#6&R]-[#8])]");
    CHECK(getSubstanceGroups(*m).empty());
  }
  SECTION("parsing v3000, v2000 compatibility version") {
    auto m = R"CTAB(query
  Mrv2108 07152116012D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.8333 4.5421 0 0
M  V30 2 C 0.5003 5.3121 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 2) -
M  V30 FIELDDISP="    0.0000    0.0000    DR    ALL  0       0" -
M  V30 QUERYTYPE=SQ QUERYOP== FIELDDATA=[#6;R]
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->hasQuery());
    CHECK(SmartsWrite::GetAtomSmarts(
              static_cast<QueryAtom *>(m->getAtomWithIdx(1))) == "[#6&R]");
    CHECK(getSubstanceGroups(*m).empty());
  }
  SECTION("parsing v2000") {
    auto m = R"CTAB(query
  Mrv2102 09032106302D          

  2  1  0  0  0  0            999 V2000
   -0.4464    2.4334    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2680    2.8459    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  STY  1   1 DAT
M  SAL   1  1   2
M  SDT   1                                                     SQ=
M  SDD   1     0.0000    0.0000    DR    ALL  0       0  
M  SED   1 [#6;R]
M  END
)CTAB"_ctab;
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->hasQuery());
    CHECK(SmartsWrite::GetAtomSmarts(
              static_cast<QueryAtom *>(m->getAtomWithIdx(1))) == "[#6&R]");
    CHECK(getSubstanceGroups(*m).empty());
  }
}

TEST_CASE("Github #4561: failure to parse CTAB with LINKNODE and SGROUP") {
  SECTION("BASICS") {
    auto mol1 = R"CTAB(
  Mrv2108 09252106182D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 7 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 12.719 -9.1518 0 0
M  V30 2 O 14.0458 -9.9326 0 0
M  V30 3 * 15.3857 -9.1735 0 0
M  V30 4 * 11.379 -9.9108 0 0
M  V30 5 C 12.7317 -7.6118 0 0
M  V30 6 C 12.2558 -6.1472 0 0
M  V30 7 C 13.7622 -6.4674 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 5 1 7 6
M  V30 6 1 6 5
M  V30 7 1 5 7
M  V30 END BOND
M  V30 LINKNODE 1 2 2 7 5 7 6
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(5 2 1 5 7 6) XBONDS=(2 2 3) BRKXYZ=(9 12.044 -10.2974 0 -
M  V30 12.044 -8.7593 0 0 0 0) BRKXYZ=(9 14.7161 -8.7854 0 14.7161 -10.3235 0 -
M  V30 0 0 0) CONNECT=HT LABEL=n
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol1);
    CHECK(getSubstanceGroups(*mol1).size() == 1);
    CHECK(mol1->hasProp(common_properties::molFileLinkNodes));
  }
}

TEST_CASE(
    "Github #4785: MDL query with aromatic bond sets aromatic flag on atoms "
    "even though they are not in an aromatic ring") {
  SECTION("benzene") {
    auto mol = R"CTAB(
  Mrv2108 12102110572D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 128.125 -103.585 0 0
M  V30 2 C 126.7913 -104.355 0 0
M  V30 3 C 126.7913 -105.895 0 0
M  V30 4 C 128.125 -106.665 0 0
M  V30 5 C 129.4587 -105.895 0 0
M  V30 6 C 129.4587 -104.355 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 4 1 2
M  V30 2 4 2 3
M  V30 3 4 3 4
M  V30 4 4 4 5
M  V30 5 4 5 6
M  V30 6 4 1 6
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(0)->getIsAromatic());
    CHECK(mol->getBondWithIdx(0)->getIsAromatic());
  }
  SECTION("non-kekulizeable") {
    auto mol = R"CTAB(
  Mrv2108 12102110572D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 128.125 -103.585 0 0
M  V30 2 C 126.7913 -104.355 0 0
M  V30 3 C 126.7913 -105.895 0 0
M  V30 4 C 128.125 -106.665 0 0
M  V30 5 C 129.4587 -105.895 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 4 1 2
M  V30 2 4 2 3
M  V30 3 4 3 4
M  V30 4 4 4 5
M  V30 5 4 5 1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(!mol);
  }
  SECTION("as reported1") {
    auto mol = R"CTAB(
  MJ201100                      

  2  1  0  0  0  0  0  0  0  0999 V2000
   -0.3538    0.6163    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0668    0.2012    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0  0  1  0
M  MRV SMA   1 [#6;a;a]
M  MRV SMA   2 [#6;a;a]
M  END)CTAB"_ctab;
    REQUIRE(mol);
  }
  SECTION("as reported2") {
    auto mol = R"CTAB(
  MJ201100                      

  2  1  0  0  0  0  0  0  0  0999 V2000
   -0.3538    0.6163    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0668    0.2012    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0  0  1  0
M  END)CTAB"_ctab;
    REQUIRE(mol);
  }
  SECTION("as reported3") {
    auto mol = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 A -0.353800 0.616300 0.000000 0
M  V30 2 A -1.066800 0.201200 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 4 1 2 TOPO=1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
  }
}
TEST_CASE("checking array bounds") {
  SECTION("XBONDS") {
    auto mb = R"CTAB(
  Mrv2108 01202214292D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 1 0 0
M  V30 BEGIN ATOM
M  V30 1 * -6.6667 7.5833 0 0
M  V30 2 C -5.333 8.3533 0 0
M  V30 3 C -3.9993 7.5833 0 0
M  V30 4 * -2.6656 8.3533 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(2 2 3) XBONDS=(20 1 3) BRKXYZ=(9 -3.9121 8.7006 0 -
M  V30 -2.9881 7.1002 0 0 0 0) BRKXYZ=(9 -5.4201 7.2361 0 -6.3441 8.8365 0 0 -
M  V30 0 0) CONNECT=HT LABEL=n
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB";
    std::unique_ptr<RWMol> mol;
    REQUIRE_THROWS_AS(mol.reset(MolBlockToMol(mb)), FileParseException);
  }
  SECTION("ATOMS") {
    auto mb = R"CTAB(
  Mrv2108 01202214292D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 1 0 0
M  V30 BEGIN ATOM
M  V30 1 * -6.6667 7.5833 0 0
M  V30 2 C -5.333 8.3533 0 0
M  V30 3 C -3.9993 7.5833 0 0
M  V30 4 * -2.6656 8.3533 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(12 2 3) XBONDS=(2 1 3) BRKXYZ=(9 -3.9121 8.7006 0 -
M  V30 -2.9881 7.1002 0 0 0 0) BRKXYZ=(9 -5.4201 7.2361 0 -6.3441 8.8365 0 0 -
M  V30 0 0) CONNECT=HT LABEL=n
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB";
    std::unique_ptr<RWMol> mol;
    REQUIRE_THROWS_AS(mol.reset(MolBlockToMol(mb)), FileParseException);
  }
}

TEST_CASE("Github #5108: Wiggly bonds don't override wedged bonds") {
  SECTION("as reported") {
    auto m = R"CTAB(
  Mrv2102 03212207042D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C 1.54 -1.54 0 0
M  V30 2 C 1.54 0 0 0
M  V30 3 O 1.54 1.54 0 0
M  V30 4 F 3.08 -0 0 0
M  V30 5 Cl 0 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 2 3 CFG=2
M  V30 3 1 2 4 CFG=1
M  V30 4 1 2 5
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
  SECTION("as reported, bond ordering changed") {
    auto m = R"CTAB(
  Mrv2102 03212207042D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C 1.54 -1.54 0 0
M  V30 2 C 1.54 0 0 0
M  V30 3 O 1.54 1.54 0 0
M  V30 4 F 3.08 -0 0 0
M  V30 5 Cl 0 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 2 4 CFG=1
M  V30 3 1 2 3 CFG=2
M  V30 4 1 2 5
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
  SECTION("assignChiralTypesFromBondDirs details") {
    auto m = R"CTAB(
  Mrv2102 03212207042D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C 1.54 -1.54 0 0
M  V30 2 C 1.54 0 0 0
M  V30 3 O 1.54 1.54 0 0
M  V30 4 F 3.08 -0 0 0
M  V30 5 Cl 0 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 2 4 CFG=1
M  V30 3 1 2 3
M  V30 4 1 2 5
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    m->getBondBetweenAtoms(1, 2)->setBondDir(Bond::BondDir::UNKNOWN);
    bool replaceExistingTags = false;
    MolOps::assignChiralTypesFromBondDirs(*m, -1, replaceExistingTags);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    replaceExistingTags = true;
    MolOps::assignChiralTypesFromBondDirs(*m, -1, replaceExistingTags);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
}

TEST_CASE(
    "Github #5152: presence of exocyclic query bonds in CTAB prevents "
    "aromaticity perception") {
  SECTION("as reported") {
    auto m = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C -2.229300 0.915100 0.000000 0
M  V30 2 C -3.562800 0.145100 0.000000 0
M  V30 3 C -3.562800 -1.395100 0.000000 0
M  V30 4 C -2.229300 -2.165100 0.000000 0
M  V30 5 C -0.895500 -1.395100 0.000000 0
M  V30 6 C -0.895500 0.145100 0.000000 0
M  V30 7 A 0.438100 0.915100 0.000000 0
M  V30 8 A 0.438100 -2.165100 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 6 6 7
M  V30 8 6 5 8
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(0)->getIsAromatic());
    CHECK(m->getBondWithIdx(0)->getIsAromatic());
  }
  SECTION("more detailed") {
    std::string molb = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C -2.229300 0.915100 0.000000 0
M  V30 2 C -3.562800 0.145100 0.000000 0
M  V30 3 C -3.562800 -1.395100 0.000000 0
M  V30 4 C -2.229300 -2.165100 0.000000 0
M  V30 5 C -0.895500 -1.395100 0.000000 0
M  V30 6 C -0.895500 0.145100 0.000000 0
M  V30 7 A 0.438100 0.915100 0.000000 0
M  V30 8 A 0.438100 -2.165100 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 6 6 7
M  V30 8 6 5 8
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB";
    std::string ptrn = "7 6 6 7";
    std::vector<std::string> alternatives = {
        "7 5 6 7",  // S/D
        "7 7 6 7",  // D/A
        "7 8 6 7",  // any
    };
    auto pos = molb.find(ptrn);
    REQUIRE(pos != std::string::npos);
    for (auto alternative : alternatives) {
      auto mb2 = molb.replace(pos, ptrn.size(), alternative);
      std::unique_ptr<RWMol> m(MolBlockToMol(mb2));
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(0)->getIsAromatic());
      CHECK(m->getBondWithIdx(0)->getIsAromatic());
    }
  }
}

TEST_CASE(
    "Github #5165: issue with V3000 SD files containing enhanced "
    "stereochemistry information") {
  SECTION("as reported") {
    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase +
                        "/Code/GraphMol/FileParsers/test_data/"
                        "mol_with_enhanced_stereo_2_And_groups.sdf";
    SDMolSupplier suppl(fName);
    std::unique_ptr<ROMol> mol{suppl.next()};
    REQUIRE(mol);
    auto groups = mol->getStereoGroups();
    REQUIRE(groups.size() == 2);
    CHECK(groups[0].getGroupType() == RDKit::StereoGroupType::STEREO_AND);
    CHECK(groups[1].getGroupType() == RDKit::StereoGroupType::STEREO_AND);
  }

  SECTION("as reported, less whitespace") {
    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/m_with_enh_stereo.sdf";
    SDMolSupplier suppl(fName);
    std::unique_ptr<ROMol> mol{suppl.next()};
    REQUIRE(mol);
    auto groups = mol->getStereoGroups();
    REQUIRE(groups.size() == 2);
    CHECK(groups[0].getGroupType() == RDKit::StereoGroupType::STEREO_AND);
    CHECK(groups[1].getGroupType() == RDKit::StereoGroupType::STEREO_AND);
  }
}

TEST_CASE("POL atoms in CTABS") {
  SECTION("V3000") {
    auto mol = R"CTAB(
  Mrv2102 05042219282D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Pol -6.25 3.375 0 0
M  V30 2 C -4.9163 4.145 0 0
M  V30 3 C -3.5826 3.375 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    std::string val;
    CHECK(mol->getAtomWithIdx(0)->getPropIfPresent(
        common_properties::dummyLabel, val));
    CHECK(val == "Pol");

    auto mb = MolToV3KMolBlock(*mol);
    CHECK(mb.find("1 Pol") != std::string::npos);
    mol->clearConformers();
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "*CC |$Pol_p;;$|");
  }
  SECTION("V2000") {
    auto mol = R"CTAB(
  Mrv2102 05042219412D          

  3  2  0  0  0  0            999 V2000
   -3.3482    1.8080    0.0000 Mod 0  0  0  0  0  0  0  0  0  0  0  0
   -2.6337    2.2205    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9193    1.8080    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    std::string val;
    CHECK(mol->getAtomWithIdx(0)->getPropIfPresent(
        common_properties::dummyLabel, val));
    CHECK(val == "Mod");
    auto mb = MolToMolBlock(*mol);
    CHECK(mb.find("0 Mod 0") != std::string::npos);
    mol->clearConformers();
    auto smi = MolToCXSmiles(*mol);
    CHECK(smi == "*CC |$Mod_p;;$|");
  }
}

TEST_CASE("PDB ACE caps bond order") {
  auto mol = R"DATA(HEADER TEST
ATOM      1  H1  ACE     1      25.950  25.179  24.582  1.00  0.00           H
ATOM      2  CH3 ACE     1      25.986  26.176  24.145  1.00  0.00           C
ATOM      3  H2  ACE     1      25.332  26.843  24.703  1.00  0.00           H
ATOM      4  H3  ACE     1      25.673  26.131  23.104  1.00  0.00           H
ATOM      5  C   ACE     1      27.405  26.691  24.218  1.00  0.00           C
ATOM      6  O   ACE     1      28.285  25.999  24.713  1.00  0.00           O
ATOM      7  N   ALA     2      27.621  27.909  23.728  1.00  0.00           N
ATOM      8  H   ALA     2      26.838  28.435  23.370  1.00  0.00           H
ATOM      9  CA  ALA     2      28.916  28.589  23.730  1.00  0.00           C
ATOM     10  HA  ALA     2      29.471  28.288  24.620  1.00  0.00           H
ATOM     11  CB  ALA     2      29.710  28.153  22.489  1.00  0.00           C
ATOM     12  HB1 ALA     2      29.172  28.440  21.584  1.00  0.00           H
ATOM     13  HB2 ALA     2      30.691  28.627  22.488  1.00  0.00           H
ATOM     14  HB3 ALA     2      29.844  27.070  22.499  1.00  0.00           H
ATOM     15  C   ALA     2      28.737  30.119  23.778  1.00  0.00           C
ATOM     16  O   ALA     2      27.675  30.634  23.429  1.00  0.00           O
ATOM     17  N   NME     3      29.784  30.841  24.197  1.00  0.00           N
ATOM     18  H   NME     3      30.622  30.348  24.461  1.00  0.00           H
ATOM     19  CH3 NME     3      29.784  32.300  24.293  1.00  0.00           C
ATOM     20 HH31 NME     3      28.951  32.628  24.918  1.00  0.00           H
ATOM     21 HH32 NME     3      30.720  32.652  24.729  1.00  0.00           H
ATOM     22 HH33 NME     3      29.663  32.734  23.299  1.00  0.00           H
TER      23      NME     3
END
)DATA"_pdb;
  REQUIRE(mol);
  // Oxygen in ACE (3rd heavy atom in mol) should be C=O, i.e. not OH
  CHECK(mol->getAtomWithIdx(2)->getTotalNumHs() == 0);
  CHECK(mol->getAtomWithIdx(2)->getTotalDegree() == 1);
  auto bond = mol->getBondBetweenAtoms(1, 2);
  REQUIRE(bond);
  CHECK(bond->getBondType() == Bond::BondType::DOUBLE);
}

TEST_CASE(
    "github #5327: MolFromMolBlock should correctly assign stereochemistry "
    "to 3D molecules") {
  SECTION("basics") {
    auto m = R"CTAB(
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.900794 -0.086835 0.009340 0
M  V30 2 C -0.552652 0.319534 0.077502 0
M  V30 3 F -0.861497 0.413307 1.437370 0
M  V30 4 Cl -0.784572 1.925710 -0.672698 0
M  V30 5 O -1.402227 -0.583223 -0.509512 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4
M  V30 4 1 2 5
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getConformer().is3D());
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
  }
  SECTION("wiggly bond") {
    auto m = R"CTAB(
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.900794 -0.086835 0.009340 0
M  V30 2 C -0.552652 0.319534 0.077502 0
M  V30 3 F -0.861497 0.413307 1.437370 0
M  V30 4 Cl -0.784572 1.925710 -0.672698 0
M  V30 5 O -1.402227 -0.583223 -0.509512 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4 CFG=2
M  V30 4 1 2 5
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getConformer().is3D());
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
  SECTION("3D as 2D") {
    // here we lie to the RDKit and tell it that a 3D conformer is 2D,
    //  the code detects that and still sets the conformer to be 3D and
    //  assigns stereo:
    auto m = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.900794 -0.086835 0.009340 0
M  V30 2 C -0.552652 0.319534 0.077502 0
M  V30 3 F -0.861497 0.413307 1.437370 0
M  V30 4 Cl -0.784572 1.925710 -0.672698 0
M  V30 5 O -1.402227 -0.583223 -0.509512 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4
M  V30 4 1 2 5
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getConformer().is3D());
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
  }
  SECTION("2D as 3D") {
    // here we lie to the RDKit and tell it that a 2D conformer is 3D,
    //  there's no chiral volume, so we don't end up with a chiral center
    auto m = R"CTAB(
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.299038 -0.750000 0.000000 0
M  V30 2 C 0.000000 -0.000000 0.000000 0
M  V30 3 F 0.750000 -1.299038 0.000000 0
M  V30 4 Cl -0.750000 1.299038 0.000000 0
M  V30 5 O 1.299038 0.750000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4
M  V30 4 1 2 5
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getConformer().is3D());
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
  SECTION("double bond") {
    auto m = R"CTAB(
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1.911935 -0.058147 -0.007384 0
M  V30 2 C 0.477913 -0.091130 -0.413392 0
M  V30 3 C -0.494810 0.079132 0.449979 0
M  V30 4 C -1.932350 0.037738 -0.006356 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getConformer().is3D());
    CHECK(Chirality::translateEZLabelToCisTrans(
              m->getBondWithIdx(1)->getStereo()) ==
          Bond::BondStereo::STEREOTRANS);
  }
  SECTION("double bond, crossed") {
    auto m = R"CTAB(
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1.911935 -0.058147 -0.007384 0
M  V30 2 C 0.477913 -0.091130 -0.413392 0
M  V30 3 C -0.494810 0.079132 0.449979 0
M  V30 4 C -1.932350 0.037738 -0.006356 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3 CFG=2
M  V30 3 1 3 4
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getConformer().is3D());
    CHECK(Chirality::translateEZLabelToCisTrans(
              m->getBondWithIdx(1)->getStereo()) ==
          Bond::BondStereo::STEREOANY);
  }
  SECTION("double bond, wiggly bond") {
    auto m = R"CTAB(
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1.911935 -0.058147 -0.007384 0
M  V30 2 C 0.477913 -0.091130 -0.413392 0
M  V30 3 C -0.494810 0.079132 0.449979 0
M  V30 4 C -1.932350 0.037738 -0.006356 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4 CFG=2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getConformer().is3D());
    CHECK(Chirality::translateEZLabelToCisTrans(
              m->getBondWithIdx(1)->getStereo()) ==
          Bond::BondStereo::STEREOANY);
  }
  SECTION("non-tetrahedral") {
    auto m = R"CTAB(
  Mrv2108 05252216313D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.7191 0.2488 -3.5085 0
M  V30 2 As -1.0558 1.9209 -2.6345 0
M  V30 3 F -0.4636 3.422 -1.7567 0
M  V30 4 O -2.808 2.4243 -2.1757 0
M  V30 5 Cl -0.1145 2.6609 -4.5048 0
M  V30 6 Br 0.2255 0.6458 -1.079 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4
M  V30 4 1 2 5
M  V30 5 1 2 6
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getConformer().is3D());
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL);
  }
  SECTION("non-tetrahedral, wiggly") {
    auto m = R"CTAB(
  Mrv2108 05252216313D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.7191 0.2488 -3.5085 0
M  V30 2 As -1.0558 1.9209 -2.6345 0
M  V30 3 F -0.4636 3.422 -1.7567 0
M  V30 4 O -2.808 2.4243 -2.1757 0
M  V30 5 Cl -0.1145 2.6609 -4.5048 0
M  V30 6 Br 0.2255 0.6458 -1.079 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4
M  V30 4 1 2 5 CFG=2
M  V30 5 1 2 6
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getConformer().is3D());
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
}

TEST_CASE("Force use of MolBlock wedges", "") {
  SECTION("basics") {
    auto m = R"CTAB(bad wedging
  ChemDraw07092209022D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.714471 0.825000 0.000000 0
M  V30 2 C -0.714471 0.000000 0.000000 0
M  V30 3 C -0.000000 -0.412500 0.000000 0
M  V30 4 C 0.714471 0.000000 0.000000 0
M  V30 5 C 0.714471 0.825000 0.000000 0
M  V30 6 C -0.000000 1.237500 0.000000 0
M  V30 7 C -0.000000 -1.237500 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4 CFG=1
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 1
M  V30 7 1 3 7
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getBondWithIdx(2)->getBondDir() == Bond::BondDir::NONE);
    reapplyMolBlockWedging(*m);
    CHECK(m->getBondWithIdx(2)->getBondDir() == Bond::BondDir::BEGINWEDGE);
  }
  SECTION("GitHub5448") {
    {
      auto m = R"CTAB(
  ChemDraw07232208492D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 12 12 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C 1.151421 0.903801 0.000000 0
M  V30 2 C 1.151421 0.078801 0.000000 0
M  V30 3 N 1.936021 -0.176200 0.000000 0
M  V30 4 C 2.420921 0.491301 0.000000 0
M  V30 5 N 1.936021 1.158699 0.000000 0
M  V30 6 C 0.436921 -0.333699 0.000000 0
M  V30 7 C -0.277478 0.078801 0.000000 0
M  V30 8 C -0.991978 -0.333699 0.000000 0
M  V30 9 C 0.436921 -1.158699 0.000000 0
M  V30 10 C -1.706442 0.078813 0.000000 0
M  V30 11 C -2.420921 -0.333674 0.000000 0
M  V30 12 F -1.706428 0.903813 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 1 5 1
M  V30 6 1 2 6
M  V30 7 1 6 7
M  V30 8 2 7 8 CFG=2
M  V30 9 1 6 9 CFG=2
M  V30 10 1 8 10
M  V30 11 1 10 11
M  V30 12 1 10 12 CFG=1
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 10)
M  V30 END COLLECTION
M  V30 END CTAB
M  END)CTAB"_ctab;
      REQUIRE(m);
      WedgeMolBonds(*m, &m->getConformer());
      CHECK(m->getBondWithIdx(10)->getBondDir() == Bond::BondDir::BEGINWEDGE);
      CHECK(m->getBondWithIdx(11)->getBondDir() == Bond::BondDir::NONE);
      reapplyMolBlockWedging(*m);
      CHECK(m->getBondWithIdx(10)->getBondDir() == Bond::BondDir::NONE);
      CHECK(m->getBondWithIdx(11)->getBondDir() == Bond::BondDir::BEGINWEDGE);
    }
  }
}

TEST_CASE(
    "GitHub Issue #5423: Parsing a Mol block/file does not clear the "
    "\"molTotValence\" property from atoms") {
  auto m = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N -3.657143 -0.742857 0.000000 0 CHG=1 VAL=4
M  V30 END ATOM
M  V30 END CTAB
M  END)CTAB"_ctab;
  REQUIRE(m);
  CHECK(!m->getAtomWithIdx(0)->hasProp(common_properties::molTotValence));
}

TEST_CASE("Github #5433: PRECONDITION error with nonsense molecule") {
  auto m = R"CTAB(
  SomeFailingMolFile

  8  8  0  0  1  0            999 V2000
   -0.7145    1.2375    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    0.0000    0.8250    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7145   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7145   -1.2375    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.6500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145   -1.2375    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  6  0  0  0
  2  3  2  0  0  0  0
  3  4  1  4  0  0  0
  3  8  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
  7  8  2  0  0  0  0
M  END)CTAB"_ctab;
  REQUIRE(m);
}