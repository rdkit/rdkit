//
//  Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <string>

#include <catch2/catch_all.hpp>

#include <RDGeneral/test.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/test_fixtures.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>

using namespace RDKit;
using namespace std;

TEST_CASE("Testing molecules which should parse.") {
  int i = 0;
  ROMol *mol, *mol2;
  string smis[] = {
      "C1CC2C1CC2",
      "c1cccn(=O)c1",
      "C",
      "CC",
      "C-C",
      "C=C",
      "[CH2+]C[CH+2]",
      "C1CC1",
      "C1CC=1",
      "C=1CC1",
      "C=C-O",
      "C1CC1",
      "C1NC1",
      "C1=CC1",
      "C1CCC1",
      "CC(C)CC",
      "CC(=O)O",
      "C1C(=O)C1",
      "C1C(N)C1",
      "CC(O)C",
      "OC=CCC",
      "CC([O-])O",
      "C1CC2C1CC2",
      "Cl/C=C/Cl",
      "Cl/C=C\\Cl",
      "Cl/C=C/Cl",
      "Cl/C=C\\Cl",
      "Cl/C=C\\\\Cl",
      "C1CC.CC1",
      "C1C(C2CC2).C2CC2C1",
      "[Na+].[Cl-].[NH4+].[Cl-]",
      "C[35Cl]",
      "C%10CC%10",
      "[H][H]",
      "[H+]",
      "C[N+](=O)[O-]",
      "N1C(=N)SC=C1",
      "[O-][N+](=O)C1=CNC(=N)S1",
      "CN(=O)=O",
      "C1=CC=C[N+]([O-])=C1",
      "C1=CC=CN(=O)=C1",
      // test whitespace tolerance:
      "  C1=CC=CN(=O)=C1",
      "C1=CC=CN(=O)=C1  ",
      "  C1=CC=CN(=O)=C1  ",
      "\tC1=CC=CN(=O)=C1\r\n",
      // test dummy atoms:
      "c1ccccc1[*]",
      "c1ccccc1[1*]",
      "S1cccc1",
      "*1ccccc1",
      "C1=CC=CC=C1",
      "*1=CC=CC=C1",
      "*1*cccc1",
      "*1**ccc1",
      // test aromatic se and te:
      "c1ccc[se]1",
      "c1ccc[te]1",
      // test zeros as ring indices, issue 2690982:
      "C0CC0",
      // test canonization error, issue 3018558:
      "C/C(/C=C2\\Sc1ccc(cc1N\\2C))=C5\\SC4=NccN4C\\5=O",
      // "the most common molecule in the universe",
      // expressed in an ugly way:
      "[HH]",
      "[2HH]",
      "[HH2-]",   // issue 3535669
      "[2HH2-]",  // issue 3535669
      // problems handling aromatic boron, issue 3480481
      "b1ccccc1",
      "C[Rf]C",  // issue 3535668
      "[C:1]",
      "[C:0]",           // issue 3525776
      "[si]1cccc[si]1",  // aromatic Si (github issue #5)
      "[asH]1cccc1",     // aromatic As (github issue #682)
      "[Db][Sg][Bh][Hs][Mt][Ds][Rg][Cn][Nh][Fl][Mc][Lv][Ts][Og]",  // new
                                                                   // elements
      "[Uun][Uuu][Uub][Uut][Uuq][Uup][Uuh][Uus][Uuo]",  // old names for new
                                                        // elements
      "['Db']['Sg']['Bh']['Hs']['Mt']['Ds']['Rg']['Cn']['Nh']['Fl']['Mc']['Lv']['"
      "Ts']['Og']",  // a biovia pathology
      "[#6]",        // feature borrowed from SMARTS
      "[12#6]",
      "C$C",  // quadruple bonds
      // extended chirality
      "C[Fe@TH](O)(Cl)F",
      "C[Fe@TH1](O)(Cl)F",
      "C[Fe@SP](O)(Cl)F",
      "C[Fe@SP1](O)(Cl)F",
      "C[Fe@TB](O)(Cl)(Br)F",
      "C[Fe@TB10](O)(Cl)(Br)F",
      "C[Fe@OH](O)(Cl)(Br)(N)F",
      "C[Fe@OH20](O)(Cl)(Br)(N)F",
      "EOS",
  };
  while (smis[i] != "EOS") {
    string smi = smis[i];
    mol = SmilesToMol(smi);
    REQUIRE(mol);
    if (mol) {
      unsigned int nAts = mol->getNumAtoms();
      REQUIRE(nAts != 0);
      smi = MolToSmiles(*mol);
      mol2 = SmilesToMol(smi);
      REQUIRE(mol2->getNumAtoms() == nAts);
      delete mol;
      delete mol2;
    }
    i++;
  }
}

TEST_CASE("Testing molecules which should fail to parse/sanitize.") {
  int i = 0;
  ROMol *mol;

  // alternate good and bad smiles here to ensure that the parser can resume
  // parsing on good input:
  string smis[] = {
      "CC=(CO)C",    "CC(=CO)C", "C1CC",
      "C1CC1",       "Ccc",      "CCC",
      "fff",  // tests the situation where the parser cannot do anything at all
      "CCC",
      "N(=O)(=O)=O",  // bad sanitization failure
      "C1CC1",
      "C=0",  // part of sf.net issue 2525792
      "C1CC1",
      "C0",  // part of sf.net issue 2525792
      "C1CC1",
      "C-0",  // part of sf.net issue 2525792
      "C1CC1",
      "C+0",  // part of sf.net issue 2525792
      "C1CC1",       "[H2H]",    "C1CC1",
      "[HH2]",       "C1CC1",    "[555555555555555555C]",
      "C1CC1",             //
      "[Fe@TD]",     "C",  //
      "[Fe@TH3]",    "C",  //
      "[Fe@SP4]",    "C",  //
      "[Fe@AL3]",    "C",  //
      "[Fe@TB21]",   "C",  //
      "[Fe@OH31]",   "C",  //
      "baz",         "C",  //
      "EOS",
  };

  // turn off the error log temporarily:
  while (smis[i] != "EOS") {
    string smi = smis[i];
    boost::logging::disable_logs("rdApp.error");
    try {
      mol = SmilesToMol(smi);
    } catch (MolSanitizeException &) {
      mol = (ROMol *)nullptr;
    }
    boost::logging::enable_logs("rdApp.error");
    if (!(i % 2)) {
      REQUIRE(!mol);
    } else {
      REQUIRE(mol);
      delete mol;
    }
    i++;
  }
}

TEST_CASE("Testing details") {
  ROMol *mol;
  Atom *a;
  std::string smi;

  // implicit/explicit H handling
  smi = "OC([OH])C[O-]";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 5);
  a = mol->getAtomWithIdx(0);
  REQUIRE(a->getValence(Atom::ValenceType::IMPLICIT) == 1);
  REQUIRE(a->getValence(Atom::ValenceType::EXPLICIT) == 1);
  REQUIRE(a->getNoImplicit() == 0);
  REQUIRE(a->getFormalCharge() == 0);
  a = mol->getAtomWithIdx(2);
  REQUIRE(a->getValence(Atom::ValenceType::IMPLICIT) == 0);
  REQUIRE(a->getValence(Atom::ValenceType::EXPLICIT) == 2);
  REQUIRE(a->getNoImplicit() == 1);
  REQUIRE(a->getFormalCharge() == 0);
  a = mol->getAtomWithIdx(4);
  REQUIRE(a->getValence(Atom::ValenceType::IMPLICIT) == 0);
  REQUIRE(a->getValence(Atom::ValenceType::EXPLICIT) == 1);
  REQUIRE(a->getNoImplicit() == 1);
  REQUIRE(a->getFormalCharge() == -1);

  delete mol;
}

TEST_CASE("Testing smiles that were previously problems") {
  ROMol *mol;
  std::string smi;

  // ring closure handling with branches/fragments
  VECT_INT_VECT rings;
  smi = "C1(CC1CC1CC1)";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  int ringCount = MolOps::findSSSR(*mol, rings);
  REQUIRE(ringCount == 2);
  REQUIRE(rings.size() == 2);
  REQUIRE(rings[0].size() == 3);
  REQUIRE(rings[1].size() == 3);

  // this is truly pathological, but both daylight
  //   and chemdraw parse it properly
  smi = "C1.C1CC1CC1";
  delete mol;
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  ringCount = MolOps::findSSSR(*mol, rings);
  REQUIRE(ringCount == 1);
  REQUIRE(rings.size() == 1);
  REQUIRE(rings[0].size() == 3);

  // here's another stupid case that we need to handle:
  delete mol;
  smi = "C1CC11CC1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  ringCount = MolOps::findSSSR(*mol, rings);
  REQUIRE(ringCount == 2);
  REQUIRE(rings.size() == 2);
  REQUIRE(rings[0].size() == 3);
  REQUIRE(rings[1].size() == 3);

  delete mol;
}

TEST_CASE("Testing basic SMILES canonicalization") {
  ROMol *mol;
  std::string smi, refSmi;

  smi = "C1OCCCC1";
  mol = SmilesToMol(smi);
  refSmi = MolToSmiles(*mol);

  delete mol;
  smi = "C1COCCC1";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol);
  REQUIRE(refSmi == smi);

  delete mol;
  smi = "O1CCCCC1";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol);
  REQUIRE(refSmi == smi);

  delete mol;
  smi = "OC=CC";
  mol = SmilesToMol(smi);
  refSmi = MolToSmiles(*mol);
  delete mol;
  smi = "CC=CO";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol);
  REQUIRE(refSmi == smi);
  delete mol;
  smi = "C(C)=CO";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol);
  REQUIRE(refSmi == smi);
  delete mol;
  smi = "C(O)=CC";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol);
  REQUIRE(refSmi == smi);

  // --- These are related to Issue 109
  delete mol;
  smi = "C([H])Cl";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getNumAtoms() == 2);
  refSmi = MolToSmiles(*mol);
  delete mol;
  smi = "CCl";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol);
  REQUIRE(refSmi == smi);
  delete mol;
  // -- Issue 131
  smi = "P#[Ga]";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getNumAtoms() == 2);
  refSmi = MolToSmiles(*mol);
  delete mol;
  mol = SmilesToMol(refSmi);
  smi = MolToSmiles(*mol);
  REQUIRE(refSmi == smi);

  delete mol;
  smi = "O=[Ba]";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getNumAtoms() == 2);
  refSmi = MolToSmiles(*mol);
  delete mol;
  mol = SmilesToMol(refSmi);
  smi = MolToSmiles(*mol);
  REQUIRE(refSmi == smi);

  // make sure empty molecules return empty SMILES:
  delete mol;
  mol = new ROMol();
  smi = MolToSmiles(*mol);
  REQUIRE(smi == "");

  delete mol;
}

TEST_CASE("Testing handling of stereochemical smiles") {
  ROMol *mol;
  std::string smi, refSmi, cip;

  smi = "F[C@](Cl)(Br)I";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  REQUIRE(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  refSmi = MolToSmiles(*mol, 1);

  delete mol;
  smi = "F[C@](Br)(I)Cl";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  REQUIRE(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "F[C@](I)(Cl)Br";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "Cl[C@](Br)(F)I";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "Cl[C@](F)(I)Br";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "I[C@](F)(Br)Cl";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "I[C@](Br)(Cl)F";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "F[C@@](Br)(Cl)I";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "F[C@@](Cl)(I)Br";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "Cl[C@@](Br)(I)F";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "Cl[C@@](F)(Br)I";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "[C@@](Cl)(F)(Br)I";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "F[C@H](Cl)Br";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  refSmi = MolToSmiles(*mol, 1);

  delete mol;
  smi = "Br[C@H](F)Cl";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  smi = MolToSmiles(*mol, 1);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "Br[C@]([H])(F)Cl";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "Br[C@](F)(Cl)[H]";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "Br[C@]1(F)(Cl).[H]1";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "Br[C@H]1Cl.F1";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "Br[C@]12Cl.F2.[H]1";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "Br[C@]21Cl.F1.[H]2";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "[C@@H](Br)(F)Cl";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  delete mol;
  smi = "[H][C@@](Br)(F)Cl";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  smi = MolToSmiles(*mol, 1);
  REQUIRE(smi == refSmi);

  // an additional set of test cases from the Chirality notes document.
  // one can never have too many tests of this stuff.
  delete mol;
  smi = "F[C@]([H])(O)C";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  delete mol;
  smi = "F[C@]1([H])OC1";
  mol = SmilesToMol(smi);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");

  delete mol;
  smi = "F[C@H](O)C";
  mol = SmilesToMol(smi);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  delete mol;
  smi = "F[C@@H]1OC1";
  mol = SmilesToMol(smi);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");

  delete mol;
  smi = "[C@](F)([H])(O)C";
  mol = SmilesToMol(smi);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  delete mol;
  smi = "[C@@]1(F)([H])OC1";
  mol = SmilesToMol(smi);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");

  delete mol;
  smi = "[C@@H](F)(O)C";
  mol = SmilesToMol(smi);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  smi = MolToSmiles(*mol, true);
  REQUIRE(smi == "C[C@@H](O)F");
  smi = MolToSmiles(*mol, true, false, 0);
  REQUIRE(smi == "[C@H](C)(O)F");

  delete mol;
  smi = "[C@@H]1(F)OC1";
  mol = SmilesToMol(smi);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  smi = MolToSmiles(*mol, true);
  REQUIRE(smi == "F[C@H]1CO1");
  smi = MolToSmiles(*mol, true, false, 0);
  REQUIRE(smi == "[C@H]1(F)CO1");

  delete mol;
  smi = "C1O[C@H]1F";
  mol = SmilesToMol(smi);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  delete mol;
  smi = "C1O[C@@]1([H])F";
  mol = SmilesToMol(smi);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");

  // -----------------------------------
  // test some double-bond containing molecules:

  //-- cis --
  delete mol;
  smi = "F\\C=C/Br";
  mol = SmilesToMol(smi);
  refSmi = MolToSmiles(*mol, 1);

  delete mol;
  mol = SmilesToMol(refSmi);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);

  delete mol;
  smi = "Br\\C=C/F";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);

  delete mol;
  smi = "Br/C=C\\F";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);

  delete mol;
  smi = "F/C=C\\Br";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);

  //-- trans --
  delete mol;
  smi = "F\\C=C\\Br";
  mol = SmilesToMol(smi);
  refSmi = MolToSmiles(*mol, 1);

  delete mol;
  mol = SmilesToMol(refSmi);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);

  delete mol;
  smi = "Br\\C=C\\F";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);

  delete mol;
  smi = "Br/C=C/F";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);

  delete mol;
  smi = "F/C=C/Br";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);

  //-- more complex --
  delete mol;
  smi = "F\\C=C(/Cl)\\Br";
  mol = SmilesToMol(smi);
  refSmi = MolToSmiles(*mol, 1);

  delete mol;
  mol = SmilesToMol(refSmi);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);

  delete mol;
  smi = "F/C=C(\\Cl)/Br";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);

  delete mol;
  smi = "F/C=C(\\Cl)Br";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);

  delete mol;
  smi = "F/C=C(Cl)/Br";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);

  //-- combine chirality with cis/trans --
  delete mol;
  smi = "F[C@H](Cl)\\C=C(/F)";
  mol = SmilesToMol(smi);
  refSmi = MolToSmiles(*mol, 1);

  delete mol;
  mol = SmilesToMol(refSmi);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);

  delete mol;
  smi = "F[C@H](Cl)/C=C(\\F)";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);
  delete mol;

  smi = "Cl[C@@H](F)/C=C(\\F)";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);
  delete mol;

  smi = "Cl[C@@H](F)\\C=C(/F)";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);
  delete mol;
}

TEST_CASE("Testing Issue 127 (chiral smiles with fused rings)") {
  ROMol *mol, *mol2;
  std::string smi, refSmi, tempStr;

  smi = "Cl[C@]12[Si]C(C2)O1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);

  // first roundtrip the non-chiral SMILES:
  refSmi = MolToSmiles(*mol);
  mol2 = SmilesToMol(refSmi);
  REQUIRE(mol2);
  tempStr = MolToSmiles(*mol2);
  REQUIRE(refSmi == tempStr);
  delete mol2;

  // now do the true SMILES:
  refSmi = MolToSmiles(*mol, 1);
  mol2 = SmilesToMol(refSmi);
  REQUIRE(mol2);
  tempStr = MolToSmiles(*mol2, 1);
  REQUIRE(refSmi == tempStr);
  delete mol2;
  delete mol;
}

TEST_CASE("Testing Issue 143 (removing chiral tags for non-chiral centers)") {
  ROMol *mol;
  std::string smi, refSmi, tempStr;

  smi = "C[C@](C)(C)C";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  refSmi = MolToSmiles(*mol, true);
  REQUIRE(refSmi == "CC(C)(C)C");
  delete mol;

  smi = "CC[C@](C)(C)C=O";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  refSmi = MolToSmiles(*mol, true);
  REQUIRE(refSmi == "CCC(C)(C)C=O");
  delete mol;
}

TEST_CASE(
    "Testing Issue 151 (Chiral centers in rings with hydrogen on them not handled correctly") {
  ROMol *mol, *mol2;
  std::string smi, refSmi, tempStr;

  smi = "C1S[C@H]1O";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  REQUIRE(mol->getAtomWithIdx(2)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  REQUIRE(mol->getAtomWithIdx(2)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);

  refSmi = MolToSmiles(*mol, true);
  REQUIRE(refSmi == "O[C@H]1CS1");
  mol2 = SmilesToMol(refSmi);
  REQUIRE(mol2);
  smi = MolToSmiles(*mol2, true);
  REQUIRE(refSmi == smi);
  delete mol;
  delete mol2;

  smi = "F[C@@H]1O[C@H](Cl)S1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  REQUIRE(mol->getAtomWithIdx(2)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  REQUIRE(mol->getAtomWithIdx(4)->getChiralTag() == Atom::CHI_UNSPECIFIED);

  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  REQUIRE(mol->getAtomWithIdx(3)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  REQUIRE(mol->getAtomWithIdx(3)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);

  refSmi = MolToSmiles(*mol, true);
  REQUIRE(refSmi == "F[C@@H]1O[C@H](Cl)S1");
  mol2 = SmilesToMol(refSmi);
  REQUIRE(mol2);
  smi = MolToSmiles(*mol2, true);
  REQUIRE(refSmi == smi);
  delete mol;
  delete mol2;

  smi = "Cl[C@@H]1S[C@@H](O1)F";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  REQUIRE(mol->getAtomWithIdx(2)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  REQUIRE(mol->getAtomWithIdx(4)->getChiralTag() == Atom::CHI_UNSPECIFIED);

  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  REQUIRE(mol->getAtomWithIdx(3)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  REQUIRE(mol->getAtomWithIdx(3)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);

  refSmi = MolToSmiles(*mol, true);
  REQUIRE(refSmi == "F[C@@H]1O[C@H](Cl)S1");
  mol2 = SmilesToMol(refSmi);
  REQUIRE(mol2);
  smi = MolToSmiles(*mol2, true);
  REQUIRE(refSmi == smi);
  delete mol;
  delete mol2;

  smi = "Cl[C@@H]1O[C@H](F)S1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  REQUIRE(mol->getAtomWithIdx(2)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  REQUIRE(mol->getAtomWithIdx(4)->getChiralTag() == Atom::CHI_UNSPECIFIED);

  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  REQUIRE(mol->getAtomWithIdx(3)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  REQUIRE(mol->getAtomWithIdx(3)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);

  refSmi = MolToSmiles(*mol, true);
  REQUIRE(refSmi == "F[C@H]1O[C@@H](Cl)S1");
  mol2 = SmilesToMol(refSmi);
  REQUIRE(mol2);
  smi = MolToSmiles(*mol2, true);
  REQUIRE(refSmi == smi);
  delete mol;
  delete mol2;
}

TEST_CASE(
    "Testing Issue 153 (Incorrect order of ring-closure bonds from SMILES)") {
  std::string code;
  ROMol *mol, *mol2;
  std::string smi, refSmi, tempStr;

  for (const bool useLegacy : {true, false}) {
    UseLegacyStereoPerceptionFixture fx(useLegacy);
    smi = "C1(O[C@H]12)S2";
    mol = SmilesToMol(smi);
    REQUIRE(mol);
    REQUIRE(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    REQUIRE(mol->getAtomWithIdx(2)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    REQUIRE(mol->getAtomWithIdx(2)->getChiralTag() ==
            Atom::CHI_TETRAHEDRAL_CCW);
    if (useLegacy) {
      MolOps::assignStereochemistry(*mol);
    } else {
      CIPLabeler::assignCIPLabels(*mol);
    }
    REQUIRE(mol->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, code);
    REQUIRE(code == "S");

    refSmi = MolToSmiles(*mol, true);
    REQUIRE(refSmi == "O1C2S[C@H]12");
    mol2 = SmilesToMol(refSmi);
    REQUIRE(mol2);
    smi = MolToSmiles(*mol2, true);
    REQUIRE(refSmi == smi);
    delete mol;
    delete mol2;

    smi = "C1(O[C@H]21)S2";
    mol = SmilesToMol(smi);
    REQUIRE(mol);
    REQUIRE(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    REQUIRE(mol->getAtomWithIdx(2)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    REQUIRE(mol->getAtomWithIdx(2)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
    if (useLegacy) {
      MolOps::assignStereochemistry(*mol);
    } else {
      CIPLabeler::assignCIPLabels(*mol);
    }
    REQUIRE(mol->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, code);
    REQUIRE(code == "R");

    refSmi = MolToSmiles(*mol, true);
    REQUIRE(refSmi == "O1C2S[C@@H]12");
    mol2 = SmilesToMol(refSmi);
    REQUIRE(mol2);
    smi = MolToSmiles(*mol2, true);
    REQUIRE(refSmi == smi);
    delete mol;
    delete mol2;
  }
}

TEST_CASE(
    "Testing Issue 157 (Symmetric molecules with multiple chiral centers badly canonicalized)") {
  std::string code;
  ROMol *mol, *mol2;
  std::string smi, refSmi, tempStr;

  smi = "O[C@](C)(Cl)[C@@](O)(Cl)C";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  REQUIRE(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
  REQUIRE(mol->getAtomWithIdx(4)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);

  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "S");

  refSmi = MolToSmiles(*mol, true);
  mol2 = SmilesToMol(refSmi);
  REQUIRE(mol2);
  smi = MolToSmiles(*mol2, true);
  REQUIRE(refSmi == smi);
  delete mol;
  delete mol2;

  smi = "Cl[C@@](C)1CC[C@@](C)(C1)Cl";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  refSmi = MolToSmiles(*mol, true);
  mol2 = SmilesToMol(refSmi);
  REQUIRE(mol2);
  smi = MolToSmiles(*mol2, true);
  REQUIRE(refSmi == smi);
  delete mol;
  delete mol2;

  smi = "[H][C@@]12CC(CO1)CN2";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, smi);
  REQUIRE(smi == "S");
  refSmi = MolToSmiles(*mol, true);
  mol2 = SmilesToMol(refSmi);
  REQUIRE(mol2);
  smi = MolToSmiles(*mol2, true);
  REQUIRE(refSmi == smi);
  delete mol;
  delete mol2;
  smi = "[H][C@@]12C[14C@@](C=C1)(C3C2C(NC3=O)=O)[H]";

  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, smi);
  REQUIRE(smi == "R");
  mol->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, smi);
  REQUIRE(smi == "S");

  refSmi = MolToSmiles(*mol, true);
  mol2 = SmilesToMol(refSmi);
  REQUIRE(mol2);
  smi = MolToSmiles(*mol2, true);
  REQUIRE(refSmi == smi);
  delete mol;
  delete mol2;
}

TEST_CASE("Testing Issue 159 (cis/trans wrong in some branched systems)") {
  ROMol *mol;
  std::string smi, refSmi, tempStr;

  smi = "C/C=C/O";
  mol = SmilesToMol(smi);
  REQUIRE(mol);

  REQUIRE(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  refSmi = MolToSmiles(*mol, 1);

  delete mol;
  smi = "C(\\C)=C/O";
  mol = SmilesToMol(smi);
  REQUIRE(mol);

  REQUIRE(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);

  delete mol;
  smi = "C(\\\\C)=C/O";
  mol = SmilesToMol(smi);
  REQUIRE(mol);

  REQUIRE(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);

  delete mol;
  smi = "C(=C/O)\\C";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(0)->getStereo() == Bond::STEREOE);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);

  delete mol;
  smi = "C(\\C/C=C/Cl)=C/O";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(4)->getStereo() == Bond::STEREOE);
  REQUIRE(mol->getBondWithIdx(2)->getStereo() == Bond::STEREOE);

  delete mol;
  smi = "O=C\\C=C/F";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(0)->getBondType() == Bond::DOUBLE);
  REQUIRE(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
  REQUIRE(mol->getBondWithIdx(2)->getStereo() == Bond::STEREOZ);

  delete mol;
  smi = "C(/C=O)=C/F";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);
  REQUIRE(mol->getBondWithIdx(2)->getStereo() == Bond::STEREOZ);

  delete mol;
  smi = "C(=C/F)/C=O";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(0)->getStereo() == Bond::STEREOZ);
  REQUIRE(mol->getBondWithIdx(3)->getStereo() == Bond::STEREONONE);

  delete mol;
  smi = "C(=O)\\C=C/Br";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(2)->getStereo() == Bond::STEREOZ);
  REQUIRE(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);

  delete mol;
  smi = "CC(=O)\\C=C/Br";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(3)->getStereo() == Bond::STEREOZ);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

  delete mol;
  smi = "C(=O)\\N=C\\Br";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(2)->getStereo() == Bond::STEREOE);
  REQUIRE(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);

  delete mol;
  smi = "CC(=O)\\N=C\\Br";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(3)->getStereo() == Bond::STEREOE);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

  delete mol;
  smi = "C(/Br)(=C/Cl)Cl";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);

  delete mol;
  smi = "C(=C/Cl)(/Br)Cl";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(0)->getStereo() == Bond::STEREOZ);

  delete mol;
  smi = "Cl\\C=C(\\Br)";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);

  delete mol;
  smi = "Cl\\C(=C\\Br)";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);

  delete mol;
  smi = "C(/C=C/C)";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  delete mol;
  smi = "C(/C)=C/C";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);

  // ---------
  // These next few molecules test propagation of bond flips:
  // ---------
  delete mol;
  smi = "Cl/C=C(/C=C/C)\\C=C\\Br";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
  REQUIRE(mol->getBondWithIdx(3)->getStereo() == Bond::STEREOE);
  REQUIRE(mol->getBondWithIdx(6)->getStereo() == Bond::STEREOE);

  delete mol;
  smi = "C(/C=C/C)(\\C=C\\Br)=C\\Cl";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  REQUIRE(mol->getBondWithIdx(4)->getStereo() == Bond::STEREOE);
  REQUIRE(mol->getBondWithIdx(6)->getStereo() == Bond::STEREOZ);

  delete mol;
  smi = "Br/C=C/C(/C=C/C)=C\\Cl";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  REQUIRE(mol->getBondWithIdx(4)->getStereo() == Bond::STEREOE);
  REQUIRE(mol->getBondWithIdx(6)->getStereo() == Bond::STEREOZ);

  delete mol;
  smi = "Cl/C=C(/C=C/C=C\\F)\\C=C\\Br";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
  REQUIRE(mol->getBondWithIdx(3)->getStereo() == Bond::STEREOE);
  REQUIRE(mol->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
  REQUIRE(mol->getBondWithIdx(8)->getStereo() == Bond::STEREOE);
  delete mol;
}

TEST_CASE("Testing Issue 175 (cis/trans wrong on ring closures)") {
  ROMol *mol;
  std::string smi, refSmi, tempStr;

  smi = "Cl\\C=C1.F/1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  delete mol;

  smi = "Cl\\C=C1CN/1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  delete mol;

  smi = "C/1=C/F.F1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(0)->getStereo() == Bond::STEREOZ);
  delete mol;
}

TEST_CASE("Testing Issue 176 (problems with 'mol BOND ring_number')") {
  ROMol *mol;
  std::string smi, refSmi, tempStr;

  smi = "C1CC1C1CC1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getNumBonds() == 7);

  delete mol;
  smi = "C1CC1C1CC-1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getNumBonds() == 7);

  delete mol;
  smi = "C1CC1C1CC=1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getNumBonds() == 7);

  delete mol;
  smi = "C1CC1C=1CC1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getNumBonds() == 7);

  delete mol;
}

TEST_CASE("Testing Issue 180: Z/E problems") {
  ROMol *mol;
  std::string smi, refSmi;

  smi = "Cl/C(=N\\O)/C(=N\\O)Br";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
  REQUIRE(mol->getBondWithIdx(4)->getStereo() == Bond::STEREOE);
  refSmi = MolToSmiles(*mol, 1);

  delete mol;
  smi = "Cl/C(/C(Br)=N\\O)=N\\O";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(3)->getStereo() == Bond::STEREOE);
  REQUIRE(mol->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);
  delete mol;
}

TEST_CASE("Testing Issue 184: Cis/Trans incorrect on ring-closure bonds") {
  ROMol *mol;
  std::string smi, refSmi;

  smi = "C1NC(Cl)C(=N\\O)/C1=N\\O";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(4)->getBondType() == Bond::DOUBLE);
  REQUIRE(mol->getBondWithIdx(4)->getStereo() == Bond::STEREOZ);
  REQUIRE(mol->getBondWithIdx(7)->getBondType() == Bond::DOUBLE);
  REQUIRE(mol->getBondWithIdx(7)->getStereo() == Bond::STEREOZ);
  refSmi = MolToSmiles(*mol, 1);
  delete mol;
  mol = SmilesToMol(refSmi);
  REQUIRE(mol);

  for (RWMol::BondIterator bondIt = mol->beginBonds();
       bondIt != mol->endBonds(); bondIt++) {
    if ((*bondIt)->getBondType() == Bond::DOUBLE) {
      REQUIRE((*bondIt)->getStereo() == Bond::STEREOZ);
    }
  }

  smi = MolToSmiles(*mol, 1);
  REQUIRE(refSmi == smi);
  delete mol;
}

TEST_CASE("Testing Issue 185: Cis/Trans incorrect on writing branches") {
  ROMol *mol;
  std::string smi, refSmi;

  // start with a simple E/Z handling case with branches:
  smi = "C(/C)=N/O";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
  refSmi = MolToSmiles(*mol, 1, 0, 0);
  CHECK(refSmi == "C(/C)=N/O");
  delete mol;
  // make sure we can round-trip:
  mol = SmilesToMol(refSmi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
  delete mol;

  // now make it more complex
  smi = "CC(=N\\O)/C=P/N";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
  REQUIRE(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  REQUIRE(mol->getBondWithIdx(4)->getBondType() == Bond::DOUBLE);
  REQUIRE(mol->getBondWithIdx(4)->getStereo() == Bond::STEREOE);
  refSmi = MolToSmiles(*mol, 1);
  delete mol;
  mol = SmilesToMol(refSmi);
  REQUIRE(mol);

  for (RWMol::BondIterator bondIt = mol->beginBonds();
       bondIt != mol->endBonds(); bondIt++) {
    if ((*bondIt)->getBondType() == Bond::DOUBLE) {
      CHECK((*bondIt)->getStereo() == Bond::STEREOE);
    }
  }
  smi = MolToSmiles(*mol, 1);
  CHECK(refSmi == smi);

  // now repeat that experiment, but this time root the SMILES so that
  // we go in a "sensible" order:
  delete mol;
  smi = "CC(=N\\O)/C=P/N";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  refSmi = MolToSmiles(*mol, true, false, 6);
  CHECK(refSmi == "N/P=C/C(C)=N/O");
  delete mol;
  mol = SmilesToMol(refSmi);
  REQUIRE(mol);
  for (RWMol::BondIterator bondIt = mol->beginBonds();
       bondIt != mol->endBonds(); bondIt++) {
    if ((*bondIt)->getBondType() == Bond::DOUBLE) {
      CHECK((*bondIt)->getStereo() == Bond::STEREOE);
    }
  }
  delete mol;
}

TEST_CASE("Testing Issue 191: Bad bond directions in a branch") {
  // Only 1 of the double bonds has stereo defined!
  constexpr const char *smi = R"SMI(C2=NNC(N=C2)=N\N=C\c1ccccc1)SMI";
  constexpr const char *refSmi = R"SMI(C(=N\N=c1nccn[nH]1)/c1ccccc1)SMI";
  ROMol *mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(7)->getBondType() == Bond::DOUBLE);
  REQUIRE(mol->getBondWithIdx(7)->getStereo() == Bond::STEREOE);
  auto tmpSmi = MolToSmiles(*mol, 1);
  CHECK(tmpSmi == refSmi);
  delete mol;
  mol = SmilesToMol(tmpSmi);
  REQUIRE(mol);

  int numE = 0;
  for (auto bond : mol->bonds()) {
    if (bond->getBondType() == Bond::DOUBLE) {
      CHECK(bond->getStereo() != Bond::STEREOZ);
      if (bond->getStereo() == Bond::STEREOE) {
        ++numE;
      }
    }
  }
  CHECK(numE == 1);
  tmpSmi = MolToSmiles(*mol, 1);
  CHECK(tmpSmi == refSmi);
  delete mol;
}

TEST_CASE("Testing Issue 256: SMILES yields incorrect structure") {
  v2::SmilesParse::SmilesParserParams ps;
  ps.sanitize = false;
  {
    auto smi = "C1CC[C+]1=1CCC1";
    auto mol = v2::SmilesParse::MolFromSmiles(smi, ps);
    REQUIRE(mol);
    auto bond = mol->getBondBetweenAtoms(3, 0);
    REQUIRE(bond);
    REQUIRE(bond->getBondType() == Bond::SINGLE);
    bond = mol->getBondBetweenAtoms(3, 6);
    REQUIRE(bond);
    REQUIRE(bond->getBondType() == Bond::DOUBLE);
  }

  {
    auto smi = "C1CC[C+]=11CCC1";
    auto mol = v2::SmilesParse::MolFromSmiles(smi, ps);
    REQUIRE(mol);
    auto bond = mol->getBondBetweenAtoms(3, 0);
    REQUIRE(bond);
    REQUIRE(bond->getBondType() == Bond::DOUBLE);
    bond = mol->getBondBetweenAtoms(3, 6);
    REQUIRE(bond);
    REQUIRE(bond->getBondType() == Bond::SINGLE);
  }
}

TEST_CASE("Testing Issue 266: kekulized SMILES output") {
  RWMol *mol;
  std::string smi;

  smi = "c1ccccc1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi = MolToSmiles(*mol);
  REQUIRE(smi == "c1ccccc1");

  MolOps::Kekulize(*mol);
  smi = MolToSmiles(*mol);
  REQUIRE(smi == "C1=CC=CC=C1");
  delete mol;

  smi = "c1ccccc1c1ccccc1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi = MolToSmiles(*mol);
  REQUIRE(smi == "c1ccc(-c2ccccc2)cc1");

  MolOps::Kekulize(*mol);
  smi = MolToSmiles(*mol);
  REQUIRE(smi == "C1=CC=C(C2=CC=CC=C2)C=C1");
  delete mol;
}

TEST_CASE("Testing rootedAtAtom functionality") {
  {
    RWMol *mol;
    std::string smi;
    smi = "CN(C)C";
    mol = SmilesToMol(smi);
    REQUIRE(mol);
    smi = MolToSmiles(*mol, false, false, -1);
    REQUIRE(smi == "CN(C)C");
    smi = MolToSmiles(*mol, false, false, 1);
    REQUIRE(smi == "N(C)(C)C");
    smi = MolToSmiles(*mol, false, false, 2);
    REQUIRE(smi == "CN(C)C");
    delete mol;
  }
  {
    // This was github issue #182:
    RWMol mol;
    std::string smi;
    smi = MolToSmiles(mol);
    REQUIRE(smi == "");
    smi = MolToSmiles(mol, false, false, 0);
    REQUIRE(smi == "");
  }
}

TEST_CASE("Testing isotope handling") {
  {
    std::string smi = "C[13C](C)(C)C";
    RWMol *mol = SmilesToMol(smi);
    REQUIRE(mol);
    REQUIRE(feq(mol->getAtomWithIdx(1)->getMass(), 13.0034));
    smi = MolToSmiles(*mol, false);
    REQUIRE(smi == "CC(C)(C)C");
    smi = MolToSmiles(*mol, true);
    REQUIRE(smi == "C[13C](C)(C)C");
    delete mol;
  }
  {
    std::string smi = "C[12C](C)(C)C";
    RWMol *mol = SmilesToMol(smi);
    REQUIRE(mol);
    REQUIRE(mol->getAtomWithIdx(1)->getMass() == 12.0);
    smi = MolToSmiles(*mol, false);
    REQUIRE(smi == "CC(C)(C)C");
    smi = MolToSmiles(*mol, true);
    REQUIRE(smi == "C[12C](C)(C)C");
    delete mol;
  }
  {
    std::string smi = "CC[U]";
    RWMol *mol = SmilesToMol(smi);
    REQUIRE(mol);
    smi = MolToSmiles(*mol, false);
    REQUIRE(smi == "C[CH2][U]");
    smi = MolToSmiles(*mol, true);
    REQUIRE(smi == "C[CH2][U]");
    delete mol;
  }
  {
    std::string smi = "CC[238U]";
    RWMol *mol = SmilesToMol(smi);
    REQUIRE(mol);
    smi = MolToSmiles(*mol, false);
    REQUIRE(smi == "C[CH2][U]");
    smi = MolToSmiles(*mol, true);
    REQUIRE(smi == "C[CH2][238U]");
    delete mol;
  }
  {
    // issue 3526814
    std::string smi = "CCCCS(=[18O])(=O)CCCCl";
    RWMol *mol = SmilesToMol(smi);
    REQUIRE(mol);
    smi = MolToSmiles(*mol, false);
    REQUIRE(smi == "CCCCS(=O)(=O)CCCCl");
    smi = MolToSmiles(*mol, true);
    REQUIRE(smi == "CCCCS(=O)(=[18O])CCCCl");
    delete mol;
  }
  {
    // issue 3526814
    std::string smi = "CCCCS(=[24O])(=O)CCCCl";
    RWMol *mol = SmilesToMol(smi);
    REQUIRE(mol);
    smi = MolToSmiles(*mol, false);
    REQUIRE(smi == "CCCCS(=O)(=O)CCCCl");
    smi = MolToSmiles(*mol, true);
    REQUIRE(smi == "CCCCS(=O)(=[24O])CCCCl");
    delete mol;
  }
  {
    // issue 3526814
    std::string smi = "CCCCS(=O)(=[24O])CCCCl";
    RWMol *mol = SmilesToMol(smi);
    REQUIRE(mol);
    smi = MolToSmiles(*mol, false);
    REQUIRE(smi == "CCCCS(=O)(=O)CCCCl");
    smi = MolToSmiles(*mol, true);
    REQUIRE(smi == "CCCCS(=O)(=[24O])CCCCl");
    delete mol;
  }
}

TEST_CASE("Testing SF.net bug 1670149") {
  RWMol *mol;
  std::string smi;

  smi = "C1[NH2+]CCC1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  REQUIRE(smi == "C1CC[NH2+]C1");

  mol->getAtomWithIdx(1)->setNumExplicitHs(0);
  mol->getAtomWithIdx(1)->setNoImplicit(false);
  mol->getAtomWithIdx(1)->updatePropertyCache();
  REQUIRE(mol->getAtomWithIdx(1)->getNumImplicitHs() == 2);
  smi = MolToSmiles(*mol, false, false, -1);
  REQUIRE(smi == "C1CC[NH2+]C1");
  delete mol;
}

TEST_CASE("Testing SF.net bug 1719046: explicit Hs in canonical smiles") {
  RWMol *mol;
  std::string smi;

  smi = "Cl[CH]1CCCCC1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  REQUIRE(smi == "ClC1CCCCC1");

  delete mol;
  smi = "Cl[C@H]1CCCCC1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  REQUIRE(smi == "ClC1CCCCC1");

  delete mol;
  smi = "Cl[C@H]1C(Br)CCCC1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  REQUIRE(smi == "ClC1CCCCC1Br");

  delete mol;
  smi = "[CH]1=[CH][CH]=[CH][CH]=[CH]1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  REQUIRE(smi == "c1ccccc1");

  delete mol;
  smi = "c1ccccn1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  REQUIRE(smi == "c1ccncc1");

  delete mol;
  smi = "C1=CNC=C1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  REQUIRE(smi == "c1cc[nH]c1");

  delete mol;
  smi = "[CH]1=[CH][NH][CH]=[CH]1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  REQUIRE(smi == "c1cc[nH]c1");

  delete mol;
  // this was Issue 35525671
  smi = "P1C=CC=C1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  REQUIRE(smi == "c1cc[pH]c1");

  delete mol;
}

TEST_CASE("Testing SF.net bug 1842174: bad bond dirs in branches") {
  RWMol *mol;
  std::string smi;

  smi = "F/C=N/Cl";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi = MolToSmiles(*mol, true, false, -1);
  CHECK(smi == "F/C=N/Cl");

  smi = MolToSmiles(*mol, true, false, 1);
  CHECK(smi == R"SMI(C(/F)=N\Cl)SMI");

  delete mol;
  smi = "C(\\C=C\\F)=C(/Cl)Br";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi = MolToSmiles(*mol, true, false, -1);
  CHECK(smi == "F/C=C/C=C(/Cl)Br");

  smi = MolToSmiles(*mol, true, false, 0);
  CHECK(smi == "C(/C=C/F)=C(\\Cl)Br");
  delete mol;

  smi = "O=NC1=NOC(=N\\O)/C1=N\\O";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi = MolToSmiles(*mol, true, false, -1);
  CHECK(smi == "O=NC1=NOC(=N\\O)/C1=N\\O");

  // ----------------------
  //  the next two examples are a pair:
  // vvvvvvvvvvvvvvvvvvvvvv
  delete mol;
  smi = "O/N=C/1COCC1=N\\O";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi = MolToSmiles(*mol, true, false, -1);
  CHECK(smi == R"SMI(O/N=C1COCC\1=N\O)SMI");

  // this time the algorithm is forced to set
  // the directionality on the ring closure bond:
  delete mol;
  smi = "O/N=C/1COC[N+]1=N\\O";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi = MolToSmiles(*mol, true, false, -1);
  CHECK(smi == R"SMI(O/N=C1COC[N+]\1=N\O)SMI");
  // ^^^^^^^^^^^^^^^^^^^^^^
  // end of the pair
  // ----------------------

  delete mol;
}

TEST_CASE(
    "Testing SF.net bug 1844617: oscillating chirality in canonical smiles") {
  RWMol *mol;
  std::string smi, smi2;
  std::string label;

  smi = "O=C1CC[C@@]2(O)[C@@H]3N(C)CC[C@]22[C@H]1OC[C@H]2CC3";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  REQUIRE(mol->getAtomWithIdx(6)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(6)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "R");
  REQUIRE(mol->getAtomWithIdx(11)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(11)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  REQUIRE(mol->getAtomWithIdx(12)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(12)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "R");
  REQUIRE(mol->getAtomWithIdx(15)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(15)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi2 = MolToSmiles(*mol, true);
  REQUIRE(smi == smi2);

  delete mol;
  smi = "O=C1CC[C@@]2(O)[C@@H]3N(C)CC[C@]22[C@H]1OC[C@H]2CC3";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  REQUIRE(mol->getAtomWithIdx(6)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(6)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "R");
  REQUIRE(mol->getAtomWithIdx(11)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(11)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  REQUIRE(mol->getAtomWithIdx(12)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(12)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "R");
  REQUIRE(mol->getAtomWithIdx(15)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(15)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  smi = MolToSmiles(*mol, true, false, 0);
  delete mol;
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi2 = MolToSmiles(*mol, true, false, 0);
  REQUIRE(smi == smi2);

  delete mol;
  smi = "O=C1CC[C@@]2(O)[C@@H]3N(CC4CC4)CC[C@]22[C@H]1OC[C@H]2CC3";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  REQUIRE(mol->getAtomWithIdx(6)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(6)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "R");
  REQUIRE(mol->getAtomWithIdx(14)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(14)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  REQUIRE(mol->getAtomWithIdx(15)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(15)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "R");
  REQUIRE(mol->getAtomWithIdx(18)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(18)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  smi2 = MolToSmiles(*mol, true);
  REQUIRE(smi == smi2);
  delete mol;
}

TEST_CASE("Testing SF.net bug 1844959: bad handling of Hs in chiral smiles") {
  RWMol *mol;
  std::string smi, smi2;
  std::string label;

  // ----------------------
  //  the next examples are a set:
  //  (this is the part that was originally working):
  // vvvvvvvvvvvvvvvvvvvvvv
  smi = "C[C@]12CNOC2.F1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "R");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "R");
  smi2 = MolToSmiles(*mol, true);
  REQUIRE(smi == smi2);

  // swap the order and make sure the chirality swaps with it:
  delete mol;
  smi = "C[C@]12CNOC1.F2";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  smi2 = MolToSmiles(*mol, true);
  REQUIRE(smi == smi2);
  delete mol;

  // now make sure it works with a reversed chiral tag:
  smi = "C[C@@]12CNOC2.F1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  smi2 = MolToSmiles(*mol, true);
  REQUIRE(smi == smi2);
  delete mol;
  smi = "C[C@@]12CNOC1.F2";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "R");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "R");
  smi2 = MolToSmiles(*mol, true);
  REQUIRE(smi == smi2);
  delete mol;
  // ^^^^^^^^^^^^^^^^^^^^^^
  // end of the set
  // ----------------------

  // ----------------------
  //  the next examples are a set:
  //  (this is the part that was originally failing):
  // vvvvvvvvvvvvvvvvvvvvvv
  smi = "C[C@]12CNOC2.[H]1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  smi2 = MolToSmiles(*mol, true);
  REQUIRE(smi == smi2);

  // swap the order and make sure the chirality swaps with it:
  delete mol;
  smi = "C[C@]12CNOC1.[H]2";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "R");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "R");
  smi2 = MolToSmiles(*mol, true);
  REQUIRE(smi == smi2);
  delete mol;

  // now make sure it works with a reversed chiral tag:
  smi = "C[C@@]12CNOC2.[H]1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "R");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "R");
  smi2 = MolToSmiles(*mol, true);
  REQUIRE(smi == smi2);
  delete mol;
  smi = "C[C@@]12CNOC1.[H]2";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  MolOps::assignStereochemistry(*mol);
  REQUIRE(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  REQUIRE(label == "S");
  smi2 = MolToSmiles(*mol, true);
  REQUIRE(smi == smi2);
  // ^^^^^^^^^^^^^^^^^^^^^^
  // end of the set
  // ----------------------

  delete mol;
}

TEST_CASE("Testing sf.net bug 1942220") {
  RWMol *m;
  std::string smi;

  smi = "[C](Cl)Br";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 3);
  REQUIRE(m->getNumAtoms(false) == 3);
  smi = MolToSmiles(*m);
  REQUIRE(smi == "Cl[C]Br");

  delete m;
  smi = "[CH2](Cl)Br";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 3);
  REQUIRE(m->getNumAtoms(false) == 5);
  smi = MolToSmiles(*m);
  REQUIRE(smi == "ClCBr");

  delete m;
  smi = "C(Cl)Br";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 3);
  REQUIRE(m->getNumAtoms(false) == 5);
  smi = MolToSmiles(*m);
  REQUIRE(smi == "ClCBr");

  delete m;
  smi = "OS(=O)=O";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  smi = MolToSmiles(*m);
  REQUIRE(smi == "O=[SH](=O)O");

  delete m;
}

TEST_CASE("Testing error reporting with ring stereochem") {
  RWMol *m;
  std::string smi;

  smi = "C[C@H]1CC[C@@H](C)CC1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 8);

  smi = MolToSmiles(*m, false);
  REQUIRE((!m->hasProp(common_properties::_ringStereoWarning)));

  delete m;
}

TEST_CASE("Testing sf.net issue 3127883 (kekulization failing)") {
  {
    ROMol *m;
    std::string smi;
    smi = "c(:c:c:1):c:c:c:1";
    m = SmilesToMol(smi);
    REQUIRE(m);
    delete m;
  }
  {
    ROMol *m;
    std::string smi;
    smi = "c1(:c(:c(:c(-C(-c2:c(:c(:c(:c(:c:2)))))=C):c(:c:1))))";
    m = SmilesToMol(smi);
    REQUIRE(m);
    delete m;
  }
}

TEST_CASE("Issue 3139534: stereochemistry in larger rings") {
  // the parsing part of this is in ../testChirality.cpp, here we look at
  // smiles generation

  {
    RWMol *m;
    std::string smiles = "C1COC/C=C\\CCC1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(4)->getStereo() == Bond::STEREOZ);

    smiles = MolToSmiles(*m, true);
    CHECK(smiles == "C1=C\\COCCCCC/1");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C1COC/C=C/CCC1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(4)->getStereo() == Bond::STEREOE);

    smiles = MolToSmiles(*m, true);
    CHECK(smiles == "C1=C/COCCCCC/1");

    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C1CC/C=C/C=C/CCC1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    smiles = MolToSmiles(*m, true, false, -1, false);
    CHECK(smiles == "C1CC/C=C/C=C/CCC1");

    smiles = MolToSmiles(*m, true);
    CHECK(smiles == "C1=C/CCCCCC/C=C/1");
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "C/1=C/C=C/CCCCCC1";
    m = SmilesToMol(smiles);
    REQUIRE(m);

    smiles = MolToSmiles(*m, true);
    CHECK(smiles == "C1=C\\CCCCCC/C=C/1");
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "C1COC/C=C/C=C/C1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(4)->getStereo() == Bond::STEREOE);

    smiles = MolToSmiles(*m, true);
    CHECK(smiles == "C1=C/CCCOC/C=C/1");

    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "C1=C/OCC/C=C\\CC\\1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    CHECK(m->getBondWithIdx(0)->getStereo() == Bond::STEREOZ);
    CHECK(m->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "C1CCCCN/C=C/1";
    m = SmilesToMol(smiles);
    REQUIRE(m);

    smiles = MolToSmiles(*m, true, false, 7, false);
    CHECK(smiles == "C1=C/NCCCCC/1");

    smiles = MolToSmiles(*m, true, false, 0, false);
    CHECK(smiles == "C1CCCCN/C=C/1");

    delete m;
  }

  {
    RWMol *m;
    // the 2 initial directed bonds are redundant (/bad ??)
    std::string smiles = "CCC/[N+]/1=C/c2ccccc2OC(=O)/C=C1/O";
    m = SmilesToMol(smiles);
    REQUIRE(m);

    REQUIRE(m->getBondWithIdx(3)->getStereo() == Bond::STEREOZ);
    REQUIRE(m->getBondWithIdx(14)->getStereo() == Bond::STEREOE);

    smiles = MolToSmiles(*m, true);
    CHECK(smiles == R"(CCC[N+]1=C/c2ccccc2OC(=O)/C=C\1O)");

    delete m;

    // 2nd pass to check stability
    m = SmilesToMol(smiles);
    REQUIRE(m);

    REQUIRE(m->getBondWithIdx(3)->getStereo() == Bond::STEREOZ);
    REQUIRE(m->getBondWithIdx(14)->getStereo() == Bond::STEREOE);

    smiles = MolToSmiles(*m, true);
    CHECK(smiles == R"(CCC[N+]1=C/c2ccccc2OC(=O)/C=C\1O)");

    delete m;
  }

  {  // Github #2023
    RWMol *m;
    // the initial directed bond is redundant
    std::string smiles = R"(CO/C1=C/C=C\C=C/C=N\1)";
    m = SmilesToMol(smiles);
    REQUIRE(m);

    REQUIRE(m->getBondWithIdx(2)->getStereo() == Bond::STEREOE);
    REQUIRE(m->getBondWithIdx(4)->getStereo() == Bond::STEREOZ);
    REQUIRE(m->getBondWithIdx(6)->getStereo() == Bond::STEREOZ);
    REQUIRE(m->getBondWithIdx(8)->getStereo() == Bond::STEREOZ);

    smiles = MolToSmiles(*m, true);
    CHECK(smiles == R"(COC1=C/C=C\C=C/C=N\1)");

    delete m;

    // 2nd pass to check stability
    m = SmilesToMol(smiles);
    REQUIRE(m);

    REQUIRE(m->getBondWithIdx(2)->getStereo() == Bond::STEREOE);
    REQUIRE(m->getBondWithIdx(4)->getStereo() == Bond::STEREOZ);
    REQUIRE(m->getBondWithIdx(6)->getStereo() == Bond::STEREOZ);
    REQUIRE(m->getBondWithIdx(8)->getStereo() == Bond::STEREOZ);

    smiles = MolToSmiles(*m, true);
    CHECK(smiles == R"(COC1=C/C=C\C=C/C=N\1)");

    delete m;
  }

  // some torture tests with natural products (thanks to James Davidson for the
  // examples)
  {
    RWMol *m;
    std::string smiles =
        "NC(=O)O[C@H]1C(/C)=C/[C@H](C)[C@@H](O)[C@@H](OC)C[C@H](C)C\\C2=C(/"
        "OC)C(=O)\\C=C(\\NC(=O)C(\\C)=C\\C=C/[C@@H]1OC)C2=O";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(30, 32)->getStereo() == Bond::STEREOE);
    REQUIRE(m->getBondBetweenAtoms(33, 34)->getStereo() == Bond::STEREOZ);
    REQUIRE(m->getBondBetweenAtoms(5, 7)->getStereo() == Bond::STEREOE);

    std::string csmiles = MolToSmiles(*m, true);

    RWMol *m2;
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      std::string nsmiles = MolToSmiles(*m, true, false, i, false);
      m2 = SmilesToMol(nsmiles);
      REQUIRE(m2);
      std::string ncsmiles = MolToSmiles(*m2, true);
      CHECK(ncsmiles == csmiles);
      delete m2;
    }
    delete m;
  }

  {
    RWMol *m;
    std::string smiles =
        "CC(O[C@@H]1C=C(C)[C@H]2[C@H]([C@H]3O[C@@H]2C/"
        "C(C)=C\\CC[C@@]3(C)OC(C)=O)[C@H]1C(OC(C)=O)(C)C)=O";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(13, 15)->getStereo() == Bond::STEREOZ);

    std::string csmiles = MolToSmiles(*m, true);

    RWMol *m2;
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      std::string nsmiles = MolToSmiles(*m, true, false, i, false);
      m2 = SmilesToMol(nsmiles);
      REQUIRE(m2);
      std::string ncsmiles = MolToSmiles(*m2, true);
      CHECK(ncsmiles == csmiles);
      delete m2;
    }
    delete m;
  }

  {
    RWMol *m;
    std::string smiles =
        "CC(O[C@@H]1C=C(C)[C@H]2[C@H]([C@H]3O[C@@H]2C/C(C)=C/"
        "CC[C@@]3(C)OC(C)=O)[C@H]1C(OC(C)=O)(C)C)=O";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(13, 15)->getStereo() == Bond::STEREOE);

    std::string csmiles = MolToSmiles(*m, true);

    RWMol *m2;
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      std::string nsmiles = MolToSmiles(*m, true, false, i, false);
      m2 = SmilesToMol(nsmiles);
      REQUIRE(m2);
      std::string ncsmiles = MolToSmiles(*m2, true);
      CHECK(ncsmiles == csmiles);
      delete m2;
    }
    delete m;
  }

  {
    RWMol *m;
    std::string smiles =
        "CC(=O)[C@@H]1CC=C(C)[C@@H]2[C@@H]3O[C@@H]([C@@H](O)C/"
        "C=C\\CC3)[C@@H]12";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(15, 16)->getStereo() == Bond::STEREOZ);

    std::string csmiles = MolToSmiles(*m, true);

    RWMol *m2;
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      std::string nsmiles = MolToSmiles(*m, true, false, i, false);
      m2 = SmilesToMol(nsmiles);
      REQUIRE(m2);
      std::string ncsmiles = MolToSmiles(*m2, true);
      CHECK(ncsmiles == csmiles);
      delete m2;
    }
    delete m;
  }
}

TEST_CASE("test adding atom-map information") {
  RWMol *m;
  std::string smiles = "[*:1]CCC([C:200])C";
  m = SmilesToMol(smiles);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(0)->hasProp(common_properties::molAtomMapNumber));

  // changed: smiles does not need to be canonical
  smiles = MolToSmiles(*m, true, false, -1, false);
  CHECK(smiles == "[*:1]CCC([C:200])C");

  delete m;
}

TEST_CASE("Issue 3145697 repeated ring labels in disconnected structures") {
  {
    RWMol *m;
    std::string smiles = "C1.C11.C1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    smiles = MolToSmiles(*m, true);
    REQUIRE(smiles == "CCC");
    delete m;

    smiles = "C1.C11.C";
    m = SmilesToMol(smiles);
    REQUIRE(!m);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C1.C11.O1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    smiles = MolToSmiles(*m, true);
    REQUIRE(smiles == "CCO");
    delete m;

    smiles = "C1.C1=1.O1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    smiles = MolToSmiles(*m, true);
    REQUIRE(smiles == "CC=O");
    delete m;

    smiles = "C1.C=11.O1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    smiles = MolToSmiles(*m, true);
    REQUIRE(smiles == "C=CO");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C1C.CC11CCC1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    smiles = MolToSmiles(*m, true);
    REQUIRE(smiles == "CCC1(C)CCC1");
    delete m;
    smiles = "C1C.CC11CCC";
    m = SmilesToMol(smiles);
    REQUIRE(!m);
  }
}

TEST_CASE("Issue 3152751 cannot roundtrip charged aromatic Se and Te") {
  {
    RWMol *m;
    std::string smiles = "c1cccc[te+]1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
    smiles = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "c1cccc[se+]1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
    smiles = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "c1ccc[te]1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
    smiles = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "c1ccc[se]1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
    smiles = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
    delete m;
  }
}

TEST_CASE("Testing use of replacement patterns in input") {
  {
    std::string smi = "C{cycloprop}C";
    std::map<std::string, std::string> repls;
    repls["{cycloprop}"] = "C1(CC1)";
    RWMol *mol = SmilesToMol(smi, 0, true, &repls);
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 5);
    REQUIRE(mol->getAtomWithIdx(1)->getDegree() == 4);
    delete mol;
  }

  {
    std::string smi = "C{cycloprop}C";
    std::map<std::string, std::string> repls;
    repls["{cycloprop}"] = "C1(C({acid})C1)";
    repls["{acid}"] = "C(=O)O";
    RWMol *mol = SmilesToMol(smi, 0, true, &repls);
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 8);
    REQUIRE(mol->getAtomWithIdx(1)->getDegree() == 4);
    delete mol;
  }
}

TEST_CASE("Testing forcing explicit bonds in the output SMILES") {
  {
    std::string smi = "CCC";
    RWMol *mol = SmilesToMol(smi);
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 3);
    smi = MolToSmiles(*mol, true);
    REQUIRE(smi == "CCC");
    smi = MolToSmiles(*mol, true, false, -1, true, true);
    REQUIRE(smi == "C-C-C");

    delete mol;
  }
  {
    std::string smi = "C1CC1";
    RWMol *mol = SmilesToMol(smi);
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 3);
    smi = MolToSmiles(*mol, true);
    REQUIRE(smi == "C1CC1");
    smi = MolToSmiles(*mol, true, false, -1, true, true);
    REQUIRE(smi == "C1-C-C-1");

    delete mol;
  }
  {
    std::string smi = "c1ccccc1";
    RWMol *mol = SmilesToMol(smi);
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 6);
    smi = MolToSmiles(*mol, true);
    REQUIRE(smi == "c1ccccc1");
    smi = MolToSmiles(*mol, true, false, -1, true, true);
    REQUIRE(smi == "c1:c:c:c:c:c:1");

    delete mol;
  }
  {
    std::string smi = "c1ccccc1c1ccccc1";
    RWMol *mol = SmilesToMol(smi);
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 12);
    smi = MolToSmiles(*mol, true);
    REQUIRE(smi == "c1ccc(-c2ccccc2)cc1");
    smi = MolToSmiles(*mol, true, false, -1, true, true);
    REQUIRE(smi == "c1:c:c:c(-c2:c:c:c:c:c:2):c:c:1");

    delete mol;
  }
}

TEST_CASE("Issue 3525799: bad smiles for r groups") {
  {
    RWMol *m;
    std::string smiles = "CC*";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    smiles = MolToSmiles(*m, true);
    REQUIRE(smiles == "*CC");
    m->getAtomWithIdx(2)->setProp(common_properties::dummyLabel, "foo");
    smiles = MolToSmiles(*m, true);
    REQUIRE(smiles == "*CC");
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "CC*";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    smiles = MolToSmiles(*m, true);
    REQUIRE(smiles == "*CC");
    m->getAtomWithIdx(2)->setProp(common_properties::smilesSymbol, "Xa");
    smiles = MolToSmiles(*m, true);
    REQUIRE(smiles == "[Xa]CC");
    delete m;
  }
}

TEST_CASE("Issue 3526810: canonical smiles failure in symmetric heterocycles") {
  {
    RWMol *m;
    std::string smiles = "C1SCCSCCCSCCSCC1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmiles1 = MolToSmiles(*m, true);
    delete m;
    std::string smiles2 = "C1CSCCSCCCSCCSC1";
    m = SmilesToMol(smiles2);
    REQUIRE(m);
    std::string csmiles2 = MolToSmiles(*m, true);
    delete m;

    REQUIRE(csmiles1 == csmiles2);
  }

  {
    RWMol *m;
    std::string smiles = "C1NCCNCCCNCCNCC1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmiles1 = MolToSmiles(*m, true);
    delete m;
    std::string smiles2 = "C1CNCCNCCCNCCNC1";
    m = SmilesToMol(smiles2);
    REQUIRE(m);
    std::string csmiles2 = MolToSmiles(*m, true);
    delete m;

    REQUIRE(csmiles1 == csmiles2);
  }

  {
    RWMol *m;
    std::string smiles = "C1CNCCCNCCNCCCNC1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmiles1 = MolToSmiles(*m, true);
    delete m;
    std::string smiles2 = "C1CCNCCCNCCNCCCN1";
    m = SmilesToMol(smiles2);
    REQUIRE(m);
    std::string csmiles2 = MolToSmiles(*m, true);
    delete m;

    REQUIRE(csmiles1 == csmiles2);
  }
}

TEST_CASE(
    "Issue 3526815: canonical smiles failure in many symmetric fragments") {
  RWMol *m;
  std::string smiles =
      "O.O.O.O.O.O.O.O.O.[Pd].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+]"
      ".[Na+].[O-]S(=O)(=O)c1cccc(c1)P(c1cccc(c1)S(=O)(=O)[O-])c1cccc(c1)S(="
      "O)(=O)[O-].[O-]S(=O)(=O)c1cccc(c1)P(c1cccc(c1)S(=O)(=O)[O-])c1cccc(c1)"
      "S(=O)(=O)[O-].[O-]S(=O)(=O)c1cccc(c1)P(c1cccc(c1)S(=O)(=O)[O-])c1cccc("
      "c1)S(=O)(=O)[O-]";
  m = SmilesToMol(smiles);
  REQUIRE(m);
  std::string csmiles1 = MolToSmiles(*m, true);
  delete m;
  std::string smiles2 =
      "O.O.O.O.O.O.O.O.O.[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+"
      "].[O-]S(c1cccc(P(c2cccc(S([O-])(=O)=O)c2)c2cccc(S([O-])(=O)=O)c2)c1)(="
      "O)=O.[Pd].[O-]S(=O)(=O)c1cccc(P(c2cccc(S([O-])(=O)=O)c2)c2cccc(S([O-])"
      "(=O)=O)c2)c1.[O-]S(=O)(=O)c1cccc(P(c2cccc(S([O-])(=O)=O)c2)c2cccc(S(["
      "O-])(=O)=O)c2)c1";
  m = SmilesToMol(smiles2);
  REQUIRE(m);
  std::string csmiles2 = MolToSmiles(*m, true);
  delete m;

  REQUIRE(csmiles1 == csmiles2);
}

TEST_CASE("Testing Fragment Smiles") {
  {
    RWMol *m;
    std::string smiles = "OCCCC";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {0, 1, 2};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse);
    REQUIRE(csmiles == "CCO");
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "OCCCCCC";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {0, 1, 2, 3};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse);
    REQUIRE(csmiles == "CCCO");
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "OC1CC1CCC";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {1, 2, 3};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse);
    REQUIRE(csmiles == "C1CC1");
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "OC1CC1CCC";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {1, 2, 3};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    int bs[] = {1, 2, 6};
    std::vector<int> bondsToUse(bs, bs + sizeof(bs) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse, &bondsToUse);
    REQUIRE(csmiles == "C1CC1");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "OC1CC1CCC";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {1, 2, 3};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    int bs[] = {1, 2};
    std::vector<int> bondsToUse(bs, bs + sizeof(bs) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse, &bondsToUse);
    REQUIRE(csmiles == "CCC");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "OC1CCCCC1N";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {1, 2, 3, 4, 5, 6};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse);
    REQUIRE(csmiles == "C1CCCCC1");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "OCCCCCCN";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {1, 2, 3, 4, 5, 6};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse);
    REQUIRE(csmiles == "CCCCCC");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "OCCCCCCN";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {1, 2, 3, 4, 5, 6};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    int bs[] = {1, 2, 3, 4, 5};
    std::vector<int> bondsToUse(bs, bs + sizeof(bs) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse, &bondsToUse);

    REQUIRE(csmiles == "CCCCCC");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "OC1CCCCC1N";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {1, 2, 3, 4, 5, 6};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    int bs[] = {1, 2, 3, 4, 5};
    std::vector<int> bondsToUse(bs, bs + sizeof(bs) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse, &bondsToUse);
    REQUIRE(csmiles == "CCCCCC");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "Oc1ccccc1N";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {1, 2, 3, 4, 5, 6};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse);
    REQUIRE(csmiles == "c1ccccc1");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "Oc1ccccc1N";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {1, 2, 3, 4, 5, 6};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    int bs[] = {1, 2, 3, 4, 5};
    std::vector<int> bondsToUse(bs, bs + sizeof(bs) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse, &bondsToUse);
    REQUIRE(csmiles == "cccccc");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "OCCCC";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {0, 1, 2};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string labels[5] = {"[A]", "[B]", "[B]", "", ""};
    std::vector<std::string> atomLabels(labels, labels + 5);
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, &atomLabels);
    REQUIRE(csmiles == "[A][B][B]");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CCCCO";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {2, 3, 4};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string labels[5] = {"", "", "[B]", "[B]", "[A]"};
    std::vector<std::string> atomLabels(labels, labels + 5);
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, &atomLabels);
    REQUIRE(csmiles == "[A][B][B]");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CCCCO";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {2, 3, 4};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string labels[5] = {"", "", "[B]", "[A]", "[B]"};
    std::vector<std::string> atomLabels(labels, labels + 5);
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, &atomLabels);
    REQUIRE(csmiles == "[B][A][B]");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CC(=O)OCC";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {0, 1, 2, 3};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, nullptr, nullptr);
    REQUIRE(csmiles == "CC(=O)O");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CC(=O)OCC";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {0, 1, 2, 3};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string labels[5] = {"-", "=", "-", "", ""};
    std::vector<std::string> bondLabels(labels, labels + 5);
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, nullptr, &bondLabels);
    REQUIRE(csmiles == "C-C(=O)-O");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CC(=O)OCC";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {0, 1, 2, 3};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string labels[5] = {"a", "b", "a", "", ""};
    std::vector<std::string> bondLabels(labels, labels + 5);
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, nullptr, &bondLabels);
    REQUIRE(csmiles == "CaC(bO)aO");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CC(=CC)CCC";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {0, 1, 2, 4};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse);
    REQUIRE(csmiles == "C=C(C)C");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CC(=CC)CCC";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {0, 1, 2, 4};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string labels[6] = {"a", "b", "", "a", "", ""};
    std::vector<std::string> bondLabels(labels, labels + 6);
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, nullptr, &bondLabels);
    REQUIRE(csmiles == "CbC(aC)aC");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CC(=CC)CCC";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {0, 1, 2, 4};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string labels[6] = {"b", "a", "", "a", "", ""};
    std::vector<std::string> bondLabels(labels, labels + 6);
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, nullptr, &bondLabels);
    REQUIRE(csmiles == "CaC(aC)bC");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CC(=CC)CCC";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {0, 1, 2, 4};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string labels[6] = {"b", "b", "", "a", "", ""};
    std::vector<std::string> bondLabels(labels, labels + 6);
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, nullptr, &bondLabels);
    REQUIRE(csmiles == "CbC(aC)bC");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "OC1CC1CC";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    int as[] = {0, 4};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse, nullptr, nullptr,
                                              nullptr, false, false, -1, false);
    REQUIRE(csmiles == "O.C");
    delete m;
  }
}

TEST_CASE("Issue 3528556: canonical smiles failure in cycle") {
  RWMol *m;
  std::string smiles = "N12.N13.C24.C35.C46.C56";
  m = SmilesToMol(smiles);
  REQUIRE(m);
  std::string csmiles1 = MolToSmiles(*m, true);
  delete m;
  std::string smiles2 = "N1NCCCC1";
  m = SmilesToMol(smiles2);
  REQUIRE(m);
  std::string csmiles2 = MolToSmiles(*m, true);
  delete m;
  REQUIRE(csmiles1 == csmiles2);
}

TEST_CASE("Issue 253: do not repeat ring closure digits on the same atom") {
  RWMol *m;
  std::string smiles = "C1CCCC1CCC1CCCCC11CCCCC1";
  m = SmilesToMol(smiles);
  REQUIRE(m);
  std::string csmiles1 = MolToSmiles(*m, true);
  REQUIRE(csmiles1 == "C1CCC2(CC1)CCCCC2CCC1CCCC1");
  delete m;
}

TEST_CASE("Issue 257: unrecognized bonds are in SMILES as ?s") {
  RWMol *m;
  std::string smiles = "CCO";
  m = SmilesToMol(smiles);
  REQUIRE(m);
  m->getBondWithIdx(1)->setBondType(Bond::UNSPECIFIED);
  std::string csmiles = MolToSmiles(*m);
  REQUIRE(csmiles == "CC~O");
  delete m;
  m = SmilesToMol(csmiles);
  REQUIRE(m);
  REQUIRE(m->getBondWithIdx(1)->getBondType() == Bond::UNSPECIFIED);
  delete m;
}

TEST_CASE("Testing Github 12: non-canonical fragment smiles") {
  RWMol *m;
  std::string smiles = "c1c(C)cccc1";
  m = SmilesToMol(smiles);
  REQUIRE(m);
  int as[] = {0, 1, 2};
  std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
  std::string csmiles1 = MolFragmentToSmiles(*m, atomsToUse);
  int as2[] = {1, 2, 3};
  std::vector<int> atomsToUse2(as2, as2 + sizeof(as2) / sizeof(int));
  std::string csmiles2 = MolFragmentToSmiles(*m, atomsToUse2);
  REQUIRE(csmiles1 == csmiles2);
  delete m;
}

TEST_CASE("Testing handling of ring stereochemistry") {
  // examples from TJ O'Donnell
  std::string inSmiles[] = {
      "N#Cc1ccc2[nH]cc([C@H]3CC[C@@H](N4CCN(c5cccc6nccnc65)CC4)CC3)c2c1",
      "N#Cc1ccc2[nH]cc([C@@H]3CC[C@H](N4CCN(c5cccc6nccnc65)CC4)CC3)c2c1",
      "N#Cc1ccc2[nH]cc([C@@H]3CC[C@@H](N4CCN(c5cccc6nccnc65)CC4)CC3)c2c1",
      "N#Cc1ccc2[nH]cc([C@H]3CC[C@H](N4CCN(c5cccc6nccnc65)CC4)CC3)c2c1",
      "N#Cc1ccc2[nH]cc([C@@H]3CC[C@@H](N4CCN(c5cccc6nccnc65)CC4)CC3)c2c1",
      "N#Cc1ccc2[nH]cc([C@H]3CC[C@H](N4CCN(c5cccc6nccnc65)CC4)CC3)c2c1",
      "O=C(N[C@@H]1CC[C@@H](CCN2CCN(c3cccc(Cl)c3Cl)CC2)CC1)c1cccs1",
      "O=C(N[C@H]1CC[C@H](CCN2CCN(c3cccc(Cl)c3Cl)CC2)CC1)c1cccs1",
      "Cn1ccc2ccc3c4[nH]c5c(cccc5CCN[C@H]5CC[C@H](O)CC5)c4c4c(c3c21)C(=O)NC4="
      "O",
      "Cn1ccc2ccc3c4[nH]c5c(cccc5CCN[C@@H]5CC[C@@H](O)CC5)c4c4c(c3c21)C(=O)"
      "NC4=O",
      "N=C(N)Nc1ccc(CNC(=O)N2CCN(C(=O)O[C@@H]3CCC[C@H](OC(=O)N4CCN(C(=O)"
      "CCCCn5ccnc5)CC4)CCC3)CC2)cc1",
      "N=C(N)Nc1ccc(CNC(=O)N2CCN(C(=O)O[C@H]3CCC[C@@H](OC(=O)N4CCN(C(=O)"
      "CCCCn5ccnc5)CC4)CCC3)CC2)cc1",
      "CC(C)c1cc(C(C)C)c(S(=O)(=O)NC[C@H]2CC[C@H](C(=O)NNC(=O)c3cc4ccccc4s3)"
      "CC2)c(C(C)C)c1",
      "CC(C)c1cc(C(C)C)c(S(=O)(=O)NC[C@@H]2CC[C@@H](C(=O)NNC(=O)"
      "c3cc4ccccc4s3)CC2)c(C(C)C)c1",
      "O=C(CCC[C@@H]1OO[C@H](CCCC(=O)c2ccccc2)OO1)c1ccccc1",
      "O=C(CCC[C@@H]1OO[C@H](CCCC(=O)c2ccccc2)OO1)c1ccccc1",
      "O=C(CCC[C@@H]1OO[C@@H](CCCC(=O)c2ccccc2)OO1)c1ccccc1",
      "O=C(CCC[C@H]1OO[C@H](CCCC(=O)c2ccccc2)OO1)c1ccccc1",
      "CCCn1c2[nH]c([C@@H]3CC[C@@H](CNC(C)=O)CC3)nc2c(=O)n(CCC)c1=O",
      "CCCn1c2[nH]c([C@H]3CC[C@H](CNC(C)=O)CC3)nc2c(=O)n(CCC)c1=O",
      "c1cc2c(cccc2N2CCN([C@H]3CC[C@@H](c4c[nH]c5ccccc54)CC3)CC2)[nH]1",
      "c1cc2c(cccc2N2CCN([C@@H]3CC[C@H](c4c[nH]c5ccccc54)CC3)CC2)[nH]1",
      "c1cc2c(cccc2N2CCN([C@H]3CC[C@@H](c4c[nH]c5ccccc54)CC3)CC2)[nH]1",
      "c1cc2c(cccc2N2CCN([C@@H]3CC[C@H](c4c[nH]c5ccccc54)CC3)CC2)[nH]1",
      "c1cc2c(cccc2N2CCN([C@@H]3CC[C@@H](c4c[nH]c5ccccc54)CC3)CC2)[nH]1",
      "c1cc2c(cccc2N2CCN([C@H]3CC[C@H](c4c[nH]c5ccccc54)CC3)CC2)[nH]1",
      "c1cc2c(cccc2N2CCN([C@@H]3CC[C@@H](c4c[nH]c5ccccc54)CC3)CC2)[nH]1",
      "c1cc2c(cccc2N2CCN([C@H]3CC[C@H](c4c[nH]c5ccccc54)CC3)CC2)[nH]1",
      "CCCCC(=O)N[C@@]1(C(=O)N[C@H](Cc2ccccc2)C(=O)N[C@@H](CCCN=C(N)N)C(=O)N["
      "C@@H](Cc2c[nH]c3ccccc23)C(=O)NCC(N)=O)CC[C@@H](c2ccc(C)cc2)CC1",
      "CCCCC(=O)N[C@]1(C(=O)N[C@H](Cc2ccccc2)C(=O)N[C@@H](CCCN=C(N)N)C(=O)N["
      "C@@H](Cc2c[nH]c3ccccc23)C(=O)NCC(N)=O)CC[C@H](c2ccc(C)cc2)CC1",
      "CC(C)Oc1ccccc1N1CCN([C@H]2CC[C@@H](NS(=O)(=O)c3cnc(Cl)c(Br)c3)CC2)CC1",
      "CC(C)Oc1ccccc1N1CCN([C@@H]2CC[C@H](NS(=O)(=O)c3cnc(Cl)c(Br)c3)CC2)CC1",
      "CC(C)Oc1ccccc1N1CCN([C@H]2CC[C@@H](NS(=O)(=O)c3cnc(Cl)c(Br)c3)CC2)CC1",
      "CC(C)Oc1ccccc1N1CCN([C@@H]2CC[C@H](NS(=O)(=O)c3cnc(Cl)c(Br)c3)CC2)CC1",
      "EOS"};
  unsigned int idx = 0;
  while (inSmiles[idx] != "EOS") {
    std::string smi1 = inSmiles[idx++];
    std::string smi2 = inSmiles[idx++];

    RWMol *m1 = SmilesToMol(smi1);
    REQUIRE(m1);
    RWMol *m2 = SmilesToMol(smi2);
    REQUIRE(m2);
    REQUIRE(m1->getNumAtoms() == m2->getNumAtoms());
    REQUIRE(m1->getNumBonds() == m2->getNumBonds());

    std::string csmiles1 = MolToSmiles(*m1, true);
    std::string csmiles2 = MolToSmiles(*m2, true);
    REQUIRE(csmiles1 == csmiles2);
    delete m1;
    delete m2;
  }
}

TEST_CASE(
    "Testing Github 45: stereochemistry information influencing non-stereo SMILES") {
  {
    RWMol *m;
    std::string smiles = "CC1CCC[13C]2(C)C1CC[14CH]2C(C)=O";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmiles1a = MolToSmiles(*m, false);
    std::string csmiles1b = MolToSmiles(*m, true);
    std::string smiles2 = "CC1CCC[C]2(C)C1CC[CH]2C(C)=O";
    delete m;
    m = SmilesToMol(smiles2);
    REQUIRE(m);
    std::string csmiles2a = MolToSmiles(*m, false);
    std::string csmiles2b = MolToSmiles(*m, true);

    REQUIRE(csmiles1a == csmiles2a);
    REQUIRE(csmiles1b != csmiles2b);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CC1CCC[C@@]2(C)C1CC[C@@H]2C(C)=O";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmiles1a = MolToSmiles(*m, false);
    std::string csmiles1b = MolToSmiles(*m, true);
    std::string smiles2 = "CC1CCC[C]2(C)C1CC[CH]2C(C)=O";
    delete m;
    m = SmilesToMol(smiles2);
    REQUIRE(m);
    std::string csmiles2a = MolToSmiles(*m, false);
    std::string csmiles2b = MolToSmiles(*m, true);

    REQUIRE(csmiles1a == csmiles2a);
    REQUIRE(csmiles1b != csmiles2b);
    delete m;
  }
}

TEST_CASE("Testing Github 206: Problems round-tripping P") {
  {
    RWMol *m;
    std::string smiles = "O=[PH3]";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmiles = MolToSmiles(*m, true);
    REQUIRE(csmiles == "O=[PH3]");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "O=P";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmiles = MolToSmiles(*m, true);
    REQUIRE(csmiles == "O=P");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "O=[PH]";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmiles = MolToSmiles(*m, true);
    REQUIRE(csmiles == "O=P");
    delete m;
  }
}

TEST_CASE(
    "Testing Github 210: flag possible stereocenters when calling assignStereochemistry()") {
  RWMol *m;
  std::string smiles = "O[C@H](F)CC(F)(Cl)I";
  m = SmilesToMol(smiles);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  REQUIRE(m->getAtomWithIdx(4)->hasProp(common_properties::_ChiralityPossible));
  delete m;
}

TEST_CASE("Testing Github 298: cannot generate smiles for ChEBI_50252") {
  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/test_data/ChEBI_50252.mol";
  RWMol *m = MolFileToMol(fName, false, false);

  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 80);
  REQUIRE(m->getNumBonds() == 210);
  m->updatePropertyCache(false);
  MolOps::fastFindRings(*m);

  std::string csmiles = MolToSmiles(*m);
  REQUIRE(csmiles != "");
  REQUIRE(csmiles.find("%100") == std::string::npos);

  delete m;
  m = SmilesToMol(csmiles, 0, false);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 80);
  REQUIRE(m->getNumBonds() == 210);
  delete m;
}

TEST_CASE(
    "Testing Github 378: SMILES parser doing the wrong thing for odd dot-disconnected construct") {
  {
    RWMol *m;
    std::string smiles = "C1.C1CO1.N1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(0, 1));
    REQUIRE(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::SINGLE);
    REQUIRE(m->getBondBetweenAtoms(3, 4));
    REQUIRE(m->getBondBetweenAtoms(3, 4)->getBondType() == Bond::SINGLE);
    REQUIRE(!m->getBondBetweenAtoms(1, 3));
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C1(O.C1)CO1.N1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(0, 2));
    REQUIRE(m->getBondBetweenAtoms(0, 2)->getBondType() == Bond::SINGLE);
    REQUIRE(m->getBondBetweenAtoms(0, 3));
    REQUIRE(m->getBondBetweenAtoms(0, 3)->getBondType() == Bond::SINGLE);
    REQUIRE(m->getBondBetweenAtoms(5, 4));
    REQUIRE(m->getBondBetweenAtoms(5, 4)->getBondType() == Bond::SINGLE);
    REQUIRE(!m->getBondBetweenAtoms(2, 3));
    delete m;
  }
}

TEST_CASE(
    "Testing Github 389: Add option to SmilesWriter to allow writing of all explicit hydrogens") {
  RWMol *m;
  std::string smiles = "CCO";
  m = SmilesToMol(smiles);
  REQUIRE(m);

  std::string csmiles = MolToSmiles(*m, true, false, -1, true, false, true);
  REQUIRE(csmiles != "");
  REQUIRE(csmiles.find("[CH3]") != std::string::npos);
  REQUIRE(csmiles.find("[CH2]") != std::string::npos);
  REQUIRE(csmiles.find("[OH]") != std::string::npos);

  delete m;
}

TEST_CASE("Testing handling of empty SMILES/SMARTS strings") {
  {
    RWMol *m;
    std::string smiles = "";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 0);

    std::string csmiles = MolToSmiles(*m);
    REQUIRE(csmiles == "");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "";
    m = SmartsToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 0);

    std::string csmiles = MolToSmarts(*m);
    REQUIRE(csmiles == "");
    delete m;
  }
}

TEST_CASE("testing smiles writing/canonicalization for modified molecules.") {
  std::string smiles = "c1ccccc1";
  ROMol *m = SmilesToMol(smiles);
  REQUIRE(m);

  m->getAtomWithIdx(0)->setAtomicNum(8);
  std::string smi = MolToSmiles(*m, true);
  REQUIRE(smi == "c1ccocc1");
  delete m;
}

TEST_CASE(
    "testing github issue 532: _smilesAtomOutputOrder incorrect for dot disconnected molecules") {
  {
    std::string smiles = "O.CO";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);

    std::string smi = MolToSmiles(*m, true);
    REQUIRE(smi == "CO.O");

    std::vector<unsigned int> atmOrder;
    REQUIRE(m->hasProp(common_properties::_smilesAtomOutputOrder));
    m->getProp(common_properties::_smilesAtomOutputOrder, atmOrder);
    REQUIRE(atmOrder.size() == 3);
    REQUIRE(atmOrder[0] == 1);
    REQUIRE(atmOrder[1] == 2);
    REQUIRE(atmOrder[2] == 0);

    delete m;
  }
  {
    std::string smiles = "CO.O";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);

    std::string smi = MolToSmiles(*m, true);
    REQUIRE(smi == "CO.O");

    std::vector<unsigned int> atmOrder;
    REQUIRE(m->hasProp(common_properties::_smilesAtomOutputOrder));
    m->getProp(common_properties::_smilesAtomOutputOrder, atmOrder);
    REQUIRE(atmOrder.size() == 3);
    REQUIRE(atmOrder[0] == 0);
    REQUIRE(atmOrder[1] == 1);
    REQUIRE(atmOrder[2] == 2);

    delete m;
  }
}

TEST_CASE(
    "testing github issue 760: reversed stereochemistry with sulfoxides and ring closures") {
  std::string smiles = "C[S@](Cl)=O";
  ROMol *m = SmilesToMol(smiles);
  REQUIRE(m);
  std::string csmi = MolToSmiles(*m, true);
  delete m;

  smiles = "C[S@]2=O.Cl2";
  m = SmilesToMol(smiles);
  REQUIRE(m);
  std::string csmi2 = MolToSmiles(*m, true);
  REQUIRE(csmi == csmi2);
  delete m;
}

TEST_CASE(
    "testing github issue 786: chiral order for ring closure after branch") {
  {
    std::string smiles = "C1CN[C@H]1O";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;

    smiles = "C1CN[C@@H](O)1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmi2 = MolToSmiles(*m, true);
    REQUIRE(csmi == csmi2);
    delete m;
  }
  {
    std::string smiles = "C1CN[C@]1(O)N";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;

    smiles = "C1CN[C@](O)(N)1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmi2 = MolToSmiles(*m, true);
    REQUIRE(csmi == csmi2);
    delete m;
  }
  {
    std::string smiles = "C1CN[C@]12(O).N2";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;

    smiles = "C1CN[C@](O)12.N2";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmi2 = MolToSmiles(*m, true);
    REQUIRE(csmi == csmi2);
    delete m;

    smiles = "C1CN[C@@]1(O)2.N2";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    csmi2 = MolToSmiles(*m, true);
    REQUIRE(csmi == csmi2);
    delete m;

    smiles = "C1CN[C@]2(O)1.N2";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    csmi2 = MolToSmiles(*m, true);
    REQUIRE(csmi == csmi2);
    delete m;
  }
  {
    std::string smiles = "C[C@]1(O)NCC1";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;

    smiles = "C[C@@](O)1NCC1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmi2 = MolToSmiles(*m, true);
    REQUIRE(csmi == csmi2);
    delete m;
  }
  {
    std::string smiles = "C[C@]1(NCC1)O";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;
    // so many pathologically ugly SMILES:
    smiles = "C[C@](NCC1)(O)1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmi2 = MolToSmiles(*m, true);
    REQUIRE(csmi == csmi2);
    delete m;
  }

  {  // Andrew's original real example:
    std::string smiles =
        "CC(C)[C@]1(N)CC[C@]2([C@@H](O2)CCC(=C)[C@H](CC[C@@](/C=C1)(C)O)O)C";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;
    smiles =
        "CC(C)[C@@](N)1CC[C@]2([C@@H](O2)CCC(=C)[C@H](CC[C@@](/C=C1)(C)O)O)C";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmi2 = MolToSmiles(*m, true);
    REQUIRE(csmi == csmi2);
    delete m;
  }
}

TEST_CASE(
    "testing github issue 1652: chiral order for ring closure after branch for the first atom in the SMILES string") {
  {
    std::string smiles = "Cl[C@](F)1CC[C@H](F)CC1";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;

    smiles = "[C@](Cl)(F)1CC[C@H](F)CC1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmi2 = MolToSmiles(*m, true);
    REQUIRE(csmi == csmi2);
    delete m;
  }
  {
    std::string smiles = "F[C@@]1(C)CCO1";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;

    smiles = "[C@@](F)1(C)CCO1";
    m = SmilesToMol(smiles);
    REQUIRE(m);
    std::string csmi2 = MolToSmiles(*m, true);
    REQUIRE(csmi == csmi2);
    delete m;
  }
}

TEST_CASE("testing dative bond support") {
  {
    std::string smiles = "CCC(=O)O->[Cu]";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);

    int dative_bond_count = 0;
    for (size_t i = 0; i < m->getNumBonds(); i++) {
      if (m->getBondWithIdx(i)->getBondType() == Bond::DATIVE) {
        dative_bond_count++;
      }
    }
    REQUIRE(dative_bond_count == 1);

    std::string out_smiles = MolToSmiles(*m, true);
    delete m;
    REQUIRE(out_smiles == "CCC(=O)[OH]->[Cu]");
  }
  {
    std::string smiles = "CCC(=O)O->[Cu]<-OC(O)CC";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);

    int dative_bond_count = 0;
    for (size_t i = 0; i < m->getNumBonds(); i++) {
      if (m->getBondWithIdx(i)->getBondType() == Bond::DATIVE) {
        dative_bond_count++;
      }
    }
    REQUIRE(dative_bond_count == 2);

    std::string out_smiles = MolToSmiles(*m, true);
    delete m;
    REQUIRE(out_smiles == "CCC(=O)[OH]->[Cu]<-[OH]C(O)CC");
  }
}

TEST_CASE(
    "GitHub Issue 1219: Stereochemistry not output to SMILES when allHsExplicit=True") {
  {
    std::string smiles = "C[C@H](F)Cl";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);

    bool doIsomericSmiles = true;
    bool doKekule = false;
    int rootedAtAtom = -1;
    bool canonical = true, allBondsExplicit = false, allHsExplicit = true;
    std::string csmi = MolToSmiles(*m, doIsomericSmiles, doKekule, rootedAtAtom,
                                   canonical, allBondsExplicit, allHsExplicit);
    REQUIRE(csmi == "[CH3][C@H]([F])[Cl]");
    delete m;
  }
  {  // another manifestation was that chiral flags were not output for atoms
     // not in the organic subset
    std::string smiles = "C[Si@H](F)Cl";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    bool doIsomericSmiles = true;
    std::string csmi = MolToSmiles(*m, doIsomericSmiles);
    REQUIRE(csmi == "C[Si@H](F)Cl");
    delete m;
  }
}

TEST_CASE("Testing the SmilesParseParams class") {
  {
    std::string smiles = "C1=CC=CC=C1[H]";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 6);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
    delete m;

    {
      SmilesParserParams params;
      m = SmilesToMol(smiles, params);
      REQUIRE(m);
      REQUIRE(m->getNumAtoms() == 6);
      REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
      delete m;
    }
    {  // no removeHs, with sanitization
      SmilesParserParams params;
      params.removeHs = false;
      m = SmilesToMol(smiles, params);
      REQUIRE(m);
      REQUIRE(m->getNumAtoms() == 7);
      REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
      delete m;
    }
    {  // removeHs, no sanitization
      SmilesParserParams params;
      params.sanitize = false;
      m = SmilesToMol(smiles, params);
      REQUIRE(m);
      REQUIRE(m->getNumAtoms() == 6);
      REQUIRE(!m->getBondWithIdx(0)->getIsAromatic());
      delete m;
    }
    {  // no removeHs, no sanitization
      SmilesParserParams params;
      params.removeHs = false;
      params.sanitize = false;
      m = SmilesToMol(smiles, params);
      REQUIRE(m);
      REQUIRE(m->getNumAtoms() == 7);
      REQUIRE(!m->getBondWithIdx(0)->getIsAromatic());
      delete m;
    }
  }

  {  // basic name parsing
    std::string smiles = "CCCC the_name";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    delete m;
    {  // it's parsed:
      SmilesParserParams params;
      params.allowCXSMILES = false;
      m = SmilesToMol(smiles, params);
      REQUIRE(m);
      REQUIRE(m->hasProp(common_properties::_Name));
      REQUIRE(m->getProp<std::string>(common_properties::_Name) == "the_name");
      delete m;
    }
    {
      SmilesParserParams params;
      params.strictCXSMILES = false;
      params.parseName = false;
      m = SmilesToMol(smiles, params);
      REQUIRE(m);
      REQUIRE(m->getNumAtoms() == 4);
      REQUIRE(!m->hasProp(common_properties::_Name));
      delete m;
    }
  }
  {  // name parsing2
    std::string smiles = "CCCC\tthe_name";
    {  // no removeHs, no sanitization
      SmilesParserParams params;
      params.parseName = true;
      RWMol *m = SmilesToMol(smiles, params);
      REQUIRE(m);
      REQUIRE(m->getNumAtoms() == 4);
      REQUIRE(m->hasProp(common_properties::_Name));
      REQUIRE(m->getProp<std::string>(common_properties::_Name) == "the_name");
      delete m;
    }
  }
  {  // name parsing3
    std::string smiles = "CCCC\t  the_name  ";
    {  // no removeHs, no sanitization
      SmilesParserParams params;
      params.parseName = true;
      RWMol *m = SmilesToMol(smiles, params);
      REQUIRE(m);
      REQUIRE(m->getNumAtoms() == 4);
      REQUIRE(m->hasProp(common_properties::_Name));
      REQUIRE(m->getProp<std::string>(common_properties::_Name) == "the_name");
      delete m;
    }
  }
}

TEST_CASE("Testing the %(....) notation for SMILES ring closure numbers") {
  {
    const char *benzenes[6] = {
        "c1ccccc1",           "c%(1)ccccc%(1)",       "c%(12)ccccc%(12)",
        "c%(123)ccccc%(123)", "c%(1234)ccccc%(1234)", "c%(99999)ccccc%(99999)"};
    for (auto &i : benzenes) {
      BOOST_LOG(rdInfoLog) << "Test: " << i << " (should be read)" << std::endl;
      ROMol *m = SmilesToMol(i);
      REQUIRE(m);
      REQUIRE(m->getNumAtoms() == 6);
      REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
      std::string benzene = MolToSmiles(*m, false, false, -1, false);
      REQUIRE(benzene == "c1ccccc1");
      delete m;
    }

    const char *not_allowed[2] = {"c%()ccccc%()", "c%(100000)ccccc%(100000)"};
    for (auto &i : not_allowed) {
      BOOST_LOG(rdInfoLog) << "Test: " << i << " (should NOT be read)"
                           << std::endl;
      ROMol *m = SmilesToMol(i);
      REQUIRE(m == (ROMol *)nullptr);
      delete m;
    }
  }
}

TEST_CASE("Testing that isomeric SMILES is now the default output") {
  std::string smi = "C[C@H](Cl)Br";
  auto m = SmilesToMol(smi);
  REQUIRE(m);
  auto csmi = MolToSmiles(*m);
  REQUIRE(csmi.find("@") != std::string::npos);
  delete m;
}

TEST_CASE("Testing constructs like [#6]") {
  std::string smi = "[#6][12#6]";
  auto m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(0)->getAtomicNum() == 6);
  REQUIRE(m->getAtomWithIdx(0)->getIsotope() == 0);
  REQUIRE(m->getAtomWithIdx(1)->getAtomicNum() == 6);
  REQUIRE(m->getAtomWithIdx(1)->getIsotope() == 12);
  delete m;
}

TEST_CASE(
    "Testing Github #1925: Atom with bond to itself is accepted by the SMILES parser.") {
  std::string smi = "C1CC111";
  RWMol *m = nullptr;
  m = SmilesToMol(smi);
  REQUIRE(!m);
}

TEST_CASE("Testing the random Generator for SMILES") {
  // it's not trivial to test this because we're using std::rand(), which does
  // not give consistent results across platforms. It's not worth adding the
  // complexity of a real RNG, so we do some hand waving in the tests
  std::srand(0xf00d);  // be sure we use it for testcase!
  const std::vector<std::string> benzenes = {"COc1ccnc(CC)c1C"};

  const std::vector<std::string> rulesmiles = {
      "COc1ccnc(CC)c1C",     "O(C)c1ccnc(CC)c1C",   "c1(OC)ccnc(CC)c1C",
      "c1c(OC)c(C)c(CC)nc1", "c1cc(OC)c(C)c(CC)n1", "n1ccc(OC)c(C)c1CC",
      "c1(CC)nccc(OC)c1C",   "C(c1nccc(OC)c1C)C",   "CCc1nccc(OC)c1C",
      "c1(C)c(OC)ccnc1CC",   "Cc1c(OC)ccnc1CC",
  };

  for (auto bz : benzenes) {
    ROMol *m = SmilesToMol(bz);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 11);
    for (unsigned int j = 0; j < m->getNumAtoms(); ++j) {
      auto rulebenzene =
          MolToSmiles(*m, true, false, j, false, false, false, false);
      REQUIRE(rulebenzene == rulesmiles[j]);
      std::set<std::string> rsmis;
      for (unsigned int iter = 0; iter < 10; ++iter) {
        auto randombenzene =
            MolToSmiles(*m, true, false, j, false, false, false, true);
        rsmis.insert(randombenzene);
      }
      // we will get dupes, but there's enough choice available here that we
      // should have gotten at least 3 unique
      REQUIRE(rsmis.size() >= 3);
    }

    // confirm that we also use random starting points:
    std::set<char> starts;
    for (unsigned int iter = 0; iter < 50; ++iter) {
      auto randombenzene =
          MolToSmiles(*m, true, false, -1, false, false, false, true);
      starts.insert(randombenzene[0]);
    }
    // we will get dupes, but there's enough choice available here that we
    // should have gotten at least 3 unique
    REQUIRE(starts.find('C') != starts.end());
    REQUIRE(starts.find('c') != starts.end());
    REQUIRE(
        (starts.find('n') != starts.end() || starts.find('O') != starts.end()));

    delete m;
  }
}

TEST_CASE(
    "Testing Github #1972: Incorrect tetrahedral stereo when reading SMILES with ring closure as last neighbor") {
  {
    std::vector<std::vector<std::string>> smiles = {
        {"[C@@]1(Cl)(F)(I).Br1", "[C@@](Br)(Cl)(F)(I)"},
        {"[C@@](Cl)(F)(I)1.Br1", "[C@@](Cl)(F)(I)Br"},
        {"[C@@](Cl)1(F)(I).Br1", "[C@@](Cl)(Br)(F)(I)"},
        {"[C@@](Cl)(F)1(I).Br1", "[C@@](Cl)(F)(Br)(I)"}};
    for (const auto &pr : smiles) {
      std::unique_ptr<ROMol> m1(SmilesToMol(pr[0]));
      std::unique_ptr<ROMol> m2(SmilesToMol(pr[1]));
      REQUIRE(m1);
      REQUIRE(m2);
      auto csmi1 = MolToSmiles(*m1);
      auto csmi2 = MolToSmiles(*m2);
      REQUIRE(csmi1 == csmi2);
    }
  }
  {  // even stupider examples
    std::vector<std::vector<std::string>> smiles = {
        {"[C@@]1(Cl)2(I).Br1.F2", "[C@@](Br)(Cl)(F)(I)"},
        {"[C@@](Cl)2(I)1.Br1.F2", "[C@@](Cl)(F)(I)Br"},
        {"[C@@]12(Cl)(I).Br1.F2", "[C@@](Br)(F)(Cl)(I)"},
        {"[C@@]21(Cl)(I).Br1.F2", "[C@@](F)(Br)(Cl)(I)"},
        {"[C@@](Cl)12(I).Br1.F2", "[C@@](Cl)(Br)(F)(I)"},
        {"[C@@](Cl)21(I).Br1.F2", "[C@@](Cl)(F)(Br)(I)"},
        {"[C@@](Cl)(I)21.Br1.F2", "[C@@](Cl)(I)(F)(Br)"},
        {"[C@@](Cl)(I)12.Br1.F2", "[C@@](Cl)(I)(Br)(F)"}};
    for (const auto &pr : smiles) {
      std::unique_ptr<ROMol> m1(SmilesToMol(pr[0]));
      std::unique_ptr<ROMol> m2(SmilesToMol(pr[1]));
      REQUIRE(m1);
      REQUIRE(m2);
      auto csmi1 = MolToSmiles(*m1);
      auto csmi2 = MolToSmiles(*m2);
      REQUIRE(csmi1 == csmi2);
    }
  }
}

TEST_CASE(
    "Testing Github #2556: Test correct parsing and fix memory leak for C1C1") {
  RWMol *m = nullptr;
  m = SmilesToMol("C1C1");
  REQUIRE(!m);
}

TEST_CASE(
    "Testing github issue #1028: Alternating canonical SMILES for ring with chiral N") {
  // note that due to the changes made for #3631, the N's originally used in
  // these tests are no longer considered to be chiral. I switched to using P
  // (and verified that P was also a problem before #1028 was fixed)

  std::string smi = "O[C@H]1CC2CCC(C1)[P@@]2C";
  const std::string ref = "C[P@]1C2CCC1C[C@H](O)C2";
  for (int i = 0; i < 3; ++i) {
    const auto mol = std::unique_ptr<ROMol>(SmilesToMol(smi));
    REQUIRE(mol);
    const std::string out = MolToSmiles(*mol);
    REQUIRE(out == ref);
    smi = out;
  }

  {
    std::string smi = "C[P@]1C[C@@H](O)C1";
    const std::string ref = smi;
    for (int i = 0; i < 3; ++i) {
      const auto mol = std::unique_ptr<ROMol>(SmilesToMol(smi));
      REQUIRE(mol);
      const std::string out = MolToSmiles(*mol);
      REQUIRE(out == ref);
      smi = out;
    }
  }
}

TEST_CASE("Testing github issue #3139: Partial bond mem leak") {
  const std::string smi = "COc(c1)cccc1C#";
  for (int i = 0; i < 3; ++i) {
    const auto mol = std::unique_ptr<ROMol>(SmilesToMol(smi));
    const auto sma = std::unique_ptr<ROMol>(SmartsToMol(smi));
  }
}

TEST_CASE("Failures/problems detected by OSS Fuzz") {
  // examples that should produce no molecule
  std::vector<std::string> failing_examples = {"C)"};
  for (auto smi : failing_examples) {
    const auto mol = std::unique_ptr<ROMol>(SmilesToMol(smi));
    // output which molecule is failing
    REQUIRE(!mol);
  }
}

TEST_CASE(
    "Testing Github Issue 3967: Double bond stereo gets flipped by SMILES reader/writer") {
  {
    auto mol = "C=c1s/c2n(c1=O)CCCCCCC\\N=2"_smiles;
    REQUIRE(mol);
    auto smi = MolToSmiles(*mol);
    CHECK(smi == "C=c1s/c2n(c1=O)CCCCCCC\\N=2");
  }
  {
    auto mol = R"SMI(C1=C\C/C=C2C3=C/C/C=C\C=C/C\3C\2\C=C/1)SMI"_smiles;
    REQUIRE(mol);
    auto smi = MolToSmiles(*mol);
    CHECK(smi == R"SMI(C1=C\C/C=C2\C3=C\C/C=C\C=C/C3C2\C=C/1)SMI");
  }
}

TEST_CASE(
    "Testing Github Issue 6349: Different SMARTS input formats lead to different SMILES outputs.") {
  auto checkSmartsToSmiles = [](const std::string &sma,
                                const std::string &refSmi) {
    std::unique_ptr<ROMol> molFromSmarts(SmartsToMol(sma));
    {
      std::string smi = MolToSmiles(*molFromSmarts);
      REQUIRE(smi == refSmi);
    }

    std::string molBlock = MolToMolBlock(*molFromSmarts);
    std::unique_ptr<ROMol> molFromBlock(
        MolBlockToMol(molBlock, /*sanitize =*/false, /*removeHs =*/false));
    {
      std::string smi = MolToSmiles(*molFromBlock);
      REQUIRE(smi == refSmi);
    }
  };
  checkSmartsToSmiles("[C]", "C");
  checkSmartsToSmiles("[C,N]", "*");
  checkSmartsToSmiles("[C,N]~[O,S]", "*~*");
  checkSmartsToSmiles("C-C(-[Cl,F,Br])-C", "*C(C)C");
}

TEST_CASE("test Parser Error Message") {
  const std::string smis[] = {
      "CC=(CO)C",
      "baz",
      "fff",
      "C+0",
      "[555555555555555555C]",
      "[Fe@TD]",
      "c%()ccccc%()",
      "c%(100000)ccccc%(100000)",
      "COc(c1)cccc1C#",
      "C)",
  };
  for (const auto &smi : smis) {
    // Test SMILES parsing
    {
      std::stringstream ss;
      rdErrorLog->SetTee(ss);

      auto mol = v2::SmilesParse::MolFromSmiles(smi);
      REQUIRE(!mol);

      rdErrorLog->ClearTee();
      auto error_msg = ss.str();
      REQUIRE(error_msg.find("check for mistakes around position") !=
              std::string::npos);
    }

    // Test SMARTS parsing
    {
      std::stringstream ss;
      rdErrorLog->SetTee(ss);

      auto mol = v2::SmilesParse::MolFromSmarts(smi);
      REQUIRE(!mol);

      rdErrorLog->ClearTee();
      auto error_msg = ss.str();
      REQUIRE(error_msg.find("check for mistakes around position") !=
              std::string::npos);
    }
  }
}
