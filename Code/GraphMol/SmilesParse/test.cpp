//
//  Copyright (C) 2003-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <iostream>
#include <string>
#include <GraphMol/RDKitBase.h>
#include "SmilesParse.h"
#include "SmilesWrite.h"
#include "SmartsWrite.h"
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/RDLog.h>
//#include <boost/log/functions.hpp>
using namespace RDKit;
using namespace std;
typedef ROMol Mol;

void testPass() {
  int i = 0;
  ROMol *mol, *mol2;
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing molecules which should parse." << std::endl;
  string smis[] = {
#if 1
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
#endif
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
    "[Db][Sg][Bh][Hs][Mt][Ds][Rg][Cn][Nh][Fl][Mc][Lv][Ts][Og]",  // new elements
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
    "EOS"
  };
  while (smis[i] != "EOS") {
    string smi = smis[i];
    BOOST_LOG(rdInfoLog) << "***: " << smi << std::endl;
    mol = SmilesToMol(smi);
    CHECK_INVARIANT(mol, smi);
    if (mol) {
      unsigned int nAts = mol->getNumAtoms();
      CHECK_INVARIANT(nAts != 0, smi.c_str());
      smi = MolToSmiles(*mol);
      // BOOST_LOG(rdInfoLog)<< "  > " << smi << std::endl;
      mol2 = SmilesToMol(smi);
      CHECK_INVARIANT(mol2->getNumAtoms() == nAts, smi.c_str())
      delete mol;
      delete mol2;
    }
    i++;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testFail() {
  int i = 0;
  Mol *mol;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing molecules which should fail to parse/sanitize." << std::endl;

  // alternate good and bad smiles here to ensure that the parser can resume
  // parsing
  // on good input:
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
      "EOS"};

  // turn off the error log temporarily:
  while (smis[i] != "EOS") {
    string smi = smis[i];
    boost::logging::disable_logs("rdApp.error");
    try {
      mol = SmilesToMol(smi);
    } catch (MolSanitizeException &) {
      mol = (Mol *)nullptr;
    }
    boost::logging::enable_logs("rdApp.error");
    if (!(i % 2)) {
      CHECK_INVARIANT(!mol, smi);
    } else {
      CHECK_INVARIANT(mol, smi);
      delete mol;
    }
    i++;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testDetails() {
  Mol *mol;
  Atom *a;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing details" << std::endl;
  // implicit/explicit H handling
  smi = "OC([OH])C[O-]";
  mol = SmilesToMol(smi);
  CHECK_INVARIANT(mol, smi);
  CHECK_INVARIANT(mol->getNumAtoms() == 5, "");
  a = mol->getAtomWithIdx(0);
  CHECK_INVARIANT(a->getImplicitValence() == 1, "");
  CHECK_INVARIANT(a->getExplicitValence() == 1, "");
  CHECK_INVARIANT(a->getNoImplicit() == 0, "");
  CHECK_INVARIANT(a->getFormalCharge() == 0, "");
  a = mol->getAtomWithIdx(2);
  CHECK_INVARIANT(a->getImplicitValence() == 0, "");
  CHECK_INVARIANT(a->getExplicitValence() == 2, "");
  CHECK_INVARIANT(a->getNoImplicit() == 1, "");
  CHECK_INVARIANT(a->getFormalCharge() == 0, "");
  a = mol->getAtomWithIdx(4);
  CHECK_INVARIANT(a->getImplicitValence() == 0, "");
  CHECK_INVARIANT(a->getExplicitValence() == 1, "");
  CHECK_INVARIANT(a->getNoImplicit() == 1, "");
  CHECK_INVARIANT(a->getFormalCharge() == -1, "");

  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testProblems() {
  Mol *mol;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing smiles that were previously problems"
                       << std::endl;

  // ring closure handling with branches/fragments
  VECT_INT_VECT rings;
  smi = "C1(CC1CC1CC1)";
  mol = SmilesToMol(smi);
  CHECK_INVARIANT(mol, smi);
  int ringCount = MolOps::findSSSR(*mol, rings);
  CHECK_INVARIANT(ringCount == 2, "");
  CHECK_INVARIANT(rings.size() == 2, "");
  CHECK_INVARIANT(rings[0].size() == 3, "");
  CHECK_INVARIANT(rings[1].size() == 3, "");

  // this is truly pathological, but both daylight
  //   and chemdraw parse it properly
  smi = "C1.C1CC1CC1";
  delete mol;
  mol = SmilesToMol(smi);
  CHECK_INVARIANT(mol, smi);
  ringCount = MolOps::findSSSR(*mol, rings);
  CHECK_INVARIANT(ringCount == 1, "");
  CHECK_INVARIANT(rings.size() == 1, "");
  CHECK_INVARIANT(rings[0].size() == 3, "");

  // here's another stupid case that we need to handle:
  delete mol;
  smi = "C1CC11CC1";
  mol = SmilesToMol(smi);
  CHECK_INVARIANT(mol, smi);
  ringCount = MolOps::findSSSR(*mol, rings);
  CHECK_INVARIANT(ringCount == 2, "");
  CHECK_INVARIANT(rings.size() == 2, "");
  CHECK_INVARIANT(rings[0].size() == 3, "");
  CHECK_INVARIANT(rings[1].size() == 3, "");

  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testBasicCanon() {
  Mol *mol;
  std::string smi, refSmi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing basic SMILES canonicalization" << std::endl;
#if 1
  smi = "C1OCCCC1";
  mol = SmilesToMol(smi);
  refSmi = MolToSmiles(*mol);

  delete mol;
  smi = "C1COCCC1";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol);
  TEST_ASSERT(refSmi == smi);

  delete mol;
  smi = "O1CCCCC1";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol);
  TEST_ASSERT(refSmi == smi);

  delete mol;
  smi = "OC=CC";
  mol = SmilesToMol(smi);
  refSmi = MolToSmiles(*mol);
  delete mol;
  smi = "CC=CO";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol);
  TEST_ASSERT(refSmi == smi);
  delete mol;
  smi = "C(C)=CO";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol);
  TEST_ASSERT(refSmi == smi);
  delete mol;
  smi = "C(O)=CC";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol);
  TEST_ASSERT(refSmi == smi);

  // --- These are related to Issue 109
  delete mol;
  smi = "C([H])Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getNumAtoms() == 2);
  refSmi = MolToSmiles(*mol);
  delete mol;
  smi = "CCl";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol);
  TEST_ASSERT(refSmi == smi);
  delete mol;
#endif
  // -- Issue 131
  smi = "P#[Ga]";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getNumAtoms() == 2);
  refSmi = MolToSmiles(*mol);
  delete mol;
  mol = SmilesToMol(refSmi);
  smi = MolToSmiles(*mol);
  TEST_ASSERT(refSmi == smi);

  delete mol;
  smi = "O=[Ba]";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getNumAtoms() == 2);
  refSmi = MolToSmiles(*mol);
  delete mol;
  mol = SmilesToMol(refSmi);
  smi = MolToSmiles(*mol);
  TEST_ASSERT(refSmi == smi);

  // make sure empty molecules return empty SMILES:
  delete mol;
  mol = new ROMol();
  smi = MolToSmiles(*mol);
  TEST_ASSERT(smi == "");

  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testLeak() {
  int i = 0;
  Mol *mol;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing a leak" << std::endl;

  smi = "C1CC1";
  for (i = 0; i < 1000000; i++) {
    mol = SmilesToMol(smi, 0, 1);
    if (!(i % 1000)) {
      BOOST_LOG(rdInfoLog) << i << std::endl;
    }

    delete mol;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testStereochem() {
  Mol *mol;
  std::string smi, refSmi, cip;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing handling of stereochemical smiles"
                       << std::endl;

  smi = "F[C@](Cl)(Br)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  refSmi = MolToSmiles(*mol, 1);

  delete mol;
  smi = "F[C@](Br)(I)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "F[C@](I)(Cl)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "Cl[C@](Br)(F)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "Cl[C@](F)(I)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "I[C@](F)(Br)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "I[C@](Br)(Cl)F";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "F[C@@](Br)(Cl)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "F[C@@](Cl)(I)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "Cl[C@@](Br)(I)F";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "Cl[C@@](F)(Br)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "[C@@](Cl)(F)(Br)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "F[C@H](Cl)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  refSmi = MolToSmiles(*mol, 1);

  delete mol;
  smi = "Br[C@H](F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  smi = MolToSmiles(*mol, 1);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "Br[C@]([H])(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "Br[C@](F)(Cl)[H]";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "Br[C@]1(F)(Cl).[H]1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "Br[C@H]1Cl.F1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "Br[C@]12Cl.F2.[H]1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "Br[C@]21Cl.F1.[H]2";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "[C@@H](Br)(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  delete mol;
  smi = "[H][C@@](Br)(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(smi == refSmi);

  // an additional set of test cases from the Chirality notes document.
  // one can never have too many tests of this stuff.
  delete mol;
  smi = "F[C@]([H])(O)C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  delete mol;
  smi = "F[C@]1([H])OC1";
  mol = SmilesToMol(smi);
  // TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete mol;
  smi = "F[C@H](O)C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  delete mol;
  smi = "F[C@@H]1OC1";
  mol = SmilesToMol(smi);
  // TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete mol;
  smi = "[C@](F)([H])(O)C";
  mol = SmilesToMol(smi);
  // TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  delete mol;
  smi = "[C@@]1(F)([H])OC1";
  mol = SmilesToMol(smi);
  // TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete mol;
  smi = "[C@@H](F)(O)C";
  mol = SmilesToMol(smi);
  // TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  smi = MolToSmiles(*mol, true);
  TEST_ASSERT(smi == "C[C@@H](O)F")
  smi = MolToSmiles(*mol, true, false, 0);
  TEST_ASSERT(smi == "[C@H](C)(O)F")

  delete mol;
  smi = "[C@@H]1(F)OC1";
  mol = SmilesToMol(smi);
  // TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  smi = MolToSmiles(*mol, true);
  TEST_ASSERT(smi == "F[C@H]1CO1")
  smi = MolToSmiles(*mol, true, false, 0);
  TEST_ASSERT(smi == "[C@H]1(F)CO1")

  delete mol;
  smi = "C1O[C@H]1F";
  mol = SmilesToMol(smi);
  // TEST_ASSERT(mol->getAtomWithIdx(2)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  delete mol;
  smi = "C1O[C@@]1([H])F";
  mol = SmilesToMol(smi);
  // TEST_ASSERT(mol->getAtomWithIdx(2)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

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
  TEST_ASSERT(refSmi == smi);

  delete mol;
  smi = "Br\\C=C/F";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);

  delete mol;
  smi = "Br/C=C\\F";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);

  delete mol;
  smi = "F/C=C\\Br";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);

  //-- trans --
  delete mol;
  smi = "F\\C=C\\Br";
  mol = SmilesToMol(smi);
  refSmi = MolToSmiles(*mol, 1);

  delete mol;
  mol = SmilesToMol(refSmi);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);

  delete mol;
  smi = "Br\\C=C\\F";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);

  delete mol;
  smi = "Br/C=C/F";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);

  delete mol;
  smi = "F/C=C/Br";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);

  //-- more complex --
  delete mol;
  smi = "F\\C=C(/Cl)\\Br";
  mol = SmilesToMol(smi);
  refSmi = MolToSmiles(*mol, 1);

  delete mol;
  mol = SmilesToMol(refSmi);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);

  delete mol;
  smi = "F/C=C(\\Cl)/Br";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);

  delete mol;
  smi = "F/C=C(\\Cl)Br";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);

  delete mol;
  smi = "F/C=C(Cl)/Br";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);

  //-- combine chirality with cis/trans --
  delete mol;
  smi = "F[C@H](Cl)\\C=C(/F)";
  mol = SmilesToMol(smi);
  refSmi = MolToSmiles(*mol, 1);

  delete mol;
  mol = SmilesToMol(refSmi);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);

  delete mol;
  smi = "F[C@H](Cl)/C=C(\\F)";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);
  delete mol;

  smi = "Cl[C@@H](F)/C=C(\\F)";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);
  delete mol;

  smi = "Cl[C@@H](F)\\C=C(/F)";
  mol = SmilesToMol(smi);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);
  delete mol;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue127() {
  Mol *mol, *mol2;
  std::string smi, refSmi, tempStr;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 127 (chiral smiles with fused rings)"
                       << std::endl;

  smi = "Cl[C@]12[Si]C(C2)O1";
  mol = SmilesToMol(smi);
  // mol->debugMol(std::cout);
  TEST_ASSERT(mol);

#if 1
  // first roundtrip the non-chiral SMILES:
  refSmi = MolToSmiles(*mol);
  mol2 = SmilesToMol(refSmi);
  TEST_ASSERT(mol2);
  tempStr = MolToSmiles(*mol2);
  TEST_ASSERT(refSmi == tempStr);
  delete mol2;
#endif

  // now do the true SMILES:
  refSmi = MolToSmiles(*mol, 1);
  mol2 = SmilesToMol(refSmi);
  // mol2->debugMol(std::cout);
  TEST_ASSERT(mol2);
  tempStr = MolToSmiles(*mol2, 1);
  // std::cout << refSmi << " : " << tempStr << std::endl;
  TEST_ASSERT(refSmi == tempStr);
  delete mol2;
  delete mol;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue143() {
  Mol *mol;
  std::string smi, refSmi, tempStr;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing Issue 143 (removing chiral tags for non-chiral centers)"
      << std::endl;

  smi = "C[C@](C)(C)C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  refSmi = MolToSmiles(*mol, true);
  TEST_ASSERT(refSmi == "CC(C)(C)C");
  delete mol;

  smi = "CC[C@](C)(C)C=O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  refSmi = MolToSmiles(*mol, true);
  TEST_ASSERT(refSmi == "CCC(C)(C)C=O");
  delete mol;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue151() {
  Mol *mol, *mol2;
  std::string smi, refSmi, tempStr;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 151 (Chiral centers in rings with "
                          "hydrogen on them not handled correctly)"
                       << std::endl;

  smi = "C1S[C@H]1O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(2)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(2)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);

  refSmi = MolToSmiles(*mol, true);
  TEST_ASSERT(refSmi == "O[C@H]1CS1");
  mol2 = SmilesToMol(refSmi);
  TEST_ASSERT(mol2);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(refSmi == smi);
  delete mol;
  delete mol2;

  smi = "F[C@@H]1O[C@H](Cl)S1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(2)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(4)->getChiralTag() == Atom::CHI_UNSPECIFIED);

  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  TEST_ASSERT(mol->getAtomWithIdx(3)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(3)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);

  refSmi = MolToSmiles(*mol, true);
  TEST_ASSERT(refSmi == "F[C@@H]1O[C@H](Cl)S1");
  mol2 = SmilesToMol(refSmi);
  TEST_ASSERT(mol2);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(refSmi == smi);
  delete mol;
  delete mol2;

  smi = "Cl[C@@H]1S[C@@H](O1)F";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(2)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(4)->getChiralTag() == Atom::CHI_UNSPECIFIED);

  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  TEST_ASSERT(mol->getAtomWithIdx(3)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(3)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);

  refSmi = MolToSmiles(*mol, true);
  TEST_ASSERT(refSmi == "F[C@@H]1O[C@H](Cl)S1");
  mol2 = SmilesToMol(refSmi);
  TEST_ASSERT(mol2);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(refSmi == smi);
  delete mol;
  delete mol2;

  smi = "Cl[C@@H]1O[C@H](F)S1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(2)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(4)->getChiralTag() == Atom::CHI_UNSPECIFIED);

  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  TEST_ASSERT(mol->getAtomWithIdx(3)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(3)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);

  refSmi = MolToSmiles(*mol, true);
  TEST_ASSERT(refSmi == "F[C@H]1O[C@@H](Cl)S1");
  mol2 = SmilesToMol(refSmi);
  TEST_ASSERT(mol2);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(refSmi == smi);
  delete mol;
  delete mol2;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue153() {
  std::string code;
  Mol *mol, *mol2;
  std::string smi, refSmi, tempStr;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing Issue 153 (Incorrect order of ring-closure bonds from SMILES)"
      << std::endl;

  smi = "C1(O[C@H]12)S2";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(2)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(2)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, code);
  TEST_ASSERT(code == "S");

  refSmi = MolToSmiles(*mol, true);
  TEST_ASSERT(refSmi == "O1C2S[C@H]12");
  mol2 = SmilesToMol(refSmi);
  TEST_ASSERT(mol2);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(refSmi == smi);
  delete mol;
  delete mol2;

  smi = "C1(O[C@H]21)S2";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(2)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(2)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, code);
  TEST_ASSERT(code == "R");

  refSmi = MolToSmiles(*mol, true);
  TEST_ASSERT(refSmi == "O1C2S[C@@H]12");
  mol2 = SmilesToMol(refSmi);
  TEST_ASSERT(mol2);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(refSmi == smi);
  delete mol;
  delete mol2;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue157() {
  std::string code;
  Mol *mol, *mol2;
  std::string smi, refSmi, tempStr;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 157 (Symmetric molecules with "
                          "multiple chiral centers badly canonicalized)"
                       << std::endl;

#if 1
  smi = "O[C@](C)(Cl)[C@@](O)(Cl)C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  TEST_ASSERT(mol->getAtomWithIdx(4)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);

  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, code);
  TEST_ASSERT(code == "R");
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, code);
  TEST_ASSERT(code == "S");

  refSmi = MolToSmiles(*mol, true);
  mol2 = SmilesToMol(refSmi);
  TEST_ASSERT(mol2);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(refSmi == smi);
  delete mol;
  delete mol2;

  smi = "Cl[C@@](C)1CC[C@@](C)(C1)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  refSmi = MolToSmiles(*mol, true);
  mol2 = SmilesToMol(refSmi);
  TEST_ASSERT(mol2);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(refSmi == smi);
  delete mol;
  delete mol2;

  BOOST_LOG(rdInfoLog) << "-**-**---------------------------------------"
                       << std::endl;
  smi = "[H][C@@]12CC(CO1)CN2";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, smi);
  TEST_ASSERT(smi == "S");
  refSmi = MolToSmiles(*mol, true);
  BOOST_LOG(rdInfoLog) << refSmi << std::endl;
  mol2 = SmilesToMol(refSmi);
  TEST_ASSERT(mol2);
  smi = MolToSmiles(*mol2, true);
  BOOST_LOG(rdInfoLog) << refSmi << std::endl;
  BOOST_LOG(rdInfoLog) << smi << std::endl;
  TEST_ASSERT(refSmi == smi);
  delete mol;
  delete mol2;
#endif
  smi = "[H][C@@]12C[14C@@](C=C1)(C3C2C(NC3=O)=O)[H]";
  // smi="C1=C[C@@H]2C[C@H]1C1C(=O)NC(=O)C21";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, smi);
  TEST_ASSERT(smi == "R");
  mol->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, smi);
  TEST_ASSERT(smi == "S");
  // mol->debugMol(std::cout);
  refSmi = MolToSmiles(*mol, true);
  mol2 = SmilesToMol(refSmi);
  TEST_ASSERT(mol2);
  smi = MolToSmiles(*mol2, true);
  BOOST_LOG(rdInfoLog) << refSmi << std::endl;
  BOOST_LOG(rdInfoLog) << smi << std::endl;
  TEST_ASSERT(refSmi == smi);
  delete mol;
  delete mol2;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue159() {
  Mol *mol;
  std::string smi, refSmi, tempStr;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing Issue 159 (cis/trans wrong in some branched systems)"
      << std::endl;

  smi = "C/C=C/O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);

  TEST_ASSERT(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  refSmi = MolToSmiles(*mol, 1);

  delete mol;
  smi = "C(\\C)=C/O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);

  TEST_ASSERT(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);

  delete mol;
  smi = "C(\\\\C)=C/O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);

  TEST_ASSERT(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);

  delete mol;
  smi = "C(=C/O)\\C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(0)->getStereo() == Bond::STEREOE);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);

  delete mol;
  smi = "C(\\C/C=C/Cl)=C/O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(4)->getStereo() == Bond::STEREOE);
  TEST_ASSERT(mol->getBondWithIdx(2)->getStereo() == Bond::STEREOE);

  delete mol;
  smi = "O=C\\C=C/F";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(0)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
  TEST_ASSERT(mol->getBondWithIdx(2)->getStereo() == Bond::STEREOZ);

  delete mol;
  smi = "C(/C=O)=C/F";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);
  TEST_ASSERT(mol->getBondWithIdx(2)->getStereo() == Bond::STEREOZ);

  delete mol;
  smi = "C(=C/F)/C=O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(0)->getStereo() == Bond::STEREOZ);
  TEST_ASSERT(mol->getBondWithIdx(3)->getStereo() == Bond::STEREONONE);

  delete mol;
  smi = "C(=O)\\C=C/Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(2)->getStereo() == Bond::STEREOZ);
  TEST_ASSERT(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);

  delete mol;
  smi = "CC(=O)\\C=C/Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(3)->getStereo() == Bond::STEREOZ);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

  delete mol;
  smi = "C(=O)\\N=C\\Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(2)->getStereo() == Bond::STEREOE);
  TEST_ASSERT(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);

  delete mol;
  smi = "CC(=O)\\N=C\\Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(3)->getStereo() == Bond::STEREOE);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

  delete mol;
  smi = "C(/Br)(=C/Cl)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);

  delete mol;
  smi = "C(=C/Cl)(/Br)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(0)->getStereo() == Bond::STEREOZ);

  delete mol;
  smi = "Cl\\C=C(\\Br)";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);

  delete mol;
  smi = "Cl\\C(=C\\Br)";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);

  delete mol;
  smi = "C(/C=C/C)";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  delete mol;
  smi = "C(/C)=C/C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);

  // ---------
  // These next few molecules test propagation of bond flips:
  // ---------
  delete mol;
  smi = "Cl/C=C(/C=C/C)\\C=C\\Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
  TEST_ASSERT(mol->getBondWithIdx(3)->getStereo() == Bond::STEREOE);
  TEST_ASSERT(mol->getBondWithIdx(6)->getStereo() == Bond::STEREOE);

  delete mol;
  smi = "C(/C=C/C)(\\C=C\\Br)=C\\Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  TEST_ASSERT(mol->getBondWithIdx(4)->getStereo() == Bond::STEREOE);
  TEST_ASSERT(mol->getBondWithIdx(6)->getStereo() == Bond::STEREOZ);

  delete mol;
  smi = "Br/C=C/C(/C=C/C)=C\\Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  TEST_ASSERT(mol->getBondWithIdx(4)->getStereo() == Bond::STEREOE);
  TEST_ASSERT(mol->getBondWithIdx(6)->getStereo() == Bond::STEREOZ);

  delete mol;
  smi = "Cl/C=C(/C=C/C=C\\F)\\C=C\\Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
  TEST_ASSERT(mol->getBondWithIdx(3)->getStereo() == Bond::STEREOE);
  TEST_ASSERT(mol->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
  TEST_ASSERT(mol->getBondWithIdx(8)->getStereo() == Bond::STEREOE);
  delete mol;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue175() {
  Mol *mol;
  std::string smi, refSmi, tempStr;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 175 (cis/trans wrong on ring closures)"
                       << std::endl;

  smi = "Cl\\C=C1.F/1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  delete mol;

  smi = "Cl\\C=C1CN/1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  delete mol;

  smi = "C/1=C/F.F1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(0)->getStereo() == Bond::STEREOZ);
  delete mol;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue176() {
  Mol *mol;
  std::string smi, refSmi, tempStr;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing Issue 176 (problems with 'mol BOND ring_number')"
      << std::endl;

  smi = "C1CC1C1CC1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumBonds() == 7);

  delete mol;
  smi = "C1CC1C1CC-1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumBonds() == 7);

  delete mol;
  smi = "C1CC1C1CC=1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumBonds() == 7);

  delete mol;
  smi = "C1CC1C=1CC1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumBonds() == 7);

  delete mol;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue180() {
  Mol *mol;
  std::string smi, refSmi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 180: Z/E problems" << std::endl;

  smi = "Cl/C(=N\\O)/C(=N\\O)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
  TEST_ASSERT(mol->getBondWithIdx(4)->getStereo() == Bond::STEREOE);
  refSmi = MolToSmiles(*mol, 1);

  delete mol;
  smi = "Cl/C(/C(Br)=N\\O)=N\\O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(3)->getStereo() == Bond::STEREOE);
  TEST_ASSERT(mol->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);
  delete mol;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue184() {
  Mol *mol;
  std::string smi, refSmi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing Issue 184: Cis/Trans incorrect on ring-closure bonds"
      << std::endl;

  smi = "C1NC(Cl)C(=N\\O)/C1=N\\O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  // mol->debugMol(std::cout);
  TEST_ASSERT(mol->getBondWithIdx(4)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(mol->getBondWithIdx(4)->getStereo() == Bond::STEREOZ);
  TEST_ASSERT(mol->getBondWithIdx(7)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(mol->getBondWithIdx(7)->getStereo() == Bond::STEREOZ);
  refSmi = MolToSmiles(*mol, 1);
  delete mol;
  mol = SmilesToMol(refSmi);
  TEST_ASSERT(mol);

  for (RWMol::BondIterator bondIt = mol->beginBonds();
       bondIt != mol->endBonds(); bondIt++) {
    if ((*bondIt)->getBondType() == Bond::DOUBLE) {
      TEST_ASSERT((*bondIt)->getStereo() == Bond::STEREOZ);
    }
  }

  smi = MolToSmiles(*mol, 1);
  TEST_ASSERT(refSmi == smi);
  delete mol;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue185() {
  Mol *mol;
  std::string smi, refSmi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing Issue 185: Cis/Trans incorrect on writing branches"
      << std::endl;

  // start with a simple E/Z handling case with branches:
  smi = "C(/C)=N/O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
  refSmi = MolToSmiles(*mol, 1, 0, 0);
  BOOST_LOG(rdInfoLog) << refSmi << std::endl;
  TEST_ASSERT(refSmi == "C(\\C)=N\\O");
  delete mol;
  // make sure we can round-trip:
  mol = SmilesToMol(refSmi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
  delete mol;

  // now make it more complex
  smi = "CC(=N\\O)/C=P/N";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  TEST_ASSERT(mol->getBondWithIdx(4)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(mol->getBondWithIdx(4)->getStereo() == Bond::STEREOE);
  refSmi = MolToSmiles(*mol, 1);
  BOOST_LOG(rdInfoLog) << refSmi << std::endl;
  delete mol;
  mol = SmilesToMol(refSmi);
  TEST_ASSERT(mol);

  for (RWMol::BondIterator bondIt = mol->beginBonds();
       bondIt != mol->endBonds(); bondIt++) {
    if ((*bondIt)->getBondType() == Bond::DOUBLE) {
      TEST_ASSERT((*bondIt)->getStereo() == Bond::STEREOE);
    }
  }
  smi = MolToSmiles(*mol, 1);
  // std::cout << "ref: " << refSmi << " -> " << smi << std::endl;
  TEST_ASSERT(refSmi == smi);

  // now repeat that experiment, but this time root the SMILES so that
  // we go in a "sensible" order:
  delete mol;
  smi = "CC(=N\\O)/C=P/N";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  refSmi = MolToSmiles(*mol, true, false, 6);
  BOOST_LOG(rdInfoLog) << refSmi << std::endl;
  TEST_ASSERT(refSmi == "N/P=C/C(C)=N/O");
  delete mol;
  mol = SmilesToMol(refSmi);
  TEST_ASSERT(mol);
  for (RWMol::BondIterator bondIt = mol->beginBonds();
       bondIt != mol->endBonds(); bondIt++) {
    if ((*bondIt)->getBondType() == Bond::DOUBLE) {
      TEST_ASSERT((*bondIt)->getStereo() == Bond::STEREOE);
    }
  }
  delete mol;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue191() {
  Mol *mol;
  std::string smi, refSmi;
  int numE = 0;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 191: Bad bond directions in a branch"
                       << std::endl;

  smi = "C2=NNC(N=C2)=N\\N=C\\c1ccccc1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(7)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(mol->getBondWithIdx(7)->getStereo() == Bond::STEREOE);
  refSmi = MolToSmiles(*mol, 1);
  delete mol;
  // std::cout << "ref: " << refSmi << std::endl;
  mol = SmilesToMol(refSmi);
  TEST_ASSERT(mol);
  // mol->debugMol(std::cout);
  numE = 0;
  for (RWMol::BondIterator bondIt = mol->beginBonds();
       bondIt != mol->endBonds(); bondIt++) {
    if ((*bondIt)->getBondType() == Bond::DOUBLE) {
      TEST_ASSERT((*bondIt)->getStereo() != Bond::STEREOZ);
      if ((*bondIt)->getStereo() == Bond::STEREOE) {
        numE++;
      }
    }
  }
  TEST_ASSERT(numE == 1);
  smi = MolToSmiles(*mol, 1);
  // std::cout << "ref: " << refSmi << " -> " << smi << std::endl;
  TEST_ASSERT(refSmi == smi);
  delete mol;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue256() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 256: SMILES yields incorrect structure"
                       << std::endl;

  v2::SmilesParse::SmilesParserParams ps;
  ps.sanitize = false;
  {
    auto smi = "C1CC[C+]1=1CCC1";
    auto mol = v2::SmilesParse::MolFromSmiles(smi, ps);
    TEST_ASSERT(mol);
    auto bond = mol->getBondBetweenAtoms(3, 0);
    TEST_ASSERT(bond)
    TEST_ASSERT(bond->getBondType() == Bond::SINGLE);
    bond = mol->getBondBetweenAtoms(3, 6);
    TEST_ASSERT(bond)
    TEST_ASSERT(bond->getBondType() == Bond::DOUBLE);
  }

  {
    auto smi = "C1CC[C+]=11CCC1";
    auto mol = v2::SmilesParse::MolFromSmiles(smi, ps);
    TEST_ASSERT(mol);
    auto bond = mol->getBondBetweenAtoms(3, 0);
    TEST_ASSERT(bond)
    TEST_ASSERT(bond->getBondType() == Bond::DOUBLE);
    bond = mol->getBondBetweenAtoms(3, 6);
    TEST_ASSERT(bond)
    TEST_ASSERT(bond->getBondType() == Bond::SINGLE);
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue266() {
  RWMol *mol;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 266: kekulized SMILES output"
                       << std::endl;

  smi = "c1ccccc1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol);
  TEST_ASSERT(smi == "c1ccccc1");

  MolOps::Kekulize(*mol);
  smi = MolToSmiles(*mol);
  TEST_ASSERT(smi == "C1=CC=CC=C1");
  delete mol;

  smi = "c1ccccc1c1ccccc1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol);
  TEST_ASSERT(smi == "c1ccc(-c2ccccc2)cc1");

  MolOps::Kekulize(*mol);
  smi = MolToSmiles(*mol);
  TEST_ASSERT(smi == "C1=CC=C(C2=CC=CC=C2)C=C1");
  delete mol;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testRootedAt() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing rootedAtAtom functionality" << std::endl;

  {
    RWMol *mol;
    std::string smi;
    smi = "CN(C)C";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    smi = MolToSmiles(*mol, false, false, -1);
    TEST_ASSERT(smi == "CN(C)C");
    smi = MolToSmiles(*mol, false, false, 1);
    TEST_ASSERT(smi == "N(C)(C)C");
    smi = MolToSmiles(*mol, false, false, 2);
    TEST_ASSERT(smi == "CN(C)C");
    delete mol;
  }
  {
    // This was github issue #182:
    RWMol mol;
    std::string smi;
    smi = MolToSmiles(mol);
    TEST_ASSERT(smi == "");
    smi = MolToSmiles(mol, false, false, 0);
    TEST_ASSERT(smi == "");
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIsotopes() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing isotope handling" << std::endl;

  {
    std::string smi = "C[13C](C)(C)C";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(feq(mol->getAtomWithIdx(1)->getMass(), 13.0034));
    smi = MolToSmiles(*mol, false);
    TEST_ASSERT(smi == "CC(C)(C)C");
    smi = MolToSmiles(*mol, true);
    TEST_ASSERT(smi == "C[13C](C)(C)C");
    delete mol;
  }
  {
    std::string smi = "C[12C](C)(C)C";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getAtomWithIdx(1)->getMass() == 12.0);
    smi = MolToSmiles(*mol, false);
    TEST_ASSERT(smi == "CC(C)(C)C");
    smi = MolToSmiles(*mol, true);
    TEST_ASSERT(smi == "C[12C](C)(C)C");
    delete mol;
  }
  {
    std::string smi = "CC[U]";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    smi = MolToSmiles(*mol, false);
    TEST_ASSERT(smi == "CC[U]");
    smi = MolToSmiles(*mol, true);
    TEST_ASSERT(smi == "CC[U]");
    delete mol;
  }
  {
    std::string smi = "CC[238U]";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    smi = MolToSmiles(*mol, false);
    TEST_ASSERT(smi == "CC[U]");
    smi = MolToSmiles(*mol, true);
    TEST_ASSERT(smi == "CC[238U]");
    delete mol;
  }
  {
    // issue 3526814
    std::string smi = "CCCCS(=[18O])(=O)CCCCl";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    smi = MolToSmiles(*mol, false);
    TEST_ASSERT(smi == "CCCCS(=O)(=O)CCCCl");
    smi = MolToSmiles(*mol, true);
    TEST_ASSERT(smi == "CCCCS(=O)(=[18O])CCCCl");
    delete mol;
  }
  {
    // issue 3526814
    std::string smi = "CCCCS(=[24O])(=O)CCCCl";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    smi = MolToSmiles(*mol, false);
    TEST_ASSERT(smi == "CCCCS(=O)(=O)CCCCl");
    smi = MolToSmiles(*mol, true);
    TEST_ASSERT(smi == "CCCCS(=O)(=[24O])CCCCl");
    delete mol;
  }
  {
    // issue 3526814
    std::string smi = "CCCCS(=O)(=[24O])CCCCl";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    smi = MolToSmiles(*mol, false);
    TEST_ASSERT(smi == "CCCCS(=O)(=O)CCCCl");
    smi = MolToSmiles(*mol, true);
    TEST_ASSERT(smi == "CCCCS(=O)(=[24O])CCCCl");
    delete mol;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testBug1670149() {
  RWMol *mol;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing SF.net bug 1670149" << std::endl;

  smi = "C1[NH2+]CCC1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  TEST_ASSERT(smi == "C1CC[NH2+]C1");

  mol->getAtomWithIdx(1)->setNumExplicitHs(0);
  mol->getAtomWithIdx(1)->setNoImplicit(false);
  mol->getAtomWithIdx(1)->updatePropertyCache();
  TEST_ASSERT(mol->getAtomWithIdx(1)->getNumImplicitHs() == 2);
  smi = MolToSmiles(*mol, false, false, -1);
  TEST_ASSERT(smi == "C1CC[NH2+]C1");
  delete mol;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testBug1719046() {
  RWMol *mol;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing SF.net bug 1719046: explicit Hs in canonical smiles"
      << std::endl;

  smi = "Cl[CH]1CCCCC1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  std::cerr << "smi: " << smi << std::endl;
  TEST_ASSERT(smi == "ClC1CCCCC1");

  delete mol;
  smi = "Cl[C@H]1CCCCC1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  std::cerr << "smi: " << smi << std::endl;
  TEST_ASSERT(smi == "ClC1CCCCC1");

  delete mol;
  smi = "Cl[C@H]1C(Br)CCCC1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  std::cerr << "smi: " << smi << std::endl;
  TEST_ASSERT(smi == "ClC1CCCCC1Br");

  delete mol;
  smi = "[CH]1=[CH][CH]=[CH][CH]=[CH]1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  std::cerr << "smi: " << smi << std::endl;
  TEST_ASSERT(smi == "c1ccccc1");

  delete mol;
  smi = "c1ccccn1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  std::cerr << "smi: " << smi << std::endl;
  TEST_ASSERT(smi == "c1ccncc1");

  delete mol;
  smi = "C1=CNC=C1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  std::cerr << "smi: " << smi << std::endl;
  TEST_ASSERT(smi == "c1cc[nH]c1");

  delete mol;
  smi = "[CH]1=[CH][NH][CH]=[CH]1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  std::cerr << "smi: " << smi << std::endl;
  TEST_ASSERT(smi == "c1cc[nH]c1");

  delete mol;
  // this was Issue 35525671
  smi = "P1C=CC=C1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol, false, false, -1);
  std::cerr << "smi: " << smi << std::endl;
  TEST_ASSERT(smi == "c1cc[pH]c1");

  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testBug1842174() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing SF.net bug 1842174: bad bond dirs in branches" << std::endl;
  RWMol *mol;
  std::string smi;

  smi = "F/C=N/Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol, true, false, -1);
  BOOST_LOG(rdInfoLog) << smi << std::endl;
  TEST_ASSERT(smi == "F/C=N/Cl");

  smi = MolToSmiles(*mol, true, false, 1);
  BOOST_LOG(rdInfoLog) << smi << std::endl;
  TEST_ASSERT(smi == "C(\\F)=N/Cl");

  delete mol;
  smi = "C(\\C=C\\F)=C(/Cl)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol, true, false, -1);
  BOOST_LOG(rdInfoLog) << smi << std::endl;
  TEST_ASSERT(smi == "F/C=C/C=C(/Cl)Br");

  smi = MolToSmiles(*mol, true, false, 0);
  BOOST_LOG(rdInfoLog) << smi << std::endl;
  TEST_ASSERT(smi == "C(/C=C/F)=C(\\Cl)Br");
  delete mol;

  smi = "O=NC1=NOC(=N\\O)/C1=N\\O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol, true, false, -1);
  BOOST_LOG(rdInfoLog) << smi << std::endl;
  TEST_ASSERT(smi == "O=NC1=NOC(=N\\O)/C1=N\\O");

  // ----------------------
  //  the next two examples are a pair:
  // vvvvvvvvvvvvvvvvvvvvvv
  delete mol;
  smi = "O/N=C/1COCC1=N\\O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol, true, false, -1);
  BOOST_LOG(rdInfoLog) << smi << std::endl;
  TEST_ASSERT(smi == "O/N=C1\\COC\\C1=N\\O");

  // this time the algorithm is forced to set
  // the directionality on the ring closure bond:
  delete mol;
  smi = "O/N=C/1COC[N+]1=N\\O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol, true, false, -1);
  BOOST_LOG(rdInfoLog) << smi << std::endl;
  TEST_ASSERT(smi == "O/N=C1\\COC\\[N+]1=N\\O");
  // ^^^^^^^^^^^^^^^^^^^^^^
  // end of the pair
  // ----------------------

  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testBug1844617() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing SF.net bug 1844617: oscillating chirality in canonical smiles"
      << std::endl;
  RWMol *mol;
  std::string smi, smi2;
  std::string label;

#if 0
  smi ="O=C1C2OCC[C@@]22C(CC1)CNCC2";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  //mol->debugMol(std::cout);
  TEST_ASSERT(mol->getAtomWithIdx(6)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(6)->getProp(common_properties::_CIPCode,label);
  TEST_ASSERT(label=="S");

  smi = MolToSmiles(*mol,true);
  BOOST_LOG(rdInfoLog) << smi << std::endl;
  delete mol;
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi2 = MolToSmiles(*mol,true);
  BOOST_LOG(rdInfoLog) << smi2 << std::endl;
  TEST_ASSERT(smi==smi2);

  delete mol;
#endif
  smi = "O=C1CC[C@@]2(O)[C@@H]3N(C)CC[C@]22[C@H]1OC[C@H]2CC3";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  // mol->debugMol(std::cout);
  MolOps::assignStereochemistry(*mol);
  // mol->debugMol(std::cout);
  TEST_ASSERT(mol->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
  TEST_ASSERT(mol->getAtomWithIdx(6)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(6)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "R");
  TEST_ASSERT(mol->getAtomWithIdx(11)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(11)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
  TEST_ASSERT(mol->getAtomWithIdx(12)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(12)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "R");
  TEST_ASSERT(mol->getAtomWithIdx(15)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(15)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
#if 1
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  // mol->debugMol(std::cout);
  smi2 = MolToSmiles(*mol, true);
  BOOST_LOG(rdInfoLog) << smi << std::endl;
  BOOST_LOG(rdInfoLog) << smi2 << std::endl;
  TEST_ASSERT(smi == smi2);
#endif

  delete mol;
  smi = "O=C1CC[C@@]2(O)[C@@H]3N(C)CC[C@]22[C@H]1OC[C@H]2CC3";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  // mol->debugMol(std::cout);
  MolOps::assignStereochemistry(*mol);
  // mol->debugMol(std::cout);
  TEST_ASSERT(mol->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
  TEST_ASSERT(mol->getAtomWithIdx(6)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(6)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "R");
  TEST_ASSERT(mol->getAtomWithIdx(11)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(11)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
  TEST_ASSERT(mol->getAtomWithIdx(12)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(12)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "R");
  TEST_ASSERT(mol->getAtomWithIdx(15)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(15)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
#if 1
  smi = MolToSmiles(*mol, true, false, 0);
  delete mol;
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  // mol->debugMol(std::cout);
  smi2 = MolToSmiles(*mol, true, false, 0);
  BOOST_LOG(rdInfoLog) << smi << std::endl;
  BOOST_LOG(rdInfoLog) << smi2 << std::endl;
  TEST_ASSERT(smi == smi2);
#endif

  delete mol;
  smi = "O=C1CC[C@@]2(O)[C@@H]3N(CC4CC4)CC[C@]22[C@H]1OC[C@H]2CC3";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  // mol->debugMol(std::cout);
  TEST_ASSERT(mol->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
  TEST_ASSERT(mol->getAtomWithIdx(6)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(6)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "R");
  TEST_ASSERT(mol->getAtomWithIdx(14)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(14)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
  TEST_ASSERT(mol->getAtomWithIdx(15)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(15)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "R");
  TEST_ASSERT(mol->getAtomWithIdx(18)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(18)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
#if 1
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  smi2 = MolToSmiles(*mol, true);
  BOOST_LOG(rdInfoLog) << smi << std::endl;
  BOOST_LOG(rdInfoLog) << smi2 << std::endl;
  TEST_ASSERT(smi == smi2);
#endif
  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testBug1844959() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing SF.net bug 1844959: bad handling of Hs in chiral smiles"
      << std::endl;
  RWMol *mol;
  std::string smi, smi2;
  std::string label;

  // ----------------------
  //  the next examples are a set:
  //  (this is the part that was originally working):
  // vvvvvvvvvvvvvvvvvvvvvv
  smi = "C[C@]12CNOC2.F1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "R");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "R");
  smi2 = MolToSmiles(*mol, true);
  TEST_ASSERT(smi == smi2);

  // swap the order and make sure the chirality swaps with it:
  delete mol;
  smi = "C[C@]12CNOC1.F2";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
  smi2 = MolToSmiles(*mol, true);
  TEST_ASSERT(smi == smi2);
  delete mol;

  // now make sure it works with a reversed chiral tag:
  smi = "C[C@@]12CNOC2.F1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
  smi2 = MolToSmiles(*mol, true);
  TEST_ASSERT(smi == smi2);
  delete mol;
  smi = "C[C@@]12CNOC1.F2";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "R");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "R");
  smi2 = MolToSmiles(*mol, true);
  TEST_ASSERT(smi == smi2);
  delete mol;
  // ^^^^^^^^^^^^^^^^^^^^^^
  // end of the set
  // ----------------------

  // ----------------------
  //  the next examples are a set:
  //  (this is the part that was originally failing):
  // vvvvvvvvvvvvvvvvvvvvvv
  BOOST_LOG(rdInfoLog) << "--------------------------------------------"
                       << std::endl;
  smi = "C[C@]12CNOC2.[H]1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  // mol->debugMol(std::cerr);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
  smi = MolToSmiles(*mol, true);
  BOOST_LOG(rdInfoLog) << smi << std::endl;
  delete mol;
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  // mol->debugMol(std::cerr);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
  smi2 = MolToSmiles(*mol, true);
  TEST_ASSERT(smi == smi2);

  // swap the order and make sure the chirality swaps with it:
  delete mol;
  smi = "C[C@]12CNOC1.[H]2";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "R");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "R");
  smi2 = MolToSmiles(*mol, true);
  TEST_ASSERT(smi == smi2);
  delete mol;

  // now make sure it works with a reversed chiral tag:
  smi = "C[C@@]12CNOC2.[H]1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "R");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "R");
  smi2 = MolToSmiles(*mol, true);
  TEST_ASSERT(smi == smi2);
  delete mol;
  smi = "C[C@@]12CNOC1.[H]2";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
  smi = MolToSmiles(*mol, true);
  delete mol;
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, label);
  TEST_ASSERT(label == "S");
  smi2 = MolToSmiles(*mol, true);
  TEST_ASSERT(smi == smi2);
  // ^^^^^^^^^^^^^^^^^^^^^^
  // end of the set
  // ----------------------

  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testBug1942220() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing sf.net bug 1942220" << std::endl;

  RWMol *m;
  std::string smi;

  smi = "[C](Cl)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 3);
  TEST_ASSERT(m->getNumAtoms(false) == 3);
  smi = MolToSmiles(*m);
  TEST_ASSERT(smi == "Cl[C]Br");

  delete m;
  smi = "[CH2](Cl)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 3);
  TEST_ASSERT(m->getNumAtoms(false) == 5);
  smi = MolToSmiles(*m);
  TEST_ASSERT(smi == "ClCBr");

  delete m;
  smi = "C(Cl)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 3);
  TEST_ASSERT(m->getNumAtoms(false) == 5);
  smi = MolToSmiles(*m);
  TEST_ASSERT(smi == "ClCBr");

  delete m;
  smi = "OS(=O)=O";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  // TEST_ASSERT(m->getNumAtoms(false)==5);
  smi = MolToSmiles(*m);
  TEST_ASSERT(smi == "O=[SH](=O)O");

  delete m;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testRingStereochemReporting() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing error reporting with ring stereochem"
                       << std::endl;

  RWMol *m;
  std::string smi;

  smi = "C[C@H]1CC[C@@H](C)CC1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 8);

  smi = MolToSmiles(*m, true);
  TEST_ASSERT(m->hasProp(common_properties::_ringStereoWarning));

  smi = MolToSmiles(*m, false);
  TEST_ASSERT((!m->hasProp(common_properties::_ringStereoWarning)));

  delete m;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testBug3127883() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue "
                          "3127883 (kekulization failing) "
                       << std::endl;
  {
    ROMol *m;
    std::string smi;
    smi = "c(:c:c:1):c:c:c:1";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    delete m;
  }

  {
    ROMol *m;
    std::string smi;
    smi = "c1(:c(:c(:c(-C(-c2:c(:c(:c(:c(:c:2)))))=C):c(:c:1))))";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testBug3139534() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Issue 3139534: stereochemistry in larger rings"
                       << std::endl;

  // the parsing part of this is in ../testChirality.cpp, here we look at
  // smiles generation

  {
    RWMol *m;
    std::string smiles = "C1COC/C=C\\CCC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(4)->getStereo() == Bond::STEREOZ);

    smiles = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << "smiles: " << smiles << std::endl;
    TEST_ASSERT(smiles == "C1=C\\COCCCCC/1");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C1COC/C=C/CCC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(4)->getStereo() == Bond::STEREOE);

    smiles = MolToSmiles(*m, true);
    TEST_ASSERT(smiles == "C1=C/COCCCCC/1");

    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C1CC/C=C/C=C/CCC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    smiles = MolToSmiles(*m, true, false, -1, false);
    // TEST_ASSERT(smiles=="C1CC/C=C/C=C/CCC1");
    TEST_ASSERT(smiles == "C1CC/C=C/C=C/CCC1");

    smiles = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << "smiles: " << smiles << std::endl;
    TEST_ASSERT(smiles == "C1=C/CCCCCC/C=C/1");
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "C/1=C/C=C/CCCCCC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    smiles = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << "smiles: " << smiles << std::endl;
    TEST_ASSERT(smiles == "C1=C\\CCCCCC/C=C/1");
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "C1COC/C=C/C=C/C1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(4)->getStereo() == Bond::STEREOE);

    smiles = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << "smiles: " << smiles << std::endl;
    TEST_ASSERT(smiles == "C1=C/CCCOC/C=C/1");

    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "C1=C/OCC/C=C\\CC\\1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "C1CCCCN/C=C/1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    smiles = MolToSmiles(*m, true, false, 7, false);
    BOOST_LOG(rdInfoLog) << "smiles: " << smiles << std::endl;
    TEST_ASSERT(smiles == "C1=C/NCCCCC/1");

    smiles = MolToSmiles(*m, true, false, 0, false);
    BOOST_LOG(rdInfoLog) << "smiles: " << smiles << std::endl;
    TEST_ASSERT(smiles == "C1CCCCN/C=C/1");

    delete m;
  }

  {
    RWMol *m;
    // the 2 initial directed bonds are redundant (/bad ??)
    std::string smiles = "CCC/[N+]/1=C/c2ccccc2OC(=O)/C=C1/O";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getBondWithIdx(3)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(14)->getStereo() == Bond::STEREOE);

    smiles = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << "smiles: " << smiles << std::endl;
    TEST_ASSERT(smiles == R"(CCC[N+]1=C/c2ccccc2OC(=O)/C=C\1O)");

    delete m;

    // 2nd pass to check stability
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getBondWithIdx(3)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(14)->getStereo() == Bond::STEREOE);

    smiles = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << "smiles: " << smiles << std::endl;
    TEST_ASSERT(smiles == R"(CCC[N+]1=C/c2ccccc2OC(=O)/C=C\1O)");

    delete m;
  }

  {  // Github #2023
    RWMol *m;
    // the initial directed bond is redundant
    std::string smiles = R"(CO/C1=C/C=C\C=C/C=N\1)";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getBondWithIdx(2)->getStereo() == Bond::STEREOE);
    TEST_ASSERT(m->getBondWithIdx(4)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(6)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(8)->getStereo() == Bond::STEREOZ);

    smiles = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << "smiles: " << smiles << std::endl;
    TEST_ASSERT(smiles == R"(COC1=C/C=C\C=C/C=N\1)");

    delete m;

    // 2nd pass to check stability
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getBondWithIdx(2)->getStereo() == Bond::STEREOE);
    TEST_ASSERT(m->getBondWithIdx(4)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(6)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(8)->getStereo() == Bond::STEREOZ);

    smiles = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << "smiles: " << smiles << std::endl;
    TEST_ASSERT(smiles == R"(COC1=C/C=C\C=C/C=N\1)");

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
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(30, 32)->getStereo() == Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(33, 34)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondBetweenAtoms(5, 7)->getStereo() == Bond::STEREOE);

    std::string csmiles = MolToSmiles(*m, true);

    RWMol *m2;
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      std::string nsmiles = MolToSmiles(*m, true, false, i, false);
      m2 = SmilesToMol(nsmiles);
      TEST_ASSERT(m2);
      std::string ncsmiles = MolToSmiles(*m2, true);
      if (ncsmiles != csmiles) {
        std::cerr << " failed in iteration: " << i << "\n"
                  << csmiles << "\n != \n"
                  << ncsmiles << "\n starting from:\n"
                  << nsmiles << "\n";
        m2->debugMol(std::cerr);
        TEST_ASSERT(ncsmiles == csmiles);
      }
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
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(13, 15)->getStereo() == Bond::STEREOZ);

    std::string csmiles = MolToSmiles(*m, true);

    RWMol *m2;
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      std::string nsmiles = MolToSmiles(*m, true, false, i, false);
      m2 = SmilesToMol(nsmiles);
      TEST_ASSERT(m2);
      std::string ncsmiles = MolToSmiles(*m2, true);
      if (ncsmiles != csmiles) {
        std::cerr << " failed in iteration: " << i << "\n"
                  << csmiles << "\n != \n"
                  << ncsmiles << "\n starting from:\n"
                  << nsmiles << "\n";
        m2->debugMol(std::cerr);
        TEST_ASSERT(ncsmiles == csmiles);
      }
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
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(13, 15)->getStereo() == Bond::STEREOE);

    std::string csmiles = MolToSmiles(*m, true);

    RWMol *m2;
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      std::string nsmiles = MolToSmiles(*m, true, false, i, false);
      m2 = SmilesToMol(nsmiles);
      TEST_ASSERT(m2);
      std::string ncsmiles = MolToSmiles(*m2, true);
      if (ncsmiles != csmiles) {
        std::cerr << " failed in iteration: " << i << "\n"
                  << csmiles << "\n != \n"
                  << ncsmiles << "\n starting from:\n"
                  << nsmiles << "\n";
        m2->debugMol(std::cerr);
        TEST_ASSERT(ncsmiles == csmiles);
      }
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
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(15, 16)->getStereo() == Bond::STEREOZ);

    std::string csmiles = MolToSmiles(*m, true);

    RWMol *m2;
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      std::string nsmiles = MolToSmiles(*m, true, false, i, false);
      m2 = SmilesToMol(nsmiles);
      TEST_ASSERT(m2);
      std::string ncsmiles = MolToSmiles(*m2, true);
      if (ncsmiles != csmiles) {
        std::cerr << " failed in iteration: " << i << "\n"
                  << csmiles << "\n != \n"
                  << ncsmiles << "\n starting from:\n"
                  << nsmiles << "\n";
        m2->debugMol(std::cerr);
        TEST_ASSERT(ncsmiles == csmiles);
      }
      delete m2;
    }
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testAtomMaps() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test adding atom-map information" << std::endl;

  {
    RWMol *m;
    std::string smiles = "[*:1]CCC([C:200])C";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(
        m->getAtomWithIdx(0)->hasProp(common_properties::molAtomMapNumber));

    // changed: smiles does not need to be canonical
    smiles = MolToSmiles(*m, true, false, -1, false);
    TEST_ASSERT(smiles == "[*:1]CCC([C:200])C");

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testBug3145697() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Issue 3145697 repeated ring labels in disconnected structures"
      << std::endl;

  {
    RWMol *m;
    std::string smiles = "C1.C11.C1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    smiles = MolToSmiles(*m, true);
    TEST_ASSERT(smiles == "CCC");
    delete m;

    smiles = "C1.C11.C";
    m = SmilesToMol(smiles);
    TEST_ASSERT(!m);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C1.C11.O1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    smiles = MolToSmiles(*m, true);
    TEST_ASSERT(smiles == "CCO");
    delete m;

    smiles = "C1.C1=1.O1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    smiles = MolToSmiles(*m, true);
    TEST_ASSERT(smiles == "CC=O");
    delete m;

    smiles = "C1.C=11.O1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    smiles = MolToSmiles(*m, true);
    TEST_ASSERT(smiles == "C=CO");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C1C.CC11CCC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    smiles = MolToSmiles(*m, true);
    TEST_ASSERT(smiles == "CCC1(C)CCC1");
    delete m;
    smiles = "C1C.CC11CCC";
    m = SmilesToMol(smiles);
    TEST_ASSERT(!m);
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testBug3152751() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Issue 3152751 cannot roundtrip charged aromatic Se and Te"
      << std::endl;

  {
    RWMol *m;
    std::string smiles = "c1cccc[te+]1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
    smiles = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "c1cccc[se+]1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
    smiles = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "c1ccc[te]1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
    smiles = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "c1ccc[se]1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
    smiles = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testReplacementPatterns() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing use of replacement patterns in input"
                       << std::endl;

  {
    std::string smi = "C{cycloprop}C";
    std::map<std::string, std::string> repls;
    repls["{cycloprop}"] = "C1(CC1)";
    RWMol *mol = SmilesToMol(smi, 0, true, &repls);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 5);
    TEST_ASSERT(mol->getAtomWithIdx(1)->getDegree() == 4);
    delete mol;
  }

  {
    std::string smi = "C{cycloprop}C";
    std::map<std::string, std::string> repls;
    repls["{cycloprop}"] = "C1(C({acid})C1)";
    repls["{acid}"] = "C(=O)O";
    RWMol *mol = SmilesToMol(smi, 0, true, &repls);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 8);
    TEST_ASSERT(mol->getAtomWithIdx(1)->getDegree() == 4);
    delete mol;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testAllBondsExplicit() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing forcing explicit bonds in the output SMILES"
                       << std::endl;

  {
    std::string smi = "CCC";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 3);
    smi = MolToSmiles(*mol, true);
    TEST_ASSERT(smi == "CCC");
    smi = MolToSmiles(*mol, true, false, -1, true, true);
    TEST_ASSERT(smi == "C-C-C");

    delete mol;
  }
  {
    std::string smi = "C1CC1";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 3);
    smi = MolToSmiles(*mol, true);
    TEST_ASSERT(smi == "C1CC1");
    smi = MolToSmiles(*mol, true, false, -1, true, true);
    TEST_ASSERT(smi == "C1-C-C-1");

    delete mol;
  }
  {
    std::string smi = "c1ccccc1";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 6);
    smi = MolToSmiles(*mol, true);
    TEST_ASSERT(smi == "c1ccccc1");
    smi = MolToSmiles(*mol, true, false, -1, true, true);
    TEST_ASSERT(smi == "c1:c:c:c:c:c:1");

    delete mol;
  }
  {
    std::string smi = "c1ccccc1c1ccccc1";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 12);
    smi = MolToSmiles(*mol, true);
    TEST_ASSERT(smi == "c1ccc(-c2ccccc2)cc1");
    smi = MolToSmiles(*mol, true, false, -1, true, true);
    TEST_ASSERT(smi == "c1:c:c:c(-c2:c:c:c:c:c:2):c:c:1");

    delete mol;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testBug3525799() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Issue 3525799: bad smiles for r groups" << std::endl;

  {
    RWMol *m;
    std::string smiles = "CC*";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    smiles = MolToSmiles(*m, true);
    TEST_ASSERT(smiles == "*CC");
    m->getAtomWithIdx(2)->setProp(common_properties::dummyLabel, "foo");
    smiles = MolToSmiles(*m, true);
    TEST_ASSERT(smiles == "*CC");
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "CC*";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    smiles = MolToSmiles(*m, true);
    TEST_ASSERT(smiles == "*CC");
    m->getAtomWithIdx(2)->setProp(common_properties::smilesSymbol, "Xa");
    smiles = MolToSmiles(*m, true);
    TEST_ASSERT(smiles == "[Xa]CC");
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testBug3526810() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Issue 3526810: canonical smiles failure in symmetric heterocycles"
      << std::endl;

  {
    RWMol *m;
    std::string smiles = "C1SCCSCCCSCCSCC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmiles1 = MolToSmiles(*m, true);
    delete m;
    std::string smiles2 = "C1CSCCSCCCSCCSC1";
    m = SmilesToMol(smiles2);
    TEST_ASSERT(m);
    std::string csmiles2 = MolToSmiles(*m, true);
    delete m;

    // std::cerr<<"csmi1: "<<csmiles1<<std::endl;
    // std::cerr<<"csmi2: "<<csmiles2<<std::endl;
    TEST_ASSERT(csmiles1 == csmiles2);
  }

  {
    RWMol *m;
    std::string smiles = "C1NCCNCCCNCCNCC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmiles1 = MolToSmiles(*m, true);
    delete m;
    std::string smiles2 = "C1CNCCNCCCNCCNC1";
    m = SmilesToMol(smiles2);
    TEST_ASSERT(m);
    std::string csmiles2 = MolToSmiles(*m, true);
    delete m;

    // std::cerr<<"csmi1: "<<csmiles1<<std::endl;
    // std::cerr<<"csmi2: "<<csmiles2<<std::endl;
    TEST_ASSERT(csmiles1 == csmiles2);
  }

  {
    RWMol *m;
    std::string smiles = "C1CNCCCNCCNCCCNC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmiles1 = MolToSmiles(*m, true);
    delete m;
    std::string smiles2 = "C1CCNCCCNCCNCCCN1";
    m = SmilesToMol(smiles2);
    TEST_ASSERT(m);
    std::string csmiles2 = MolToSmiles(*m, true);
    delete m;

    // std::cerr<<"csmi1: "<<csmiles1<<std::endl;
    // std::cerr<<"csmi2: "<<csmiles2<<std::endl;
    TEST_ASSERT(csmiles1 == csmiles2);
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testBug3526815() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Issue 3526815: canonical smiles failure in many symmetric fragments"
      << std::endl;

  {
    RWMol *m;
    std::string smiles =
        "O.O.O.O.O.O.O.O.O.[Pd].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+]"
        ".[Na+].[O-]S(=O)(=O)c1cccc(c1)P(c1cccc(c1)S(=O)(=O)[O-])c1cccc(c1)S(="
        "O)(=O)[O-].[O-]S(=O)(=O)c1cccc(c1)P(c1cccc(c1)S(=O)(=O)[O-])c1cccc(c1)"
        "S(=O)(=O)[O-].[O-]S(=O)(=O)c1cccc(c1)P(c1cccc(c1)S(=O)(=O)[O-])c1cccc("
        "c1)S(=O)(=O)[O-]";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmiles1 = MolToSmiles(*m, true);
    delete m;
    std::string smiles2 =
        "O.O.O.O.O.O.O.O.O.[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+].[Na+"
        "].[O-]S(c1cccc(P(c2cccc(S([O-])(=O)=O)c2)c2cccc(S([O-])(=O)=O)c2)c1)(="
        "O)=O.[Pd].[O-]S(=O)(=O)c1cccc(P(c2cccc(S([O-])(=O)=O)c2)c2cccc(S([O-])"
        "(=O)=O)c2)c1.[O-]S(=O)(=O)c1cccc(P(c2cccc(S([O-])(=O)=O)c2)c2cccc(S(["
        "O-])(=O)=O)c2)c1";
    m = SmilesToMol(smiles2);
    TEST_ASSERT(m);
    std::string csmiles2 = MolToSmiles(*m, true);
    delete m;

    // std::cerr<<"csmi1: "<<csmiles1<<std::endl;
    // std::cerr<<"csmi2: "<<csmiles2<<std::endl;
    TEST_ASSERT(csmiles1 == csmiles2);
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testFragmentSmiles() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Fragment Smiles" << std::endl;
  {
    RWMol *m;
    std::string smiles = "OCCCC";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {0, 1, 2};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse);
    TEST_ASSERT(csmiles == "CCO");
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "OCCCCCC";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {0, 1, 2, 3};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse);
    TEST_ASSERT(csmiles == "CCCO");
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "OC1CC1CCC";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {1, 2, 3};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse);
    TEST_ASSERT(csmiles == "C1CC1");
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "OC1CC1CCC";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {1, 2, 3};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    int bs[] = {1, 2, 6};
    std::vector<int> bondsToUse(bs, bs + sizeof(bs) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse, &bondsToUse);
    TEST_ASSERT(csmiles == "C1CC1");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "OC1CC1CCC";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {1, 2, 3};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    int bs[] = {1, 2};
    std::vector<int> bondsToUse(bs, bs + sizeof(bs) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse, &bondsToUse);
    TEST_ASSERT(csmiles == "CCC");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "OC1CCCCC1N";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {1, 2, 3, 4, 5, 6};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse);
    TEST_ASSERT(csmiles == "C1CCCCC1");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "OCCCCCCN";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {1, 2, 3, 4, 5, 6};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse);
    TEST_ASSERT(csmiles == "CCCCCC");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "OCCCCCCN";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {1, 2, 3, 4, 5, 6};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    int bs[] = {1, 2, 3, 4, 5};
    std::vector<int> bondsToUse(bs, bs + sizeof(bs) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse, &bondsToUse);

    TEST_ASSERT(csmiles == "CCCCCC");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "OC1CCCCC1N";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {1, 2, 3, 4, 5, 6};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    int bs[] = {1, 2, 3, 4, 5};
    std::vector<int> bondsToUse(bs, bs + sizeof(bs) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse, &bondsToUse);
    TEST_ASSERT(csmiles == "CCCCCC");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "Oc1ccccc1N";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {1, 2, 3, 4, 5, 6};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse);
    TEST_ASSERT(csmiles == "c1ccccc1");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "Oc1ccccc1N";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {1, 2, 3, 4, 5, 6};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    int bs[] = {1, 2, 3, 4, 5};
    std::vector<int> bondsToUse(bs, bs + sizeof(bs) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse, &bondsToUse);
    TEST_ASSERT(csmiles == "cccccc");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "OCCCC";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {0, 1, 2};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string labels[5] = {"[A]", "[B]", "[B]", "", ""};
    std::vector<std::string> atomLabels(labels, labels + 5);
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, &atomLabels);
    TEST_ASSERT(csmiles == "[A][B][B]");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CCCCO";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {2, 3, 4};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string labels[5] = {"", "", "[B]", "[B]", "[A]"};
    std::vector<std::string> atomLabels(labels, labels + 5);
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, &atomLabels);
    TEST_ASSERT(csmiles == "[A][B][B]");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CCCCO";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {2, 3, 4};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string labels[5] = {"", "", "[B]", "[A]", "[B]"};
    std::vector<std::string> atomLabels(labels, labels + 5);
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, &atomLabels);
    TEST_ASSERT(csmiles == "[B][A][B]");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CC(=O)OCC";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {0, 1, 2, 3};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, nullptr, nullptr);
    TEST_ASSERT(csmiles == "CC(=O)O");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CC(=O)OCC";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {0, 1, 2, 3};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string labels[5] = {"-", "=", "-", "", ""};
    std::vector<std::string> bondLabels(labels, labels + 5);
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, nullptr, &bondLabels);
    TEST_ASSERT(csmiles == "C-C(=O)-O");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CC(=O)OCC";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {0, 1, 2, 3};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string labels[5] = {"a", "b", "a", "", ""};
    std::vector<std::string> bondLabels(labels, labels + 5);
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, nullptr, &bondLabels);
    TEST_ASSERT(csmiles == "CaC(bO)aO");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CC(=CC)CCC";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {0, 1, 2, 4};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse);
    TEST_ASSERT(csmiles == "C=C(C)C");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CC(=CC)CCC";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {0, 1, 2, 4};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string labels[6] = {"a", "b", "", "a", "", ""};
    std::vector<std::string> bondLabels(labels, labels + 6);
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, nullptr, &bondLabels);
    std::cerr << csmiles << std::endl;
    TEST_ASSERT(csmiles == "CbC(aC)aC");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CC(=CC)CCC";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {0, 1, 2, 4};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string labels[6] = {"b", "a", "", "a", "", ""};
    std::vector<std::string> bondLabels(labels, labels + 6);
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, nullptr, &bondLabels);
    std::cerr << csmiles << std::endl;
    TEST_ASSERT(csmiles == "CaC(bC)aC" || csmiles == "CbC(ac)ac");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CC(=CC)CCC";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {0, 1, 2, 4};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string labels[6] = {"b", "b", "", "a", "", ""};
    std::vector<std::string> bondLabels(labels, labels + 6);
    std::string csmiles =
        MolFragmentToSmiles(*m, atomsToUse, nullptr, nullptr, &bondLabels);
    std::cerr << csmiles << std::endl;
    TEST_ASSERT(csmiles == "CaC(bC)bC" || csmiles == "CbC(bC)aC");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "OC1CC1CC";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {0, 4};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles = MolFragmentToSmiles(*m, atomsToUse, nullptr, nullptr,
                                              nullptr, false, false, -1, false);
    std::cerr << csmiles << std::endl;
    TEST_ASSERT(csmiles == "O.C");
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testBug3528556() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Issue 3528556: canonical smiles failure in cycle"
                       << std::endl;

  {
    RWMol *m;
    std::string smiles = "N12.N13.C24.C35.C46.C56";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmiles1 = MolToSmiles(*m, true);
    delete m;
    std::string smiles2 = "N1NCCCC1";
    m = SmilesToMol(smiles2);
    TEST_ASSERT(m);
    std::string csmiles2 = MolToSmiles(*m, true);
    delete m;
    TEST_ASSERT(csmiles1 == csmiles2);
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testBug253() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "do not repeat ring closure digits on the same atom"
                       << std::endl;

  {
    RWMol *m;
    std::string smiles = "C1CCCC1CCC1CCCCC11CCCCC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmiles1 = MolToSmiles(*m, true);
    std::cerr << "--" << csmiles1 << std::endl;
    TEST_ASSERT(csmiles1 == "C1CCC2(CC1)CCCCC2CCC1CCCC1");
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testBug257() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Issue 257: unrecognized bonds are in SMILES as ?s"
                       << std::endl;

  {
    RWMol *m;
    std::string smiles = "CCO";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    m->getBondWithIdx(1)->setBondType(Bond::UNSPECIFIED);
    std::string csmiles = MolToSmiles(*m);
    TEST_ASSERT(csmiles == "CC~O");
    delete m;
    m = SmilesToMol(csmiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::UNSPECIFIED);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub12() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github 12: non-canonical fragment smiles"
                       << std::endl;
  {
    RWMol *m;
    std::string smiles = "c1c(C)cccc1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int as[] = {0, 1, 2};
    std::vector<int> atomsToUse(as, as + sizeof(as) / sizeof(int));
    std::string csmiles1 = MolFragmentToSmiles(*m, atomsToUse);
    int as2[] = {1, 2, 3};
    std::vector<int> atomsToUse2(as2, as2 + sizeof(as2) / sizeof(int));
    std::string csmiles2 = MolFragmentToSmiles(*m, atomsToUse2);
    TEST_ASSERT(csmiles1 == csmiles2);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testRingStereochem() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing handling of ring stereochemistry"
                       << std::endl;

  {
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
      ;
      TEST_ASSERT(m1);
      RWMol *m2 = SmilesToMol(smi2);
      ;
      TEST_ASSERT(m2);
      TEST_ASSERT(m1->getNumAtoms() == m2->getNumAtoms());
      TEST_ASSERT(m1->getNumBonds() == m2->getNumBonds());

      std::string csmiles1 = MolToSmiles(*m1, true);
      std::string csmiles2 = MolToSmiles(*m2, true);
      if (csmiles1 != csmiles2) {
        std::cerr << "---------\n" << csmiles1 << "\n" << csmiles2 << std::endl;
      }
      TEST_ASSERT(csmiles1 == csmiles2);
      delete m1;
      delete m2;
    }
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub45() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github 45: stereochemistry information "
                          "influencing non-stereo SMILES"
                       << std::endl;
  {
    RWMol *m;
    std::string smiles = "CC1CCC[13C]2(C)C1CC[14CH]2C(C)=O";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmiles1a = MolToSmiles(*m, false);
    std::string csmiles1b = MolToSmiles(*m, true);
    std::string smiles2 = "CC1CCC[C]2(C)C1CC[CH]2C(C)=O";
    delete m;
    m = SmilesToMol(smiles2);
    TEST_ASSERT(m);
    std::string csmiles2a = MolToSmiles(*m, false);
    std::string csmiles2b = MolToSmiles(*m, true);

    TEST_ASSERT(csmiles1a == csmiles2a);
    TEST_ASSERT(csmiles1b != csmiles2b);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "CC1CCC[C@@]2(C)C1CC[C@@H]2C(C)=O";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmiles1a = MolToSmiles(*m, false);
    std::string csmiles1b = MolToSmiles(*m, true);
    std::string smiles2 = "CC1CCC[C]2(C)C1CC[CH]2C(C)=O";
    delete m;
    m = SmilesToMol(smiles2);
    TEST_ASSERT(m);
    std::string csmiles2a = MolToSmiles(*m, false);
    std::string csmiles2b = MolToSmiles(*m, true);

    TEST_ASSERT(csmiles1a == csmiles2a);
    TEST_ASSERT(csmiles1b != csmiles2b);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub206() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github 206: Problems round-tripping P"
                       << std::endl;
  {
    RWMol *m;
    std::string smiles = "O=[PH3]";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmiles = MolToSmiles(*m, true);
    TEST_ASSERT(csmiles == "O=[PH3]");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "O=P";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmiles = MolToSmiles(*m, true);
    TEST_ASSERT(csmiles == "O=P");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "O=[PH]";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmiles = MolToSmiles(*m, true);
    TEST_ASSERT(csmiles == "O=P");
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}
void testGithub210() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github 210: flag possible stereocenters "
                          "when calling assignStereochemistry()"
                       << std::endl;
  {
    RWMol *m;
    std::string smiles = "O[C@H](F)CC(F)(Cl)I";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(
        m->getAtomWithIdx(4)->hasProp(common_properties::_ChiralityPossible));
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}
void testGithub298() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing Github 298: cannot generate smiles for ChEBI_50252"
      << std::endl;
  {
    std::string rdbase = getenv("RDBASE");
    std::string fName = rdbase + "/Code/GraphMol/test_data/ChEBI_50252.mol";
    RWMol *m = MolFileToMol(fName, false, false);

    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 80);
    TEST_ASSERT(m->getNumBonds() == 210);
    m->updatePropertyCache(false);
    MolOps::fastFindRings(*m);

    std::string csmiles = MolToSmiles(*m);
    TEST_ASSERT(csmiles != "");
    TEST_ASSERT(csmiles.find("%100") == std::string::npos);

    delete m;
    m = SmilesToMol(csmiles, 0, false);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 80);
    TEST_ASSERT(m->getNumBonds() == 210);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub378() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github 378: SMILES parser doing the wrong "
                          "thing for odd dot-disconnected construct"
                       << std::endl;
  {
    RWMol *m;
    std::string smiles = "C1.C1CO1.N1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1));
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondBetweenAtoms(3, 4));
    TEST_ASSERT(m->getBondBetweenAtoms(3, 4)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(!m->getBondBetweenAtoms(1, 3));
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C1(O.C1)CO1.N1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 2));
    TEST_ASSERT(m->getBondBetweenAtoms(0, 2)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 3));
    TEST_ASSERT(m->getBondBetweenAtoms(0, 3)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondBetweenAtoms(5, 4));
    TEST_ASSERT(m->getBondBetweenAtoms(5, 4)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(!m->getBondBetweenAtoms(2, 3));
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub389() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github 389: Add option to SmilesWriter to "
                          "allow writing of all explicit hydrogens"
                       << std::endl;
  {
    RWMol *m;
    std::string smiles = "CCO";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    std::string csmiles = MolToSmiles(*m, true, false, -1, true, false, true);
    TEST_ASSERT(csmiles != "");
    TEST_ASSERT(csmiles.find("[CH3]") != std::string::npos);
    TEST_ASSERT(csmiles.find("[CH2]") != std::string::npos);
    TEST_ASSERT(csmiles.find("[OH]") != std::string::npos);

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testEmptyStrings() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing handling of empty SMILES/SMARTS strings"
                       << std::endl;
  {
    RWMol *m;
    std::string smiles = "";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 0);

    std::string csmiles = MolToSmiles(*m);
    TEST_ASSERT(csmiles == "");
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "";
    m = SmartsToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 0);

    std::string csmiles = MolToSmarts(*m);
    TEST_ASSERT(csmiles == "");
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testSmilesWriteForModifiedMolecules() {
  BOOST_LOG(rdInfoLog)
      << "testing smiles writing/canonicalization for modified molecules."
      << std::endl;
  {
    std::string smiles = "c1ccccc1";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    m->getAtomWithIdx(0)->setAtomicNum(8);
    std::string smi = MolToSmiles(*m, true);
    // std::cerr<< smi <<std::endl;
    TEST_ASSERT(smi == "c1ccocc1");
    delete m;
  }
}

void testGithub532() {
  BOOST_LOG(rdInfoLog) << "testing github issue 532: _smilesAtomOutputOrder "
                          "incorrect for dot disconnected molecules"
                       << std::endl;
  {
    std::string smiles = "O.CO";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    std::string smi = MolToSmiles(*m, true);
    TEST_ASSERT(smi == "CO.O");

    std::vector<unsigned int> atmOrder;
    TEST_ASSERT(m->hasProp(common_properties::_smilesAtomOutputOrder));
    m->getProp(common_properties::_smilesAtomOutputOrder, atmOrder);
    TEST_ASSERT(atmOrder.size() == 3);
    TEST_ASSERT(atmOrder[0] == 1);
    TEST_ASSERT(atmOrder[1] == 2);
    TEST_ASSERT(atmOrder[2] == 0);

    delete m;
  }
  {
    std::string smiles = "CO.O";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    std::string smi = MolToSmiles(*m, true);
    TEST_ASSERT(smi == "CO.O");

    std::vector<unsigned int> atmOrder;
    TEST_ASSERT(m->hasProp(common_properties::_smilesAtomOutputOrder));
    m->getProp(common_properties::_smilesAtomOutputOrder, atmOrder);
    TEST_ASSERT(atmOrder.size() == 3);
    TEST_ASSERT(atmOrder[0] == 0);
    TEST_ASSERT(atmOrder[1] == 1);
    TEST_ASSERT(atmOrder[2] == 2);

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub760() {
  BOOST_LOG(rdInfoLog) << "testing github issue 760: reversed stereochemistry "
                          "with sulfoxides and ring closures"
                       << std::endl;
  {
    std::string smiles = "C[S@](Cl)=O";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;

    smiles = "C[S@]2=O.Cl2";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi2 = MolToSmiles(*m, true);
    TEST_ASSERT(csmi == csmi2);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub786() {
  BOOST_LOG(rdInfoLog) << "testing github issue 786: chiral order for "
                          "ring closure after branch"
                       << std::endl;
  {
    std::string smiles = "C1CN[C@H]1O";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;

    smiles = "C1CN[C@@H](O)1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi2 = MolToSmiles(*m, true);
    TEST_ASSERT(csmi == csmi2);
    delete m;
  }
  {
    std::string smiles = "C1CN[C@]1(O)N";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;

    smiles = "C1CN[C@](O)(N)1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi2 = MolToSmiles(*m, true);
    TEST_ASSERT(csmi == csmi2);
    delete m;
  }
  {
    std::string smiles = "C1CN[C@]12(O).N2";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;

    smiles = "C1CN[C@](O)12.N2";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi2 = MolToSmiles(*m, true);
    TEST_ASSERT(csmi == csmi2);
    delete m;

    smiles = "C1CN[C@@]1(O)2.N2";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    csmi2 = MolToSmiles(*m, true);
    TEST_ASSERT(csmi == csmi2);
    delete m;

    smiles = "C1CN[C@]2(O)1.N2";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    csmi2 = MolToSmiles(*m, true);
    TEST_ASSERT(csmi == csmi2);
    delete m;
  }
  {
    std::string smiles = "C[C@]1(O)NCC1";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;

    smiles = "C[C@@](O)1NCC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi2 = MolToSmiles(*m, true);
    TEST_ASSERT(csmi == csmi2);
    delete m;
  }
  {
    std::string smiles = "C[C@]1(NCC1)O";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;
    // so many pathologically ugly SMILES:
    smiles = "C[C@](NCC1)(O)1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi2 = MolToSmiles(*m, true);
    TEST_ASSERT(csmi == csmi2);
    delete m;
  }

  {  // Andrew's original real example:
    std::string smiles =
        "CC(C)[C@]1(N)CC[C@]2([C@@H](O2)CCC(=C)[C@H](CC[C@@](/C=C1)(C)O)O)C";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;
    smiles =
        "CC(C)[C@@](N)1CC[C@]2([C@@H](O2)CCC(=C)[C@H](CC[C@@](/C=C1)(C)O)O)C";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi2 = MolToSmiles(*m, true);
    std::cerr << csmi << " " << csmi2 << std::endl;
    TEST_ASSERT(csmi == csmi2);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub1652() {
  BOOST_LOG(rdInfoLog)
      << "testing github issue 1652: chiral order for "
         "ring closure after branch for the first atom in the SMILES string"
      << std::endl;
  {
    std::string smiles = "Cl[C@](F)1CC[C@H](F)CC1";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;

    smiles = "[C@](Cl)(F)1CC[C@H](F)CC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi2 = MolToSmiles(*m, true);
    TEST_ASSERT(csmi == csmi2);
    delete m;
  }
  {
    std::string smiles = "F[C@@]1(C)CCO1";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;

    smiles = "[C@@](F)1(C)CCO1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string csmi2 = MolToSmiles(*m, true);
    TEST_ASSERT(csmi == csmi2);
    delete m;
  }
}

void testDativeBonds() {
  BOOST_LOG(rdInfoLog) << "testing dative bond support" << std::endl;
  {
    std::string smiles = "CCC(=O)O->[Cu]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    int dative_bond_count = 0;
    for (size_t i = 0; i < m->getNumBonds(); i++) {
      if (m->getBondWithIdx(i)->getBondType() == Bond::DATIVE) {
        dative_bond_count++;
      }
    }
    TEST_ASSERT(dative_bond_count == 1);

    std::string out_smiles = MolToSmiles(*m, true);
    delete m;
    TEST_ASSERT(out_smiles == smiles);
  }
  {
    std::string smiles = "CCC(=O)O->[Cu]<-OC(O)CC";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    int dative_bond_count = 0;
    for (size_t i = 0; i < m->getNumBonds(); i++) {
      if (m->getBondWithIdx(i)->getBondType() == Bond::DATIVE) {
        dative_bond_count++;
      }
    }
    TEST_ASSERT(dative_bond_count == 2);

    std::string out_smiles = MolToSmiles(*m, true);
    delete m;
    TEST_ASSERT(out_smiles == smiles);
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub1219() {
  BOOST_LOG(rdInfoLog)
      << "Stereochemistry not output to SMILES when allHsExplicit=True"
      << std::endl;
  {
    std::string smiles = "C[C@H](F)Cl";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    bool doIsomericSmiles = true;
    bool doKekule = false;
    int rootedAtAtom = -1;
    bool canonical = true, allBondsExplicit = false, allHsExplicit = true;
    std::string csmi = MolToSmiles(*m, doIsomericSmiles, doKekule, rootedAtAtom,
                                   canonical, allBondsExplicit, allHsExplicit);
    TEST_ASSERT(csmi == "[CH3][C@H]([F])[Cl]");
    delete m;
  }
  {  // another manifestation was that chiral flags were not output for atoms
     // not in the organic subset
    std::string smiles = "C[Si@H](F)Cl";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    bool doIsomericSmiles = true;
    std::string csmi = MolToSmiles(*m, doIsomericSmiles);
    TEST_ASSERT(csmi == "C[Si@H](F)Cl");
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testSmilesParseParams() {
  BOOST_LOG(rdInfoLog) << "Testing the SmilesParseParams class" << std::endl;
  {
    std::string smiles = "C1=CC=CC=C1[H]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 6);
    TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
    delete m;

    {
      SmilesParserParams params;
      m = SmilesToMol(smiles, params);
      TEST_ASSERT(m);
      TEST_ASSERT(m->getNumAtoms() == 6);
      TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
      delete m;
    }
    {  // no removeHs, with sanitization
      SmilesParserParams params;
      params.removeHs = false;
      m = SmilesToMol(smiles, params);
      TEST_ASSERT(m);
      TEST_ASSERT(m->getNumAtoms() == 7);
      TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
      delete m;
    }
    {  // removeHs, no sanitization
      SmilesParserParams params;
      params.sanitize = false;
      m = SmilesToMol(smiles, params);
      TEST_ASSERT(m);
      TEST_ASSERT(m->getNumAtoms() == 6);
      TEST_ASSERT(!m->getBondWithIdx(0)->getIsAromatic());
      delete m;
    }
    {  // no removeHs, no sanitization
      SmilesParserParams params;
      params.removeHs = false;
      params.sanitize = false;
      m = SmilesToMol(smiles, params);
      TEST_ASSERT(m);
      TEST_ASSERT(m->getNumAtoms() == 7);
      TEST_ASSERT(!m->getBondWithIdx(0)->getIsAromatic());
      delete m;
    }
  }

  {  // basic name parsing
    std::string smiles = "CCCC the_name";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    delete m;
    {  // it's parsed:
      SmilesParserParams params;
      params.allowCXSMILES = false;
      m = SmilesToMol(smiles, params);
      TEST_ASSERT(m);
      TEST_ASSERT(m->hasProp(common_properties::_Name));
      TEST_ASSERT(m->getProp<std::string>(common_properties::_Name) ==
                  "the_name");
      delete m;
    }
    {
      SmilesParserParams params;
      params.strictCXSMILES = false;
      params.parseName = false;
      m = SmilesToMol(smiles, params);
      TEST_ASSERT(m);
      TEST_ASSERT(m->getNumAtoms() == 4);
      TEST_ASSERT(!m->hasProp(common_properties::_Name));
      delete m;
    }
  }
  {  // name parsing2
    std::string smiles = "CCCC\tthe_name";
    {  // no removeHs, no sanitization
      SmilesParserParams params;
      params.parseName = true;
      RWMol *m = SmilesToMol(smiles, params);
      TEST_ASSERT(m);
      TEST_ASSERT(m->getNumAtoms() == 4);
      TEST_ASSERT(m->hasProp(common_properties::_Name));
      TEST_ASSERT(m->getProp<std::string>(common_properties::_Name) ==
                  "the_name");
      delete m;
    }
  }
  {  // name parsing3
    std::string smiles = "CCCC\t  the_name  ";
    {  // no removeHs, no sanitization
      SmilesParserParams params;
      params.parseName = true;
      RWMol *m = SmilesToMol(smiles, params);
      TEST_ASSERT(m);
      TEST_ASSERT(m->getNumAtoms() == 4);
      TEST_ASSERT(m->hasProp(common_properties::_Name));
      TEST_ASSERT(m->getProp<std::string>(common_properties::_Name) ==
                  "the_name");
      delete m;
    }
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testRingClosureNumberWithBrackets() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
  BOOST_LOG(rdInfoLog)
      << "Testing the %(....) notation for SMILES ring closure numbers\n"
      << std::endl;
  {
    const char *benzenes[6] = {
        "c1ccccc1",           "c%(1)ccccc%(1)",       "c%(12)ccccc%(12)",
        "c%(123)ccccc%(123)", "c%(1234)ccccc%(1234)", "c%(99999)ccccc%(99999)"};
    for (auto &i : benzenes) {
      BOOST_LOG(rdInfoLog) << "Test: " << i << " (should be read)" << std::endl;
      ROMol *m = SmilesToMol(i);
      TEST_ASSERT(m);
      TEST_ASSERT(m->getNumAtoms() == 6);
      TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
      std::string benzene = MolToSmiles(*m, false, false, -1, false);
      TEST_ASSERT(benzene == "c1ccccc1");
      delete m;
    }

    const char *not_allowed[2] = {"c%()ccccc%()", "c%(100000)ccccc%(100000)"};
    for (auto &i : not_allowed) {
      BOOST_LOG(rdInfoLog) << "Test: " << i << " (should NOT be read)"
                           << std::endl;
      ROMol *m = SmilesToMol(i);
      TEST_ASSERT(m == (ROMol *)nullptr);
      delete m;
    }
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIsomericSmilesIsDefault() {
  BOOST_LOG(rdInfoLog)
      << "Testing that isomeric SMILES is now the default output" << std::endl;
  {
    std::string smi = "C[C@H](Cl)Br";
    auto m = SmilesToMol(smi);
    TEST_ASSERT(m)
    auto csmi = MolToSmiles(*m);
    TEST_ASSERT(csmi.find("@") != std::string::npos);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testHashAtomExtension() {
  BOOST_LOG(rdInfoLog) << "Testing constructs like [#6]" << std::endl;
  {
    std::string smi = "[#6][12#6]";
    auto m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getAtomicNum() == 6);
    TEST_ASSERT(m->getAtomWithIdx(0)->getIsotope() == 0);
    TEST_ASSERT(m->getAtomWithIdx(1)->getAtomicNum() == 6);
    TEST_ASSERT(m->getAtomWithIdx(1)->getIsotope() == 12);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub1925() {
  BOOST_LOG(rdInfoLog) << "Testing Github #1925: Atom with bond to itself is "
                          "accepted by the SMILES parser."
                       << std::endl;
  {
    std::string smi = "C1CC111";
    RWMol *m = nullptr;
    m = SmilesToMol(smi);
    TEST_ASSERT(!m);
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testdoRandomSmileGeneration() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Testing the random Generator for SMILES"
                       << std::endl;
  {
    // it's not trivial to test this because we're using std::rand(), which does
    // not give consistent results across platforms. It's not worth adding the
    // complexity of a real RNG, so we do some hand waving in the tests
    std::srand(0xf00d);  // be sure we use it for testcase!
    const std::vector<std::string> benzenes = {"COc1ccnc(CC)c1C"};

    const std::vector<std::string> rulesmiles = {
        "COc1ccnc(CC)c1C",     "O(C)c1ccnc(CC)c1C",   "c1(OC)ccnc(CC)c1C",
        "c1c(OC)c(C)c(CC)nc1", "c1cc(OC)c(C)c(CC)n1", "n1ccc(OC)c(C)c1CC",
        "c1(CC)nccc(OC)c1C",   "C(c1nccc(OC)c1C)C",   "CCc1nccc(OC)c1C",
        "c1(C)c(OC)ccnc1CC",   "Cc1c(OC)ccnc1CC"};

    for (auto bz : benzenes) {
      ROMol *m = SmilesToMol(bz);
      TEST_ASSERT(m);
      TEST_ASSERT(m->getNumAtoms() == 11);
      for (unsigned int j = 0; j < m->getNumAtoms(); ++j) {
        auto rulebenzene =
            MolToSmiles(*m, true, false, j, false, false, false, false);
        // BOOST_LOG(rdInfoLog) << "rule :" << rulebenzene << std::endl;
        // std::cout << "\"" << rulebenzene << "\", ";
        TEST_ASSERT(rulebenzene == rulesmiles[j]);
        std::set<std::string> rsmis;
        for (unsigned int iter = 0; iter < 10; ++iter) {
          auto randombenzene =
              MolToSmiles(*m, true, false, j, false, false, false, true);
          // BOOST_LOG(rdInfoLog) << "random :" << j << " " << iter << " "
          //                      << randombenzene << std::endl;
          rsmis.insert(randombenzene);
        }
        // we will get dupes, but there's enough choice available here that we
        // should have gotten at least 3 unique
        TEST_ASSERT(rsmis.size() >= 3);
      }
      // std::cout << std::endl;

      // confirm that we also use random starting points:
      std::set<char> starts;
      for (unsigned int iter = 0; iter < 50; ++iter) {
        auto randombenzene =
            MolToSmiles(*m, true, false, -1, false, false, false, true);
        // BOOST_LOG(rdInfoLog) << "random :" << j << " " << iter << " "
        //                      << randombenzene << std::endl;
        starts.insert(randombenzene[0]);
      }
      // we will get dupes, but there's enough choice available here that we
      // should have gotten at least 3 unique
      TEST_ASSERT(starts.find('C') != starts.end());
      TEST_ASSERT(starts.find('c') != starts.end());
      TEST_ASSERT(starts.find('n') != starts.end() ||
                  starts.find('O') != starts.end());

      delete m;
    }
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub1972() {
  BOOST_LOG(rdInfoLog)
      << "Testing Github #1972: Incorrect tetrahedral stereo when reading "
         "SMILES with ring closure as last neighbor"
      << std::endl;
  {
    std::vector<std::vector<std::string>> smiles = {
        {"[C@@]1(Cl)(F)(I).Br1", "[C@@](Br)(Cl)(F)(I)"},
        {"[C@@](Cl)(F)(I)1.Br1", "[C@@](Cl)(F)(I)Br"},
        {"[C@@](Cl)1(F)(I).Br1", "[C@@](Cl)(Br)(F)(I)"},
        {"[C@@](Cl)(F)1(I).Br1", "[C@@](Cl)(F)(Br)(I)"}};
    for (const auto &pr : smiles) {
      // std::cerr << "--------------------------" << std::endl;
      // std::cerr << pr[0] << " " << pr[1] << std::endl;
      std::unique_ptr<ROMol> m1(SmilesToMol(pr[0]));
      // std::cerr << "------------" << std::endl;
      std::unique_ptr<ROMol> m2(SmilesToMol(pr[1]));
      TEST_ASSERT(m1);
      TEST_ASSERT(m2);
      // m1->debugMol(std::cerr);
      // std::cerr << "------------" << std::endl;
      // m2->debugMol(std::cerr);
      auto csmi1 = MolToSmiles(*m1);
      auto csmi2 = MolToSmiles(*m2);
      // std::cerr << ">>> " << (csmi1 == csmi2) << " " << csmi1 << " " << csmi2
      //           << std::endl;
      TEST_ASSERT(csmi1 == csmi2);
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
      // std::cerr << "--------------------------" << std::endl;
      // std::cerr << pr[0] << " " << pr[1] << std::endl;
      std::unique_ptr<ROMol> m1(SmilesToMol(pr[0]));
      // std::cerr << "------------" << std::endl;
      std::unique_ptr<ROMol> m2(SmilesToMol(pr[1]));
      TEST_ASSERT(m1);
      TEST_ASSERT(m2);
      // m1->debugMol(std::cerr);
      // std::cerr << "------------" << std::endl;
      // m2->debugMol(std::cerr);
      auto csmi1 = MolToSmiles(*m1);
      auto csmi2 = MolToSmiles(*m2);
      // std::cerr << ">>> " << (csmi1 == csmi2) << " " << csmi1 << " " << csmi2
      //           << std::endl;
      TEST_ASSERT(csmi1 == csmi2);
    }
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub2556() {
  BOOST_LOG(rdInfoLog) << "Testing Github #2556: Test correct parsing and fix "
                          "memory leak for C1C1"
                       << std::endl;
  RWMol *m = nullptr;
  m = SmilesToMol("C1C1");
  TEST_ASSERT(!m);
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub1028() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github issue #1028: Alternating canonical "
                          "SMILES for ring with chiral N"
                       << std::endl;
  // note that due to the changes made for #3631, the N's originally used in
  // these tests are no longer considered to be chiral. I switched to using P
  // (and verified that P was also a problem before #1028 was fixed)
  {
    std::string smi = "O[C@H]1CC2CCC(C1)[P@@]2C";
    const std::string ref = "C[P@]1C2CCC1C[C@H](O)C2";
    for (int i = 0; i < 3; ++i) {
      const auto mol = std::unique_ptr<ROMol>(SmilesToMol(smi));
      TEST_ASSERT(mol);
      const std::string out = MolToSmiles(*mol);
      TEST_ASSERT(out == ref);
      smi = out;
    }

    {
      std::string smi = "C[P@]1C[C@@H](O)C1";
      const std::string ref = smi;
      for (int i = 0; i < 3; ++i) {
        const auto mol = std::unique_ptr<ROMol>(SmilesToMol(smi));
        TEST_ASSERT(mol);
        const std::string out = MolToSmiles(*mol);
        TEST_ASSERT(out == ref);
        smi = out;
      }
    }
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub3139() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github issue #3139: Partial bond mem leak"
                       << std::endl;

  {
    const std::string smi = "COc(c1)cccc1C#";
    for (int i = 0; i < 3; ++i) {
      const auto mol = std::unique_ptr<ROMol>(SmilesToMol(smi));
      const auto sma = std::unique_ptr<ROMol>(SmartsToMol(smi));
    }
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testOSSFuzzFailures() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Failures/problems detected by OSS Fuzz" << std::endl;

  {  // examples that should produce no molecule
    std::vector<std::string> failing_examples = {"C)"};
    for (auto smi : failing_examples) {
      const auto mol = std::unique_ptr<ROMol>(SmilesToMol(smi));
      // output which molecule is failing
      if (mol) {
        std::cerr << "  Should have failed: " << smi << std::endl;
        TEST_ASSERT(!mol);
      }
    }
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub3967() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github Issue 3967: Double bond stereo gets "
                          "flipped by SMILES reader/writer"
                       << std::endl;

  {
    auto mol = "C=c1s/c2n(c1=O)CCCCCCC\\N=2"_smiles;
    TEST_ASSERT(mol);
    auto smi = MolToSmiles(*mol);
    std::cerr << smi << std::endl;
    TEST_ASSERT(smi == "C=c1s/c2n(c1=O)CCCCCCC\\N=2");
  }
  {
    auto mol = R"SMI(C1=C\C/C=C2C3=C/C/C=C\C=C/C\3C\2\C=C/1)SMI"_smiles;
    TEST_ASSERT(mol);
    auto smi = MolToSmiles(*mol);
    std::cerr << smi << std::endl;
    TEST_ASSERT(smi == R"SMI(C1=C\C/C=C2C3=C\C/C=C\C=C/C/3C\2\C=C/1)SMI");
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub6349() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing Github Issue 6349: Different SMARTS input formats lead to different SMILES outputs."
      << std::endl;

  auto checkSmartsToSmiles = [](const std::string &sma,
                                const std::string &refSmi) {
    std::unique_ptr<ROMol> molFromSmarts(SmartsToMol(sma));
    {
      std::string smi = MolToSmiles(*molFromSmarts);
      TEST_ASSERT(smi == refSmi);
    }

    std::string molBlock = MolToMolBlock(*molFromSmarts);
    std::unique_ptr<ROMol> molFromBlock(
        MolBlockToMol(molBlock, /*sanitize =*/false, /*removeHs =*/false));
    {
      std::string smi = MolToSmiles(*molFromBlock);
      TEST_ASSERT(smi == refSmi);
    }
  };
  checkSmartsToSmiles("[C]", "C");
  checkSmartsToSmiles("[C,N]", "*");
  checkSmartsToSmiles("[C,N]~[O,S]", "*~*");
  checkSmartsToSmiles("C-C(-[Cl,F,Br])-C", "*C(C)C");

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();
// boost::logging::enable_logs("rdApp.debug");
#if 1
  testPass();
  testFail();

  testDetails();
  testProblems();
  // testLeak();
  testBasicCanon();
  testIssue127();
  testIssue143();
  testIssue151();
  testIssue153();
  testIssue175();
  testIssue176();
  testIssue180();
  testIssue159();
  testIssue184();
  testIssue185();
  testIssue191();
  testIssue256();
  testIssue266();
  testRootedAt();
  testIsotopes();
  testBug1670149();
  testBug1842174();
  testBug1844959();
  testIssue157();
  testStereochem();
  testBug1942220();
  testBug3127883();
  testAtomMaps();
  testBug3145697();
  testBug3152751();
  testReplacementPatterns();
  testAllBondsExplicit();
  testBug3139534();
  testBug3526815();
  testBug3525799();
  testBug3526810();
  testBug3528556();
  testBug253();
  testBug257();
  testRingStereochem();
  testGithub45();
  testGithub206();
  testGithub210();
  testGithub378();
  testGithub389();
  testBug1719046();
  testBug1844617();

#if 1  // POSTPONED during canonicalization rewrite
  // testGithub298();
  testFragmentSmiles();
  testGithub12();
#endif
  testSmilesWriteForModifiedMolecules();
  testGithub532();
  testGithub786();
  testGithub760();
  testDativeBonds();
  testGithub1219();
  testSmilesParseParams();
  testRingClosureNumberWithBrackets();
  testGithub1652();
  testIsomericSmilesIsDefault();
  testHashAtomExtension();
  testGithub1925();
  testGithub1972();
  testGithub2556();
  testdoRandomSmileGeneration();
  testGithub1028();
  testGithub3139();
  testGithub3967();
  testGithub6349();
#endif
  testOSSFuzzFailures();
}
