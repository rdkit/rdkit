//
//  Copyright (C) 2006-2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <string>
#include <iostream>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/Exceptions.h>

using namespace RDKit;

void testDeleteSubstruct() {
  ROMol *mol1 = nullptr, *mol2 = nullptr, *matcher1 = nullptr,
        *matcher2 = nullptr, *matcher3 = nullptr;
  std::string smi, sma;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing deleteSubstruct" << std::endl;

  // a lot of the seemingly repetitive stuff is here for Issue96
  smi = "CCC(=O).C=O";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  TEST_ASSERT(mol1->getNumAtoms() == 6)
  sma = "C=O";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  mol2 = deleteSubstructs(*mol1, *matcher1, 0);
  delete mol1;
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 2)
  mol1 = deleteSubstructs(*mol2, *matcher1, 0);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 2)
  delete matcher1;
  delete mol1;
  delete mol2;

  sma = "[Cl;H1&X1,-]";
  matcher1 = SmartsToMol(sma);
  sma = "[Na+]";
  matcher2 = SmartsToMol(sma);
  sma = "[O;H2,H1&-,X0&-2]";
  matcher3 = SmartsToMol(sma);
  mol1 = SmilesToMol("CCO.Cl");
  TEST_ASSERT(mol1);
  TEST_ASSERT(mol1->getNumAtoms() == 4);
  mol2 = deleteSubstructs(*mol1, *matcher1, true);
  delete mol1;
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 3);
  mol1 = deleteSubstructs(*mol2, *matcher2, true);
  delete mol1;
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 3);
  mol1 = deleteSubstructs(*mol2, *matcher3, true);
  delete mol1;
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 3);
  delete mol2;
  delete matcher3;

  mol1 = SmilesToMol("CC(=O)[O-].[Na+]");
  TEST_ASSERT(mol1);
  TEST_ASSERT(mol1->getNumAtoms() == 5);

  delete matcher1;
  matcher1 = SmartsToMol("[Cl;H1&X1,-]");
  delete matcher2;
  matcher2 = SmartsToMol("[Na+]");

  mol2 = deleteSubstructs(*mol1, *matcher1, true);
  delete mol1;
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 5);
  mol1 = deleteSubstructs(*mol2, *matcher2, true);
  delete mol2;
  TEST_ASSERT(mol1);
  TEST_ASSERT(mol1->getNumAtoms() == 4);
  mol2 = deleteSubstructs(*mol1, *matcher1, true);
  TEST_ASSERT(mol1);
  TEST_ASSERT(mol1->getNumAtoms() == 4);
  delete mol2;
  delete mol1;
  delete matcher1;
  delete matcher2;

  // test chirality
  smi = "CCO[C@H](N)(P)";
  mol1 = SmilesToMol(smi);
  matcher1 = SmartsToMol("O[C@H](N)(P)");
  mol2 = deleteSubstructs(*mol1, *matcher1, false, true);
  std::string smi1 = MolToSmiles(*mol2, true);
  std::cerr << "1 smi: " << smi1 << std::endl;
  TEST_ASSERT(smi1 == "CC");
  delete mol2;

  // still matches with non-chiral
  mol2 = deleteSubstructs(*mol1, *matcher1, false, false);
  smi1 = MolToSmiles(*mol2, true);
  std::cerr << "1 smi: " << smi1 << std::endl;
  TEST_ASSERT(smi1 == "CC");
  delete mol2;
  delete matcher1;

  matcher1 = SmartsToMol("O[C@@H](N)(P)");
  mol2 = deleteSubstructs(*mol1, *matcher1, false, true);
  smi1 = MolToSmiles(*mol2, true);
  std::cerr << "2 smi: " << smi1 << std::endl;
  TEST_ASSERT(smi1 == "CCO[C@H](N)P");
  delete mol2;

  // still matches with non-chiral
  mol2 = deleteSubstructs(*mol1, *matcher1, false, false);
  smi1 = MolToSmiles(*mol2, true);
  std::cerr << "2 smi: " << smi1 << std::endl;
  TEST_ASSERT(smi1 == "CC");
  delete mol2;
  delete mol1;
  delete matcher1;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testReplaceSubstructs() {
  ROMol *mol1 = nullptr, *matcher1 = nullptr, *frag = nullptr;
  std::string smi, sma;
  std::vector<ROMOL_SPTR> vect;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing replaceSubstruct" << std::endl;

  smi = "CCCC";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);

  sma = "C(=O)O";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  smi = "N";
  frag = SmilesToMol(smi);
  TEST_ASSERT(frag);

  vect = replaceSubstructs(*mol1, *matcher1, *frag);
  TEST_ASSERT(vect.size() == 1);
  TEST_ASSERT(mol1->getNumAtoms() == 4);
  TEST_ASSERT(vect[0]->getNumAtoms() == 4);

  delete mol1;
  smi = "CCCC(=O)O";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  vect = replaceSubstructs(*mol1, *matcher1, *frag);
  TEST_ASSERT(vect.size() == 1);
  TEST_ASSERT(mol1->getNumAtoms() == 6);
  TEST_ASSERT(vect[0]->getNumAtoms() == 4);

  vect = replaceSubstructs(*mol1, *matcher1, *frag, true);
  TEST_ASSERT(vect.size() == 1);
  TEST_ASSERT(mol1->getNumAtoms() == 6);
  TEST_ASSERT(vect[0]->getNumAtoms() == 4);

  delete mol1;
  smi = "OC(=O)CCCC(=O)O";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  vect = replaceSubstructs(*mol1, *matcher1, *frag);
  TEST_ASSERT(vect.size() == 2);
  TEST_ASSERT(mol1->getNumAtoms() == 9);
  TEST_ASSERT(vect[0]->getNumAtoms() == 7);
  TEST_ASSERT(vect[1]->getNumAtoms() == 7);

  // use replaceAll when a single replacement is available:
  vect = replaceSubstructs(*mol1, *matcher1, *frag, true);
  TEST_ASSERT(vect.size() == 1);
  TEST_ASSERT(mol1->getNumAtoms() == 9);
  TEST_ASSERT(vect[0]->getNumAtoms() == 5);

  delete matcher1;
  sma = "C=O";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  delete mol1;
  smi = "CC(=O)C";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  vect = replaceSubstructs(*mol1, *matcher1, *frag);
  TEST_ASSERT(vect.size() == 1);
  TEST_ASSERT(mol1->getNumAtoms() == 4);
  TEST_ASSERT(vect[0]->getNumAtoms() == 3);

  delete mol1;
  smi = "C1C(=O)C1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  vect = replaceSubstructs(*mol1, *matcher1, *frag);
  TEST_ASSERT(vect.size() == 1);
  TEST_ASSERT(mol1->getNumAtoms() == 4);
  TEST_ASSERT(vect[0]->getNumAtoms() == 3);
  TEST_ASSERT(vect[0]->getNumBonds() == 3);

  // test chirality
  delete mol1;
  delete matcher1;
  sma = "P([C@H](N)O)";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  // test stereo matching
  smi = "CCP([C@H](N)O)CC";
  mol1 = SmilesToMol(smi);
  vect = replaceSubstructs(*mol1, *matcher1, *frag, false, 0, true);
  TEST_ASSERT(vect.size() == 1);
  std::string smi1 = MolToSmiles(*vect[0], true);
  TEST_ASSERT(smi1 == "CCNCC");

  // should also match when matching non-chiral
  vect = replaceSubstructs(*mol1, *matcher1, *frag, false, 0, false);
  TEST_ASSERT(vect.size() == 1);
  smi1 = MolToSmiles(*vect[0], true);
  TEST_ASSERT(smi1 == "CCNCC");
  delete mol1;

  smi = "CCP([C@@H](N)O)CC";
  mol1 = SmilesToMol(smi);
  vect = replaceSubstructs(*mol1, *matcher1, *frag, false, 0, true);
  TEST_ASSERT(vect.size() == 1);
  smi1 = MolToSmiles(*vect[0], true);
  std::cerr << "replaceSub smi1:" << smi1 << std::endl;
  ;
  // no change
  TEST_ASSERT(smi1 == "CCP(CC)[C@@H](N)O");

  // should also match when matching non-chiral
  vect = replaceSubstructs(*mol1, *matcher1, *frag, false, 0, false);
  TEST_ASSERT(vect.size() == 1);
  smi1 = MolToSmiles(*vect[0], true);
  std::cerr << "replaceSub smi1:" << smi1 << std::endl;
  delete mol1;
  delete frag;
  delete matcher1;

  // no change
  TEST_ASSERT(smi1 == "CCNCC");

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testReplaceSubstructs2() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "More testing of replaceSubstruct" << std::endl;

  {
    std::string smi = "CCC(=O)O";
    ROMol *mol1 = SmilesToMol(smi);
    TEST_ASSERT(mol1);

    std::string sma = "[$(O-C(=O))]";
    ROMol *matcher1 = SmartsToMol(sma);
    TEST_ASSERT(matcher1);

    std::string repl = "NC";
    ROMol *frag = SmilesToMol(repl);
    TEST_ASSERT(frag);

    std::vector<ROMOL_SPTR> vect = replaceSubstructs(*mol1, *matcher1, *frag);
    TEST_ASSERT(vect.size() == 1);
    TEST_ASSERT(mol1->getNumAtoms() == 5);

    TEST_ASSERT(vect[0]->getNumAtoms() == 6);
    std::string csmi1 = MolToSmiles(*vect[0], true);

    repl = "CN";
    delete frag;
    frag = SmilesToMol(repl);
    TEST_ASSERT(frag);
    vect = replaceSubstructs(*mol1, *matcher1, *frag, true, 1);
    TEST_ASSERT(vect.size() == 1);
    TEST_ASSERT(mol1->getNumAtoms() == 5);
    TEST_ASSERT(vect[0]->getNumAtoms() == 6);
    std::string csmi2 = MolToSmiles(*vect[0], true);
    TEST_ASSERT(csmi2 == csmi1);
    delete mol1;
    delete frag;
    delete matcher1;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testReplaceSidechains() {
  ROMol *mol1 = nullptr, *mol2 = nullptr, *matcher1 = nullptr;
  std::string smi, sma;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing replaceSidechains" << std::endl;

  smi = "ClC1CC(F)C1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);

  sma = "C1CCC1";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  mol2 = replaceSidechains(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 6);
  delete mol2;
  delete mol1;

  smi = "ClC1C(F)C1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  mol2 = replaceSidechains(*mol1, *matcher1);
  TEST_ASSERT(!mol2);
  delete mol1;
  delete mol2;

  delete matcher1;
  sma = "C=O";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  smi = "CC=O";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  mol2 = replaceSidechains(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 3);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(smi == "[1*]C=O");
  delete mol1;
  delete mol2;

  smi = "CC(C)=O";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  mol2 = replaceSidechains(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 4);
  smi = MolToSmiles(*mol2, true);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi == "[1*]C([2*])=O" || smi == "[2*]C([1*])=O");

  delete mol1;
  delete mol2;
  delete matcher1;

  sma = "P([C@H](N)O)";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  smi = "CCP([C@H](N)O)CC";
  mol1 = SmilesToMol(smi);
  mol2 = replaceSidechains(*mol1, *matcher1, true);
  TEST_ASSERT(mol2);
  smi = MolToSmiles(*mol2, true);
  std::cerr << "sidechains chiral=true smi1;;; " << smi << std::endl;
  delete mol2;

  mol2 = replaceSidechains(*mol1, *matcher1, true);
  TEST_ASSERT(mol2);
  smi = MolToSmiles(*mol2, true);
  std::cerr << "sidechains chiral=false smi1;;; " << smi << std::endl;

  delete mol1;
  delete mol2;
  delete matcher1;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testReplaceCore() {
  ROMol *mol1 = nullptr, *mol2 = nullptr, *matcher1 = nullptr;
  std::string smi, sma;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing replaceCore" << std::endl;

  smi = "ClC1CC(F)C1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  sma = "C1CC1";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(!mol2);
  delete mol1;
  delete mol2;
  delete matcher1;

  smi = "ClC1CC(F)C1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  sma = "C1CCC1";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 4);
  smi = MolToSmiles(*mol2, true);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi == "[1*]Cl.[2*]F" || smi == "[2*]Cl.[1*]F");
  delete mol1;
  delete mol2;
  delete matcher1;

  smi = "CCC=O";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  sma = "C=O";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 3);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(smi == "[1*]CC");
  delete mol1;
  delete mol2;
  delete matcher1;

  smi = "C1C(=O)CC1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  sma = "C=O";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 5);
  smi = MolToSmiles(*mol2, true);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi == "[1*]CCC[2*]" || smi == "[1*]CCC[2*]");
  delete mol1;
  delete mol2;
  delete matcher1;

  smi = "CNC";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  sma = "N";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 4);
  smi = MolToSmiles(*mol2, true);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi == "[1*]C.[2*]C" || smi == "[2*]C.[1*]C");
  delete mol1;
  delete mol2;
  delete matcher1;

  smi = "OC1CCC1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  sma = "[CH2][CH2][CH2]";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(mol2->getNumAtoms() == 4);
  smi = MolToSmiles(*mol2, true);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi == "[2*]C([1*])O" || smi == "[1*]C([2*])O");
  delete mol1;
  delete mol2;
  delete matcher1;

  smi = "C/C=C/CN/C=C/O";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  sma = "N";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 9);
  smi = MolToSmiles(*mol2, true);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi == "[1*]C/C=C/C.[2*]/C=C/O" ||
              smi == "[1*]/C=C/O.[2*]C/C=C/C");
  delete mol1;
  delete mol2;
  delete matcher1;

  smi = "C[C@](F)(Cl)N";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  sma = "N";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 5);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(smi == "[1*][C@@](C)(F)Cl");
  delete mol1;
  delete mol2;

  smi = "C[C@](F)(N)Cl";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 5);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(smi == "[1*][C@](C)(F)Cl");
  delete mol1;
  delete mol2;

  smi = "N[C@](C)(F)Cl";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 5);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(smi == "[1*][C@](C)(F)Cl");
  delete mol1;
  delete mol2;
  delete matcher1;

  smi = "C[C@@](F)(Cl)N";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  sma = "N";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 5);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(smi == "[1*][C@](C)(F)Cl");
  delete mol1;
  delete mol2;

  smi = "C[C@@](F)(N)Cl";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 5);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(smi == "[1*][C@@](C)(F)Cl");
  smi = "C[C@@](N)(Cl)F";
  delete mol1;
  delete mol2;

  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 5);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(smi == "[1*][C@@](C)(F)Cl");
  delete mol1;
  delete mol2;
  delete matcher1;

  smi = "C[C@H](F)N";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  sma = "N";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 4);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(smi == "[1*][C@H](C)F");
  delete mol1;
  delete mol2;

  smi = "C[C@H](N)F";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 4);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(smi == "[1*][C@@H](C)F");
  delete mol1;
  delete mol2;

  smi = "N[C@H](C)F";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 4);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(smi == "[1*][C@H](C)F");
  delete mol1;
  delete mol2;

  smi = "F[C@H](C)N";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 4);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(smi == "[1*][C@@H](C)F");
  delete mol1;
  delete mol2;

  smi = "N[C@H]1CCCO1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 6);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(smi == "[1*][C@H]1CCCO1");
  delete mol1;
  delete mol2;
  delete matcher1;

  smi = "ClC1CC(F)C1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  sma = "*C1CC(*)C1";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  mol2 = replaceCore(*mol1, *matcher1, false);
  TEST_ASSERT(mol2);
  smi = MolToSmiles(*mol2, true);
  TEST_ASSERT(mol2->getNumAtoms() == 4);
  smi = MolToSmiles(*mol2, true);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi == "[1*]F.[2*]Cl" || smi == "[1*]Cl.[2*]F");

  delete mol1;
  delete mol2;
  delete matcher1;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

struct CoreTest {
  const char *smiles;
  const char *smarts;
  bool replaceDummies;
  bool labelByIndex;
  bool requireDummyMatch;
  bool useChirality;
  const char *expected;
};

void testReplaceCore2() {
  const CoreTest tests[] = {
      {"C1O[C@@]1(OC)NC", "C1O[C@]1(*)*", false, false, false, false,
       "[1*]OC.[2*]NC"},
      {"C1O[C@@]1(OC)NC", "C1O[C@]1(*)*", false, false, false, true,
       "[1*]NC.[2*]OC"},
      {"C1O[C@@]1(OC)NC", "C1O[C@]1(*)*", false, true, false, false,
       "[3*]OC.[4*]NC"},
      {"C1O[C@@]1(OC)NC", "C1O[C@]1(*)*", false, true, false, true,
       "[3*]NC.[4*]OC"},
      {"C1O[C@@]1(OC)NC", "C1O[C@]1(*)*", true, false, false, false,
       "[1*]C.[2*]C"},
      {"C1O[C@@]1(OC)NC", "C1O[C@]1(*)*", true, false, false, true,
       "[1*]C.[2*]C"},
      {"C1O[C@@]1(OC)NC", "C1O[C@]1(*)*", true, true, false, false,
       "[3*]C.[4*]C"},
      {"C1O[C@@]1(OC)NC", "C1O[C@]1(*)*", true, true, false, true,
       "[3*]C.[4*]C"},

      {"C1O[C@]1(OC)NC", "C1O[C@]1(*)*", false, false, false, false,
       "[1*]OC.[2*]NC"},
      {"C1O[C@]1(OC)NC", "C1O[C@]1(*)*", false, false, false, true,
       "[1*]OC.[2*]NC"},
      {"C1O[C@]1(OC)NC", "C1O[C@]1(*)*", false, true, false, false,
       "[3*]OC.[4*]NC"},
      {"C1O[C@]1(OC)NC", "C1O[C@]1(*)*", false, true, false, true,
       "[3*]OC.[4*]NC"},

      {"C1O[C@]1(OC)NC", "C1O[C@]1(*)*", true, false, false, false,
       "[1*]C.[2*]C"},
      {"C1O[C@]1(OC)NC", "C1O[C@]1(*)*", true, false, false, true,
       "[1*]C.[2*]C"},
      {"C1O[C@]1(OC)NC", "C1O[C@]1(*)*", true, true, false, false,
       "[3*]C.[4*]C"},
      {"C1O[C@]1(OC)NC", "C1O[C@]1(*)*", true, true, false, true,
       "[3*]C.[4*]C"},

      {"C1O[C@@]1(OC)NCC", "C1O[C@]1(*)*", true, false, false, false,
       "[1*]C.[2*]CC"},
      {"C1O[C@@]1(OC)NCC", "C1O[C@]1(*)*", true, false, false, true,
       "[1*]CC.[2*]C"},

      {"C1O[C@@]1(OC)NCC", "C1O[C@]1(*)*", true, true, false, false,
       "[3*]C.[4*]CC"},
      {"C1O[C@@]1(OC)NCC", "C1O[C@]1(*)*", true, true, false, true,
       "[3*]CC.[4*]C"},

      {"C1O[C@]1(OC)NCC", "C1O[C@]1(*)*", true, false, false, false,
       "[1*]C.[2*]CC"},
      {"C1O[C@]1(OC)NCC", "C1O[C@]1(*)*", true, false, false, true,
       "[1*]C.[2*]CC"},
      {"C1O[C@]1(OC)NCC", "C1O[C@]1(*)*", true, true, false, false,
       "[3*]C.[4*]CC"},
      {"C1O[C@]1(OC)NCC", "C1O[C@]1(*)*", true, true, false, true,
       "[3*]C.[4*]CC"},

      {"CNOC", "CONC", false, true, false, false, ""},
      {"PCNOCS", "CONC", false, true, false, false, "*S.[3*]P"},
      {"PCNOCS", "CONC", false, false, false, false, "[1*]S.[2*]P"},

      {"PCONCS", "CONC", false, true, false, false, "*P.[3*]S"},
      {"PCONCS", "CONC", false, false, false, false, "[1*]P.[2*]S"}

  };
  size_t num_tests = sizeof(tests) / sizeof(CoreTest);
  for (size_t i = 0; i < num_tests; ++i) {
    ROMOL_SPTR mol(SmilesToMol(tests[i].smiles));
    ROMOL_SPTR smarts(SmartsToMol(tests[i].smarts));
    ROMOL_SPTR res(replaceCore(*mol.get(), *smarts.get(),
                               tests[i].replaceDummies, tests[i].labelByIndex,
                               tests[i].requireDummyMatch,
                               tests[i].useChirality));
    if (tests[i].expected) {
      TEST_ASSERT(res.get());
      std::string smi = MolToSmiles(*res.get(), true);
      if (smi != tests[i].expected) {
        std::cerr << i << " " << tests[i].smiles << " " << tests[i].smarts
                  << " " << (int)tests[i].replaceDummies << " "
                  << (int)tests[i].labelByIndex << " "
                  << (int)tests[i].requireDummyMatch << " "
                  << (int)tests[i].useChirality
                  << " expected:" << tests[i].expected << " got => " << smi
                  << std::endl;
      }
      TEST_ASSERT(smi == tests[i].expected);
    } else {
      TEST_ASSERT(!res.get());
    }

    MatchVectType matchV;

    // do the substructure matching and get the atoms that match the query
    const bool recursionPossible = true;
    bool matchFound = SubstructMatch(*mol.get(), *smarts.get(), matchV,
                                     recursionPossible, tests[i].useChirality);
    TEST_ASSERT(matchFound);
    res = ROMOL_SPTR(replaceCore(*mol.get(), *smarts.get(), matchV,
                                 tests[i].replaceDummies, tests[i].labelByIndex,
                                 tests[i].requireDummyMatch));

    if (tests[i].expected) {
      TEST_ASSERT(res.get());
      std::string smi = MolToSmiles(*res.get(), true);
      if (smi != tests[i].expected) {
        std::cerr << i << " " << tests[i].smiles << " " << tests[i].smarts
                  << " " << (int)tests[i].replaceDummies << " "
                  << (int)tests[i].labelByIndex << " "
                  << (int)tests[i].requireDummyMatch << " "
                  << (int)tests[i].useChirality
                  << " expected:" << tests[i].expected << " got => " << smi
                  << std::endl;
      }
      TEST_ASSERT(smi == tests[i].expected);
    } else {
      TEST_ASSERT(!res.get());
    }
  }
}

void testReplaceCoreLabels() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing replaceCore with labels" << std::endl;

  {
    std::string sma = "n1cocc1";
    ROMol *matcher = SmartsToMol(sma);
    TEST_ASSERT(matcher);

    std::string smi = "n1c(CC)oc(C)c1";
    ROMol *mol1 = SmilesToMol(smi);
    TEST_ASSERT(mol1);

    ROMol *mol2 = replaceCore(*mol1, *matcher, true, true);
    TEST_ASSERT(mol2);
    TEST_ASSERT(mol2->getNumAtoms() == 5);
    smi = MolToSmiles(*mol2, true);
    TEST_ASSERT(smi == "[1*]CC.[3*]C");

    delete mol1;
    delete mol2;
    delete matcher;
  }

  {
    std::string sma = "n1cocc1";
    ROMol *matcher = SmartsToMol(sma);
    TEST_ASSERT(matcher);

    std::string smi = "n1c(C)oc(CC)c1";
    ROMol *mol1 = SmilesToMol(smi);
    TEST_ASSERT(mol1);

    ROMol *mol2 = replaceCore(*mol1, *matcher, true, true);
    TEST_ASSERT(mol2);
    smi = MolToSmiles(*mol2, true);
    TEST_ASSERT(smi == "[1*]C.[3*]CC");

    delete mol1;
    delete mol2;
    delete matcher;
  }

  {
    std::string sma = "n1cocc1";
    ROMol *matcher = SmartsToMol(sma);
    TEST_ASSERT(matcher);

    std::string smi = "CCC1=C(OC)N=C(C)O1";
    ROMol *mol1 = SmilesToMol(smi);
    TEST_ASSERT(mol1);

    ROMol *mol2 = replaceCore(*mol1, *matcher, true, true);
    TEST_ASSERT(mol2);
    smi = MolToSmiles(*mol2, true);
    TEST_ASSERT(smi == "[1*]C.[3*]CC.[4*]OC");

    delete mol1;
    delete mol2;
    delete matcher;
  }

  {
    std::string sma = "n1cocc1";
    ROMol *matcher = SmartsToMol(sma);
    TEST_ASSERT(matcher);

    std::string smi = "C1=C(OC)N=CO1";
    ROMol *mol1 = SmilesToMol(smi);
    TEST_ASSERT(mol1);

    ROMol *mol2 = replaceCore(*mol1, *matcher, true, true);
    TEST_ASSERT(mol2);
    smi = MolToSmiles(*mol2, true);
    TEST_ASSERT(smi == "[4*]OC");

    delete mol1;
    delete mol2;
    delete matcher;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testReplaceCoreCrash() {
  ROMol *mol1 = nullptr, *mol2 = nullptr, *matcher1 = nullptr;
  std::string smi, sma;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing a former crash in replaceCore" << std::endl;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ChemTransforms/testData/oldcrash.mol";
  mol1 = MolFileToMol(pathName);
  TEST_ASSERT(mol1);

  sma = "N";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 6);
  smi = MolToSmiles(*mol2, true);
  delete mol1;
  delete mol2;
  delete matcher1;

  // there's no way to guarantee the order here:
  TEST_ASSERT(smi == "[1*]CC.[2*]CC" || smi == "[2*]CC.[1*]CC");
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testReplaceCorePositions() {
  ROMol *mol1 = nullptr, *mol2 = nullptr, *matcher1 = nullptr;
  std::string smi, sma;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing that atom positions are correctly copied by replaceCore"
      << std::endl;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ChemTransforms/testData/core_location.mol";
  mol1 = MolFileToMol(pathName);
  TEST_ASSERT(mol1);
  TEST_ASSERT(mol1->getNumAtoms() == 8);

  sma = "c1ccccc1";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  mol2 = replaceCore(*mol1, *matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms() == 4);

  RDGeom::Point3D op, np;
  op = mol1->getConformer().getAtomPos(6);
  np = mol2->getConformer().getAtomPos(0);
  TEST_ASSERT(feq(op.x, np.x));
  TEST_ASSERT(feq(op.y, np.y));
  TEST_ASSERT(feq(op.z, np.z));
  op = mol1->getConformer().getAtomPos(7);
  np = mol2->getConformer().getAtomPos(1);
  TEST_ASSERT(feq(op.x, np.x));
  TEST_ASSERT(feq(op.y, np.y));
  TEST_ASSERT(feq(op.z, np.z));
  op = mol1->getConformer().getAtomPos(0);
  np = mol2->getConformer().getAtomPos(2);
  TEST_ASSERT(feq(op.x, np.x));
  TEST_ASSERT(feq(op.y, np.y));
  TEST_ASSERT(feq(op.z, np.z));
  op = mol1->getConformer().getAtomPos(4);
  np = mol2->getConformer().getAtomPos(3);
  TEST_ASSERT(feq(op.x, np.x));
  TEST_ASSERT(feq(op.y, np.y));
  TEST_ASSERT(feq(op.z, np.z));
  delete mol1;
  delete mol2;
  delete matcher1;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testReplaceCoreMatchVectMultipleMappingToCore() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing replaceCore with a matchvect with multiple "
                          "dummy target atoms mapping to single core atom"
                       << std::endl;

  std::unique_ptr<RWMol> core, mol;

  // 2 dummies in core map to single N in molecule
  core = "C1([*:1])C([*:2])CC1"_smiles;
  mol = "C1CC2NC12"_smiles;
  MatchVectType mapping{{0, 4}, {1, 3}, {2, 2}, {3, 3}, {4, 1}, {5, 0}};

  ROMOL_SPTR result(replaceCore(*mol.get(), *core.get(), mapping, false));
  std::string smi = MolToSmiles(*result.get(), true);
  TEST_ASSERT("[1*]N[2*]" == smi);

  // 3 dummies in core map to single N in molecule
  core = "C1([*:1])C([*:2])C([*:3])C1"_smiles;
  mol = "C1C2C3N2C13"_smiles;
  mapping.clear();
  mapping.insert(mapping.end(),
                 {{0, 4}, {1, 3}, {2, 2}, {3, 3}, {4, 1}, {5, 3}, {6, 0}});
  result.reset(replaceCore(*mol.get(), *core.get(), mapping, false));
  smi = MolToSmiles(*result.get(), true);
  TEST_ASSERT("[1*]N([2*])[3*]" == smi);
}

void testReplaceCoreMatchVect() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing replaceCore with a matchvect" << std::endl;

  std::map<std::string, int> expected;
  expected["[1*]N.[2*]O.[3*]P"] = 1;
  expected["[1*]N.[2*]P.[3*]O"] = 1;
  expected["[1*]O.[2*]N.[3*]P"] = 1;
  expected["[1*]O.[2*]P.[3*]N"] = 1;
  expected["[1*]P.[2*]N.[3*]O"] = 1;
  expected["[1*]P.[2*]O.[3*]N"] = 1;

  const char *smiles = "NC1C(O)C1P";
  const char *smarts = "*C1C(*)C1*";
  ROMOL_SPTR mol(SmilesToMol(smiles));
  ROMOL_SPTR query(SmartsToMol(smarts));

  std::vector<MatchVectType> matches;
  const bool uniquify = false;
  unsigned int matchCount =
      SubstructMatch(*mol.get(), *query.get(), matches, uniquify);
  TEST_ASSERT(matchCount);

  std::map<std::string, int> results;
  const bool replaceDummies = false;
  for (unsigned int i = 0; i < matchCount; ++i) {
    ROMOL_SPTR res(
        replaceCore(*mol.get(), *query.get(), matches[i], replaceDummies));
    std::string smi = MolToSmiles(*res.get(), true);
    results[smi] = results[smi] + 1;
  }

  TEST_ASSERT(expected == results);

  {
    MatchVectType fake;
    fake.push_back(std::make_pair(1, 100));
    try {
      ROMOL_SPTR res(
          replaceCore(*mol.get(), *query.get(), fake, replaceDummies));
      TEST_ASSERT(0);
    } catch (...) {
      BOOST_LOG(rdInfoLog) << "-- caught error";
    }
  }

  {
    MatchVectType fake;
    fake.push_back(std::make_pair(100, 1));
    try {
      ROMOL_SPTR res(
          replaceCore(*mol.get(), *query.get(), fake, replaceDummies));
      TEST_ASSERT(0);
    } catch (...) {
      BOOST_LOG(rdInfoLog) << "-- caught error";
    }
  }

  {
    MatchVectType fake;
    fake.push_back(std::make_pair(1, -10));
    try {
      ROMOL_SPTR res(
          replaceCore(*mol.get(), *query.get(), fake, replaceDummies));
      TEST_ASSERT(0);
    } catch (...) {
      BOOST_LOG(rdInfoLog) << "-- caught error";
    }
  }

  {
    MatchVectType fake;
    fake.push_back(std::make_pair(-10, 1));
    try {
      ROMOL_SPTR res(
          replaceCore(*mol.get(), *query.get(), fake, replaceDummies));
      TEST_ASSERT(0);
    } catch (...) {
      BOOST_LOG(rdInfoLog) << "-- caught error";
    }
  }
}

void testMurckoDecomp() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing murcko decomposition" << std::endl;

  const char *testMolecules[][2] = {
      {"C1CC1CCC1CC(C)C1", "C1CC1CCC1CCC1"},
      {"CC1CCC1", "C1CCC1"},
      {"NCNCC2CC2C1CC1O", "C1CC1C1CC1"},
      {"OC2C(C)C21C(N)C1C", "C2CC12CC1"},  // Spiro
      {"C1CC1C(=O)OC", "C1CC1"},           // Carbonyl outside scaffold
      {"C1CC1C=C", "C1CC1"},               // Double bond outside scaffold
      {"C1CC1C=CC1CC1C=CNNCO", "C1CC1C=CC1CC1"},  // Double bond in scaffold
      {"CC1CC1C(N)C1C(N)C1", "C1CC1CC1CC1"},
      {"C1CC1S(=O)C1CC1C=CNNCO", "C1CC1S(=O)C1CC1"},  // S=O group in scaffold
      {"O=SCNC1CC1S(=O)C1CC1C=CNNCO",
       "C1CC1S(=O)C1CC1"},  // S=O group outside scaffold
      {"C1CC1S(=O)(=O)C1CC1C=CNNCO",
       "C1CC1S(=O)(=O)C1CC1"},  // SO2 group in scaffold
      {"O=S(CNCNC)(=O)CNC1CC1S(=O)(=O)C1CC1C=CNNCO",
       "C1CC1S(=O)(=O)C1CC1"},                    // SO2 group outside scaffold
      {"C1CC1C=NO", "C1CC1"},                     // Hydroxamide
      {"C1CC1C(C(C)C)=NC1CC1", "C1CC1C=NC1CC1"},  // Hydroxamide
      {"C1CC1C#N", "C1CC1"},                      // Cyano group
      {"C1CC1C#CNC", "C1CC1"},                    // Acetylene group
      {"O=C1N(C)C(=O)N1C#CNC", "O=C1NC(=O)N1"},   // Acetylene group
      {"[O-][N+](=O)c1cc(ccc1Cl)NS(=O)(=O)Cc2ccccc2",
       "c1ccccc1NS(=O)(=O)Cc2ccccc2"},
      {"Cn1cccc1", "c1ccc[nH]1"},
      {"C1CC1[CH](C)C1CC1", "C1CC1CC1CC1"},
      {"CC(C)c1c(Cl)cc(F)c(Nc2ccc(C)cc2CC(=O)O)c1F", "c1ccc(Nc2ccccc2)cc1"},
      {"C1CC1C[C@](Cl)(F)C1CCC1", "C1CC1CCC1CCC1"},
      {"C1CC1.C1CC1", "C1CC1.C1CC1"},  // this was a crash at one point
      {"CCC", ""},
      {"EOS", "EOS"}};
  unsigned int i = 0;
  while (1) {
    std::string smi = testMolecules[i][0];
    std::string tgt = testMolecules[i][1];
    ++i;
    if (smi == "EOS") {
      break;
    }
    ROMol *mol = SmilesToMol(smi);
    ROMol *nMol = MurckoDecompose(*mol);
    TEST_ASSERT(nMol);
    delete mol;
    if (tgt != "") {
      TEST_ASSERT(nMol->getNumAtoms());
      MolOps::sanitizeMol(static_cast<RWMol &>(*nMol));
      smi = MolToSmiles(*nMol, true);
      mol = SmilesToMol(tgt);
      TEST_ASSERT(mol);
      tgt = MolToSmiles(*mol, true);
      delete mol;
      TEST_ASSERT(smi == tgt);
    } else {
      TEST_ASSERT(nMol->getNumAtoms() == 0);
    }
    delete nMol;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testReplaceCoreRequireDummies() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing replaceCore requiring dummy atoms"
                       << std::endl;

  {
    std::string sma = "n1c(*)oc(*)c1";
    ROMol *matcher = SmartsToMol(sma);
    TEST_ASSERT(matcher);

    std::string smi = "n1c(CC)oc(C)c1";
    ROMol *mol1 = SmilesToMol(smi);
    TEST_ASSERT(mol1);

    bool replaceDummies = false;
    bool labelByIndex = true;
    bool requireDummyMatch = false;
    ROMol *mol2 = replaceCore(*mol1, *matcher, replaceDummies, labelByIndex,
                              requireDummyMatch);
    TEST_ASSERT(mol2);
    TEST_ASSERT(mol2->getNumAtoms() == 5);
    smi = MolToSmiles(*mol2, true);
    // If dummies existed and were replaced
    //  use THEIR atom indices
    TEST_ASSERT(smi == "[2*]CC.[5*]C");

    delete mol1;
    delete mol2;

    smi = "n1c(CC)oc(C)c1C";
    mol1 = SmilesToMol(smi);
    TEST_ASSERT(mol1);

    replaceDummies = false;
    labelByIndex = true;
    requireDummyMatch = true;
    mol2 = replaceCore(*mol1, *matcher, replaceDummies, labelByIndex,
                       requireDummyMatch);

    TEST_ASSERT(!mol2);
    delete mol1;

    smi = "n1c(CC)occ1C";
    mol1 = SmilesToMol(smi);
    TEST_ASSERT(mol1);

    replaceDummies = false;
    labelByIndex = true;
    requireDummyMatch = true;
    mol2 = replaceCore(*mol1, *matcher, replaceDummies, labelByIndex,
                       requireDummyMatch);

    TEST_ASSERT(!mol2);
    delete mol1;

    delete matcher;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue3453144() {
  ROMol *mol1 = nullptr, *matcher1 = nullptr, *replacement = nullptr;
  std::string smi, sma;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing that atom positions are correctly copied by replaceSubstructs"
      << std::endl;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ChemTransforms/testData/ethanol.mol";
  mol1 = MolFileToMol(pathName);
  TEST_ASSERT(mol1);
  TEST_ASSERT(mol1->getNumAtoms() == 3);

  sma = "O";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  smi = "N";
  replacement = SmilesToMol(smi);
  TEST_ASSERT(replacement);

  std::vector<ROMOL_SPTR> rs;
  rs = replaceSubstructs(*mol1, *matcher1, *replacement);
  TEST_ASSERT(rs.size() == 1);

  TEST_ASSERT(rs[0]->getNumAtoms() == 3);
  TEST_ASSERT(rs[0]->getAtomWithIdx(2)->getAtomicNum() == 7);

  TEST_ASSERT(rs[0]->getNumConformers() == mol1->getNumConformers());

  RDGeom::Point3D op, np;
  op = mol1->getConformer().getAtomPos(0);
  np = rs[0]->getConformer().getAtomPos(0);
  TEST_ASSERT(feq(op.x, np.x));
  TEST_ASSERT(feq(op.y, np.y));
  TEST_ASSERT(feq(op.z, np.z));
  op = mol1->getConformer().getAtomPos(1);
  np = rs[0]->getConformer().getAtomPos(1);
  TEST_ASSERT(feq(op.x, np.x));
  TEST_ASSERT(feq(op.y, np.y));
  TEST_ASSERT(feq(op.z, np.z));
  delete mol1;
  delete matcher1;
  delete replacement;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue3537675() {
  std::string smi, sma;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing handling of substituted non-carbon aromatic atoms"
      << std::endl;

  {
    std::string smi = "c1cccp1C";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    ROMol *nMol = MurckoDecompose(*mol);
    TEST_ASSERT(nMol->getAtomWithIdx(4)->getNumExplicitHs());
    delete mol;
    delete nMol;
  }
  {
    std::string smi = "c1cccn1C";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    ROMol *nMol = MurckoDecompose(*mol);
    TEST_ASSERT(nMol->getAtomWithIdx(4)->getNumExplicitHs());
    delete mol;
    delete nMol;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testCombineMols() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing combination of molecules" << std::endl;

  {
    std::string smi1 = "C1CCC1";
    std::string smi2 = "C1CC1";
    ROMol *mol1 = SmilesToMol(smi1);
    ROMol *mol2 = SmilesToMol(smi2);

    ROMol *mol3 = combineMols(*mol1, *mol2);
    TEST_ASSERT(mol3->getNumAtoms() ==
                (mol1->getNumAtoms() + mol2->getNumAtoms()));
    MolOps::findSSSR(*mol3);
    TEST_ASSERT(mol3->getRingInfo()->numRings() == 2);
    delete mol1;
    delete mol2;
    delete mol3;
  }
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/ChemTransforms/testData/ethanol.mol";
    ROMol *mol1 = MolFileToMol(pathName);
    TEST_ASSERT(mol1);
    ROMol *mol2 = MolFileToMol(pathName);
    TEST_ASSERT(mol2);
    TEST_ASSERT(mol1->getNumAtoms() == 3);
    TEST_ASSERT(mol2->getNumAtoms() == 3);
    ROMol *mol3 = combineMols(*mol1, *mol2);
    TEST_ASSERT(mol3->getNumAtoms() ==
                (mol1->getNumAtoms() + mol2->getNumAtoms()));
    TEST_ASSERT(mol3->getNumConformers() == 1);
    TEST_ASSERT(feq(mol3->getConformer().getAtomPos(0).x,
                    mol3->getConformer().getAtomPos(3).x));
    TEST_ASSERT(feq(mol3->getConformer().getAtomPos(0).y,
                    mol3->getConformer().getAtomPos(3).y));
    TEST_ASSERT(feq(mol3->getConformer().getAtomPos(0).z,
                    mol3->getConformer().getAtomPos(3).z));

    delete mol3;
    mol3 = combineMols(*mol1, *mol2, RDGeom::Point3D(1, 0, 0));
    TEST_ASSERT(mol3->getNumAtoms() ==
                (mol1->getNumAtoms() + mol2->getNumAtoms()));
    TEST_ASSERT(mol3->getNumConformers() == 1);
    TEST_ASSERT(!feq(mol3->getConformer().getAtomPos(0).x,
                     mol3->getConformer().getAtomPos(3).x));
    TEST_ASSERT(feq(mol3->getConformer().getAtomPos(0).y,
                    mol3->getConformer().getAtomPos(3).y));
    TEST_ASSERT(feq(mol3->getConformer().getAtomPos(0).z,
                    mol3->getConformer().getAtomPos(3).z));
    delete mol1;
    delete mol2;
    delete mol3;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testAddRecursiveQueries() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing adding recursive queries to a molecule"
                       << std::endl;

  {
    std::string smi1 = "CC";
    ROMol *mol1 = SmilesToMol(smi1);

    std::string smi2 = "CO";
    ROMOL_SPTR q1(SmilesToMol(smi2));
    std::map<std::string, ROMOL_SPTR> mp;
    mp["foo"] = q1;

    TEST_ASSERT(!mol1->getAtomWithIdx(0)->hasQuery());
    addRecursiveQueries(*mol1, mp, "replaceme");
    TEST_ASSERT(!mol1->getAtomWithIdx(0)->hasQuery());
    mol1->getAtomWithIdx(0)->setProp("replaceme", "foo");
    addRecursiveQueries(*mol1, mp, "replaceme");
    TEST_ASSERT(mol1->getAtomWithIdx(0)->hasQuery());
    TEST_ASSERT(mol1->getAtomWithIdx(0)->getQuery()->getDescription() ==
                "AtomAnd");
    TEST_ASSERT(!mol1->getAtomWithIdx(1)->hasQuery());

    mol1->getAtomWithIdx(0)->setProp("replaceme", "bar");
    bool ok = false;
    try {
      addRecursiveQueries(*mol1, mp, "replaceme");
    } catch (KeyErrorException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
    delete mol1;
  }

  {
    std::string smi1 = "CC";
    ROMol *mol1 = SmilesToMol(smi1);

    std::string smi2 = "CO";
    ROMOL_SPTR q1(SmilesToMol(smi2));
    std::map<std::string, ROMOL_SPTR> mp;
    mp["foo"] = q1;

    std::vector<std::pair<unsigned int, std::string>> labels;

    mol1->getAtomWithIdx(0)->setProp("replaceme", "foo");
    addRecursiveQueries(*mol1, mp, "replaceme", &labels);
    TEST_ASSERT(mol1->getAtomWithIdx(0)->hasQuery());
    TEST_ASSERT(labels.size() == 1);

    delete mol1;
  }

  {
    std::string smi1 = "CC";
    ROMol *mol1 = SmilesToMol(smi1);

    std::string smi2 = "CO";
    ROMOL_SPTR q1(SmilesToMol(smi2));
    std::map<std::string, ROMOL_SPTR> mp;
    mp["foo"] = q1;

    smi2 = "OC";
    ROMOL_SPTR q2(SmilesToMol(smi2));
    mp["bar"] = q2;

    std::vector<std::pair<unsigned int, std::string>> labels;

    mol1->getAtomWithIdx(0)->setProp("replaceme", "foo");
    mol1->getAtomWithIdx(1)->setProp("replaceme", "bar");
    addRecursiveQueries(*mol1, mp, "replaceme", &labels);
    TEST_ASSERT(mol1->getAtomWithIdx(0)->hasQuery());
    TEST_ASSERT(mol1->getAtomWithIdx(1)->hasQuery());
    TEST_ASSERT(mol1->getAtomWithIdx(0)->getQuery()->getDescription() ==
                "AtomAnd");
    TEST_ASSERT(mol1->getAtomWithIdx(1)->getQuery()->getDescription() ==
                "AtomAnd");
    TEST_ASSERT(labels.size() == 2);

    delete mol1;
  }

  {
    std::string smi1 = "CC";
    ROMol *mol1 = SmilesToMol(smi1);

    std::string smi2 = "CO";
    ROMOL_SPTR q1(SmilesToMol(smi2));
    std::map<std::string, ROMOL_SPTR> mp;
    mp["foo"] = q1;

    smi2 = "CN";
    ROMOL_SPTR q2(SmilesToMol(smi2));
    mp["bar"] = q2;

    mol1->getAtomWithIdx(0)->setProp("replaceme", "foo,bar");
    addRecursiveQueries(*mol1, mp, "replaceme");
    TEST_ASSERT(mol1->getAtomWithIdx(0)->hasQuery());
    TEST_ASSERT(!mol1->getAtomWithIdx(1)->hasQuery());
    TEST_ASSERT(mol1->getAtomWithIdx(0)->getQuery()->getDescription() ==
                "AtomAnd");

    MatchVectType mv;
    std::string msmi = "CCC";
    ROMol *mmol = SmilesToMol(msmi);
    TEST_ASSERT(mmol);
    TEST_ASSERT(!SubstructMatch(*mmol, *mol1, mv));
    delete mmol;

    msmi = "CCO";
    mmol = SmilesToMol(msmi);
    TEST_ASSERT(mmol);
    TEST_ASSERT(SubstructMatch(*mmol, *mol1, mv));
    delete mmol;

    msmi = "CCN";
    mmol = SmilesToMol(msmi);
    TEST_ASSERT(mmol);
    TEST_ASSERT(SubstructMatch(*mmol, *mol1, mv));
    delete mmol;

    delete mol1;
  }

  // Tests smarts as props
  {
    std::string smi1 = "CC";
    ROMol *mol1 = SmilesToMol(smi1);

    std::map<std::string, ROMOL_SPTR> mp;
    mol1->getAtomWithIdx(0)->setProp("replaceme", "CO");
    addRecursiveQueries(*mol1, mp, "replaceme");
    TEST_ASSERT(mol1->getAtomWithIdx(0)->hasQuery());
    TEST_ASSERT(!mol1->getAtomWithIdx(1)->hasQuery());
    TEST_ASSERT(mol1->getAtomWithIdx(0)->getQuery()->getDescription() ==
                "AtomAnd");

    MatchVectType mv;
    std::string msmi = "CCC";
    ROMol *mmol = SmilesToMol(msmi);
    TEST_ASSERT(mmol);
    TEST_ASSERT(!SubstructMatch(*mmol, *mol1, mv));
    delete mmol;

    msmi = "CCO";
    mmol = SmilesToMol(msmi);
    TEST_ASSERT(mmol);
    TEST_ASSERT(SubstructMatch(*mmol, *mol1, mv));
    delete mmol;

    delete mol1;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testParseQueryDefFile() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing parsing query definition file" << std::endl;

  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/ChemTransforms/testData/query_file1.txt";
    std::map<std::string, ROMOL_SPTR> qdefs;
    parseQueryDefFile(pathName, qdefs, false);
    TEST_ASSERT(!qdefs.empty());
    TEST_ASSERT(qdefs.size() == 7);
    TEST_ASSERT(qdefs.find("AcidChloride") != qdefs.end());
    TEST_ASSERT(qdefs.find("AcideChloride") == qdefs.end());
    TEST_ASSERT(qdefs.find("AcidChloride.Aliphatic") != qdefs.end());
    TEST_ASSERT(qdefs.find("CarboxylicAcid.AlphaAmino") != qdefs.end());

    std::string msmi = "CCC(=O)Cl";
    ROMol *mmol = SmilesToMol(msmi);
    TEST_ASSERT(mmol);
    MatchVectType mv;
    TEST_ASSERT(SubstructMatch(*mmol, *(qdefs["AcidChloride"]), mv));
    delete mmol;
  }
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/ChemTransforms/testData/query_file2.txt";
    std::map<std::string, ROMOL_SPTR> qdefs;
    parseQueryDefFile(pathName, qdefs);
    TEST_ASSERT(qdefs.empty());
  }
  /*{
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/ChemTransforms/testData/query_file1.txt";
    std::map<std::string, ROMOL_SPTR> qdefs;
    parseQueryDefFile(pathName, qdefs);
    TEST_ASSERT(!qdefs.empty());
    TEST_ASSERT(qdefs.size()==7);
    TEST_ASSERT(qdefs.find("acidchloride")!=qdefs.end());
    TEST_ASSERT(qdefs.find("AcidChloride")==qdefs.end());
  }*/
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue275() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing sf.net issue 275: Murcko decomposition with chiral atoms"
      << std::endl;

  {
    std::string smi = "CCCCC[C@H]1CC[C@H](C(=O)O)CC1";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    ROMol *nMol = MurckoDecompose(*mol);
    smi = MolToSmiles(*nMol, true);
    TEST_ASSERT(smi == "C1CCCCC1");
    delete mol;
    delete nMol;
  }

  {
    std::string smi = "CCCCC[C@H]1CC[C@H](C(=O)O)CC1";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    ROMol *nMol = MurckoDecompose(*mol);
    smi = MolToSmiles(*nMol, false);
    TEST_ASSERT(smi == "C1CCCCC1");
    delete mol;
    delete nMol;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testFragmentOnBonds() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing fragmentOnBonds" << std::endl;

  {
    std::string smi = "OCCCN";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 5);
    unsigned int indices[] = {0, 3};
    std::vector<unsigned int> bindices(
        indices, indices + (sizeof(indices) / sizeof(indices[0])));
    ROMol *nmol = MolFragmenter::fragmentOnBonds(*mol, bindices, false);
    TEST_ASSERT(nmol);
    TEST_ASSERT(nmol->getNumAtoms() == 5);
    delete mol;
    delete nmol;
  }
  {
    std::string smi = "OCCCN";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 5);
    unsigned int indices[] = {0, 1};
    std::vector<unsigned int> bindices(
        indices, indices + (sizeof(indices) / sizeof(indices[0])));
    std::vector<unsigned int> cutsPerAtom(mol->getNumAtoms());
    ROMol *nmol = MolFragmenter::fragmentOnBonds(*mol, bindices, false, nullptr,
                                                 nullptr, &cutsPerAtom);
    TEST_ASSERT(nmol);
    TEST_ASSERT(nmol->getNumAtoms() == 5);
    TEST_ASSERT(cutsPerAtom[0] == 1);
    TEST_ASSERT(cutsPerAtom[1] == 2);
    TEST_ASSERT(cutsPerAtom[2] == 1);
    TEST_ASSERT(cutsPerAtom[3] == 0);
    delete mol;
    delete nmol;
  }
  {
    std::string smi = "OCCCN";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 5);
    TEST_ASSERT(mol->getBondBetweenAtoms(0, 1));
    TEST_ASSERT(mol->getBondBetweenAtoms(3, 4));
    unsigned int indices[] = {0, 3};
    std::vector<unsigned int> bindices(
        indices, indices + (sizeof(indices) / sizeof(indices[0])));
    ROMol *nmol = MolFragmenter::fragmentOnBonds(*mol, bindices);
    TEST_ASSERT(nmol);
    TEST_ASSERT(nmol->getNumAtoms() == 9);
    TEST_ASSERT(!nmol->getBondBetweenAtoms(0, 1));
    TEST_ASSERT(!nmol->getBondBetweenAtoms(3, 4));
    TEST_ASSERT(nmol->getAtomWithIdx(5)->getAtomicNum() == 0);
    TEST_ASSERT(nmol->getAtomWithIdx(5)->getIsotope() == 0);
    TEST_ASSERT(nmol->getBondBetweenAtoms(1, 5));
    TEST_ASSERT(nmol->getAtomWithIdx(6)->getAtomicNum() == 0);
    TEST_ASSERT(nmol->getAtomWithIdx(6)->getIsotope() == 1);
    TEST_ASSERT(nmol->getBondBetweenAtoms(0, 6));

    TEST_ASSERT(nmol->getAtomWithIdx(7)->getAtomicNum() == 0);
    TEST_ASSERT(nmol->getAtomWithIdx(7)->getIsotope() == 3);
    TEST_ASSERT(nmol->getBondBetweenAtoms(4, 7));
    TEST_ASSERT(nmol->getAtomWithIdx(8)->getAtomicNum() == 0);
    TEST_ASSERT(nmol->getAtomWithIdx(8)->getIsotope() == 4);
    TEST_ASSERT(nmol->getBondBetweenAtoms(3, 8));
    delete mol;
    delete nmol;
  }
  {
    std::string smi = "OCCCN";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 5);
    TEST_ASSERT(mol->getBondBetweenAtoms(0, 1));
    TEST_ASSERT(mol->getBondBetweenAtoms(3, 4));
    unsigned int indices[] = {0, 3};
    std::vector<unsigned int> bindices(
        indices, indices + (sizeof(indices) / sizeof(indices[0])));
    std::vector<std::pair<unsigned int, unsigned int>> dummyLabels(2);
    dummyLabels[0] = std::make_pair(10, 11);
    dummyLabels[1] = std::make_pair(100, 110);
    ROMol *nmol =
        MolFragmenter::fragmentOnBonds(*mol, bindices, true, &dummyLabels);
    TEST_ASSERT(nmol);
    TEST_ASSERT(nmol->getNumAtoms() == 9);
    TEST_ASSERT(!nmol->getBondBetweenAtoms(0, 1));
    TEST_ASSERT(!nmol->getBondBetweenAtoms(3, 4));
    TEST_ASSERT(nmol->getAtomWithIdx(5)->getAtomicNum() == 0);
    TEST_ASSERT(nmol->getAtomWithIdx(5)->getIsotope() == 10);
    TEST_ASSERT(nmol->getBondBetweenAtoms(1, 5));
    TEST_ASSERT(nmol->getAtomWithIdx(6)->getAtomicNum() == 0);
    TEST_ASSERT(nmol->getAtomWithIdx(6)->getIsotope() == 11);
    TEST_ASSERT(nmol->getBondBetweenAtoms(0, 6));

    TEST_ASSERT(nmol->getAtomWithIdx(7)->getAtomicNum() == 0);
    TEST_ASSERT(nmol->getAtomWithIdx(7)->getIsotope() == 100);
    TEST_ASSERT(nmol->getBondBetweenAtoms(4, 7));
    TEST_ASSERT(nmol->getAtomWithIdx(8)->getAtomicNum() == 0);
    TEST_ASSERT(nmol->getAtomWithIdx(8)->getIsotope() == 110);
    TEST_ASSERT(nmol->getBondBetweenAtoms(3, 8));
    delete mol;
    delete nmol;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testFragmentOnBRICSBonds() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing fragmentOnBRICSBonds" << std::endl;

  {
    std::string smi = "c1ccccc1OC";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 8);
    TEST_ASSERT(mol->getBondBetweenAtoms(5, 6));
    TEST_ASSERT(mol->getBondBetweenAtoms(6, 7));

    std::vector<MolFragmenter::FragmenterBondType> fbts;
    MolFragmenter::constructBRICSBondTypes(fbts);
    ROMol *nmol = MolFragmenter::fragmentOnBonds(*mol, fbts);
    TEST_ASSERT(nmol);

    TEST_ASSERT(nmol->getNumAtoms() == 10);
    TEST_ASSERT(!nmol->getBondBetweenAtoms(5, 6));
    TEST_ASSERT(nmol->getBondBetweenAtoms(6, 7));

    smi = MolToSmiles(*nmol, true);
    TEST_ASSERT(smi == "[16*]c1ccccc1.[3*]OC");

    TEST_ASSERT(nmol->getAtomWithIdx(8)->getAtomicNum() == 0);
    TEST_ASSERT(nmol->getAtomWithIdx(8)->getIsotope() == 3);
    TEST_ASSERT(nmol->getBondBetweenAtoms(6, 8));
    TEST_ASSERT(nmol->getAtomWithIdx(9)->getAtomicNum() == 0);
    TEST_ASSERT(nmol->getAtomWithIdx(9)->getIsotope() == 16);
    TEST_ASSERT(nmol->getBondBetweenAtoms(5, 9));

    delete mol;
    delete nmol;
  }

  {
    std::string smi = "c1ccccc1";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 6);

    std::vector<MolFragmenter::FragmenterBondType> fbts;
    MolFragmenter::constructBRICSBondTypes(fbts);
    ROMol *nmol = MolFragmenter::fragmentOnBonds(*mol, fbts);
    TEST_ASSERT(nmol);

    TEST_ASSERT(nmol->getNumAtoms() == 6);

    smi = MolToSmiles(*nmol, true);
    TEST_ASSERT(smi == "c1ccccc1");

    delete mol;
    delete nmol;
  }

  {
    std::string smi = "OC(C)=CC";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 5);

    std::vector<MolFragmenter::FragmenterBondType> fbts;
    MolFragmenter::constructBRICSBondTypes(fbts);
    ROMol *nmol = MolFragmenter::fragmentOnBonds(*mol, fbts);
    TEST_ASSERT(nmol);

    TEST_ASSERT(nmol->getNumAtoms() == 7);

    smi = MolToSmiles(*nmol, true);
    TEST_ASSERT(smi == "[7*]=C(C)O.[7*]=CC");

    delete mol;
    delete nmol;
  }

  {
    std::string smi = "c1ccccc1OC";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 8);
    ROMol *nmol = MolFragmenter::fragmentOnBRICSBonds(*mol);
    TEST_ASSERT(nmol);
    TEST_ASSERT(nmol->getNumAtoms() == 10);
    smi = MolToSmiles(*nmol, true);
    TEST_ASSERT(smi == "[16*]c1ccccc1.[3*]OC");

    delete mol;
    delete nmol;
  }
  {
    std::string smi = "OC(C)=CC";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 5);

    ROMol *nmol = MolFragmenter::fragmentOnBRICSBonds(*mol);
    TEST_ASSERT(nmol);

    TEST_ASSERT(nmol->getNumAtoms() == 7);
    smi = MolToSmiles(*nmol, true);
    TEST_ASSERT(smi == "[7*]=C(C)O.[7*]=CC");

    delete mol;
    delete nmol;
  }

  {
    std::string smi = "CCCOCCC(=O)c1ccccc1";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 14);

    ROMol *nmol = MolFragmenter::fragmentOnBRICSBonds(*mol);
    TEST_ASSERT(nmol);

    TEST_ASSERT(nmol->getNumAtoms() == 20);
    smi = MolToSmiles(*nmol, true);
    TEST_ASSERT(smi == "[16*]c1ccccc1.[3*]O[3*].[4*]CCC.[4*]CCC([6*])=O");
    MolOps::sanitizeMol(static_cast<RWMol &>(*nmol));
    smi = MolToSmiles(*nmol, true);
    TEST_ASSERT(smi == "[16*]c1ccccc1.[3*]O[3*].[4*]CCC.[4*]CCC([6*])=O");

    delete mol;
    delete nmol;
  }

  {
    std::string smi = "Cl.CC(=O)O[C@]1(c2ccccc2)CCN(C)[C@H]2CCCC[C@@H]21";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 22);
    ROMol *nmol = MolFragmenter::fragmentOnBRICSBonds(*mol);
    TEST_ASSERT(nmol);

    TEST_ASSERT(nmol->getNumAtoms() == 28);
    smi = MolToSmiles(*nmol, true);
    // std::cerr<<smi<<std::endl;
    TEST_ASSERT(smi ==
                "Cl.[1*]C(C)=O.[15*]C1([15*])CCN(C)[C@H]2CCCC[C@@H]21.[16*]"
                "c1ccccc1.[3*]O[3*]");
    MolOps::sanitizeMol(static_cast<RWMol &>(*nmol));
    smi = MolToSmiles(*nmol, true);
    TEST_ASSERT(smi ==
                "Cl.[1*]C(C)=O.[15*]C1([15*])CCN(C)[C@H]2CCCC[C@@H]21.[16*]"
                "c1ccccc1.[3*]O[3*]");

    delete mol;
    delete nmol;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void benchFragmentOnBRICSBonds() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing fragmentOnBRICSBonds" << std::endl;
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Regress/Data/mols.1000.sdf";
    SDMolSupplier suppl(pathName);
    while (!suppl.atEnd()) {
      ROMol *m = suppl.next();
      ROMol *nmol = MolFragmenter::fragmentOnBRICSBonds(*m);
      delete m;
      delete nmol;
    }
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testFragmentOnSomeBonds() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing fragmentOnSomeBonds" << std::endl;

  {
    std::string smi = "OCCCCN";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms());
    unsigned int indices[] = {0, 2, 4};
    std::vector<unsigned int> bindices(
        indices, indices + (sizeof(indices) / sizeof(indices[0])));
    std::vector<ROMOL_SPTR> frags;
    MolFragmenter::fragmentOnSomeBonds(*mol, bindices, frags, 2);
    TEST_ASSERT(frags.size() == 3);
    std::vector<std::vector<int>> fragMap;

    TEST_ASSERT(MolOps::getMolFrags(*frags[0], fragMap) == 3);
    TEST_ASSERT(fragMap.size() == 3);
    TEST_ASSERT(fragMap[0].size() == 2);
    TEST_ASSERT(fragMap[1].size() == 4);
    TEST_ASSERT(fragMap[2].size() == 4);

    TEST_ASSERT(MolOps::getMolFrags(*frags[1], fragMap) == 3);
    TEST_ASSERT(fragMap.size() == 3);
    TEST_ASSERT(fragMap[0].size() == 2);
    TEST_ASSERT(fragMap[1].size() == 6);
    TEST_ASSERT(fragMap[2].size() == 2);

    TEST_ASSERT(MolOps::getMolFrags(*frags[2], fragMap) == 3);
    TEST_ASSERT(fragMap.size() == 3);
    TEST_ASSERT(fragMap[0].size() == 4);
    TEST_ASSERT(fragMap[1].size() == 4);
    TEST_ASSERT(fragMap[2].size() == 2);

    delete mol;
  }

  {
    std::string smi = "OCCCCN";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms());
    unsigned int indices[] = {0, 2, 4};
    std::vector<unsigned int> bindices(
        indices, indices + (sizeof(indices) / sizeof(indices[0])));
    std::vector<ROMOL_SPTR> frags;
    std::vector<std::vector<unsigned int>> cpa;
    MolFragmenter::fragmentOnSomeBonds(*mol, bindices, frags, 2, false, nullptr,
                                       nullptr, &cpa);
    TEST_ASSERT(frags.size() == 3);
    TEST_ASSERT(cpa.size() == 3);
    TEST_ASSERT(cpa[0].size() == mol->getNumAtoms());
    TEST_ASSERT(cpa[1].size() == mol->getNumAtoms());
    TEST_ASSERT(cpa[2].size() == mol->getNumAtoms());

    delete mol;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithubIssue429() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github issue 429: fragmentOnSomeBonds() "
                          "should update implicit H count on aromatic "
                          "heteroatoms when addDummies is False"
                       << std::endl;

  {
    std::string smi = "c1cccn1CC1CC1";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms());
    unsigned int indices[] = {4, 5};
    std::vector<unsigned int> bindices(
        indices, indices + (sizeof(indices) / sizeof(indices[0])));
    std::vector<ROMOL_SPTR> frags;
    MolFragmenter::fragmentOnSomeBonds(*mol, bindices, frags, 1, false);
    TEST_ASSERT(frags.size() == 2);
    std::vector<std::vector<int>> fragMap;

    for (auto romol : frags) {
      auto *rwmol = (RWMol *)(romol.get());
      MolOps::sanitizeMol(*rwmol);
    }

    // we actually changed fragmentOnBonds(), check that too:
    auto *nmol = (RWMol *)MolFragmenter::fragmentOnBonds(*mol, bindices, false);
    MolOps::sanitizeMol(*nmol);
    delete nmol;

    delete mol;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}
void testGithubIssue430() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github issue 430: fragmentOnSomeBonds() "
                          "crashes if bond list is empty"
                       << std::endl;

  {
    std::string smi = "c1cccn1CC1CC1";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms());
    std::vector<unsigned int> bindices;
    std::vector<ROMOL_SPTR> frags;
    MolFragmenter::fragmentOnSomeBonds(*mol, bindices, frags, 1, false);
    TEST_ASSERT(frags.size() == 0);
    delete mol;
  }

  {
    std::string smi = "c1cccn1CC1CC1";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms());
    unsigned int indices[] = {4, 5};
    std::vector<unsigned int> bindices(
        indices, indices + (sizeof(indices) / sizeof(indices[0])));
    std::vector<ROMOL_SPTR> frags;
    MolFragmenter::fragmentOnSomeBonds(*mol, bindices, frags, 0, false);
    TEST_ASSERT(frags.size() == 0);
    delete mol;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithubIssue511() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github issue 511: FragmentOnBonds() "
                          "producing incorrect chirality"
                       << std::endl;

  {  // start with an example with no ring bonds... just to dodge that
     // complication
    std::string smi = "CO[C@](CC)(C)N";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 7)
    std::vector<unsigned int> bindices;
    bindices.push_back(1);
    ROMol *frags = MolFragmenter::fragmentOnBonds(*mol, bindices, true);
    TEST_ASSERT(frags->getNumAtoms() == 9);

    std::string csmi1 = MolToSmiles(*mol, true);
    // std::cerr << csmi1 << std::endl;

    TEST_ASSERT(csmi1 == "CC[C@@](C)(N)OC");
    std::string csmi2 = MolToSmiles(*frags, true);
    // std::cerr << csmi2 << std::endl;
    TEST_ASSERT(csmi2 == "[1*][C@](C)(N)CC.[2*]OC");

    delete mol;
    delete frags;
  }

  {  // Steve's original bug report
    std::string smi = "CO[C@]1(C)NCC1";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 7)
    std::vector<unsigned int> bindices;
    bindices.push_back(1);
    ROMol *frags = MolFragmenter::fragmentOnBonds(*mol, bindices, true);
    TEST_ASSERT(frags->getNumAtoms() == 9);

    std::string csmi1 = MolToSmiles(*mol, true);
    TEST_ASSERT(csmi1 == "CO[C@@]1(C)CCN1");
    std::string csmi2 = MolToSmiles(*frags, true);
    // std::cerr << csmi2 << std::endl;
    TEST_ASSERT(csmi2 == "[1*][C@@]1(C)CCN1.[2*]OC");

    delete mol;
    delete frags;
  }

  {  // further complicated by adding stereochem on both sides
    std::string smi = "C[C@](O)(F)[C@@]1(C)CCN1";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 9)
    std::vector<unsigned int> bindices;
    bindices.push_back(3);
    ROMol *frags = MolFragmenter::fragmentOnBonds(*mol, bindices, true);
    TEST_ASSERT(frags->getNumAtoms() == 11);

    std::string csmi1 = MolToSmiles(*mol, true);
    // std::cerr << csmi1 << std::endl;
    TEST_ASSERT(csmi1 == "C[C@](O)(F)[C@@]1(C)CCN1");
    std::string csmi2 = MolToSmiles(*frags, true);
    // std::cerr << csmi2 << std::endl;
    TEST_ASSERT(csmi2 == "[1*][C@@]1(C)CCN1.[4*][C@@](C)(O)F");

    delete mol;
    delete frags;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub1734() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github issue 1734: BreakBRICSBonds() not "
                          "preserving stereochemistry"
                       << std::endl;

  {  // no implicit H on chiral center
    std::string smi = "c1ccccc1[C@]([I])(C)NC";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 11)
    ROMol *nmol = MolFragmenter::fragmentOnBRICSBonds(*mol);
    TEST_ASSERT(nmol);
    TEST_ASSERT(nmol->getNumAtoms() == 15);
    smi = MolToSmiles(*nmol, true);
    TEST_ASSERT(smi == "[16*]c1ccccc1.[4*][C@]([8*])(C)I.[5*]NC");
    delete mol;
    delete nmol;
  }
  {  // The original example
    std::string smi = "c1ccccc1[C@H](C)NC";
    RWMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 10)
    ROMol *nmol = MolFragmenter::fragmentOnBRICSBonds(*mol);
    TEST_ASSERT(nmol);
    TEST_ASSERT(nmol->getNumAtoms() == 14);
    smi = MolToSmiles(*nmol, true);
    TEST_ASSERT(smi == "[16*]c1ccccc1.[4*][C@H]([8*])C.[5*]NC");
    delete mol;
    delete nmol;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub3206() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing GitHub #3206: Queries generated from "
                          "PreprocessReaction cannot be translated to SMARTS"
                       << std::endl;

  {
    auto mol1 = "CC"_smiles;

    std::map<std::string, ROMOL_SPTR> mp;
    mp["foo"] = ROMOL_SPTR(SmilesToMol("CO"));
    mp["bar"] = ROMOL_SPTR(SmilesToMol("CN"));
    mp["baz"] = ROMOL_SPTR(SmilesToMol("CF"));

    TEST_ASSERT(!mol1->getAtomWithIdx(0)->hasQuery());
    addRecursiveQueries(*mol1, mp, "replaceme");
    TEST_ASSERT(!mol1->getAtomWithIdx(0)->hasQuery());
    mol1->getAtomWithIdx(0)->setProp("replaceme", "foo,bar,baz");
    addRecursiveQueries(*mol1, mp, "replaceme");
    TEST_ASSERT(mol1->getAtomWithIdx(0)->hasQuery());
    TEST_ASSERT(mol1->getAtomWithIdx(0)->getQuery()->getDescription() ==
                "AtomAnd");
    TEST_ASSERT(!mol1->getAtomWithIdx(1)->hasQuery());

    auto sma = MolToSmarts(*mol1);
    TEST_ASSERT(sma == "[#6;$([#6]-[#8]),$([#6]-[#7]),$([#6]-[#9])]-[#6]");
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub4019() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing GitHub #4019: dummy atoms should not be marked "
      << "as aromatic, not have explicit Hs, not "
      << "be bonded through an aromatic bond and "
      << "not be bonded with each other" << std::endl;
  {
    auto mol = "c1ncccc1n1ccc2ccccc12"_smiles;
    ROMOL_SPTR core(SmartsToMol("n1ccc2ccccc12"));
    ROMOL_SPTR molNoSidechain(replaceSidechains(*mol, *core));
    bool hasDummy = false;
    for (auto a : molNoSidechain->atoms()) {
      if (a->getAtomicNum() == 0) {
        hasDummy = true;
        TEST_ASSERT(!a->getIsAromatic());
        a->setAtomicNum(1);
        a->setIsotope(0);
      }
    }
    TEST_ASSERT(hasDummy);
    molNoSidechain.reset(MolOps::removeHs(*molNoSidechain));
    TEST_ASSERT(MolToSmiles(*molNoSidechain) == "c1ccc2[nH]ccc2c1");
  }
  {
    auto mol = "c1ccc2[nH]ccc2c1"_smiles;
    ROMOL_SPTR core(SmartsToMol("c1ccccc1"));
    ROMOL_SPTR molNoSidechain(replaceSidechains(*mol, *core));
    unsigned int nDummies = 0;
    for (auto a : molNoSidechain->atoms()) {
      if (a->getAtomicNum() == 0) {
        ++nDummies;
        TEST_ASSERT(!a->getIsAromatic());
        TEST_ASSERT(!a->getNumExplicitHs());
      }
    }
    TEST_ASSERT(nDummies == 2);
    for (auto b : molNoSidechain->bonds()) {
      if (b->getBeginAtom()->getAtomicNum() == 0 ||
          b->getEndAtom()->getAtomicNum() == 0) {
        TEST_ASSERT(!b->getIsAromatic());
        TEST_ASSERT(b->getBondType() == Bond::SINGLE);
      }
    }
  }
  {
    auto mol = "c1ccccc1C1CC1"_smiles;
    ROMOL_SPTR core(SmartsToMol("c1ccccc1C"));
    ROMOL_SPTR molNoSidechain(replaceSidechains(*mol, *core));
    std::vector<unsigned int> dummies;
    for (auto a : molNoSidechain->atoms()) {
      if (a->getAtomicNum() == 0) {
        dummies.push_back(a->getIdx());
        TEST_ASSERT(!a->getIsAromatic());
        TEST_ASSERT(!a->getNumExplicitHs());
      }
    }
    TEST_ASSERT(dummies.size() == 2);
    TEST_ASSERT(
        !molNoSidechain->getBondBetweenAtoms(dummies.front(), dummies.back()));
  }
}

int main() {
  RDLog::InitLogs();

  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing Chemical Transforms \n";

#if 1
  testDeleteSubstruct();
  testReplaceSubstructs();
  testReplaceSubstructs2();
  testReplaceSidechains();
  testReplaceCore();
  testReplaceCoreLabels();
  testReplaceCoreCrash();
  testReplaceCorePositions();
  testReplaceCoreMatchVect();
  testReplaceCoreMatchVectMultipleMappingToCore();

  testMurckoDecomp();
  testReplaceCoreRequireDummies();

  testIssue3453144();
  testIssue3537675();

  testCombineMols();
  testAddRecursiveQueries();
  testParseQueryDefFile();
  testIssue275();

  testFragmentOnBonds();
  testFragmentOnBRICSBonds();
  testFragmentOnSomeBonds();
  // benchFragmentOnBRICSBonds();
  testGithubIssue429();
  testGithubIssue430();
  testGithubIssue511();
  testReplaceCore2();
#endif
  testGithub1734();
  testGithub3206();
  testGithub4019();
  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  return (0);
}
