//  $Id$
// 
//   Copyright (C) 2007-2010 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
// There are chirality test cases spread all over the place. Many of the
// tests here are repeats, but it's good to have everything in one place.
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
//#include <boost/log/functions.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Canon.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>

#include <iostream>

using namespace RDKit;
using namespace std;

void testMol1(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "CIP codes from a mol file (1)" << std::endl;
  std::string rdbase = getenv("RDBASE");
  RWMol *m;
  std::string fName,smi;
  std::string cip;

  // start with SMILES:
  BOOST_LOG(rdInfoLog) << " >>>>>>>>>>>>> smiles 1 <<<<<<<<<<<<<< " << std::endl;
  smi="O[C@@H](N)I";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  smi="[C@H](O)(N)I";
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  MolOps::removeStereochemistry(*m);
  TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp("_CIPCode"));

  BOOST_LOG(rdInfoLog) << " >>>>>>>>>>>>> mol file <<<<<<<<<<<<<< " << std::endl;
  delete m;
  fName = rdbase+"/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir1a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete m;
  fName = rdbase+"/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir1b.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  MolOps::removeStereochemistry(*m);
  TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp("_CIPCode"));

  delete m;
  fName = rdbase+"/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir2a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  MolOps::removeStereochemistry(*m);
  TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp("_CIPCode"));

  delete m;
  fName = rdbase+"/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir2b.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  MolOps::removeStereochemistry(*m);
  TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp("_CIPCode"));

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
};


void testRoundTrip(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "CIP codes from a mol->smiles conversion (1)" << std::endl;
  std::string rdbase = getenv("RDBASE");
  RWMol *m;
  std::string fName,smi,smi2;
  std::string cip;

  // start with SMILES:
  smi="O[C@@H](N)I";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  smi=MolToSmiles(*m,true);
  std::cerr<<"smi: "<<smi<<std::endl;
  delete m;
  m=SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  smi2=MolToSmiles(*m,true);
  TEST_ASSERT(smi==smi2);

  
  smi="[C@H](O)(N)I";
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
#if 1
  smi=MolToSmiles(*m,true);
  BOOST_LOG(rdInfoLog)<<"smiout: "<<smi<<std::endl;
  TEST_ASSERT(smi=="N[C@H](O)I");
  delete m;
  m=SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  smi2=MolToSmiles(*m,true);
  TEST_ASSERT(smi==smi2);
#endif

  BOOST_LOG(rdInfoLog) << " >>>>>>>>>>>>> mol file <<<<<<<<<<<<<< " << std::endl;
  delete m;
  fName = rdbase+"/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir1a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
#if 1
  smi=MolToSmiles(*m,true);
  delete m;
  m=SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  smi2=MolToSmiles(*m,true);
  TEST_ASSERT(smi==smi2);
#endif

  delete m;
  fName = rdbase+"/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir1b.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
#if 1
  smi=MolToSmiles(*m,true);
  delete m;
  m=SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  smi2=MolToSmiles(*m,true);
  TEST_ASSERT(smi==smi2);
#endif

  delete m;
  fName = rdbase+"/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir2a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
#if 1
  smi=MolToSmiles(*m,true);
  delete m;
  m=SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  smi2=MolToSmiles(*m,true);
  TEST_ASSERT(smi==smi2);
#endif

  delete m;
  fName = rdbase+"/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir2b.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
#if 1
  smi=MolToSmiles(*m,true);
  delete m;
  m=SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  smi2=MolToSmiles(*m,true);
  TEST_ASSERT(smi==smi2);
#endif

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
};


void testMol2(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "CIP codes from a mol file (2)" << std::endl;
  std::string rdbase = getenv("RDBASE");
  RWMol *m;
  std::string fName,smi;
  std::string cip;

  // start with SMILES:
  BOOST_LOG(rdInfoLog) << " >>>>>>>>>>>>> smiles 1 <<<<<<<<<<<<<< " << std::endl;
  smi="[C@]1(SC[C@@]([H])(F)[C@]1(Br)O)([I])[H]";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==9);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
  m->getAtomWithIdx(3)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  TEST_ASSERT(m->getAtomWithIdx(5)->hasProp("_CIPCode"));
  m->getAtomWithIdx(5)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  // same molecule, H combined with the first atom (reproduces 
  // exact situation in upcoming mol file)
  BOOST_LOG(rdInfoLog) << " >>>>>>>>>>>>> smiles 2 <<<<<<<<<<<<<< " << std::endl;
  delete m;
  smi="[C@@H]1(SC[C@@]([H])(F)[C@]1(Br)O)([I])";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==9);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
  m->getAtomWithIdx(3)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  TEST_ASSERT(m->getAtomWithIdx(5)->hasProp("_CIPCode"));
  m->getAtomWithIdx(5)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  
  delete m;
  BOOST_LOG(rdInfoLog) << " >>>>>>>>>>>>> mol file <<<<<<<<<<<<<< " << std::endl;
  fName = rdbase+"/Code/GraphMol/FileParsers/test_data/Issue142b.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==9);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
  m->getAtomWithIdx(3)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
};

void testSmiles1(){
  ROMol *mol;
  std::string smi,cip;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "CIP codes from SMILES" << std::endl;

  smi = "F[C@](Cl)(Br)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_UNSPECIFIED);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "F[C@](Br)(I)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_UNSPECIFIED);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "F[C@](I)(Cl)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "Cl[C@](Br)(F)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "Cl[C@](F)(I)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "I[C@](F)(Br)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "I[C@](Br)(Cl)F";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "F[C@@](Br)(Cl)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "F[C@@](Cl)(I)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "Cl[C@@](Br)(I)F";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "Cl[C@@](F)(Br)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "[C@@](Cl)(F)(Br)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "F[C@H](Cl)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "Br[C@H](F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "Br[C@]([H])(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "Br[C@](F)(Cl)[H]";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "Br[C@]1(F)(Cl).[H]1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "Br[C@H]1Cl.F1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "Br[C@]12Cl.F2.[H]1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "Br[C@]21Cl.F1.[H]2";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "Br[C@]12Cl.F1.[H]2";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "[C@@](C)(Br)(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "[C@@]([H])(Br)(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "[C@@H](Br)(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "[H][C@@](Br)(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


void testChiralityCleanup(){
  ROMol *mol,*mol2;
  std::string smi,cip;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "chirality cleanup" << std::endl;

  smi = "F[C@H+](Cl)(Br)I";
  mol = SmilesToMol(smi,false,false);
  mol2 = MolOps::removeHs(*mol,false,false);
  delete mol;
  mol=mol2;
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol,true);
  TEST_ASSERT(!mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_UNSPECIFIED);
  delete mol;

  smi = "F[C@+](C)(Cl)(Br)I";
  mol = SmilesToMol(smi,false,false);
  mol2 = MolOps::removeHs(*mol,false,false);
  delete mol;
  mol=mol2;
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol,true);
  TEST_ASSERT(!mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_UNSPECIFIED);
  delete mol;


  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


void testRingStereochemistry(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test ring stereochemistry " << std::endl;
  // NOTE: this test is for correctness, not canonicality
  {
    std::string smi = "B[C@H]1CC[C@H](C)CC1";
    RWMol *m = SmilesToMol(smi);
    std::string smi1=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi1=="B[C@H]1CC[C@H](C)CC1");
    delete m;
#if 0
    smi="B[C@@H]1CC[C@@H](C)CC1";
    m = SmilesToMol(smi);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi2<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi2==smi);
    delete m;
#endif
  }

  {
    std::string smi = "C1[C@@H](B)CC[C@H](C)C1";
    RWMol *m = SmilesToMol(smi);
    std::string smi1=MolToSmiles(*m,true);
    smi = "B[C@@H]1CC[C@@H](C)CC1";
    BOOST_LOG(rdInfoLog)<<" : "<<smi<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi1==smi);
    delete m;
#if 0
    smi="C1[C@H](B)CC[C@@H](C)C1";
    m = SmilesToMol(smi);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi2<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi2==smi1);
    delete m;
#endif
  }

  {
    std::string smi = "C[C@H]1CC[C@H](F)CC1";
    RWMol *m = SmilesToMol(smi);
    std::string smi1=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi1=="C[C@H]1CC[C@H](F)CC1");
    delete m;
#if 0
    smi="C[C@@H]1CC[C@@H](F)CC1";
    m = SmilesToMol(smi);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi2<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi2==smi1);
    delete m;
#endif
  }

  {
    std::string smi = "F[C@H]1CC[C@H](C)CC1";
    RWMol *m = SmilesToMol(smi);
    std::string smi1=MolToSmiles(*m,true);
    delete m;
#if 0
    smi="F[C@@H]1CC[C@@H](C)CC1";
    m = SmilesToMol(smi);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi2<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi2==smi1);
    delete m;
#endif
  }

  {
    std::string smi = "F[C@H]1CC[C@](C)(C)CC1";
    RWMol *m = SmilesToMol(smi);
    std::string smi1=MolToSmiles(*m,true);
    delete m;
    smi="FC1CCC(C)(C)CC1";
    m = SmilesToMol(smi);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi2<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi2==smi1);
    delete m;
  }

  {
    std::string smi = "C1C[C@H]2CC[C@@H]1CC2";
    RWMol *m = SmilesToMol(smi);
    std::string smi1=MolToSmiles(*m,true);
    delete m;
    smi="C1CC2CCC1CC2";
    m = SmilesToMol(smi);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi2<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi2==smi1);
    delete m;
  }

  {
    std::string smi = "C[C@]12CC[C@](C)(CC1)CC2";
    RWMol *m = SmilesToMol(smi);
    std::string smi1=MolToSmiles(*m,true);
    delete m;
    smi="CC12CCC(C)(CC1)CC2";
    m = SmilesToMol(smi);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi2<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi2==smi1);
    delete m;
  }

  {
    // make sure we aren't removing stereochem that should still be there
    std::string smi = "C[C@@]12CC[C@@](C)(NC1)OC2";
    RWMol *m = SmilesToMol(smi);
    std::string smi1=MolToSmiles(*m,true);
    delete m;
    smi="CC12CCC(C)(NC1)OC2";
    m = SmilesToMol(smi);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi2<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi2!=smi1);
    delete m;
  }

#if 0
  // FIX : these tests do not pass
  {
    std::string smi = "C[C@H]1CC[C@H](C)CC1";
    RWMol *m = SmilesToMol(smi);
    std::string smi1=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi1==smi);
    delete m;
  }

  {
    std::string smi = "C1[C@@H](C)CC[C@H](C)C1";
    RWMol *m = SmilesToMol(smi);
    m->debugMol(std::cerr);
    std::string smi1=MolToSmiles(*m,true);
    smi = "C[C@H]1CC[C@H](C)CC1";
    BOOST_LOG(rdInfoLog)<<" : "<<smi<<" "<<smi1<<std::endl;
    m->debugMol(std::cerr);
    TEST_ASSERT(smi1==smi);
    delete m;
  }
#endif


  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testChiralityFrom3D(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "chirality perception from 3D coordinates: " << std::endl;

  std::string rdbase = getenv("RDBASE");
  RWMol *m;
  std::string fName,smi;
  std::string cip;

  fName = rdbase+"/Code/GraphMol/test_data/chi3d_r1.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==5);

  MolOps::assignChiralTypesFrom3D(*m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");


  delete m;
  fName = rdbase+"/Code/GraphMol/test_data/chi3d_s1.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==5);

  MolOps::assignChiralTypesFrom3D(*m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");


  delete m;
  fName = rdbase+"/Code/GraphMol/test_data/chi3d_r2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);

  MolOps::assignChiralTypesFrom3D(*m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");


  delete m;
  fName = rdbase+"/Code/GraphMol/test_data/chi3d_s2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);

  MolOps::assignChiralTypesFrom3D(*m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete m;
  fName = rdbase+"/Code/GraphMol/test_data/chi3d_r1_bad.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==5);

  // this molecule starts out with incorrect stereochemistry (e.g. the bond wedging
  // does not match the 3D structure. Start by verifying that the start position is bad:
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");
  // now assign the stereochem based on the 3D structure and check that we get it
  // right:
  MolOps::assignChiralTypesFrom3D(*m,-1,true);
  MolOps::assignStereochemistry(*m,true,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


void testIterativeChirality(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "iterative chirality (sf.net issue 1931470): " << std::endl;

  std::string rdbase = getenv("RDBASE");

  // unless otherwise noted, the R/S and Z/E assignments here 
  // match Marvin and ChemDraw.
#if 1  
  {  // atom-chirality -> atom-chirality
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi1a.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==9);
 
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");

    TEST_ASSERT(m->getAtomWithIdx(5)->hasProp("_CIPCode"));
    m->getAtomWithIdx(5)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
    m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    std::string smi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
    
    delete m;
  }

  { // atom-chirality -> atom-chirality
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi1b.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==9);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");

    TEST_ASSERT(m->getAtomWithIdx(5)->hasProp("_CIPCode"));
    m->getAtomWithIdx(5)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
    m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");

    std::string smi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
    
    delete m;
  }

  { // atom-chirality -> atom-chirality
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi1c.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==9);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");

    TEST_ASSERT(m->getAtomWithIdx(5)->hasProp("_CIPCode"));
    m->getAtomWithIdx(5)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");

    TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp("_CIPCode"));

#if 1 // this fails due to sf.net bug 1896935
    std::string smi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
#endif
    
    delete m;
  }

  { // atom-chirality -> atom-chirality
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi1d.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==9);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");

    TEST_ASSERT(m->getAtomWithIdx(5)->hasProp("_CIPCode"));
    m->getAtomWithIdx(5)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");

    TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp("_CIPCode"));

#if 1 // this fails due to sf.net bug 1896935
    std::string smi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
#endif    
    delete m;
  }

  { // atom-chirality -> atom-chirality
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi1e.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==9);

    TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp("_CIPCode"));
    TEST_ASSERT(!m->getAtomWithIdx(2)->hasProp("_CIPCode"));
    TEST_ASSERT(!m->getAtomWithIdx(4)->hasProp("_CIPCode"));

#if 0 // this fails due to sf.net bug 1896935
    std::cerr<<"m pre -----"<<std::endl;
    m->debugMol(std::cerr);
    std::cerr<<"-----"<<std::endl;
    std::string smi1=MolToSmiles(*m,true);
    std::cerr<<"m post -----"<<std::endl;
    m->debugMol(std::cerr);
    std::cerr<<"-----"<<std::endl;
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::cerr<<"m2 pre -----"<<std::endl;
    m->debugMol(std::cerr);
    std::cerr<<"-----"<<std::endl;
    std::string smi2=MolToSmiles(*m,true);
    std::cerr<<"m post -----"<<std::endl;
    m->debugMol(std::cerr);
    std::cerr<<"-----"<<std::endl;
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
#endif    
    delete m;
  }

  { // bond-stereochem -> atom-chirality
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi2a.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==8);

    TEST_ASSERT(m->getBondBetweenAtoms(2,5)->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(0,4)->getStereo()==Bond::STEREOZ);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");  // this value is from ChemDraw, Marvin doesn't tag it.

    std::string smi1=MolToSmiles(*m,true);

    MolOps::removeStereochemistry(*m);
    TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    TEST_ASSERT(m->getBondBetweenAtoms(0,4)->getStereo()==Bond::STEREONONE);
    TEST_ASSERT(m->getBondBetweenAtoms(2,5)->getStereo()==Bond::STEREONONE);
    
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);

    delete m;
  }

  {  // bond-stereochem -> atom-chirality
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi2b.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==8);

    TEST_ASSERT(m->getBondBetweenAtoms(2,5)->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(0,4)->getStereo()==Bond::STEREOZ);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");  // this value is from ChemDraw, Marvin doesn't tag it.

    std::string smi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
    
    delete m;
  }

  {  // bond-stereochem -> atom-chirality
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi2c.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==8);

    TEST_ASSERT(m->getBondBetweenAtoms(2,5)->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(0,4)->getStereo()==Bond::STEREOE);

    TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp("_CIPCode"));

    std::string smi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
    
    delete m;
  }

  {  // bond-stereochem -> atom-chirality
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi2d.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==8);

    TEST_ASSERT(m->getBondBetweenAtoms(2,5)->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(0,4)->getStereo()==Bond::STEREOANY);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");  // this value is from ChemDraw, Marvin doesn't tag it.

    std::string smi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
    
    delete m;
  }

  {  // bond-stereochem -> atom-chirality
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi2e.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==8);

    TEST_ASSERT(m->getBondBetweenAtoms(2,5)->getStereo()==Bond::STEREOANY);
    TEST_ASSERT(m->getBondBetweenAtoms(0,4)->getStereo()==Bond::STEREOZ);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");  // this value is from ChemDraw, Marvin doesn't tag it.

    std::string smi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
    
    delete m;
  }

  {  // atom chirality -> bond stereochemistry
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi3a.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==11);

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
    m->getAtomWithIdx(3)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");  

    TEST_ASSERT(m->getAtomWithIdx(7)->hasProp("_CIPCode"));
    m->getAtomWithIdx(7)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");  


    TEST_ASSERT(m->getBondBetweenAtoms(1,2)->getStereo()==Bond::STEREOZ); // this value is from ChemDraw, Marvin doesn't tag it.

    std::string smi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
    
    delete m;
  }

  {  // atom chirality -> bond stereochemistry
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi3b.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==11);

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
    m->getAtomWithIdx(3)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");  

    TEST_ASSERT(m->getAtomWithIdx(7)->hasProp("_CIPCode"));
    m->getAtomWithIdx(7)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");  


    TEST_ASSERT(m->getBondBetweenAtoms(1,2)->getStereo()==Bond::STEREOE); // this value is from ChemDraw, Marvin doesn't tag it.

    std::string smi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
    
    delete m;
  }

  {  // atom chirality -> bond stereochemistry
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi3c.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==11);

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
    m->getAtomWithIdx(3)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");  

    TEST_ASSERT(m->getAtomWithIdx(7)->hasProp("_CIPCode"));
    m->getAtomWithIdx(7)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");  


    TEST_ASSERT(m->getBondBetweenAtoms(1,2)->getStereo()==Bond::STEREONONE);

    std::string smi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
    
    delete m;
  }
#endif

  {  // bond stereochemistry -> bond stereochemistry
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi4a.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==9);

    TEST_ASSERT(m->getBondBetweenAtoms(4,5)->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(3,7)->getStereo()==Bond::STEREOZ);

    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getStereo()==Bond::STEREOE);
    
    std::string smi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
    
    delete m;
  }

  {  // bond stereochemistry -> bond stereochemistry
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi4b.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==9);

    TEST_ASSERT(m->getBondBetweenAtoms(4,5)->getStereo()==Bond::STEREOZ);
    TEST_ASSERT(m->getBondBetweenAtoms(3,7)->getStereo()==Bond::STEREOE);

    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getStereo()==Bond::STEREOZ);
    
    std::string smi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
    
    delete m;
  }

  {  // bond stereochemistry -> bond stereochemistry
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi4c.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==9);

    TEST_ASSERT(m->getBondBetweenAtoms(4,5)->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(3,7)->getStereo()==Bond::STEREOE);

    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getStereo()==Bond::STEREONONE);
    
    std::string smi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
    
    delete m;
  }

  {  // bond stereochemistry -> bond stereochemistry
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi4d.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==9);

    TEST_ASSERT(m->getBondBetweenAtoms(4,5)->getStereo()==Bond::STEREOZ);
    TEST_ASSERT(m->getBondBetweenAtoms(3,7)->getStereo()==Bond::STEREOZ);

    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getStereo()==Bond::STEREONONE);
    
    std::string smi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
    
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


void testBondDirRemoval(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testing that the removal of bond directions is correct: " << std::endl;

  std::string rdbase = getenv("RDBASE");


  {  
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/stereoOrder1.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==7);
 
    TEST_ASSERT(m->getBondBetweenAtoms(1,2)->getStereo()==Bond::STEREOZ);
    TEST_ASSERT(m->getBondBetweenAtoms(4,5)->getStereo()==Bond::STEREOE);

    // on input all the single bonds are in the same direction:
    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getBondDir()==m->getBondBetweenAtoms(1,4)->getBondDir());
    TEST_ASSERT(m->getBondBetweenAtoms(2,3)->getBondDir()==m->getBondBetweenAtoms(1,4)->getBondDir());
    TEST_ASSERT(m->getBondBetweenAtoms(5,6)->getBondDir()==m->getBondBetweenAtoms(1,4)->getBondDir());

    std::string smi1=MolToSmiles(*m,true);

    // check removal of redundant bond direction information:
    std::vector<int> ranks(m->getNumAtoms(),0);
    MolOps::rankAtoms(*m,ranks);
    std::vector<Canon::AtomColors> colors(m->getNumAtoms());
    Canon::MolStack stack;
    Canon::canonicalizeFragment(*m,0,colors,ranks,stack);

    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getBondDir()==Bond::NONE);
    TEST_ASSERT(m->getBondBetweenAtoms(2,3)->getBondDir()==m->getBondBetweenAtoms(1,4)->getBondDir());
    TEST_ASSERT(m->getBondBetweenAtoms(5,6)->getBondDir()==m->getBondBetweenAtoms(1,4)->getBondDir());

    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);

    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
    
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


void testIssue2705543(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Issue 2705543: " << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;
  RWMol *m;
  std::string cip;

  {
    fName = rdbase+"/Code/GraphMol/test_data/Issue2705543.1h.mol";
    m = MolFileToMol(fName,true,false);
    TEST_ASSERT(m);

    MolOps::assignChiralTypesFrom3D(*m);
    MolOps::assignStereochemistry(*m,true);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
    m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");

    TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
    m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
    m->getAtomWithIdx(3)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(4)->hasProp("_CIPCode"));
    m->getAtomWithIdx(4)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");

    delete m;
  }
  {
    fName = rdbase+"/Code/GraphMol/test_data/Issue2705543.2h.mol";
    m = MolFileToMol(fName,true,false);
    TEST_ASSERT(m);

    MolOps::assignChiralTypesFrom3D(*m);
    MolOps::assignStereochemistry(*m,true);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
    m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");

    TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
    m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
    m->getAtomWithIdx(3)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(4)->hasProp("_CIPCode"));
    m->getAtomWithIdx(4)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");

    delete m;
  }
  {
    fName = rdbase+"/Code/GraphMol/test_data/Issue2705543.1.mol";
    m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==13);

    MolOps::assignChiralTypesFrom3D(*m);
    MolOps::assignStereochemistry(*m,true);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
    m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");

    TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
    TEST_ASSERT(m->getAtomWithIdx(2)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
    m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
    m->getAtomWithIdx(3)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(4)->hasProp("_CIPCode"));
    TEST_ASSERT(m->getAtomWithIdx(4)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
    m->getAtomWithIdx(4)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");

    delete m;
  }
  {
    fName = rdbase+"/Code/GraphMol/test_data/Issue2705543.2.mol";
    m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==13);

    MolOps::assignChiralTypesFrom3D(*m);
    MolOps::assignStereochemistry(*m,true);

#if 0
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      if(m->getAtomWithIdx(i)->hasProp("_CIPCode")){
        m->getAtomWithIdx(i)->getProp("_CIPCode",cip);
        std::cerr<<"  >> "<<i<<" "<<cip<<std::endl;
      }
    }
#endif    
    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
    m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");

    TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
    TEST_ASSERT(m->getAtomWithIdx(2)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
    m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
    m->getAtomWithIdx(3)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(4)->hasProp("_CIPCode"));
    TEST_ASSERT(m->getAtomWithIdx(4)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
    m->getAtomWithIdx(4)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");

    delete m;
  }


  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


void testIssue2762917(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Issue 2762917: chirality swap on addHs()" << std::endl;

  std::string rdbase = getenv("RDBASE");

  {
    RWMol *m;
    std::string cip;
    std::string smiles="[C@@H](C)(Cl)O";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);

    MolOps::assignStereochemistry(*m,true);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
    m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
    MolOps::assignStereochemistry(*m,true);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
    m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    delete m;
  }

  {
    RWMol *m;
    std::string cip;
    std::string smiles="CCC.[C@@H](C)(Cl)O";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);

    MolOps::assignStereochemistry(*m,true);

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
    m->getAtomWithIdx(3)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    MolOps::addHs(*m);

    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
    MolOps::assignStereochemistry(*m,true);

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
    m->getAtomWithIdx(3)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    delete m;
  }

  {
    RWMol *m;
    std::string cip;
    std::string smiles="[C@@H]([C@H](C)O)(C)O";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);

    MolOps::assignStereochemistry(*m,true);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
    m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    MolOps::addHs(*m);

    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);

    MolOps::assignStereochemistry(*m,true);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
    m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    delete m;
  }


  {
    RWMol *m;
    std::string cip;
    std::string smiles="C1CC.[C@@H]1(Cl)O";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);

    MolOps::assignStereochemistry(*m,true);

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
    m->getAtomWithIdx(3)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    MolOps::addHs(*m);

    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
    MolOps::assignStereochemistry(*m,true);

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
    m->getAtomWithIdx(3)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


void testIssue3009911(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Issue 3009911: bad atom priorities" << std::endl;
  
  {
    RWMol *m;
    std::string smiles="F[C@](O)(c1ccccc1)C(=C)CO";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int *ranks;
    ranks = new int[m->getNumAtoms()];
    MolOps::assignStereochemistry(*m,true);
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      int rank;
      TEST_ASSERT(m->getAtomWithIdx(i)->hasProp("_CIPRank"))
      m->getAtomWithIdx(i)->getProp("_CIPRank",rank);
      ranks[i]=rank;
    }
    // basics:
    TEST_ASSERT(ranks[0]>ranks[1]);
    TEST_ASSERT(ranks[2]>ranks[1]);
    TEST_ASSERT(ranks[0]>ranks[2]);
    // now the key point:
    TEST_ASSERT(ranks[3]<ranks[9]);

    std::string cip;
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    delete m;
    delete [] ranks;
  }
  {
    RWMol *m;
    std::string smiles="COC(C)(OC)[C@](O)(F)C(C)=O";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int *ranks;
    ranks = new int[m->getNumAtoms()];
    MolOps::assignStereochemistry(*m,true);
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      int rank;
      TEST_ASSERT(m->getAtomWithIdx(i)->hasProp("_CIPRank"))
      m->getAtomWithIdx(i)->getProp("_CIPRank",rank);
      ranks[i]=rank;
    }
    // basics:
    TEST_ASSERT(ranks[8]>ranks[7]);
    TEST_ASSERT(ranks[7]>ranks[9]);
    TEST_ASSERT(ranks[7]>ranks[2]);
    // FIX: these are the key points, but at the moment they are not handled correctly
    // due to a weakness in the CIP-ranking algorithm.
    //TEST_ASSERT(ranks[2]>ranks[9]);
    std::string cip;
    TEST_ASSERT(m->getAtomWithIdx(6)->hasProp("_CIPCode"));
    m->getAtomWithIdx(6)->getProp("_CIPCode",cip);
    //TEST_ASSERT(cip=="R");

    delete m;
    delete [] ranks;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue3139534(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Issue 3139534: stereochemistry in larger rings" << std::endl;

  // the smiles generation part of this is in SmilesParse/test.cpp


  // tests that the creation and assignment are correct:
  {
    RWMol *m;
    std::string smiles="C1COCC/C=C\\CC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(5)->getStereo()==Bond::STEREOZ);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles="C1COCC/C=C/CC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(5)->getStereo()==Bond::STEREOE);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles="C/1=C/OCCC=CCC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo()==Bond::STEREOZ);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles="C1=C/OCCC=CCC\\1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo()==Bond::STEREOZ);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles="C\\1=C/OCCC=CCC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo()==Bond::STEREOE);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles="C1=C/OCCC=CCC/1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo()==Bond::STEREOE);
    delete m;
  }

  {
    RWMol *m;
    std::string smiles="C/1=C/OCC/C=C\\CC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo()==Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(5)->getStereo()==Bond::STEREOZ);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles="C\\1=C/OCC/C=C\\CC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m->getBondWithIdx(5)->getStereo()==Bond::STEREOZ);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles="C1=C/OCC/C=C\\CC\\1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo()==Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(5)->getStereo()==Bond::STEREOZ);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles="C1=C/OCC/C=C\\CC/1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m->getBondWithIdx(5)->getStereo()==Bond::STEREOZ);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testFindChiralAtoms(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test findChiralAtoms." << std::endl;

  {
    // by default the chirality possible flag is not assigned:
    RWMol *m;
    std::string smiles="F[C@H](Cl)C(Cl)(Br)C(F)(F)F";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    MolOps::assignStereochemistry(*m,true,true);
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    TEST_ASSERT(!(m->getAtomWithIdx(3)->hasProp("_CIPCode")));
    TEST_ASSERT(!(m->getAtomWithIdx(6)->hasProp("_CIPCode")));
    TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp("_ChiralityPossible"));
    TEST_ASSERT(!m->getAtomWithIdx(3)->hasProp("_ChiralityPossible"));
    TEST_ASSERT(!(m->getAtomWithIdx(6)->hasProp("_ChiralityPossible")));

    // but we can force it:
    MolOps::assignStereochemistry(*m,true,true,true);
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    TEST_ASSERT(!(m->getAtomWithIdx(3)->hasProp("_CIPCode")));
    TEST_ASSERT(!(m->getAtomWithIdx(6)->hasProp("_CIPCode")));
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_ChiralityPossible"));
    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_ChiralityPossible"));
    TEST_ASSERT(!(m->getAtomWithIdx(6)->hasProp("_ChiralityPossible")));
    
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue3453172(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Issue 3453172: stereochemistry at three-coordinate S and Se" << std::endl;

  {
    RWMol *m;
    std::string smiles="C=[S@](F)Br";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles="C[S@+](F)Br";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles="C=[Se@](F)Br";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles="C[Se@+](F)Br";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles="C=[S@](Br)Br";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_UNSPECIFIED);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles="C[S@+](Br)Br";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_UNSPECIFIED);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles="C=[Se@](Br)Br";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_UNSPECIFIED);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles="C[Se@+](Br)Br";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_UNSPECIFIED);
    delete m;
  }

  {
    // this was issue 254
    RWMol *m;
    std::string smiles="O=[S@](c1ccccc1)C";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP3);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


void testGithub87(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github issue 87: removal of bond wedging" << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    std::string fName = rdbase+"/Code/GraphMol/test_data/github87.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==5);
    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
    WedgeMolBonds(*m,&m->getConformer());
    MolOps::assignStereochemistry(*m,true,true);
    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getBondDir()==Bond::BEGINWEDGE);
    m->getAtomWithIdx(0)->setChiralTag(Atom::CHI_UNSPECIFIED);
    MolOps::assignStereochemistry(*m,true,true);
    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getBondDir()==Bond::NONE);
    
    delete m;
  }
  {
    std::string fName = rdbase+"/Code/GraphMol/test_data/github87.2.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==5);
    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
    WedgeMolBonds(*m,&m->getConformer());
    MolOps::assignStereochemistry(*m,true,true);
    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getBondDir()==Bond::BEGINDASH);
    m->getAtomWithIdx(0)->setChiralTag(Atom::CHI_UNSPECIFIED);
    MolOps::assignStereochemistry(*m,true,true);
    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getBondDir()==Bond::NONE);
    
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub90(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github issue 90: isotopes and chirality" << std::endl;

  {
    std::string smi="C[C@@H](F)[13CH3]";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==4);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
    MolOps::assignStereochemistry(*m,true,true);

    std::string cip;
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    delete m;
  }
  {
    std::string smi="[13CH3][C@@H](F)C";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==4);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
    MolOps::assignStereochemistry(*m,true,true);

    std::string cip;
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");
    delete m;
  }
  {
    std::string smi="[CH3][C@@H](F)C";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==4);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_UNSPECIFIED);
    MolOps::assignStereochemistry(*m,true,true);
    TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    delete m;
  }
  {
    std::string smi="C\\C([13CH3])=C(/C)[13CH3]";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==6);
    TEST_ASSERT(m->getBondWithIdx(2)->getStereo()==Bond::STEREOZ);
    delete m;
  }
  {
    std::string smi="C\\C([CH3])=C(/C)[13CH3]";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==6);
    TEST_ASSERT(m->getBondWithIdx(2)->getStereo()==Bond::STEREONONE);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}



int main(){
  RDLog::InitLogs();
  //boost::logging::enable_logs("rdApp.debug");

#if 1
  testSmiles1();
  testMol1();
  testMol2();
  testRoundTrip();
  testChiralityCleanup();
  testChiralityFrom3D();
  testIterativeChirality();
  testBondDirRemoval();
  testIssue2705543();
  testIssue2762917();
  testIssue3009911();
  testIssue3139534();
  testFindChiralAtoms();
  testIssue3453172();
  testRingStereochemistry();
  testGithub87();
#endif
  testGithub90();
  return 0;
}


