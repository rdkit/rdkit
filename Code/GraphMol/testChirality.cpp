//  $Id$
// 
//   Copyright (C) 2007, 2008 Greg Landrum
//
//   @@ All Rights Reserved  @@
//
//
// There are chirality test cases spread all over the place. Many of the
// tests here are repeats, but it's good to have everything in one place.
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
//#include <boost/log/functions.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>

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
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  smi="[C@H](O)(N)I";
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  BOOST_LOG(rdInfoLog) << " >>>>>>>>>>>>> mol file <<<<<<<<<<<<<< " << std::endl;
  delete m;
  fName = rdbase+"/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir1a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete m;
  fName = rdbase+"/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir1b.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete m;
  fName = rdbase+"/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir2a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete m;
  fName = rdbase+"/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir2b.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

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
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  smi=MolToSmiles(*m,true);
  delete m;
  m=SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
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
  MolOps::assignAtomChiralCodes(*m);
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
  MolOps::assignAtomChiralCodes(*m);
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
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
#if 1
  smi=MolToSmiles(*m,true);
  delete m;
  m=SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
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
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
#if 1
  smi=MolToSmiles(*m,true);
  delete m;
  m=SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
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
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
#if 1
  smi=MolToSmiles(*m,true);
  delete m;
  m=SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
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
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
#if 1
  smi=MolToSmiles(*m,true);
  delete m;
  m=SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
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
  MolOps::assignAtomChiralCodes(*m);
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
  MolOps::assignAtomChiralCodes(*m);
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
  MolOps::assignAtomChiralCodes(*m);
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
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "F[C@](Br)(I)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_UNSPECIFIED);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "F[C@](I)(Cl)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "Cl[C@](Br)(F)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "Cl[C@](F)(I)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "I[C@](F)(Br)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "I[C@](Br)(Cl)F";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "F[C@@](Br)(Cl)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "F[C@@](Cl)(I)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "Cl[C@@](Br)(I)F";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "Cl[C@@](F)(Br)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "[C@@](Cl)(F)(Br)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "F[C@H](Cl)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "Br[C@H](F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "Br[C@]([H])(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "Br[C@](F)(Cl)[H]";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "Br[C@]1(F)(Cl).[H]1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "Br[C@H]1Cl.F1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "Br[C@]12Cl.F2.[H]1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "Br[C@]21Cl.F1.[H]2";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "Br[C@]12Cl.F1.[H]2";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete mol;
  smi = "[C@@](C)(Br)(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "[C@@]([H])(Br)(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "[C@@H](Br)(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp("_CIPCode"));
  mol->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete mol;
  smi = "[H][C@@](Br)(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignAtomChiralCodes(*mol);
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
  MolOps::assignAtomChiralCodes(*mol,true);
  TEST_ASSERT(!mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_UNSPECIFIED);
  delete mol;

  smi = "F[C@+](C)(Cl)(Br)I";
  mol = SmilesToMol(smi,false,false);
  mol2 = MolOps::removeHs(*mol,false,false);
  delete mol;
  mol=mol2;
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignAtomChiralCodes(*mol,true);
  TEST_ASSERT(!mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_UNSPECIFIED);
  delete mol;

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
  TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
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
  TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
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
    delete m;
  }

  {  // bond-stereochem -> atom-chirality
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi2d.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==8);

    TEST_ASSERT(m->getBondBetweenAtoms(2,5)->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(0,4)->getStereo()==Bond::STEREONONE);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");  // this value is from ChemDraw, Marvin doesn't tag it.

    delete m;
  }

  {  // bond-stereochem -> atom-chirality
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi2e.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==8);

    TEST_ASSERT(m->getBondBetweenAtoms(2,5)->getStereo()==Bond::STEREONONE);
    TEST_ASSERT(m->getBondBetweenAtoms(0,4)->getStereo()==Bond::STEREOZ);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
    m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");  // this value is from ChemDraw, Marvin doesn't tag it.

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
    delete m;
  }


  {  // bond stereochemistry -> bond stereochemistry
    std::string cip;

    std::string fName = rdbase+"/Code/GraphMol/test_data/iChi4a.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==9);

    TEST_ASSERT(m->getBondBetweenAtoms(4,5)->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(3,7)->getStereo()==Bond::STEREOZ);

    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getStereo()==Bond::STEREOE);
    
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

    //TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getStereo()==Bond::STEREOZ);
    
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
#endif
  testIterativeChirality();
  
  return 0;
}


