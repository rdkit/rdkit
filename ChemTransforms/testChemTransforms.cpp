// $Id$
//
//  Copyright (C) 2006 Greg Landrum
//
//   @@ All Rights Reserved  @@
//
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <string>
#include <iostream>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>

using namespace RDKit;

void testDeleteSubstruct() 
{
  int i = 0;
  ROMol *mol1=0,*mol2=0,*matcher1=0,*matcher2=0,*matcher3=0;
  std::string smi,sma;
  
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing deleteSubstruct" << std::endl;

  // a lot of the seemingly repetitive stuff is here for Issue96
  smi = "CCC(=O).C=O";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  TEST_ASSERT(mol1->getNumAtoms()==6)
  sma = "C=O";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  
  mol2 = deleteSubstructs(*mol1,*matcher1,0);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==2)
  mol2 = deleteSubstructs(*mol2,*matcher1,0);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==2)
  
  delete matcher1;
  sma = "[Cl;H1&X1,-]";
  matcher1 = SmartsToMol(sma);
  sma = "[Na+]";
  matcher2 = SmartsToMol(sma);
  sma = "[O;H2,H1&-,X0&-2]";
  matcher3 = SmartsToMol(sma);
  delete mol1;
  mol1 = SmilesToMol("CCO.Cl");
  TEST_ASSERT(mol1);
  TEST_ASSERT(mol1->getNumAtoms()==4);

  delete mol2;
  mol2 = deleteSubstructs(*mol1,*matcher1,true);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==3);
  mol2 = deleteSubstructs(*mol2,*matcher2,true);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==3);
  mol2 = deleteSubstructs(*mol2,*matcher3,true);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==3);

  delete mol1;
  mol1 = SmilesToMol("CC(=O)[O-].[Na+]");
  TEST_ASSERT(mol1);
  TEST_ASSERT(mol1->getNumAtoms()==5);

  delete matcher1;
  matcher1 = SmartsToMol("[Cl;H1&X1,-]");
  delete matcher2;
  matcher2 = SmartsToMol("[Na+]");

  mol2 = deleteSubstructs(*mol1,*matcher1,true);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==5);

  mol2 = deleteSubstructs(*mol2,*matcher2,true);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==4);

  mol2 = deleteSubstructs(*mol2,*matcher1,true);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==4);

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testReplaceSubstructs() 
{
  int i = 0;
  ROMol *mol1=0,*mol2=0,*matcher1=0,*frag=0;
  std::string smi,sma;
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

  vect = replaceSubstructs(*mol1,*matcher1,*frag);
  TEST_ASSERT(vect.size()==1);
  TEST_ASSERT(mol1->getNumAtoms()==4);
  TEST_ASSERT(vect[0]->getNumAtoms()==4);
  

  delete mol1;
  smi = "CCCC(=O)O";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  vect = replaceSubstructs(*mol1,*matcher1,*frag);
  TEST_ASSERT(vect.size()==1);
  TEST_ASSERT(mol1->getNumAtoms()==6);
  TEST_ASSERT(vect[0]->getNumAtoms()==4);
  
  vect = replaceSubstructs(*mol1,*matcher1,*frag,true);
  TEST_ASSERT(vect.size()==1);
  TEST_ASSERT(mol1->getNumAtoms()==6);
  TEST_ASSERT(vect[0]->getNumAtoms()==4);
  
  delete mol1;
  smi = "OC(=O)CCCC(=O)O";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  vect = replaceSubstructs(*mol1,*matcher1,*frag);
  TEST_ASSERT(vect.size()==2);
  TEST_ASSERT(mol1->getNumAtoms()==9);
  TEST_ASSERT(vect[0]->getNumAtoms()==7);
  TEST_ASSERT(vect[1]->getNumAtoms()==7);

  // use replaceAll when a single replacement is available:
  vect = replaceSubstructs(*mol1,*matcher1,*frag,true);
  TEST_ASSERT(vect.size()==1);
  TEST_ASSERT(mol1->getNumAtoms()==9);
  TEST_ASSERT(vect[0]->getNumAtoms()==5);
  
  delete matcher1;
  sma = "C=O";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  delete mol1;
  smi = "CC(=O)C";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  vect = replaceSubstructs(*mol1,*matcher1,*frag);
  TEST_ASSERT(vect.size()==1);
  TEST_ASSERT(mol1->getNumAtoms()==4);
  TEST_ASSERT(vect[0]->getNumAtoms()==3);

  delete mol1;
  smi = "C1C(=O)C1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  vect = replaceSubstructs(*mol1,*matcher1,*frag);
  TEST_ASSERT(vect.size()==1);
  TEST_ASSERT(mol1->getNumAtoms()==4);
  TEST_ASSERT(vect[0]->getNumAtoms()==3);
  TEST_ASSERT(vect[0]->getNumBonds()==3);

  
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}




int main() { 
  RDLog::InitLogs();
    
  BOOST_LOG(rdInfoLog) << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing Chemical Transforms \n";

#if 1
  testDeleteSubstruct();
  testReplaceSubstructs();
#endif

  BOOST_LOG(rdInfoLog) << "*******************************************************\n";
  return(0);
}

