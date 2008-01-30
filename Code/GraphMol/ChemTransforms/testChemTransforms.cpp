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
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/FileParsers/FileParsers.h>

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


void testReplaceSidechains() 
{
  int i = 0;
  ROMol *mol1=0,*mol2=0,*matcher1=0;
  std::string smi,sma;
  
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing replaceSidechains" << std::endl;

  smi = "ClC1CC(F)C1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);

  sma = "C1CCC1";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  mol2 = replaceSidechains(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==6);
  
  delete mol1;
  smi = "ClC1C(F)C1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  mol2 = replaceSidechains(*mol1,*matcher1);
  TEST_ASSERT(!mol2);

  delete matcher1;  
  sma = "C=O";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  delete mol1;
  smi = "CC=O";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);

  delete mol2;
  mol2 = replaceSidechains(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==3);
  smi = MolToSmiles(*mol2);
  TEST_ASSERT(smi=="[Xa]C=O");
  
  delete mol1;
  smi = "CC(C)=O";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);

  delete mol2;
  mol2 = replaceSidechains(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==4);
  smi = MolToSmiles(*mol2);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi=="[Xa]C([Xb])=O"||smi=="[Xb]C([Xa])=O");

  delete mol1;
  delete mol2;
  delete matcher1;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testReplaceCore() 
{
  int i = 0;
  ROMol *mol1=0,*mol2=0,*matcher1=0;
  std::string smi,sma;
  
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing replaceCore" << std::endl;

  smi = "ClC1CC(F)C1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);

  sma = "C1CCC1";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==4);
  smi = MolToSmiles(*mol2);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi=="[Xa]F.Cl[Xb]"||smi=="[Xb]F.Cl[Xa]");

  delete mol1;
  smi = "CCC=O";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);

  delete matcher1;
  sma = "C=O";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==3);
  smi = MolToSmiles(*mol2);
  TEST_ASSERT(smi=="[Xa]CC");
  
  delete mol1;
  smi = "C1C(=O)CC1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);

  delete matcher1;
  sma = "C=O";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==5);
  smi = MolToSmiles(*mol2);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi=="[Xa]CCC[Xb]"||smi=="[Xb]CCC[Xa]");


  delete mol1;
  smi = "CNC";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);

  delete matcher1;
  sma = "N";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==4);
  smi = MolToSmiles(*mol2);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi=="[Xa]C.[Xb]C"||smi=="[Xb]C.[Xa]C");

  delete mol1;
  smi = "OC1CCC1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);

  delete matcher1;
  sma = "[CH2][CH2][CH2]";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  smi = MolToSmiles(*mol2);
  TEST_ASSERT(mol2->getNumAtoms()==4);
  smi = MolToSmiles(*mol2);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi=="[Xb]C([Xa])O"||smi=="[Xa]C([Xb])O");

  delete mol1;
  smi = "C/C=C/CN/C=C/O";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);

  delete matcher1;
  sma = "N";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==9);
  smi = MolToSmiles(*mol2,true);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi=="[Xa]C/C=C/C.O/C=C/[Xb]"||smi=="[Xb]C/C=C/C.O/C=C/[Xa]"||
              smi=="[Xa]/C=C/O.C/C=C/C[Xb]"||smi=="[Xb]/C=C/O.C/C=C/C[Xa]");

  delete mol1;
  smi = "C[C@](F)(Cl)N";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete matcher1;
  sma = "N";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==5);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[Xa][C@@](C)(F)Cl");
  delete mol1;
  smi = "C[C@](F)(N)Cl";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==5);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[Xa][C@](C)(F)Cl");
  delete mol1;
  smi = "N[C@](C)(F)Cl";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==5);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[Xa][C@](C)(F)Cl");

  delete mol1;
  smi = "C[C@@](F)(Cl)N";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete matcher1;
  sma = "N";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==5);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[Xa][C@](C)(F)Cl");
  delete mol1;
  smi = "C[C@@](F)(N)Cl";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==5);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[Xa][C@@](C)(F)Cl");
  smi = "C[C@@](N)(Cl)F";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==5);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[Xa][C@@](C)(F)Cl");


  delete mol1;
  smi = "C[C@H](F)N";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete matcher1;
  sma = "N";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==4);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[Xa][C@H](C)F");
  delete mol1;
  smi = "C[C@H](N)F";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==4);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[Xa][C@@H](C)F");

  delete mol1;
  smi = "N[C@H](C)F";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==4);
  smi = MolToSmiles(*mol2,true);
  BOOST_LOG(rdInfoLog)<<"smiles: "<<smi <<std::endl;
  TEST_ASSERT(smi=="[Xa][C@H](C)F");
  delete mol1;
  smi = "F[C@H](C)N";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==4);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[Xa][C@@H](C)F");
  
  delete mol1;
  smi = "N[C@H]1CCCO1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==6);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[Xa][C@@H]1OCCC1");

  delete mol1;
  smi = "ClC1CC(F)C1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);

  delete matcher1;
  sma = "[*]C1CC([*])C1";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  mol2 = replaceCore(*mol1,*matcher1,false);
  TEST_ASSERT(mol2);
  smi = MolToSmiles(*mol2);
  BOOST_LOG(rdInfoLog)<<"smiles: "<<smi <<std::endl;
  TEST_ASSERT(mol2->getNumAtoms()==4);
  smi = MolToSmiles(*mol2);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi=="[Xa]F.Cl[Xb]"||smi=="[Xb]F.Cl[Xa]");

  delete mol1;
  delete mol2;
  delete matcher1;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testReplaceCoreCrash() 
{
  int i = 0;
  ROMol *mol1=0,*mol2=0,*matcher1=0;
  std::string smi,sma;
  
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing a former crash in replaceCore" << std::endl;

  mol1 = MolFileToMol("testData/oldcrash.mol");
  TEST_ASSERT(mol1);

  sma = "N";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==6);
  smi = MolToSmiles(*mol2);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi=="[Xa]CC.[Xb]CC"||smi=="[Xb]CC.[Xa]CC");

}

int main() { 
  RDLog::InitLogs();
    
  BOOST_LOG(rdInfoLog) << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing Chemical Transforms \n";

#if 1
  testDeleteSubstruct();
  testReplaceSubstructs();
  testReplaceSidechains();
#endif
  testReplaceCore();
  testReplaceCoreCrash();

  BOOST_LOG(rdInfoLog) << "*******************************************************\n";
  return(0);
}

