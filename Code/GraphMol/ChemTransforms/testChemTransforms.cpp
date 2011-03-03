// $Id$
//
//  Copyright (C) 2006-2011 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
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
  ROMol *mol1=0,*matcher1=0,*frag=0;
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
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[1*]C=O");
  
  delete mol1;
  smi = "CC(C)=O";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);

  delete mol2;
  mol2 = replaceSidechains(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==4);
  smi = MolToSmiles(*mol2,true);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi=="[1*]C([2*])=O"||smi=="[2*]C([1*])=O");

  delete mol1;
  delete mol2;
  delete matcher1;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testReplaceCore() 
{
  ROMol *mol1=0,*mol2=0,*matcher1=0;
  std::string smi,sma;
  
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing replaceCore" << std::endl;

  smi = "ClC1CC(F)C1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  sma = "C1CC1";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(!mol2);

  smi = "ClC1CC(F)C1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);

  sma = "C1CCC1";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==4);
  smi = MolToSmiles(*mol2,true);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi=="[1*]Cl.[2*]F"||smi=="[2*]Cl.[1*]F");

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
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[1*]CC");
  
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
  smi = MolToSmiles(*mol2,true);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi=="[1*]CCC[2*]"||smi=="[1*]CCC[2*]");


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
  smi = MolToSmiles(*mol2,true);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi=="[1*]C.[2*]C"||smi=="[2*]C.[1*]C");

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
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(mol2->getNumAtoms()==4);
  smi = MolToSmiles(*mol2,true);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi=="[2*]C([1*])O"||smi=="[1*]C([2*])O");

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
  TEST_ASSERT(smi=="[1*]C/C=C/C.[2*]/C=C/O"||smi=="[1*]/C=C/O.[2*]C/C=C/C");

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
  TEST_ASSERT(smi=="[1*][C@@](C)(F)Cl");
  delete mol1;
  smi = "C[C@](F)(N)Cl";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==5);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[1*][C@](C)(F)Cl");
  delete mol1;
  smi = "N[C@](C)(F)Cl";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==5);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[1*][C@](C)(F)Cl");

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
  TEST_ASSERT(smi=="[1*][C@](C)(F)Cl");
  delete mol1;
  smi = "C[C@@](F)(N)Cl";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==5);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[1*][C@@](C)(F)Cl");
  smi = "C[C@@](N)(Cl)F";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==5);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[1*][C@@](C)(F)Cl");


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
  TEST_ASSERT(smi=="[1*][C@H](C)F");
  delete mol1;
  smi = "C[C@H](N)F";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==4);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[1*][C@@H](C)F");

  delete mol1;
  smi = "N[C@H](C)F";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==4);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[1*][C@H](C)F");
  delete mol1;
  smi = "F[C@H](C)N";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==4);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[1*][C@@H](C)F");
  
  delete mol1;
  smi = "N[C@H]1CCCO1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  delete mol2;
  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==6);
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(smi=="[1*][C@@H]1OCCC1");

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
  smi = MolToSmiles(*mol2,true);
  TEST_ASSERT(mol2->getNumAtoms()==4);
  smi = MolToSmiles(*mol2,true);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi=="[1*]F.[2*]Cl"||smi=="[1*]Cl.[2*]F");

  delete mol1;
  delete mol2;
  delete matcher1;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testReplaceCoreLabels() 
{
  
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing replaceCore with labels" << std::endl;

  {
    std::string sma = "n1cocc1";
    ROMol *matcher = SmartsToMol(sma);
    TEST_ASSERT(matcher);

    std::string smi = "n1c(CC)oc(C)c1";
    ROMol *mol1 = SmilesToMol(smi);
    TEST_ASSERT(mol1);

    ROMol *mol2 = replaceCore(*mol1,*matcher,true,true);
    TEST_ASSERT(mol2);
    TEST_ASSERT(mol2->getNumAtoms()==5);
    smi = MolToSmiles(*mol2,true);
    TEST_ASSERT(smi=="[1*]CC.[3*]C");

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

    ROMol *mol2 = replaceCore(*mol1,*matcher,true,true);
    TEST_ASSERT(mol2);
    smi = MolToSmiles(*mol2,true);
    TEST_ASSERT(smi=="[1*]C.[3*]CC");

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

    ROMol *mol2 = replaceCore(*mol1,*matcher,true,true);
    TEST_ASSERT(mol2);
    smi = MolToSmiles(*mol2,true);
    TEST_ASSERT(smi=="[1*]C.[3*]CC.[4*]OC");

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

    ROMol *mol2 = replaceCore(*mol1,*matcher,true,true);
    TEST_ASSERT(mol2);
    smi = MolToSmiles(*mol2,true);
    TEST_ASSERT(smi=="[4*]OC");

    delete mol1;
    delete mol2;
    delete matcher;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}


void testReplaceCoreCrash() 
{
  ROMol *mol1=0,*mol2=0,*matcher1=0;
  std::string smi,sma;
  
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing a former crash in replaceCore" << std::endl;

  std::string pathName=getenv("RDBASE");
  pathName += "/Code/GraphMol/ChemTransforms/testData/oldcrash.mol";
  mol1 = MolFileToMol(pathName);
  TEST_ASSERT(mol1);

  sma = "N";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==6);
  smi = MolToSmiles(*mol2,true);
  // there's no way to guarantee the order here:
  TEST_ASSERT(smi=="[1*]CC.[2*]CC"||smi=="[2*]CC.[1*]CC");
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;

}

void testReplaceCorePositions() 
{
  ROMol *mol1=0,*mol2=0,*matcher1=0;
  std::string smi,sma;
  
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing that atom positions are correctly copied by replaceCore" << std::endl;

  std::string pathName=getenv("RDBASE");
  pathName += "/Code/GraphMol/ChemTransforms/testData/core_location.mol";
  mol1 = MolFileToMol(pathName);
  TEST_ASSERT(mol1);
  TEST_ASSERT(mol1->getNumAtoms()==8);

  sma = "c1ccccc1";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  mol2 = replaceCore(*mol1,*matcher1);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getNumAtoms()==4);

  RDGeom::Point3D op,np;
  op = mol1->getConformer().getAtomPos(6);
  np = mol2->getConformer().getAtomPos(0);  
  TEST_ASSERT(feq(op.x,np.x));
  TEST_ASSERT(feq(op.y,np.y));
  TEST_ASSERT(feq(op.z,np.z));
  op = mol1->getConformer().getAtomPos(7);
  np = mol2->getConformer().getAtomPos(1);  
  TEST_ASSERT(feq(op.x,np.x));
  TEST_ASSERT(feq(op.y,np.y));
  TEST_ASSERT(feq(op.z,np.z));
  op = mol1->getConformer().getAtomPos(0);
  np = mol2->getConformer().getAtomPos(2);  
  TEST_ASSERT(feq(op.x,np.x));
  TEST_ASSERT(feq(op.y,np.y));
  TEST_ASSERT(feq(op.z,np.z));
  op = mol1->getConformer().getAtomPos(4);
  np = mol2->getConformer().getAtomPos(3);  
  TEST_ASSERT(feq(op.x,np.x));
  TEST_ASSERT(feq(op.y,np.y));
  TEST_ASSERT(feq(op.z,np.z));


  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}


void testMurckoDecomp() 
{
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing murcko decomposition" << std::endl;


  const char *testMolecules[][2]=
    {
      {"C1CC1CCC1CC(C)C1", "C1CC1CCC1CCC1"},
      {"CC1CCC1", "C1CCC1"},
      {"NCNCC2CC2C1CC1O", "C1CC1C1CC1"},
      {"OC2C(C)C21C(N)C1C", "C2CC12CC1"}, // Spiro
      {"C1CC1C(=O)OC", "C1CC1"},          // Carbonyl outside scaffold
      {"C1CC1C=C", "C1CC1"},              // Double bond outside scaffold
      {"C1CC1C=CC1CC1C=CNNCO", "C1CC1C=CC1CC1"},              // Double bond in scaffold
      {"CC1CC1C(N)C1C(N)C1", "C1CC1CC1CC1"},
      {"C1CC1S(=O)C1CC1C=CNNCO", "C1CC1S(=O)C1CC1"}, // S=O group in scaffold
      {"O=SCNC1CC1S(=O)C1CC1C=CNNCO", "C1CC1S(=O)C1CC1"}, // S=O group outside scaffold
      {"C1CC1S(=O)(=O)C1CC1C=CNNCO", "C1CC1S(=O)(=O)C1CC1"}, // SO2 group in scaffold
      {"O=S(CNCNC)(=O)CNC1CC1S(=O)(=O)C1CC1C=CNNCO", "C1CC1S(=O)(=O)C1CC1"}, // SO2 group outside scaffold
      {"C1CC1C=NO","C1CC1"}, //Hydroxamide
      {"C1CC1C(C(C)C)=NC1CC1","C1CC1C=NC1CC1"}, //Hydroxamide
      {"C1CC1C#N","C1CC1"}, //Cyano group
      {"C1CC1C#CNC","C1CC1"}, //Acetylene group
      {"O=C1N(C)C(=O)N1C#CNC","O=C1NC(=O)N1"}, //Acetylene group
      {"[O-][N+](=O)c1cc(ccc1Cl)NS(=O)(=O)Cc2ccccc2","c1ccccc1NS(=O)(=O)Cc2ccccc2"},
      {"Cn1cccc1", "c1ccc[nH]1"},
      {"C1CC1[CH](C)C1CC1", "C1CC1CC1CC1"},
      {"CC(C)c1c(Cl)cc(F)c(Nc2ccc(C)cc2CC(=O)O)c1F","c1ccc(Nc2ccccc2)cc1"},
      {"C1CC1C[C@](Cl)(F)C1CCC1","C1CC1CCC1CCC1"},
      {"EOS","EOS"}
    };
  unsigned int i=0;
  while(1){
    std::string smi=testMolecules[i][0];
    std::string tgt=testMolecules[i][1];
    ++i;
    if(smi=="EOS") break;
    ROMol *mol=SmilesToMol(smi);
    ROMol *nMol=MurckoDecompose(*mol);
    MolOps::sanitizeMol(static_cast<RWMol &>(*nMol));
    TEST_ASSERT(nMol);
    TEST_ASSERT(nMol->getNumAtoms());
    delete mol;
    smi = MolToSmiles(*nMol,true);
    mol = SmilesToMol(tgt);
    TEST_ASSERT(mol);
    tgt=MolToSmiles(*mol,true);
    delete mol;
    std::cerr<<smi<<" "<<tgt<<std::endl;
    TEST_ASSERT(smi==tgt);
    delete nMol;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
  
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
  testReplaceCoreLabels();
  testReplaceCoreCrash();
  testReplaceCorePositions();

  testMurckoDecomp();

  BOOST_LOG(rdInfoLog) << "*******************************************************\n";
  return(0);
}

