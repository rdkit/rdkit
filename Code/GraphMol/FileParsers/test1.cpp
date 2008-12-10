// $Id$
//
//  Copyright (C) 2002-2008 Greg Landrum and Rational Discovery LLC
//       All Rights Reserved
//

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h> 
#include "FileParsers.h"
#include "MolFileStereochem.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>

#include <string>

using namespace RDKit;

void test1(){
  BOOST_LOG(rdInfoLog) << "testing atom query parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/list-query.mol";
  RWMol *m = MolFileToMol(fName,false);
  //MolOps::sanitizeMol(*m);
  TEST_ASSERT(m);

  TEST_ASSERT(m->getNumAtoms()==6);
  std::string smi = MolToSmiles(*m);
  TEST_ASSERT(smi=="C1=CC=CC=C1");  

  smi = MolToSmarts(*m);
  TEST_ASSERT(smi=="[#6]1=[#6]-[#6]=[#6]-[#6]=[#6,#7,#15]-1");


  smi = "C1=CC=CC=C1";
  RWMol *m2 = SmilesToMol(smi,false,false);
  MatchVectType mv;
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==6);
  // sanitize it, which will aromatize the bonds... we will not match:
  MolOps::sanitizeMol(*m2);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);

  delete m2;
  smi = "N1=CC=CC=C1";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==6);
  delete m2;
  smi = "S1=CC=CC=C1";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;
  smi = "P1=CC=CC=C1";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==6);
  delete m2;
  
  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/not-list-query.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC(=C)C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==4);
  delete m2;
  
  smi = "CC(=O)C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  delete m2;
  
  smi = "CC(=N)C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  delete m2;
  
  smi = "CC(=O)C(=C)C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==4);
  delete m2;
  
  smi = "C(=C)C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  delete m2;

  // make sure new-style atom lists override old-style atom lists:
  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/conflicting-list-query.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC(=C)C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==4);
  delete m2;
  
  smi = "CC(=O)C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  delete m2;
  
  // longer list queries, this was issue 2413431:
  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/list-query-long.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(14)->hasQuery());

  smi = "C1COC2=CC3=CC4=C(C=CC=C4)C=C3C=C2C1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  delete m2;
  smi = "C1C[Se]C2=CC3=CC4=C(C=CC=C4)C=C3C=C2C1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  delete m2;
  smi = "C1C[Te]C2=CC3=CC4=C(C=CC=C4)C=C3C=C2C1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  delete m2;
  smi = "C1C[As]C2=CC3=CC4=C(C=CC=C4)C=C3C=C2C1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  delete m2;

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test2(){
  BOOST_LOG(rdInfoLog) << "testing bond query parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/bond-query.mol";

  RWMol *m = MolFileToMol(fName);

  TEST_ASSERT(m->getNumAtoms()==5);
  std::string smi = MolToSmiles(*m);
  TEST_ASSERT(smi=="C=CC?CC");
  smi = "C1=CC=CC=C1";
  RWMol *m2 = SmilesToMol(smi,false,false);
  MatchVectType mv;
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==5);
  // sanitize it (making bonds aromatic) ... we will not match:
  MolOps::sanitizeMol(*m2);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);

  delete m2;
  smi = "C=CC=CC=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==5);
  
  delete m2;
  smi = "C=CCCC=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==5);
  
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/RingBondQuery.mol";
  delete m;
  m = MolFileToMol(fName);
  TEST_ASSERT(m->getNumAtoms()==5);
  delete m2;
  smi = "C1CCC1C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  delete m2;
  smi = "C1CC2C1C2";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==5);

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ChainBondQuery.mol";
  delete m;
  m = MolFileToMol(fName);
  TEST_ASSERT(m->getNumAtoms()==5);
  delete m2;
  smi = "C1CCC1C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==5);
  delete m2;
  smi = "C1CC2C1C2";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));

  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


void test4(){
  // basic writing test
  BOOST_LOG(rdInfoLog) << " ----------> Test4 "<< std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";
  std::string fName = rdbase + "test_data/mol1.mol";

  RWMol *m = MolFileToMol(fName);
  TEST_ASSERT(m);
  std::string smi=MolToSmiles(*m);
  CHECK_INVARIANT(smi=="c1cc[cH-]c1",smi);

  TEST_ASSERT(m->getConformer().is3D()==false);

  std::string molBlock = MolToMolBlock(*m);
  delete m;
  m = MolBlockToMol(molBlock);
  TEST_ASSERT(m);
  smi = MolToSmiles(*m);
  CHECK_INVARIANT(smi=="c1cc[cH-]c1",smi);
  TEST_ASSERT(m->getConformer().is3D()==false);

  BOOST_LOG(rdInfoLog) << " Finished <---------- "<< std::endl;
}

void test5(){
  // formerly problematic molecules
  BOOST_LOG(rdInfoLog) << " ----------> Test5 "<< std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";
  std::string fName = rdbase + "test_data/issue123.mol";

  RWMol *m = MolFileToMol(fName);
  CHECK_INVARIANT(m,"");
  TEST_ASSERT(m->getNumAtoms()==23);
  TEST_ASSERT(m->getConformer().is3D()==true);
  std::string molBlock = MolToMolBlock(*m);
  delete m;
  m = MolBlockToMol(molBlock);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==23);
  TEST_ASSERT(m->getConformer().is3D()==true);

  delete m;
  // now try without removing the Hs:
  m = MolFileToMol(fName,true,false);
  CHECK_INVARIANT(m,"");
  TEST_ASSERT(m->getNumAtoms()==39);

}  

void test6(){
  BOOST_LOG(rdInfoLog) << "testing chirality parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";
  std::string fName = rdbase + "test_data/chiral1.mol";

  RWMol *m;
  std::string smi,cip;

  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==5);
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@](F)(Cl)Br");
  delete m;

  fName = rdbase+"test_data/chiral1a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==5);
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@](F)(Cl)Br");
  delete m;

  fName = rdbase+"test_data/chiral2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==5);
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@@](F)(Cl)Br");
  delete m;

  fName = rdbase+"test_data/chiral2a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==5);
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@@](F)(Cl)Br");
  delete m;

  fName = rdbase+"test_data/chiral3.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");  
#if 1
  smi = MolToSmiles(*m,true);
  //BOOST_LOG(rdInfoLog) << " smi: " << smi << std::endl;
  TEST_ASSERT(smi=="C[C@H](F)Cl");
#endif
  delete m;

  fName = rdbase+"test_data/chiral3a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");  
#if 1
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@H](F)Cl");
#endif
  delete m;

  fName = rdbase+"test_data/chiral4.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");  
#if 1
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@@H](F)Cl");
#endif
  delete m;

  fName = rdbase+"test_data/chiral4a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");  
#if 1
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@@H](F)Cl");
#endif
  delete m;

  fName = rdbase+"test_data/chiral5.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==5);
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp("_CIPCode"));
#if 1
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="CC(C)(Cl)Br");
#endif
  delete m;


  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test7(){
  BOOST_LOG(rdInfoLog) << "testing roundtrip chirality parsing" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";

  RWMol *m,*m2;
  std::string fName;
  std::string smi,molBlock,smi2,cip;

  fName = rdbase+"test_data/chiral1.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==5);
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@](F)(Cl)Br");
  molBlock=MolToMolBlock(*m);
  m2=MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  smi2 = MolToSmiles(*m2,true);
  TEST_ASSERT(smi==smi2);
  delete m;
  delete m2;

  fName = rdbase+"test_data/chiral2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==5);
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@@](F)(Cl)Br");
  molBlock=MolToMolBlock(*m);
  m2=MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  smi2 = MolToSmiles(*m2,true);
  TEST_ASSERT(smi==smi2);
  delete m;
  delete m2;
  fName = rdbase+"test_data/chiral3.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");  
#if 1
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@H](F)Cl");
  molBlock=MolToMolBlock(*m);
  //BOOST_LOG(rdInfoLog) << molBlock << std::endl;
  m2=MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  MolOps::assignAtomChiralCodes(*m2);
  TEST_ASSERT(m2->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m2->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");  
  smi2 = MolToSmiles(*m2,true);
  TEST_ASSERT(smi==smi2);
  delete m2;
#endif
  delete m;

  fName = rdbase+"test_data/chiral4.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");  
#if 1
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@@H](F)Cl");
  molBlock=MolToMolBlock(*m);
  //BOOST_LOG(rdInfoLog) << molBlock << std::endl;
  m2=MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  MolOps::assignAtomChiralCodes(*m2);
  TEST_ASSERT(m2->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m2->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");  
  //smi2 = MolToSmiles(*m2,true);
  //TEST_ASSERT(smi==smi2);
  delete m2;
#endif
  delete m;

  fName = rdbase+"test_data/Issue142d.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  MolOps::assignAtomChiralCodes(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");  
#if 1
  smi = MolToSmiles(*m,true);
  m2=SmilesToMol(smi);
  smi2 = MolToSmiles(*m2,true);
  if(smi!=smi2){
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi==smi2);
  delete m2;
  //BOOST_LOG(rdInfoLog) << "SMI: "<< smi << std::endl;
  molBlock=MolToMolBlock(*m);
  //BOOST_LOG(rdInfoLog) << molBlock << std::endl;
  m2=MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  smi2 = MolToSmiles(*m2,true);
  if(smi!=smi2){
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi==smi2);
  delete m2;
#endif

  delete m;
  fName = rdbase+"test_data/Issue142b.mol";
  m = MolFileToMol(fName);
  //BOOST_LOG(rdInfoLog) << m->getNumAtoms() << "\n";
  //BOOST_LOG(rdInfoLog) << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << std::endl;
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
#if 1
  smi = MolToSmiles(*m,true);
  m2=SmilesToMol(smi);
  smi2 = MolToSmiles(*m2,true);
  if(smi!=smi2){
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi==smi2);
  delete m2;
  //BOOST_LOG(rdInfoLog) << "SMI: "<< smi << std::endl;
  BOOST_LOG(rdInfoLog) << m->getNumAtoms() << " " << m->getConformer().getNumAtoms() << "\n";
  molBlock=MolToMolBlock(*m);
  //BOOST_LOG(rdInfoLog) << molBlock << std::endl;
  m2=MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  smi2 = MolToSmiles(*m2,true);
  if(smi!=smi2){
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi==smi2);
  delete m2;
#endif
  delete m;

  fName = rdbase+"test_data/issue142a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==28);
#if 1
  smi = MolToSmiles(*m,true);
  m2=SmilesToMol(smi);
  smi2 = MolToSmiles(*m2,true);
  if(smi!=smi2){
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi==smi2);
  delete m2;

  molBlock=MolToMolBlock(*m);
  //BOOST_LOG(rdInfoLog) << molBlock << std::endl;
  m2=MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  smi2 = MolToSmiles(*m2,true);
  if(smi!=smi2){
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi==smi2);
  delete m2;
#endif
  delete m;


  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test8(){
  BOOST_LOG(rdInfoLog) << "testing reading without sanitization" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";

  RWMol *m;
  std::string fName;
  std::string smi,molBlock,smi2;

  // in this case the test means to not remove Hs:
  fName = rdbase+"test_data/unsanitary.mol";
  m = MolFileToMol(fName,false);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==6);

  delete m;
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==3);

  delete m;
  fName = rdbase+"test_data/unsanitary2.mol";
  m = MolFileToMol(fName,false);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==9);
  
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


void testIssue145(){
  BOOST_LOG(rdInfoLog) << "testing Issue145:\n Mol parsing: molecule yields non-canonical smiles from mol block" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";

  RWMol *m,*m2;
  std::string fName;
  std::string smi,molBlock,smi2;

  fName = rdbase+"test_data/issue145.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==19);
  smi = MolToSmiles(*m,true);
  m2=SmilesToMol(smi);
  smi2 = MolToSmiles(*m2,true);
  if(smi!=smi2){
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi==smi2);
  delete m2;

  molBlock=MolToMolBlock(*m);
  m2=MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  smi2 = MolToSmiles(*m2,true);
  if(smi!=smi2){
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi==smi2);
  delete m;
  delete m2;


  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue148(){
  BOOST_LOG(rdInfoLog) << "testing Issue148:\n Mol files containing mis-drawn nitro groups not properly parsed" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";

  RWMol *m;
  std::string fName;
  std::string smi,molBlock,smi2;

  fName = rdbase+"test_data/issue148.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==9);
  TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge()==1);
  

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue180(){
  BOOST_LOG(rdInfoLog) << "testing Issue180: bad Z/E assignments" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";

  RWMol *m;
  std::string fName;
  std::string code;
  Bond *bond;

  fName = rdbase+"test_data/Issue180.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  bond = m->getBondWithIdx(2);
  TEST_ASSERT(bond->getBondType()==Bond::DOUBLE);
  TEST_ASSERT(bond->getStereo()==Bond::STEREOZ);

  bond = m->getBondWithIdx(5);
  TEST_ASSERT(bond->getBondType()==Bond::DOUBLE);
  TEST_ASSERT(bond->getStereo()==Bond::STEREOE);
  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue264(){
  BOOST_LOG(rdInfoLog) << "testing Issue264: bad stereochemistry from mol files" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";

  RWMol *m1,*m2;
  std::string smi1,smi2;
  std::string fName;


  fName = rdbase+"test_data/Issue264-1.mol";
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);

  fName = rdbase+"test_data/Issue264-2.mol";
  m2 = MolFileToMol(fName);
  TEST_ASSERT(m2);

  smi1 = MolToSmiles(*m1,false);
  smi2 = MolToSmiles(*m2,false);
  TEST_ASSERT(smi1==smi2);
  
  smi1 = MolToSmiles(*m1,true);
  smi2 = MolToSmiles(*m2,true);
  //BOOST_LOG(rdInfoLog) << smi1 << std::endl;
  //BOOST_LOG(rdInfoLog) << smi2 << std::endl;
  TEST_ASSERT(smi1!=smi2);
  
  delete m1;
  delete m2;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue399(){
  BOOST_LOG(rdInfoLog) << "testing Issue399: bond wedging cleanup" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  RWMol *m1;
  std::string smi1,smi2;
  std::string fName;

  fName = rdbase+"Issue399a.mol";
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);
#if 1    
  smi1 = MolToSmiles(*m1,true);
  TEST_ASSERT(smi1=="C[C@H]1CO1");
#endif  
  MolOps::assignAtomChiralCodes(*m1);
  TEST_ASSERT(m1->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m1->getAtomWithIdx(1)->getProp("_CIPCode",smi2);
  TEST_ASSERT(smi2=="S");
#if 1
  WedgeMolBonds(*m1,&m1->getConformer());
  TEST_ASSERT(m1->getBondWithIdx(0)->getBondDir()==Bond::BEGINWEDGE);  
  TEST_ASSERT(m1->getBondWithIdx(1)->getBondDir()==Bond::NONE);  
  TEST_ASSERT(m1->getBondWithIdx(2)->getBondDir()==Bond::NONE);  
  TEST_ASSERT(m1->getBondWithIdx(3)->getBondDir()==Bond::NONE);  
#endif
  
  delete m1;


  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileChgLines(){
  BOOST_LOG(rdInfoLog) << "testing handling of charge lines" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";


  // SF.Net Issue1603923: problems with multiple chg lines
  {
    RWMol *m1;
    std::string fName;
    fName = rdbase+"MolFileChgBug.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getAtomWithIdx(24)->getFormalCharge()==-1);
    TEST_ASSERT(m1->getAtomWithIdx(25)->getFormalCharge()==-1);
    delete m1;
  }

  // many charges in one molecule:
  {
    RWMol *m1;
    std::string fName;
    fName = rdbase+"manycharges.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getFormalCharge()==-1);
    TEST_ASSERT(m1->getAtomWithIdx(13)->getFormalCharge()==-1);

    std::string molBlock = MolToMolBlock(*m1);
    //std::cerr<<molBlock<<std::endl;
    delete m1;
    m1 = MolBlockToMol(molBlock);
    //m1->debugMol(std::cerr);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getFormalCharge()==-1);
    TEST_ASSERT(m1->getAtomWithIdx(13)->getFormalCharge()==-1);
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


void testDblBondStereochem(){
  BOOST_LOG(rdInfoLog) << "testing basic double bond stereochemistry" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    RWMol *m1;
    std::string fName=rdbase+"simple_z.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getBondWithIdx(0)->getStereo()==Bond::STEREOZ);
    delete m1;
  }

  {
    RWMol *m1;
    std::string fName=rdbase+"simple_e.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getBondWithIdx(0)->getStereo()==Bond::STEREOE);
    delete m1;
  }

  {
    RWMol *m1;
    std::string fName=rdbase+"simple_either.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getBondWithIdx(0)->getStereo()==Bond::STEREONONE);
    TEST_ASSERT(m1->getBondWithIdx(0)->getBondDir()==Bond::EITHERDOUBLE);
    delete m1;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}



void testSymmetricDblBondStereochem(){
  // this was sf.net issue 1718794:
  // http://sourceforge.net/tracker/index.php?func=detail&aid=1718794&group_id=160139&atid=814650)
  BOOST_LOG(rdInfoLog) << "testing double bonds with symmetric substituents" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  RWMol *m1;
  std::string fName,smi;

  fName = rdbase+"cistrans.1a.mol";
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);
  TEST_ASSERT(m1->getBondWithIdx(0)->getStereo()==Bond::STEREOE);
  smi = MolToSmiles(*m1,true);
  TEST_ASSERT(smi=="C/C=C/Cl");

  fName = rdbase+"cistrans.2a.mol";
  delete m1;
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);
  TEST_ASSERT(m1->getBondWithIdx(0)->getStereo()==Bond::STEREOZ);

  smi = MolToSmiles(*m1,true);
  TEST_ASSERT(smi=="C/C=C\\Cl");

  fName = rdbase+"cistrans.1.mol";
  delete m1;
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);
  TEST_ASSERT(m1->getBondWithIdx(0)->getStereo()==Bond::STEREOE);

  smi = MolToSmiles(*m1,true);
  TEST_ASSERT(smi=="C/C=C/C");

  fName = rdbase+"cistrans.2.mol";
  delete m1;
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);
  TEST_ASSERT(m1->getBondWithIdx(0)->getStereo()==Bond::STEREOZ);

  smi = MolToSmiles(*m1,true);
  TEST_ASSERT(smi=="C/C=C\\C");

  fName = rdbase+"cistrans.3.mol";
  delete m1;
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);
  TEST_ASSERT(m1->getBondWithIdx(0)->getStereo()==Bond::STEREONONE);

  smi = MolToSmiles(*m1,true);
  TEST_ASSERT(smi=="CC=CC");


  delete m1;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testRingDblBondStereochem(){
  // this was sf.net issue 1725068:
  // http://sourceforge.net/tracker/index.php?func=detail&aid=1725068&group_id=160139&atid=814650
  BOOST_LOG(rdInfoLog) << "testing double bonds in rings with stereochem specifications" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  RWMol *m1;
  std::string fName,smi;

  fName = rdbase+"badringstereochem3.mol";
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);

  smi = MolToSmiles(*m1,true);
  TEST_ASSERT(smi.find("/",0)==std::string::npos);
  TEST_ASSERT(smi.find("\\",0)==std::string::npos);
  delete m1;

  fName = rdbase+"badringstereochem2.mol";
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);

  smi = MolToSmiles(*m1,true);
  TEST_ASSERT(smi.find("/",0)==std::string::npos);
  TEST_ASSERT(smi.find("\\",0)==std::string::npos);
  delete m1;

  fName = rdbase+"badringstereochem.mol";
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);

  smi = MolToSmiles(*m1,true);
  TEST_ASSERT(smi.find("/",0)==std::string::npos);
  TEST_ASSERT(smi.find("\\",0)==std::string::npos);
  delete m1;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}



void testMolFileRGroups(){
  BOOST_LOG(rdInfoLog) << "testing mol file R-group parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/rgroups1.mol";
  RWMol *m = MolFileToMol(fName);
  TEST_ASSERT(m);
  unsigned int idx;
  
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_MolFileRLabel"));
  m->getAtomWithIdx(3)->getProp("_MolFileRLabel",idx);
  TEST_ASSERT(idx==2);
  
  TEST_ASSERT(m->getAtomWithIdx(4)->hasProp("_MolFileRLabel"));
  m->getAtomWithIdx(4)->getProp("_MolFileRLabel",idx);
  TEST_ASSERT(idx==1);
  
  RWMol *m2;
  MatchVectType mv;

  std::string smi;
  smi = "C1C(O)C1C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==5);
  delete m2;
  
  smi = "C1CC(O)C1C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  delete m2;
  
  smi = "C1C(CO)C1CC";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==5);
  delete m2;
  
  smi = "CC(=O)CC";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  delete m2;

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/rgroups2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_MolFileRLabel"));
  m->getAtomWithIdx(3)->getProp("_MolFileRLabel",idx);
  TEST_ASSERT(idx==1);
  
  TEST_ASSERT(m->getAtomWithIdx(4)->hasProp("_MolFileRLabel"));
  m->getAtomWithIdx(4)->getProp("_MolFileRLabel",idx);
  TEST_ASSERT(idx==1);

  smi = "C1C(O)C1C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==5);
  delete m2;
  
  smi = "C1CC(O)C1C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  delete m2;
  
  smi = "C1C(CO)C1CC";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==5);
  delete m2;
  
  smi = "CC(=O)CC";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  delete m2;

  
  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileDegreeQueries(){
  BOOST_LOG(rdInfoLog) << "testing mol file degree queries" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/subst1.mol";
  RWMol *m = MolFileToMol(fName);
  TEST_ASSERT(m);
    
  RWMol *m2;
  MatchVectType mv;
  std::string smi;

  smi = "CC(=O)O";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==4);
  delete m2;
  
  smi = "CC(=O)C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==4);
  delete m2;
  
  smi = "CC(=O)";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  delete m2;
  
  
  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/subst2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  
  smi = "CC(=O)O";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==4);
  delete m2;
  
  smi = "CC(=O)C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==4);
  delete m2;
  
  smi = "CC(O)C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==4);
  delete m2;
  
  smi = "CC(O)";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  delete m2;

  smi = "CC(O)(C)C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  delete m2;

  
  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/subst3.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  
  smi = "CC(=O)O";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==3);
  delete m2;
  
  smi = "CC(=O)C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==3);
  delete m2;
  
  smi = "CC(O)C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  
  smi = "CC(O)";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  delete m2;

  smi = "CC(O)(C)C";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  delete m2;

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/subst4.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  
  smi = "CC(=O)O";
  m2 = SmilesToMol(smi,false,false);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  delete m2;

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileRBCQueries(){
  BOOST_LOG(rdInfoLog) << "testing mol file ring-bond count queries" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;
  RWMol *m;
  RWMol *m2;
  MatchVectType mv;
  std::string smi;

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  
  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;
  smi = "C1CCC1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==2);
  delete m2;
  smi = "C12C3C4C1C5C2C3C45";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_3.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;
  smi = "C1CCC1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;
  smi = "C12C3C4C1C5C2C3C45";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==2);
  delete m2;

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_0.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==2);
  delete m2;
  smi = "C1CCC1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;
  

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_4.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;
  smi = "C1CCC1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;
  smi = "C12C3C4C1C5C2C3C45";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;

  smi = "C1CS234C5CC2CC13CC4C5";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==2);
  delete m2;
  smi = "C1C2CC3CC4CC1S234";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==2);


  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_star.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==2);
  delete m2;
  smi = "C1CCC1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;
  
  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_star2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;
  smi = "C1CCC1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==4);
  delete m2;
  smi = "C1CC2C1CC2";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==4);
  delete m2;
  smi = "C12C3C4C1C5C2C3C45";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_star3.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;
  smi = "C1CCC1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==4);
  delete m2;
  smi = "C1CC2C1CC2";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;
  smi = "C12C3C4C1C5C2C3C45";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


void testMolFileUnsaturationQueries(){
  BOOST_LOG(rdInfoLog) << "testing mol file unsaturation queries" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;
  RWMol *m;
  RWMol *m2;
  MatchVectType mv;
  std::string smi;

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/unsaturation.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CO";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;

  smi = "C=O";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==2);
  delete m2;

  smi = "CCO";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==0);
  delete m2;

  smi = "C=CO";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==2);
  delete m2;
  
  smi = "C#CO";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2,*m,mv));
  TEST_ASSERT(mv.size()==2);
  delete m2;
  
  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileQueryToSmarts(){
  BOOST_LOG(rdInfoLog) << "testing mol file queries -> SMARTS " << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;
  RWMol *m;
  std::string sma;

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  sma = MolToSmarts(*m,true);
  TEST_ASSERT(sma=="[#6&$(*(@*)@*)&!$(*(@*)(@*)@*)]-[#6]")
  
  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_3.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  sma = MolToSmarts(*m,true);
  TEST_ASSERT(sma=="[#6&$(*(@*)(@*)@*)&!$(*(@*)(@*)(@*)@*)]-[#6]")

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_0.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  sma = MolToSmarts(*m,true);
  TEST_ASSERT(sma=="[#6&!$(*@*)]-[#6]")

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_4.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  sma = MolToSmarts(*m,true);
  TEST_ASSERT(sma=="[#16&$(*(@*)(@*)(@*)@*)&!$(*(@*)(@*)(@*)(@*)@*)]-[#6]")

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_star.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  sma = MolToSmarts(*m,true);
  TEST_ASSERT(sma=="[#6&!$(*@*)]-[#6]")
  
  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_star2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  sma = MolToSmarts(*m,true);
  TEST_ASSERT(sma.find("[#6&$(*(@*)@*)&!$(*(@*)(@*)@*)]")!=std::string::npos);

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/unsaturation.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  sma = MolToSmarts(*m,true);
  TEST_ASSERT(sma=="[#6;$(*=,:,#*)]~[#8]")

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMissingFiles(){
  BOOST_LOG(rdInfoLog) << "testing handling of missing files" << std::endl;

  std::string fName;
  bool ok;
  RWMol *m;

  fName = "bogus_file.mol";
  ok=false;
  try{
    m = MolFileToMol(fName);
  } catch (BadFileException &e){
    ok=true;
  }
  TEST_ASSERT(ok);

  ok=false;
  try{
    m = TPLFileToMol(fName);
  } catch (BadFileException &e){
    ok=true;
  }
  TEST_ASSERT(ok);

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}



void testIssue1965035(){
  BOOST_LOG(rdInfoLog) << "testing issue Issue1965035: problems with WedgeMolBonds " << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;
  RWMol *m;
  std::string sma;

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/Issue1965035.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  // the mol file parser removes bond wedging info:
  TEST_ASSERT(m->getBondWithIdx(4)->getBondDir()==Bond::NONE);
  // but a chiral tag is assigned:
  TEST_ASSERT(m->getAtomWithIdx(2)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);

  WedgeMolBonds(*m,&m->getConformer());
  TEST_ASSERT(m->getBondWithIdx(4)->getBondDir()==Bond::BEGINDASH);
  
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testRadicals(){
  BOOST_LOG(rdInfoLog) << "testing handling of radicals " << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;
  RWMol *m;
  std::string smiles;

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/radical.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==0);
  TEST_ASSERT(m->getAtomWithIdx(1)->getNumRadicalElectrons()==1);
  
  std::string molBlock = MolToMolBlock(*m);
  delete m;
  m = MolBlockToMol(molBlock);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==0);
  TEST_ASSERT(m->getAtomWithIdx(1)->getNumRadicalElectrons()==1);
  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;

}

void testBadBondOrders(){
  BOOST_LOG(rdInfoLog) << "testing handling of bogus bond orders (issue 2337369)" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;
  RWMol *m;

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/bondorder0.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getBondType()==Bond::UNSPECIFIED);
  TEST_ASSERT(!m->getBondBetweenAtoms(0,1)->hasQuery());
  delete m;

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/bondorder9.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getBondBetweenAtoms(0,1)->hasQuery());
  TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getQuery()->getDescription()=="BondNull");
  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;

}

int main(int argc,char *argv[]){
  RDLog::InitLogs();
#if 1
  test1();
  test2();
  test4();
  test5();
  test6();
  testIssue145();
  testIssue148();
  test7();
  test8();
#endif
  testIssue180();
  testIssue264();
  testIssue399();
  testMolFileChgLines();
  testDblBondStereochem();
  testSymmetricDblBondStereochem();
  testRingDblBondStereochem();
  testMolFileRGroups();
  testMolFileDegreeQueries();
  testMolFileRBCQueries();
  testMolFileUnsaturationQueries();
  testMolFileQueryToSmarts();
  testMissingFiles();
  testIssue1965035();
  testRadicals();
  testBadBondOrders();
  //testCrash();
  return 0;
}
