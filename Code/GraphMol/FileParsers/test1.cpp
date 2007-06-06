// $Id$
//
//  Copyright (C) 2002-2005 Greg Landrum and Rational Discovery LLC
//       All Rights Reserved
//

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include "FileParsers.h"
#include "MolFileStereochem.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <string>

using namespace RDKit;

void test1(){
  BOOST_LOG(rdInfoLog) << "testing atom query parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/list-query.mol";
  RWMol *m = MolFileToMol(fName);
  //MolOps::sanitizeMol(*m);
  TEST_ASSERT(m);

  TEST_ASSERT(m->getNumAtoms()==6);
  std::string smi = MolToSmiles(*m);
  TEST_ASSERT(smi=="[Du]1=CC=CC=C1");
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
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test2(){
  BOOST_LOG(rdInfoLog) << "testing bond query parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/bond-query.mol";

  RWMol *m = MolFileToMol(fName);

  TEST_ASSERT(m->getNumAtoms()==5);
  std::string smi = MolToSmiles(*m);
  TEST_ASSERT(smi=="CC?CC=C");
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
  
  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


void test3(){
  // basic SD Parsing
  BOOST_LOG(rdInfoLog) << " ----------> Test3 "<< std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";
  std::string fName = rdbase + "test_data/esters.sdf";


  RWMOL_SPTR_VECT mols = SDFileToMols(fName);
  CHECK_INVARIANT(mols.size()==6,"");

  fName = rdbase+"test_data/earlyEOF.sdf";
  mols = SDFileToMols(fName);
  CHECK_INVARIANT(mols.size()==6,"");

  BOOST_LOG(rdInfoLog) << " Finished <---------- "<< std::endl;
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

  std::string molBlock = MolToMolBlock(*m);
  delete m;
  m = MolBlockToMol(molBlock);
  TEST_ASSERT(m);
  smi = MolToSmiles(*m);
  CHECK_INVARIANT(smi=="c1cc[cH-]c1",smi);

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
}  


void test6(){
  BOOST_LOG(rdInfoLog) << "testing chirality parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";
  std::string fName = rdbase + "test_data/chiral1.mol";

  RWMol *m;
  std::string smi;

#if 1
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
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@H](F)Cl");
  delete m;

  fName = rdbase+"test_data/chiral3a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@H](F)Cl");
  delete m;

#endif

  fName = rdbase+"test_data/chiral4.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@@H](F)Cl");
  delete m;

  fName = rdbase+"test_data/chiral4a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@@H](F)Cl");
  delete m;

  fName = rdbase+"test_data/chiral5.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==5);
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="CC(C)(Cl)Br");
  delete m;


  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test7(){
  BOOST_LOG(rdInfoLog) << "testing roundtrip chirality parsing" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";

  RWMol *m,*m2;
  std::string fName;
  std::string smi,molBlock,smi2;

#if 1
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
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@H](F)Cl");
  molBlock=MolToMolBlock(*m);
  //BOOST_LOG(rdInfoLog) << molBlock << std::endl;
  m2=MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  smi2 = MolToSmiles(*m2,true);
  TEST_ASSERT(smi==smi2);
  delete m;
  delete m2;

  fName = rdbase+"test_data/chiral4.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  smi = MolToSmiles(*m,true);
  TEST_ASSERT(smi=="C[C@@H](F)Cl");
  molBlock=MolToMolBlock(*m);
  //BOOST_LOG(rdInfoLog) << molBlock << std::endl;
  m2=MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  smi2 = MolToSmiles(*m2,true);
  TEST_ASSERT(smi==smi2);
  delete m;
  delete m2;
#endif
  fName = rdbase+"test_data/Issue142d.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
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
  BOOST_LOG(rdInfoLog) << molBlock << std::endl;
  m2=MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  smi2 = MolToSmiles(*m2,true);
  if(smi!=smi2){
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi==smi2);
  delete m;
  delete m2;
  fName = rdbase+"test_data/Issue142b.mol";
  m = MolFileToMol(fName);
  BOOST_LOG(rdInfoLog) << m->getNumAtoms() << "\n";
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==9);
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
  BOOST_LOG(rdInfoLog) << molBlock << std::endl;
  m2=MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  smi2 = MolToSmiles(*m2,true);
  if(smi!=smi2){
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi==smi2);
  delete m;
  delete m2;

  fName = rdbase+"test_data/issue142a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==28);
  smi = MolToSmiles(*m,true);
  m2=SmilesToMol(smi);
  smi2 = MolToSmiles(*m2,true);
  if(smi!=smi2){
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi==smi2);
  delete m2;

  molBlock=MolToMolBlock(*m);
  BOOST_LOG(rdInfoLog) << molBlock << std::endl;
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
  BOOST_LOG(rdInfoLog) << smi1 << std::endl;
  BOOST_LOG(rdInfoLog) << smi2 << std::endl;
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
    
  smi1 = MolToSmiles(*m1,true);
  TEST_ASSERT(smi1=="C[C@H]1CO1");
  
  MolOps::assignAtomChiralCodes(*m1);
  TEST_ASSERT(m1->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m1->getAtomWithIdx(1)->getProp("_CIPCode",smi2);
  TEST_ASSERT(smi2=="S");
  WedgeMolBonds(*m1);
  TEST_ASSERT(m1->getBondWithIdx(0)->getBondDir()==Bond::BEGINWEDGE);  
  TEST_ASSERT(m1->getBondWithIdx(1)->getBondDir()==Bond::NONE);  
  TEST_ASSERT(m1->getBondWithIdx(2)->getBondDir()==Bond::NONE);  
  TEST_ASSERT(m1->getBondWithIdx(3)->getBondDir()==Bond::NONE);  

  
  delete m1;


  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileChgLines(){
  BOOST_LOG(rdInfoLog) << "testing SF.Net Issue1603923: problems with multiple chg lines" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  RWMol *m1;
  std::string fName;

  fName = rdbase+"MolFileChgBug.mol";
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);
  TEST_ASSERT(m1->getAtomWithIdx(24)->getFormalCharge()==-1);
  TEST_ASSERT(m1->getAtomWithIdx(25)->getFormalCharge()==-1);
  

  delete m1;


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

  smi = MolToSmiles(*m1,true);
  TEST_ASSERT(smi=="C/C=C/Cl");

  fName = rdbase+"cistrans.2a.mol";
  delete m1;
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);

  smi = MolToSmiles(*m1,true);
  TEST_ASSERT(smi=="C/C=C\\Cl");

  fName = rdbase+"cistrans.1.mol";
  delete m1;
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);

  smi = MolToSmiles(*m1,true);
  TEST_ASSERT(smi=="C/C=C/C");

  fName = rdbase+"cistrans.2.mol";
  delete m1;
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);

  smi = MolToSmiles(*m1,true);
  TEST_ASSERT(smi=="C/C=C\\C");

  fName = rdbase+"cistrans.3.mol";
  delete m1;
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);

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

int main(int argc,char *argv[]){
  RDLog::InitLogs();
#if 0
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  testIssue145();
  testIssue148();
  test7();
  test8();
  testIssue180();
  testIssue264();
  testIssue399();
  testMolFileChgLines();
  testSymmetricDblBondStereochem();
#endif
  testRingDblBondStereochem();
  //testCrash();
  return 0;
}
