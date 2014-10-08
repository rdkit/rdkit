//  $Id$
// 
//   Copyright (C) 2002-2013 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <iostream>
#include <map>
#include <algorithm>
#include <boost/foreach.hpp>


#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

using namespace RDKit;
using namespace std;
RWMol _t;
typedef class ROMol Mol;

void test1(){
  string smi;
  Mol *m;
  INT_VECT iv;
  unsigned int count;
  std::vector<ROMOL_SPTR> frags;

  smi = "CCCC(=O)O";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m,"smiles parse failed");
  count = MolOps::getMolFrags(*m,iv);
  CHECK_INVARIANT(count==1,"bad frag count");
  frags = MolOps::getMolFrags(*m);
  CHECK_INVARIANT(frags.size()==1,"bad frag count");
  TEST_ASSERT(frags[0]->getNumAtoms()==6);
  delete m;

  smi = "CCCC(=O)[O-].[Na+]";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m,"smiles parse failed");
  count = MolOps::getMolFrags(*m,iv);
  CHECK_INVARIANT(count==2,"bad frag count");
  frags = MolOps::getMolFrags(*m);
  CHECK_INVARIANT(frags.size()==2,"bad frag count");
  TEST_ASSERT(frags[0]->getNumAtoms()==6);
  TEST_ASSERT(frags[1]->getNumAtoms()==1);
  frags = MolOps::getMolFrags(*m,true,&iv);
  CHECK_INVARIANT(frags.size()==2,"bad frag count");
  TEST_ASSERT(frags[0]->getNumAtoms()==6);
  TEST_ASSERT(frags[1]->getNumAtoms()==1);
  TEST_ASSERT(iv.size()==7);
  TEST_ASSERT(iv[0]==0)
    TEST_ASSERT(iv[6]==1)
    delete m;

  smi = "CCCC(=O)[O-].[Na+].[NH4+].[Cl-]";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m,"smiles parse failed");
  count = MolOps::getMolFrags(*m,iv);
  CHECK_INVARIANT(count==4,"bad frag count");
  frags = MolOps::getMolFrags(*m);
  CHECK_INVARIANT(frags.size()==4,"bad frag count");
  TEST_ASSERT(frags[0]->getNumAtoms()==6);
  TEST_ASSERT(frags[1]->getNumAtoms()==1);
  TEST_ASSERT(frags[2]->getNumAtoms()==1);
  TEST_ASSERT(frags[3]->getNumAtoms()==1);
  delete m;

};

void test2(){
  string smi;
  Mol *m;
  INT_VECT iv;
  int count;
  smi = "CCCC(=O)O";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m,"smiles parse failed");
  count = MolOps::getMolFrags(*m,iv);
  CHECK_INVARIANT(count==1,"bad frag count");
  CHECK_INVARIANT(iv[0]==0,"bad frag membership");
  CHECK_INVARIANT(iv[1]==0,"bad frag membership");
  CHECK_INVARIANT(iv[4]==0,"bad frag membership");
  CHECK_INVARIANT(iv[5]==0,"bad frag membership");
  delete m;

  smi = "CCCC(=O)[O-].[Na+]";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m,"smiles parse failed");
  count = MolOps::getMolFrags(*m,iv);
  CHECK_INVARIANT(count==2,"bad frag count");
  CHECK_INVARIANT(iv[0]==0,"bad frag membership");
  CHECK_INVARIANT(iv[1]==0,"bad frag membership");
  CHECK_INVARIANT(iv[4]==0,"bad frag membership");
  CHECK_INVARIANT(iv[5]==0,"bad frag membership");
  CHECK_INVARIANT(iv[6]==1,"bad frag membership");
  delete m;

};



void test3(){
  string smi;
  Mol *m;
  VECT_INT_VECT sssr;
  INT_VECT rings;
  int count;

  smi = "C1CC1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getRingInfo()->numRings()==1);
  count = MolOps::findSSSR(*m,sssr);
  TEST_ASSERT(count==1);
  TEST_ASSERT(sssr[0].size()==3);
  TEST_ASSERT(m->getRingInfo()->numRings()==1);
  
  for(unsigned int i=0;i<m->getNumAtoms();i++){
    TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(i,3));
    TEST_ASSERT(!m->getRingInfo()->isAtomInRingOfSize(i,4));
    TEST_ASSERT(m->getRingInfo()->numAtomRings(i)==1);
  }
  for(unsigned int i=0;i<m->getNumBonds();i++){
    TEST_ASSERT(m->getRingInfo()->isBondInRingOfSize(i,3));
    TEST_ASSERT(!m->getRingInfo()->isBondInRingOfSize(i,4));
    TEST_ASSERT(m->getRingInfo()->numBondRings(i)==1);
  }
  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;
  
  smi = "C1CCC1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getRingInfo()->numRings()==1);
  count = MolOps::findSSSR(*m,sssr);
  TEST_ASSERT(count==1);
  TEST_ASSERT(sssr[0].size()==4);
  TEST_ASSERT(m->getRingInfo()->numRings()==1);
  for(unsigned int i=0;i<m->getNumAtoms();i++){
    TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(i,4));
    TEST_ASSERT(!m->getRingInfo()->isAtomInRingOfSize(i,3));
    TEST_ASSERT(m->getRingInfo()->numAtomRings(i)==1);
  }
  TEST_ASSERT(m->getRingInfo()->isBondInRingOfSize(0,4));
  TEST_ASSERT(m->getRingInfo()->numBondRings(0)==1);

  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;


  smi = "C1CCCCCC1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getRingInfo()->numRings()==1);
  count = MolOps::findSSSR(*m,sssr);
  TEST_ASSERT(count==1);
  TEST_ASSERT(sssr[0].size()==7);
  TEST_ASSERT(m->getRingInfo()->numRings()==1);
  for(unsigned int i=0;i<m->getNumAtoms();i++){
    TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(i,7));
    TEST_ASSERT(m->getRingInfo()->numAtomRings(i)==1);
  }
  TEST_ASSERT(m->getRingInfo()->isBondInRingOfSize(0,7));
  TEST_ASSERT(m->getRingInfo()->numBondRings(0)==1);

  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;



  smi = "C1C(CCC)CC(C(C)CCC(CC))CCC1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  count = MolOps::findSSSR(*m,sssr);
  TEST_ASSERT(count==1);
  TEST_ASSERT(sssr[0].size()==7);
  TEST_ASSERT(m->getRingInfo()->numAtomRings(0)==1);
  TEST_ASSERT(m->getRingInfo()->numBondRings(m->getBondBetweenAtoms(0,1)->getIdx()));
  TEST_ASSERT(!m->getRingInfo()->numBondRings(m->getBondBetweenAtoms(1,2)->getIdx()));
  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;

  smi = "CC1CCC1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  count = MolOps::findSSSR(*m,sssr);
  TEST_ASSERT(count==1);
  TEST_ASSERT(sssr[0].size()==4);
  TEST_ASSERT(!m->getBondBetweenAtoms(0,1)->hasProp("ringMembership"));
  TEST_ASSERT(!m->getRingInfo()->numBondRings(m->getBondBetweenAtoms(0,1)->getIdx()));
  TEST_ASSERT(m->getRingInfo()->numBondRings(m->getBondBetweenAtoms(1,2)->getIdx()));
  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;
  
  smi = "CC1C(C2)CCC2C1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  count = MolOps::findSSSR(*m, sssr);
  TEST_ASSERT(count==2)
    TEST_ASSERT(sssr[0].size()==5);
  TEST_ASSERT(sssr[1].size()==5);
  TEST_ASSERT(!m->getRingInfo()->isAtomInRingOfSize(0,5));
  TEST_ASSERT(m->getRingInfo()->numAtomRings(0)==0);
  TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(1,5));
  TEST_ASSERT(m->getRingInfo()->numAtomRings(1)==1);
  TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(2,5));
  TEST_ASSERT(m->getRingInfo()->numAtomRings(2)==2);
  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;

  smi = "C(C1C2C3C41)(C2C35)C45"; // cubane
  //smi = "C1(C2C3C4C5C6C72)C3C4C5C6C71"; // from Figureras paper
  //smi = "C17C5C4C3C2C1C6C2C3C4C5C67"; 
  // we cannot use the sanitzation code, because that finds *symmetric*
  // rings, which will break this case:
  m = SmilesToMol(smi,0,0); 
  int bfs = MolOps::findSSSR(*m);
  TEST_ASSERT(bfs==5);
  BOOST_LOG(rdInfoLog) << "BFSR: " << bfs << "\n";
  VECT_INT_VECT bfrs;
  bfrs.resize(0);
  bfs = MolOps::symmetrizeSSSR(*m, bfrs);
  TEST_ASSERT(bfs==6);
  BOOST_LOG(rdInfoLog) << "BFSR: " << bfs << "\n";
  //VECT_INT_VECT_I ri;
  //for (ri == bfrs.begin(); ri != bfrs.end(); ri++) {
  for (unsigned int ri = 0; ri < bfrs.size(); ri++) {
    INT_VECT_I mi;
    INT_VECT bring = bfrs[ri];
    BOOST_LOG(rdInfoLog) << "( ";
    //for (mi = (*ri).begin(); mi != (*ri).end(); mi++) {
    for (mi = bring.begin(); mi != bring.end(); mi++) {
      BOOST_LOG(rdInfoLog) << " " << (*mi);
    }
    BOOST_LOG(rdInfoLog) << ")\n";
  }
  BOOST_LOG(rdInfoLog) << smi << "\n";
  
  delete m;
  

  smi = "C1CC2C1CCC2";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  count = MolOps::findSSSR(*m,sssr);
  TEST_ASSERT(count==2);
  TEST_ASSERT(sssr[0].size()==4);
  TEST_ASSERT(sssr[1].size()==5);
  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;

  smi = "C12=C3C=CC=C1C=CC2=CC=C3";
  BOOST_LOG(rdInfoLog) << "\n" << smi << "\n";
  m = SmilesToMol(smi,0,0);
  TEST_ASSERT(m);
  count = MolOps::findSSSR(*m,sssr);
  TEST_ASSERT(count==3);
  TEST_ASSERT(sssr[0].size()==6);
  TEST_ASSERT(sssr[1].size()==5);
  TEST_ASSERT(sssr[2].size()==6);
  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;


  smi = "C1(O)C(O)C(O)C1O";
  m = SmilesToMol(smi,0,0);
  TEST_ASSERT(m);
  count = MolOps::findSSSR(*m,sssr);
  TEST_ASSERT(count==1);
  TEST_ASSERT(sssr[0].size()==4);
  for(unsigned i=0;i<m->getNumAtoms();i++){
    if(!(i%2)){
      TEST_ASSERT(m->getRingInfo()->numAtomRings(i)==1);
      TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(i,4));
    } else {
      TEST_ASSERT(m->getRingInfo()->numAtomRings(i)==0);
    }
  }
  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;

  // this molecule is from issue 134
  // it should come up with three rings
  smi = "SC(C3C1CC(C3)CC(C2S)(O)C1)2S";
  m = SmilesToMol(smi,0,0);
  TEST_ASSERT(m);
  count = MolOps::findSSSR(*m,sssr);
  TEST_ASSERT(count==3);
  TEST_ASSERT(sssr[0].size()==5);
  TEST_ASSERT(sssr[1].size()==6);
  TEST_ASSERT(sssr[2].size()==6);

  // this yet another painful case
  smi = "CC1=CC=C(C=C1)S(=O)(=O)O[CH]2[CH]3CO[CH](O3)[CH]4OC(C)(C)O[CH]24";
  m = SmilesToMol(smi,0,0);
  TEST_ASSERT(m);
  count = MolOps::findSSSR(*m,sssr);
  TEST_ASSERT(count==4);
  TEST_ASSERT(sssr[0].size()==6);
  TEST_ASSERT(sssr[1].size()==5);
  TEST_ASSERT(sssr[2].size()==5);
  TEST_ASSERT(sssr[3].size()==6);

  smi = "C1CC2C1C2";
  m = SmilesToMol(smi,0,0);
  TEST_ASSERT(m);
  count = MolOps::findSSSR(*m,sssr);
  TEST_ASSERT(count==2);
  TEST_ASSERT(sssr[0].size()==4);
  TEST_ASSERT(sssr[1].size()==3);

  TEST_ASSERT(m->getRingInfo()->numAtomRings(0)==1);
  TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(0,4));

  TEST_ASSERT(m->getRingInfo()->numAtomRings(1)==1);
  TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(1,4));

  TEST_ASSERT(m->getRingInfo()->numAtomRings(2)==2);
  TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(2,3));
  TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(2,4));

  TEST_ASSERT(m->getRingInfo()->numAtomRings(3)==2);
  TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(3,4));
  TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(3,3));

  TEST_ASSERT(m->getRingInfo()->numAtomRings(4)==1);
  TEST_ASSERT(!m->getRingInfo()->isAtomInRingOfSize(4,4));
  TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(4,3));
  delete m;

  // This is a test of Issue 217
  smi = "C=C1C2CC1C2";
  m = SmilesToMol(smi,0,0);
  TEST_ASSERT(m);
  count = MolOps::findSSSR(*m,sssr);
  TEST_ASSERT(count==2);
  TEST_ASSERT(sssr[0].size()==4);
  TEST_ASSERT(sssr[1].size()==4);
  count = MolOps::symmetrizeSSSR(*m,sssr);
  TEST_ASSERT(count==3);
  TEST_ASSERT(sssr[0].size()==4);
  TEST_ASSERT(sssr[1].size()==4);
  TEST_ASSERT(sssr[2].size()==4);

  TEST_ASSERT(m->getRingInfo()->numAtomRings(0)==0);
  TEST_ASSERT(m->getRingInfo()->numAtomRings(1)==2);
  TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(1,4));
  TEST_ASSERT(m->getRingInfo()->numAtomRings(2)==3);
  TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(2,4));
  TEST_ASSERT(m->getRingInfo()->numAtomRings(3)==2);
  TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(3,4));
  TEST_ASSERT(m->getRingInfo()->numAtomRings(4)==3);
  TEST_ASSERT(m->getRingInfo()->isAtomInRingOfSize(4,4));
  delete m;




}



void test4(){
  string smi;
  Mol *m;
  INT_VECT iv;
  VECT_INT_VECT sssr;
  smi = "CC";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  double *adjMat = MolOps::getAdjacencyMatrix(*m);
  TEST_ASSERT(adjMat);
  TEST_ASSERT(adjMat[0]==0);
  TEST_ASSERT(adjMat[1]==1);
  TEST_ASSERT(adjMat[2]==1);
  TEST_ASSERT(adjMat[3]==0);
  adjMat = MolOps::getAdjacencyMatrix(*m);
  TEST_ASSERT(adjMat);
  TEST_ASSERT(adjMat[0]==0);
  TEST_ASSERT(adjMat[1]==1);
  TEST_ASSERT(adjMat[2]==1);
  TEST_ASSERT(adjMat[3]==0);
}

void test5(){
  string smi;
  Mol *m;
  VECT_INT_VECT sssr;

  int count; 
  smi = "C1C4C5C3C(=O)C2C5C1C2C34";
  m = SmilesToMol(smi,0,0);
  count = MolOps::findSSSR(*m,sssr);
  BOOST_LOG(rdInfoLog) << "Count: " << count << "\n";
  CHECK_INVARIANT(count==5,"");

  smi = "C1C(C2)CCC2C1";
  m = SmilesToMol(smi);
  count = MolOps::findSSSR(*m,sssr);
  CHECK_INVARIANT(count==2,"");
}

/*
  void test6(){
  string smi;
  Mol *m;
  VECT_INT_VECT sssr;

  int c1,c2;
  smi = "C1(Cl)C(Cl)C1Cl";
  m = SmilesToMol(smi);
  INT_SET ringAtoms,ringBonds;
  //boost::tie(c1,c2) = MolOps::findRingAtomsAndBonds(*m,ringAtoms,ringBonds);
  
  CHECK_INVARIANT(c1==3,"bad nRingAtoms");
  CHECK_INVARIANT(ringAtoms.count(0)==1,"bad RingAtoms");
  CHECK_INVARIANT(ringAtoms.count(1)==0,"bad RingAtoms");
  CHECK_INVARIANT(ringAtoms.count(2)==1,"bad RingAtoms");
  CHECK_INVARIANT(ringAtoms.count(3)==0,"bad RingAtoms");
  CHECK_INVARIANT(ringAtoms.count(4)==1,"bad RingAtoms");
  CHECK_INVARIANT(ringAtoms.count(5)==0,"bad RingAtoms");
  
  CHECK_INVARIANT(c2==3,"bad nRingBonds");
  CHECK_INVARIANT(ringBonds.count(0)==0,"");
  CHECK_INVARIANT(ringBonds.count(1)==1,"");
  CHECK_INVARIANT(ringBonds.count(2)==0,"");
  CHECK_INVARIANT(ringBonds.count(3)==1,"");
  CHECK_INVARIANT(ringBonds.count(4)==0,"");
  CHECK_INVARIANT(ringBonds.count(5)==1,"");

  
  }
*/

void test7(){
#if 0
  string smi;
  Mol *m;
  INT_VECT tree;
#if 1
  smi = "C(CO)OCC";
  m = SmilesToMol(smi);
  MolOps::findSpanningTree(*m,tree);
  CHECK_INVARIANT(tree.size()==5,"bad mst");
  delete m;

  smi = "C1CC1";
  m = SmilesToMol(smi);
  MolOps::findSpanningTree(*m,tree);
  CHECK_INVARIANT(tree.size()==2,"bad mst");
  delete m;

  smi = "C1C=C1";
  m = SmilesToMol(smi);
  MolOps::findSpanningTree(*m,tree);
  CHECK_INVARIANT(tree.size()==2,"bad mst");
  CHECK_INVARIANT(std::find(tree.begin(),tree.end(),1)==tree.end(),"bogus idx in mst");
  delete m;
#endif
  
  smi = "C1C=CC=CC=1";
  m = SmilesToMol(smi);
  MolOps::findSpanningTree(*m,tree);
  CHECK_INVARIANT(tree.size()==5,"bad mst");
  delete m;


  smi = "C1C(=CC1)";
  m = SmilesToMol(smi);
  MolOps::findSpanningTree(*m,tree);
  CHECK_INVARIANT(tree.size()==3,"bad mst");
  delete m;


  smi = "C1C(C=C1)";
  m = SmilesToMol(smi);
  MolOps::findSpanningTree(*m,tree);
  CHECK_INVARIANT(tree.size()==3,"bad mst");
  delete m;

  smi = "C1C(C2)CCC2C1";
  m = SmilesToMol(smi);
  MolOps::findSpanningTree(*m,tree);
  CHECK_INVARIANT(tree.size()==6,"bad mst");
  delete m;

  smi = "C1C2CC3CCCCC3CC2CCC1";
  m = SmilesToMol(smi);
  MolOps::findSpanningTree(*m,tree);
  CHECK_INVARIANT(tree.size()==13,"bad mst");
  delete m;
#endif
}

void test8()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Hydrogen Ops" << std::endl;
  ROMol *m,*m2,*m3;
  INT_VECT tree;

  std::string smi = "CCC";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m,"");
  CHECK_INVARIANT(m->getNumAtoms()==3,"");

  //BOOST_LOG(rdInfoLog) << "1" << std::endl;
  m2 = MolOps::addHs(*m);
  CHECK_INVARIANT(m2->getNumAtoms()==11,"");

  smi = "CC(=O)[OH]";
  delete m;
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m,"");
  CHECK_INVARIANT(m->getNumAtoms()==4,"");

  //BOOST_LOG(rdInfoLog) << "2" << std::endl;
  delete m2;
  m2 = MolOps::addHs(*m,true);
  CHECK_INVARIANT(m2->getNumAtoms()==5,"");

  
  //BOOST_LOG(rdInfoLog) << "3" << std::endl;
  m3 = MolOps::addHs(*m2,false);
  CHECK_INVARIANT(m3->getNumAtoms()==8,"");

  //BOOST_LOG(rdInfoLog) << "4" << std::endl;
  delete m2;
  m2 = MolOps::addHs(*m,false);
  CHECK_INVARIANT(m2->getNumAtoms()==8,"");
  delete m3;
  // remove all
  //BOOST_LOG(rdInfoLog) << "5" << std::endl;
  m3 = MolOps::removeHs(*m2,false);
  CHECK_INVARIANT(m3->getNumAtoms()==4,"");
  delete m3;
  // remove only implicit
  //BOOST_LOG(rdInfoLog) << "6" << std::endl;
  m3 = MolOps::removeHs(*m2,true);
  CHECK_INVARIANT(m3->getNumAtoms()==5,"");

  //BOOST_LOG(rdInfoLog) << "7" << std::endl;
  // remove all after removing only implicit
  MolOps::removeHs(static_cast<RWMol &>(*m3),false);
  CHECK_INVARIANT(m3->getNumAtoms()==4,"");

  // this test is also done in the same order in the python tests:
  delete m;
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m,"");
  CHECK_INVARIANT(m->getNumAtoms()==4,"");

  delete m2;
  m2 = MolOps::addHs(*m,true);
  CHECK_INVARIANT(m2->getNumAtoms()==5,"");
  //BOOST_LOG(rdInfoLog) << "8" << std::endl;
  m3 = MolOps::removeHs(*m2,true);
  CHECK_INVARIANT(m3->getNumAtoms()==5,"");
  delete m3;
  //BOOST_LOG(rdInfoLog) << "9" << std::endl;
  m3 = MolOps::removeHs(*m2,false);
  CHECK_INVARIANT(m3->getNumAtoms()==4,"");

  delete m2;
  //BOOST_LOG(rdInfoLog) << "10" << std::endl;
  m2 = MolOps::addHs(*m,false);
  CHECK_INVARIANT(m2->getNumAtoms()==8,"");
  delete m3;
  //BOOST_LOG(rdInfoLog) << "11" << std::endl;
  m3 = MolOps::removeHs(*m2,true);
  CHECK_INVARIANT(m3->getNumAtoms()==5,"");
  delete m3;
  //BOOST_LOG(rdInfoLog) << "12" << std::endl;
  m3 = MolOps::removeHs(*m2,false);
  CHECK_INVARIANT(m3->getNumAtoms()==4,"");


  // related to RDTrack Issues 109 and 110:
  smi = "C1C=C([C@H](N)C(=O)N[C@@]2([H])[C@]3([H])SC(C)(C)[C@@H](C(=O)O)N3C(=O)2)C=CC=1";
  delete m;
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m,"");
  CHECK_INVARIANT(m->getNumAtoms()==24,"");
  delete m3;
  //BOOST_LOG(rdInfoLog) << "13" << std::endl;
  m3 = MolOps::removeHs(*m,false);
  CHECK_INVARIANT(m3->getNumAtoms()==24,"");


  // RDTrack Issue 130:
  delete m;
  smi = "[H][N+]([H])([H])[H]";
  m = SmilesToMol(smi,false,false);
  CHECK_INVARIANT(m,"");
  CHECK_INVARIANT(m->getNumAtoms()==5,"");
  delete m2;
  //BOOST_LOG(rdInfoLog) << "14" << std::endl;
  m2 = MolOps::removeHs(*m,0,false);
  CHECK_INVARIANT(m2->getNumAtoms()==1,"");
  delete m;
  smi = "[H][N+]([H])([H])[H]";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m,"");
  CHECK_INVARIANT(m->getNumAtoms()==1,"");


  delete m;
  smi = "[H][H]";
  m = SmilesToMol(smi,false,false);
  CHECK_INVARIANT(m,"");
  CHECK_INVARIANT(m->getNumAtoms()==2,"");
  delete m2;
  //BOOST_LOG(rdInfoLog) << "15" << std::endl;
  m2 = MolOps::removeHs(*m,0,false);
  CHECK_INVARIANT(m2->getNumAtoms()==2,"");

  std::string sma;
  delete m;
  smi = "CC";
  m = SmartsToMol(smi);
  MolOps::sanitizeMol(*((RWMol *)m));
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==2);
  sma = SmartsWrite::GetAtomSmarts(static_cast<const QueryAtom *>(m->getAtomWithIdx(0)));
  TEST_ASSERT(sma=="C");

  delete m2;
  //BOOST_LOG(rdInfoLog) << "16" << std::endl;
  m2 = MolOps::addHs(*m);
  TEST_ASSERT(m2->getNumAtoms()==8);
  sma = SmartsWrite::GetAtomSmarts(static_cast<const QueryAtom *>(m2->getAtomWithIdx(0)));
  TEST_ASSERT(sma=="C");

  delete m;
  //BOOST_LOG(rdInfoLog) << "17" << std::endl;
  m = MolOps::mergeQueryHs(*m2);
  TEST_ASSERT(m->getNumAtoms()==2);
  sma = SmartsWrite::GetAtomSmarts(static_cast<const QueryAtom *>(m->getAtomWithIdx(0)));
  //BOOST_LOG(rdInfoLog) << "sma: " << sma<<std::endl;
  // this was sf.net issue 3415204:
  TEST_ASSERT(sma=="[C&!H0&!H1&!H2]");

  // RDTrack Issue 1228:
  delete m;
  smi = "c1c[nH]cc1";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m,"");
  CHECK_INVARIANT(m->getNumAtoms()==5,"");
  delete m2;
  //BOOST_LOG(rdInfoLog) << "18" << std::endl;
  m2 = MolOps::addHs(*m,false,false);
  CHECK_INVARIANT(m2->getNumAtoms()==10,"");
  delete m;
  //BOOST_LOG(rdInfoLog) << "19" << std::endl;
  m = MolOps::removeHs(*m2);
  CHECK_INVARIANT(m,"");
  CHECK_INVARIANT(m->getNumAtoms()==5,"");


  // labelling:
  delete m;
  smi = "c1cn([H])cc1";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m,"");
  CHECK_INVARIANT(m->getNumAtoms()==5,"");
  delete m2;
  //BOOST_LOG(rdInfoLog) << "19" << std::endl;
  m2 = MolOps::removeHs(*m);
  CHECK_INVARIANT(m,"");
  CHECK_INVARIANT(m->getNumAtoms()==5,"");

  delete m;
  smi = "c1cn([2H])cc1";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m,"");
  CHECK_INVARIANT(m->getNumAtoms()==6,"");
  delete m2;
  //BOOST_LOG(rdInfoLog) << "19" << std::endl;
  m2 = MolOps::removeHs(*m);
  CHECK_INVARIANT(m,"");
  CHECK_INVARIANT(m->getNumAtoms()==6,"");
  delete m;
  delete m2;

  smi = "CC[H]";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==3);
  m2 = MolOps::mergeQueryHs(*m);
  TEST_ASSERT(m2->getNumAtoms()==2);
  TEST_ASSERT(!m->getAtomWithIdx(1)->hasQuery());
  TEST_ASSERT(m2->getAtomWithIdx(1)->hasQuery());
  delete m;
  delete m2;

  // sf.net issue 3415206
  smi = "CO[H]";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==3);
  m2 = MolOps::mergeQueryHs(*m);
  TEST_ASSERT(m2->getNumAtoms()==2);
  sma = SmartsWrite::GetAtomSmarts(static_cast<const QueryAtom *>(m2->getAtomWithIdx(1)));
  //BOOST_LOG(rdInfoLog) << "sma: " << sma<<std::endl;
  TEST_ASSERT(sma=="[#8&!H0]");
  delete m;
  delete m2;

  smi = "CN([H])[H]";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  m2 = MolOps::mergeQueryHs(*m);
  TEST_ASSERT(m2->getNumAtoms()==2);
  sma = SmartsWrite::GetAtomSmarts(static_cast<const QueryAtom *>(m2->getAtomWithIdx(1)));
  //BOOST_LOG(rdInfoLog) << "sma: " << sma<<std::endl;
  TEST_ASSERT(sma=="[#7&!H0&!H1]");
  delete m;
  delete m2;

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}


void test9()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Distance Matrix Operations" << std::endl;
  ROMol *m;
  std::string smi = "CC=C";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==3);
  double *dMat;
  dMat = MolOps::getDistanceMat(*m,false,false);
  TEST_ASSERT(dMat);
  TEST_ASSERT(dMat[0]==0.0);
  TEST_ASSERT(dMat[1]==1.0);
  TEST_ASSERT(dMat[2]==2.0);
  TEST_ASSERT(dMat[3]==1.0);
  TEST_ASSERT(dMat[4]==0.0);
  TEST_ASSERT(dMat[5]==1.0);
  TEST_ASSERT(dMat[6]==2.0);
  TEST_ASSERT(dMat[7]==1.0);
  TEST_ASSERT(dMat[8]==0.0);

  dMat = MolOps::getDistanceMat(*m,false,false);
  TEST_ASSERT(dMat);
  TEST_ASSERT(dMat[0]==0.0);
  TEST_ASSERT(dMat[1]==1.0);
  TEST_ASSERT(dMat[2]==2.0);
  TEST_ASSERT(dMat[3]==1.0);
  TEST_ASSERT(dMat[4]==0.0);
  TEST_ASSERT(dMat[5]==1.0);
  TEST_ASSERT(dMat[6]==2.0);
  TEST_ASSERT(dMat[7]==1.0);
  TEST_ASSERT(dMat[8]==0.0);

  // test Issue328:
  dMat = MolOps::getDistanceMat(*m,true,false);
  TEST_ASSERT(dMat);
  TEST_ASSERT(dMat[0]==0.0);
  TEST_ASSERT(dMat[1]==1.0);
  TEST_ASSERT(dMat[2]==1.5);
  TEST_ASSERT(dMat[3]==1.0);
  TEST_ASSERT(dMat[4]==0.0);
  TEST_ASSERT(dMat[5]==0.5);
  TEST_ASSERT(dMat[6]==1.5);
  TEST_ASSERT(dMat[7]==0.5);
  TEST_ASSERT(dMat[8]==0.0);


  dMat = MolOps::getDistanceMat(*m,false,false);
  TEST_ASSERT(dMat);
  TEST_ASSERT(dMat[0]==0.0);
  TEST_ASSERT(dMat[1]==1.0);
  TEST_ASSERT(dMat[2]==2.0);
  TEST_ASSERT(dMat[3]==1.0);
  TEST_ASSERT(dMat[4]==0.0);
  TEST_ASSERT(dMat[5]==1.0);
  TEST_ASSERT(dMat[6]==2.0);
  TEST_ASSERT(dMat[7]==1.0);
  TEST_ASSERT(dMat[8]==0.0);
  
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test10()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Atom Ranking" << std::endl;
  ROMol *m;
  std::string smi = "FC(Cl)(Br)C";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);

  INT_VECT ranks;
  ranks.resize(m->getNumAtoms());
  Chirality::assignAtomCIPRanks(*m,ranks);

  int cip1,cip2;
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPRank"));
  m->getAtomWithIdx(0)->getProp("_CIPRank",cip1);
  TEST_ASSERT(cip1==ranks[0]);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPRank"));
  m->getAtomWithIdx(2)->getProp("_CIPRank",cip2);
  TEST_ASSERT(cip2==ranks[2]);
  TEST_ASSERT(cip1<cip2);
  TEST_ASSERT(m->getAtomWithIdx(4)->hasProp("_CIPRank"));
  m->getAtomWithIdx(4)->getProp("_CIPRank",cip2);
  TEST_ASSERT(cip1>cip2);
  m->getAtomWithIdx(2)->getProp("_CIPRank",cip1);
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPRank"));
  m->getAtomWithIdx(3)->getProp("_CIPRank",cip2);
  TEST_ASSERT(cip1<cip2);
  
  delete m;
  smi = "FC(Cl)(Br)C(F)(F)F";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  ranks.resize(m->getNumAtoms());
  Chirality::assignAtomCIPRanks(*m,ranks);
  for(unsigned int i=0;i<m->getNumAtoms();i++){
    int cip;
    TEST_ASSERT(m->getAtomWithIdx(i)->hasProp("_CIPRank"));
    m->getAtomWithIdx(i)->getProp("_CIPRank",cip);
  }
  
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}


void test11()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing CIP chirality" << std::endl;
  ROMol *m;
  std::string cip;
  std::string smi = "F[C@]([C@])(Cl)Br";

  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  // make sure the cleanup worked:
  TEST_ASSERT(m->getAtomWithIdx(2)->getChiralTag()==Atom::CHI_UNSPECIFIED);

  TEST_ASSERT(!(m->getAtomWithIdx(0)->hasProp("_CIPCode")));
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  TEST_ASSERT(!(m->getAtomWithIdx(2)->hasProp("_CIPCode")));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");


  delete m;
  smi = "F[C@H](C)C";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(!(m->getAtomWithIdx(0)->hasProp("_CIPCode")));
  TEST_ASSERT(!(m->getAtomWithIdx(1)->hasProp("_CIPCode")));
  TEST_ASSERT(!(m->getAtomWithIdx(2)->hasProp("_CIPCode")));
  // test Issue 194:
  TEST_ASSERT(m->getAtomWithIdx(1)->getNumExplicitHs()==0);


  delete m;
  smi = "F[C@]1CC(Cl)C1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(!(m->getAtomWithIdx(1)->hasProp("_CIPCode")));

  delete m;
  smi = "F[C@H]1C(Cl)CC1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));

  delete m;
  smi = "F[C@@](C)(Cl)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(!(m->getAtomWithIdx(0)->hasProp("_CIPCode")));
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  TEST_ASSERT(!(m->getAtomWithIdx(2)->hasProp("_CIPCode")));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete m;
  smi = "F[C@](Br)(C)Cl";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  delete m;
  smi = "F[C@](Cl)(Br)C";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete m;
  smi = "FC(F)(F)[C@](Br)(F)C(Cl)(Cl)Cl";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(4)->hasProp("_CIPCode"));
  m->getAtomWithIdx(4)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete m;
  smi = "C[C@](C=C)(F)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete m;
  smi = "CC[C@](C=C)(F)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");
  delete m;
  smi = "[CH2-][C@](C)(F)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete m;
  smi = "F[C@]([H])(Cl)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete m;
  smi = "[C@H](Cl)(F)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete m;
  smi = "[C@]([H])(Cl)(F)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete m;
  smi = "F[C@H](Cl)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  
  delete m;
  smi = "CC[C@H](C=C)C";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete m;
  smi = "OC[C@H](C=C)C";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete m;
  smi = "CC[C@H](C=C)O";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete m;
  smi = "OC[C@H](C=C)O";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete m;
  smi = "C[C@H]1C[C@H](C=C1)N";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
  m->getAtomWithIdx(3)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  // a couple random molecules from the BBB data set
  delete m;
  smi = "OC[C@H]1C[C@@H](N2C=NC3=C2N=C(N)N=C3NC4CC4)C=C1";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");
  TEST_ASSERT(m->getAtomWithIdx(4)->hasProp("_CIPCode"));
  m->getAtomWithIdx(4)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");



  delete m;
  smi = "N[C@H]1O[C@@H](SC1)CO";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
  m->getAtomWithIdx(3)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");


  delete m;
  smi = "C1(N([C@H]2O[C@H](CO)SC2)C=CC(N)=N1)=O";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");
  TEST_ASSERT(m->getAtomWithIdx(4)->hasProp("_CIPCode"));
  m->getAtomWithIdx(4)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  // this is Issue 152:
  smi = "C1[C@H](N)C[C@H](C)C=1";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  TEST_ASSERT(m->getAtomWithIdx(4)->hasProp("_CIPCode"));
  m->getAtomWithIdx(4)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  
  // -----------------------------------------------
  // these are related to Issue 397:
  smi = "C(=O)[C@@H](C)N";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  smi = "C(=O)[C@@H](CO)N";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  smi = "C(O)[C@@H](C)N";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  // -----------------------------------------------


  // NOTE: This test gives correct results according to the current
  // CIP ranking procedure, but the results are actually incorrect.
  // This arises because of the inclusion of hybridization in the
  // chiral atom invariants
  // (see the note in Chirality.cpp:buildCIPInvariants())
  smi = "[H][C@@](O)(C=C)C(C)CC";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  smi = "[H][C@@](O)(C=C)C(C)CO";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  smi = "[H][C@@]12C[C@@](NC1)(OC2)[H]";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  smi = "[H][C@@]12C[C@@](C=C1)(CC2)[H]";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  smi = "[H][C@@]12O[C@@](CC1)(C3C2C(NC3=O)=O)[H]";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  smi = "[H][C@@]12O[C@@](C=C1)(C3C2C(NC3=O)=O)[H]";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");


  // -----------------------------------------------

  
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}


void test12()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing double bond stereochemistry" << std::endl;
  ROMol *m;
  RWMol *m2;
  std::string smi = "F\\C=C/Cl";
  std::string refSmi;

  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
  TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
  TEST_ASSERT(m->getBondWithIdx(2)->getStereo() == Bond::STEREONONE);

  delete m;
  smi = "F/C=CCl";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

  delete m;
  smi = "F/C=C/Cl";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getBondWithIdx(1)->getStereo()==Bond::STEREOE);

  delete m;
  smi = "F/C=C(/Br)Cl";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getBondWithIdx(1)->getStereo()==Bond::STEREOE);


  delete m;
  smi = "F/C=C(/Cl)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getBondWithIdx(1)->getStereo()==Bond::STEREOZ);


  delete m;
  smi = "F/C(Br)=C/Cl";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getBondWithIdx(2)->getStereo()==Bond::STEREOZ);

  delete m;
  smi = "F/C=C(/Cl)Cl";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

  // build a molecule from scratch to test problems
  // around Issue 180. The molecule corresponds to SMILES
  // F/C=C(/Br)C
  delete m;
  m2 = new RWMol();
  m2->addAtom(new Atom(9),true,true);
  m2->addAtom(new Atom(6),true,true);
  m2->addAtom(new Atom(6),true,true);
  m2->addAtom(new Atom(35),true,true);
  m2->addAtom(new Atom(6),true,true);
  m2->addBond(0,1,Bond::SINGLE);
  m2->addBond(1,2,Bond::DOUBLE);
  m2->addBond(2,3,Bond::SINGLE);
  m2->addBond(2,4,Bond::SINGLE);
  m2->getBondWithIdx(0)->setBondDir(Bond::ENDUPRIGHT);
  m2->getBondWithIdx(2)->setBondDir(Bond::ENDUPRIGHT);
  m2->getBondWithIdx(3)->setBondDir(Bond::ENDDOWNRIGHT);
  MolOps::sanitizeMol(*m2);
  MolOps::assignStereochemistry(*m2);
  TEST_ASSERT(m2->getBondWithIdx(1)->getStereo()==Bond::STEREOE);


  m2->getBondWithIdx(0)->setBondDir(Bond::ENDDOWNRIGHT);
  MolOps::assignStereochemistry(*m2,true,true);
  TEST_ASSERT(m2->getBondWithIdx(1)->getStereo()==Bond::STEREOZ);

  delete m2;
  m2 = new RWMol();
  m2->addAtom(new Atom(9),true,true);
  m2->addAtom(new Atom(6),true,true);
  m2->addAtom(new Atom(6),true,true);
  m2->addAtom(new Atom(35),true,true);
  m2->addAtom(new Atom(6),true,true);
  m2->addBond(1,0,Bond::SINGLE);
  m2->addBond(1,2,Bond::DOUBLE);
  m2->addBond(2,3,Bond::SINGLE);
  m2->addBond(2,4,Bond::SINGLE);
  m2->getBondWithIdx(0)->setBondDir(Bond::ENDUPRIGHT);
  m2->getBondWithIdx(2)->setBondDir(Bond::ENDUPRIGHT);
  m2->getBondWithIdx(3)->setBondDir(Bond::ENDDOWNRIGHT);
  MolOps::sanitizeMol(*m2);
  MolOps::assignStereochemistry(*m2);
  TEST_ASSERT(m2->getBondWithIdx(1)->getStereo()==Bond::STEREOZ);

  m2->getBondWithIdx(0)->setBondDir(Bond::ENDDOWNRIGHT);
  MolOps::assignStereochemistry(*m2,true,true);
  TEST_ASSERT(m2->getBondWithIdx(1)->getStereo()==Bond::STEREOE);


  // ----------------------
  // test Issue 174:
  delete m2;
  smi = "O\\N=C\\C=N/O";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(m2->getBondWithIdx(1)->getStereo()==Bond::STEREOE);
  TEST_ASSERT(m2->getBondWithIdx(3)->getStereo()==Bond::STEREOZ);
  refSmi = MolToSmiles(*m2,1);
  m = SmilesToMol(refSmi);
  TEST_ASSERT(m);
  smi = MolToSmiles(*m,1);
  TEST_ASSERT(refSmi==smi);

  delete m;
  delete m2;

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}


void testIssue183()
{
  // ----------------------
  // test "unsetting" of redundant bond directions:

  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Issue 183\n" << std::endl;
  RWMol *m,*m2;
  std::string smi;
  std::string refSmi;

  smi = "Cl\\C(C)=C(\\C(F)=C(/F)C)/C(C)=C(\\F)C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(m2->getBondWithIdx(2)->getStereo()==Bond::STEREOE);
  TEST_ASSERT(m2->getBondWithIdx(5)->getStereo()==Bond::STEREOE);
  TEST_ASSERT(m2->getBondWithIdx(10)->getStereo()==Bond::STEREOZ);

  m2->debugMol(std::cerr);  
  refSmi = MolToSmiles(*m2,1);
  BOOST_LOG(rdInfoLog) << "ref: " << refSmi << std::endl;
  m = SmilesToMol(refSmi);
  TEST_ASSERT(m);
  smi = MolToSmiles(*m,1);
  BOOST_LOG(rdInfoLog) << "smi: " << smi << std::endl;
  TEST_ASSERT(refSmi==smi);

  int nEs=0,nZs=0,nDbl=0;
  for(RWMol::BondIterator bondIt=m->beginBonds();
      bondIt!=m->endBonds();
      bondIt++){
    if((*bondIt)->getBondType()==Bond::DOUBLE){
      nDbl++;
      if((*bondIt)->getStereo()==Bond::STEREOE) nEs++;
      else if((*bondIt)->getStereo()==Bond::STEREOZ) nZs++;
    }
  }
  //BOOST_LOG(rdInfoLog) << ">> " << nDbl << " " << nEs << " " << nZs << std::endl;
  TEST_ASSERT(nDbl==3);
  TEST_ASSERT(nEs==2);
  TEST_ASSERT(nZs==1);
  delete m;
  delete m2;
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}


void testIssue188()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Issue 188: bad CIP rankings" << std::endl;
  ROMol *m;
  std::string smi;
  int cip1,cip2,cip3;

  smi = "OC[C@H](C=C)C";
  m = SmilesToMol(smi);
  INT_VECT ranks;
  ranks.resize(m->getNumAtoms());
  Chirality::assignAtomCIPRanks(*m,ranks);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPRank"));
  m->getAtomWithIdx(1)->getProp("_CIPRank",cip1);
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPRank"));
  m->getAtomWithIdx(3)->getProp("_CIPRank",cip2);
  TEST_ASSERT(cip1>cip2);
  TEST_ASSERT(m->getAtomWithIdx(5)->hasProp("_CIPRank"));
  m->getAtomWithIdx(5)->getProp("_CIPRank",cip3);
  TEST_ASSERT(cip1>cip3);
  TEST_ASSERT(cip2>cip3);

  delete m;
  smi = "CC(=N\\N)/C=N/N";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  ranks.resize(m->getNumAtoms());
  Chirality::assignAtomCIPRanks(*m,ranks);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPRank"));
  m->getAtomWithIdx(0)->getProp("_CIPRank",cip1);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPRank"));
  m->getAtomWithIdx(1)->getProp("_CIPRank",cip2);
  TEST_ASSERT(cip2>cip1);
  TEST_ASSERT(m->getAtomWithIdx(4)->hasProp("_CIPRank"));
  m->getAtomWithIdx(4)->getProp("_CIPRank",cip3);
  TEST_ASSERT(cip3>cip1);
  TEST_ASSERT(cip2>cip3);
  delete m;
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}



void testIssue189()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Issue 189: BondDirs not getting properly cleared." << std::endl;
  ROMol *m;
  std::string smi,refSmi;
  int count;

  smi = "C(=S)/N=c(/n1C)scc1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);

  TEST_ASSERT(m->getBondWithIdx(2)->getBondType()==Bond::DOUBLE);
  TEST_ASSERT(m->getBondWithIdx(2)->getStereo()==Bond::STEREOZ);

  count = 0;
  for(unsigned int i=0;i<m->getNumBonds();i++){
    if(m->getBondWithIdx(i)->getBondDir()!=Bond::NONE){
      count++;
    }
  }
  TEST_ASSERT(count==2);

  refSmi = MolToSmiles(*m,1);
  count = 0;
  for(unsigned int i=0;i<m->getNumBonds();i++){
    if(m->getBondWithIdx(i)->getBondDir()!=Bond::NONE){
      count++;
    }
  }
  TEST_ASSERT(count==2);

  delete m;
  m = SmilesToMol(refSmi);
  TEST_ASSERT(m);
  count = 0;
  for(unsigned int i=0;i<m->getNumBonds();i++){
    if(m->getBondWithIdx(i)->getBondDir()!=Bond::NONE){
      count++;
    }
  }
  TEST_ASSERT(count==2);
  smi = MolToSmiles(*m,1);
  count = 0;
  for(unsigned int i=0;i<m->getNumBonds();i++){
    if(m->getBondWithIdx(i)->getBondDir()!=Bond::NONE){
      count++;
    }
  }
  TEST_ASSERT(count==2);
  TEST_ASSERT(smi==refSmi);

  
  delete m;
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testIssue190()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Issue 190: BondDirs incorrectly cleared." << std::endl;
  ROMol *m;
  std::string smi,refSmi;
  int count;

  smi ="O\\N=C\\NC(\\C)=N/OC";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  TEST_ASSERT(m->getBondWithIdx(5)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(m->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
  refSmi = MolToSmiles(*m,1);

  count = 0;
  for(unsigned int i=0;i<m->getNumBonds();i++){
    if(m->getBondWithIdx(i)->getBondDir()!=Bond::NONE){
      count++;
    }
  }
  TEST_ASSERT(count==4);

  delete m;
  m = SmilesToMol(refSmi);
  TEST_ASSERT(m);
  count = 0;
  for(unsigned int i=0;i<m->getNumBonds();i++){
    if(m->getBondWithIdx(i)->getBondDir()!=Bond::NONE){
      count++;
    }
  }
  TEST_ASSERT(count==4);
  smi = MolToSmiles(*m,1);
  count = 0;
  for(unsigned int i=0;i<m->getNumBonds();i++){
    if(m->getBondWithIdx(i)->getBondDir()!=Bond::NONE){
      count++;
    }
  }
  TEST_ASSERT(count==4);
  TEST_ASSERT(smi==refSmi);



  delete m;
  smi ="O\\N=C\\CC(\\C)=N/OC";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  TEST_ASSERT(m->getBondWithIdx(5)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(m->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
  refSmi = MolToSmiles(*m,1);

  count = 0;
  for(unsigned int i=0;i<m->getNumBonds();i++){
    if(m->getBondWithIdx(i)->getBondDir()!=Bond::NONE){
      count++;
    }
  }
  TEST_ASSERT(count==4);
  delete m;
  m = SmilesToMol(refSmi);
  TEST_ASSERT(m);
  count = 0;
  for(unsigned int i=0;i<m->getNumBonds();i++){
    if(m->getBondWithIdx(i)->getBondDir()!=Bond::NONE){
      count++;
    }
  }
  TEST_ASSERT(count==4);
  smi = MolToSmiles(*m,1);
  count = 0;
  for(unsigned int i=0;i<m->getNumBonds();i++){
    if(m->getBondWithIdx(i)->getBondDir()!=Bond::NONE){
      count++;
    }
  }
  TEST_ASSERT(count==4);
  TEST_ASSERT(smi==refSmi);



  delete m;
  smi ="O\\N=C\\C(=O)C(\\C)=N/OC";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  TEST_ASSERT(m->getBondWithIdx(6)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(m->getBondWithIdx(6)->getStereo() == Bond::STEREOZ);
  refSmi = MolToSmiles(*m,1);

  count = 0;
  for(unsigned int i=0;i<m->getNumBonds();i++){
    if(m->getBondWithIdx(i)->getBondDir()!=Bond::NONE){
      count++;
    }
  }
  TEST_ASSERT(count==4);
  delete m;
  m = SmilesToMol(refSmi);
  TEST_ASSERT(m);
  count = 0;
  for(unsigned int i=0;i<m->getNumBonds();i++){
    if(m->getBondWithIdx(i)->getBondDir()!=Bond::NONE){
      count++;
    }
  }
  TEST_ASSERT(count==4);
  smi = MolToSmiles(*m,1);
  count = 0;
  for(unsigned int i=0;i<m->getNumBonds();i++){
    if(m->getBondWithIdx(i)->getBondDir()!=Bond::NONE){
      count++;
    }
  }
  TEST_ASSERT(count==4);
  TEST_ASSERT(smi==refSmi);

  
  delete m;
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testShortestPath() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing shortest path code. This should finish very quickly." << std::endl;
  {
    std::string smi ="CC(OC1C(CCCC3)C3C(CCCC2)C2C1OC(C)=O)=O";
    ROMol *m = SmilesToMol(smi);

    INT_LIST path = MolOps::getShortestPath(*m, 1, 20);
    CHECK_INVARIANT(path.size() == 7, "");
    INT_LIST_CI pi = path.begin();
    CHECK_INVARIANT((*pi) == 1, ""); pi++;
    CHECK_INVARIANT((*pi) == 2, ""); pi++;
    CHECK_INVARIANT((*pi) == 3, ""); pi++;
    CHECK_INVARIANT((*pi) == 16, ""); pi++;
    CHECK_INVARIANT((*pi) == 17, ""); pi++;
    CHECK_INVARIANT((*pi) == 18, ""); pi++;
    CHECK_INVARIANT((*pi) == 20, ""); pi++;
    delete m;
  }
  {
    // issue 2219400
    std::string smi ="CC.C";
    ROMol *m = SmilesToMol(smi);

    INT_LIST path = MolOps::getShortestPath(*m, 0, 1);
    std::cerr<<"path: "<<path.size()<<std::endl;
    CHECK_INVARIANT(path.size() == 2, "");
    INT_LIST_CI pi = path.begin();
    CHECK_INVARIANT((*pi) == 0, ""); pi++;
    CHECK_INVARIANT((*pi) == 1, "");

    path = MolOps::getShortestPath(*m, 1, 2);
    CHECK_INVARIANT(path.size() == 0, "");

    path = MolOps::getShortestPath(*m, 0, 2);
    CHECK_INVARIANT(path.size() == 0, "");
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}


void testIssue210()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Issue 210" << std::endl;
  ROMol *m,*m2;

  std::string smi = "C1CC1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==3);
  TEST_ASSERT(m->getRingInfo()->isInitialized());

  m2 = MolOps::addHs(*m);
  TEST_ASSERT(m2->getNumAtoms()==9);
  TEST_ASSERT(m2->getRingInfo()->isInitialized());

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
  delete m;
  delete m2;
}


void testIssue211()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Issue 211" << std::endl;
  ROMol *m;

  std::string smi = "P(c1ccccc1)(c1ccccc1)c1ccccc1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==19);

  const Atom *at=m->getAtomWithIdx(0);
  TEST_ASSERT(at->getHybridization()==Atom::SP3);

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testIssue212()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Issue 212" << std::endl;
  ROMol *m,*m2;
  std::string smi,mb;
  smi = "C";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==1);
  Conformer *conf = new Conformer(1);
  m->addConformer(conf);
  conf->setAtomPos(0,RDGeom::Point3D(0,0,0));
  m2 = MolOps::addHs(*m,false,true);
  TEST_ASSERT(m2->getNumAtoms()==5);

  try{
    mb = MolToMolBlock(*m2);
  } catch (...) {
    TEST_ASSERT(0); //,"MolToMolBlock() failed");
  }

  delete m;
  delete m2;
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testAddHsCoords()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing AddHs with coordinates" << std::endl;
  ROMol *m,*m2;
  RDGeom::Point3D v;
  double bondLength = PeriodicTable::getTable()->getRb0(1) +
    PeriodicTable::getTable()->getRb0(6);
  double tetDist=2.*sin((109.471/2.)*M_PI/180)*bondLength;
  double sp2Dist=2.*sin(60.*M_PI/180)*bondLength;

  std::string smi;

  smi = "C";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==1);
  Conformer *conf = new Conformer(1);
  m->addConformer(conf);
  conf->setAtomPos(0,RDGeom::Point3D(0,0,0));

  m2 = MolOps::addHs(*m,false,true);
  const Conformer *conf2 = &(m2->getConformer());
  TEST_ASSERT(m2->getNumAtoms()==5);
  TEST_ASSERT(feq((conf2->getAtomPos(0) - conf2->getAtomPos(1)).length(),
                  bondLength));
  TEST_ASSERT(feq((conf2->getAtomPos(0) - conf2->getAtomPos(2)).length(),
                  bondLength));
  TEST_ASSERT(feq((conf2->getAtomPos(0) - conf2->getAtomPos(3)).length(),
                  bondLength));
  TEST_ASSERT(feq((conf2->getAtomPos(0) - conf2->getAtomPos(4)).length(),
                  bondLength));
  TEST_ASSERT(feq((conf2->getAtomPos(1) - conf2->getAtomPos(2)).length(),
                  tetDist));
  TEST_ASSERT(feq((conf2->getAtomPos(1) - conf2->getAtomPos(3)).length(),
                  tetDist));
  TEST_ASSERT(feq((conf2->getAtomPos(1) - conf2->getAtomPos(4)).length(),
                  tetDist));
  TEST_ASSERT(feq((conf2->getAtomPos(2) - conf2->getAtomPos(3)).length(),
                  tetDist));
  TEST_ASSERT(feq((conf2->getAtomPos(2) - conf2->getAtomPos(4)).length(),
                  tetDist));
  TEST_ASSERT(feq((conf2->getAtomPos(3) - conf2->getAtomPos(4)).length(),
                  tetDist));          
  delete m;
  delete m2;
  
  smi = "CC";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==2);
  conf = new Conformer(2);
  m->addConformer(conf);
  conf->setAtomPos(0,RDGeom::Point3D(0,0,0));
  conf->setAtomPos(1,RDGeom::Point3D(1.54,0,0));

  m2 = MolOps::addHs(*m,false,true);
  conf2 = &(m2->getConformer());
  TEST_ASSERT(m2->getNumAtoms()==8);
  TEST_ASSERT(feq((conf2->getAtomPos(0) - conf2->getAtomPos(2)).length(),
                  bondLength));
  TEST_ASSERT(feq((conf2->getAtomPos(0) - conf2->getAtomPos(3)).length(),
                  bondLength));
  TEST_ASSERT(feq((conf2->getAtomPos(0) - conf2->getAtomPos(4)).length(),
                  bondLength));
  TEST_ASSERT(feq((conf2->getAtomPos(2) - conf2->getAtomPos(3)).length(),
                  tetDist));
  TEST_ASSERT(feq((conf2->getAtomPos(2) - conf2->getAtomPos(4)).length(),
                  tetDist));
  TEST_ASSERT(feq((conf2->getAtomPos(3) - conf2->getAtomPos(4)).length(),
                  tetDist));
  TEST_ASSERT(feq((conf2->getAtomPos(1) - conf2->getAtomPos(5)).length(),
                  bondLength));
  TEST_ASSERT(feq((conf2->getAtomPos(1) - conf2->getAtomPos(6)).length(),
                  bondLength));
  TEST_ASSERT(feq((conf2->getAtomPos(1) - conf2->getAtomPos(7)).length(),
                  bondLength));
  TEST_ASSERT(feq((conf2->getAtomPos(5) - conf2->getAtomPos(6)).length(),
                  tetDist));
  TEST_ASSERT(feq((conf2->getAtomPos(5) - conf2->getAtomPos(7)).length(),
                  tetDist));
  TEST_ASSERT(feq((conf2->getAtomPos(6) - conf2->getAtomPos(7)).length(),
                  tetDist));

  delete m;
  delete m2;
  
  smi = "C=C";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==2);
  conf = new Conformer(2);
  m->addConformer(conf);
  conf->setAtomPos(0,RDGeom::Point3D(0,0,0));
  conf->setAtomPos(1,RDGeom::Point3D(1.3,0,0));

  m2 = MolOps::addHs(*m,false,true);
  
  conf2 = &(m2->getConformer());
  
  TEST_ASSERT(m2->getNumAtoms()==6);
  
  TEST_ASSERT(feq((conf2->getAtomPos(0) - conf2->getAtomPos(2)).length(),
                  bondLength));
  TEST_ASSERT(feq((conf2->getAtomPos(0) - conf2->getAtomPos(3)).length(),
                  bondLength));
  TEST_ASSERT(feq((conf2->getAtomPos(1) - conf2->getAtomPos(4)).length(),
                  bondLength));
  TEST_ASSERT(feq((conf2->getAtomPos(1) - conf2->getAtomPos(5)).length(),
                  bondLength));
  TEST_ASSERT(feq((conf2->getAtomPos(2) - conf2->getAtomPos(3)).length(),
                  sp2Dist));
  TEST_ASSERT(feq((conf2->getAtomPos(4) - conf2->getAtomPos(5)).length(),
                  sp2Dist));
  delete m;
  delete m2;
  
  smi = "C#C";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==2);
  conf = new Conformer(2);
  m->addConformer(conf);
  conf->setAtomPos(0,RDGeom::Point3D(0,0,0));
  conf->setAtomPos(1,RDGeom::Point3D(1.2,0,0));

  m2 = MolOps::addHs(*m,false,true);
  conf2 = &(m2->getConformer());
  TEST_ASSERT(m2->getNumAtoms()==4);
  TEST_ASSERT(feq((conf2->getAtomPos(0) - conf2->getAtomPos(2)).length(),
                  bondLength));
  TEST_ASSERT(feq((conf2->getAtomPos(1) - conf2->getAtomPos(3)).length(),
                  bondLength));
  TEST_ASSERT(feq((conf2->getAtomPos(0) - conf2->getAtomPos(3)).length(),
                  bondLength+1.2));
  TEST_ASSERT(feq((conf2->getAtomPos(1) - conf2->getAtomPos(2)).length(),
                  bondLength+1.2));

  delete m;
  delete m2;

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testSanitOps()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Sanitization special cases" << std::endl;
  ROMol *m;
  std::string smi,pathName;

  smi = "CN(=O)=O";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  TEST_ASSERT(m->getAtomWithIdx(1)->getFormalCharge()==1);
  delete m;

  smi = "C[N+](=O)[O-]";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  TEST_ASSERT(m->getAtomWithIdx(1)->getFormalCharge()==1);
  delete m;

  smi = "Cl(=O)(=O)(=O)[O-]";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==5);
  TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge()==3);
  delete m;

  pathName=getenv("RDBASE");
  pathName += "/Code/GraphMol/test_data/";
  m = MolFileToMol(pathName+"perchlorate1.mol");
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==51);
  TEST_ASSERT(m->getAtomWithIdx(7)->getFormalCharge()==3);
  delete m;

  smi = "CN=N#N";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  TEST_ASSERT(m->getAtomWithIdx(1)->getFormalCharge()==0);
  TEST_ASSERT(m->getAtomWithIdx(2)->getFormalCharge()==1);
  TEST_ASSERT(m->getAtomWithIdx(3)->getFormalCharge()==-1);
  TEST_ASSERT(m->getBondBetweenAtoms(2,3)->getBondType()==Bond::DOUBLE);
  delete m;

  smi = "N#N=NC";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  TEST_ASSERT(m->getAtomWithIdx(2)->getFormalCharge()==0);
  TEST_ASSERT(m->getAtomWithIdx(1)->getFormalCharge()==1);
  TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge()==-1);
  TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getBondType()==Bond::DOUBLE);
  delete m;

  smi = "N#N";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==2);
  TEST_ASSERT(m->getAtomWithIdx(1)->getFormalCharge()==0);
  TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge()==0);
  TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getBondType()==Bond::TRIPLE);
  delete m;
  
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testAddConformers() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Add Confomers" << std::endl;

  std::string smi = "CC";
  ROMol *m = SmilesToMol(smi);
  unsigned int i;
  for (i = 0; i < 5; i++) {
    Conformer *conf = new Conformer(2);
    conf->setAtomPos(0, RDGeom::Point3D(0.0, 0.0, 0.0));
    conf->setAtomPos(1, RDGeom::Point3D(1.5, 0.0, 0.0));
    m->addConformer(conf, true);
  }
  CHECK_INVARIANT(m->getNumConformers() == 5, "");
  
  ROMol *m2 = MolOps::addHs(*m,false,true);
  CHECK_INVARIANT(m2->getNumConformers() == 5, "");
  //const ROMol::CONF_SPTR_LIST &confs = m2->getConformers();
  ROMol::ConstConformerIterator ci;
  i = 0;
  for (ci = m2->beginConformers(); ci != m2->endConformers(); ci++) {
    CHECK_INVARIANT((*ci)->getNumAtoms() == 8, "");
    CHECK_INVARIANT((*ci)->getId() == i, "");
    const ROMol *mn = &((*ci)->getOwningMol());
    CHECK_INVARIANT(mn->getNumAtoms() == 8, "");
    i++;
  }
  //std::cout << m2->getNumAtoms() << " " << m2->getNumConformers() << "\n";
  delete m;
  delete m2;
  BOOST_LOG(rdInfoLog) << "Finished \n ";
  
}

void testIssue252() {
  // lets check if we can sanitize C60
  std::string smi = "C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3=C4C4=C5C5=C9C6=C7C6=C7C8=C1C1=C8C%10=C%10C%11=C2C2=C3C3=C4C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1=C%12C5=C%11C4=C3C3=C5C(=C81)C%10=C23";
  ROMol *mol = SmilesToMol(smi);
  for(ROMol::BondIterator it=mol->beginBonds();
      it!=mol->endBonds();it++){
    TEST_ASSERT((*it)->getIsAromatic());
  }
  std::string asmi = MolToSmiles(*mol);
  // check if we can do it in the aromatic form
  ROMol *nmol = SmilesToMol(asmi);
  for(ROMol::BondIterator it=nmol->beginBonds();
      it!=nmol->endBonds();it++){
    TEST_ASSERT((*it)->getIsAromatic());
  }

  std::string nsmi = MolToSmiles(*nmol);
  delete mol;
  delete nmol;
  // This is a check for Issue253
  CHECK_INVARIANT(asmi == nsmi, "");

}

void testIssue276() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Issue 276" << std::endl;
  std::string smi = "CP1(C)=CC=CN=C1C";
  ROMol *mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  // as of this writing, I'm not 100% sure what the right answer is here,
  // but the hybridization definitely should *not* be SP2:
  TEST_ASSERT(mol->getAtomWithIdx(1)->getHybridization()>Atom::SP2);
  delete mol;

  BOOST_LOG(rdInfoLog) << "Finished \n ";
}

void testHsAndAromaticity() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Additional Aromaticity Cases" << std::endl;
  std::string smi;
  ROMol *mol;

  smi = "[CH]1-[CH]-[CH]-[CH]-[CH]-[CH]-1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  //std::cerr << mol->getAtomWithIdx(0)->getHybridization() << std::endl;
  TEST_ASSERT(mol->getAtomWithIdx(0)->getHybridization()==Atom::SP3);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getImplicitValence()==0);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getNumImplicitHs()==0);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getNumRadicalElectrons()==1);
  TEST_ASSERT(!mol->getAtomWithIdx(0)->getIsAromatic());
  TEST_ASSERT(!mol->getBondBetweenAtoms(0,1)->getIsAromatic());

  smi = "C1=CC(=C)C(=C)C=C1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getHybridization()==Atom::SP2);
  TEST_ASSERT(mol->getAtomWithIdx(2)->getHybridization()==Atom::SP2);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getIsAromatic());
  TEST_ASSERT(mol->getBondBetweenAtoms(0,1)->getIsAromatic());

  delete mol;

  BOOST_LOG(rdInfoLog) << "Finished \n ";
}


void testSFIssue1694023()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 1694023 (bad chiral smiles after removing Hs) " << std::endl;
  ROMol *m;

  std::string smi;

  smi = "[C@@]([H])(F)(Cl)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);

  smi = "[C@@](F)([H])(Cl)Br";
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);

  smi = "[C@@](F)(Cl)([H])Br";
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);

  smi = "[C@@](F)(Cl)(Br)[H]";
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  
  smi = "[H][C@@](F)(Cl)Br";
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  
  smi = "F[C@@]([H])(Cl)Br";
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  
  smi = "F[C@@](Cl)([H])Br";
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  
  smi = "F[C@@](Cl)(Br)[H]";
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==4);
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW);
  
  smi = "C1CO[C@@H]1Cl";
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==5);
  TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);

  smi = "C1CO[C@]1([H])Cl";
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==5);
  TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  
  smi = "C1CO[C@@]1(Cl)[H]";
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms()==5);
  TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW);
  
  
  delete m;
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testSFIssue1719053()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 1719053 (Ring stereochemistry incorrectly removed) " << std::endl;
  ROMol *m;

  std::string smi;

  smi = "C[C@@H]1CCCCC1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_UNSPECIFIED);

  delete m;
  smi = "C[C@@H]1CC[C@@H](C)CC1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(m->getAtomWithIdx(4)->getChiralTag()!=Atom::CHI_UNSPECIFIED);

  delete m;
  smi = "C[C@@H]1C(C)CCCC1C";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_UNSPECIFIED);

  delete m;
  smi = "C[C@@H]1[C@H](C)CCC[C@H]1C";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  // this is a truly symmetric case, so the stereochem should be removed:
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(m->getAtomWithIdx(2)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(m->getAtomWithIdx(7)->getChiralTag()!=Atom::CHI_UNSPECIFIED);

  delete m;
  smi = "C[C@@H]1C=C[C@@H](C)C=C1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(m->getAtomWithIdx(4)->getChiralTag()!=Atom::CHI_UNSPECIFIED);

  delete m;
  smi = "C[N@@]1C=C[C@@H](C)C=C1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(m->getAtomWithIdx(4)->getChiralTag()==Atom::CHI_UNSPECIFIED);

  delete m;
  smi = "C[N@@]1CC[C@@H](C)CC1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(m->getAtomWithIdx(4)->getChiralTag()!=Atom::CHI_UNSPECIFIED);

  
  delete m;
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testSFIssue1811276()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 1811276 (kekulization failing) " << std::endl;
  ROMol *m;

  std::string smi;

  smi = "[O-]N1C=C[N+](=O)C=C1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  smi = MolToSmiles(*m);
  TEST_ASSERT(smi=="O=[n+]1ccn([O-])cc1");
  delete m;

  smi="o1ccc(=O)cc1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  smi = MolToSmiles(*m);
  TEST_ASSERT(smi=="O=c1ccocc1");
    
  smi="O=[n+]1ccocc1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  smi = MolToSmiles(*m);
  TEST_ASSERT(smi=="O=[n+]1ccocc1");

  smi="O=[n+]1ccn([O-])cc1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  smi = MolToSmiles(*m);
  TEST_ASSERT(smi=="O=[n+]1ccn([O-])cc1");

  smi="O=n1ccccc1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  smi = MolToSmiles(*m);
  TEST_ASSERT(smi=="[O-][n+]1ccccc1");

  delete m;
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testSFIssue1836576()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 1836576 (sanitization crash) " << std::endl;
  RWMol *m;

  std::string smi;
  bool ok;
  
  // the original form of the test runs foul of the rules for explicit
  // valence on B:
  smi = "[BH]123[BH]45[BH]167[BH]289[BH]312[BH]838[BH]966[Co]74479%10%11%12[CH]633[BH]811[CH]345[BH]21[BH]1234[BH]75[BH]911[BH]226[BH]%1011[BH]227[BH]633[BH]44[BH]322[CH]%1145[CH]%12271";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);

  unsigned int opThatFailed;
  ok=false;
  try {
    MolOps::sanitizeMol(*m,opThatFailed);
  } catch (MolSanitizeException &vee){
    ok=true;
  }
  TEST_ASSERT(ok);
  TEST_ASSERT(opThatFailed==MolOps::SANITIZE_PROPERTIES);

  // this molecule shows a known bug related to ring
  // ring finding in a molecule where all atoms are 4 connected.
  smi = "C123C45C11C44C55C22C33C14C523";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  
  ok=false;
  try {
    MolOps::sanitizeMol(*m,opThatFailed);
  } catch (ValueErrorException &vee){
    ok=true;
  }
  TEST_ASSERT(ok);
  TEST_ASSERT(opThatFailed==MolOps::SANITIZE_SYMMRINGS);

  delete m;

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}


void testChiralityAndRemoveHs()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing impact of removeHs on chirality" << std::endl;
  ROMol *m,*m2;

  std::string smi,code;

  smi = "F[C@]([H])(Cl)Br";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  m2=MolOps::removeHs(*m);
  TEST_ASSERT(m2);
  MolOps::assignStereochemistry(*m2,true,true);
  TEST_ASSERT(m2->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m2->getAtomWithIdx(1)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  delete m;
  delete m2;

  smi = "F[C@H](Cl)Br";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  m2=MolOps::removeHs(*m);
  TEST_ASSERT(m2);
  MolOps::assignStereochemistry(*m2,true,true);
  TEST_ASSERT(m2->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m2->getAtomWithIdx(1)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  delete m;
  delete m2;

  smi = "[C@]([H])(Cl)(F)Br";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true,true);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  m2=MolOps::removeHs(*m);
  TEST_ASSERT(m2);
  MolOps::assignStereochemistry(*m2,true,true);
  TEST_ASSERT(m2->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m2->getAtomWithIdx(0)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  delete m;
  delete m2;

  smi = "[C@H](Cl)(F)Br";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true,true);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  m2=MolOps::removeHs(*m);
  TEST_ASSERT(m2);
  MolOps::assignStereochemistry(*m2,true,true);
  TEST_ASSERT(m2->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m2->getAtomWithIdx(0)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  delete m;
  delete m2;

  smi = "[H]1.F[C@]1(Cl)Br";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true,true);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  m2=MolOps::removeHs(*m);
  TEST_ASSERT(m2);
  MolOps::assignStereochemistry(*m2,true,true);
  TEST_ASSERT(m2->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m2->getAtomWithIdx(1)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  delete m;
  delete m2;

  smi = "F[C@]1(Cl)Br.[H]1";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  m2=MolOps::removeHs(*m);
  TEST_ASSERT(m2);
  MolOps::assignStereochemistry(*m2,true,true);
  TEST_ASSERT(m2->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m2->getAtomWithIdx(1)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  delete m;
  delete m2;

  smi = "[H]1.[C@]1(Cl)(F)Br";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  m2=MolOps::removeHs(*m);
  TEST_ASSERT(m2);
  MolOps::assignStereochemistry(*m2,true,true);
  TEST_ASSERT(m2->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m2->getAtomWithIdx(0)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  delete m;
  delete m2;

  smi = "[C@]1(Cl)(F)Br.[H]1";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true,true);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  m2=MolOps::removeHs(*m);
  TEST_ASSERT(m2);
  MolOps::assignStereochemistry(*m2,true,true);
  TEST_ASSERT(m2->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m2->getAtomWithIdx(0)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  delete m;
  delete m2;
  

  smi = "Cl1.F2.Br3.[C@H]123";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true,true);
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
  m->getAtomWithIdx(3)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  delete m;

  smi = "[C@H]123.Cl1.F2.Br3";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true,true);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPCode"));
  m->getAtomWithIdx(0)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  delete m;

  smi = "F2.Cl1.Br3.[C@H]123";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true,true);
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
  m->getAtomWithIdx(3)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  delete m;


  smi = "Cl2.F1.Br3.[C@H]213";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  MolOps::assignStereochemistry(*m,true,true);
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp("_CIPCode"));
  m->getAtomWithIdx(3)->getProp("_CIPCode",code);
  TEST_ASSERT(code=="R");
  delete m;

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}


void testSFIssue1894348()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing SFIssue1894348 (impact of removeHs on bond stereo atoms" << std::endl;
  RWMol *m,*m2;

  std::string smi;

  smi = "Cl/C([H])=C/Cl";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  MolOps::sanitizeMol(*m);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getBondWithIdx(2)->getStereoAtoms().size()==2);
  TEST_ASSERT(m->getBondWithIdx(2)->getStereoAtoms()[0]==0);
  TEST_ASSERT(m->getBondWithIdx(2)->getStereoAtoms()[1]==4);
  m2=static_cast<RWMol *>(MolOps::removeHs(static_cast<const ROMol &>(*m)));
  TEST_ASSERT(m->getBondWithIdx(2)->getStereoAtoms().size()==2);
  TEST_ASSERT(m->getBondWithIdx(2)->getStereoAtoms()[0]==0);
  TEST_ASSERT(m->getBondWithIdx(2)->getStereoAtoms()[1]==4);
  TEST_ASSERT(m2->getBondWithIdx(1)->getStereoAtoms().size()==2);
  TEST_ASSERT(m2->getBondWithIdx(1)->getStereoAtoms()[0]==0);
  TEST_ASSERT(m2->getBondWithIdx(1)->getStereoAtoms()[1]==3);

  delete m;
  delete m2;

  smi = "Cl/C([H])=C/Cl";
  m = SmilesToMol(smi,false,false);
  TEST_ASSERT(m);
  MolOps::sanitizeMol(*m);
  TEST_ASSERT(m->getBondWithIdx(2)->getStereoAtoms().size()==0);
  m2=static_cast<RWMol *>(MolOps::removeHs(static_cast<const ROMol &>(*m)));
  // if we don't assign stereocodes in the original we shouldn't have them here:
  TEST_ASSERT(m2->getBondWithIdx(1)->getStereoAtoms().size()==0);
  delete m;
  delete m2;

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testAromaticityEdges()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing some aromaticity edge cases " << std::endl;
  RWMol *m;

  std::string smi;

  // ------
  // this was sf.net bug 1934360
  smi = "C1=C=NC=N1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(!m->getAtomWithIdx(0)->getIsAromatic());
  TEST_ASSERT(!m->getBondWithIdx(0)->getIsAromatic());
  delete m;

  smi = "C1=CNC=N1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(0)->getIsAromatic());
  TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
  delete m;

  smi = "C=[C+]1=CNC=N1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(!m->getAtomWithIdx(1)->getIsAromatic());
  TEST_ASSERT(!m->getBondWithIdx(1)->getIsAromatic());
  delete m;

  // ------
  // this was sf.net bug 1940646
  smi = "C1#CC=C1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(!m->getAtomWithIdx(0)->getIsAromatic());
  TEST_ASSERT(!m->getBondWithIdx(0)->getIsAromatic());
  delete m;
  smi = "C1#CC=CC=C1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(0)->getIsAromatic());
  TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
  delete m;

  // ------
  // this was sf.net bug 2091839
  
  smi = "c1cccc[c]1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(0)->getIsAromatic());
  TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
  delete m;
  
  smi = "C1=CC=CC=[C]1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(0)->getIsAromatic());
  TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
  delete m;

  smi = "c1cccc[n+]1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(0)->getIsAromatic());
  TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
  delete m;

#if 0 // currently fails
  smi = "[n]1cccc1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(0)->getIsAromatic());
  TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==1);
  TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
  delete m;
#endif
  
  smi = "[n]1ccccc1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(0)->getIsAromatic());
  TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==0);
  TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
  delete m;

  smi = "[H]n1cccc1";
  m = SmilesToMol(smi,0,0);
  TEST_ASSERT(m);
  MolOps::sanitizeMol(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->getIsAromatic());
  TEST_ASSERT(m->getAtomWithIdx(1)->getNumRadicalElectrons()==0);
  TEST_ASSERT(!m->getAtomWithIdx(0)->getIsAromatic());
  TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==0);
  delete m;

  smi = "[H]";
  m = SmilesToMol(smi,0,0);
  TEST_ASSERT(m);
  MolOps::sanitizeMol(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==1);
  delete m;
  
  // ------
  // this was sf.net bug 2787221.
  smi = "O=C1C(=O)C=C1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(1)->getIsAromatic());
  TEST_ASSERT(m->getBondBetweenAtoms(1,2)->getIsAromatic());
  delete m;


  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testSFIssue1942657()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 1942657 " << std::endl;
  RWMol *m;

  std::string smi;

  smi = "C[C](C)(C)(C)C";
  try{
    m = SmilesToMol(smi);
  } catch (MolSanitizeException &e){
    m=0;
  }
  TEST_ASSERT(!m);

  smi = "C[CH](C)(C)C";
  try{
    m = SmilesToMol(smi);
  } catch (MolSanitizeException &e){
    m=0;
  }
  TEST_ASSERT(!m);

  smi = "C[C](=C)(C)C";
  try{
    m = SmilesToMol(smi);
  } catch (MolSanitizeException &e){
    m=0;
  }
  TEST_ASSERT(!m);

  smi = "C[Si](=C)(=C)=C";
  try{
    m = SmilesToMol(smi);
  } catch (MolSanitizeException &e){
    m=0;
  }
  TEST_ASSERT(!m);


  
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}


void testSFIssue1968608()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 198608 " << std::endl;
  RWMol *m;

  std::string smi;

  smi = "C1CC1CC1CC1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m->getRingInfo()->minAtomRingSize(0)==3);
  TEST_ASSERT(m->getRingInfo()->minAtomRingSize(3)==0);
  TEST_ASSERT(m->getRingInfo()->minBondRingSize(0)==3);
  TEST_ASSERT(m->getRingInfo()->minBondRingSize(3)==0);

  
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testHybridization()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing hybridization assignment " << std::endl;

  {
    RWMol *m;
    std::string smi="CCC";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getHybridization()==Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="CNC";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getHybridization()==Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="COC";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getHybridization()==Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="C[C-2]C";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getHybridization()==Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="C[CH-]C";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getHybridization()==Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="C[CH]C";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getHybridization()==Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="C[C]C";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getHybridization()==Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="C[C-]C";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getHybridization()==Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="C[CH+]C";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP2);
    TEST_ASSERT(m->getAtomWithIdx(0)->getHybridization()==Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="CC=C";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP2);
    TEST_ASSERT(m->getAtomWithIdx(0)->getHybridization()==Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="CN=C";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP2);
    TEST_ASSERT(m->getAtomWithIdx(0)->getHybridization()==Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="C[C-]=C";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP2);
    TEST_ASSERT(m->getAtomWithIdx(0)->getHybridization()==Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="C[C]=C";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP2);
    TEST_ASSERT(m->getAtomWithIdx(0)->getHybridization()==Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="C[N+]=C";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP2);
    TEST_ASSERT(m->getAtomWithIdx(0)->getHybridization()==Atom::SP3);
    delete m;
  }


  
  {
    RWMol *m;
    std::string smi="C#C";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="C#[C-]";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="C#[C]";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="C[O]";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi="C[N-]";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP3);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testSFNetIssue2196817() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 2196817: handling of aromatic dummies" << std::endl;

  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName+"dummyArom.mol");
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getAtomicNum()==0);
    TEST_ASSERT(m->getAtomWithIdx(0)->getIsAromatic()==true);

    MolOps::Kekulize(*m);
    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getBondType()==Bond::SINGLE);
    TEST_ASSERT(m->getBondBetweenAtoms(0,4)->getBondType()==Bond::SINGLE);
      
    delete m;
  }
  
  {
    std::string smi="*1cncc1";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getAtomicNum()==0);
    TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic()==true);
    delete m;
  }
  
  {
    std::string smi="*1C=NC=C1";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getAtomicNum()==0);
    TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic()==true);
    delete m;
  }

  {
    // case where all must be ignored:
    std::string smi="c1*ccc1-c1*ccc1-c1*ccc1";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    delete m;
  }
  
  {
    std::string smi="c1*[nH]*c1";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    smi=MolToSmiles(*m);
    TEST_ASSERT(smi=="[*]1cc[*][nH]1");
    delete m;
    smi="c1***c1";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    smi=MolToSmiles(*m);
    TEST_ASSERT(smi=="[*]1:[*]cc[*]:1");
    delete m;
    smi="c:1:*:*:*:*1";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    smi=MolToSmiles(*m);
    TEST_ASSERT(smi=="[*]1:[*]:[*]c[*]:1");

    delete m;
    smi="*:1:*:*:*:*:1";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    smi=MolToSmiles(*m);
    TEST_ASSERT(smi=="[*]1:[*]:[*]:[*]:[*]:1");
    delete m;
  }
  
  {
    std::string smi="c1*[nH]cc1-c1*[nH]cc1-c1*ccc1";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    delete m;
    smi="c1*[nH]cc1-c1*ccc1-c1*[nH]cc1";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    delete m;
    smi="c1*ccc1-c1*[nH]cc1-c1*[nH1]cc1";
    m = SmilesToMol(smi);
    TEST_ASSERT(m);
    delete m;
  }
  
  {
    std::string smi="c1*[nH]cc1-c1*[nH]cc1-c1*[nH]cc1";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    delete m;
  }
  
  {
    std::string smi="c1ccc(C2CC(n4cc[*]c4=C2))cc1";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getIsAromatic());
    TEST_ASSERT(m->getBondBetweenAtoms(0,14)->getIsAromatic());
    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getBondType()==Bond::AROMATIC);
    TEST_ASSERT(m->getBondBetweenAtoms(0,14)->getBondType()==Bond::AROMATIC);
    MolOps::Kekulize(*m);
    TEST_ASSERT(!m->getBondBetweenAtoms(0,1)->getIsAromatic());
    TEST_ASSERT(!m->getBondBetweenAtoms(0,14)->getIsAromatic());
    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getBondType()==Bond::DOUBLE||
                m->getBondBetweenAtoms(0,14)->getBondType()==Bond::DOUBLE);
    MolOps::setAromaticity(*m);
    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getIsAromatic());
    TEST_ASSERT(m->getBondBetweenAtoms(0,14)->getIsAromatic());
    TEST_ASSERT(m->getBondBetweenAtoms(0,1)->getBondType()==Bond::AROMATIC);
    TEST_ASSERT(m->getBondBetweenAtoms(0,14)->getBondType()==Bond::AROMATIC);
    delete m;
  }
  
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}


void testSFNetIssue2208994() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 2208994 : kekulization error" << std::endl;

  {
    std::string smi="Cn1ccc(=O)n1C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getIsAromatic()==true);
    TEST_ASSERT(m->getBondWithIdx(1)->getIsAromatic()==true);

    delete m;
  }
  
  {
    std::string smi="c:1:c:c:c:c:c1";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getIsAromatic()==true);
    TEST_ASSERT(m->getBondWithIdx(1)->getIsAromatic()==true);

    delete m;
  }
  
  {
    std::string smi="c1:c:c:c:c:c:1";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getIsAromatic()==true);
    TEST_ASSERT(m->getBondWithIdx(1)->getIsAromatic()==true);

    delete m;
  }
  
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testSFNetIssue2313979() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 2313979: aromaticity assignment hangs " << std::endl;
  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    SDMolSupplier suppl(pathName+"Issue2313979.sdf",false);

    while(!suppl.atEnd()){
      ROMol *m=suppl.next();
      TEST_ASSERT(m);
      std::string nm;
      m->getProp("_Name",nm);
      BOOST_LOG(rdInfoLog) << "   Doing molecule: "<<nm<< std::endl;
      
      BOOST_LOG(rdInfoLog) << "     This should finish in a few seconds.  >>>" << std::endl;
      MolOps::sanitizeMol(*(RWMol *)m);
      delete m;
      BOOST_LOG(rdInfoLog) << "   <<< Done." << std::endl;
    }

  }
  
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}



void testSFNetIssue2316677() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 2316677 : canonicalization error" << std::endl;
  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName+"Issue2316677.mol");
    TEST_ASSERT(m);
    std::string smi=MolToSmiles(*m,true);
    std::cerr<<"smi: "<<smi<<std::endl;
    TEST_ASSERT(smi=="Cc1ccc(S(=O)(=O)/N=C2\\CC(=N\\C(C)(C)C)/C2=N\\C(C)(C)C)cc1");
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testSanitizeNonringAromatics() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 2830244: make sure that non-ring aromatic atoms generate errors:" << std::endl;
  {
    std::string smi="c-C";
    
    RWMol *m = SmilesToMol(smi,0,false);
    bool ok=false;
    try {
      MolOps::Kekulize(*m);
    } catch (MolSanitizeException &vee){
      ok=true;
    }
    TEST_ASSERT(ok);
    delete m;
  }
  {
    std::string smi="c-C";
    
    RWMol *m = SmilesToMol(smi,0,false);
    bool ok=false;
    unsigned int opThatFailed;
    try {
      MolOps::sanitizeMol(*m,opThatFailed);
    } catch (MolSanitizeException &vee){
      ok=true;
    }
    TEST_ASSERT(ok);
    TEST_ASSERT(opThatFailed==MolOps::SANITIZE_KEKULIZE);
    delete m;
  }
  
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}


void testSFNetIssue2951221() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 2951221 : hydrogens added with bad coordinates" << std::endl;
  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    ROMol *m = MolFileToMol(pathName+"Issue2951221.1.mol");
    TEST_ASSERT(m);
    ROMol *m2 = MolOps::addHs(*m,false,true);
    TEST_ASSERT(m2);
    delete m;
    TEST_ASSERT(m2->getNumAtoms(false)==12);
    RDGeom::Point3D coords[4];
    coords[0]= m2->getConformer().getAtomPos(2);
    coords[1]= m2->getConformer().getAtomPos(0);
    coords[2]= m2->getConformer().getAtomPos(1);
    coords[3]= m2->getConformer().getAtomPos(9);
    double dot=(coords[3]-coords[0]).dotProduct((coords[1]-coords[0]).crossProduct(coords[2]-coords[0]));
    TEST_ASSERT(dot>1.0);
  }

  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    ROMol *m = MolFileToMol(pathName+"Issue2951221.2.mol");
    TEST_ASSERT(m);
    ROMol *m2 = MolOps::addHs(*m,false,true);
    TEST_ASSERT(m2);
    delete m;
    TEST_ASSERT(m2->getNumAtoms(false)==5);
    MolOps::assignChiralTypesFrom3D(*m2);
    MolOps::assignStereochemistry(*m2,true,true);
    TEST_ASSERT(m2->getAtomWithIdx(1)->hasProp("_CIPCode"));
    std::string cip;
    m2->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");
  }
  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    ROMol *m = MolFileToMol(pathName+"Issue2951221.3.mol");
    TEST_ASSERT(m);
    ROMol *m2 = MolOps::addHs(*m,false,true);
    TEST_ASSERT(m2);
    delete m;
    TEST_ASSERT(m2->getNumAtoms(false)==5);
    MolOps::assignChiralTypesFrom3D(*m2);
    MolOps::assignStereochemistry(*m2,true,true);
    TEST_ASSERT(m2->getAtomWithIdx(1)->hasProp("_CIPCode"));
    std::string cip;
    m2->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testSFNetIssue2952255() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 2952255 : bad assignment of radicals to early elements" << std::endl;
  {
    std::string smi="[C](C)(C)C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==1);
    delete m;
  }
  {
    std::string smi="[C](C)C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==2);
    delete m;
  }
  {
    std::string smi="[CH](C)C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==1);
    delete m;
  }
  {
    std::string smi="[CH+](C)C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==0);
    delete m;
  }
  {
    std::string smi="[C-](C)C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==1);
    delete m;
  }
  {
    std::string smi="[C+](C)(C)(C)C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==1);
    delete m;
  }
  {
    std::string smi="C(C)(C)(C)C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==0);
    delete m;
  }
  {
    std::string smi="[N](C)C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==1);
    delete m;
  }
  {
    std::string smi="[N+](C)(C)C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==1);
    delete m;
  }
  {
    std::string smi="[Cl]";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==1);
    delete m;
  }
  {
    std::string smi="[Cl-]";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==0);
    delete m;
  }
  {
    std::string smi="[Cl]C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==0);
    delete m;
  }
  {
    std::string smi="[Na]";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==1);
    delete m;
  }
  {
    std::string smi="[Na+]";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==0);
    delete m;
  }
  {
    std::string smi="[Na]C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==0);
    delete m;
  }
  {
    std::string smi="[Mg+]C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==0);
    delete m;
  }
  {
    std::string smi="[Mg]C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==1);
    delete m;
  }
  {
    std::string smi="[Mg+]";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==1);
    delete m;
  }
  {
    std::string smi="[Mg+2]";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons()==0);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testSFNetIssue3185548() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 3185548 : problems with SSSR code" << std::endl;

  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    BOOST_LOG(rdInfoLog) << "  Starting file read 1" << std::endl;
    RWMol *m = MolFileToMol(pathName+"Issue3185548.mol");
    BOOST_LOG(rdInfoLog) << "  finished" << std::endl;
    TEST_ASSERT(m);
  }

  
  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    BOOST_LOG(rdInfoLog) << "  Starting file read 2" << std::endl;
    RWMol *m = MolFileToMol(pathName+"Issue3185548.2.mol");
    BOOST_LOG(rdInfoLog) << "  finished" << std::endl;
    TEST_ASSERT(m);

    m->getRingInfo()->reset();
    unsigned int nsssr;
    VECT_INT_VECT sssrs;
    nsssr=MolOps::findSSSR(*m,sssrs);
    TEST_ASSERT(nsssr=48);
    nsssr=MolOps::symmetrizeSSSR(*m,sssrs);
    TEST_ASSERT(nsssr=56);    
  }

  
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testSFNetIssue3349243(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 3349243" << std::endl;
  {
    std::string smi="c1cccc[n+]1";
    RWMol *m=SmilesToMol(smi);
    TEST_ASSERT(m);
    MolOps::Kekulize(*m);
    // just finishing is good
    TEST_ASSERT(m->getBondWithIdx(0)->getBondType()!=Bond::AROMATIC);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testFastFindRings(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing fast find rings" << std::endl;
  {
    std::string smi="CCC";
    RWMol *m=SmilesToMol(smi,0,0);
    TEST_ASSERT(m);
    MolOps::fastFindRings(*m);
    TEST_ASSERT(m->getRingInfo());
    TEST_ASSERT(m->getRingInfo()->isInitialized());
    TEST_ASSERT(m->getRingInfo()->numRings()==0);
    delete m;
  }
  {
    std::string smi="C1CC1";
    RWMol *m=SmilesToMol(smi,0,0);
    TEST_ASSERT(m);
    MolOps::fastFindRings(*m);
    TEST_ASSERT(m->getRingInfo());
    TEST_ASSERT(m->getRingInfo()->isInitialized());
    TEST_ASSERT(m->getRingInfo()->numRings()==1);
    delete m;
  }

  {
    std::string smi="CC1CC1";
    RWMol *m=SmilesToMol(smi,0,0);
    TEST_ASSERT(m);
    MolOps::fastFindRings(*m);
    TEST_ASSERT(m->getRingInfo());
    TEST_ASSERT(m->getRingInfo()->isInitialized());
    TEST_ASSERT(m->getRingInfo()->numRings()==1);
    delete m;
  }

  {
    std::string smi="C1CC1.C1CC1";
    RWMol *m=SmilesToMol(smi,0,0);
    TEST_ASSERT(m);
    MolOps::fastFindRings(*m);
    TEST_ASSERT(m->getRingInfo());
    TEST_ASSERT(m->getRingInfo()->isInitialized());
    TEST_ASSERT(m->getRingInfo()->numRings()==2);
    delete m;
  }
  {
    std::string smi="C1C(C)C1";
    RWMol *m=SmilesToMol(smi,0,0);
    TEST_ASSERT(m);
    MolOps::fastFindRings(*m);
    TEST_ASSERT(m->getRingInfo());
    TEST_ASSERT(m->getRingInfo()->isInitialized());
    TEST_ASSERT(m->getRingInfo()->numRings()==1);
    delete m;
  }
  {
    std::string smi="c1c(=O)nc2[nH]cnn2c1O";
    RWMol *m=SmilesToMol(smi,0,0);
    TEST_ASSERT(m);
    MolOps::fastFindRings(*m);
    TEST_ASSERT(m->getRingInfo());
    TEST_ASSERT(m->getRingInfo()->isInitialized());
    TEST_ASSERT(m->getRingInfo()->numRings()==2);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}


void testSFNetIssue3487473(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 3487473" << std::endl;
  {
    std::string smi="C*C";
    RWMol *m=SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::UNSPECIFIED);
    delete m;
  }

  {
    std::string smi="C*C";
    RWMol *m = SmartsToMol(smi);
    TEST_ASSERT(m);
    m->updatePropertyCache(false);
    MolOps::setConjugation(*m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::UNSPECIFIED);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testSFNetIssue3480481(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 3480481" << std::endl;
  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName+"Issue3480481.mol");
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getIsAromatic()==true);
    TEST_ASSERT(m->getAtomWithIdx(0)->getExplicitValence()==4);
    TEST_ASSERT(m->getAtomWithIdx(0)->getImplicitValence()==0);
    TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge()==-1);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}


void aamatchtest(std::string smi1,std::string smi2,bool shouldMatch,int idx1,int idx2){
  RWMol *m1=SmilesToMol(smi1);
  RWMol *m2=SmilesToMol(smi2);
  TEST_ASSERT(m1);
  TEST_ASSERT(m2);
  //std::cerr<<"   "<<smi1<<" "<<smi2<<std::endl;
  TEST_ASSERT(m2->getAtomWithIdx(idx2)->Match(m1->getAtomWithIdx(idx1))==shouldMatch);
  delete m1;
  delete m2;
}

void testAtomAtomMatch(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Atom-Atom matching behavior" << std::endl;
  /* Here's what we're testing:

     | Molecule | Query   | Match |
     | CCO      | CCO     | Yes   |
     | CC[O-]   | CCO     | Yes   |
     | CCO      | CC[O-]  | No    |
     | CC[O-]   | CC[O-]  | Yes   |
     | CC[O-]   | CC[OH]  | Yes   |
     | CCOC     | CC[OH]  | Yes   |
     | CCOC     | CCO     | Yes   |
     | CCC      | CCC     | Yes   |
     | CC[14C]  | CCC     | Yes   |
     | CCC      | CC[14C] | No    |
     | CC[14C]  | CC[14C] | Yes   |
     | OCO      | C       | Yes   |
     | OCO      | [CH2]   | Yes   |
     | OCO      | [CH3]   | Yes   |
     | O[CH2]O  | C       | Yes   |
     | O[CH2]O  | [CH2]   | Yes   |
     | OCO      | [CH]    | Yes   |

     This is a large superset of issue 3495370

  */

  // note that in some cases here we have to be fairly careful about Hs on 
  // the query to make sure that it doesn't have radicals (radical handling
  // added to fix github #165
  aamatchtest("CCO","O",true,2,0);
  aamatchtest("CC[O-]","O",true,2,0);
  aamatchtest("CCO","[OH-]",false,2,0);
  aamatchtest("CC[O-]","[OH-]",true,2,0);
  aamatchtest("CC[O-]","[OH2]",true,2,0);
  aamatchtest("CCOC","[OH2]",true,2,0);
  aamatchtest("CCOC","O",true,2,0);
  aamatchtest("CCC","C",true,2,0);
  aamatchtest("CC[14C]","C",true,2,0);
  aamatchtest("CCC","[14CH4]",false,2,0);
  aamatchtest("CC[14C]","[14CH4]",true,2,0);
  aamatchtest("CC[13C]","[14CH4]",false,2,0);
  aamatchtest("OCO","C",true,1,0);
  aamatchtest("OCO","[CH4]",true,1,0);
  aamatchtest("O[CH2]O","C",true,1,0);
  aamatchtest("O[CH2]O","[CH4]",true,1,0);
  aamatchtest("OCO","[CH2]",false,1,0); // doesn't match due to radical count
  aamatchtest("O[CH2]O","[CH2]",false,1,0); // doesn't match due to radical count
  aamatchtest("O[CH]O","[CH3]",true,1,0);
  aamatchtest("O[CH]O","[CH2]",false,1,0); // doesn't match due to radical count
  aamatchtest("CC","*",false,1,0);
  aamatchtest("C*","*",true,1,0);
  aamatchtest("C[1*]","*",true,1,0);
  aamatchtest("C[1*]","[1*]",true,1,0);
  aamatchtest("C*","[1*]",true,1,0);
  aamatchtest("C[2*]","[1*]",false,1,0);

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}


void testSFNetIssue3525076(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 3525076" << std::endl;
  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName+"Issue3525076.sdf");
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(18)->getIsAromatic()==false);
    TEST_ASSERT(m->getBondWithIdx(18)->getBondType()==Bond::SINGLE);
    MolOps::Kekulize(*m);
    TEST_ASSERT(m->getBondWithIdx(18)->getIsAromatic()==false);
    TEST_ASSERT(m->getBondWithIdx(18)->getBondType()==Bond::SINGLE);
    MolOps::sanitizeMol(*m);
    TEST_ASSERT(m->getBondWithIdx(18)->getIsAromatic()==false);
    TEST_ASSERT(m->getBondWithIdx(18)->getBondType()==Bond::SINGLE);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testBasicCanon(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing canonicalization basics" << std::endl;
  // these are all cases that were problematic at one time or another during
  // the canonicalization rewrite.
#if 1
  {
    std::string smi="FC1C(=C/Cl)\\C1";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(2,3)->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(2,3)->getStereo()==Bond::STEREOZ);

    //std::cerr<<"-------------\n";
    std::string csmi1=MolToSmiles(*m,true);
    //std::cerr<<csmi1<<std::endl;
    RWMol *m2 = SmilesToMol(csmi1);
    TEST_ASSERT(m2);

    MatchVectType mv;
    TEST_ASSERT(SubstructMatch(*m,*m2,mv));
    std::map<int,int> mmap;
    for(MatchVectType::const_iterator mvit=mv.begin();mvit!=mv.end();++mvit){
      mmap[mvit->second]=mvit->first;
    }
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[2],mmap[3])->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[2],mmap[3])->getStereo()==Bond::STEREOZ);

    
    //std::cerr<<"-------------\n";
    std::string csmi2=MolToSmiles(*m2,true);
    //std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    TEST_ASSERT(csmi1==csmi2);
    delete m;
    delete m2;
  }



  {
    std::string smi="CC1(C)C2CCC1(C)C(=O)/C2=C\\C(N=N/c1ccccc1)=N/Nc1ccccc1";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(10,11)->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(10,11)->getStereo()==Bond::STEREOZ);
    TEST_ASSERT(m->getBondBetweenAtoms(12,21)->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(12,21)->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(13,14)->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(13,14)->getStereo()==Bond::STEREONONE);
    //std::cerr<<"-------------\n";
    std::string csmi1=MolToSmiles(*m,true);

    RWMol *m2 = SmilesToMol(csmi1);
    TEST_ASSERT(m2);

    MatchVectType mv;
    TEST_ASSERT(SubstructMatch(*m,*m2,mv));
    std::map<int,int> mmap;
    for(MatchVectType::const_iterator mvit=mv.begin();mvit!=mv.end();++mvit){
      mmap[mvit->second]=mvit->first;
    }
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[10],mmap[11])->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[10],mmap[11])->getStereo()==Bond::STEREOZ);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[12],mmap[21])->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[12],mmap[21])->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[13],mmap[14])->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[13],mmap[14])->getStereo()==Bond::STEREONONE);
    
    //std::cerr<<"-------------\n";
    std::string csmi2=MolToSmiles(*m2,true);
    //std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    TEST_ASSERT(csmi1==csmi2);
    delete m;
    delete m2;
  }


  {
    std::string smi="COc1ccc(OC)c2[nH]c(=O)cc(C)c21";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    //std::cerr<<"-------------\n";
    std::string csmi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(csmi1);
    TEST_ASSERT(m);
    //std::cerr<<"-------------\n";
    std::string csmi2=MolToSmiles(*m,true);
    //std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    TEST_ASSERT(csmi1==csmi2);
    delete m;
  }
  {
    std::string smi="COc1cc(C)c(C(=O)[O-])cc1OC";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    //std::cerr<<"-------------\n";
    std::string csmi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(csmi1);
    TEST_ASSERT(m);
    //std::cerr<<"-------------\n";
    std::string csmi2=MolToSmiles(*m,true);
    //std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    TEST_ASSERT(csmi1==csmi2);
    delete m;
  }
  {
    std::string smi="COc1ccc(C(=O)OC(c2ccc(OC)cc2)C(C)O)cc1";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    //std::cerr<<"-------------\n";
    std::string csmi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(csmi1);
    TEST_ASSERT(m);
    //std::cerr<<"-------------\n";
    std::string csmi2=MolToSmiles(*m,true);
    //std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    TEST_ASSERT(csmi1==csmi2);
    delete m;
  }
  {
    std::string smi="CC(C)C1CCC(C)=CC1=NNC(N)=O";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    //std::cerr<<"-------------\n";
    std::string csmi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(csmi1);
    TEST_ASSERT(m);
    //std::cerr<<"-------------\n";
    std::string csmi2=MolToSmiles(*m,true);
    //std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    TEST_ASSERT(csmi1==csmi2);
    delete m;
  }
  {
    std::string smi="COCCNC(=O)c1ccccc1N1C(=O)C2(C)c3[nH]c4ccccc4c3CCN2C1=O";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    //std::cerr<<"-------------\n";
    std::string csmi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(csmi1);
    TEST_ASSERT(m);
    //std::cerr<<"-------------\n";
    std::string csmi2=MolToSmiles(*m,true);
    //std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    TEST_ASSERT(csmi1==csmi2);
    delete m;
  }
  {
    std::string smi="Cc1c(Br)cc(Br)cc1C(F)(F)F";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    //std::cerr<<"-------------\n";
    std::string csmi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(csmi1);
    TEST_ASSERT(m);
    //std::cerr<<"-------------\n";
    std::string csmi2=MolToSmiles(*m,true);
    //std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    TEST_ASSERT(csmi1==csmi2);
    delete m;
  }
  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName+"zinc4235774a.mol");
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(1,2)->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(1,2)->getStereo()==Bond::STEREOZ);
    TEST_ASSERT(m->getBondBetweenAtoms(7,8)->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(7,8)->getStereo()==Bond::STEREOZ);
    std::string smi=MolToSmiles(*m,true);
    //std::cerr<<"SMILES: "<<smi<<std::endl;
    RWMol *m2 = SmilesToMol(smi);
    MatchVectType mv;
    TEST_ASSERT(SubstructMatch(*m,*m2,mv));
    std::map<int,int> mmap;
    for(MatchVectType::const_iterator mvit=mv.begin();mvit!=mv.end();++mvit){
      mmap[mvit->second]=mvit->first;
    }
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[1],mmap[2])->getBondType()==Bond::DOUBLE);    
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[7],mmap[8])->getBondType()==Bond::DOUBLE);        
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[1],mmap[2])->getStereo()==Bond::STEREOZ);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[7],mmap[8])->getStereo()==Bond::STEREOZ);
  }
  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName+"zinc4235774.mol");
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(4,5)->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(4,5)->getStereo()==Bond::STEREOZ);
    TEST_ASSERT(m->getBondBetweenAtoms(14,15)->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(14,15)->getStereo()==Bond::STEREOZ);
    std::string smi=MolToSmiles(*m,true);
    RWMol *m2 = SmilesToMol(smi);
    MatchVectType mv;
    TEST_ASSERT(SubstructMatch(*m,*m2,mv));
    std::map<int,int> mmap;
    for(MatchVectType::const_iterator mvit=mv.begin();mvit!=mv.end();++mvit){
      mmap[mvit->second]=mvit->first;
    }
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[4],mmap[5])->getBondType()==Bond::DOUBLE);    
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[14],mmap[15])->getBondType()==Bond::DOUBLE);        
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[4],mmap[5])->getStereo()==Bond::STEREOZ);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[14],mmap[15])->getStereo()==Bond::STEREOZ);
  }
#endif

  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName+"zinc3850436piece.mol");
    TEST_ASSERT(m);
    std::string csmi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(csmi1);
    //std::cerr<<"-------------\n";
    std::string csmi2=MolToSmiles(*m,true);
    //std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    TEST_ASSERT(csmi1==csmi2);
    delete m;
  }
  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName+"zinc13403961piece.mol");
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(1,2)->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(1,2)->getStereo()==Bond::STEREOZ);
    TEST_ASSERT(m->getBondBetweenAtoms(3,7)->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(3,7)->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(4,5)->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(4,5)->getStereo()==Bond::STEREOE);
    std::string csmi1=MolToSmiles(*m,true);
    //std::cerr<<"SMI1: "<<csmi1<<std::endl;
    RWMol *m2 = SmilesToMol(csmi1);
    MatchVectType mv;
    TEST_ASSERT(SubstructMatch(*m,*m2,mv));
    std::map<int,int> mmap;
    for(MatchVectType::const_iterator mvit=mv.begin();mvit!=mv.end();++mvit){
      mmap[mvit->second]=mvit->first;
    }

    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[1],mmap[2])->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[1],mmap[2])->getStereo()==Bond::STEREOZ);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[3],mmap[7])->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[3],mmap[7])->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[4],mmap[5])->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[4],mmap[5])->getStereo()==Bond::STEREOE);

    //std::cerr<<"-------------\n";
    std::string csmi2=MolToSmiles(*m2,true);
    //std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    TEST_ASSERT(csmi1==csmi2);
    delete m;
  }
  {
    std::string smi="C\\N=c1/s/c(=N\\Cl)/c/1=N/F";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    //std::cerr<<"-------------\n";
    std::string csmi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(csmi1);
    TEST_ASSERT(m);
    //std::cerr<<"-------------\n";
    std::string csmi2=MolToSmiles(*m,true);
    //std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    TEST_ASSERT(csmi1==csmi2);
    delete m;
  }
  {
    std::string smi="Cc1ccc(S(=O)(=O)/N=c2sc(=N\\C(C)(C)C)/c\\2=N/C(C)(C)C)cc1";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    //std::cerr<<"-------------\n";
    std::string csmi1=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(csmi1);
    TEST_ASSERT(m);
    //std::cerr<<"-------------\n";
    std::string csmi2=MolToSmiles(*m,true);
    //std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    TEST_ASSERT(csmi1==csmi2);
    delete m;
  }
  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName+"zinc6624278.mol");
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(21,13)->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(21,13)->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(5,12)->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(5,12)->getStereo()==Bond::STEREOZ);

    std::string csmi1=MolToSmiles(*m,true);
    //std::cerr<<"SMI1: "<<csmi1<<std::endl;
    RWMol *m2 = SmilesToMol(csmi1);
    MatchVectType mv;
    TEST_ASSERT(SubstructMatch(*m,*m2,mv));
    std::map<int,int> mmap;
    for(MatchVectType::const_iterator mvit=mv.begin();mvit!=mv.end();++mvit){
      mmap[mvit->second]=mvit->first;
    }

    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[21],mmap[13])->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[21],mmap[13])->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[5],mmap[12])->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[5],mmap[12])->getStereo()==Bond::STEREOZ);

    //std::cerr<<"-------------\n";
    std::string csmi2=MolToSmiles(*m2,true);
    //std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    TEST_ASSERT(csmi1==csmi2);
    delete m2;

    std::string tsmi=MolToSmiles(*m,true,false,7,false);
    //std::cerr<<"-------------\n";
    //std::cerr<<"T:\n"<<tsmi<<"\n-------------\n"<<std::endl;
    m2 = SmilesToMol(tsmi);
    TEST_ASSERT(SubstructMatch(*m,*m2,mv));
    mmap.clear();
    for(MatchVectType::const_iterator mvit=mv.begin();mvit!=mv.end();++mvit){
      mmap[mvit->second]=mvit->first;
    }
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[21],mmap[13])->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[21],mmap[13])->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[5],mmap[12])->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[5],mmap[12])->getStereo()==Bond::STEREOZ);
    
    //std::cerr<<"-------------\n";
    csmi2=MolToSmiles(*m2,true);
    //std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    TEST_ASSERT(csmi1==csmi2);
    delete m2;
    
    delete m;

    
  }
  {
    std::string smi="F/C=C/C=C(C)/C=C/Cl";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(1,2)->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(1,2)->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(3,4)->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(3,4)->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(6,7)->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(6,7)->getStereo()==Bond::STEREOE);



    std::string tsmi=MolToSmiles(*m,true,false,3,false);
    //std::cerr<<"-------------\n";
    //std::cerr<<"T:\n"<<tsmi<<"\n-------------\n"<<std::endl;
    RWMol *m2 = SmilesToMol(tsmi);
    MatchVectType mv;
    TEST_ASSERT(SubstructMatch(*m,*m2,mv));
    std::map<int,int> mmap;
    mmap.clear();
    for(MatchVectType::const_iterator mvit=mv.begin();mvit!=mv.end();++mvit){
      mmap[mvit->second]=mvit->first;
    }
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[1],mmap[2])->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[1],mmap[2])->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[3],mmap[4])->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[3],mmap[4])->getStereo()==Bond::STEREOE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[6],mmap[7])->getBondType()==Bond::DOUBLE);
    TEST_ASSERT(m2->getBondBetweenAtoms(mmap[6],mmap[7])->getStereo()==Bond::STEREOE);

    
    std::string csmi1=MolToSmiles(*m,true);
    //std::cerr<<"SMI1: "<<csmi1<<std::endl;
    //std::cerr<<"-------------\n";
    std::string csmi2=MolToSmiles(*m2,true);
    //std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    TEST_ASSERT(csmi1==csmi2);

    delete m2;
    delete m;
  }

  {
    // this was issue 3528556
    std::string smi="N12.N13.C24.C35.C46.C56";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    std::string csmi=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(csmi);
    TEST_ASSERT(m);
    smi=MolToSmiles(*m,true);
    TEST_ASSERT(csmi==smi);
  }
  {
    // this was issue 3526831
    std::string smi="CO/N=C/C(=C(\\O)/c1ccc(Cl)cc1)/C=N\\OC";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    std::string csmi=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(csmi);
    TEST_ASSERT(m);
    smi=MolToSmiles(*m,true);
    TEST_ASSERT(csmi==smi);
  }

  
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}


void testSFNetIssue3549146() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 3549146: problems after mergeQueryHs" << std::endl;

  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    ROMol *m = MolFileToMol(pathName+"Issue3549146.mol",true,false);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==16);
    ROMol *m2 = MolOps::mergeQueryHs(*m);
    TEST_ASSERT(m2);
    TEST_ASSERT(m2->getNumAtoms()==13);
    TEST_ASSERT(!(m2->getRingInfo()->isInitialized()));
    delete m;
    delete m2;
  }
  {
    std::string smi="CCC.C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==4);
    TEST_ASSERT((m->getRingInfo()->isInitialized()));
    m->addBond(1,3,Bond::SINGLE);
    TEST_ASSERT((m->getRingInfo()->isInitialized()));
    m->addBond(0,2,Bond::SINGLE);
    TEST_ASSERT(!(m->getRingInfo()->isInitialized()));

    delete m;
  }
  {
    std::string smi="C1CC1C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==4);
    TEST_ASSERT((m->getRingInfo()->isInitialized()));
    m->removeBond(2,3);
    TEST_ASSERT(!(m->getRingInfo()->isInitialized()));
    delete m;
  }
  {
    std::string smi="C1CC1C";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==4);
    TEST_ASSERT((m->getRingInfo()->isInitialized()));
    m->removeAtom(3);
    TEST_ASSERT(!(m->getRingInfo()->isInitialized()));
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}


void testSFNetIssue249() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 249: finding rings consumes all memory" << std::endl;

  {
    std::string smi="Cc1cc2cc(c1)C(=O)NCc1cc-3cc(CNC(=O)c4cc(C)cc(c4)C(=O)NCc4cc(cc(CNC2=O)c4O)-c2cc4CNC(=O)c5cc(C)cc(c5)C(=O)NCc5cc-3cc(CNC(=O)c3cc(C)cc(c3)C(=O)NCc(c2)c4O)c5O)c1O";
    ROMol *m = SmilesToMol(smi,0,0);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==88);
    m->updatePropertyCache(false);
    std::cerr<<"starting ring finding"<<std::endl;
    MolOps::findSSSR(*m);
    std::cerr<<"done"<<std::endl;    
    
    delete m;
  }

  {
    std::string smi="CCCOc1c2CNC(=O)c3cc(cc(c3)C(=O)NCc3cc4cc(CNC(=O)c5cc(C(=O)NCc1cc(c2)c1cc2CNC(=O)c6cc(cc(c6)C(=O)NCc6cc4cc(CNC(=O)c4cc(C(=O)NCc(c1)c2OCCC)cc(c4)C(=O)NC(COCCC(=O)O)(COCCC(=O)O)COCCC(=O)O)c6OCCC)C(=O)NC(COCCC(=O)O)(COCCC(=O)O)COCCC(=O)O)cc(c5)C(=O)NC(COCCC(=O)O)(COCCC(=O)O)COCCC(=O)O)c3OCCC)C(=O)NC(COCCC(=O)O)(COCCC(=O)O)COCCC(=O)O";
    ROMol *m = SmilesToMol(smi,0,0);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==196);
    m->updatePropertyCache(false);
    std::cerr<<"starting ring finding"<<std::endl;
    MolOps::findSSSR(*m);
    std::cerr<<"done"<<std::endl;    
    
    delete m;
  }

  {
    std::string smi="CCn1nnc(c1)CN=C(C1CC2C3CCC4C5C3C3C6C2C2C1C1CCC7C(C1)C1C8C9C7C(C(=O)O)C(C(=O)O)C7C9C(C9C8C8C%10C1C1C(C2C2C6C6C%11C3C(C5)C3C(C(=O)O)C5C%12CC9C9C8C8C(C%10)C%10C%13C(C%14C(C2C1C(=O)O)C6C1C2C%11C3C3C5C(C5C%12C9C(C8CC%10)CC5)C(CC3C2C(C(C1C%14CC%13C(=NCc1nnn(c1)CC)O)C(=O)O)C(=O)O)C(=NCc1nnn(c1)CC)O)C(=O)O)C(=O)O)CC1C(C4C(CC71)C(=NCc1nnn(c1)CC)O)C(=O)O)O";
    ROMol *m = SmilesToMol(smi,0,0);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==167);
    m->updatePropertyCache(false);
    std::cerr<<"starting ring finding"<<std::endl;
    MolOps::findSSSR(*m);
    std::cerr<<"done"<<std::endl;    
    delete m;
  }

  {
    std::string smi="C/C=N/c1ccc(cc1)c1cc2cc(c1)c1ccc(cc1)N=Cc1ccc(cc1)c1cc(cc(c1)c1ccc(cc1)C)c1ccc(cc1)C=Nc1ccc(cc1)c1cc(cc(c1)c1ccc(cc1)/N=C/C)c1ccc(cc1)N=Cc1ccc(cc1)c1cc(c3ccc(C=Nc4ccc(c5cc6c7ccc(N=Cc8ccc(c9cc(c%10ccc(C=Nc%11ccc2cc%11)cc%10)cc(c9)c2ccc(cc2)C=Nc2ccc(cc2)c2cc(cc(c2)c2ccc(cc2)N=Cc2ccc(cc2)c2cc(cc(c2)c2ccc(cc2)C)c2ccc(cc2)C=Nc2ccc(cc2)c2cc(c9ccc(N=Cc%10ccc(c%11cc(c%12ccc(C=Nc%13ccc(c(c6)c5)cc%13)cc%12)cc(c%11)c5ccc(cc5)C)cc%10)cc9)cc(c2)c2ccc(cc2)/N=C/C)c2ccc(cc2)/N=C/C)cc8)cc7)cc4)cc3)cc(c1)c1ccc(cc1)C";
    RWMol *m = SmilesToMol(smi,0,0);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==278);
    std::cerr<<"starting sanitization"<<std::endl;
    MolOps::sanitizeMol(*m);
    std::cerr<<"done"<<std::endl;    
    delete m;
  }

  {
    std::string smi="COc1cc2ccc1n1c(=O)n(c3ccc(cc3OC)c3ccc(c(c3)OC)n3c(=O)n(c4ccc(cc4OC)c4ccc(c(c4)OC)n4c(=O)n(c5ccc(cc5OC)c5ccc(n6c(=O)n(c7ccc(c8ccc(n9c(=O)n(c%10ccc(c%11ccc(n%12c(=O)n(c%13ccc2cc%13OC)c(=O)n(c%12=O)c2ccc(cc2OC)C)c(OC)c%11)cc%10OC)c(=O)n(c9=O)c2ccc(cc2OC)C)c(OC)c8)cc7OC)c(=O)n(c6=O)c2ccc(cc2OC)C)c(c5)OC)c(=O)n(c4=O)c2ccc(cc2OC)C)c(=O)n(c3=O)c2ccc(cc2OC)C)c(=O)n(c1=O)c1ccc(cc1OC)C";
    RWMol *m = SmilesToMol(smi,0,0);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==204);
    std::cerr<<"starting sanitization"<<std::endl;
    MolOps::sanitizeMol(*m);
    std::cerr<<"done"<<std::endl;    
    delete m;
  }


  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testSFNetIssue256() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 256: bad atom counts" << std::endl;

  {
    std::string smi="*CC[H]";
    ROMol *m = SmilesToMol(smi,0,0);
    TEST_ASSERT(m);
    m->updatePropertyCache(false);
    TEST_ASSERT(m->getNumAtoms()==4);
    TEST_ASSERT(m->getNumAtoms(false)==8);
    TEST_ASSERT(m->getNumHeavyAtoms()==2);    
    delete m;
  }
  {
    std::string smi="*CC[2H]";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==4);
    TEST_ASSERT(m->getNumAtoms(false)==8);
    TEST_ASSERT(m->getNumHeavyAtoms()==2);    
    delete m;
  }


  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testSFNetIssue266() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 266: ring finding error" << std::endl;

  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName+"Issue266.mol",false);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==19);
    TEST_ASSERT(m->getNumBonds()==25);
    std::cerr<<"starting ring finding"<<std::endl;
    MolOps::findSSSR(*m);
    std::cerr<<"done"<<std::endl;    
    TEST_ASSERT(m->getRingInfo()->numRings()==(m->getNumBonds()-m->getNumAtoms()+1));
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testSFNetIssue272() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing sf.net issue 272: removing two-coordinate Hs" << std::endl;

  {
    std::string smi="C[H-]C";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==3);
    delete m;
  }
  {
    std::string smi="C[H].C";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==2);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}



void testGitHubIssue8()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Github issue 8 (impact of removeAtom on bond stereo atoms)" << std::endl;
  {
    std::string smi= "Cl/C=C/Cl";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    MolOps::assignStereochemistry(*m);
    TEST_ASSERT(m->getBondWithIdx(1)->getStereoAtoms().size()==2);
    m->removeAtom((unsigned int)0);
    TEST_ASSERT(m->getBondWithIdx(1)->getStereoAtoms().size()==0);
    delete m;
  }
  {
    std::string smi= "CC/C=C/Cl";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    MolOps::assignStereochemistry(*m);
    INT_VECT &sas=m->getBondWithIdx(2)->getStereoAtoms();
    TEST_ASSERT(sas.size()==2);
    TEST_ASSERT(std::find(sas.begin(),sas.end(),1)!=sas.end());
    m->removeAtom((unsigned int)0);
    TEST_ASSERT(m->getBondWithIdx(1)->getStereoAtoms().size()==2);
    TEST_ASSERT(std::find(sas.begin(),sas.end(),0)!=sas.end());
    TEST_ASSERT(std::find(sas.begin(),sas.end(),1)==sas.end());
    delete m;
  }
  {
    std::string smi= "C/C=C/CC";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    MolOps::assignStereochemistry(*m);
    INT_VECT &sas=m->getBondWithIdx(1)->getStereoAtoms();
    TEST_ASSERT(sas.size()==2);
    TEST_ASSERT(std::find(sas.begin(),sas.end(),0)!=sas.end());
    TEST_ASSERT(std::find(sas.begin(),sas.end(),3)!=sas.end());
    m->removeAtom((unsigned int)4);
    TEST_ASSERT(m->getBondWithIdx(1)->getStereoAtoms().size()==2);
    TEST_ASSERT(std::find(sas.begin(),sas.end(),0)!=sas.end());
    TEST_ASSERT(std::find(sas.begin(),sas.end(),3)!=sas.end());
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGitHubIssue42()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Github issue 42 (impact of removeAtom on atom stereochem)" << std::endl;
  {
    std::string smi= "CCN1CCN(c2cc3[nH]c(C(=O)[C@@]4(CC)CC[C@](C)(O)CC4)nc3cc2Cl)CC1";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    int indices[]={29, 28, 27, 26, 25, 24, 8, 7, 6, 5, 4, 3, 2, 1, 0,-1};
    for(unsigned int i=0;indices[i]>-1;++i){
      m->removeAtom((unsigned int)indices[i]);
    }
    smi=MolToSmiles(*m,true);
    std::cerr<<"smiles: "<<smi<<std::endl;
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGitHubIssue65()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Github issue 65 (kekulization of boron-containing aromatic rings)" << std::endl;
  {
    std::string smi= "C[B-]1=CC=CC=C1";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getIsAromatic());
    TEST_ASSERT(m->getBondWithIdx(1)->getIsAromatic());

    MolOps::Kekulize(*m);

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGitHubIssue72()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Github issue 72 (problems with bad benzothiazolium structure)" << std::endl;
  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName+"github72.mol");
    TEST_ASSERT(m);
    TEST_ASSERT(!m->getBondBetweenAtoms(0,8)->getIsAromatic());
    TEST_ASSERT(m->getBondBetweenAtoms(1,6)->getIsAromatic());
    delete m;
  }

  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName+"github72.2.mol");
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(0,8)->getIsAromatic());
    TEST_ASSERT(m->getBondBetweenAtoms(1,6)->getIsAromatic());
    delete m;
  }
  {
    std::string pathName=getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName+"github72.3.mol");
    TEST_ASSERT(m);
    TEST_ASSERT(!m->getBondBetweenAtoms(0,8)->getIsAromatic());
    TEST_ASSERT(m->getBondBetweenAtoms(1,6)->getIsAromatic());

    std::string smi=MolToSmiles(*m,true);
    delete m;
    m = SmilesToMol(smi);
    TEST_ASSERT(m);

    
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

namespace{
  void _renumberTest(const ROMol *m){
    PRECONDITION(m,"no molecule");
    std::vector<unsigned int> idxV(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i) idxV[i]=i;

    std::string refSmi=MolToSmiles(*m,true);
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      std::vector<unsigned int> nVect(idxV);
      std::random_shuffle(nVect.begin(),nVect.end());
      //std::copy(nVect.begin(),nVect.end(),std::ostream_iterator<int>(std::cerr,", "));
      //std::cerr<<std::endl;

      ROMol *nm=MolOps::renumberAtoms(*m,nVect);
      TEST_ASSERT(nm);
      TEST_ASSERT(nm->getNumAtoms()==m->getNumAtoms());
      TEST_ASSERT(nm->getNumBonds()==m->getNumBonds());

      // checking the SSS is a test for Github #317
      MatchVectType mv; 
      TEST_ASSERT(SubstructMatch(*m,*nm,mv));
      TEST_ASSERT(mv.size()==nm->getNumAtoms());
      
      std::string nSmi=MolToSmiles(*nm,true);
      TEST_ASSERT(nSmi==refSmi);
      delete nm;
    }
  }
}

void testRenumberAtoms()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing renumbering atoms" << std::endl;
  {
    std::string smiles="CC1CCCC(C)C1C";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    _renumberTest(m);
    delete m;
  }
  {
    std::string smiles="C[C@H]1C[C@H](F)C1";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    _renumberTest(m);
    delete m;
  }
  {
    std::string smiles="C[C@H]1CC[C@H](C/C=C/[C@H](F)Cl)CC1";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    _renumberTest(m);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}
void testGithubIssue141()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing github issue 141: Kekulization of molecule with aromatic N leaves the explicit H there." << std::endl;
  {
    std::string smiles="N1C=CC=C1";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    MolOps::Kekulize(*m,true);
    m->updatePropertyCache(true);
    TEST_ASSERT(!m->getAtomWithIdx(0)->getIsAromatic())
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumImplicitHs()==1)
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumExplicitHs()==0)
    
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testZBO()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Testing ZBO basics" << std::endl;
  {
    RWMol *m=new RWMol();

    m->addAtom(new Atom(26));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addBond(1,2,Bond::AROMATIC);
    m->addBond(2,3,Bond::AROMATIC);
    m->addBond(3,4,Bond::AROMATIC);
    m->addBond(4,5,Bond::AROMATIC);
    m->addBond(5,6,Bond::AROMATIC);
    m->addBond(1,6,Bond::AROMATIC);

    m->addBond(1,0,Bond::ZERO);
    m->addBond(2,0,Bond::ZERO);
    m->addBond(3,0,Bond::ZERO);
    m->addBond(4,0,Bond::ZERO);
    m->addBond(5,0,Bond::ZERO);
    m->addBond(6,0,Bond::ZERO);

    MolOps::sanitizeMol(*m);

    TEST_ASSERT(m->getRingInfo()->numAtomRings(0)==0);
    TEST_ASSERT(m->getRingInfo()->numAtomRings(1)==1);

    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization()==Atom::SP2);
    TEST_ASSERT(m->getAtomWithIdx(2)->getHybridization()==Atom::SP2);
    TEST_ASSERT(m->getAtomWithIdx(3)->getHybridization()==Atom::SP2);
    TEST_ASSERT(m->getAtomWithIdx(4)->getHybridization()==Atom::SP2);
    TEST_ASSERT(m->getAtomWithIdx(5)->getHybridization()==Atom::SP2);
    TEST_ASSERT(m->getAtomWithIdx(6)->getHybridization()==Atom::SP2);

    TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
    TEST_ASSERT(m->getBondWithIdx(1)->getIsAromatic());
    TEST_ASSERT(m->getBondWithIdx(2)->getIsAromatic());
    TEST_ASSERT(m->getBondWithIdx(3)->getIsAromatic());
    TEST_ASSERT(m->getBondWithIdx(4)->getIsAromatic());
    TEST_ASSERT(m->getBondWithIdx(5)->getIsAromatic());
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testMolAssignment()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing operator= on molecules" << std::endl;
  {
    std::string smi= "CCN1CCN(c2cc3[nH]c(C(=O)[C@@]4(CC)CC[C@](C)(O)CC4)nc3cc2Cl)CC1";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    std::string csmi=MolToSmiles(*m,true);

    RWMol m2=*m;
    std::string nsmi=MolToSmiles(m2,true);
    TEST_ASSERT(nsmi==csmi);

    RWMol *m3 = SmilesToMol("C2CC2[C@H](F)Cl");
    TEST_ASSERT(m3);
    *m3 = *m;
    nsmi=MolToSmiles(*m3,true);
    TEST_ASSERT(nsmi==csmi);
    delete m3;
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGithubIssue190()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing github issue 190: Don't merge Hs onto dummy atoms." << std::endl;
  {
    std::string smiles="*[H]";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==2);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

namespace {
  int getAtNum(const ROMol &m,const Atom *at){
    return at->getAtomicNum();
  }
  std::string getSymbol(const ROMol &m,const Atom *at){
    return at->getSymbol();
  }
}
void testMolFragsWithQuery()
{
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing getMolFragsWithQuery()." << std::endl;
  {
    std::string smiles="C1CCC1ONNC";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==8);
    std::map<int,boost::shared_ptr<ROMol> > res=MolOps::getMolFragsWithQuery(*m,getAtNum);
    TEST_ASSERT(res.size()==3);
    TEST_ASSERT(res.find(6)!=res.end());
    TEST_ASSERT(res.find(7)!=res.end());
    TEST_ASSERT(res.find(8)!=res.end());
    TEST_ASSERT(res.find(5)==res.end());
    TEST_ASSERT(res[6]->getNumAtoms()==5);
    TEST_ASSERT(res[6]->getNumBonds()==4);
    TEST_ASSERT(res[7]->getNumAtoms()==2);
    TEST_ASSERT(res[7]->getNumBonds()==1);
    TEST_ASSERT(res[8]->getNumAtoms()==1);
    TEST_ASSERT(res[8]->getNumBonds()==0);
    delete m;
  }
  {
    std::string smiles="C1CCC1ONNC";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==8);
    std::vector<int> keep;
    keep.push_back(6);
    keep.push_back(8);
    std::map<int,boost::shared_ptr<ROMol> > res=MolOps::getMolFragsWithQuery(*m,getAtNum,true,&keep);
    TEST_ASSERT(res.size()==2);
    TEST_ASSERT(res.find(6)!=res.end());
    TEST_ASSERT(res.find(7)==res.end());
    TEST_ASSERT(res.find(8)!=res.end());
    TEST_ASSERT(res[6]->getNumAtoms()==5);
    TEST_ASSERT(res[6]->getNumBonds()==4);
    TEST_ASSERT(res[8]->getNumAtoms()==1);
    TEST_ASSERT(res[8]->getNumBonds()==0);
    delete m;
  }
  {
    std::string smiles="C1CCC1ONNC";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==8);
    std::vector<int> keep;
    keep.push_back(6);
    keep.push_back(8);
    std::map<int,boost::shared_ptr<ROMol> > res=MolOps::getMolFragsWithQuery(*m,getAtNum,true,&keep,true);
    TEST_ASSERT(res.size()==1);
    TEST_ASSERT(res.find(6)==res.end());
    TEST_ASSERT(res.find(7)!=res.end());
    TEST_ASSERT(res.find(8)==res.end());
    TEST_ASSERT(res[7]->getNumAtoms()==2);
    TEST_ASSERT(res[7]->getNumBonds()==1);

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}


int main(){
  RDLog::InitLogs();
  //boost::logging::enable_logs("rdApp.debug");

#if 1
  test1();
  test2();
  test3();
  test4();
  test5();
  //test6();
  test7();
  test8();
  test9();
  test10();
  test11();
  test12();
  testIssue183();
  testIssue188();
  testIssue189();
  testShortestPath();
  testIssue190();
  testIssue211();
  testIssue210();
  testIssue212();
  testAddHsCoords();
  testAddConformers();
  testSanitOps();
  testIssue252();
  testIssue276();
  testHsAndAromaticity();
  testSFIssue1694023();
  testSFIssue1719053();
  testSFIssue1811276();
  testSFIssue1836576();
  testChiralityAndRemoveHs();
  testSFIssue1894348();
  testSFIssue1942657();
  testSFIssue1968608();
  testHybridization();
  testAromaticityEdges();
  testSFNetIssue2196817();
  testSFNetIssue2208994();
  testSFNetIssue2313979();
  testSFNetIssue2316677();
  testSFNetIssue2951221();
  testSFNetIssue2952255();
  testSFNetIssue3185548();
  testSFNetIssue3349243();
  testFastFindRings();
  testSanitizeNonringAromatics();  
  testSFNetIssue3349243();
  testBasicCanon();
  testSFNetIssue3549146();
  testSFNetIssue249();
  testSFNetIssue256();
  testSFNetIssue266();
  testSFNetIssue266();
  testSFNetIssue272();
  testGitHubIssue8();
  testGitHubIssue42();
  testGitHubIssue65();
  testGitHubIssue72();
  testRenumberAtoms();
  testGithubIssue141();
#endif
  testZBO();
  testMolAssignment();

  testAtomAtomMatch();
  testGithubIssue190();
  testMolFragsWithQuery();
  return 0;
}


