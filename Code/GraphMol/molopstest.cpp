//  $Id$
// 
//   Copyright (C) 2002-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include <iostream>

using namespace RDKit;
using namespace std;
RWMol _t;
typedef class ROMol Mol;

void test1(){
  string smi;
  Mol *m;
  INT_VECT iv;
  int count;
  smi = "CCCC(=O)O";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m,"smiles parse failed");
  count = MolOps::getMolFrags(*m,iv);
  CHECK_INVARIANT(count==1,"bad frag count");
  delete m;

  smi = "CCCC(=O)[O-].[Na+]";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m,"smiles parse failed");
  count = MolOps::getMolFrags(*m,iv);
  CHECK_INVARIANT(count==2,"bad frag count");
  delete m;

  smi = "CCCC(=O)[O-].[Na+].[NH4+].[Cl-]";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m,"smiles parse failed");
  count = MolOps::getMolFrags(*m,iv);
  CHECK_INVARIANT(count==4,"bad frag count");
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
  delete m2;
  //BOOST_LOG(rdInfoLog) << "7" << std::endl;
  // remove all after removing only implicit
  m2 = MolOps::removeHs(*m3,false);
  CHECK_INVARIANT(m2->getNumAtoms()==4,"");

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
  sma = GetAtomSmarts(static_cast<const QueryAtom *>(m->getAtomWithIdx(0)));
  TEST_ASSERT(sma=="C");

  delete m2;
  //BOOST_LOG(rdInfoLog) << "16" << std::endl;
  m2 = MolOps::addHs(*m);
  TEST_ASSERT(m2->getNumAtoms()==8);
  sma = GetAtomSmarts(static_cast<const QueryAtom *>(m2->getAtomWithIdx(0)));
  TEST_ASSERT(sma=="C");

  delete m;
  //BOOST_LOG(rdInfoLog) << "17" << std::endl;
  m = MolOps::mergeQueryHs(*m2);
  TEST_ASSERT(m->getNumAtoms()==2);
  sma = GetAtomSmarts(static_cast<const QueryAtom *>(m->getAtomWithIdx(0)));
  TEST_ASSERT(sma=="[C&!H0]");
  


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

  MolOps::assignAtomChiralCodes(*m);

  int cip1,cip2;
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("_CIPRank"));
  m->getAtomWithIdx(0)->getProp("_CIPRank",cip1);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPRank"));
  m->getAtomWithIdx(2)->getProp("_CIPRank",cip2);
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
  MolOps::assignAtomChiralCodes(*m);
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

#if 1
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignAtomChiralCodes(*m,true);
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
  MolOps::assignAtomChiralCodes(*m,true);
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
  MolOps::assignAtomChiralCodes(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()==Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(!(m->getAtomWithIdx(1)->hasProp("_CIPCode")));

  delete m;
  smi = "F[C@]1C(Cl)CC1";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignAtomChiralCodes(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));


  
  delete m;
  smi = "F[C@@](C)(Cl)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignAtomChiralCodes(*m,true);
  TEST_ASSERT(!(m->getAtomWithIdx(0)->hasProp("_CIPCode")));
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  TEST_ASSERT(!(m->getAtomWithIdx(2)->hasProp("_CIPCode")));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete m;
  smi = "F[C@](Br)(C)Cl";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignAtomChiralCodes(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  delete m;
  smi = "F[C@](Cl)(Br)C";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignAtomChiralCodes(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete m;
  smi = "FC(F)(F)[C@](Br)(F)C(Cl)(Cl)Cl";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignAtomChiralCodes(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(4)->hasProp("_CIPCode"));
  m->getAtomWithIdx(4)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete m;
  smi = "C[C@](C=C)(F)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignAtomChiralCodes(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete m;
  smi = "CC[C@](C=C)(F)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignAtomChiralCodes(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");
  delete m;
  smi = "[CH2-][C@](C)(F)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignAtomChiralCodes(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete m;
  smi = "F[C@]([H])(Cl)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignAtomChiralCodes(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  delete m;
  smi = "F[C@H](Cl)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignAtomChiralCodes(*m,true);
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
  MolOps::assignAtomChiralCodes(*m,true);
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
  MolOps::assignAtomChiralCodes(*m,true);
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
  MolOps::assignAtomChiralCodes(*m,true);
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
  MolOps::assignAtomChiralCodes(*m,true);
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
  MolOps::assignAtomChiralCodes(*m,true);
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
  MolOps::assignAtomChiralCodes(*m,true);
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
  MolOps::assignAtomChiralCodes(*m,true);
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
  MolOps::assignAtomChiralCodes(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");
  TEST_ASSERT(m->getAtomWithIdx(4)->hasProp("_CIPCode"));
  m->getAtomWithIdx(4)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  // this is Issue 152:
  delete m;
  smi = "C1[C@H](N)C[C@H](C)C=1";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignAtomChiralCodes(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("_CIPCode"));
  m->getAtomWithIdx(1)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
  TEST_ASSERT(m->getAtomWithIdx(4)->hasProp("_CIPCode"));
  m->getAtomWithIdx(4)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="S");

  delete m;
#endif
  
  // -----------------------------------------------
  // these are related to Issue 397:
  smi = "C(=O)[C@@H](C)N";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignAtomChiralCodes(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  smi = "C(=O)[C@@H](CO)N";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignAtomChiralCodes(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");

  smi = "C(O)[C@@H](C)N";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------" << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  MolOps::assignAtomChiralCodes(*m,true);
  TEST_ASSERT(m->getAtomWithIdx(2)->hasProp("_CIPCode"));
  m->getAtomWithIdx(2)->getProp("_CIPCode",cip);
  TEST_ASSERT(cip=="R");
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
  //MolOps::assignBondStereoCodes(*m);
  TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
  TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
  TEST_ASSERT(m->getBondWithIdx(2)->getStereo() == Bond::STEREONONE);

  delete m;
  smi = "F/C=CCl";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  //MolOps::assignBondStereoCodes(*m);
  TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

  delete m;
  smi = "F/C=C/Cl";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  //MolOps::assignBondStereoCodes(*m);
  TEST_ASSERT(m->getBondWithIdx(1)->getStereo()==Bond::STEREOE);

  delete m;
  smi = "F/C=C(/Br)Cl";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  //MolOps::assignBondStereoCodes(*m);
  TEST_ASSERT(m->getBondWithIdx(1)->getStereo()==Bond::STEREOE);


  delete m;
  smi = "F/C=C(/Cl)Br";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  //MolOps::assignBondStereoCodes(*m);
  TEST_ASSERT(m->getBondWithIdx(1)->getStereo()==Bond::STEREOZ);


  delete m;
  smi = "F/C(Br)=C/Cl";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  //MolOps::assignBondStereoCodes(*m);
  TEST_ASSERT(m->getBondWithIdx(2)->getStereo()==Bond::STEREOZ);

  delete m;
  smi = "F/C=C(/Cl)Cl";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  //MolOps::assignBondStereoCodes(*m);
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
  MolOps::assignBondStereoCodes(*m2);
  TEST_ASSERT(m2->getBondWithIdx(1)->getStereo()==Bond::STEREOE);


  m2->getBondWithIdx(0)->setBondDir(Bond::ENDDOWNRIGHT);
  MolOps::assignBondStereoCodes(*m2,true,true);
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
  MolOps::assignBondStereoCodes(*m2);
  TEST_ASSERT(m2->getBondWithIdx(1)->getStereo()==Bond::STEREOZ);

  m2->getBondWithIdx(0)->setBondDir(Bond::ENDDOWNRIGHT);
  MolOps::assignBondStereoCodes(*m2,true,true);
  TEST_ASSERT(m2->getBondWithIdx(1)->getStereo()==Bond::STEREOE);


  // ----------------------
  // test Issue 174:
  delete m2;
  smi = "O\\N=C\\C=N/O";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  //MolOps::assignBondStereoCodes(*m2);
  TEST_ASSERT(m2->getBondWithIdx(1)->getStereo()==Bond::STEREOE);
  TEST_ASSERT(m2->getBondWithIdx(3)->getStereo()==Bond::STEREOZ);
  refSmi = MolToSmiles(*m2,1);
  m = SmilesToMol(refSmi);
  TEST_ASSERT(m);
  //MolOps::assignBondStereoCodes(*m);
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
  TEST_ASSERT(m2->getBondWithIdx(4)->getStereo()==Bond::STEREOE);
  TEST_ASSERT(m2->getBondWithIdx(10)->getStereo()==Bond::STEREOZ);

  refSmi = MolToSmiles(*m2,1);
  //BOOST_LOG(rdInfoLog) << "ref: " << refSmi << std::endl;
  m = SmilesToMol(refSmi);
  TEST_ASSERT(m);
  smi = MolToSmiles(*m,1);
  //BOOST_LOG(rdInfoLog) << "smi: " << smi << std::endl;
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
  std::string smi ="CC(OC1C(CCCC3)C3C(CCCC2)C2C1OC(C)=O)=O";
  ROMol *m = SmilesToMol(smi);

  INT_LIST path = MolOps::getShortestPath(*m, 1, 20);
  CHECK_INVARIANT(path.size() == 5, "");
  INT_LIST_CI pi = path.begin();
  CHECK_INVARIANT((*pi) == 2, ""); pi++;
  CHECK_INVARIANT((*pi) == 3, ""); pi++;
  CHECK_INVARIANT((*pi) == 16, ""); pi++;
  CHECK_INVARIANT((*pi) == 17, ""); pi++;
  CHECK_INVARIANT((*pi) == 18, "");
  delete m;
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
  
  
  
  
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testAddConformers() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Add Confomers" << std::endl;

  std::string smi = "CC";
  ROMol *m = SmilesToMol(smi);
  int i;
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

  
  unsigned int na = mol->getNumAtoms();
  unsigned int nb = mol->getNumBonds();
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

int main(){
  RDLog::InitLogs();
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
#endif
  
  return 0;
}


