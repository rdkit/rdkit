// $Id$
//
//  Copyright (C) 2003-2011 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>

#include <fstream>
#include "SmilesParse.h"
#include "SmilesWrite.h"
#include "SmartsWrite.h"
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <RDGeneral/RDLog.h>
//#include <boost/log/functions.hpp>
using namespace RDKit;
using namespace std;
typedef ROMol Mol;

#define LOCAL_TEST_ALL 1
#if LOCAL_TEST_ALL
void testPass(){
  int i = 0;
  Mol *mol;
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing patterns which should parse." << std::endl;
  string smis[]={
#if 1
    "C",
    "CC",
    "C-C",
    "C=C",
    "[CH2+]C[CH+2]",
#endif
    "C1CC=1",
    "C=1CC1",
    "Ccc",
    "C=C-O",
    "C1CC1",
    "C1NC1",
    "C1=CC1",
    "C1CCC1",
    "CC(C)CC",
    "CC(=O)O",
    "C1C(=O)C1",
    "C1C(N)C1",
    "CC(O)C",
    "OC=CCC",
    "CC([O-])O",
    "C1CC2C1CC2",
    "Cl/C=C/Cl",
    "Cl/C=C\\Cl",
    "Cl/C=C/Cl",
    "Cl/C=C\\Cl",
    "C1CC.CC1",
    "C1C(C2CC2).C2CC2C1",
    "[Na+].[Cl-].[NH4+].[Cl-]",
    "[$(CO)]CO",
    "[!$(*#*)]",
    "[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]",
    "[C;!$(C-[OH])]=O",
    "[CH2;X4][N,O,S,F,Cl,Br,I]",
    "C=O",
    "[C;!$(C-[OH])]=O",
    "[#6]-!:[#6]",
    "[C^3]",
    "[se]",
    "[te]",
    // test zeros as ring indices, issue 2690982:
    "C0CC0", 
    // these used to fail before Roger Sayle's SMARTS patch:
    "[HO]",
    "[13]",
    "[+]",
    "[C:1]",
    "[C:0]",  // issue 3525776
    "[$([#0].[#0])]", // sf.net issue 281
    //"CC((C)C)", // github #102
    "c1ccccb1", // github 220
    "EOS"};
  while( smis[i] != "EOS" ){
    string smi = smis[i];
    mol = SmartsToMol(smi);
    CHECK_INVARIANT(mol,smi);
    if (mol) {
      int nAts = mol->getNumAtoms();
      CHECK_INVARIANT(nAts!=0,smi.c_str());
      delete mol;
    }
    i++;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testFail(){
  int i = 0;
  Mol *mol;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing patterns which should fail to parse." << std::endl;
  BOOST_LOG(rdInfoLog) << "\tExpect Parse error messages" << std::endl;

  // alternate good and bad smiles here to ensure that the parser can resume parsing
  // on good input:
  string smis[]={
    "CC=(CO)C",
    "CC(=CO)C",
    "C1CC",
    "C1CC1",
    "fff",
    "C1CC1",
    "C=0", // part of sf.net issue 2525792
    "C1CC1",
    "C0", // part of sf.net issue 2525792
    "C1CC1",
    "C-0", // part of sf.net issue 2525792
    "C1CC1",
    "C+0", // part of sf.net issue 2525792
    "C1CC1",
    "[HQ]",
    "C1CC1",
    "EOS"};
  while( smis[i] != "EOS" ){
    string smi = smis[i];
    boost::logging::disable_logs("rdApp.error");
    mol = SmartsToMol(smi);
    boost::logging::enable_logs("rdApp.error");
    if(!(i%2)) {
      CHECK_INVARIANT(!mol,smi);
    }
    else{
      CHECK_INVARIANT(mol,smi);
    }
    i++;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;

}

std::vector< MatchVectType > _checkMatches(std::string smarts, std::string smiles, 
                   unsigned int nMatches, unsigned int lenFirst, bool addHs=false) {
  // utility function that will find the matches between a smarts and smiles 
  // if they match the expected values
  //  smarts : smarts string
  //  smiles : smiles string
  //  nMatches : expected number of matches
  //  lenFirst : length of the first match
  //
  // Return the list of all matches just in case want to do aditional testing
  ROMol *mol, *mol2, *matcher, *matcher2;
  bool matches;
  unsigned int matchCount;
  std::string pickle;
  MatchVectType mV;
  std::vector< MatchVectType > mVV;

  matcher = SmartsToMol(smarts);
  CHECK_INVARIANT(matcher,smarts);
  // we will at the same time test the serialization:
  MolPickler::pickleMol(matcher,pickle);
  matcher2 = new ROMol();
  MolPickler::molFromPickle(pickle,matcher2);
  CHECK_INVARIANT(matcher2,smarts);

  //std::cerr << "\tSMA: " << smarts << " -> " << MolToSmarts(*matcher) << std::endl;;
  
  mol = SmilesToMol(smiles);
  CHECK_INVARIANT(mol,smiles);
  if (addHs) {
    mol2 = MolOps::addHs(*mol);
    delete mol;
    mol = mol2;
  }
  MolOps::findSSSR(*mol);
  
  matches = SubstructMatch(*mol,*matcher,mV);
  CHECK_INVARIANT(matches, smarts + " " + smiles);
  CHECK_INVARIANT(mV.size()==lenFirst, smarts + " " + smiles);
  matchCount = SubstructMatch(*mol,*matcher,mVV,true);
  CHECK_INVARIANT(matchCount==nMatches, smarts + " " + smiles);
  CHECK_INVARIANT(mVV[0].size()==lenFirst, smarts + " " + smiles);
  delete matcher;
  matcher=0;

  matches = SubstructMatch(*mol,*matcher2,mV);
  CHECK_INVARIANT(matches, smarts + " " + smiles);
  CHECK_INVARIANT(mV.size()==lenFirst, smarts + " " + smiles);
  matchCount = SubstructMatch(*mol,*matcher2,mVV,true);
  CHECK_INVARIANT(matchCount==nMatches, smarts + " " + smiles);
  CHECK_INVARIANT(mVV[0].size()==lenFirst, smarts + " " + smiles);
  delete matcher2;
  matcher2=0;

  delete mol;
  
  return mVV;
}

void _checkNoMatches(std::string smarts, std::string smiles, bool addHs=false) {
  ROMol *mol,*matcher,*matcher2;
  std::string pickle;
  bool matches;
  MatchVectType mV;

  matcher = SmartsToMol(smarts);
  CHECK_INVARIANT(matcher,smarts);
  // we will at the same time test the serialization:
  MolPickler::pickleMol(matcher,pickle);
  matcher2 = new ROMol();
  MolPickler::molFromPickle(pickle,matcher2);
  CHECK_INVARIANT(matcher2,smarts);

  mol = SmilesToMol(smiles);
  CHECK_INVARIANT(mol,smiles);
  if (addHs) {
    ROMol *mol2 = MolOps::addHs(*mol);
    delete mol;
    mol = mol2;
  }
  MolOps::findSSSR(*mol);

  matches = SubstructMatch(*mol,*matcher,mV);
  CHECK_INVARIANT(!matches,smarts+"|"+smiles);
  matches = SubstructMatch(*mol,*matcher2,mV);
  CHECK_INVARIANT(!matches,smarts+"|"+smiles);
  delete mol;
  delete matcher;
  delete matcher2;
}

void testMatches(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing matching" << std::endl;

  _checkMatches("C#,=O", "CC(=O)O", 1, 2);
  _checkNoMatches("C#,=O", "CC(O)O");

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}
    
void testMatches2(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing matching2" << std::endl;

  _checkMatches("C#*", "C#CO", 1, 2);
  
  _checkMatches("[D3]O", "C1C(O)CCC1O", 2, 2);
  _checkNoMatches("[D3]O", "C1CCCC1");

  _checkMatches("[h3]C", "C1C(C)CCC1C", 2, 2);
  _checkNoMatches("[h3]C", "C1CCCC1");

  _checkMatches("[D1;h2]C", "C1C(C)CCC1N", 1,2);

  _checkMatches("[D1,h2][#6]", "c1c(C)cccc1CC", 3, 2);

  _checkMatches("[R2]","C1CC2C1CC2", 2, 1);

  _checkMatches("[r4]", "C1CC2CCCC12", 4, 1);

  _checkMatches("[!r4]", "C1CC2CCCC12", 3, 1);

  _checkMatches("C@C", "C1CC1CC", 3, 2);
  _checkNoMatches("C@C", "CCCCC");

  _checkMatches("C!@C", "C1CC1CC", 2, 2);

  _checkMatches("C!@C", "CCCCC", 4, 2);

  _checkMatches("[x2]", "C1CCC1", 4, 1);
  _checkNoMatches("[x1]", "C1CCC1");
  _checkNoMatches("[x3]", "C1CCC1");
  _checkMatches("[x2]", "C1CC2C1CC2", 4, 1);
  _checkMatches("[x3]", "C1CC2C1CC2", 2, 1);
  _checkNoMatches("[x0]", "C1CCC1");
  _checkMatches("[x0]", "CC1CCC1", 1, 1);

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}


void testMatches3(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing matching3" << std::endl;

  _checkMatches("[#6]!@[!#6]", "C1CC1CO", 1, 2);

  _checkMatches("[#6]!@[!#6]", "CCOCC", 2, 2);

  _checkMatches("C[O-]", "CC(=O)[O-]", 1, 2);
  _checkNoMatches("C[O-]", "CC(=O)[S-]");

  _checkMatches("C[-]", "CC(=O)[O-]", 1, 2);

  _checkMatches("C[-]", "CC(=O)[S-]", 1, 2);

  _checkMatches("[CH3][CH]", "C(C)(C)CC(=O)[O-]", 2, 2);

  _checkMatches("[c;H]", "c1ccccc1", 6, 1);
  _checkNoMatches("[c;H]", "C1CCCCC1");

  _checkMatches("[#6]([#7])[#6]", "c1ncccc1", 2, 3);

  _checkMatches("[c;H]", "c1c[c-]ccc1", 5, 1);

  _checkMatches("C~O", "C(=O)[O-]", 2, 2);


  // -----
  // This block is connected to SF-Issue 1538280
  //   http://sourceforge.net/tracker/index.php?func=detail&aid=1538280&group_id=160139&atid=814650
  _checkMatches("[R]", "c1ccccc1", 6, 1);

  _checkMatches("[r]", "c1ccccc1", 6, 1);

  _checkMatches("[R;!O]", "c1ccccc1", 6, 1);

  _checkMatches("[c]", "c1ccccc1", 6, 1);

  _checkMatches("[R;c]", "c1ccccc1", 6, 1);

  _checkMatches("[c;#6]", "c1ccccc1", 6, 1);

  _checkMatches("[R;R]", "c1ccccc1", 6, 1);

  _checkMatches("[#6;R]", "c1ccccc1", 6, 1);

  _checkMatches("[c;R]", "c1ccccc1", 6, 1);

  _checkMatches("[!O;R]", "c1ccccc1", 6, 1);

  _checkMatches("[!C;R]", "c1ccccc1", 6, 1);

  _checkMatches("[!C;R]", "C1COC1", 1, 1);


  // -----
  // This block is connected to SF-Issue 1836223
  //   http://sourceforge.net/tracker/index.php?func=detail&aid=1836223&group_id=160139&atid=814650
  _checkMatches("[r6]", "C1CCCCC1", 6, 1);
  _checkNoMatches("[CR2r6]", "C1CCCC2C1CC2");
  _checkMatches("[CR2r4]", "C1CCCC2C1CC2",2,1);

  // -----
  // This block is connected to SF-Issue 1912895
  //   http://sourceforge.net/tracker/index.php?func=detail&aid=1912895&group_id=160139&atid=814650
  _checkMatches("[$(Sc)]", "Sc1ccccc1", 1, 1);
  _checkMatches("[Sc]", "[Sc]C", 1, 1);
  _checkMatches("[Sc]", "[Sc]C", 1, 1);
  _checkMatches("[$([Sc])]", "[Sc]C",1,1);
  _checkNoMatches("[$(Sc)]", "[Sc]C");
  _checkNoMatches("[$([Sc])]", "Sc1ccccc1");


  // -----
  // This block is connected to SF-Issue 1968930
  //   http://sourceforge.net/tracker/index.php?func=detail&aid=1968930&group_id=160139&atid=814650
  _checkNoMatches("C[16*,32*]", "CC(=O)[O-]");
  _checkMatches("C[16*,32*]", "CC(=O)[16O-]",1,2);
  _checkNoMatches("C[16*,32*]", "CC(=O)[S-]");
  _checkMatches("C[16*,32*]", "CC(=O)[32S-]", 1, 2);

  // -----
  // This block is connected to GitHub #60
  //
  _checkMatches("[#7h1]","c1cnc[nH]1", 1,1);
  _checkNoMatches("[#7h1]","c1cnc[nH]1",true);

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testMatches4(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing matching4" << std::endl;

  _checkMatches("[X3]", "C(C)(C)=CC(C)(C)N", 3, 1);
  
  _checkMatches("[D3]", "C(C)(C)=CC(C)(C)N", 1, 1);
  
  _checkMatches("[v3]", "C(C)(C)=CC(C)(C)N", 1, 1);
  
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testMatches5(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing matching 5" << std::endl;

  _checkMatches("[$(CO)]", "COCCC", 2, 1);
  
  _checkMatches("[$(CO)]C", "COCCC", 1, 2);
  
  std::vector< MatchVectType > mVV = _checkMatches("[C;D3](C)C", "CC(C)CCOC(C)C", 4, 3);
  CHECK_INVARIANT(mVV[0][0].second == 1, "");
  CHECK_INVARIANT(mVV[0][1].second==0||mVV[0][2].second==0,"");
  CHECK_INVARIANT(mVV[0][1].second==2||mVV[0][2].second==2,"");
  
  mVV = _checkMatches("[$(C(C)C)]O", "CC(C)CCOC(C)C", 1, 2);
  CHECK_INVARIANT(mVV[0][0].second==6,"");
  CHECK_INVARIANT(mVV[0][1].second==5,"");

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}
void testMatches6(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing matching 6" << std::endl;

  _checkMatches("[C^3]", "CCC", 3, 1);
  _checkMatches("[C^3]", "CC=C", 1, 1);
  _checkMatches("[C^2]", "CC=C", 2, 1);
  _checkNoMatches("[C^2]", "CCC");
  
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testProblems() 
{
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing former problems" << std::endl;

#if 1
  _checkMatches("[$(C(=O)O)]", "CC(=O)OC", 1, 1);

  _checkMatches("[$(C(=O)[O,N])]", "CC(=O)OC", 1, 1);

  _checkNoMatches("[$(C(=O)O)]", "CC(=O)NC");
  _checkMatches("[$(C(=O)[O,N])]", "CC(=O)NC", 1, 1);

  //
  //  Problems discovered while getting the Crippen stuff working
  //
  _checkNoMatches("[N+0]C", "C[N+](=O)[O-]");

  _checkMatches("[N+1]C", "C[N+](=O)[O-]", 1, 2);

  // Issue 71
  _checkMatches("CCC", "CCC", 1, 3);
  _checkMatches("CCC", "C1CC1", 1, 3);

  // this pattern caused crashes at one point
  _checkNoMatches("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]", "N#CCC#N");
  
  // ensure that the recursive queries get property reset
  _checkMatches("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]", "OCCC#N", 1, 2);
  
  _checkNoMatches("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]", "N#CCC#N");

  _checkMatches("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]", "OCCC#N", 1, 2);

  // these weren't being parsed at one point, make sure they match
  _checkMatches("N(=,-C)", "N=C", 1, 2);
  _checkMatches("N(=,-C)", "N-C", 1, 2);
  _checkNoMatches("N(=,-C)", "N#C");
  
  _checkMatches("[$(O),Cl]", "C(=O)O", 2, 1);
  _checkMatches("[Cl,$(O)]", "C(=O)O", 2, 1);
  
  // tests precedence of &
  _checkNoMatches("[N&v3;H1,H2]", "CCC");
  _checkMatches("[N&v3;H1,H2]", "CNC", 1, 1);

  // nested recursion (yick!)
  _checkNoMatches("[O]-[!$(*=O)]", "CC(=O)O");

  //BOOST_LOG(rdInfoLog) << "-*-*-*-*-*-*-*-*-" << std::endl;
#endif
  _checkNoMatches("[$([O]-[!$(*=O)])]", "CC(=O)O");
  
  // ISSUE 78
  _checkNoMatches("[$(C1CC1)]", "C1CCC1");
  
  // The next one got fixed too quickly to get an Issue id:
  _checkMatches("[!$(a(:a!:*):a!:*)]", "c1(C)cccc(C)c1CC", 9, 1);
  
  // Issue 99:
  std::string sma = "[r10,r11]";
  ROMol *matcher = SmartsToMol(sma);
  CHECK_INVARIANT(matcher,sma);

  // Issue 65:
  _checkMatches("C[C;H]", "FC(F)(F)C(O)O", 1, 2);
  _checkNoMatches("CC[H]", "FC(F)(F)C(O)O");
  _checkMatches("CC[H]", "FC(F)(F)C(O)O", 1, 3, true);
  
  _checkNoMatches("C[C;H]", "FC(F)(F)C(O)(O)O");
  _checkNoMatches("CC[H]", "FC(F)(F)C(O)(O)O");

  _checkNoMatches("C[C;H]", "FC(F)(F)CO");
  _checkNoMatches("CC[H]", "FC(F)(F)CO");
  _checkMatches("CC[H]", "FC(F)(F)CO", 2, 3, true);
  
  _checkMatches("*-[N;H2,H1&-1,-2]", "CC(N)C", 1, 2);
  _checkMatches("*-[N;H2,H1&-1,-2]", "CC([NH2])C", 1, 2);
  _checkMatches("*-[N;H2,H1&-1,-2]", "CC([NH-])C", 1, 2);
  _checkMatches("*-[N;H2,H1&-1,-2]", "CC([N-2])C", 1, 2);
  _checkNoMatches("*-[N;H2,H1&-1,-2]", "CC(=N)C");
  _checkNoMatches("*-[N;H2,H1&-1,-2]", "CC(NC)C");
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}


void testFrags(){
  int i = 0;
  Mol *mol;
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing fragment patterns." << std::endl;
  string smis[]={
    "C=O",
    "[C!$(C-[OH])]=O",
    "[C!$(C=O)]-[OH]",
    "c[OH]",
    "O(-[#6])-C",
    "[NH2]",
    "[NH1,nH1]",
    "[NH0,nH0]",
    "n",
    "[Nv3](=C)-[#6]",
    "C#N",
    "[#9,#17,#35,#53]",
    "[SX2](-[#6])-C",
    "[SH]",
    "C=[SX1]",
    "C=N-O",
    "[N!$(N=O)](-O)-C",
    "[N!$(N-O)]=O",
    "[NX3]-[NX3]",
    "C=N-[NX3]",
    "N(=O)(O)[#6]",
    "[#6]-N=N-[#6]",
    "[N+]#N",
    "[#6]-N=[N+]=[N-]",
    "S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])(-[#6])-[#6]",
    "N-S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])-[#6]",
    "[NH2]-S(=,-[OX1;+0;-1])(=,-[OX1;+0;-1])-[#6]",
    "C(=O)-N",
    "C(=O)-[NH2]",
    "C(=N)(-N)-[!#7]",
    "N=C=O",
    "N=C=S",
    "S-C#N",
    "C(=O)O-C",
    "C-C(=O)[O;H1,-]",
    "c-C(=O)[O;H1,-]",
    "[#6]C(=O)[O;H,-1]",
    "C1C(=O)NC(=O)NC1=O",
    "C(=O)(-N)-N",
    "N(-C(=O))-C=O",
    "C#[CH]",
    "n1cncc1",
    "o1cccc1",
    "s1cccc1",
    "c1scnc1",
    "c1ocnc1",
    "n1ccccc1",
    "N1CCCCC1",
    "N1CCNCC1",
    "O1CCNCC1",
    "N1C(=O)CC1",
    "[NX4]",
    "[nH]",
    "C(=N)(N)N",
    "c1nnnn1",
    "O1CC1",
    "EOS"};
  while( smis[i] != "EOS" ){
    string smi = smis[i];
    mol = SmartsToMol(smi);
    CHECK_INVARIANT(mol,smi);
    if (mol) {
      int nAts = mol->getNumAtoms();
      CHECK_INVARIANT(nAts!=0,smi.c_str());
      delete mol;
    }
    i++;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testSmartsWrite() {
  int i = 0;
  Mol *mol;
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Smarts Writer." << std::endl;
  string smis[]={
    "[v3]",
    "n1cncc1",
    "c[OH]",
    "S(=,-[O])",
    "CC",
    "C=O",
    "[$(C=O)]",
    "[C!$(C=O)]",
    "[C!$(C-[OH])]=O",
    "[C!$(C=O)]-[OH]",
    "O(-[#6])-C",
    "[NH2]",
    "[NH1,nH1]",
    "[NH0,nH0]",
    "n",
    "[Nv3](=C)-[#6]",
    "C#N",
    "[#9,#17,#35,#53]",
    "[SX2](-[#6])-C",
    "[SH]",
    "C=[SX1]",
    "C=N-O",
    "[N!$(N=O)](-O)-C",
    "[N!$(N-O)]=O",
    "[NX3]-[NX3]",
    "C=N-[NX3]",
    "N(=O)(O)[#6]",
    "[#6]-N=N-[#6]",
    "[N+]#N",
    "[#6]-N=[N+]=[N-]",
    "S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])(-[#6])-[#6]",
    "N-S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])-[#6]",
    "[NH2]-S(=,-[OX1;+0;-1])(=,-[OX1;+0;-1])-[#6]",
    "C(=O)-N",
    "C(=O)-[NH2]",
    "C(=N)(-N)-[!#7]",
    "N=C=O",
    "N=C=S",
    "S-C#N",
    "C(=O)O-C",
    "C-C(=O)[O;H1,-]",
    "c-C(=O)[O;H1,-]",
    "[#6]C(=O)[O;H,-1]",
    "C1C(=O)NC(=O)NC1=O",
    "C(=O)(-N)-N",
    "N(-C(=O))-C=O",
    "C#[CH]",
    "n1cncc1",
    "o1cccc1",
    "s1cccc1",
    "c1scnc1",
    "c1ocnc1",
    "n1ccccc1",
    "N1CCCCC1",
    "N1CCNCC1",
    "O1CCNCC1",
    "N1C(=O)CC1",
    "[NX4]",
    "[nH]",
    "C(=N)(N)N",
    "c1nnnn1",
    "O1CC1",
    "[C^3]",
    "[$([NH2][CX4]),$([$([NH]([CX4])[CX4]);!$([NH]([CX4])[CX4][O,N]);!$([NH]([CX4])[CX4][O,N])]),$([ND3]([CX4])([CX4])[CX4])]", // this was sf.net issue 3496800
    "[r6][r6]", // this was sf.net issue 3496799
    "EOS"};
  
  std::vector<std::string> diffSmi;

  while( smis[i] != "EOS" ){
    std::string smi = smis[i];
    mol = SmartsToMol(smi);
    CHECK_INVARIANT(mol,smi);
    std::string nsma = MolToSmarts(*mol);
    if (smi != nsma) {
      diffSmi.push_back(smi);
    }
    i++;
    delete mol;
  }

  // for the smarts that come out different from the writer verify if they are 
  // functionally the same
  
  std::string smiles[] = {
    "c1c(O)cccc1",
    "O=C(O)C[N+]#N",
    "CN=[N+]=[N-]",
    "NS(O)(=O)C",
    "NS(=O)(=O)C",
    "CC1=CC(=O)C=CC1=O",
    "OC1=C(Cl)C=C(C=C1[N+]([O-])=O)[N+]([O-])=O",
    "NC1=CC2=C(C=C1)C(=O)C3=C(C=CC=C3)C2=O",
    "NC1=CC2=C(C=C1)C(=O)C3=C(C=CC=C3)C2=O",
    "OC(=O)C1=C(C=CC=C1)C2=C3C=CC(=O)C(=C3OC4=C2C=CC(=C4Br)O)Br",
    "CN(C)C1=C(Cl)C(=O)C2=C(C=CC=C2)C1=O",
    "CC(=NO)C(C)=NO",
    "CC1=NN(C(=O)C1)C2=CC=CC=C2",
    "NC1=CC=NC2=C1C=CC(=C2)Cl",
    "BrN1C(=O)CCC1=O",
    "CCCCSCC",
    "CC(=O)NC1=NC2=C(C=C1)C(=CC=N2)O",
    "CC(=O)NC1=NC2=C(C=C1)C(=CC=N2)O",
    "NN=C(C1=CC=CC=C1)C2=CC=CC=C2",
    "OC(=O)[CH](CC1=CC=CC=C1)C2=CC=CC=C2",
    "O=S(CC1=CC=CC=C1)CC2=CC=CC=C2",
    "O=S(=O)(CC1=CC=CC=C1)CC2=CC=CC=C2",
    "O=S(=O)(CC1=CC=CC=C1)CC2=CC=CC=C2",
    "CCOC(=O)C1=CC=C(C=C1)S(=O)(=O)N(CC)CC",
    "CN1C2=C(SC3=C1C=CC=C3)C=CC=C2",
    "CC1=CC=C(S)C=C1",
    "SC1=C2C=CC=C(S)C2=CC=C1",
    "CC(=S)N1CCCCC1",
    "CC(C)CSC(C)=S",
    "NNP(=S)(NN)C1=CC=CC=C1",
    "NNC(=S)NNC1=CC=CC=C1",
    "CCCCOC(N)=O",
    "CNCC(N)=O",
    "CC(C)(O)C#C",
    "CC(C)C[C](C)(O)C#C",
    "EOS"};

  std::vector<std::string>::const_iterator dsmi;
  i = 0;
  MatchVectType mV1, mV2;
  bool mts1, mts2;
  while (smiles[i] != "EOS") {
    ROMol *nmol = SmilesToMol(smiles[i]);
    TEST_ASSERT(nmol);
    for (dsmi = diffSmi.begin(); dsmi != diffSmi.end(); dsmi++) {
      ROMol *m1 = SmartsToMol(*dsmi);
      TEST_ASSERT(m1)
      std::string wsma = MolToSmarts(*m1);
      ROMol *m2 = SmartsToMol(wsma);
      TEST_ASSERT(m2)

      mts1 = SubstructMatch(*nmol,*m1,mV1);
      mts2 = SubstructMatch(*nmol,*m2,mV2);

      CHECK_INVARIANT(mts1 == mts2, (*dsmi)+std::string(" ")+wsma);
      CHECK_INVARIANT(mV1.size() == mV2.size(), "");

      delete m1;
      delete m2;
    }
    i++;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}



void testIssue196() 
{
  ROMol *mol1=0,*matcher1=0;
  std::string smi,sma;
  
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 196: Smarts handling of 'aa' incorrect" << std::endl;

  smi = "c1ccccc1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  sma = "aa";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);

  MatchVectType mV1;
  bool mts1;
  mts1 = SubstructMatch(*mol1,*matcher1,mV1);
  TEST_ASSERT(mts1);
  TEST_ASSERT(mV1.size()==2);
  
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}


void testIssue254(){
  ROMol *mol1,*mol2;
  ROMol *matcher1,*matcher2;
  std::string smi,sma;
  MatchVectType mV;
  bool mts;


  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 254: Bad handling of unspecified bonds in SMARTS" << std::endl;

  smi ="n1nnnnc1c1nnnnn1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  smi ="c1ccccc1";
  mol2 = SmilesToMol(smi);
  TEST_ASSERT(mol2);

  sma = "cc";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  sma = "[#6][#6]";
  matcher2 = SmartsToMol(sma);
  TEST_ASSERT(matcher2);

  mts=SubstructMatch(*mol1,*matcher1,mV);
  TEST_ASSERT(mts);
  mts=SubstructMatch(*mol2,*matcher1,mV);
  TEST_ASSERT(mts);

  mts=SubstructMatch(*mol1,*matcher2,mV);
  TEST_ASSERT(mts);
  mts=SubstructMatch(*mol2,*matcher2,mV);
  TEST_ASSERT(mts);
    
  delete mol1;
  
  smi ="C1CCCCC1";
  mol1 = SmilesToMol(smi);
  TEST_ASSERT(mol1);
  mts=SubstructMatch(*mol1,*matcher1,mV);
  TEST_ASSERT(!mts);
  mts=SubstructMatch(*mol1,*matcher2,mV);
  TEST_ASSERT(mts);

  delete matcher1;
  sma = "[#6]1[#6][#6][#6][#6][#6]1";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  mts=SubstructMatch(*mol1,*matcher1,mV);
  TEST_ASSERT(mts);
  mts=SubstructMatch(*mol2,*matcher1,mV);
  TEST_ASSERT(mts);

  
  delete mol1;
  delete mol2;
  delete matcher1;
  delete matcher2;
  
  
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue255(){
  ROMol *matcher1;
  std::string sma;
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 255: Core leaks in smarts parsing.  Watch memory consumption." << std::endl;

  for(int i=0;i<10000;i++){
#if 1
    sma = "CC";
    matcher1 = SmartsToMol(sma);
    TEST_ASSERT(matcher1);
    delete matcher1;

    sma = "C-C";
    matcher1 = SmartsToMol(sma);
    TEST_ASSERT(matcher1);
    delete matcher1;

    sma = "C=C";
    matcher1 = SmartsToMol(sma);
    TEST_ASSERT(matcher1);
    delete matcher1;

    sma = "C=,#C";
    matcher1 = SmartsToMol(sma);
    TEST_ASSERT(matcher1);
    delete matcher1;
#endif
    sma = "C-1CC-1";
    //matcher1 = SmartsToMol(sma);
    matcher1 = SmilesToMol(sma);
    TEST_ASSERT(matcher1);
    delete matcher1;
#if 1

    sma = "C-1CC1";
    matcher1 = SmartsToMol(sma);
    TEST_ASSERT(matcher1);
    delete matcher1;
#endif
  }
  
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}



void testIssue330(){
  ROMol *matcher1;
  std::string sma,wsma;
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 330: problems writing some recursive smarts." << std::endl;
  sma = "[$(C=O)]";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  wsma=MolToSmarts(*matcher1);
  BOOST_LOG(rdInfoLog) << "sma: " << wsma << std::endl;

  delete matcher1;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

#endif
void testIssue351(){
  ROMol *matcher1;
  std::string sma;
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 351:" << std::endl;


  sma = "[$(C),S&v2]";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  delete matcher1;

  // this was failing:
  //std::cerr << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << std::endl;
  sma = "[$([C]),S&v2]";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  delete matcher1;
  
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testAtomMap(){
  ROMol *matcher1;
  std::string sma;
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing atom map assignment:" << std::endl;

  int mapNum;
  
  sma = "[C:10]CC";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  TEST_ASSERT(matcher1->getAtomWithIdx(0)->hasProp("molAtomMapNumber"));
  matcher1->getAtomWithIdx(0)->getProp("molAtomMapNumber",mapNum);
  TEST_ASSERT(mapNum==10);
  delete matcher1;

  sma = "[CH3:10]CC";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  TEST_ASSERT(matcher1->getAtomWithIdx(0)->hasProp("molAtomMapNumber"));
  matcher1->getAtomWithIdx(0)->getProp("molAtomMapNumber",mapNum);
  TEST_ASSERT(mapNum==10);
  delete matcher1;

  sma = "[C:10H3]CC";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(!matcher1);

  sma = "[C:10:3]ON";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(!matcher1);

  sma ="C-C";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  sma = MolToSmiles(*matcher1);
  TEST_ASSERT(sma=="CC");

  // test writing the atom map numbers:
  delete matcher1;
  sma="[C:4]";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  sma=MolToSmarts(*matcher1);
  TEST_ASSERT(sma=="[C:4]");

  delete matcher1;
  sma="[C;$(C=O):2]-[O:3]";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  sma=MolToSmarts(*matcher1);
  TEST_ASSERT(sma=="[C&$(C=O):2]-[O:3]");
  
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

#if 1
void testIssue1804420(){
  ROMol *matcher1;
  std::string sma;
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing issue 1804420, missing assignment of atoms maps" << std::endl;

  sma = "[N;D3:1]";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  TEST_ASSERT(matcher1->getAtomWithIdx(0)->hasProp("molAtomMapNumber"));
  delete matcher1;

  sma = "[N,O;D3:1]";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  TEST_ASSERT(matcher1->getAtomWithIdx(0)->hasProp("molAtomMapNumber"));
  delete matcher1;

  sma = "[N&R;X3:1]";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  TEST_ASSERT(matcher1->getAtomWithIdx(0)->hasProp("molAtomMapNumber"));
  delete matcher1;

  sma = "[NH0&R;D3,X3:1]";
  matcher1 = SmartsToMol(sma);
  TEST_ASSERT(matcher1);
  TEST_ASSERT(matcher1->getAtomWithIdx(0)->hasProp("molAtomMapNumber"));
  delete matcher1;

    
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}
#endif

void testSmartsSmiles(){
  RWMol *mol;
  std::string sma,smi;
  
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing cleaner SMARTS -> SMILES " << std::endl;

  smi ="c1ccccc1";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol);
  TEST_ASSERT(smi=="c1ccccc1");

  delete mol;
  smi ="C1CCCCC1";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  smi = MolToSmiles(*mol);
  TEST_ASSERT(smi=="C1CCCCC1");

  
  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testSmilesSmarts(){
  RWMol *mol;
  std::string sma,smi;
  
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing SMILES -> SMARTS" << std::endl;

  smi ="CC";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  sma = MolToSmarts(*mol);
  TEST_ASSERT(sma=="[#6]-[#6]");
  delete mol;

  smi ="C[Si]";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  sma = MolToSmarts(*mol);
  TEST_ASSERT(sma=="[#6]-[Si]");
  delete mol;

  smi ="[CH2-]C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  sma = MolToSmarts(*mol);
  TEST_ASSERT(sma=="[#6H2-]-[#6]");
  delete mol;

  smi ="[CH-2]C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  sma = MolToSmarts(*mol);
  TEST_ASSERT(sma=="[#6H-2]-[#6]");
  delete mol;

  smi ="[CH4+]C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  sma = MolToSmarts(*mol);
  TEST_ASSERT(sma=="[#6H4+]-[#6]");
  delete mol;

  smi ="[CH5+2]C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  sma = MolToSmarts(*mol);
  TEST_ASSERT(sma=="[#6H5+2]-[#6]");
  delete mol;

  smi ="c1ccccc1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  sma = MolToSmarts(*mol);
  TEST_ASSERT(sma=="[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1");
  delete mol;


  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue1914154(){
  RWMol *mol;
  std::string sma;
  
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 1914154: problems with generating smarts for recursive queries" << std::endl;

  sma ="[$(C);$(O)]";
  mol = SmartsToMol(sma);
  TEST_ASSERT(mol);
  sma = MolToSmarts(*mol);
  BOOST_LOG(rdInfoLog)<<sma<<std::endl;
  TEST_ASSERT(sma=="[$(C)&$(O)]");
  delete mol;

  sma ="[$(C),$(O)]";
  mol = SmartsToMol(sma);
  TEST_ASSERT(mol);
  sma = MolToSmarts(*mol);
  BOOST_LOG(rdInfoLog)<<sma<<std::endl;
  TEST_ASSERT(sma=="[$(C),$(O)]");
  delete mol;

  sma ="[!$(C);$(O)]";
  mol = SmartsToMol(sma);
  TEST_ASSERT(mol);
  sma = MolToSmarts(*mol);
  BOOST_LOG(rdInfoLog)<<sma<<std::endl;
  TEST_ASSERT(sma=="[!$(C)&$(O)]");
  delete mol;

  sma ="[C;$(C-O);$(C=O)]";
  mol = SmartsToMol(sma);
  TEST_ASSERT(mol);
  sma = MolToSmarts(*mol);
  BOOST_LOG(rdInfoLog)<<sma<<std::endl;
  TEST_ASSERT(sma=="[C&$(C-O)&$(C=O)]");
  delete mol;

  sma ="[$(C=O),$(C-O);$(C-N)]";
  mol = SmartsToMol(sma);
  TEST_ASSERT(mol);
  sma = MolToSmarts(*mol);
  BOOST_LOG(rdInfoLog)<<sma<<std::endl;
  TEST_ASSERT(sma=="[$(C=O),$(C-O);$(C-N)]");
  delete mol;

  sma ="[$(C-N);$(C=O),$(C-O)]";
  mol = SmartsToMol(sma);
  TEST_ASSERT(mol);
  sma = MolToSmarts(*mol);
  BOOST_LOG(rdInfoLog)<<sma<<std::endl;
  TEST_ASSERT(sma=="[$(C-N);$(C=O),$(C-O)]");
  delete mol;

  sma ="[$(C-N)&$(C=O),$(C-O)]";
  mol = SmartsToMol(sma);
  TEST_ASSERT(mol);
  sma = MolToSmarts(*mol);
  BOOST_LOG(rdInfoLog)<<sma<<std::endl;
  TEST_ASSERT(sma=="[$(C-N)&$(C=O),$(C-O)]");
  delete mol;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testMiscSmartsWriting(){
  RWMol *mol;
  std::string sma;
  
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing miscellaneous bits of SMARTS writing" << std::endl;

  sma ="[13C]";
  mol = SmartsToMol(sma);
  TEST_ASSERT(mol);
  sma = MolToSmarts(*mol);
  TEST_ASSERT(sma=="[C&13*]");
  delete mol;

  sma ="[C]";
  mol = SmartsToMol(sma);
  TEST_ASSERT(mol);
  sma = MolToSmarts(*mol);
  TEST_ASSERT(sma=="C");
  delete mol;

  mol = new RWMol();
  sma = MolToSmarts(*mol);
  TEST_ASSERT(sma=="");
  delete mol;


  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testSmartsStereochem(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing handling (or lack thereof) of stereochem in smarts (sf.net issue 2738320)" << std::endl;

  _checkMatches("C/C=C/C", "CC=CC", 1, 4);
  _checkMatches("C/C=C/C", "C/C=C/C", 1, 4);
  _checkMatches("C/C=C/C", "C\\C=C\\C", 1, 4);
  _checkMatches("C/C=C/C", "C/C=C\\C", 1, 4);
  
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}


void testIssue2884178_part1(){
  
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 2884178 part1: SubstructMatch not returning correct number of matches" << std::endl;

  {
    // part one of the problem: number of unique matches incorrect
    RWMol *patt,*mol;
    std::string sma;
    sma ="*~1~*~*~*~*~*~*~*~*~*~*~*~*~*1";
    patt = SmartsToMol(sma);
    TEST_ASSERT(patt);
    sma ="CC1CCC3C4CCCCC4CCC3C1";
    mol = SmilesToMol(sma);
    TEST_ASSERT(mol);
    std::vector< MatchVectType > mVV;
    unsigned int count=SubstructMatch(*mol,*patt,mVV,true);
    TEST_ASSERT(count==1);
    TEST_ASSERT(mVV.size()==1);
    delete patt;
    delete mol;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue2884178_part2(){
  
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 2884178 part2: SubstructMatch not returning correct number of matches" << std::endl;

  {
    RWMol *patt,*mol;
    std::string sma;
    sma ="C~1~C~C~1";
    patt = SmartsToMol(sma);
    unsigned int count;

    TEST_ASSERT(patt);
    TEST_ASSERT(patt->getBondBetweenAtoms(0,1));
    TEST_ASSERT(patt->getBondBetweenAtoms(0,1)->hasQuery());
    TEST_ASSERT(patt->getBondBetweenAtoms(0,1)->getQuery()->getDescription()=="BondNull");
    TEST_ASSERT(patt->getBondBetweenAtoms(0,2));
    TEST_ASSERT(patt->getBondBetweenAtoms(0,2)->hasQuery());
    TEST_ASSERT(patt->getBondBetweenAtoms(0,2)->getQuery()->getDescription()=="BondNull");

    sma ="C1CC1";
    mol = SmilesToMol(sma);
    TEST_ASSERT(mol);
    std::vector< MatchVectType > mVV;
    count=SubstructMatch(*mol,*patt,mVV,false);
    TEST_ASSERT(count==6);
    TEST_ASSERT(mVV.size()==6);

    delete mol;
    sma ="C1C=C1";
    mol = SmilesToMol(sma);
    TEST_ASSERT(mol);
    count=SubstructMatch(*mol,*patt,mVV,false);
    TEST_ASSERT(count==6);
    TEST_ASSERT(mVV.size()==6);
    
    delete patt;
    delete mol;
  }
  {
    RWMol *patt,*mol;
    std::string sma;
    sma ="C~1~C~C1";
    unsigned int count;

    patt = SmartsToMol(sma);
    TEST_ASSERT(patt);
    TEST_ASSERT(patt->getBondBetweenAtoms(0,1));
    TEST_ASSERT(patt->getBondBetweenAtoms(0,1)->hasQuery());
    TEST_ASSERT(patt->getBondBetweenAtoms(0,1)->getQuery()->getDescription()=="BondNull");
    TEST_ASSERT(patt->getBondBetweenAtoms(0,2));
    TEST_ASSERT(patt->getBondBetweenAtoms(0,2)->hasQuery());
    TEST_ASSERT(patt->getBondBetweenAtoms(0,2)->getQuery()->getDescription()=="BondNull");

    sma ="C1CC1";
    mol = SmilesToMol(sma);
    TEST_ASSERT(mol);
    std::vector< MatchVectType > mVV;
    count=SubstructMatch(*mol,*patt,mVV,false);
    TEST_ASSERT(count==6);
    TEST_ASSERT(mVV.size()==6);

    delete mol;
    sma ="C1C=C1";
    mol = SmilesToMol(sma);
    TEST_ASSERT(mol);
    count=SubstructMatch(*mol,*patt,mVV,false);
    TEST_ASSERT(count==6);
    TEST_ASSERT(mVV.size()==6);
    
    delete patt;
    delete mol;
  }
  {
    RWMol *patt,*mol;
    std::string sma;
    unsigned int count;

    sma ="C1~C~C~1";
    patt = SmartsToMol(sma);
    TEST_ASSERT(patt);
    TEST_ASSERT(patt->getBondBetweenAtoms(0,1));
    TEST_ASSERT(patt->getBondBetweenAtoms(0,1)->hasQuery());
    TEST_ASSERT(patt->getBondBetweenAtoms(0,1)->getQuery()->getDescription()=="BondNull");
    TEST_ASSERT(patt->getBondBetweenAtoms(0,2));
    TEST_ASSERT(patt->getBondBetweenAtoms(0,2)->hasQuery());
    TEST_ASSERT(patt->getBondBetweenAtoms(0,2)->getQuery()->getDescription()=="BondNull");

    sma ="C1CC1";
    mol = SmilesToMol(sma);
    TEST_ASSERT(mol);
    std::vector< MatchVectType > mVV;
    count=SubstructMatch(*mol,*patt,mVV,false);
    TEST_ASSERT(count==6);
    TEST_ASSERT(mVV.size()==6);

    delete mol;
    sma ="C1C=C1";
    mol = SmilesToMol(sma);
    TEST_ASSERT(mol);
    count=SubstructMatch(*mol,*patt,mVV,false);
    TEST_ASSERT(count==6);
    TEST_ASSERT(mVV.size()==6);
    
    delete patt;
    delete mol;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue3000399(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 3000399: incorrect behavior of X queries" << std::endl;

  {
    std::vector< MatchVectType > mVV;
    std::string sma ="[C;X4]";
    RWMol *patt = SmartsToMol(sma);
    TEST_ASSERT(patt);

    sma ="C(N1C)(CCC1)C";
    RWMol *mol = SmilesToMol(sma);
    TEST_ASSERT(mol);
    unsigned int count=SubstructMatch(*mol,*patt,mVV,false);
    TEST_ASSERT(count==6);
    TEST_ASSERT(mVV.size()==6);
    delete mol;

    sma ="[C@H](N1C)(CCC1)C";
    mol = SmilesToMol(sma);
    TEST_ASSERT(mol);
    count=SubstructMatch(*mol,*patt,mVV,false);
    TEST_ASSERT(count==6);
    TEST_ASSERT(mVV.size()==6);
    delete mol;

    delete patt;
  }    

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testRecursiveSerialNumbers(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing serial numbers in recursive SMARTS queries" << std::endl;

  {
    std::vector< MatchVectType > mVV;
    std::string sma ="[$(CO)_1]O[$(CO)_1]";
    RWMol *patt = SmartsToMol(sma);
    TEST_ASSERT(patt);

    sma ="CCOCC";
    RWMol *mol = SmilesToMol(sma);
    TEST_ASSERT(mol);
    unsigned int count=SubstructMatch(*mol,*patt,mVV,false);
    TEST_ASSERT(count==2);
    TEST_ASSERT(mVV.size()==2);
    TEST_ASSERT(mVV[0].size()==3);
    delete mol;

    delete patt;
  }    

  {
    std::vector< MatchVectType > mVV;
    std::string sma ="[$(C(C)O)]O[$(C(C)O)]";
    RWMol *patt = SmartsToMol(sma);
    TEST_ASSERT(patt);

    sma ="CCOCC";
    RWMol *mol = SmilesToMol(sma);
    TEST_ASSERT(mol);
    unsigned int count=SubstructMatch(*mol,*patt,mVV,false);
    TEST_ASSERT(count==2);
    TEST_ASSERT(mVV.size()==2);
    TEST_ASSERT(mVV[0].size()==3);
    delete mol;

    delete patt;
  }    

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}


void testReplacementPatterns(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing use of replacement patterns in input" << std::endl;

  {
    std::vector< MatchVectType > mVV;
    std::string sma ="[{CO}]O[{CO}]";
    std::map<std::string,std::string> repls;
    repls["{CO}"]="$(CO)";
    RWMol *patt = SmartsToMol(sma,0,true,&repls);
    TEST_ASSERT(patt);

    sma ="CCOCC";
    RWMol *mol = SmilesToMol(sma);
    TEST_ASSERT(mol);
    unsigned int count=SubstructMatch(*mol,*patt,mVV,false);
    TEST_ASSERT(count==2);
    TEST_ASSERT(mVV.size()==2);
    TEST_ASSERT(mVV[0].size()==3);
    delete mol;
    delete patt;
  }    

  {
    std::vector< MatchVectType > mVV;
    std::string sma ="[{Q}]O[{Q}]";
    std::map<std::string,std::string> repls;
    repls["{Q}"]="$(C(C)O)";
    RWMol *patt = SmartsToMol(sma,0,true,&repls);
    TEST_ASSERT(patt);

    sma ="CCOCC";
    RWMol *mol = SmilesToMol(sma);
    TEST_ASSERT(mol);
    unsigned int count=SubstructMatch(*mol,*patt,mVV,false);
    TEST_ASSERT(count==2);
    TEST_ASSERT(mVV.size()==2);
    TEST_ASSERT(mVV[0].size()==3);
    delete mol;
    delete patt;
  }    

  {
    std::vector< MatchVectType > mVV;
    std::string sma ="[{MyC}]O";
    std::map<std::string,std::string> repls;
    repls["{CO}"]="$(CO)";
    repls["{C.4}"]="$([C;D4])";
    repls["{MyC}"]="$([{C.4};{CO}])";
    RWMol *patt = SmartsToMol(sma,0,true,&repls);
    TEST_ASSERT(patt);

    sma ="COC(C)(F)(Cl)";
    RWMol *mol = SmilesToMol(sma);
    TEST_ASSERT(mol);
    unsigned int count=SubstructMatch(*mol,*patt,mVV,false);
    TEST_ASSERT(count==1);
    TEST_ASSERT(mVV.size()==1);
    TEST_ASSERT(mVV[0].size()==2);
    delete mol;
    delete patt;
  }    
  
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub313(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github #313: problems with 'h' in SMARTS" << std::endl;

  {
    // basics: does it parse correctly and generate the right results?
    _checkMatches("N[Ch]","FC(Cl)Nc1ccccc1", 1,2);
    _checkMatches("N[Ch1]","CNC(F)c1ccccc1", 1,2);
    _checkMatches("N[Ch]","CNCc1ccccc1", 2,2);
    _checkMatches("N[Ch]","CNc1ccccc1", 1,2);
    _checkNoMatches("N[Ch]","FC(Cl)(O)Nc1ccccc1", false);
  }

  {
    // next: can we write it?
    std::string sma="[h]";
    ROMol *matcher = SmartsToMol(sma);
    TEST_ASSERT(matcher);
    sma=MolToSmarts(*matcher);
    TEST_ASSERT(sma=="[h]");
    delete matcher;

    sma="[h1]";
    matcher = SmartsToMol(sma);
    TEST_ASSERT(matcher);
    sma=MolToSmarts(*matcher);
    TEST_ASSERT(sma=="[h1]");
    delete matcher;    
  }
  
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub314(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github #314: problems with 'x' in SMARTS" << std::endl;

  {
    // basics: does it parse correctly and generate the right results?
    _checkMatches("[x]","C1CC1", 3,1);
    _checkNoMatches("[x]","CCC");
    _checkMatches("[x2]","C1CC1", 3,1);
    _checkNoMatches("[x3]","C1CC1");
    _checkNoMatches("[x3]","CC1CC1");
  }

  {
    // next: can we write it?
    std::string sma="[x]";
    ROMol *matcher = SmartsToMol(sma);
    TEST_ASSERT(matcher);
    sma=MolToSmarts(*matcher);
    TEST_ASSERT(sma=="[x]");
    delete matcher;

    sma="[x1]";
    matcher = SmartsToMol(sma);
    TEST_ASSERT(matcher);
    sma=MolToSmarts(*matcher);
    TEST_ASSERT(sma=="[x1]");
    delete matcher;    
  }
  
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub378(){
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github 378: SMILES parser doing the wrong thing for odd dot-disconnected construct" << std::endl;
  {
    RWMol *m;
    std::string smiles="C1.C1CO1.N1";
    m = SmartsToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(0,1));
    TEST_ASSERT(m->getBondBetweenAtoms(3,4));
    TEST_ASSERT(!m->getBondBetweenAtoms(1,3));
    delete m;
  }
  {
    RWMol *m;
    std::string smiles="C1(O.C1)CO1.N1";
    m = SmartsToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(0,2));
    TEST_ASSERT(m->getBondBetweenAtoms(0,3));
    TEST_ASSERT(m->getBondBetweenAtoms(5,4));
    TEST_ASSERT(!m->getBondBetweenAtoms(2,3));
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}
  

int
main(int argc, char *argv[])
{
  RDLog::InitLogs();
#if 1
  testPass();
  testFail();
  testMatches();
  testMatches2();
  testMatches3();
  testMatches4();
  testMatches5();
  testMatches6();
  testSmartsWrite();
  testFrags();
  testProblems();
  testIssue196();
  testIssue254();
  testIssue255();
  testIssue330();
  testIssue351();
  testAtomMap();
  testSmartsSmiles();
  testSmilesSmarts();
  testIssue1914154();
  testMiscSmartsWriting();
  testIssue1804420();
  testSmartsStereochem();
  testIssue2884178_part1();
  testIssue2884178_part2();
  testIssue3000399();
  testRecursiveSerialNumbers();
  testReplacementPatterns();
#endif
  testGithub313();
  testGithub314();
  testGithub378();
  
  return 0;
}
