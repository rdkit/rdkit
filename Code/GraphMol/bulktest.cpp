//  $Id$
// 
//   Copyright (C) 2002-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/RDLog.h>

#include <iostream>
#include <stdlib.h>

using namespace RDKit;
using namespace std;

void test1(){
  char smi[512],id[256];
  char *fName = "Wrap/test_data/rtecs_smiles.5000.txt";
  FILE *inF = fopen(fName,"r");
  int n = 2;
  Mol *m;
  string canonSmi;
  vector<string> smis;
  n = 2;
  while(n==2){
    n=fscanf(inF,"%s %s",id,smi);
    if(n==2 && id[0]!='#') smis.push_back(smi);
  }
  vector<string>::const_iterator smilIt;
  for(int i=0;i<5;i++){
    BOOST_LOG(rdInfoLog) << "******* Pass " << i+1 << endl;
    for(smilIt=smis.begin();smilIt!=smis.end();smilIt++){
      m = SmilesToMol(*smilIt,false);
      //m = SmilesToMol(*(smis.begin()),false);
      if(!m){
	BOOST_LOG(rdErrorLog) << "failure:" << *smilIt << endl;
      }
#if 0
      CHECK_INVARIANT(m,"Smiles parse failed:");
      canonSmi = MolToSmiles(*m);
      string pickle;
      MolPickler::pickleMol(*m,pickle);
      Mol newM(pickle);
      string trySmi = MolToSmiles(newM);
      if(canonSmi!=trySmi){
	BOOST_LOG(rdErrorLog) << "failure:" << *smilIt << endl;
	BOOST_LOG(rdErrorLog) << canonSmi << endl;
	BOOST_LOG(rdErrorLog) << " != " << endl;
	BOOST_LOG(rdErrorLog) << trySmi << endl;
      }
      CHECK_INVARIANT(canonSmi==trySmi,(char *)(string("pickling/depickling failed: ")+trySmi).c_str());
#endif
      delete m;
    }
  }
};


int main(){
  RDLog::InitLogs()
  test1();

  return 0;
}
