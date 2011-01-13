// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;

std::string getSmiles(std::string line, std::string &smi) {
  int x = line.find("\t");
  if (x == -1) {
    x = line.find(" ");
  }

  smi = line.substr(0,x);
  std::string rem = line.substr(x+1, line.length() - x);
  x = rem.find(" ");
  std::string name;
  name = rem.substr(0, x);
  
  return name;
}

void switchAtoms(RWMol *mol, int idx) {
 
  ROMol::ADJ_ITER nbrIdx,endNbrs;
  boost::tie(nbrIdx,endNbrs) = mol->getAtomNeighbors(mol->getAtomWithIdx(idx));
  INT_VECT neighs;
  std::map<int, Bond*> nbBnd;
  std::map<int, Bond*>::const_iterator nbi;

  // before removing the atom keep track of the bonds and neighboring atoms
  while(nbrIdx != endNbrs) {
    Bond *bnd = new Bond(*(mol->getBondBetweenAtoms(idx, (*nbrIdx))));
    nbBnd[*nbrIdx] = bnd;
    nbrIdx++;
  }
                         
  Atom *oatom = mol->getAtomWithIdx(idx);
  Atom *natom = new Atom(*oatom);
  
  // remove teh aroignal atom
  mol->removeAtom(idx);
  //delete oatom;
  
  // now add it to the end of the molecule
  int nid = mol->addAtom(natom);

  // then reconnect to the old neighbors with the same bonds we stored
  // keep in mind that if a old neighbor index is greater than that of the 
  // removed atom, the id changed by -1
  
  int nbid;
  for (nbi = nbBnd.begin(); nbi != nbBnd.end(); nbi++) {
    nbid = nbi->first;
    if (nbid > idx) {
      nbid--;
    }
    Bond *bnd = new Bond(*(nbi->second));
    bnd->setBeginAtomIdx(nid);
    bnd->setEndAtomIdx(nbid);
    mol->addBond(bnd, true);
  }

  for (nbi = nbBnd.begin(); nbi != nbBnd.end(); nbi++) {
    delete nbi->second;
  }
  nbBnd.clear();
}
    
  
  
void permuteAtoms(RWMol *mol) {
  int nats = mol->getNumAtoms();
  INT_VECT picks;
  INT_VECT_CI pri;
  int i, aid;
  for (i = 0; i < 5; i++) {
    aid = rand()%nats;
    if (std::find(picks.begin(), picks.end(), aid) == picks.end()) {
      picks.push_back(aid);
    }
  }
  
  for (pri = picks.begin(); pri != picks.end(); pri++) {
    switchAtoms(mol, *pri);
  }
  if(mol->hasProp("SubstructGraphPtr")) mol->clearProp("SubstructGraphPtr");
}
  

int main(int argc, char *argv[]) {
  RDLog::InitLogs();
  std::string fname;
  if (argc > 1) {
    fname = argv[1];
  }
  else {
    BOOST_LOG(rdErrorLog) << "Pass in the list of smiles\n";
  }
  
  std::ifstream inStream(fname.c_str());
  const int MAX_LINE_LEN = 512;
  char inLine[MAX_LINE_LEN];
  std::string tmpstr;
  std::string smi;
  inStream.getline(inLine, MAX_LINE_LEN,'\n');
  tmpstr = inLine;
  //MolOps molop;
  int lineCount=0;
  while (tmpstr.size() > 0) {
    lineCount++;
    if(!(lineCount%100)){
      BOOST_LOG(rdErrorLog) << "Doing: " << lineCount << std::endl;
    }
    if(tmpstr[0] != '#'){
      std::string name = getSmiles(tmpstr, smi);
#ifdef VERBOSE_CANON
      BOOST_LOG(rdInfoLog) << "\n-----------------\nSTART: " << smi << std::endl;
#endif      
      //std::cout << " \n";
      RWMol *om = SmilesToMol(smi, 0, 0);
      RWMol *m = SmilesToMol(smi, 0, 0);
      if(m && om){
	try {
	  MolOps::sanitizeMol(*m);
	  MolOps::sanitizeMol(*om);
	  for(int iter=0;iter<5;iter++){
	    std::string preSmi = MolToSmiles(*om);
#if 1
	    permuteAtoms(m);
	    std::string postSmi = MolToSmiles(*m);
#else
	    BOOST_LOG(rdInfoLog) << " ----------------- 2 --------------" << std::endl;
	    std::string postSmi = MolToSmiles(*SmilesToMol(preSmi));
#endif      
	    std::vector<MatchVectType> fgpMatches;
	    int nmat = SubstructMatch(*om, *m, fgpMatches);

	    if (preSmi != postSmi) {
	      BOOST_LOG(rdInfoLog) << lineCount << " " << name << " " << smi << " "
			<< preSmi << " " << postSmi << std::endl;
	      iter =100;
	      std::cout << "\n\n";
	      m->debugMol(std::cout);
	      std::cout << "-*-*-*- MAP -*-*-*-" << std::endl;
	      MatchVectType::const_iterator matchI;
	      for(matchI=fgpMatches[0].begin();
		  matchI!=fgpMatches[0].end();
		  matchI++){
		BOOST_LOG(rdInfoLog) << "\t" << matchI->second << "\t" << matchI->first << std::endl;
	      }
	      break;
	    }
      
	    if (nmat < 1) {
	      BOOST_LOG(rdInfoLog) << lineCount << " Not the same mol: " << name << "\n";
	      exit(-1);
	    }
	  }
	  delete om;
	  delete m;
	}
	catch (MolSanitizeException){
	  BOOST_LOG(rdErrorLog) << smi << "\n";
	  delete om;
	  delete m; 
	}
      } else {
	BOOST_LOG(rdErrorLog) << lineCount << " Parse Failed: " << smi << std::endl;
      }
    }
    inStream.getline(inLine, MAX_LINE_LEN,'\n');
    tmpstr = inLine;
  }
  BOOST_LOG(rdErrorLog) << "FINISHED after " << lineCount << " lines."<< std::endl;
  return 0;

}
      
      
