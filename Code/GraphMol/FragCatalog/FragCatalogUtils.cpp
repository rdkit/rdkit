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

// REVIEW: move this to a GraphMol/FuncGroups directory

#include <RDGeneral/BadFileException.h>
#include "FragCatalogUtils.h"
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <fstream>
#include <string>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
#include <boost/algorithm/string.hpp>

namespace RDKit {
  // local utility namespace
  namespace {
    ROMol * getSmarts(const std::string &tmpStr) {
      ROMol *mol = 0;
      if (tmpStr.length() == 0) {
	// empty line
	return mol;
      }
      if (tmpStr.substr(0,2) == "//") {
	// comment line
	return mol;
      }

      boost::char_separator<char> tabSep("\t");
      tokenizer tokens(tmpStr,tabSep);
      tokenizer::iterator token=tokens.begin();
    
      // name of the functional groups
      std::string name = *token;
      boost::erase_all(name," ");
      ++token;

      // grab the smarts:
      std::string smarts = *token;
      boost::erase_all(smarts," ");
      ++token;

      mol = SmartsToMol(smarts);
      CHECK_INVARIANT(mol,smarts);
      mol->setProp("_Name", name);
      mol->setProp("_fragSMARTS",smarts);
      return mol;
    }    
  } // end of local utility namespace
    
  MOL_SPTR_VECT readFuncGroups(std::istream &inStream,int nToRead) {
    MOL_SPTR_VECT funcGroups;
    funcGroups.clear();
    if (inStream.bad()) {
      throw BadFileException("Bad stream contents.");
    }
    const int MAX_LINE_LEN = 512;
    char inLine[MAX_LINE_LEN];
    std::string tmpstr;
    int nRead=0;
    while (!inStream.eof() && (nToRead<0 || nRead<nToRead)) {
      inStream.getline(inLine, MAX_LINE_LEN,'\n');
      tmpstr = inLine;
      // parse the molecule on this line (if there is one)
      ROMol *mol = getSmarts(tmpstr);
      if (mol) {
	funcGroups.push_back(ROMOL_SPTR(mol));
	nRead++;
      }
    }
    return funcGroups;
  }

  MOL_SPTR_VECT readFuncGroups(std::string fileName) {
    std::ifstream inStream(fileName.c_str());
    if ((!inStream) || (inStream.bad()) ) {
      std::ostringstream errout;
      errout << "Bad input file " << fileName;
      throw BadFileException(errout.str());
    }
    MOL_SPTR_VECT funcGroups;
    funcGroups=readFuncGroups(inStream);
    return funcGroups;
  }

  MatchVectType findFuncGroupsOnMol(const ROMol &mol, 
				    const FragCatParams *params,
				    INT_VECT &fgBonds) {
    PRECONDITION(params,"bad params");

    fgBonds.clear();
    
    std::pair<int, int> amat;
    MatchVectType aidFgrps;
    std::vector<MatchVectType> fgpMatches;
    std::vector<MatchVectType>::const_iterator mati;
    MatchVectType::const_iterator mi;
    int aid;
    //const ROMol *fgrp;

    INT_VECT_CI bi;
    aidFgrps.clear();
    
    int fid = 0;
    const MOL_SPTR_VECT &fgrps = params->getFuncGroups();
    MOL_SPTR_VECT::const_iterator fgci;
    
    for (fgci = fgrps.begin(); fgci != fgrps.end(); fgci++) {
      const ROMol *fgrp = fgci->get();
      std::string fname;
      (*fgci)->getProp("_Name", fname);
      //std::cout << "Groups number: " << fname << "\n";
      //(*fgci)->debugMol(std::cout);
      //mol->debugMol(std::cout);
      // match this functional group onto the molecule
      SubstructMatch(mol, *fgrp, fgpMatches);

      // loop over all the matches we get for this fgroup
      for (mati = fgpMatches.begin(); mati != fgpMatches.end(); mati++) {
	//FIX: we will assume that the first atom in fgrp is always the connection
	// atom
	amat = mati->front();
	aid = amat.second; //FIX: is this correct - the second entry in the pair is the atom ID from mol

	// grab the list of atom Ids from mol that match the functional group 
	INT_VECT bondIds, maids;
	for (mi = mati->begin(); mi != mati->end(); mi++) {
	  maids.push_back(mi->second);
	}

	// create a list of bond IDs from these atom ID 
	// these are the bond in mol that are part of portion that matches the 
	// functional group
	bondIds = Subgraphs::bondListFromAtomList(mol, maids);
	
	// now check if all these bonds have been covered as part of larger 
	// functional group that was dealt with earlier
	// FIX: obviously we assume here that the function groups in params 
	// come in decreasing order of size.
	bool allDone = true;
	for (bi = bondIds.begin(); bi != bondIds.end(); bi++) {
	  if (std::find(fgBonds.begin(), fgBonds.end(), (*bi)) == fgBonds.end()) {
	    allDone = false;
	    fgBonds.push_back(*bi);
	  }
	}
	
	if (!allDone) {
	  // this functional group mapping onto mol is not part of a larger func
	  // group mapping so record it
	  aidFgrps.push_back(std::pair<int, int>(aid, fid));
	}
      }
      fid++;

    }

    
    return aidFgrps;
  }

  ROMol *prepareMol(const ROMol &mol, const FragCatParams *fparams,
		    MatchVectType &aToFmap) {
    PRECONDITION(fparams,"");
    
    // get a mapping of the functional groups onto the molecule
    INT_VECT fgBonds;
    MatchVectType aidToFid = findFuncGroupsOnMol(mol, fparams, fgBonds);
    
    // get the core piece of molecule (i.e. the part of the molecule 
    // without the functional groups). This basically the part of the molecule
    // that does not contain the function group bonds given by "fgBonds"
    INT_VECT cBonds;
    int bid, nbds = mol.getNumBonds();

    for (bid = 0; bid < nbds; bid++) {
      if (std::find(fgBonds.begin(), fgBonds.end(), bid) == fgBonds.end()) {
	cBonds.push_back(bid);
      }
    }

    INT_MAP_INT aIdxMap; // a map from atom id in mol to the new atoms id in coreMol
    
    ROMol *coreMol = Subgraphs::pathToSubmol(mol, cBonds, false, aIdxMap);
    
    // now map the functional groups on mol to coreMol using aIdxMap
    MatchVectType::iterator mati;
    
    int newID;
    for (mati = aidToFid.begin(); mati != aidToFid.end(); mati++) {
      newID = aIdxMap[mati->first];
      aToFmap.push_back(std::pair<int, int>(newID, mati->second));
    }
    
    return coreMol;
  }
}

