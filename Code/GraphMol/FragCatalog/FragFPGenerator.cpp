// $Id$
//
// Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <Catalogs/Catalog.h>
#include "FragFPGenerator.h"
#include "FragCatalogEntry.h"
#include "FragCatParams.h"
#include "FragCatalogUtils.h"

#include <RDGeneral/types.h>
#include <DataStructs/BitVects.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Subgraphs/Subgraphs.h>

namespace RDKit {
  
  ExplicitBitVect *FragFPGenerator::getFPForMol(const ROMol &mol,
						const FragCatalog &fcat) {
    INT_PATH_LIST_MAP allPaths;

    ExplicitBitVect *fp = new ExplicitBitVect(fcat.getFPLength());
    const FragCatParams *fparams = fcat.getCatalogParams();

    // prepare the molecule for fingerprinting
    // i.e. find functional groups, remove them from the mol etc.
    MatchVectType newAidToFid;
    ROMol *coreMol = prepareMol(mol, fparams, newAidToFid);

    computeFP(*coreMol, fcat, newAidToFid, fp);

    // this was Issue294;
    delete coreMol;
    return fp;
  }

  void FragFPGenerator::computeFP(const ROMol &mol, const FragCatalog &fcat,
				 const MatchVectType &aidToFid, ExplicitBitVect *fp) {
    PRECONDITION(fp, "Bad ExplicitBitVect - FingerPrint");

    const FragCatParams *fparams = fcat.getCatalogParams();

    int uLen = fparams->getUpperFragLength();
    double tol = fparams->getTolerance();

    DOUBLE_INT_MAP mapkm1, mapk;

    // get all the paths in the molecule mapped by their order
    INT_PATH_LIST_MAP allPathsMap = findAllSubgraphsOfLengthsMtoN(mol, 1, uLen);

    // first deal with order 1 stuff
    PATH_LIST_CI pi;
    INT_VECT_CI eti;
    const INT_VECT &o1entries = fcat.getEntriesOfOrder(1);
    const FragCatalogEntry *entry;
    int bitId;
    double invar;
    for (pi = allPathsMap[1].begin(); pi != allPathsMap[1].end(); pi++) {
      //std::cout << "-*-*-* Fragment *-*-*-*-" << std::endl;
      FragCatalogEntry *nent = new FragCatalogEntry(&mol, (*pi), aidToFid);
      nent->setDescription(fparams);
      invar = computeIntVectPrimesProduct(*pi);
      // ok here is the plan - initialize the entry for this path in mapkm1 to -1
      // which will be overwritten to the correct entry id in the catalog if we find 
      // a match. This -1 initialization will be useful when we move onto higher order stuff
      mapkm1[invar] = -1;
      for (eti = o1entries.begin(); eti != o1entries.end(); eti++) {
	entry = fcat.getEntryWithIdx(*eti);
	if (nent->match(entry, tol)) {
	  bitId = entry->getBitId();
	  if (bitId >= 0) {
	    fp->setBit(bitId);
	  }
	  mapkm1[invar] = (*eti);
	  delete nent;
	  break;
	}
      }
    }

    // now deal with the higher order stuff. 
    // - for each higher order path obtain the invariants for order k-1 subpaths
    // - if entries for this invariant values exist in mapkm1, it means that we have
    //   seen this subpath before
    // - if the value for this entry is -1 it mean that we did not find a match for this subpath
    //   in the catalog i.e. we wont find a match for this order k path either
    // - if the entries for all order k-1 subpaths have non-negative values - search the intersection
    //   of the down entries for these catalog entries to find a match for the order k path
    INT_PATH_LIST_MAP_CI ordi;
    double sinvar;
    int entId;

    for (ordi = allPathsMap.begin(); ordi != allPathsMap.end(); ordi++) {
      if (ordi->first < 2) {
	continue; //ignore order 1 paths - we are done with them
      }
      mapk.clear();
      
      for (pi = ordi->second.begin(); pi != ordi->second.end(); pi++) {
	invar = computeIntVectPrimesProduct(*pi);
	mapk[invar] = -1;
	FragCatalogEntry *nent = new FragCatalogEntry(&mol, (*pi), aidToFid);
	nent->setDescription(fparams);
	//std::cout << "Testing 2nd order fragment: " << nent->getDescription() << std::endl;
	int scnt = 0;
	INT_VECT intersect, tmpVect;
	INT_VECT_CI iti;
	PATH_TYPE::const_iterator pii;
	
	// loop over the subpaths (order (k-1) ) (by ignoring one bond
	// at a time from consideration) and find out which entries
	// int eh catalog they correspond to and make an interestion
	// of the down entries (i.e. order k entries that contain
	// these order k-1 entries. - we can baiscally limit our
	// search for an isomorphic entry in the catalog of the order
	// k path from the molecule to this intersection list
	for (pii = (*pi).begin(); pii != (*pi).end(); pii++) {
	  sinvar = invar/firstThousandPrimes[*pii];
	  
	  // here is a check for "did we see this path before ?" 
	  // this should take care of disconnected subpaths (since the
	  // catalog should have only connected subgraphs)
	  if (mapkm1.find(sinvar) == mapkm1.end()) {
	    continue;
	  }
	  
	  if (mapkm1[sinvar] == -1) { 
	    // we have seen this order k-1 pat before but we couldn't find a match
	    // for it in the catalog
	    // which mean that we wont find a match for the order k path either in
	    // the catalog
	    intersect.clear(); // the intersection list we built so far is useless
	    break;
	  }

	  entId = mapkm1[sinvar];
	  if (scnt == 0) {
	    intersect = fcat.getDownEntryList(entId);
	    scnt++;
	  }
	  else {
	    tmpVect = intersect;
	    Intersect(fcat.getDownEntryList(entId), tmpVect, intersect);
	    scnt++;
	  }
	}
	// now search through the intersection list to check if we already have a isomorphic
	// entry in the catalog
	for (iti = intersect.begin(); iti != intersect.end(); iti++) {
	  entry = fcat.getEntryWithIdx(*iti);
	  if (nent->match(entry, tol) ) {
	    mapk[invar] = (*iti);
	    bitId = entry->getBitId();
	    if (bitId >= 0) {
	      fp->setBit(bitId);
	    }
	    delete nent;
	    break;
	  }
	}
      }
      
      // overwrite mapkm1 with mapk before we move on to order k+1
      mapkm1 = mapk;
    } // end of loop over path order
  
	  
#if 0
    /**************************************************************
     * Another way of dealing with higher order paths - didn't work 
     * too well for compound with lots of matches

    PATH_LIST_CI km11, km12, tmpM;
    PATH_TYPE oKpath;
    INT_VECT dEnt1, dEnt2, intersect;
    INT_VECT_CI iti;
    double sinvar1, sinvar2;
    int i, j, eid1, eid2;
    int ord = 2;
    for (ord = 2; ord <= uLen; ord++) {
      std::cout << "Order: " << ord << " " << km1matches.size() << "\n";
      kmatches.clear();
      mapk.clear();
      for (km11 = km1matches.begin(); km11 != km1matches.end(); km11++) {
	tmpM = km11;
	tmpM++;
	
	for (km12 = tmpM;  km12 != km1matches.end(); km12++) {
	  Union((*km11), (*km12), oKpath);
	  if (oKpath.size() == ord) {
	    FragCatalogEntry *nent = new FragCatalogEntry(mol, oKpath, aidToFid);
	    sinvar1 = computeIntVectPrimesProduct(*km11);
	    sinvar2 = computeIntVectPrimesProduct(*km12);
	    
	    eid1 = mapkm1[sinvar1];
	    eid2 = mapkm1[sinvar2];
	    Intersect(fcat.getDownEntryList(eid1),
		      fcat.getDownEntryList(eid2),
		      intersect);
	    // now search through the intersection list to check if we have a isomorphic
	    // entry in the catalog
	    for (iti = intersect.begin(); iti != intersect.end(); iti++) {
	      entry = fcat.getEntryWithIdx(*iti);
	      if (nent->match(entry, tol) ) {
		bitId = entry->getBitId();
		if (bitId >= 0) {
		  fp->setBit(bitId);
		}
		invar = computeIntVectPrimesProduct(oKpath);
		kmatches.push_back(oKpath);
		mapk[invar] = (*iti);
	      }
	    }// end of for loop over intersect entries
	  } // end of if block (if we found a order path)
	}
      }
      km1matches = kmatches;
      mapkm1 = mapk;
    } // end of loop over higher order paths order
    **************************/
#endif    
  }
}
	    
      
      
    
    
