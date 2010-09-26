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
#include <Catalogs/Catalog.h>
#include "FragCatGenerator.h"
#include "FragCatalogEntry.h"
#include "FragCatParams.h"
#include "FragCatalogUtils.h"

#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace RDKit {
  unsigned int addOrder1Paths(PATH_LIST &paths, const ROMol &mol, 
                              FragCatalog *fcat, DOUBLE_INT_MAP &mapkm1,
                              const MatchVectType &aidToFid) {
    PRECONDITION(fcat,"");
    bool found;
    const FragCatalogEntry *entry;
    //INT_VECT o1entries;
    const FragCatParams *fparams = fcat->getCatalogParams();

    unsigned int lLen = fparams->getLowerFragLength();
    unsigned int uLen = fparams->getUpperFragLength();
    CHECK_INVARIANT(lLen<=uLen,"");
    unsigned int n01 = 0;
    double tol = fparams->getTolerance();

    PATH_LIST_CI pi;
    INT_VECT_CI eti;
    double invar;
    int vid;
    for (pi = paths.begin(); pi != paths.end(); pi++) {
      FragCatalogEntry *nent = new FragCatalogEntry(&mol, (*pi), aidToFid);
      // loop over each order 1 path
      found = false;
      const INT_VECT &o1entries = fcat->getEntriesOfOrder(1);
      for (eti = o1entries.begin(); eti != o1entries.end(); eti++) {
        // loop over all the order 1 entries all ready present in the catalog
        entry = fcat->getEntryWithIdx(*eti);
        if (nent->match(entry, tol)) {
          found = true;
          invar = computeIntVectPrimesProduct(*pi);
          mapkm1[invar] = (*eti);
          delete nent;
          break;
        }
      }
      if (!found) {
        bool updateSigL = false; 
        if ((lLen <= 1) && (uLen >= 1)) {
          // if order 1 subgraphs are of interest in fingerprinting 
          // asign a bit to this fragment when we add to the catalog and update 
          // fingerprint len
          updateSigL = true;
        }
        nent->setDescription(fparams);
        vid = fcat->addEntry(nent, updateSigL);
        invar = computeIntVectPrimesProduct(*pi);
        mapkm1[invar] = vid;
        n01++;
      }
    }
    return n01;
  }

  unsigned int addHigherOrderPaths(const INT_PATH_LIST_MAP &allPaths,
                                   const ROMol &mol, FragCatalog *fcat,
                                   DOUBLE_INT_MAP &mapkm1,
                                   const MatchVectType &aidToFid) {

    PRECONDITION(fcat,"");

    // This works something like this
    // - for each path of order k in the mol
    //    - we find all connected subpaths
    //      of order (k-1) 
    //    - find the entries in the catalog that correspond to each of these
    //      order (k-1) paths (using mapkm1 - remember that this maps the invariant
    //      of a path to the entry ID in the catalog graph)
    //    - Find the intersection of the down entries of these order (k-1) 
    //    - check if order k path we are testing matches any of the order k entries in this intersection
    //    - if we find a match move onto the next order k path
    //    - if we do not find a match 
    //       - create an entry for the order k path and add it to the catalog
    //       - also add out edges from each of the entries corresponding to the order k-1 
    //         subgraphs to this path.

    PATH_LIST paths;
    PATH_LIST_CI pi;
    bool found;
    double invar, sinvar;
    int entId;
    DOUBLE_INT_MAP mapk;
    int mEntId, vid;
    const FragCatParams *fparams = fcat->getCatalogParams();

    unsigned int lLen = fparams->getLowerFragLength();
    unsigned int uLen = fparams->getUpperFragLength();
    double tol = fparams->getTolerance();
    unsigned int nrem = 0; // counter for number of fragments added to the catalog
    const FragCatalogEntry *entry;
    
    INT_PATH_LIST_MAP_CI ordi;
    for (ordi = allPaths.begin(); ordi != allPaths.end(); ordi++) {
      if (ordi->first < 2) {
        continue;
      }
      mapk.clear();
      
      for (pi = (*ordi).second.begin(); pi != (*ordi).second.end(); pi++) {
        found = false;
    
        FragCatalogEntry *nent = new FragCatalogEntry(&mol, (*pi), aidToFid);
        nent->setDescription(fparams);
        
        unsigned int scnt = 0;
        INT_VECT intersect, tmpVect;
        INT_VECT_CI iti;
        invar = computeIntVectPrimesProduct(*pi);
        DOUBLE_VECT sinvarV;
        DOUBLE_VECT_CI sci;
    
        // loop over the subpaths (order (k-1) ) (by ignoring one bond
        // at a time from consideration) and find out which entries int eh catalog they correspond to
        // and make an interestion of the down entries (i.e. order k entries that contain these order k-1
        // entries. - we can baiscally limit our search for an isomorphic entry in the 
        // catalog of the order k path from the molecule to this intersection list
        PATH_TYPE::const_iterator pii;
        for (pii = pi->begin(); pii != pi->end(); pii++) {
          sinvar = invar/firstThousandPrimes[*pii];
    
          // here is a check for "did we see this path before ?" 
          // this should also take care of disconnected subpaths (since the
          // catalog should have only connected subgraphs)
          if (mapkm1.find(sinvar) == mapkm1.end()) {
            continue;
          }
          
          // push this sinvar onto a vector 
          // we need them to add edges int he catalog graph
          sinvarV.push_back(sinvar);
    
          entId = mapkm1[sinvar];
          if (scnt == 0) {
            intersect = fcat->getDownEntryList(entId);
            scnt++;
          } else {
            tmpVect = intersect;
            Intersect(fcat->getDownEntryList(entId), tmpVect, intersect);
            scnt++;
          }
        }
    
        // now search through the intersection list to check if we already have a isomorphic
        // entry in the catalog
        for (iti = intersect.begin(); iti != intersect.end(); iti++) {
          entry = fcat->getEntryWithIdx(*iti);
          if (nent->match(entry, tol) ) {
            found = true;
            mEntId = (*iti);
            delete nent;
            break;
          }
        }
            
        if (found) {
          // update the mapk so that the next time we see this path (when 
          // dealing with order k+1 path we know which entry in the catalog
          // to look at
          mapk[invar] = mEntId;
        } else {
          // we have never seen this subgraph before add it to the catalog
          unsigned int ordr = nent->getOrder();
          bool updateSigL = false;
          if ((ordr >= lLen) && (ordr <= uLen)) {
            // if this order subgraphs are of interest in fingerprinting 
            // asign a bit to this fragment when we add to the catalog and update 
            // fingerprint len
            updateSigL = true;
          }
          
          vid = fcat->addEntry(nent, updateSigL);
          mapk[invar] = vid;
          nrem++; // increment the fragment counter
          // loop over the entries corresponding to the subpaths and 
          // add connections to them
          for (sci = sinvarV.begin(); sci != sinvarV.end(); sci++) {
            entId = mapkm1[*sci];
            fcat->addEdge(entId, vid);
          }
        } // end of never seen this order k subgraph
      } // end of loop over order k paths in mol
      // overwrite mapkm1 with mapk before we move on to order k+1
      mapkm1 = mapk;
    } // end of loop over path order
    return nrem;
  }



  unsigned int FragCatGenerator::addFragsFromMol(const ROMol &mol, FragCatalog *fcat) {
    PRECONDITION(fcat,"");

    INT_PATH_LIST_MAP allPaths;
    allPaths.clear();

    DOUBLE_INT_MAP mapkm1;

    mapkm1.clear();
    const FragCatParams *fparams = fcat->getCatalogParams();
    
    unsigned int lLen = fparams->getLowerFragLength();
    unsigned int uLen = fparams->getUpperFragLength();
    CHECK_INVARIANT(lLen<=uLen,"");

    // prepare the molecule to add to the catalog
    // i.e. find functional groups, remove them from the mol etc.
    MatchVectType newAidToFid;
    INT_VECT fgBonds;
    
    
    ROMol *coreMol = prepareMol(mol, fparams, newAidToFid);
    //mol->debugMol(std::cout);
    allPaths = findAllSubgraphsOfLengthsMtoN(*coreMol, 1, uLen);

    // deal with order 1 paths
    unsigned int nO1Pths = addOrder1Paths(allPaths[1], *coreMol, fcat, mapkm1, newAidToFid);
    
    // now deal with the higher order paths
    unsigned int nremPths = addHigherOrderPaths(allPaths, *coreMol, fcat, mapkm1, newAidToFid);
    
    delete coreMol;
    return (nO1Pths + nremPths);
  }
}
	    
	      
	    
