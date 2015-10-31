//
//  Copyright (c) 2015 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/GraphMol.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>

#include <vector>
#include <algorithm> 

namespace RDKit {
  namespace MolOps {
    ROMol *adjustQueryProperties(const ROMol &mol,const AdjustQueryParameters *params){
      RWMol *res = new RWMol(mol);
      try{
        adjustQueryProperties(*res,params);
      } catch (MolSanitizeException &se){
        delete res;
        throw se;
      }
      return static_cast<ROMol *>(res);
    }
    void adjustQueryProperties(RWMol &mol,const AdjustQueryParameters *inParams){
      AdjustQueryParameters params;
      if(inParams){
        params = *inParams;
      }
      const RingInfo *ringInfo=mol.getRingInfo();
      if( !ringInfo->isInitialized() ){
        MolOps::symmetrizeSSSR(mol);
      }

      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        Atom *at = mol.getAtomWithIdx(i);
        // pull properties we need from the atom here, once we
        // create a query atom they may no longer be valid.
        unsigned int nRings = ringInfo->numAtomRings(i);
        int atomicNum = at->getAtomicNum();
        if( params.adjustDegree && 
            !( (params.adjustDegreeFlags & ADJUST_RINGSONLY) && !nRings ) && 
            !( (params.adjustDegreeFlags & ADJUST_IGNOREDUMMIES) && !atomicNum) 
            ){
          QueryAtom *qa;
          if(!at->hasQuery()){
            qa = new QueryAtom(*at);
            mol.replaceAtom(i,qa);
            delete qa;
            qa = static_cast<QueryAtom *>(mol.getAtomWithIdx(i));
            at = static_cast<Atom *>(qa);
          } else {
            qa = static_cast<QueryAtom *>(at);
          }
          qa->expandQuery(makeAtomExplicitDegreeQuery(qa->getDegree()));
        } // end of adjust degree
        if( params.adjustRingCount && 
            !( (params.adjustRingCountFlags & ADJUST_RINGSONLY) && !nRings ) && 
            !( (params.adjustRingCountFlags & ADJUST_IGNOREDUMMIES) && !atomicNum )
            ){
          QueryAtom *qa;
          if(!at->hasQuery()){
            qa = new QueryAtom(*at);
            mol.replaceAtom(i,qa);
            delete qa;
            qa = static_cast<QueryAtom *>(mol.getAtomWithIdx(i));
            at = static_cast<Atom *>(qa);
          } else {
            qa = static_cast<QueryAtom *>(at);
          }
          qa->expandQuery(makeAtomInNRingsQuery(nRings));
        } // end of adjust ring count
        if(params.makeDummiesQueries && 
           atomicNum == 0 && !at->hasQuery() &&
           !at->getIsotope() ){
            QueryAtom *qa = new QueryAtom();
            qa->setQuery(makeAtomNullQuery());
            mol.replaceAtom(i,qa);
            delete qa;
            at = mol.getAtomWithIdx(i);
        } // end of makeDummiesQueries
      } // end of loop over atoms
    }
  } // end of MolOps namespace
} // end of RDKit namespace
