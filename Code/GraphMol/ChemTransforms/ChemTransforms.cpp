// $Id$
//
//  Copyright (C) 2006 Greg Landrum
//
//   @@ All Rights Reserved  @@
//
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <vector>
#include <algorithm>

namespace RDKit{

  ROMol *deleteSubstructs(const ROMol &mol, const ROMol &query,
                          bool onlyFrags) {
    RWMol *res = static_cast<RWMol*>(new ROMol(mol,false));
    std::vector<MatchVectType> fgpMatches;
    std::vector<MatchVectType>::const_iterator mati;
    std::pair<int, int> amat;
    VECT_INT_VECT matches; // all matches onto the molecule - list of list of atom ids
    MatchVectType::const_iterator mi;
    // do the substructure matching and get the atoms that match the query
    SubstructMatch(*res, query, fgpMatches);

    // if didn't find any matches nothing to be done here
    // simply return a copy of the molecule
    if (fgpMatches.size() == 0) {
      return res;
    }

    for (mati = fgpMatches.begin(); mati != fgpMatches.end(); mati++) {
      INT_VECT match; // each match onto the molecule - list of atoms ids
      for (mi = mati->begin(); mi != mati->end(); mi++) {
        match.push_back(mi->second);
      }
      matches.push_back(match);
    }

    // now loop over the list of matches and check if we can delete any of them
    INT_VECT delList;
    
    VECT_INT_VECT_I mxi, fi;
    if (onlyFrags) {
      VECT_INT_VECT frags;
      
      unsigned int nfrags = MolOps::getMolFrags(*res, frags);
      for (fi = frags.begin(); fi != frags.end(); fi++) {
        std::sort(fi->begin(), fi->end());
        for (mxi = matches.begin(); mxi != matches.end(); mxi++) {
          std::sort(mxi->begin(), mxi->end());
          if ((*fi) == (*mxi) ) {
            INT_VECT tmp; 
            Union((*mxi), delList, tmp);
            delList = tmp;
            break;
          } // end of if we found a matching fragment
        } // endof loop over matches
      } // end of loop over fragments
    } // end of if onlyFrags
    else {
      // in this case we want to delete any matches we find
      // simply loop over the matches and collect the atoms that need to 
      // be removed
      for (mxi = matches.begin(); mxi != matches.end(); mxi++) {
        INT_VECT tmp; 
        Union((*mxi), delList, tmp);
        delList = tmp;
      }
    }

    // now loop over the union list and delete the atoms
    // Will do this in the decreasing order of the atomIds
    // this is so that the AtomIds ids in the "delList" are
    // not invalidated by a previous removal (removing atom number i changes 
    // the atom indices only atoms with indices >i )
    std::sort(delList.begin(), delList.end());

    INT_VECT_RI dri;
    for (dri = delList.rbegin(); dri != delList.rend(); dri++) {
      res->removeAtom(*dri);
    }
    // if we removed any atoms, clear the computed properties:
    if(delList.size()){
      res->clearComputedProps(true);
      // update our properties, but allow unhappiness:
      res->updatePropertyCache(false);
    }
    return res;
  }

  std::vector<ROMOL_SPTR>
  replaceSubstructs(const ROMol &mol, const ROMol &query,const ROMol &replacement,
                    bool replaceAll) {
    std::vector<ROMOL_SPTR> res;
    std::vector<MatchVectType> fgpMatches;

    // do the substructure matching and get the atoms that match the query
    SubstructMatch(mol, query, fgpMatches);

    // if we didn't find any matches, there's nothing to be done here
    // simply return a list with a copy of the starting molecule
    if (fgpMatches.size() == 0) {
      res.push_back(ROMOL_SPTR(new ROMol(mol,false)));
      return res;
    }

    INT_VECT delList;

    // now loop over the list of matches and replace them:
    for (std::vector<MatchVectType>::const_iterator mati = fgpMatches.begin();
         mati != fgpMatches.end(); mati++) {
      INT_VECT match; // each match onto the molecule - list of atoms ids
      for (MatchVectType::const_iterator mi = mati->begin();
           mi != mati->end(); mi++) {
        match.push_back(mi->second);
      }
      INT_VECT sortMatch = match;
      std::sort(sortMatch.begin(),sortMatch.end());

      if( !replaceAll || !res.size() ) {
        res.push_back(ROMOL_SPTR(new ROMol(mol,false)));
      }
      RWMol *newMol = static_cast<RWMol *>(res.rbegin()->get());

      // we need a tab to the orig number of atoms because the
      // new molecule will start numbered above this:
      int numOrigAtoms=newMol->getNumAtoms();
      
      // Add the atoms and bonds from the replacement:
      newMol->insertMol(replacement);

      // loop over the central atom's (the first atom in match) bonds
      // and duplicate any that connect to the remainder of the molecule:
      Atom *origAtom=newMol->getAtomWithIdx(match[0]);
      ROMol::ADJ_ITER nbrIdx,endNbrs;
      boost::tie(nbrIdx,endNbrs) = newMol->getAtomNeighbors(origAtom);
      while(nbrIdx!=endNbrs){
        // we don't want to duplicate any "intra-match" bonds:
        if(!std::binary_search(sortMatch.begin(),sortMatch.end(),int(*nbrIdx))){
          Bond *oBond=newMol->getBondBetweenAtoms(match[0],*nbrIdx);
          CHECK_INVARIANT(oBond,"required bond not found");
          newMol->addBond(numOrigAtoms,*nbrIdx,oBond->getBondType());
        }
        nbrIdx++;
      }
      
      if(replaceAll){
        // we'll accumulate  a list of atoms to be removed:
        INT_VECT tmp; 
        Union(sortMatch, delList, tmp);
        delList = tmp;
      } else {
        // just delete the atoms now:
        for (INT_VECT_RI dri = sortMatch.rbegin(); dri != sortMatch.rend(); dri++) {
          newMol->removeAtom(*dri);
        }
      }
    }

    if(replaceAll){
      // remove the atoms from the delList:
      std::sort(delList.begin(),delList.end());
      RWMol *newMol = static_cast<RWMol *>(res[0].get());
      for (INT_VECT_RI dri = delList.rbegin(); dri != delList.rend(); dri++) {
        newMol->removeAtom(*dri);
      }
    }

    // clear computed props and do basic updates on the
    // the resulting molecules, but allow unhappiness:
    for(std::vector<ROMOL_SPTR>::iterator resI=res.begin();
        resI!=res.end();resI++){
      (*resI)->clearComputedProps(true);
      (*resI)->updatePropertyCache(false);
    }

    return res;
  }

}
