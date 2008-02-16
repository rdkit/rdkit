// $Id$
//
//  Copyright (C) 2006-22007 Greg Landrum
//
//   @@ All Rights Reserved  @@
//
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/dynamic_bitset.hpp>
#include <vector>
#include <algorithm>

namespace RDKit{

  ROMol *deleteSubstructs(const ROMol &mol, const ROMol &query,
                          bool onlyFrags) {
    RWMol *res = static_cast<RWMol*>(new ROMol(mol,false));
    std::vector<MatchVectType> fgpMatches;
    std::vector<MatchVectType>::const_iterator mati;
    std::pair<int, int> amat;
    VECT_INT_VECT matches; // all matches on the molecule - list of list of atom ids
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
      
      MolOps::getMolFrags(*res, frags);
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
      res[0]->clearComputedProps(false);
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

  ROMol *replaceSidechains(const ROMol &mol, const ROMol &coreQuery){
    MatchVectType matchV;

    // do the substructure matching and get the atoms that match the query
    bool matchFound=SubstructMatch(mol, coreQuery, matchV);

    // if we didn't find any matches, there's nothing to be done here
    // simply return null to indicate the problem
    if (!matchFound || matchV.size()==0){
      return 0;
    }

    boost::dynamic_bitset<> matchingIndices(mol.getNumAtoms());
    for(MatchVectType::const_iterator mvit=matchV.begin();
        mvit!=matchV.end();mvit++){
      matchingIndices[mvit->second] = 1;
    }
    std::vector<Atom *> keepList;

    RWMol *newMol = static_cast<RWMol *>(new ROMol(mol));
    unsigned int nDummies=0;
    for(MatchVectType::const_iterator mvit=matchV.begin();
        mvit!=matchV.end();mvit++){
      keepList.push_back(newMol->getAtomWithIdx(mvit->second));
      // if the atom in the molecule has higher degree than the atom in the
      // core, we have an attachment point:
      if( newMol->getAtomWithIdx(mvit->second)->getDegree() > 
          coreQuery.getAtomWithIdx(mvit->first)->getDegree()){ 
        ROMol::ADJ_ITER nbrIdx,endNbrs;
        boost::tie(nbrIdx,endNbrs) = newMol->getAtomNeighbors(newMol->getAtomWithIdx(mvit->second));
        while(nbrIdx!=endNbrs){
          if(!matchingIndices[*nbrIdx]){
            // this neighbor isn't in the match, convert it to a dummy atom and save it
            Atom *at=newMol->getAtomWithIdx(*nbrIdx);
            at->setAtomicNum(0);
            std::string label="X";
            label+=static_cast<char>(static_cast<short>('a')+nDummies);
            at->setProp("dummyLabel",label);
            keepList.push_back(at);
            ++nDummies;
          }
          nbrIdx++;
        }  
      }
    }
    std::vector<Atom *> delList;
    for(RWMol::AtomIterator atIt=newMol->beginAtoms();atIt!=newMol->endAtoms();atIt++){
      Atom *tmp=*atIt;
      if(std::find(keepList.begin(),keepList.end(),tmp)==keepList.end()){
        delList.push_back(tmp);
      } 
    }    
    for(std::vector<Atom *>::const_iterator delIt=delList.begin();delIt!=delList.end();delIt++){
      newMol->removeAtom(*delIt);
    }

    // clear computed props and do basic updates on the
    // the resulting molecule, but allow unhappiness:
    newMol->clearComputedProps(true);
    newMol->updatePropertyCache(false);
    return static_cast<ROMol *>(newMol);
  }

  ROMol *replaceCore(const ROMol &mol, const ROMol &coreQuery, bool replaceDummies){
    MatchVectType matchV;

    // do the substructure matching and get the atoms that match the query
    bool matchFound=SubstructMatch(mol, coreQuery, matchV);

    // if we didn't find any matches, there's nothing to be done here
    // simply return null to indicate the problem
    if (!matchFound || matchV.size()==0){
      return 0;
    }

    unsigned int origNumAtoms=mol.getNumAtoms();
    boost::dynamic_bitset<> matchingIndices(origNumAtoms);
    for(MatchVectType::const_iterator mvit=matchV.begin();
        mvit!=matchV.end();mvit++){
      if(replaceDummies || coreQuery.getAtomWithIdx(mvit->first)->getAtomicNum()>0){
        matchingIndices[mvit->second] = 1;
      }
    }

    RWMol *newMol = static_cast<RWMol *>(new ROMol(mol));
    std::vector<Atom *> keepList;
    unsigned int nDummies=0;
    for(unsigned int i=0;i<origNumAtoms;++i){
      if(!matchingIndices[i]){
        Atom *sidechainAtom=newMol->getAtomWithIdx(i);
        // we're keeping the sidechain atoms:
        keepList.push_back(sidechainAtom);

        // loop over our neighbors and see if any are in the match:
        std::list<unsigned int> nbrList;
        ROMol::ADJ_ITER nbrIter,endNbrs;
        boost::tie(nbrIter,endNbrs) = newMol->getAtomNeighbors(sidechainAtom);
        while(nbrIter!=endNbrs && (*nbrIter) < origNumAtoms){
          // we need to add bonds and atoms to the molecule while looping
          // over neighbors. This invalidates iterators, so collect a list
          // of our neighbors now:
          nbrList.push_back(*nbrIter);
          ++nbrIter;
        }
        unsigned int whichNbr=0;
        std::list<Bond *> newBonds;
        for(std::list<unsigned int>::const_iterator lIter=nbrList.begin();
            lIter!=nbrList.end();++lIter){
          unsigned int nbrIdx=*lIter;
          Bond *connectingBond=newMol->getBondBetweenAtoms(i,nbrIdx);
          if(matchingIndices[nbrIdx]){
            bool removedPrecedingAtom=false;
            Atom *newAt=new Atom(0);
            std::string label="X";
            label+=static_cast<char>(static_cast<short>('a')+nDummies);
            ++nDummies;
            newAt->setProp("dummyLabel",label);
            newMol->addAtom(newAt,false,true);
            keepList.push_back(newAt);
            Bond *bnd=connectingBond->copy();
            if(bnd->getBeginAtomIdx()==i){
              bnd->setEndAtomIdx(newAt->getIdx());
              removedPrecedingAtom=false;
            } else {
              bnd->setBeginAtomIdx(newAt->getIdx());
              if(nbrIdx<i) removedPrecedingAtom=true;
            }
            newBonds.push_back(bnd);

            // we may be changing the bond ordering at the atom.
            // e.g. replacing the N in C[C@](Cl)(N)F gives an atom ordering of C[C?](Cl)(F)[X]
            // so we need the SMILES C[C@@](Cl)(F)[X] to maintain the appropriate chirality 
            // check for these cases and adjust our chirality flags as appropriate.
            //
            if(sidechainAtom->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW
               || sidechainAtom->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW ){
              bool switchIt=false;
              switch(newMol->getAtomDegree(sidechainAtom)){
              case 4:
                //     start:         ordering:        swap?
                //   N[C@](F)(Cl)C -> F[C@@](Cl)(C)X    yes 
                //   F[C@](N)(Cl)C -> F[C@](Cl)(C)X     no
                //   F[C@](Cl)(N)C -> F[C@@](Cl)(C)X    yes
                //   F[C@](Cl)(C)N -> F[C@](Cl)(C)X     no
                if(!(whichNbr%2)) switchIt=true;
                break;  
              case 3:
                // things are different in the degree three case because of the implicit H:
                //     start:         ordering:     swap?
                //   N[C@H](F)C -> [C@H](F)(C)X     yes
                //   [C@H](N)(F)C -> [C@H](F)(C)X   no
                //   F[C@H](N)C -> F[C@@H](C)X      yes        
                //   F[C@H](C)N -> F[C@H](C)X       no
                if(whichNbr==1 || (whichNbr==0&&removedPrecedingAtom) ) switchIt=true;
                break;  
              }
              if(switchIt){
                sidechainAtom->invertChirality();
              }
            }
          }
          ++whichNbr;
        }
        // add the bonds now, after we've finished the loop over neighbors:
        for(std::list<Bond *>::iterator bi=newBonds.begin();bi!=newBonds.end();++bi){
          newMol->addBond(*bi,true);
        }
      }
    }
    std::vector<Atom *> delList;
    for(RWMol::AtomIterator atIt=newMol->beginAtoms();atIt!=newMol->endAtoms();atIt++){
      Atom *tmp=*atIt;
      if(std::find(keepList.begin(),keepList.end(),tmp)==keepList.end()){
        delList.push_back(tmp);
      } 
    }    
    for(std::vector<Atom *>::const_iterator delIt=delList.begin();delIt!=delList.end();delIt++){
      newMol->removeAtom(*delIt);
    }

    // clear computed props and do basic updates on
    // the resulting molecule, but allow unhappiness:
    newMol->clearComputedProps(true);
    newMol->updatePropertyCache(false);

    return static_cast<ROMol *>(newMol);

    
  }


}  // end of namespace RDKit
