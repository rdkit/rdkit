// $Id$
//
//  Copyright (C) 2006-2012 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitQueries.h>
#include <RDBoost/Exceptions.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>

namespace RDKit{
  namespace {
    void updateSubMolConfs(const ROMol &mol,RWMol &res,
                           boost::dynamic_bitset<> &removedAtoms){
      // update conformer information:
      res.clearConformers();
      for(ROMol::ConstConformerIterator citer=mol.beginConformers();
          citer!=mol.endConformers();++citer){
        Conformer *newConf=new Conformer(res.getNumAtoms());
        newConf->setId((*citer)->getId());
        int aIdx=0;
        for(unsigned int i=0;i<mol.getNumAtoms();++i){
          if(!removedAtoms[i]){
            newConf->setAtomPos(aIdx,(*citer)->getAtomPos(i));
            ++aIdx;
          }
        }
        res.addConformer(newConf,false);
      }
    }

  }
  
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

    boost::dynamic_bitset<> removedAtoms(mol.getNumAtoms());
    for (INT_VECT_RI dri = delList.rbegin(); dri != delList.rend(); dri++) {
      removedAtoms.set(*dri);
      res->removeAtom(*dri);
    }
    // if we removed any atoms, clear the computed properties:
    if(delList.size()){
      updateSubMolConfs(mol,*res,removedAtoms);

      res->clearComputedProps(true);
      // update our properties, but allow unhappiness:
      res->updatePropertyCache(false);
    }
    return res;
  }

  std::vector<ROMOL_SPTR>
  replaceSubstructs(const ROMol &mol, const ROMol &query,const ROMol &replacement,
                    bool replaceAll,unsigned int replacementConnectionPoint) {
    PRECONDITION(replacementConnectionPoint<replacement.getNumAtoms(),"bad replacementConnectionPoint");
    std::vector<ROMOL_SPTR> res;
    std::vector<MatchVectType> fgpMatches;

    boost::dynamic_bitset<> removedAtoms(mol.getNumAtoms());

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
          newMol->addBond(numOrigAtoms+replacementConnectionPoint,
                          *nbrIdx,oBond->getBondType());
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
          removedAtoms.set(*dri);
          newMol->removeAtom(*dri);
        }
      }
    }

    if(delList.size()){
      if(replaceAll){
        // remove the atoms from the delList:
        std::sort(delList.begin(),delList.end());
        RWMol *newMol = static_cast<RWMol *>(res[0].get());
        for (INT_VECT_RI dri = delList.rbegin(); dri != delList.rend(); dri++) {
          removedAtoms.set(*dri);
          newMol->removeAtom(*dri);
        }
      }
        
      // clear conformers and computed props and do basic updates 
      // on the the resulting molecules, but allow unhappiness:
      for(std::vector<ROMOL_SPTR>::iterator resI=res.begin();
          resI!=res.end();resI++){
        updateSubMolConfs(mol,*(RWMol *)(*resI).get(),removedAtoms);
        (*resI)->clearComputedProps(true);
        (*resI)->updatePropertyCache(false);
      }
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
            ++nDummies;
            at->setIsotope(nDummies);
            keepList.push_back(at);
          }
          nbrIdx++;
        }  
      }
    }
    std::vector<Atom *> delList;
    boost::dynamic_bitset<> removedAtoms(mol.getNumAtoms());
    for(RWMol::AtomIterator atIt=newMol->beginAtoms();atIt!=newMol->endAtoms();atIt++){
      Atom *tmp=*atIt;
      if(std::find(keepList.begin(),keepList.end(),tmp)==keepList.end()){
        delList.push_back(tmp);
        removedAtoms.set(tmp->getIdx());
      } 
    }    
    for(std::vector<Atom *>::const_iterator delIt=delList.begin();delIt!=delList.end();delIt++){
      newMol->removeAtom(*delIt);
    }

    updateSubMolConfs(mol,*newMol,removedAtoms);
    
    // clear computed props and do basic updates on the
    // the resulting molecule, but allow unhappiness:
    newMol->clearComputedProps(true);
    newMol->updatePropertyCache(false);
    return static_cast<ROMol *>(newMol);
  }

  ROMol *replaceCore(const ROMol &mol, const ROMol &coreQuery, bool replaceDummies,
                     bool labelByIndex,bool requireDummyMatch){
    MatchVectType matchV;

    // do the substructure matching and get the atoms that match the query
    bool matchFound=SubstructMatch(mol, coreQuery, matchV);
    
    // if we didn't find any matches, there's nothing to be done here
    // simply return null to indicate the problem
    if (!matchFound || matchV.size()==0){
      return 0;
    }

    unsigned int origNumAtoms=mol.getNumAtoms();
    std::vector<int> matchingIndices(origNumAtoms,-1);
    for(MatchVectType::const_iterator mvit=matchV.begin();
        mvit!=matchV.end();mvit++){
      if(replaceDummies || coreQuery.getAtomWithIdx(mvit->first)->getAtomicNum()>0){
        matchingIndices[mvit->second] = mvit->first;
      }
    }

    RWMol *newMol = static_cast<RWMol *>(new ROMol(mol));
    std::vector<Atom *> keepList;
    std::map< int,Atom *> dummyAtomMap;
    unsigned int nDummies=0;
    for(unsigned int i=0;i<origNumAtoms;++i){
      if(matchingIndices[i]==-1){
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
          if(matchingIndices[nbrIdx]>-1){
            // we've matched an atom in the core.
            if(requireDummyMatch &&
               coreQuery.getAtomWithIdx(matchingIndices[nbrIdx])->getAtomicNum()!=0){
              delete newMol;
              return NULL;
            }
            Atom *newAt=new Atom(0);
            ++nDummies;
            
            if(!labelByIndex){
              newAt->setIsotope(nDummies);
            } else {
              newAt->setIsotope(matchingIndices[nbrIdx]);
            }
            newMol->addAtom(newAt,false,true);
            dummyAtomMap[nbrIdx]=newAt;
            keepList.push_back(newAt);
            Bond *bnd=connectingBond->copy();
            if(bnd->getBeginAtomIdx()==i){
              bnd->setEndAtomIdx(newAt->getIdx());
            } else {
              bnd->setBeginAtomIdx(newAt->getIdx());
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
                if(whichNbr==1 ) switchIt=true;
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
    boost::dynamic_bitset<> removedAtoms(mol.getNumAtoms());
    for(RWMol::AtomIterator atIt=newMol->beginAtoms();atIt!=newMol->endAtoms();atIt++){
      Atom *tmp=*atIt;
      if(std::find(keepList.begin(),keepList.end(),tmp)==keepList.end()){
        delList.push_back(tmp);
        removedAtoms.set(tmp->getIdx());
      } 
    }    
    for(std::vector<Atom *>::const_iterator delIt=delList.begin();delIt!=delList.end();delIt++){
      newMol->removeAtom(*delIt);
    }

    updateSubMolConfs(mol,*newMol,removedAtoms);

    // make a guess at the position of the dummy atoms showing the attachment point:
    for(ROMol::ConstConformerIterator citer=mol.beginConformers();
        citer!=mol.endConformers();++citer){
      Conformer &newConf=newMol->getConformer((*citer)->getId());
      for(std::map<int,Atom *>::const_iterator iter=dummyAtomMap.begin();iter!=dummyAtomMap.end();++iter){
        newConf.setAtomPos(iter->second->getIdx(),(*citer)->getAtomPos(iter->first));
      }
    }
    // clear computed props and do basic updates on
    // the resulting molecule, but allow unhappiness:
    newMol->clearComputedProps(true);
    newMol->updatePropertyCache(false);

    return static_cast<ROMol *>(newMol);

    
  }

  ROMol *MurckoDecompose(const ROMol &mol){
    RWMol *res=new RWMol(mol);
    unsigned int nAtoms=res->getNumAtoms();
    if(!nAtoms) return res;

    // start by getting the shortest paths matrix:
    MolOps::getDistanceMat(mol,false,false,true);
    boost::shared_array<int> pathMat;
    mol.getProp("DistanceMatrix_Paths",pathMat);

    boost::dynamic_bitset<> keepAtoms(nAtoms);
    const RingInfo *ringInfo=res->getRingInfo();
    for(unsigned int i=0;i<nAtoms;++i){
      if(ringInfo->numAtomRings(i)) keepAtoms[i]=1;
    }
    const VECT_INT_VECT &rings=ringInfo->atomRings();
    //std::cerr<<"  rings: "<<rings.size()<<std::endl;
    // now find the shortest paths between each ring system and mark the atoms
    // along each as being keepers:
    for(VECT_INT_VECT::const_iterator ringsItI=rings.begin();
        ringsItI!=rings.end();++ringsItI){
      for(VECT_INT_VECT::const_iterator ringsItJ=ringsItI+1;
          ringsItJ!=rings.end();++ringsItJ){
        int atomI=(*ringsItI)[0];
        int atomJ=(*ringsItJ)[0];
        //std::cerr<<atomI<<" -> "<<atomJ<<": ";
        while(atomI!=atomJ){
          keepAtoms[atomI]=1;
          atomI = pathMat[atomJ*nAtoms+atomI];
          // test for the disconnected case:
          if(atomI<0) break;
          //std::cerr<<atomI<<" ";
        }
        //std::cerr<<std::endl;
      }
    }

    boost::dynamic_bitset<> removedAtoms(nAtoms);
    std::vector<Atom *> atomsToRemove;
    for(unsigned int i=0;i<nAtoms;++i){
      if(!keepAtoms[i]){
        Atom *atom=res->getAtomWithIdx(i);
        bool removeIt=true;

        // check if the atom has a neighboring keeper:
        ROMol::ADJ_ITER nbrIdx,endNbrs;
        boost::tie(nbrIdx,endNbrs) = res->getAtomNeighbors(atom);
        while(nbrIdx!=endNbrs){
          const ATOM_SPTR nbr=(*res)[*nbrIdx];
          if(keepAtoms[nbr->getIdx()]){
            if(res->getBondBetweenAtoms(atom->getIdx(),nbr->getIdx())->getBondType()==Bond::DOUBLE){
              removeIt=false;
              break;
            } else if(nbr->getIsAromatic() && nbr->getAtomicNum()!=6 ){
              // fix aromatic heteroatoms:
              nbr->setNumExplicitHs(1);
            } else if(nbr->getNoImplicit() || nbr->getChiralTag()!=Atom::CHI_UNSPECIFIED){
              nbr->setNoImplicit(false);
              nbr->setNumExplicitHs(0);
              nbr->setChiralTag(Atom::CHI_UNSPECIFIED);
            }
          }
          ++nbrIdx;
        }
        
        if(removeIt){
          atomsToRemove.push_back(atom);
          removedAtoms.set(atom->getIdx());
        }
      }
    }

    for(std::vector<Atom *>::const_iterator atomIt=atomsToRemove.begin();
        atomIt!=atomsToRemove.end();++atomIt){
      res->removeAtom(*atomIt);
    }
    updateSubMolConfs(mol,*res,removedAtoms);
    res->clearComputedProps();
    
    return (ROMol *)res;
  }

  ROMol *combineMols(const ROMol &mol1,const ROMol &mol2,
                     RDGeom::Point3D offset){
    RWMol *res=new RWMol(mol1);
    int nAtoms1=res->getNumAtoms();
    res->insertMol(mol2);

    // copy over coordinates
    if(mol1.getNumConformers() && mol2.getNumConformers()){
      if(mol1.getNumConformers() != mol2.getNumConformers()){
        BOOST_LOG(rdWarningLog)<<"combineMols: molecules have unequal numbers of conformers"<<std::endl;
      }
      for(ROMol::ConformerIterator conf1It=res->beginConformers();
          conf1It!=res->endConformers();++conf1It){
        Conformer *conf1=(*conf1It).get();
        try{
          const Conformer *conf2=&mol2.getConformer(conf1->getId());
          for(unsigned int i=0;i<mol2.getNumAtoms();++i){
            conf1->setAtomPos(i+nAtoms1,conf2->getAtomPos(i)+offset);
          }
        } catch (ConformerException &ce) {
          BOOST_LOG(rdWarningLog)<<"combineMols: conformer id "<<conf1->getId()<<" not found in mol2";
        }
      }
    }
    res->clearComputedProps(true);
    return (ROMol *)res;
  }    

  void addRecursiveQueries(ROMol &mol,const std::map<std::string,ROMOL_SPTR> &queries,std::string propName,
                             std::vector<std::pair<unsigned int, std::string> > *reactantLabels){
    std::string delim=",";
    boost::char_separator<char> sep(delim.c_str());
    if (reactantLabels!=NULL) {
      (*reactantLabels).resize(0);
    }

    ROMol::VERTEX_ITER atBegin,atEnd;
    boost::tie(atBegin,atEnd) = mol.getVertices();
    while(atBegin!=atEnd){
      Atom *at=mol[*atBegin].get();
      ++atBegin;
      if(!at->hasProp(propName)) continue;
      std::string pval;
      at->getProp(propName,pval);
      boost::algorithm::to_lower(pval);
      if (reactantLabels!=NULL) {
        std::pair<unsigned int, std::string> label (at->getIdx(), pval);
        (*reactantLabels).push_back(label);
      }

      QueryAtom::QUERYATOM_QUERY *qToAdd;
      if(pval.find(delim)!=std::string::npos){
        boost::tokenizer<boost::char_separator<char> > tokens(pval,sep);
        boost::tokenizer<boost::char_separator<char> >::iterator token;
        qToAdd = new ATOM_OR_QUERY();
        for(token=tokens.begin();token!=tokens.end();++token){
          std::map<std::string,ROMOL_SPTR>::const_iterator iter=queries.find(*token);
          if(iter==queries.end()) {
            throw KeyErrorException(pval);
          }
          RecursiveStructureQuery *tqp= new RecursiveStructureQuery(new ROMol(*(iter->second)));
          boost::shared_ptr<RecursiveStructureQuery> nq(tqp);
          qToAdd->addChild(nq);
        }
      } else {
        std::map<std::string,ROMOL_SPTR>::const_iterator iter=queries.find(pval);
        if(iter==queries.end()) {
          throw KeyErrorException(pval);
        }
        qToAdd= new RecursiveStructureQuery(new ROMol(*(iter->second)));
      }
      if(!at->hasQuery()){
        QueryAtom qAt(*at);
        static_cast<RWMol &>(mol).replaceAtom(at->getIdx(),&qAt);
        at = mol.getAtomWithIdx(at->getIdx());
      }
      at->expandQuery(qToAdd,Queries::COMPOSITE_AND);
    }
  }

  void parseQueryDefFile(std::istream *inStream,std::map<std::string,ROMOL_SPTR> &queryDefs,
                         bool standardize,std::string delimiter,std::string comment,
                         unsigned int nameColumn,unsigned int smartsColumn){
    PRECONDITION(inStream,"no stream");
    queryDefs.clear();

    boost::char_separator<char> sep(delimiter.c_str());
    unsigned int line=0;
    std::string tempStr;
    while(!inStream->eof()){
      line++;
      tempStr=getLine(inStream);
      if(tempStr=="" || tempStr.find(comment)==0 ) continue;
      boost::tokenizer<boost::char_separator<char> > tokens(tempStr,sep);
      unsigned int tpos;
      boost::tokenizer<boost::char_separator<char> >::iterator token;
      std::string qname="";
      std::string qval="";
      for(token=tokens.begin(),tpos=0;token!=tokens.end();++token,++tpos){
        if(tpos==nameColumn){
          qname=*token;
        } else if(tpos==smartsColumn){
          qval=*token;
        }
      }
      boost::trim_if(qname,boost::is_any_of(" \t"));
      boost::trim_if(qval,boost::is_any_of(" \t"));
      if(qname=="" || qval==""){
        continue;
      }
      RWMol *m=0;
      try{
        m = SmartsToMol(qval);
      } catch (...) {
        m=0;
      }
      if(!m){
        BOOST_LOG(rdWarningLog)<<"cannot convert SMARTS "<<qval<<" to molecule at line "<<line<<std::endl;
        continue;
      }
      ROMOL_SPTR msptr(m);
      if (standardize) {
          boost::algorithm::to_lower(qname);
      }
      queryDefs[qname]=msptr;
    }
  }
  void parseQueryDefFile(std::string filename,std::map<std::string,ROMOL_SPTR> &queryDefs,
                         bool standardize,std::string delimiter,std::string comment,
                         unsigned int nameColumn,unsigned int smartsColumn){
    std::ifstream inStream(filename.c_str());
    if (!inStream || (inStream.bad()) ) {
      std::ostringstream errout;
      errout << "Bad input file " << filename;
      throw BadFileException(errout.str());
    }
    parseQueryDefFile(&inStream,queryDefs,standardize,delimiter,comment,nameColumn,smartsColumn);
  }
  void parseQueryDefText(const std::string &queryDefText,
                         std::map<std::string,ROMOL_SPTR> &queryDefs,
                         bool standardize,std::string delimiter,std::string comment,
                         unsigned int nameColumn,unsigned int smartsColumn){
    std::stringstream inStream(queryDefText);
    parseQueryDefFile(&inStream,queryDefs,standardize,delimiter,comment,nameColumn,smartsColumn);
  }
}  // end of namespace RDKit
