// $Id$
//
//  Copyright (C) 2003-2009 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>

#include "Subgraphs.h"
#include "SubgraphUtils.h"

#include <RDGeneral/utils.h>
#include <RDBoost/Exceptions.h>

#include <iostream>
#include <cstring>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>

namespace RDKit {
namespace Subgraphs {
  void getNbrsList(const ROMol &mol, bool useHs, INT_INT_VECT_MAP &nbrs) {
    nbrs.clear();
    int nAtoms = mol.getNumAtoms();

    // create a list of neighbors for each bond
    // The python version (subgraph.py) did not take care of atoms being
    // hydrogens (and if a user would like to ignore them)
    // a list of list of bonds is no longer appropriate here
    // some bond IDs may be missing to index the list on.
    // so using an associative container.

    for (int i = 0; i < nAtoms; i++) {
      const Atom *atom = mol.getAtomWithIdx(i);
      // if are at a hydrogen and we are not interested in bonds connecting to them
      // move on
      if( useHs || atom->getAtomicNum()!=1 ){
        ROMol::OEDGE_ITER bIt1,end;
        boost::tie(bIt1,end) = mol.getAtomBonds(atom);
        while(bIt1!=end){
          const BOND_SPTR bond1 = mol[*bIt1];
          // if this bond connect to a hydrogen and we are not interested
          // in it ignore 
          if( useHs || bond1->getOtherAtom(atom)->getAtomicNum() != 1 ){
            int bid1 = bond1->getIdx();
            if (nbrs.find(bid1) == nbrs.end()) {
              INT_VECT nlst;
              nbrs[bid1] = nlst;
            }
            ROMol::OEDGE_ITER bIt2= mol.getAtomBonds(atom).first;
            while(bIt2 != end){
              const BOND_SPTR bond2 = mol[*bIt2];
              int bid2 = bond2->getIdx();
              if (bid1 != bid2 &&
                  (useHs || bond2->getOtherAtom(atom)->getAtomicNum() != 1 ) ){
                nbrs[bid1].push_back(bid2); //FIX: pathListType should probably be container of pointers ??
              }
              ++bIt2;
            }
          }
          ++bIt1;
        }
      }
    }
  }

  // Everything passed here by reference                               
  void recurseWalk(INT_INT_VECT_MAP &nbrs, // neighbors for each bond
                   PATH_TYPE &spath, // the current path to be build upon
                   INT_VECT &cands, // neighbors of current path
                   unsigned int targetLen, // the maximum subgraph len we are interested in
                   boost::dynamic_bitset<> forbidden, // bonds that have been covered already
                   // we don't want reference passing for forbidden, 
                   // it gets altered through the processand we want  
                   // fresh start everytime we buble back up to "FindAllSubGraphs"
                   PATH_LIST &res // the final list of subgraphs 
                   ) 
  {
    // end case for recursion
    if (spath.size() == targetLen) {
      res.push_back(spath);
      return;
    }

    // if the path is already bigger than target length don't do anything
    if (spath.size() > targetLen) {
      return;
    }
  
    // we  have the cndidates that can be used to add to the existing path
    // try extending the subgraphs
    while (cands.size() != 0) { 
      int next = cands.back(); // start with the last one in the candidate list
      cands.pop_back();
      //cands.erase(remove(cands.begin(), cands.end(), next), cands.end());
      if (!forbidden[next]){
        // this bond should not appear in the later subgraphs
        forbidden[next]=1;
      
        // update a local stack before the next recursive call
        INT_VECT tstack = cands;
        for (INT_VECT::iterator bid=nbrs[next].begin(); bid != nbrs[next].end(); bid++) {
          if (!forbidden[*bid]){
            tstack.push_back(*bid);
          }
        }
      
        PATH_TYPE tpath = spath;
        tpath.push_back(next);

        recurseWalk(nbrs,tpath, tstack, targetLen, forbidden, res);
      }
    }
  } 


  // Everything passed here by reference                               
  void recurseWalkRange(INT_INT_VECT_MAP &nbrs, // neighbors for each bond
                        PATH_TYPE &spath, // the current path to be build upon
                        INT_VECT &cands, // neighbors of current path
                        unsigned int lowerLen, // lower limit of the subgraph lengths we are interested in
                        unsigned int upperLen, // the maximum subgraph len we are interested in
                        boost::dynamic_bitset<> forbidden, // bonds that have been covered already
                        // we don't want reference passing for forbidden, 
                        // it gets altered through the processand we want  
                        // fresh start everytime we buble back up to "FindAllSubGraphs"
                        INT_PATH_LIST_MAP &res // the final list of subgraphs 
                        ) 
  {
    unsigned int nsize = spath.size();
    if ((nsize >= lowerLen) && (nsize <= upperLen)) {
      //if (res.find(nsize) == res.end()) {
      //  PATH_LIST ordern;
      //  res[nsize] = ordern;
      //}
      res[nsize].push_back(spath);
    }
  
    // end case for recursion
    if (nsize == upperLen) {
      return;
    }

    // if the path is already bigger than desired size
    if (nsize > upperLen) {
      return;
    }
  
    // we  have the cndidates that can be used to add to the existing path
    // try extending the subgraphs
    while (cands.size() != 0) { 
      int next = cands.back(); // start with the last one in the candidate list
      cands.pop_back();
      //cands.erase(remove(cands.begin(), cands.end(), next), cands.end());
      if (!forbidden[next]){
        // this bond should not appear in the later subgraphs
        forbidden[next]=1;
      
        // update a local stack before the next recursive call
        INT_VECT tstack = cands;
        for (INT_VECT::iterator bid=nbrs[next].begin(); bid != nbrs[next].end(); bid++) {
          if (!forbidden[*bid]){
            tstack.push_back(*bid);
          }
        }
      
        PATH_TYPE tpath = spath;
        tpath.push_back(next);

        recurseWalkRange(nbrs,tpath, tstack, lowerLen, upperLen, forbidden, res);
      }
    }
  }

  void dumpVIV(VECT_INT_VECT v){
    VECT_INT_VECT::iterator i;
    INT_VECT::iterator j;
    for(i=v.begin();i!=v.end();i++){
      for(j=i->begin();j!=i->end();j++){
        std::cout << *j << " ";
      }
      std::cout << std::endl;
    }
  }
  
  PATH_LIST
  extendPaths(int *adjMat,unsigned int dim,const PATH_LIST &paths,int allowRingClosures=-1)
  {
    PRECONDITION(adjMat,"no matrix");
    //
    //  extend each of the currently active paths by adding
    //   a single adjacent index to the end of each
    //
    PATH_LIST res;
    PATH_LIST::const_iterator path;
    for(path=paths.begin(); path!=paths.end(); path++){
      unsigned int endIdx = (*path)[path->size()-1];
      unsigned int iTab = endIdx*dim;
      for(unsigned int otherIdx = 0; otherIdx < dim; otherIdx++){
        if( adjMat[iTab+otherIdx] == 1){
          // test 1: make sure the new atom is not already
          //   in the path
          PATH_TYPE::const_iterator loc;
          loc = std::find(path->begin(),path->end(),static_cast<int>(otherIdx));
          // The two conditions for adding the atom are:
          //   1) it's not there already
          //   2) it's there, but ring closures are allowed and this
          //      will be the last addition to the path.
          if ( loc == path->end() ){
            // the easy case
            //PATH_TYPE newPath=*path;
            //newPath.push_back(otherIdx);
            //res.push_back(newPath);
            res.push_back(*path);
            res.rbegin()->push_back(otherIdx);
          } else if (allowRingClosures>2 &&
                     static_cast<int>(path->size())==allowRingClosures-1) {
            // We *might* be adding the atom, but we need to make sure
            // that we're not just duplicating the second to last
            // element of the path:
            PATH_TYPE::const_reverse_iterator rIt=path->rbegin();
            rIt++;
            if( *rIt != static_cast<int>(otherIdx) ){
              //PATH_TYPE newPath=*path;
              //newPath.push_back(otherIdx);
              //res.push_back(newPath);
              res.push_back(*path);
              res.rbegin()->push_back(otherIdx);
            }
          }
        }
      }
    }
    return res;
  }

  INT_PATH_LIST_MAP
  pathFinderHelper(int *adjMat,unsigned int dim,unsigned int minLen,
                   unsigned int maxLen,int rootedAtAtom)
  {
    PRECONDITION(adjMat,"no matrix");
    PRECONDITION(minLen<=maxLen,"bad lengths provided");
    // finds all paths of length N using an adjacency matrix,
    //  which is constructed elsewhere
    INT_PATH_LIST_MAP res;
    PATH_LIST paths;
    paths.clear();

    if(rootedAtAtom<0){
      // start a path at each possible index
      for(unsigned int i=0;i<dim;i++){
        PATH_TYPE tPath;
        tPath.push_back(i);
        paths.push_back(tPath);
      }
    } else if(rootedAtAtom<static_cast<int>(dim)) {
      // only start a path at the atom of interest:
      PATH_TYPE tPath;
      tPath.push_back(rootedAtAtom);
      paths.push_back(tPath);
    } else {
      return res;
    }

    // and build them up one index at a time:
    for(unsigned int length=1;length<maxLen;length++){
      // extend each path:
      if(length>=minLen) res[length]=paths;
      paths = extendPaths(adjMat,dim,paths,maxLen);
    }
    res[maxLen]=paths;
    
    return res;
  }
} // end of Subgraphs namespace

  PATH_LIST findAllSubgraphsOfLengthN (const ROMol &mol, unsigned int targetLen,
                                       bool useHs,int rootedAtAtom){
    /*********************************************
      FIX: Lots of issues here:
      - pathListType is defined as a container of "pathType", should it be a container
      of "pointers to pathtype"
      - to make few things clear it might be useful to typdef a "subgraphListType"
      even if it is exactly same as the "pathListType", just to not confuse between
      path vs. subgraph definitions
      - To make it consistent with the python version of this function in "subgraph.py"
      it return a "list of paths" instead of a "list of list of paths" (see 
      "GetPathsUpTolength" in "molgraphs.cpp")
    ****************************************************************************/
    boost::dynamic_bitset<> forbidden(mol.getNumBonds());
  
    // this should be the only dependence on mol object:
    INT_INT_VECT_MAP nbrs;
    Subgraphs::getNbrsList(mol, useHs,nbrs); 
  
    // Start path at each bond
    PATH_LIST res;

    // start paths at each bond:
    for (INT_INT_VECT_MAP::iterator nbi = nbrs.begin();
         nbi != nbrs.end(); ++nbi) {
      // don't come back to this bond in the later subgraphs
      int i = (*nbi).first;

      // if we're only returning paths rooted at a particular atom, check now
      // that this bond involves that atom:
      if(rootedAtAtom>=0 &&
         mol.getBondWithIdx(i)->getBeginAtomIdx()!=static_cast<unsigned int>(rootedAtAtom) &&
         mol.getBondWithIdx(i)->getEndAtomIdx()!=static_cast<unsigned int>(rootedAtAtom) ){
        continue;
      }

      if (forbidden[i]){
        continue;
      }
      forbidden[i]=1;
      
      // start the recursive path building with the current bond
      PATH_TYPE spath;
      spath.clear();
      spath.push_back(i);
      
      // neighbors of this bond are the next candidates
      INT_VECT cands = nbrs[i];
      
      // now call the recursive function
      // little bit different from the python version 
      // the result list of paths is passed as a reference, instead of on the fly 
      // appending 
      Subgraphs::recurseWalk(nbrs, spath, cands, targetLen, forbidden, res);
    }
    nbrs.clear();
    return res;
  }


  INT_PATH_LIST_MAP findAllSubgraphsOfLengthsMtoN(const ROMol &mol, unsigned int lowerLen,
                                                  unsigned int upperLen, bool useHs,int rootedAtAtom){
    PRECONDITION(lowerLen <= upperLen, "");
    boost::dynamic_bitset<> forbidden(mol.getNumBonds());


    INT_INT_VECT_MAP nbrs;
    Subgraphs::getNbrsList(mol, useHs,nbrs);
  
    // Start path at each bond
    INT_PATH_LIST_MAP res;
    for (unsigned int idx = lowerLen; idx <= upperLen; idx++) {
      PATH_LIST ordern;
      res[idx] = ordern;
    }

    // start paths at each bond:
    for (INT_INT_VECT_MAP::iterator nbi = nbrs.begin();
         nbi != nbrs.end(); nbi++) {
      int i = (*nbi).first;

      // if we're only returning paths rooted at a particular atom, check now
      // that this bond involves that atom:
      if(rootedAtAtom>=0 &&
         mol.getBondWithIdx(i)->getBeginAtomIdx()!=static_cast<unsigned int>(rootedAtAtom) &&
         mol.getBondWithIdx(i)->getEndAtomIdx()!=static_cast<unsigned int>(rootedAtAtom) ){
        continue;
      }

      // don't come back to this bond in the later subgraphs
      if (forbidden[i]){
        continue;
      }
      forbidden[i]=1;
      
      // start the recursive path building with the current bond
      PATH_TYPE spath;
      spath.clear();
      spath.push_back(i);

      // neighbors of this bond are the next candidates
      INT_VECT cands = nbrs[i];
      
      // now call the recursive function
      // little bit different from the python version 
      // the result list of paths is passed as a reference, instead of on the fly 
      // appending 
      Subgraphs::recurseWalkRange(nbrs, spath, cands, lowerLen, upperLen, forbidden, res);
    }
    nbrs.clear();
    return res; //FIX : need some verbose testing code here
  }
  
  PATH_LIST findUniqueSubgraphsOfLengthN (const ROMol &mol, unsigned int targetLen,
                                          bool useHs,bool useBO,int rootedAtAtom) 
  {
    // start by finding all subgraphs, then uniquify
    PATH_LIST allSubgraphs=findAllSubgraphsOfLengthN(mol,targetLen,useHs,rootedAtAtom);
    PATH_LIST res = Subgraphs::uniquifyPaths(mol,allSubgraphs,useBO);
    return res;
  }

  // ----------------------------------------------
  //
  //  You may find yourself wondering: "what's the difference between 
  //  a subgraph and path?"  Well, let me tell you: there's a big
  //  diffference! 
  //
  //  Subgraphs are potentially branched, whereas paths (in our
  //  terminology at least) cannot be.  So, the following graph:
  //
  //            C--0--C--1--C--3--C
  //                  |
  //                  2
  //                  |
  //                  C
  //  has 3 subgraphs of length 3: (0,1,2),(0,1,3),(2,1,3)
  //  but only 2 paths of length 3: (0,1,3),(2,1,3)
  //
  //
  // ----------------------------------------------
  //
  //  Args:
  //
  //    adjMat: the adjacency matrix
  //    dim: number of rows in the adjacency matrix
  //    paths: paths to be extended
  //    allowRingClosures: if > 2, paths will be allowed where the final
  //        point is a duplicate (closes rings).  In this case,
  //        allowRingClosures should equal the target path length.
  //
  INT_PATH_LIST_MAP 
  findAllPathsOfLengthsMtoN(const ROMol &mol,unsigned int lowerLen,
                           unsigned int upperLen,bool useBonds,
                           bool useHs,int rootedAtAtom) {
    //
    //  We can't be clever here and just use the bond adjacency matrix
    //  to solve this problem when useBonds is true.  This is because
    //  the bond adjacency matrices for the molecules C1CC1 and CC(C)C
    //  are indistinguishable.  In the second case, t-butane (and
    //  anything else with a T junction), we'll get some subgraphs mixed
    //  in with the paths.  So we have to construct paths of atoms and
    //  then convert them into bond paths.
    //
    PRECONDITION(lowerLen <= upperLen, "");

    int *adjMat,dim;
    dim = mol.getNumAtoms();
    adjMat = new int[dim*dim];
    memset((void *)adjMat,0,dim*dim*sizeof(int));
  
    // generate the adjacency matrix by hand by looping over the bonds
    ROMol::ConstBondIterator bondIt;
    for(bondIt=mol.beginBonds();bondIt!=mol.endBonds();bondIt++){
      Atom *beg=(*bondIt)->getBeginAtom();
      Atom *end=(*bondIt)->getEndAtom();
      // check for H, which we might be skipping
      if(useHs || (beg->getAtomicNum()!=1 && end->getAtomicNum()!=1)){
        adjMat[beg->getIdx()*dim+end->getIdx()] = 1;
        adjMat[end->getIdx()*dim+beg->getIdx()] = 1;
      }
    }

    // if we're using bonds, we'll need to find paths of length N+1,
    // then convert them
    if(useBonds) {
      ++lowerLen;
      ++upperLen;
    }

    // find the paths themselves
    INT_PATH_LIST_MAP atomPaths=Subgraphs::pathFinderHelper(adjMat,dim,lowerLen,upperLen,
                                                            rootedAtAtom);

    // clean up the adjacency matrix
    delete [] adjMat;

    INT_PATH_LIST_MAP res;

    //
    //--------------------------------------------------------
    // loop through all the paths we have and make sure that there are
    // no duplicates (duplicate = contains identical bond indices)
    //
    //  We need to use the bond paths for this duplicate finding
    //  because, in rings, there can be many paths which share atom
    //  indices but which have different bond compositions. For example,
    //  there is only one "atom unique" path of length 5 bonds (6 atoms)
    //  through a 6-ring, but there are six bond paths.
    //
    if(!useBonds && lowerLen>=1){
      res[1]=atomPaths[1];
    }
    if(useBonds || upperLen>1){
      for(unsigned int i=lowerLen;i<=upperLen;++i){
        if(i<=1){
          continue;
        }

        std::vector< boost::dynamic_bitset<> > invars;

        for(PATH_LIST::const_iterator vivI=atomPaths[i].begin();
            vivI!=atomPaths[i].end();++vivI){
          boost::dynamic_bitset<> invar(mol.getNumBonds());
          const PATH_TYPE &resi=*vivI;
          PATH_TYPE locV;
          locV.reserve(i);
          for(unsigned int j=0;j<i-1;j++){
            const Bond *bond=mol.getBondBetweenAtoms(resi[j],resi[j+1]);
            locV.push_back(bond->getIdx());
            invar.set(bond->getIdx());
          }
          if(std::find(invars.begin(),invars.end(),invar)==invars.end()){
            invars.push_back(invar);
            if(useBonds){
              res[i-1].push_back(locV);
            } else {
              res[i].push_back(resi);
            }
          }
        }
      }
    }
    return res;
  }
  PATH_LIST
  findAllPathsOfLengthN(const ROMol &mol,unsigned int targetLen,bool useBonds,
                        bool useHs,int rootedAtAtom) {
    return findAllPathsOfLengthsMtoN(mol,targetLen,targetLen,useBonds,useHs,rootedAtAtom)[targetLen];
  }

  PATH_TYPE findAtomEnvironmentOfRadiusN(const ROMol &mol,unsigned int radius,
                                         unsigned int rootedAtAtom,bool useHs){
    if(rootedAtAtom>=mol.getNumAtoms()) throw ValueErrorException("bad atom index");

    PATH_TYPE res;
    std::list< std::pair<int,int> > nbrStack;
    ROMol::OEDGE_ITER beg,end;
    boost::tie(beg,end) = mol.getAtomBonds(mol.getAtomWithIdx(rootedAtAtom));
    while(beg!=end){
      BOND_SPTR bond=mol[*beg];
      if(useHs || mol.getAtomWithIdx(bond->getOtherAtomIdx(rootedAtAtom))->getAtomicNum()!=1){
        nbrStack.push_back(std::make_pair(rootedAtAtom,bond->getIdx()));
      }
      ++beg;
    }
    boost::dynamic_bitset<> bondsIn(mol.getNumBonds());
    unsigned int i;
    for(i=0;i<radius;++i){
      if(nbrStack.empty()){
        break;
      }

      std::list< std::pair<int,int> > nextLayer;
      while(!nbrStack.empty()){
        int bondIdx,startAtom;
        boost::tie(startAtom,bondIdx)=nbrStack.front();
        nbrStack.pop_front();
        if(!bondsIn.test(bondIdx)){
          bondsIn.set(bondIdx);
          res.push_back(bondIdx);

          // add the next set of neighbors:
          int oAtom=mol.getBondWithIdx(bondIdx)->getOtherAtomIdx(startAtom);
          boost::tie(beg,end) = mol.getAtomBonds(mol.getAtomWithIdx(oAtom));
          while(beg!=end){
            BOND_SPTR bond=mol[*beg];
            if(!bondsIn.test(bond->getIdx())){
              if(useHs || mol.getAtomWithIdx(bond->getOtherAtomIdx(oAtom))->getAtomicNum()!=1){
                nextLayer.push_back(std::make_pair(oAtom,bond->getIdx()));
              }
            }
            ++beg;
          }
        }
      }
      nbrStack=nextLayer;
    }
    if(i!=radius){
      // this happens when there are no paths with the requested radius.
      // return nothing in this case:
      res.clear();
      res.resize(0);
    }

    return res;
  }

} // end of RDKit namespace
