// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "SubgraphUtils.h"
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <iostream>
#include <algorithm>
#include <map>

namespace RDKit {
  namespace Subgraphs {
ROMol *PathToSubmol(const ROMol &mol, const PATH_TYPE &path,
                    bool useQuery) {
  INT_MAP_INT aIdxMap;
  return PathToSubmol(mol, path, useQuery, aIdxMap);
}
 
ROMol *PathToSubmol(const ROMol &mol, const PATH_TYPE &path, 
                    bool useQuery,
                    INT_MAP_INT &atomIdxMap) {
  RWMol *subMol=new RWMol();
  PATH_TYPE::const_iterator pathIter;
  //std::map<int,int> atomIdxMap;
  atomIdxMap.clear();

  if (useQuery) {
    // have to do this in two different blocks because of issues with variable scopes.
    for(pathIter=path.begin(); pathIter!=path.end(); ++pathIter){
      QueryBond *bond;
      bond = new QueryBond(*(mol.getBondWithIdx(*pathIter)));
    
      int begIdx=bond->getBeginAtomIdx();
      int endIdx=bond->getEndAtomIdx();
      
      if(atomIdxMap.find(begIdx)==atomIdxMap.end()){
        QueryAtom *atom = new QueryAtom(*(mol.getAtomWithIdx(begIdx)));
        int newAtomIdx=subMol->addAtom(atom,false,true);
        atomIdxMap[begIdx] = newAtomIdx;
      }
      begIdx = atomIdxMap.find(begIdx)->second;
      if(atomIdxMap.find(endIdx)==atomIdxMap.end()){
        QueryAtom *atom = new QueryAtom(*(mol.getAtomWithIdx(endIdx)));
        int newAtomIdx=subMol->addAtom(atom,false,true);
        atomIdxMap[endIdx] = newAtomIdx;
      }
      endIdx = atomIdxMap.find(endIdx)->second;
      
      bond->setOwningMol(subMol);
      bond->setBeginAtomIdx(begIdx);
      bond->setEndAtomIdx(endIdx);
      subMol->addBond(bond,true);
    }
  }
  else {
    for(pathIter=path.begin(); pathIter!=path.end(); ++pathIter){
      Bond *bond;
      bond=mol.getBondWithIdx(*pathIter)->copy();
      
      int begIdx=bond->getBeginAtomIdx();
      int endIdx=bond->getEndAtomIdx();
      
      if(atomIdxMap.find(begIdx)==atomIdxMap.end()){
        Atom *atom = mol.getAtomWithIdx(begIdx)->copy();
        int newAtomIdx=subMol->addAtom(atom,false,true);
        atomIdxMap[begIdx] = newAtomIdx;
      }
      begIdx = atomIdxMap.find(begIdx)->second;
      if(atomIdxMap.find(endIdx)==atomIdxMap.end()){
        Atom *atom = mol.getAtomWithIdx(endIdx)->copy();
        int newAtomIdx=subMol->addAtom(atom,false,true);
        atomIdxMap[endIdx] = newAtomIdx;
      }
      endIdx = atomIdxMap.find(endIdx)->second;
      
      bond->setOwningMol(subMol);
      bond->setBeginAtomIdx(begIdx);
      bond->setEndAtomIdx(endIdx);
      subMol->addBond(bond,true);
    }
  }
  return subMol;
}

  PATH_TYPE bondListFromAtomList(const ROMol &mol, const PATH_TYPE &atomIds) {
    PATH_TYPE bids;
    unsigned int natms = atomIds.size();
    if (natms <= 1) {
      return bids; //FIX: should probably throw an exception
    }
    for (unsigned int i = 0; i < natms; i++) {
      for (unsigned int j = i+1; j < natms; j++) {
        const Bond *bnd = mol.getBondBetweenAtoms(atomIds[i], atomIds[j]);
        if (bnd) {
          int bid = bnd->getIdx();
          bids.push_back(bid);
        }
      }
    }
    return bids;
  }

PathDiscrimTuple
CalcPathDiscriminators(const ROMol &mol, const PATH_TYPE &path, bool useBO) {

  //---------------------------------------------------
  // start by building a molecule consisting of only the atoms and
  // bonds in the path being considered.
  //
  ROMol *subMol=PathToSubmol(mol,path);

  PathDiscrimTuple res = MolOps::computeDiscriminators(*subMol, useBO);
  
  delete subMol;
  return res;
}

  bool operator==(const PathDiscrimTuple &t1,const PathDiscrimTuple &t2){
#ifndef WIN32
    if( !feq(t1.get<0>(),t2.get<0>()) ){
      return false;
    } else if( !feq(t1.get<1>(),t2.get<1>()) ) {
      return false;
    } else if( !feq(t1.get<2>(),t2.get<2>()) ) {
      return false;
    }
#else
    if( !feq(boost::tuples::get<0>(t1),boost::tuples::get<0>(t2)) ){
      return false;
    } else if( !feq(boost::tuples::get<1>(t1),boost::tuples::get<1>(t2)) ) {
      return false;
    } else if( !feq(boost::tuples::get<2>(t1),boost::tuples::get<2>(t2)) ) {
      return false;
    }
#endif
    return true;
  }
  bool operator!=(const PathDiscrimTuple &t1,const PathDiscrimTuple &t2){
    return !(t1==t2);
  }

  //
  // This is intended for use on either subgraphs or paths.
  //  The entries in PATH_LIST should refer to bonds though (not
  //  atoms)
  //
  PATH_LIST uniquifyPaths (const ROMol &mol, const PATH_LIST &allPaths,
                           bool useBO,double tol) {
    PATH_LIST res;
    std::vector<PathDiscrimTuple> discrimsSeen;
    for(PATH_LIST::const_iterator path=allPaths.begin();
        path!=allPaths.end();++path){
      PathDiscrimTuple discrims = CalcPathDiscriminators(mol,*path,useBO);
      bool found=false;
      for(std::vector<PathDiscrimTuple>::iterator discrimIt=discrimsSeen.begin();
          discrimIt!=discrimsSeen.end();
          ++discrimIt){
        if( feq(boost::tuples::get<0>(discrims),boost::tuples::get<0>(*discrimIt),tol) &&
            feq(boost::tuples::get<1>(discrims),boost::tuples::get<1>(*discrimIt),tol) &&
            feq(boost::tuples::get<2>(discrims),boost::tuples::get<2>(*discrimIt),tol)){
          found=true;
          break;
        }
      }
      if(!found){
        discrimsSeen.push_back(discrims);
        res.push_back(*path);
      }
    }
    return res;
  }
  } // end of namespace Subgraphs
} // end of namespace RDKit
