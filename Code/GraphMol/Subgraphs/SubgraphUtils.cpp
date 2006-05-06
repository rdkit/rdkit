// $Id: SubgraphUtils.cpp 4975 2006-02-18 00:52:04Z glandrum $
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
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

//lapack ++ includes
#include <lafnames.h>
#include <lapack.h>
#include <symd.h>
#include <lavd.h>
#include <laslv.h>
//#include <lapack++.h>

//
//  This is required because whoever wrote lapack++ must have decided
//  that no one would ever want to actually use the standard C++
//  throw.  Duh.
//
#ifdef throw
#undef throw
#endif

namespace RDKit {

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
    for(pathIter=path.begin(); pathIter!=path.end(); pathIter++){
      QueryAtom *atom;
      QueryBond *bond;
      bond = new QueryBond(*(mol.getBondWithIdx(*pathIter)));
    
      int begIdx,endIdx,newAtomIdx;
      begIdx=bond->getBeginAtomIdx();
      endIdx=bond->getEndAtomIdx();
      
      if(atomIdxMap.find(begIdx)==atomIdxMap.end()){
	atom = new QueryAtom(*(mol.getAtomWithIdx(begIdx)));
	newAtomIdx=subMol->addAtom(atom,false,true);
	atomIdxMap[begIdx] = newAtomIdx;
      }
      begIdx = atomIdxMap.find(begIdx)->second;
      if(atomIdxMap.find(endIdx)==atomIdxMap.end()){
	atom = new QueryAtom(*(mol.getAtomWithIdx(endIdx)));
	newAtomIdx=subMol->addAtom(atom,false,true);
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
    for(pathIter=path.begin(); pathIter!=path.end(); pathIter++){
      Atom *atom;
      Bond *bond;
      bond=mol.getBondWithIdx(*pathIter)->copy();
      
      int begIdx,endIdx,newAtomIdx;
      begIdx=bond->getBeginAtomIdx();
      endIdx=bond->getEndAtomIdx();
      
      if(atomIdxMap.find(begIdx)==atomIdxMap.end()){
	atom = mol.getAtomWithIdx(begIdx)->copy();
	newAtomIdx=subMol->addAtom(atom,false,true);
	atomIdxMap[begIdx] = newAtomIdx;
      }
      begIdx = atomIdxMap.find(begIdx)->second;
      if(atomIdxMap.find(endIdx)==atomIdxMap.end()){
	atom = mol.getAtomWithIdx(endIdx)->copy();
	newAtomIdx=subMol->addAtom(atom,false,true);
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

  PATH_TYPE bondListFromAtomList(const ROMol &mol, const PATH_TYPE atomIds) {
    PATH_TYPE bids;
    int natms = atomIds.size();
    const Bond *bnd;
    if (natms <= 1) {
      return bids; //FIX: should probably throw an exception
    }
    
    int i, j, bid;
    for (i = 0; i < natms; i++) {
      for (j = i+1; j < natms; j++) {
	bnd = mol.getBondBetweenAtoms(atomIds[i], atomIds[j]);
	if (bnd) {
	  bid = bnd->getIdx();
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
  //delete [] dMat;
  //return boost::make_tuple(J,ev1,ev2); 
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
    PATH_LIST::const_iterator path;
    std::vector<PathDiscrimTuple> discrimsSeen;
    bool found;
    for(path=allPaths.begin();path!=allPaths.end();path++){
      PathDiscrimTuple discrims = CalcPathDiscriminators(mol,*path,useBO);
      std::vector<PathDiscrimTuple>::iterator discrimIt;
      found=0;
      
      for(discrimIt=discrimsSeen.begin();
	  discrimIt!=discrimsSeen.end();
	  discrimIt++){
	if( feq(boost::tuples::get<0>(discrims),boost::tuples::get<0>(*discrimIt),tol) &&
	    feq(boost::tuples::get<1>(discrims),boost::tuples::get<1>(*discrimIt),tol) &&
	    feq(boost::tuples::get<2>(discrims),boost::tuples::get<2>(*discrimIt),tol)){
	  found=1;
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

} // end of namespace
