// $Id$
//
//  Copyright (C) 2003-2009 Greg Landrum and Rational Discovery LLC
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

#ifdef RDK_USELAPACKPP
//lapack ++ includes
#include <lafnames.h>
#include <lapack.h>
#include <symd.h>
#include <lavd.h>
#include <laslv.h>
#else
// uBLAS and boost.bindings includes
#include <boost/numeric/bindings/traits/ublas_matrix.hpp> 
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/bindings/lapack/syev.hpp>
#include <boost/numeric/ublas/io.hpp> 
namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;
#endif


namespace RDKit {
  namespace Subgraphs {

    DiscrimTuple computeDiscriminators(double *distMat, unsigned int nb, unsigned int na) {
      PRECONDITION(distMat,"bogus distance matrix");
      double ev1 = 0.0;
      double ev2 = 0.0;
#ifdef RDK_USELAPACKPP
      LaSymmMatDouble A(distMat, na, na);
      LaVectorDouble eigs(na);
      LaEigSolve(A, eigs);
#else
      ublas::matrix<double> A(na,na);
      ublas::vector<double> eigs(na);
      for(unsigned int i=0;i<na;++i){
        for(unsigned int j=i;j<na;++j){
          A(i,j)=distMat[i*na+j];
        }
      }
      lapack::syev('N','L',A,eigs);
#endif
      if (na > 1) {
        ev1 = eigs(0);
      }
      if (na > 2) {
        ev2 = eigs(na-1); 
      }
      double J = MolOps::computeBalabanJ(distMat, nb, na);
    
      return boost::make_tuple(J,ev1,ev2);
    }

    DiscrimTuple computeDiscriminators(const ROMol &mol, 
                                       bool useBO,
                                       bool force) {
      DiscrimTuple res;
      if ((mol.hasProp("Discrims")) && (!force)) {
        mol.getProp("Discrims", res);
      }
      else {
        unsigned int nAts = mol.getNumAtoms();
        double *dMat;
        unsigned int nb = mol.getNumBonds();
        dMat = MolOps::getDistanceMat(mol,useBO,true);
        //  Our discriminators (based upon eigenvalues of the distance matrix
        //  and Balaban indices) are good, but is a case they don't properly
        //  handle by default.  These two fragments:
        //    C-ccc and c-ccc
        //  give the same discriminators because their distance matrices are
        //  identical.  We'll work around this by adding 0.5 to the diagonal
        //  elements of the distance matrix corresponding to aromatic atoms:
        ROMol::ConstAromaticAtomIterator atomIt;
        for(atomIt=mol.beginAromaticAtoms();
            atomIt!=mol.endAromaticAtoms();
            atomIt++){
          unsigned int idx=(*atomIt)->getIdx();
          dMat[idx*nAts+idx] += 0.5;
        }
#if 0
        BOOST_LOG(rdDebugLog)<< "--------------------" << std::endl;
        for(int i=0;i<nAts;i++){
          for(int j=0;j<nAts;j++){
            BOOST_LOG(rdDebugLog)<< "\t" << std::setprecision(4) << dMat[i*nAts+j];
          }
          BOOST_LOG(rdDebugLog)<< std::endl;
        }
#endif
      
        res = computeDiscriminators(dMat, nb, nAts);
        mol.setProp("Discrims", res, true);
      }
      return res;
    }


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

DiscrimTuple
CalcPathDiscriminators(const ROMol &mol, const PATH_TYPE &path, bool useBO) {

  //---------------------------------------------------
  // start by building a molecule consisting of only the atoms and
  // bonds in the path being considered.
  //
  ROMol *subMol=PathToSubmol(mol,path);

  DiscrimTuple res = computeDiscriminators(*subMol, useBO);
  
  delete subMol;
  return res;
}

  bool operator==(const DiscrimTuple &t1,const DiscrimTuple &t2){
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
  bool operator!=(const DiscrimTuple &t1,const DiscrimTuple &t2){
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
    std::vector<DiscrimTuple> discrimsSeen;
    for(PATH_LIST::const_iterator path=allPaths.begin();
        path!=allPaths.end();++path){
      DiscrimTuple discrims = CalcPathDiscriminators(mol,*path,useBO);
      bool found=false;
      for(std::vector<DiscrimTuple>::iterator discrimIt=discrimsSeen.begin();
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
