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

#include <RDGeneral/hash/hash.hpp>

namespace RDKit {
  namespace Subgraphs {
    ROMol *pathToSubmol(const ROMol &mol, const PATH_TYPE &path,
                        bool useQuery) {
      INT_MAP_INT aIdxMap;
      return pathToSubmol(mol, path, useQuery, aIdxMap);
    }
 
    ROMol *pathToSubmol(const ROMol &mol, const PATH_TYPE &path, 
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

    using boost::uint32_t;
    using boost::int32_t;
    DiscrimTuple
    calcPathDiscriminators(const ROMol &mol, const PATH_TYPE &path, bool useBO,
                            std::vector<boost::uint32_t> *extraInvars) {
      if(extraInvars) CHECK_INVARIANT(extraInvars->size()==mol.getNumAtoms(),"bad extra invars");
      DiscrimTuple res;
      std::vector<int32_t> atomsUsed(mol.getNumAtoms(),-1);
      std::vector<const Atom *> atoms;
      std::vector<uint32_t> pathDegrees;
      for(PATH_TYPE::const_iterator pathIter=path.begin(); pathIter!=path.end(); pathIter++){
        const Bond *bond=mol.getBondWithIdx(*pathIter);
        if(atomsUsed[bond->getBeginAtomIdx()]<0){
          atomsUsed[bond->getBeginAtomIdx()]=atoms.size();
          atoms.push_back(bond->getBeginAtom());
          pathDegrees.push_back(1);
        } else {
          pathDegrees[atomsUsed[bond->getBeginAtomIdx()]]+=1;
        }
        if(atomsUsed[bond->getEndAtomIdx()]<0){
          atomsUsed[bond->getEndAtomIdx()]=atoms.size();
          atoms.push_back(bond->getEndAtom());
          pathDegrees.push_back(1);
        } else {
          pathDegrees[atomsUsed[bond->getEndAtomIdx()]]+=1;
        }
      }

      unsigned int nAtoms=atoms.size();
      std::vector<uint32_t> invars(nAtoms);
      for(unsigned int i=0;i<nAtoms;++i){
        const Atom *atom=atoms[i];
        uint32_t invar=atom->getAtomicNum();
        gboost::hash_combine(invar,pathDegrees[i]);
        gboost::hash_combine(invar,atom->getFormalCharge());
        int deltaMass = static_cast<int>(atom->getMass() -
                                         PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
        gboost::hash_combine(invar,deltaMass);
        if(atom->getIsAromatic()){
          gboost::hash_combine(invar,1);
        }
        if(extraInvars){
          gboost::hash_combine(invar,(*extraInvars)[atom->getIdx()]);
        }
        invars[i] = invar;
      }
#if 0
      std::cerr<<"-------------------"<<std::endl;;
      std::cerr<<"      ats: >>>";
      std::copy(atomsUsed.begin(),atomsUsed.end(),std::ostream_iterator<int>(std::cerr," "));
      std::cerr<<std::endl;
#endif
      unsigned int nCycles=path.size()/2+1;
      gboost::hash<std::vector<uint32_t> > vectHasher;

      for(unsigned int cycle=0;cycle<nCycles;++cycle){
#if 0
        std::cerr<<"             round: >>>"<<" "<<cycle<<": ";
        std::copy(invars.begin(),invars.end(),std::ostream_iterator<int>(std::cerr," "));
        std::cerr<<std::endl;
#endif
        std::vector< std::vector<uint32_t> > locInvars(nAtoms);
        for(PATH_TYPE::const_iterator pathIter=path.begin(); pathIter!=path.end(); ++pathIter){
          const Bond *bond=mol.getBondWithIdx(*pathIter);
          uint32_t v1=invars[atomsUsed[bond->getBeginAtomIdx()]];
          uint32_t v2=invars[atomsUsed[bond->getEndAtomIdx()]];
          if(useBO){
            gboost::hash_combine(v1,static_cast<uint32_t>(bond->getBondType()));
            gboost::hash_combine(v2,static_cast<uint32_t>(bond->getBondType()));
          }
          locInvars[atomsUsed[bond->getBeginAtomIdx()]].push_back(v2);
          locInvars[atomsUsed[bond->getEndAtomIdx()]].push_back(v1);
        }
        for(unsigned int i=0;i<nAtoms;++i){
          std::sort(locInvars[i].begin(),locInvars[i].end());
          invars[i]=vectHasher(locInvars[i]);
        }
      }  

#if 0
      std::cerr<<"      invars: >>>";
      std::copy(invars.begin(),invars.end(),std::ostream_iterator<int>(std::cerr," "));
      std::cerr<<std::endl;
#endif
      std::sort(invars.begin(),invars.end());
      uint32_t pathInvar=vectHasher(invars);
#if 0
      std::cerr<<"      >>>";
      std::copy(path.begin(),path.end(),std::ostream_iterator<int>(std::cerr," "));
      std::cerr<<std::endl;
      std::cerr<<"         >>>"<<pathInvar<<","<<path.size()<<","<<nAtoms<<std::endl;
#endif
      return boost::make_tuple(pathInvar,path.size(),nAtoms);
    }

    //
    // This is intended for use on either subgraphs or paths.
    //  The entries in PATH_LIST should refer to bonds though (not
    //  atoms)
    //
    PATH_LIST uniquifyPaths (const ROMol &mol, const PATH_LIST &allPaths,
                             bool useBO){
      PATH_LIST res;
      std::vector<DiscrimTuple> discrimsSeen;
      for(PATH_LIST::const_iterator path=allPaths.begin();
          path!=allPaths.end();++path){
        DiscrimTuple discrims = calcPathDiscriminators(mol,*path,useBO);
        if(std::find(discrimsSeen.begin(),discrimsSeen.end(),discrims)==discrimsSeen.end()){
          discrimsSeen.push_back(discrims);
          res.push_back(*path);
        }
      }
      return res;
    }
  } // end of namespace Subgraphs
} // end of namespace RDKit
