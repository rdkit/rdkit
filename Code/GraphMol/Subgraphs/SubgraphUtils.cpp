// $Id$
//
//  Copyright (C) 2003-2013 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SubgraphUtils.h"
#include "Subgraphs.h"
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <iostream>
#include <algorithm>
#include <map>
#include <boost/tuple/tuple_comparison.hpp>
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
      if(mol.getNumConformers()){
        // copy coordinates over:
        for(ROMol::ConstConformerIterator confIt=mol.beginConformers();
            confIt!=mol.endConformers();++confIt){
          Conformer *conf=new Conformer(subMol->getNumAtoms());
          conf->set3D((*confIt)->is3D());
          for(INT_MAP_INT::const_iterator mapIt=atomIdxMap.begin();
              mapIt!=atomIdxMap.end();++mapIt){
            conf->setAtomPos(mapIt->second,(*confIt)->getAtomPos(mapIt->first));
          }
          conf->setId((*confIt)->getId());
          subMol->addConformer(conf,false);
        }
      }
      // clear computed properties
      subMol->clearComputedProps(true);

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

      // Start by collecting the atoms in the path and their degrees
      std::vector<int32_t> atomsUsed(mol.getNumAtoms(),-1); // map from atom index->path index
      std::vector<const Atom *> atoms; // to contain the atoms in the path
      std::vector<uint32_t> pathDegrees; // degrees of each atom *in the path*
      for(PATH_TYPE::const_iterator pathIter=path.begin();
	  pathIter!=path.end(); ++pathIter){
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

      // Calculate the atomic invariants
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

      // now do the Morgan iterations:
      // the most number of cycles we need for the atoms on the edges
      // to feel each other is pathSize/2
      // EFF: it may be worth revisiting this at some point to see
      // if the iteration count can be even smaller (and if it
      // makes a difference in runtime)
      unsigned int nCycles=path.size()/2+1;
      gboost::hash<std::vector<uint32_t> > vectHasher;
      for(unsigned int cycle=0;cycle<nCycles;++cycle){
	// let each atom feel it's neighbors:
        std::vector< std::vector<uint32_t> > locInvars(nAtoms);
        for(PATH_TYPE::const_iterator pathIter=path.begin();
	    pathIter!=path.end(); ++pathIter){
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
	// we need to sort by the neighbor invariants to be order
	// independent:
        for(unsigned int i=0;i<nAtoms;++i){
          std::sort(locInvars[i].begin(),locInvars[i].end());
          invars[i]=vectHasher(locInvars[i]);
        }
      }  

      // again, a sort for order independence:
      std::sort(invars.begin(),invars.end());
      uint32_t pathInvar=vectHasher(invars);

      // also include the path size (bond count) and number of atoms
      // in the discriminator
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
