//
//  Copyright (C) 2014 Greg Landrum
//  Adapted from pseudo-code from Roger Sayle
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "new_canon.h"
#include <GraphMol/RDKitBase.h>
#include <boost/cstdint.hpp>
#include <boost/foreach.hpp>
#include <cstring>
#include <iostream>
#include <cassert>

namespace RDKit {
  namespace Canon{
    void CreateSinglePartition(unsigned int nAtoms,
                               int *order,
                               int *count,
                               canon_atom *atoms) {
      for( unsigned int i=0; i<nAtoms; i++ ) {
        atoms[i].index = 0;
        order[i] = i;
        count[i] = 0;
      }
      count[0] = nAtoms;
    }

    void ActivatePartitions(unsigned int nAtoms,
                            int *order,
                            int *count,
                            int &activeset,int *next,
                            int *changed) {
      unsigned int i,j;
      activeset = -1;
      for( i=0; i<nAtoms; i++ )
        next[i] = -2;

      i = 0;
      do {
        j = order[i];
        if( count[j] > 1 ) {
          next[j] = activeset;
          activeset = j;
          i += count[j];
        } else i++;
      } while( i < nAtoms );

      for( i=0; i<nAtoms; i++ ){
        j = order[i];
        int flag=1;
#define SKIP_NODE_CHANGED_OPTIMIZATION 1
#ifndef SKIP_NODE_CHANGED_OPTIMIZATION
        if(count[j]){
          flag=(next[j]!=-2);
        }
#endif
        changed[j]=flag;
      }
    }

    template <typename T>
    void rankWithFunctor(T &ftor,
                         bool breakTies,
                         int *order){
      const ROMol &mol=*ftor.dp_mol;
      canon_atom *atoms=ftor.dp_atoms;
      unsigned int nAts=mol.getNumAtoms();
      int *count=(int *)malloc(nAts*sizeof(int));
      int activeset;
      int *next=(int *)malloc(nAts*sizeof(int));
      int *changed=(int *)malloc(nAts*sizeof(int));
      CreateSinglePartition(nAts,order,count,atoms);
      ActivatePartitions(nAts,order,count,activeset,next,changed);
      RefinePartitions(mol,atoms,ftor,false,order,count,activeset,next,changed);
      //#define VERBOSE_CANON 1
#ifdef VERBOSE_CANON
      std::cerr<<"1--------"<<std::endl;
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
      }
#endif
      ftor.df_useNbrs=true;
      ActivatePartitions(nAts,order,count,activeset,next,changed);
#ifdef VERBOSE_CANON
      std::cerr<<"1a--------"<<std::endl;
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
      }
#endif
      RefinePartitions(mol,atoms,ftor,true,order,count,activeset,next,changed);
#ifdef VERBOSE_CANON
      std::cerr<<"2--------"<<std::endl;
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
      }
#endif
      if(breakTies){
        BreakTies(mol,atoms,ftor,true,order,count,activeset,next,changed);
#ifdef VERBOSE_CANON
        std::cerr<<"3--------"<<std::endl;
        for(unsigned int i=0;i<mol.getNumAtoms();++i){
          std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
        }
#endif
      }
      free(count); free(next);
    }

    namespace {
      bool hasRingNbr(const ROMol &mol,const Atom *at) {
        ROMol::ADJ_ITER beg,end;
        boost::tie(beg,end) = mol.getAtomNeighbors(at);
        while(beg!=end){
          const ATOM_SPTR nbr=mol[*beg];
          ++beg;
          if((nbr->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW ||
              nbr->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW) &&
             nbr->hasProp("_ringStereoAtoms")){
            return true;
          }
        }
        return false;
      }
    }

    void rankMolAtoms(const ROMol &mol,std::vector<unsigned int> &res,
                      bool breakTies,
                      bool includeChirality,bool includeIsotopes) {
      std::vector<Canon::canon_atom> atoms(mol.getNumAtoms());
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        atoms[i].atom=mol.getAtomWithIdx(i);
        atoms[i].index=i;
        atoms[i].degree=mol.getAtomWithIdx(i)->getDegree();
        atoms[i].hasRingNbr=hasRingNbr(mol,atoms[i].atom);
        atoms[i].p_symbol=NULL;
      }
      AtomCompareFunctor ftor(&atoms.front(),mol);
      ftor.df_useIsotopes=includeIsotopes;
      ftor.df_useChirality=includeChirality;

      int *order=(int *)malloc(mol.getNumAtoms()*sizeof(int));
      rankWithFunctor(ftor,breakTies,order);

      res.resize(mol.getNumAtoms());
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        res[order[i]]=atoms[order[i]].index;
      }
      free(order); 
    } // end of rankMolAtoms()

    void rankFragmentAtoms(const ROMol &mol,std::vector<unsigned int> &res,
                           const boost::dynamic_bitset<> &atomsInPlay,
                           const boost::dynamic_bitset<> &bondsInPlay,
                           const std::vector<std::string> *atomSymbols,
                           bool breakTies,
                           bool includeChirality,bool includeIsotopes) {
      PRECONDITION(atomsInPlay.size()==mol.getNumAtoms(),"bad atomsInPlay size");
      PRECONDITION(bondsInPlay.size()==mol.getNumBonds(),"bad bondsInPlay size");
      PRECONDITION(!atomSymbols || atomSymbols->size()==mol.getNumAtoms(),"bad atomSymbols size");

      std::vector<Canon::canon_atom> atoms(mol.getNumAtoms());
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        atoms[i].atom=mol.getAtomWithIdx(i);
        atoms[i].index=i;
        atoms[i].degree=0;
        atoms[i].hasRingNbr=hasRingNbr(mol,atoms[i].atom);
        if(atomSymbols){
          atoms[i].p_symbol=&(*atomSymbols)[i];
        } else {
          atoms[i].p_symbol=0;
        }
      }
      for(ROMol::ConstBondIterator bI=mol.beginBonds();
          bI!=mol.endBonds();++bI){
        if(!bondsInPlay[(*bI)->getIdx()])
          continue;
        atoms[(*bI)->getBeginAtomIdx()].degree++;
        atoms[(*bI)->getEndAtomIdx()].degree++;
      }
      AtomCompareFunctor ftor(&atoms.front(),mol,
                              &atomsInPlay,&bondsInPlay);
      ftor.df_useIsotopes=includeIsotopes;
      ftor.df_useChirality=includeChirality;

      int *order=(int *)malloc(mol.getNumAtoms()*sizeof(int));
      rankWithFunctor(ftor,breakTies,order);

      res.resize(mol.getNumAtoms());
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        res[order[i]]=atoms[order[i]].index;
      }
      free(order); 
    } // end of rankFragmentAtoms()

    
    void chiralRankMolAtoms(const ROMol &mol,std::vector<unsigned int> &res){
      std::vector<Canon::canon_atom> atoms(mol.getNumAtoms());
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        atoms[i].atom=mol.getAtomWithIdx(i);
        atoms[i].index=i;
      }
      ChiralAtomCompareFunctor ftor(&atoms.front(),mol);

      int *order=(int *)malloc(mol.getNumAtoms()*sizeof(int));
      rankWithFunctor(ftor,false,order);

      res.resize(mol.getNumAtoms());
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        res[order[i]]=atoms[order[i]].index;
      }
      free(order); 
    } // end of rankMolAtoms()
  } // end of Canon namespace
} // end of RDKit namespace
