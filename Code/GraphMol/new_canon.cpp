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
//#define SKIP_NODE_CHANGED_OPTIMIZATION 0
//#ifndef SKIP_NODE_CHANGED_OPTIMIZATION
//        if(count[j]){
//          std::cout << "j " << j << std::endl;
//          flag=(next[j]!=-2);
//        }
//#endif
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
      char *touched=(char *)malloc(nAts*sizeof(char));
      memset(touched,0,nAts*sizeof(char));
      memset(changed,1,nAts*sizeof(int));
      CreateSinglePartition(nAts,order,count,atoms);
      ActivatePartitions(nAts,order,count,activeset,next,changed);
      RefinePartitions(mol,atoms,ftor,false,order,count,activeset,next,changed,touched);
      //#define VERBOSE_CANON 1
#ifdef VERBOSE_CANON
      std::cerr<<"1--------"<<std::endl;
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        std::cerr<<order[i]+1<<" "<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
      }
#endif
      ftor.df_useNbrs=true;
      ActivatePartitions(nAts,order,count,activeset,next,changed);
#ifdef VERBOSE_CANON
      std::cerr<<"1a--------"<<std::endl;
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        std::cerr<<order[i]+1<<" "<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
      }
#endif
      RefinePartitions(mol,atoms,ftor,true,order,count,activeset,next,changed,touched);
#ifdef VERBOSE_CANON
      std::cerr<<"2--------"<<std::endl;
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        std::cerr<<order[i]+1<<" "<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
      }
#endif
      if(breakTies){
        BreakTies(mol,atoms,ftor,true,order,count,activeset,next,changed,touched);
#ifdef VERBOSE_CANON
        std::cerr<<"3--------"<<std::endl;
        for(unsigned int i=0;i<mol.getNumAtoms();++i){
          std::cerr<<order[i]+1<<" "<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
        }
#endif
      }
      free(count); free(next); free(touched); free(changed);
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

      unsigned int getNumRingMember(const ROMol &mol,unsigned int idx){
        if(mol.getRingInfo()->isInitialized()){
          return mol.getRingInfo()->numAtomRings(idx);
        }
        return 0;
      }

      void getNbrs(const ROMol &mol,const Atom *at, int *ids){
        ROMol::ADJ_ITER beg,end;
        boost::tie(beg,end) = mol.getAtomNeighbors(at);
        unsigned int idx=0;
        while(beg!=end){
          const ATOM_SPTR nbr=(mol)[*beg];
          ++beg;
          ids[idx] = nbr->getIdx();
          ++idx;
        }
      }

      void getBonds(const ROMol &mol,const Atom *at,bool includeChirality,std::vector<bondholder> &nbrs){

        ROMol::OEDGE_ITER beg,end;
        boost::tie(beg,end) = mol.getAtomBonds(at);
        while(beg!=end){
          const BOND_SPTR bond=(mol)[*beg];
          ++beg;
          Bond::BondStereo stereo=Bond::STEREONONE;
          if(includeChirality){
            stereo=bond->getStereo();
          }
          unsigned int idx = bond->getOtherAtomIdx(at->getIdx());
          bondholder bh(bondholder(bond->getBondType(),stereo,idx,idx));
          nbrs.insert(std::lower_bound(nbrs.begin(),nbrs.end(),bh),1,bh);
        }
        std::reverse(nbrs.begin(),nbrs.end());
      }

      void getChiralBonds(const ROMol &mol,const Atom *at,std::vector<bondholder> &nbrs){

        ROMol::OEDGE_ITER beg,end;
        boost::tie(beg,end) = mol.getAtomBonds(at);
        while(beg!=end){
          const BOND_SPTR bond=(mol)[*beg];
          ++beg;
          unsigned int nbrIdx = bond->getOtherAtomIdx(at->getIdx());
          const Atom* nbr = mol.getAtomWithIdx(nbrIdx);
          unsigned int degreeNbr = nbr->getDegree();
          unsigned int nReps=1;
          unsigned int stereo=0;
          switch(bond->getStereo()){
          case Bond::STEREOZ:
            stereo=1;
            break;
          case Bond::STEREOE:
            stereo=2;
            break;
          default:
            stereo=0;
          }
          if(bond->getBondType() == Bond::DOUBLE &&
              nbr->getAtomicNum()==15 &&
              (degreeNbr==4 || degreeNbr==3) ) {
            // a special case for chiral phophorous compounds
            // (this was leading to incorrect assignment of
            // R/S labels ):
            nReps=1;
            // general justification of this is:
            // Paragraph 2.2. in the 1966 article is "Valence-Bond Conventions:
            // Multiple-Bond Unsaturation and Aromaticity". It contains several
            // conventions of which convention (b) is the one applying here:
            // "(b) Contibutions by d orbitals to bonds of quadriligant atoms are
            // neglected."
            // FIX: this applies to more than just P
          } else {
            nReps = static_cast<unsigned int>(floor(2.*bond->getBondTypeAsDouble()));
          }
          unsigned int symclass = nbr->getAtomicNum()*ATNUM_CLASS_OFFSET+nbrIdx+1;
          bondholder bh(bondholder(Bond::SINGLE,stereo,nbrIdx,symclass));
          std::vector<bondholder>::iterator iPos=std::lower_bound(nbrs.begin(),nbrs.end(),bh);
          nbrs.insert(iPos,nReps,bh);
        }
        std::reverse(nbrs.begin(),nbrs.end());

        if(!at->needsUpdatePropertyCache()){
          for(unsigned int ii=0;ii<at->getTotalNumHs();++ii){
            nbrs.push_back(bondholder(Bond::SINGLE,Bond::STEREONONE,ATNUM_CLASS_OFFSET,ATNUM_CLASS_OFFSET));
            nbrs.push_back(bondholder(Bond::SINGLE,Bond::STEREONONE,ATNUM_CLASS_OFFSET,ATNUM_CLASS_OFFSET));
          }
        }
      }

    }

    void initCanonAtoms(const ROMol &mol,std::vector<Canon::canon_atom> &atoms,
        bool includeChirality){
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        atoms[i].atom=mol.getAtomWithIdx(i);
        atoms[i].index=i;
        atoms[i].p_symbol=NULL;
        atoms[i].degree=atoms[i].atom->getDegree();
        atoms[i].totalNumHs=atoms[i].atom->getTotalNumHs();
        atoms[i].numRingMember=getNumRingMember(mol,i);
        atoms[i].isRingStereoAtom=(atoms[i].atom->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW ||
            atoms[i].atom->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW) &&
            atoms[i].atom->hasProp("_ringStereoAtoms");
        atoms[i].nbrIds=(int *)malloc(atoms[i].degree*sizeof(int));
        getNbrs(mol, atoms[i].atom,atoms[i].nbrIds);
        atoms[i].bonds.reserve(atoms[i].degree);
        getBonds(mol,atoms[i].atom,includeChirality,atoms[i].bonds);
        /* this could be realized using the neighbors above */
        atoms[i].hasRingNbr=hasRingNbr(mol,atoms[i].atom);
      }
    }

    void initFragmentCanonAtoms(const ROMol &mol,std::vector<Canon::canon_atom> &atoms,
        bool includeChirality, const std::vector<std::string> *atomSymbols){
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        atoms[i].atom=mol.getAtomWithIdx(i);
        atoms[i].index=i;
        if(atomSymbols){
          atoms[i].p_symbol=&(*atomSymbols)[i];
        }
        else{
          atoms[i].p_symbol=0;
        }
        atoms[i].degree=0;
        atoms[i].totalNumHs=atoms[i].atom->getTotalNumHs();
        atoms[i].numRingMember=getNumRingMember(mol,i);
        atoms[i].isRingStereoAtom=(atoms[i].atom->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW ||
            atoms[i].atom->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW) &&
            atoms[i].atom->hasProp("_ringStereoAtoms");
        atoms[i].nbrIds=(int *)malloc(atoms[i].degree*sizeof(int));
        getNbrs(mol, atoms[i].atom,atoms[i].nbrIds);
        atoms[i].bonds.reserve(atoms[i].degree);
        getBonds(mol,atoms[i].atom,includeChirality,atoms[i].bonds);
        /* this could be realized using the neighbors above */
        atoms[i].hasRingNbr=hasRingNbr(mol,atoms[i].atom);
      }
    }

    void initChiralCanonAtoms(const ROMol &mol,std::vector<Canon::canon_atom> &atoms){
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        atoms[i].atom=mol.getAtomWithIdx(i);
        atoms[i].index=i;
        atoms[i].degree=atoms[i].atom->getDegree();
        atoms[i].nbrIds=(int *)malloc(atoms[i].degree*sizeof(int));
        getNbrs(mol, atoms[i].atom,atoms[i].nbrIds);
        getChiralBonds(mol,atoms[i].atom,atoms[i].bonds);
      }
    }

    void rankMolAtoms(const ROMol &mol,std::vector<unsigned int> &res,
                      bool breakTies,
                      bool includeChirality,bool includeIsotopes) {
      std::vector<Canon::canon_atom> atoms(mol.getNumAtoms());
      initCanonAtoms(mol,atoms,includeChirality);
      AtomCompareFunctor ftor(&atoms.front(),mol);
      ftor.df_useIsotopes=includeIsotopes;
      ftor.df_useChirality=includeChirality;
      ftor.df_useChiralityRings=includeChirality;

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
      initFragmentCanonAtoms(mol,atoms,includeChirality, atomSymbols);
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
      initChiralCanonAtoms(mol,atoms);
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
