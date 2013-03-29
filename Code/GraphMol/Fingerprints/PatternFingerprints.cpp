// $Id$
//
//  Copyright (C) 2013 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include "Fingerprints.h"
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Invariant.h>
#include <boost/random.hpp>
#include <limits.h>
#include <boost/cstdint.hpp>
#include <RDGeneral/hash/hash.hpp>
#include <RDGeneral/types.h>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>

//#define VERBOSE_FINGERPRINTING 1
//#define REPORT_FP_STATS 1
#ifdef REPORT_FP_STATS
#include <GraphMol/SmilesParse/SmilesWrite.h>
#endif

namespace RDKit{
  const char *pqs[]={ "[*]~[*]",
                      "[*]~[*]~[*]",
                      "[R]~1~[R]~[R]~1",
                      //"[*]~[*]~[*]~[*]",
                      "[*]~[*](~[*])~[*]",
                      //"[*]~[R]~1[R]~[R]~1",
                      "[R]~1[R]~[R]~[R]~1",
                      //"[*]~[*]~[*]~[*]~[*]",
                      "[*]~[*]~[*](~[*])~[*]",
                      //"[*]~[R]~1[R]~[R]~1~[*]",
                      "[R]~1~[R]~[R]~[R]~[R]~1",
                      "[R]~1~[R]~[R]~[R]~[R]~[R]~1",
                      "[R2]~[R1]~[R2]",
                      "[R2]~[R1]~[R1]~[R2]",
                      "[*]!@[R]~[R]!@[*]",
                      "[*]!@[R]~[R]~[R]!@[*]",

#if 0
                      "[*]~[*](~[*])(~[*])~[*]",
                      "[*]~[*]~[*]~[*]~[*]~[*]",
                      "[*]~[*]~[*]~[*](~[*])~[*]",
                      "[*]~[*]~[*](~[*])~[*]~[*]",
                      "[*]~[*]~[*](~[*])(~[*])~[*]",
                      "[*]~[*](~[*])~[*](~[*])~[*]",
                      "[*]~[R]~1[R]~[R]~1(~[*])~[*]",
                      "[*]~[R]~1[R](~[*])~[R]~1[*]",
                      "[*]~[R]~1[R]~[R](~[*])~[R]~1",
                      "[*]~[R]~1[R]~[R]~[R]~1[*]",
                      "[*]~[R]~1[R]~[R]~[R]~[R]~1",
                      "[*]~[R]~1(~[*])~[R]~[R]~[R]~1",
                      "[*]~[*]~[*]~[*]~[*]~[*]~[*]",
                      "[*]~[*]~[*]~[*]~[*](~[*])~[*]",
                      "[*]~[*]~[*]~[*](~[*])~[*]~[*]",
                      "[*]~[*]~[*]~[*](~[*])(~[*])~[*]",
                      "[*]~[*]~[*](~[*])~[*](~[*])~[*]",
                      "[*]~[*](~[*])~[*]~[*](~[*])~[*]",
                      "[*]~[*](~[*])~[*](~[*])(~[*])~[*]",
#endif
                      ""};

  namespace detail {
    void getAtomNumbers(const Atom *a,std::vector<int> &atomNums){
      atomNums.clear();
      if( !a->hasQuery()){
        atomNums.push_back(a->getAtomicNum());
        return;
      }
      // negated things are always complex:
      if( a->getQuery()->getNegation()) return;
      std::string descr=a->getQuery()->getDescription();
      if(descr=="AtomAtomicNum"){
        atomNums.push_back(static_cast<ATOM_EQUALS_QUERY *>(a->getQuery())->getVal());
      } else if(descr=="AtomXor"){
        return;
      } else if(descr=="AtomAnd"){
        Queries::Query<int,Atom const *,true>::CHILD_VECT_CI childIt=a->getQuery()->beginChildren();
        if( (*childIt)->getDescription()=="AtomAtomicNum" &&
            ((*(childIt+1))->getDescription()=="AtomIsAliphatic" ||
             (*(childIt+1))->getDescription()=="AtomIsAromatic") &&
            (childIt+2)==a->getQuery()->endChildren()){
          atomNums.push_back(static_cast<ATOM_EQUALS_QUERY *>((*childIt).get())->getVal());
          return;
        }
      } else if(descr=="AtomOr"){
        Queries::Query<int,Atom const *,true>::CHILD_VECT_CI childIt=a->getQuery()->beginChildren();
        while(childIt !=a->getQuery()->endChildren()){
          if( (*childIt)->getDescription()=="AtomAtomicNum" ){
            atomNums.push_back(static_cast<ATOM_EQUALS_QUERY *>((*childIt).get())->getVal());
          } else if((*childIt)->getDescription()=="AtomAnd"){
            Queries::Query<int,Atom const *,true>::CHILD_VECT_CI childIt2=(*childIt)->beginChildren();
            if( (*childIt2)->getDescription()=="AtomAtomicNum" &&
                ((*(childIt2+1))->getDescription()=="AtomIsAliphatic" ||
                 (*(childIt2+1))->getDescription()=="AtomIsAromatic") &&
                (childIt2+2)==(*childIt)->endChildren()){
              atomNums.push_back(static_cast<ATOM_EQUALS_QUERY *>((*childIt2).get())->getVal());
            } else {
              atomNums.clear();
              return;
            }
          } else {
            atomNums.clear();
            return;
          }
          ++childIt;
        }
      }
      return;
    }
  }    


  // caller owns the result, it must be deleted
  ExplicitBitVect *PatternFingerprintMol(const ROMol &mol,
                                         unsigned int fpSize,
                                         std::vector<unsigned int> *atomCounts,
                                         ExplicitBitVect *setOnlyBits){
    PRECONDITION(fpSize!=0,"fpSize==0");
    PRECONDITION(!atomCounts || atomCounts->size()>=mol.getNumAtoms(),"bad atomCounts size");
    PRECONDITION(!setOnlyBits || setOnlyBits->getNumBits()==fpSize,"bad setOnlyBits size");

    static std::vector<ROMOL_SPTR> patts;
    // FIX: need a mutex here to be threadsafe
    if(patts.size()==0){
      unsigned int idx=0;
      while(1){
        std::string pq=pqs[idx];
        if(pq=="") break;
        idx++;
        RWMol *tm;
        try {
          tm = SmartsToMol(pq);
        }catch (...) {
          tm=NULL;
        }
        if(!tm) continue;
        patts.push_back(ROMOL_SPTR(static_cast<ROMol *>(tm)));
      }
    }
    if(!mol.getRingInfo()->isInitialized()){
      MolOps::findSSSR(mol);
    }

    boost::dynamic_bitset<> isQueryAtom(mol.getNumAtoms()),isQueryBond(mol.getNumBonds());
    ROMol::VERTEX_ITER firstA,lastA;
    boost::tie(firstA,lastA) = mol.getVertices();  
    while(firstA!=lastA){
      const Atom *at=mol[*firstA].get();
      if(Fingerprints::detail::isComplexQuery(at)) isQueryAtom.set(at->getIdx());
      ++firstA;
    }
    ROMol::EDGE_ITER firstB,lastB;
    boost::tie(firstB,lastB) = mol.getEdges();
    while(firstB!=lastB){
      const Bond *bond = mol[*firstB].get();
      if( Fingerprints::detail::isComplexQuery(bond) ){
        isQueryBond.set(bond->getIdx());
      }
      ++firstB;
    }
    
    ExplicitBitVect *res = new ExplicitBitVect(fpSize);
    unsigned int pIdx=0;
    BOOST_FOREACH(ROMOL_SPTR patt,patts){
      ++pIdx;
      std::vector<MatchVectType> matches;
      // uniquify matches?
      //   time for 10K molecules w/ uniquify: 5.24s
      //   time for 10K molecules w/o uniquify: 4.87s
      SubstructMatch(mol,*(patt.get()),matches,false); 
      boost::uint32_t mIdx=pIdx+patt->getNumAtoms()+patt->getNumBonds();
      BOOST_FOREACH(MatchVectType &mv,matches){
#ifdef VERBOSE_FINGERPRINTING
        std::cerr<<"\nPatt: "<<pIdx<<" | ";
#endif          
        // collect bits counting the number of occurances of the pattern:
        gboost::hash_combine(mIdx,0xBEEF);
        res->setBit(mIdx%fpSize);

        bool isQuery=false;
        boost::uint32_t bitId=pIdx;
        std::vector<unsigned int> amap(mv.size(),0);
        BOOST_FOREACH(MatchVectType::value_type &p,mv){
#ifdef VERBOSE_FINGERPRINTING
          std::cerr<<p.second<<" ";
#endif
          if(isQueryAtom[p.second]){
            isQuery=true;
#ifdef VERBOSE_FINGERPRINTING
            std::cerr<<"atom query.";
#endif
            break;
          }
          gboost::hash_combine(bitId,mol.getAtomWithIdx(p.second)->getAtomicNum());
          amap[p.first]=p.second;
        }
        if(isQuery) continue;
        ROMol::EDGE_ITER firstB,lastB;
        boost::tie(firstB,lastB) = patt->getEdges();
        while(firstB!=lastB){
          BOND_SPTR pbond = (*patt.get())[*firstB];
          ++firstB;
          if(isQueryBond[pbond->getIdx()]){
            isQuery=true;
#ifdef VERBOSE_FINGERPRINTING
            std::cerr<<"bond query: "<<pbond->getIdx();
#endif
            break;
          }
          const Bond *mbond=mol.getBondBetweenAtoms(amap[pbond->getBeginAtomIdx()],
                                                    amap[pbond->getEndAtomIdx()]);
          gboost::hash_combine(bitId,(boost::uint32_t)mbond->getBondType());
        }
        if(!isQuery){
#ifdef VERBOSE_FINGERPRINTING
          std::cerr<<" set: "<<bitId<<" "<<bitId%fpSize;
#endif
          res->setBit(bitId%fpSize);
        }
      }
    }
    return res;
  }
}
