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

#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>


//#define VERBOSE_FINGERPRINTING 1

  namespace {
    class ss_matcher {
    public:
      ss_matcher() {};
      ss_matcher(const std::string &pattern){
        RDKit::RWMol *p=RDKit::SmartsToMol(pattern);
        TEST_ASSERT(p);
        m_matcher.reset(p);
      };

      //const RDKit::ROMOL_SPTR &getMatcher() const { return m_matcher; };
      const RDKit::ROMol *getMatcher() const { return m_matcher.get(); };
    private:
      RDKit::ROMOL_SPTR m_matcher;
    };
  }

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
                      //"[R2]~[R1]~[R2]", Github #151: can't have ring counts in an SSS pattern
                      //"[R2]~[R1]~[R1]~[R2]",  Github #151: can't have ring counts in an SSS pattern
                      "[R](@[R])(@[R])~[R]~[R](@[R])(@[R])",
                      "[R](@[R])(@[R])~[R]@[R]~[R](@[R])(@[R])",
                      
                      //"[*]!@[R]~[R]!@[*]",  Github #151: can't have !@ in an SSS pattern
                      //"[*]!@[R]~[R]~[R]!@[*]", Github #151: can't have !@ in an SSS pattern
                      "[*]~[R](@[R])@[R](@[R])~[*]",
                      "[*]~[R](@[R])@[R]@[R](@[R])~[*]",
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
  typedef boost::flyweight<boost::flyweights::key_value<std::string,ss_matcher>,boost::flyweights::no_tracking > pattern_flyweight;

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

  namespace {
    bool isPatternComplexQuery(const Bond *b){
      if( !b->hasQuery()) return false;
      // negated things are always complex:
      if( b->getQuery()->getNegation()) return true;
      std::string descr=b->getQuery()->getDescription();
      //std::cerr<<"   !!!!!! "<<b->getIdx()<<" "<<b->getBeginAtomIdx()<<"-"<<b->getEndAtomIdx()<<" "<<descr<<std::endl;
      if(descr=="BondOrder") return false;
      return true;
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

    std::vector<const ROMol *> patts;
    patts.reserve(10);
    unsigned int idx=0;
    while(1){
      std::string pq=pqs[idx];
      if(pq=="") break;
      ++idx;
      const ROMol *matcher=pattern_flyweight(pq).get().getMatcher();
      CHECK_INVARIANT(matcher,"bad smarts");
      patts.push_back(matcher);
    }

    if(!mol.getRingInfo()->isInitialized()){
      MolOps::fastFindRings(mol);
    }

    boost::dynamic_bitset<> isQueryAtom(mol.getNumAtoms()),isQueryBond(mol.getNumBonds());
    ROMol::VERTEX_ITER firstA,lastA;
    boost::tie(firstA,lastA) = mol.getVertices();  
    while(firstA!=lastA){
      const Atom *at=mol[*firstA].get();
      if(Fingerprints::detail::isComplexQuery(at)){
        isQueryAtom.set(at->getIdx());
        //std::cerr<<"   complex atom: "<<at->getIdx()<<std::endl;
      }
      ++firstA;
    }
    ROMol::EDGE_ITER firstB,lastB;
    boost::tie(firstB,lastB) = mol.getEdges();
    while(firstB!=lastB){
      const Bond *bond = mol[*firstB].get();
      //if( Fingerprints::detail::isComplexQuery(bond) ){
      if( isPatternComplexQuery(bond) ){
        isQueryBond.set(bond->getIdx());
        //std::cerr<<"   complex bond: "<<bond->getIdx()<<std::endl;
      }
      ++firstB;
    }
    
    ExplicitBitVect *res = new ExplicitBitVect(fpSize);
    unsigned int pIdx=0;
    BOOST_FOREACH(const ROMol *patt,patts){
      ++pIdx;
      std::vector<MatchVectType> matches;
      // uniquify matches?
      //   time for 10K molecules w/ uniquify: 5.24s
      //   time for 10K molecules w/o uniquify: 4.87s
      SubstructMatch(mol,*patt,matches,false); 
      boost::uint32_t mIdx=pIdx+patt->getNumAtoms()+patt->getNumBonds();
      BOOST_FOREACH(MatchVectType &mv,matches){
#ifdef VERBOSE_FINGERPRINTING
        std::cerr<<"\nPatt: "<<pIdx<<" | ";
#endif          
        // collect bits counting the number of occurances of the pattern:
        gboost::hash_combine(mIdx,0xBEEF);
        res->setBit(mIdx%fpSize);
#ifdef VERBOSE_FINGERPRINTING
        std::cerr<<"count: "<<mIdx%fpSize<<" | ";
#endif          

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
#ifdef VERBOSE_FINGERPRINTING
        std::cerr<<" bs:|| ";
#endif
        while(!isQuery && firstB!=lastB){
          BOND_SPTR pbond = (*patt)[*firstB];
          ++firstB;
          const Bond *mbond=mol.getBondBetweenAtoms(amap[pbond->getBeginAtomIdx()],
                                                    amap[pbond->getEndAtomIdx()]);

          if(isQueryBond[mbond->getIdx()]){
            isQuery=true;
#ifdef VERBOSE_FINGERPRINTING
            std::cerr<<"bond query: "<<mbond->getIdx();
#endif
            break;
          }
          // makes sure aromatic bonds and single bonds from SMARTS always hash the same:
          //if(!mbond->getIsAromatic() && mbond->getBondType()!=Bond::SINGLE &&
          //   mbond->getBondType()!=Bond::AROMATIC){
          if(!mbond->getIsAromatic()){
            gboost::hash_combine(bitId,(boost::uint32_t)mbond->getBondType());
#ifdef VERBOSE_FINGERPRINTING
            std::cerr<<mbond->getBondType()<<" ";
#endif
          } else {
            gboost::hash_combine(bitId,(boost::uint32_t)Bond::AROMATIC);
#ifdef VERBOSE_FINGERPRINTING
            std::cerr<<Bond::AROMATIC<<" ";
#endif
          }
            //} else {
            //  gboost::hash_combine(bitId,(boost::uint32_t)Bond::SINGLE);
            //          }

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
