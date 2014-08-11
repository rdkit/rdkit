// $Id$
//
//  Copyright (c) 2009-2010, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met: 
//
//     * Redistributions of source code must retain the above copyright 
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following 
//       disclaimer in the documentation and/or other materials provided 
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
//       nor the names of its contributors may be used to endorse or promote 
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  Created by Greg Landrum, July 2008
//
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <RDGeneral/hash/hash.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>


#include <boost/dynamic_bitset.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/foreach.hpp>
#include <algorithm>

#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>

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
  namespace MorganFingerprints {
    using boost::uint32_t;
    using boost::int32_t;
    typedef boost::tuple<boost::dynamic_bitset<>,uint32_t,unsigned int> AccumTuple;

    // Definitions for feature points adapted from:
    // Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998)
    const char *smartsPatterns[6]={
      "[$([N;!H0;v3,v4&+1]),\
$([O,S;H1;+0]),\
n&H1&+0]", // Donor
      "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),\
$([O,S;H0;v2]),\
$([O,S;-]),\
$([N;v3;!$(N-*=[O,N,P,S])]),\
n&H0&+0,\
$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]", // Acceptor
      "[a]", //Aromatic
      "[F,Cl,Br,I]",//Halogen
      "[#7;+,\
$([N;H2&+0][$([C,a]);!$([C,a](=O))]),\
$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),\
$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]", // Basic
      "[$([C,S](=[O,S,P])-[O;H1,-1])]" //Acidic
    };
    std::vector<std::string> defaultFeatureSmarts(smartsPatterns,smartsPatterns+6);
    typedef boost::flyweight<boost::flyweights::key_value<std::string,ss_matcher>,boost::flyweights::no_tracking > pattern_flyweight;
    void getFeatureInvariants(const ROMol &mol,
                              std::vector<uint32_t> &invars,
                              std::vector<const ROMol *> *patterns){
      unsigned int nAtoms=mol.getNumAtoms();
      PRECONDITION(invars.size()>=nAtoms,"vector too small");

      std::vector<const ROMol *> featureMatchers;
      if(!patterns){
        featureMatchers.reserve(defaultFeatureSmarts.size());
        for(std::vector<std::string>::const_iterator smaIt=defaultFeatureSmarts.begin();
            smaIt!=defaultFeatureSmarts.end();++smaIt){
          const ROMol *matcher=pattern_flyweight(*smaIt).get().getMatcher();
          CHECK_INVARIANT(matcher,"bad smarts");
          featureMatchers.push_back(matcher);
        }
        patterns=&featureMatchers;
      }
      std::fill(invars.begin(),invars.end(),0);
      for(unsigned int i=0;i<patterns->size();++i){
        unsigned int mask=1<<i;
        std::vector<MatchVectType> matchVect;
        // to maintain thread safety, we have to copy the pattern
        // molecules:
        SubstructMatch(mol,ROMol(*(*patterns)[i],true),matchVect);
        for(std::vector<MatchVectType>::const_iterator mvIt=matchVect.begin();
            mvIt!=matchVect.end();++mvIt){
          for(MatchVectType::const_iterator mIt=mvIt->begin();
              mIt!=mvIt->end();++mIt){
            invars[mIt->second]|=mask;
          }
        }
      }
    } // end of getFeatureInvariants()

    void getConnectivityInvariants(const ROMol &mol,
                                   std::vector<uint32_t> &invars,
                                   bool includeRingMembership){
      unsigned int nAtoms=mol.getNumAtoms();
      PRECONDITION(invars.size()>=nAtoms,"vector too small");
      gboost::hash<std::vector<uint32_t> > vectHasher;
      for(unsigned int i=0;i<nAtoms;++i){
        Atom const *atom = mol.getAtomWithIdx(i);
        std::vector<uint32_t> components;
        components.push_back(atom->getAtomicNum());
        components.push_back(atom->getTotalDegree());
        components.push_back(atom->getTotalNumHs());
        components.push_back(atom->getFormalCharge());
        int deltaMass = static_cast<int>(atom->getMass() -
                                         PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
        components.push_back(deltaMass);

        if(includeRingMembership && 
           atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx())){
          components.push_back(1);
        }
        invars[i]=vectHasher(components);
      }
    } // end of getConnectivityInvariants()

    uint32_t updateElement(SparseIntVect<uint32_t> &v,unsigned int elem, bool counting){
      uint32_t bit=elem%v.getLength();
      if (counting) {
          v.setVal(bit,v.getVal(bit)+1);
      } else {
          v.setVal(bit,1);
      }
      return bit;
    }
    uint32_t updateElement(ExplicitBitVect &v,unsigned int elem, bool counting=false){
      uint32_t bit=elem%v.getNumBits();
      v.setBit(bit);
      return bit;
    }

    template <typename T>
    void calcFingerprint(const ROMol &mol,
                         unsigned int radius,
                         std::vector<uint32_t> *invariants,
                         const std::vector<uint32_t> *fromAtoms,
                         bool useChirality,
                         bool useBondTypes,
                         bool useCounts,
                         bool onlyNonzeroInvariants,
                         BitInfoMap *atomsSettingBits,
                         T &res){
      unsigned int nAtoms=mol.getNumAtoms();
      bool owner=false;
      if(!invariants){
        invariants = new std::vector<uint32_t>(nAtoms);
        owner=true;
        getConnectivityInvariants(mol,*invariants);
      }
      // Make a copy of the invariants:
      std::vector<uint32_t> invariantCpy(nAtoms);
      std::copy(invariants->begin(),invariants->end(),invariantCpy.begin());

      // add the round 0 invariants to the result:
      for(unsigned int i=0;i<nAtoms;++i){
        if(!fromAtoms ||
           std::find(fromAtoms->begin(),fromAtoms->end(),i)!=fromAtoms->end()){
          if(!onlyNonzeroInvariants || (*invariants)[i]){
            uint32_t bit=updateElement(res,(*invariants)[i], useCounts);
            if(atomsSettingBits) (*atomsSettingBits)[bit].push_back(std::make_pair(i,0));
          }
        }
      }

      boost::dynamic_bitset<> chiralAtoms(nAtoms);
        
      // these are the neighborhoods that have already been added
      // to the fingerprint
      std::vector< boost::dynamic_bitset<> > neighborhoods;
      // these are the environments around each atom:
      std::vector< boost::dynamic_bitset<> > atomNeighborhoods(nAtoms,
                                                               boost::dynamic_bitset<>(mol.getNumBonds()));
      boost::dynamic_bitset<> deadAtoms(nAtoms);

      boost::dynamic_bitset<> includeAtoms(nAtoms);
      if(fromAtoms){
        BOOST_FOREACH(uint32_t idx,*fromAtoms){
          includeAtoms.set(idx,1);
        }
      } else {
        includeAtoms.set();
      }

      std::vector<unsigned int> atomOrder(nAtoms);
      if( onlyNonzeroInvariants ){
        std::vector< std::pair<int32_t,uint32_t> > ordering;
        for(unsigned int i=0;i<nAtoms;++i){
          if(!(*invariants)[i])
            ordering.push_back(std::make_pair(1,i));
          else 
            ordering.push_back(std::make_pair(0,i));
        }
        std::sort(ordering.begin(),ordering.end());
        for(unsigned int i=0;i<nAtoms;++i){
          atomOrder[i]=ordering[i].second;
        }
      } else {
        for(unsigned int i=0;i<nAtoms;++i){
          atomOrder[i]=i;
        }
      }
      // now do our subsequent rounds:
      for(unsigned int layer=0;layer<radius;++layer){
        std::vector<uint32_t> roundInvariants(nAtoms);
        std::vector< boost::dynamic_bitset<> > roundAtomNeighborhoods=atomNeighborhoods;
        std::vector< AccumTuple > neighborhoodsThisRound;
          
        BOOST_FOREACH(unsigned int atomIdx,atomOrder){
          if(!deadAtoms[atomIdx]){
            std::vector< std::pair<int32_t,uint32_t> > nbrs;
            ROMol::OEDGE_ITER beg,end;
            boost::tie(beg,end) = mol.getAtomBonds(mol.getAtomWithIdx(atomIdx));
            while(beg!=end){
              const BOND_SPTR bond=mol[*beg];
              roundAtomNeighborhoods[atomIdx][bond->getIdx()]=1;

              unsigned int oIdx=bond->getOtherAtomIdx(atomIdx);
              roundAtomNeighborhoods[atomIdx] |= atomNeighborhoods[oIdx];

              if(useBondTypes){
                nbrs.push_back(std::make_pair(static_cast<int32_t>(bond->getBondType()),
                                              (*invariants)[oIdx]));
              } else {
                nbrs.push_back(std::make_pair(static_cast<int32_t>(1),
                                              (*invariants)[oIdx]));
              }

              ++beg;
            }

            // sort the neighbor list:
            std::sort(nbrs.begin(),nbrs.end());
            // and now calculate the new invariant and test if the atom is newly
            // "chiral"
            boost::uint32_t invar=layer;
            gboost::hash_combine(invar,(*invariants)[atomIdx]);
            bool looksChiral = (mol.getAtomWithIdx(atomIdx)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
            for(std::vector< std::pair<int32_t,uint32_t> >::const_iterator it=nbrs.begin();
                it!=nbrs.end();++it){
              // add the contribution to the new invariant:
              gboost::hash_combine(invar, *it);

              //std::cerr<<"     "<<atomIdx<<": "<<it->first<<" "<<it->second<<" -> "<<invar<<std::endl;
                
              // update our "chirality":
              if(useChirality && looksChiral && chiralAtoms[atomIdx]){
                if(it->first != static_cast<int32_t>(Bond::SINGLE)){
                  looksChiral=false;
                } else if(it!=nbrs.begin() && it->second == (it-1)->second) {
                  looksChiral=false;
                }
              }

            }
            if(useChirality && looksChiral){
              chiralAtoms[atomIdx]=1;
              // add an extra value to the invariant to reflect chirality:
              Atom const *tAt=mol.getAtomWithIdx(atomIdx);
              std::string cip="";
              if(tAt->hasProp("_CIPCode")){
                tAt->getProp("_CIPCode",cip);
              }
              if(cip=="R"){
                gboost::hash_combine(invar, 3);
              } else if(cip=="S"){
                gboost::hash_combine(invar, 2);
              } else {
                gboost::hash_combine(invar, 1);
              }
            }
            roundInvariants[atomIdx]=static_cast<uint32_t>(invar);
            neighborhoodsThisRound.push_back(boost::make_tuple(roundAtomNeighborhoods[atomIdx],
                                                               static_cast<uint32_t>(invar),
                                                               atomIdx));
            if(std::find(neighborhoods.begin(),neighborhoods.end(),
                         roundAtomNeighborhoods[atomIdx])!=neighborhoods.end()){
              // we have seen this exact environment before, this atom
              // is now out of consideration:
              deadAtoms[atomIdx]=1;
              //std::cerr<<"   atom: "<< atomIdx <<" is dead."<<std::endl;
            }
          }
        }
        std::sort(neighborhoodsThisRound.begin(),neighborhoodsThisRound.end());
        for(std::vector< AccumTuple >::const_iterator iter=neighborhoodsThisRound.begin();
            iter!=neighborhoodsThisRound.end();++iter){
          // if we haven't seen this exact environment before, update the fingerprint:
          if(std::find(neighborhoods.begin(),neighborhoods.end(),
                       iter->get<0>())==neighborhoods.end()){
            if(!onlyNonzeroInvariants || invariantCpy[iter->get<2>()]){
              if(includeAtoms[iter->get<2>()]){
                uint32_t bit=updateElement(res,iter->get<1>(), useCounts);
                if(atomsSettingBits) (*atomsSettingBits)[bit].push_back(std::make_pair(iter->get<2>(),
                                                                                       layer+1));
              }
              if(!fromAtoms || std::find(fromAtoms->begin(),fromAtoms->end(),
                                         iter->get<2>())!=fromAtoms->end()){
                neighborhoods.push_back(iter->get<0>());
              }
            }
            //std::cerr<<" layer: "<<layer<<" atom: "<<iter->get<2>()<<" " <<iter->get<0>()<< " " << iter->get<1>() << " " << deadAtoms[iter->get<2>()]<<std::endl;
          } else {
            // we have seen this exact environment before, this atom
            // is now out of consideration:
            //std::cerr<<"   atom: "<< iter->get<2>()<<" is dead."<<std::endl;
            deadAtoms[iter->get<2>()]=1;
          }
        }

        // the invariants from this round become the global invariants:
        std::copy(roundInvariants.begin(),roundInvariants.end(),invariants->begin());

        atomNeighborhoods=roundAtomNeighborhoods;
      }

      if(owner) delete invariants;
    }
      
    SparseIntVect<uint32_t> *
    getFingerprint(const ROMol &mol,
                   unsigned int radius,
                   std::vector<uint32_t> *invariants,
                   const std::vector<uint32_t> *fromAtoms,
                   bool useChirality,bool useBondTypes,
                   bool useCounts,
                   bool onlyNonzeroInvariants,
                   BitInfoMap *atomsSettingBits){
      SparseIntVect<uint32_t> *res;
      res = new SparseIntVect<uint32_t>(std::numeric_limits<uint32_t>::max());
      calcFingerprint(mol,radius,invariants,fromAtoms,useChirality,useBondTypes,useCounts,
                      onlyNonzeroInvariants,atomsSettingBits,*res);
      return res;
    }
    SparseIntVect<uint32_t> *
    getHashedFingerprint(const ROMol &mol,
                         unsigned int radius,
                         unsigned int nBits,
                         std::vector<uint32_t> *invariants,
                         const std::vector<uint32_t> *fromAtoms,
                         bool useChirality,bool useBondTypes,
                         bool onlyNonzeroInvariants,
                         BitInfoMap *atomsSettingBits){
      SparseIntVect<uint32_t> *res;
      res = new SparseIntVect<uint32_t>(nBits);
      calcFingerprint(mol,radius,invariants,fromAtoms,useChirality,useBondTypes,true,
                      onlyNonzeroInvariants,atomsSettingBits,*res);
      return res;
    }

    ExplicitBitVect *
    getFingerprintAsBitVect(const ROMol &mol,
                            unsigned int radius,
                            unsigned int nBits,
                            std::vector<uint32_t> *invariants,
                            const std::vector<uint32_t> *fromAtoms,
                            bool useChirality,bool useBondTypes,
                            bool onlyNonzeroInvariants,
                            BitInfoMap *atomsSettingBits){
      ExplicitBitVect *res=new ExplicitBitVect(nBits);
      calcFingerprint(mol,radius,invariants,fromAtoms,useChirality,useBondTypes,false,
                      onlyNonzeroInvariants,atomsSettingBits,*res);
      return res;
    }


  } // end of namespace MorganFingerprints
} // end of namespace RDKit
