// $Id$
//
//  Copyright (c) 2009, Novartis Institutes for BioMedical Research Inc.
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
#include <boost/functional/hash.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/foreach.hpp>
#include <algorithm>

namespace RDKit{
  namespace MorganFingerprints {
    using boost::uint32_t;
    using boost::int32_t;
    typedef boost::tuple<boost::dynamic_bitset<>,uint32_t,unsigned int> AccumTuple;

    void getConnectivityInvariants(const ROMol &mol,
                                   std::vector<uint32_t> &invars,
                                   bool includeRingMembership){
      unsigned int nAtoms=mol.getNumAtoms();
      PRECONDITION(invars.size()>=nAtoms,"vector too small");
      boost::hash<std::vector<uint32_t> > vectHasher;
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
      
    SparseIntVect<uint32_t> *
    getFingerprint(const ROMol &mol,
                   unsigned int radius,
                   std::vector<uint32_t> *invariants,
                   const std::vector<uint32_t> *fromAtoms){
      unsigned int nAtoms=mol.getNumAtoms();
      bool owner=false;
      if(!invariants){
        invariants = new std::vector<uint32_t>(nAtoms);
        owner=true;
        getConnectivityInvariants(mol,*invariants);
      }
      SparseIntVect<uint32_t> *res;
      res = new SparseIntVect<uint32_t>(std::numeric_limits<uint32_t>::max());

      // add the round 0 invariants to the result:
      for(unsigned int i=0;i<nAtoms;++i){
        if(!fromAtoms ||
           std::find(fromAtoms->begin(),fromAtoms->end(),i)!=fromAtoms->end()){
          res->setVal((*invariants)[i],res->getVal((*invariants)[i])+1);
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

      // now do our subsequent rounds:
      for(unsigned int layer=0;layer<radius;++layer){
        std::vector<uint32_t> roundInvariants(nAtoms);
        std::vector< boost::dynamic_bitset<> > roundAtomNeighborhoods=atomNeighborhoods;
        std::vector< AccumTuple > neighborhoodsThisRound;
          
        for(unsigned int atomIdx=0;atomIdx<nAtoms;++atomIdx){
          if(!deadAtoms[atomIdx]){
            std::vector< std::pair<int32_t,uint32_t> > nbrs;
            ROMol::OEDGE_ITER beg,end;
            boost::tie(beg,end) = mol.getAtomBonds(mol.getAtomWithIdx(atomIdx));
            while(beg!=end){
              const BOND_SPTR bond=mol[*beg];
              roundAtomNeighborhoods[atomIdx][bond->getIdx()]=1;

              unsigned int oIdx=bond->getOtherAtomIdx(atomIdx);
              roundAtomNeighborhoods[atomIdx] |= atomNeighborhoods[oIdx];

              nbrs.push_back(std::make_pair(static_cast<int32_t>(bond->getBondType()),
                                            (*invariants)[oIdx]));

              ++beg;
            }

            // sort the neighbor list:
            std::sort(nbrs.begin(),nbrs.end());
            // and now calculate the new invariant and test if the atom is newly
            // "chiral"
            std::size_t invar=layer;
            boost::hash_combine(invar,(*invariants)[atomIdx]);
            bool looksChiral = (mol.getAtomWithIdx(atomIdx)->getChiralTag()!=Atom::CHI_UNSPECIFIED);
            for(std::vector< std::pair<int32_t,uint32_t> >::const_iterator it=nbrs.begin();
                it!=nbrs.end();++it){
              // add the contribution to the new invariant:
              boost::hash_combine(invar, *it);

              //std::cerr<<"     "<<atomIdx<<": "<<it->first<<" "<<it->second<<" -> "<<invar<<std::endl;
                
              // update our "chirality":
              if(looksChiral && chiralAtoms[atomIdx]){
                if(it->first != static_cast<int32_t>(Bond::SINGLE)){
                  looksChiral=false;
                } else if(it!=nbrs.begin() && it->second == (it-1)->second) {
                  looksChiral=false;
                }
              }

            }
            if(looksChiral){
              chiralAtoms[atomIdx]=1;
              // add an extra value to the invariant to reflect chirality:
              boost::hash_combine(invar, 1);
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
            }
          }
        }
        std::sort(neighborhoodsThisRound.begin(),neighborhoodsThisRound.end());
        for(std::vector< AccumTuple >::const_iterator iter=neighborhoodsThisRound.begin();
            iter!=neighborhoodsThisRound.end();++iter){
          // if we haven't seen this exact environment before, update the fingerprint:
          if(std::find(neighborhoods.begin(),neighborhoods.end(),
                       iter->get<0>())==neighborhoods.end()){
            if(includeAtoms[iter->get<2>()]){
              res->setVal(iter->get<1>(),res->getVal(iter->get<1>())+1);
            }
            neighborhoods.push_back(iter->get<0>());
            //std::cerr<<" layer: "<<layer<<" atom: "<<iter->get<2>()<<" " <<iter->get<0>()<< " " << iter->get<1>() << " " << deadAtoms[iter->get<2>()]<<std::endl;
          } else {
            // we have seen this exact environment before, this atom
            // is now out of consideration:
            deadAtoms[iter->get<2>()]=1;
          }
        }

        // the invariants from this round become the global invariants:
        std::copy(roundInvariants.begin(),roundInvariants.end(),invariants->begin());

        atomNeighborhoods=roundAtomNeighborhoods;
      }

      if(owner) delete invariants;
      return res;
    }
  } // end of namespace MorganFingerprints
} // end of namespace RDKit
