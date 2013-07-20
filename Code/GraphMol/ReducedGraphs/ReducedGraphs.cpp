//
//  Copyright (C) 2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "ReducedGraphs.h"

#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>


namespace RDKit{

  namespace {
    // FIX: this is duplicated here and in the MorganFingerprints code
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


  namespace ReducedGraphs {
    // Definitions for feature points adapted from:
    // Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998)
    const char *smartsPatterns[4]={
      "[$([N;!H0;v3,v4&+1]),\
$([O,S;H1;+0]),\
n&H1&+0]", // Donor
      "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),\
$([O,S;H0;v2]),\
$([O,S;-]),\
$([N;v3;!$(N-*=[O,N,P,S])]),\
n&H0&+0,\
$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]", // Acceptor
      "[#7;+,\
$([N;H2&+0][$([C,a]);!$([C,a](=O))]),\
$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),\
$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]", // Positive
      "[$([C,S](=[O,S,P])-[O;H1,-1])]" // Negative
    };
    std::vector<std::string> defaultFeatureSmarts(smartsPatterns,smartsPatterns+4);
    typedef boost::flyweight<boost::flyweights::key_value<std::string,ss_matcher>,boost::flyweights::no_tracking > pattern_flyweight;
    
    void getErGAtomTypes(const ROMol &mol,
                      std::vector<boost::dynamic_bitset<> > &types,
                      std::vector<const ROMol *> *patterns=0){
      unsigned int nAtoms=mol.getNumAtoms();

      if(!patterns){
        std::vector<const ROMol *> featureMatchers;
        featureMatchers.reserve(defaultFeatureSmarts.size());
        for(std::vector<std::string>::const_iterator smaIt=defaultFeatureSmarts.begin();
            smaIt!=defaultFeatureSmarts.end();++smaIt){
          const ROMol *matcher=pattern_flyweight(*smaIt).get().getMatcher();
          CHECK_INVARIANT(matcher,"bad smarts");
          featureMatchers.push_back(matcher);
        }
        patterns=&featureMatchers;
      }

      types.resize(patterns->size());
    
      for(unsigned int i=0;i<patterns->size();++i){
        types[i].resize(nAtoms);
        types[i].reset();
        unsigned int mask=1<<i;
        std::vector<MatchVectType> matchVect;
        // to maintain thread safety, we have to copy the pattern
        // molecules:
        SubstructMatch(mol,ROMol(*(*patterns)[i],true),matchVect);
        for(std::vector<MatchVectType>::const_iterator mvIt=matchVect.begin();
            mvIt!=matchVect.end();++mvIt){
          types[i].set((*mvIt)[0].second);
        }
      }
    } // end of getAtomTypes;
  } // end of namespace ReducedGraphs

  ROMol *createMolExtendedReducedGraph(const ROMol &mol,
                                       std::vector<boost::dynamic_bitset<> > *atomTypes
                                       ){
    std::vector<boost::dynamic_bitset<> > *latomTypes=0;    
    if(!atomTypes){
      latomTypes = new std::vector<boost::dynamic_bitset<> >();
      atomTypes = latomTypes;
      ReducedGraphs::getErGAtomTypes(mol,*atomTypes);
    }
    RWMol *res = new RWMol(mol);

    const int aromaticFlag = atomTypes->size();
    const int aliphaticFlag = atomTypes->size()+1;

    for(ROMol::AtomIterator atIt=res->beginAtoms();atIt!=res->endAtoms();++atIt){
      std::list<int> tv;
      for(unsigned int i=0;i<atomTypes->size();++i){
        if((*atomTypes)[i][(*atIt)->getIdx()]) tv.push_back(i);
      }
      (*atIt)->setProp("_ErGAtomTypes",tv,true);
    }

    // start by adding dummies at the ring centroids
    BOOST_FOREACH(const INT_VECT &ring,mol.getRingInfo()->atomRings()){
      if(ring.size()<8){
        int nIdx=res->addAtom(new Atom(0),false,true)-1;
        int nAromatic=0,nSP2=0;
        BOOST_FOREACH(int idx,ring){
          res->addBond(idx,nIdx,Bond::SINGLE);
          if(mol.getAtomWithIdx(idx)->getIsAromatic()){
            ++nAromatic;
          } else if(mol.getAtomWithIdx(idx)->getHybridization()==Atom::SP2) {
            ++nSP2;
          }
        }
        std::list<int> tv;
        if(nAromatic>=2 || nSP2 >= ring.size()/2) tv.push_back(aromaticFlag);
        else tv.push_back(aliphaticFlag);
        res->getAtomWithIdx(nIdx)->setProp("_ErGAtomTypes",tv,true);
      }
    }

    // now remove any degree-two ring atoms that have no features:
    for(unsigned int i=mol.getNumAtoms()-1;i>=0;++i){
      if(mol.getRingInfo()->numAtomRings(i) && mol.getAtomWithIdx(i)->getDegree()==2){
        std::list<int> li;
        mol.getAtomWithIdx(i)->getProp("_ErGAtomTypes",li);
        if(li.empty()){
          res->removeAtom(i);
        }
      }
    }

    // FIX: still need to do the "highly fused rings" simplification for things like adamantane

    
    
    
    if(latomTypes) delete latomTypes;
    return res;
  }
  
} // end of namespace RDKit
