// $Id$
//
//  Copyright (C) 2011-2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <vector>
#include <string>

  namespace {
    class ss_matcher {
    public:
      ss_matcher(const std::string &pattern) : m_pattern(pattern){
        m_needCopies=(pattern.find_first_of("$")!=std::string::npos);
        RDKit::RWMol *p=RDKit::SmartsToMol(pattern);
        m_matcher=p;
        POSTCONDITION(m_matcher,"no matcher");
      };
      const RDKit::ROMol *getMatcher() const { return m_matcher; };
      unsigned int countMatches(const RDKit::ROMol &mol) const {
        PRECONDITION(m_matcher,"no matcher");
        std::vector<RDKit::MatchVectType> matches;
        // This is an ugly one. Recursive queries aren't thread safe.
        // Unfortunately we have to take a performance hit here in order
        // to guarantee thread safety
        if(m_needCopies){
	  const RDKit::ROMol nm(*(m_matcher),true);
          RDKit::SubstructMatch(mol,nm,matches);
        } else {
          const RDKit::ROMol &nm=*m_matcher;
          RDKit::SubstructMatch(mol,nm,matches);
        }
        return matches.size();
      }
      ~ss_matcher() { delete m_matcher; };
    private:
      ss_matcher() : m_pattern(""), m_needCopies(false), m_matcher(0) {};
      std::string m_pattern;
      bool m_needCopies;
      const RDKit::ROMol *m_matcher;
    };
  }


typedef boost::flyweight<boost::flyweights::key_value<std::string,ss_matcher>,boost::flyweights::no_tracking > pattern_flyweight;
#define SMARTSCOUNTFUNC(nm,pattern,vers)                             \
const std::string nm ## Version  =vers; \
unsigned int calc##nm(const RDKit::ROMol &mol){        \
  pattern_flyweight m(pattern);                        \
  return m.get().countMatches(mol);                    \
}                                                      \
extern int no_such_variable

namespace RDKit{
  namespace Descriptors {
    unsigned int calcLipinskiHBA(const ROMol &mol){
      unsigned int res=0;
      for(ROMol::ConstAtomIterator iter=mol.beginAtoms();
          iter!=mol.endAtoms();++iter){
        if((*iter)->getAtomicNum()==7 || (*iter)->getAtomicNum()==8) ++res;
      }
      return res;
    }
    unsigned int calcLipinskiHBD(const ROMol &mol){
      unsigned int res=0;
      for(ROMol::ConstAtomIterator iter=mol.beginAtoms();
          iter!=mol.endAtoms();++iter){
        if( ((*iter)->getAtomicNum()==7 || (*iter)->getAtomicNum()==8) ) {
          res += (*iter)->getTotalNumHs(true);
        }
      }
      return res;
    }

    const std::string NumRotatableBondsVersion="2.0.0";
    unsigned int calcNumRotatableBonds(const ROMol &mol,bool strict){
      if(strict){
        std::string strict_pattern="[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-!@[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])]";
        pattern_flyweight m(strict_pattern);
        return m.get().countMatches(mol);
      } else {
        std::string pattern="[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]";
        pattern_flyweight m(pattern);
        return m.get().countMatches(mol);
      }
    }

    //SMARTSCOUNTFUNC(NumHBD, "[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]","2.0.1" ) ;
    SMARTSCOUNTFUNC(NumHBD, "[N&!H0&v3,N&!H0&+1&v4,O&H1&+0,S&H1&+0,n&H1&+0]","2.0.1" ) ;
    SMARTSCOUNTFUNC(NumHBA, "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0,o,s;+0])]","2.0.1") ;
    SMARTSCOUNTFUNC(NumHeteroatoms,"[!#6;!#1]","1.0.1") ;
    SMARTSCOUNTFUNC(NumAmideBonds,"C(=[O;!R])N","1.0.0") ;
    
    const std::string NumRingsVersion="1.0.1";
    unsigned int calcNumRings(const ROMol &mol){
      return mol.getRingInfo()->numRings();
    }

    const std::string FractionCSP3Version="1.0.0";
    double calcFractionCSP3(const ROMol &mol){
      unsigned int nCSP3=0;
      unsigned int nC=0;
      ROMol::VERTEX_ITER atBegin,atEnd;
      boost::tie(atBegin,atEnd) = mol.getVertices();  
      while(atBegin!=atEnd){
        ATOM_SPTR at=mol[*atBegin];
        if(at->getAtomicNum()==6){
          ++nC;
          if(at->getTotalDegree()==4){
            ++nCSP3;
          }
        }
        ++atBegin;
      }
      if(!nC) return 0;
      return static_cast<double>(nCSP3)/nC;
    }

    const std::string NumHeterocyclesVersion="1.0.0";
    unsigned int calcNumHeterocycles(const ROMol &mol){
      unsigned int res=0;
      BOOST_FOREACH(const INT_VECT &iv,mol.getRingInfo()->atomRings()){
        BOOST_FOREACH(int i,iv){
          if(mol.getAtomWithIdx(i)->getAtomicNum()!=6){
            ++res;
            break;
          }
        }
      }
      return res;
    }
    const std::string NumAromaticRingsVersion="1.0.0";
    unsigned int calcNumAromaticRings(const ROMol &mol){
      unsigned int res=0;
      BOOST_FOREACH(const INT_VECT &iv,mol.getRingInfo()->bondRings()){
        ++res;
        BOOST_FOREACH(int i,iv){
          if(!mol.getBondWithIdx(i)->getIsAromatic()){
            --res;
            break;
          }
        }
      }
      return res;
    }
    const std::string NumSaturatedRingsVersion="1.0.0";
    unsigned int calcNumSaturatedRings(const ROMol &mol){
      unsigned int res=0;
      BOOST_FOREACH(const INT_VECT &iv,mol.getRingInfo()->bondRings()){
        ++res;
        BOOST_FOREACH(int i,iv){
          if(mol.getBondWithIdx(i)->getBondType()!=Bond::SINGLE || mol.getBondWithIdx(i)->getIsAromatic()){
            --res;
            break;
          }
        }
      }
      return res;
    }
    const std::string NumAliphaticRingsVersion="1.0.0";
    unsigned int calcNumAliphaticRings(const ROMol &mol){
      unsigned int res=0;
      BOOST_FOREACH(const INT_VECT &iv,mol.getRingInfo()->bondRings()){
        BOOST_FOREACH(int i,iv){
          if(!mol.getBondWithIdx(i)->getIsAromatic()){
            ++res;
            break;
          }
        }
      }
      return res;
    }
    const std::string NumAromaticHeterocyclesVersion="1.0.0";
    unsigned int calcNumAromaticHeterocycles(const ROMol &mol){
      unsigned int res=0;
      BOOST_FOREACH(const INT_VECT &iv,mol.getRingInfo()->bondRings()){
        bool countIt=false;
        BOOST_FOREACH(int i,iv){
          if(!mol.getBondWithIdx(i)->getIsAromatic()){
            countIt=false;
            break;
          }
          // we're checking each atom twice, which is kind of doofy, but this
          // function is hopefully not going to be a big time sink.
          if(!countIt &&
             (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum()!=6 ||
              mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum()!=6) ){
            countIt=true;
          }
        }
        if(countIt) ++res;
      }
      return res;
    }
    const std::string NumAromaticCarbocyclesVersion="1.0.0";
    unsigned int calcNumAromaticCarbocycles(const ROMol &mol){
      unsigned int res=0;
      BOOST_FOREACH(const INT_VECT &iv,mol.getRingInfo()->bondRings()){
        bool countIt=true;
        BOOST_FOREACH(int i,iv){
          if(!mol.getBondWithIdx(i)->getIsAromatic()){
            countIt=false;
            break;
          }
          // we're checking each atom twice, which is kind of doofy, but this
          // function is hopefully not going to be a big time sync.
          if(mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum()!=6 ||
             mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum()!=6 ){
            countIt=false;
            break;
          }
        }
        if(countIt) ++res;
      }
      return res;
    }
    const std::string NumAliphaticHeterocyclesVersion="1.0.0";
    unsigned int calcNumAliphaticHeterocycles(const ROMol &mol){
      unsigned int res=0;
      BOOST_FOREACH(const INT_VECT &iv,mol.getRingInfo()->bondRings()){
        bool hasAliph=false;
        bool hasHetero=false;
        BOOST_FOREACH(int i,iv){
          if(!mol.getBondWithIdx(i)->getIsAromatic()){
            hasAliph=true;
          }
          // we're checking each atom twice, which is kind of doofy, but this
          // function is hopefully not going to be a big time sink.
          if(!hasHetero &&
             (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum()!=6 ||
              mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum()!=6) ){
            hasHetero=true;
          }
        }
        if(hasHetero&&hasAliph) ++res;
      }
      return res;
    }
    const std::string NumAliphaticCarbocyclesVersion="1.0.0";
    unsigned int calcNumAliphaticCarbocycles(const ROMol &mol){
      unsigned int res=0;
      BOOST_FOREACH(const INT_VECT &iv,mol.getRingInfo()->bondRings()){
        bool hasAliph=false;
        bool hasHetero=false;
        BOOST_FOREACH(int i,iv){
          if(!mol.getBondWithIdx(i)->getIsAromatic()){
            hasAliph=true;
          }
          // we're checking each atom twice, which is kind of doofy, but this
          // function is hopefully not going to be a big time sync.
          if(mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum()!=6 ||
             mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum()!=6 ){
            hasHetero=true;
            break;
          }
        }
        if(hasAliph&&!hasHetero) ++res;
      }
      return res;
    }
    const std::string NumSaturatedHeterocyclesVersion="1.0.0";
    unsigned int calcNumSaturatedHeterocycles(const ROMol &mol){
      unsigned int res=0;
      BOOST_FOREACH(const INT_VECT &iv,mol.getRingInfo()->bondRings()){
        bool countIt=false;
        BOOST_FOREACH(int i,iv){
          if(mol.getBondWithIdx(i)->getBondType()!=Bond::SINGLE || mol.getBondWithIdx(i)->getIsAromatic()){
            countIt=false;
            break;
          }
          // we're checking each atom twice, which is kind of doofy, but this
          // function is hopefully not going to be a big time sync.
          if(!countIt &&
             (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum()!=6 ||
              mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum()!=6) ){
            countIt=true;
          }
        }
        if(countIt) ++res;
      }
      return res;
    }
    const std::string NumSaturatedCarbocyclesVersion="1.0.0";
    unsigned int calcNumSaturatedCarbocycles(const ROMol &mol){
      unsigned int res=0;
      BOOST_FOREACH(const INT_VECT &iv,mol.getRingInfo()->bondRings()){
        bool countIt=true;
        BOOST_FOREACH(int i,iv){
          if(mol.getBondWithIdx(i)->getBondType()!=Bond::SINGLE || mol.getBondWithIdx(i)->getIsAromatic()){
            countIt=false;
            break;
          }
          // we're checking each atom twice, which is kind of doofy, but this
          // function is hopefully not going to be a big time sync.
          if(mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum()!=6 ||
             mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum()!=6 ){
            countIt=false;
            break;
          }
        }
        if(countIt) ++res;
      }
      return res;
    }
  } // end of namespace Descriptors
} // end of namespace RDKit
