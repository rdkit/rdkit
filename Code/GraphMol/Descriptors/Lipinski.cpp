// $Id$
//
//  Copyright (C) 2011 Greg Landrum
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
        unsigned int res=0;
        // This is an ugly one. Recursive queries aren't thread safe.
        // Unfortunately we have to take a performance hit here in order
        // to guarantee thread safety
        if(m_needCopies){
	  const RDKit::ROMol nm(*(m_matcher),true);
          res=RDKit::SubstructMatch(mol,nm,matches);
        } else {
          const RDKit::ROMol &nm=*m_matcher;
          res=RDKit::SubstructMatch(mol,nm,matches);
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
          res += (*iter)->getTotalNumHs();
        }
      }
      return res;
    }

    SMARTSCOUNTFUNC(NumRotatableBonds, "[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]", "1.0.1" ) ;
    //SMARTSCOUNTFUNC(NumHBD, "[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]","2.0.1" ) ;
    SMARTSCOUNTFUNC(NumHBD, "[N&!H0&v3,N&!H0&+1&v4,O&H1&+0,S&H1&+0,n&H1&+0]","2.0.1" ) ;
    SMARTSCOUNTFUNC(NumHBA, "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0,o,s;+0])]","2.0.1") ;
    SMARTSCOUNTFUNC(NumHeteroatoms,"[!#6;!#1]","1.0.1") ;
    SMARTSCOUNTFUNC(NumAmideBonds,"C(=[O;!R])N","1.0.0") ;
    
    const std::string NumRingsVersion="1.0.1";
    unsigned int calcNumRings(const ROMol &mol){
      return mol.getRingInfo()->numRings();
    }

  } // end of namespace Descriptors
} // end of namespace RDKit
