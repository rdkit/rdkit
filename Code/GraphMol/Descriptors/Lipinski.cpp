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


typedef boost::flyweight<boost::flyweights::key_value<std::string,ss_matcher>,boost::flyweights::no_tracking > pattern_flyweight;
#define SMARTSCOUNTFUNC(nm,pattern,vers)                             \
const std::string nm ## Version  =vers; \
unsigned int calc##nm(const RDKit::ROMol &mol){        \
  pattern_flyweight m(pattern);                        \
  const ROMol *matcher=m.get().getMatcher();           \
  TEST_ASSERT(matcher);                                \
  std::vector< MatchVectType > matches;                \
  int res=0;                                           \
  if(std::string(pattern).find_first_of("$")!=std::string::npos){       \
    const ROMol nm(*(matcher),true);                   \
    res=SubstructMatch(mol,nm,matches);                \
  } else {                                             \
    const ROMol &nm=*matcher;                          \
    res=SubstructMatch(mol,nm,matches);                \
  }                                                    \
  return static_cast<unsigned int>(res);               \
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
