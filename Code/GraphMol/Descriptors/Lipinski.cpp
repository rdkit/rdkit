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
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <vector>
#include <string>

#define SMARTSCOUNTFUNC(nm,pattern,vers)                             \
const std::string nm ## Version  =vers; \
unsigned int calc##nm(const RDKit::ROMol &mol){        \
  static ROMOL_SPTR matcher_##nm;\
  if(!matcher_##nm.get()){\
    matcher_##nm.reset(SmartsToMol(pattern));\
  }\
  TEST_ASSERT(matcher_##nm.get());        \
  std::vector< MatchVectType > matches; \
  int res=SubstructMatch(mol,*(matcher_##nm.get()),matches);    \
  return static_cast<unsigned int>(res);\
}\
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
        if( ((*iter)->getAtomicNum()==7 || (*iter)->getAtomicNum()==8) &&
            ((*iter)->getNumImplicitHs() || (*iter)->getNumExplicitHs() || (*iter)->getTotalNumHs(true)) ) ++res;
      }
      return res;
    }

    SMARTSCOUNTFUNC(NumRotatableBonds, "[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]", "1.0.1" ) ;
    SMARTSCOUNTFUNC(NumHBD, "[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]","2.0.1" ) ;
    SMARTSCOUNTFUNC(NumHBA, "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0,o,s;+0])]","2.0.1") ;
    SMARTSCOUNTFUNC(NumHeteroatoms,"[!#6;!#1]","1.0.1") ;
    
    const std::string NumRingsVersion="1.0.1";
    unsigned int calcNumRings(const ROMol &mol){
      return mol.getRingInfo()->numRings();
    }

  } // end of namespace Descriptors
} // end of namespace RDKit
