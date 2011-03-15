//
//  Copyright (C) 2007-2011 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file Lipinski.h

  \brief Contains Lipinski and Lipinski-like descriptors. Use MolDescriptors.h in client code.

*/
#ifndef __RD_LIPINSKI_H__
#define __RD_LIPINSKI_H__
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <vector>
#include <string>


#define SMARTSCOUNTFUNC(nm,pattern,vers)                             \
static const std::string nm ## Version  =vers; \
static const ROMOL_SPTR matcher_##nm=ROMOL_SPTR(SmartsToMol(pattern)); \
template <class T> \
unsigned int calc##nm(const T &mol){ \
  TEST_ASSERT(matcher_##nm); \
  std::vector< MatchVectType > matches; \
  int res=SubstructMatch(mol,*(matcher_##nm.get()),matches);    \
  return static_cast<unsigned int>(res);\
}\
extern int no_such_variable

namespace RDKit{
  namespace Descriptors {

    static const std::string lipinskiHBAVersion="1.0.0";
    //! calculates the standard Lipinski HBA definition
    /*!  
      \param mol        the molecule of interest

      \return the number of Ns and Os in the molecule

    */
    template <class T>
    unsigned int calcLipinskiHBA(const T &mol){
      unsigned int res=0;
      for(typename T::ConstAtomIterator iter=mol.beginAtoms();
          iter!=mol.endAtoms();++iter){
        if((*iter)->getAtomicNum()==7 || (*iter)->getAtomicNum()==8) ++res;
      }
      return res;
    }


    static const std::string lipinskiHBDVersion="1.0.0";
    //! calculates the standard Lipinski HBA definition
    /*!  
      \param mol        the molecule of interest

      \return the number of Ns and Os with at least one H in the molecule

    */
    template <class T>
    unsigned int calcLipinskiHBD(const T &mol){
      unsigned int res=0;
      for(typename T::ConstAtomIterator iter=mol.beginAtoms();
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

    static const std::string NumRingsVersion="1.0.1";
    template <class T>
    unsigned int calcNumRings(const T &mol){
      return mol.getRingInfo()->numRings();
    }

  } // end of namespace Descriptors
} //end of namespace RDKit

#endif
