// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define NO_IMPORT_ARRAY
#include <boost/python.hpp>
#include <string>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/types.h>

namespace python = boost::python;
namespace RDKit{
  
#define QAFUNC1(_funcname_,_func_,_typ_) QueryAtom *_funcname_(_typ_ val,bool negate) { \
    QueryAtom *res=new QueryAtom();\
    res->setQuery(_func_(val));\
    if(negate) res->getQuery()->setNegation(true);\
    return res;\
  }
#define QAFUNC2(_funcname_,_func_,_typ_) QueryAtom *_funcname_(bool negate) { \
    QueryAtom *res=new QueryAtom();\
    res->setQuery(_func_());\
    if(negate) res->getQuery()->setNegation(true);\
    return res;\
  }
  QAFUNC1(AtomNumEqualsQueryAtom,makeAtomNumQuery,int);
  QAFUNC1(ExplicitValenceEqualsQueryAtom,makeAtomExplicitValenceQuery,int);
  QAFUNC1(TotalValenceEqualsQueryAtom,makeAtomTotalValenceQuery,int);
  QAFUNC1(ExplicitDegreeEqualsQueryAtom,makeAtomExplicitDegreeQuery,int);
  QAFUNC1(TotalDegreeEqualsQueryAtom,makeAtomTotalDegreeQuery,int);

  QAFUNC1(HCountEqualsQueryAtom,makeAtomHCountQuery,int);
  QAFUNC1(MassEqualsQueryAtom,makeAtomMassQuery,int);
  QAFUNC1(IsotopeEqualsQueryAtom,makeAtomIsotopeQuery,int);
  QAFUNC1(FormalChargeEqualsQueryAtom,makeAtomFormalChargeQuery,int);

  QAFUNC1(HybridizationEqualsQueryAtom,makeAtomHybridizationQuery,int);
  QAFUNC1(InNRingsEqualsQueryAtom,makeAtomInNRingsQuery,int);
  QAFUNC1(MinRingSizeEqualsQueryAtom,makeAtomMinRingSizeQuery,int);
  QAFUNC1(RingBondCountEqualsQueryAtom,makeAtomRingBondCountQuery,int);

  QAFUNC2(IsUnsaturatedQueryAtom,makeAtomUnsaturatedQuery,int);
  QAFUNC2(IsAromaticQueryAtom,makeAtomAromaticQuery,int);
  QAFUNC2(IsAliphaticQueryAtom,makeAtomAliphaticQuery,int);
  QAFUNC2(IsInRingQueryAtom,makeAtomInRingQuery,int);




struct queries_wrapper {
  static void wrap(){
#define QADEF1(_funcname_)     python::def(# _funcname_,_funcname_,\
                                          (python::arg("val"),python::arg("negate")=false), \
                                         python::return_value_policy<python::manage_new_object>());
#define QADEF2(_funcname_)     python::def(# _funcname_,_funcname_,\
                                          (python::arg("negate")=false), \
                                         python::return_value_policy<python::manage_new_object>());

    QADEF1(AtomNumEqualsQueryAtom);
    QADEF1(ExplicitValenceEqualsQueryAtom);
    QADEF1(TotalValenceEqualsQueryAtom);
    QADEF1(ExplicitDegreeEqualsQueryAtom);
    QADEF1(TotalDegreeEqualsQueryAtom);
    QADEF1(HCountEqualsQueryAtom);
    QADEF1(MassEqualsQueryAtom);
    QADEF1(IsotopeEqualsQueryAtom);
    QADEF1(FormalChargeEqualsQueryAtom);
    QADEF1(HybridizationEqualsQueryAtom);
    QADEF1(InNRingsEqualsQueryAtom);
    QADEF1(MinRingSizeEqualsQueryAtom);
    QADEF1(RingBondCountEqualsQueryAtom);

    QADEF2(IsUnsaturatedQueryAtom);
    QADEF2(IsAromaticQueryAtom);
    QADEF2(IsAliphaticQueryAtom);
    QADEF2(IsInRingQueryAtom);

    
  };
};
}// end of namespace


void wrap_queries() {
  RDKit::queries_wrapper::wrap();
}


