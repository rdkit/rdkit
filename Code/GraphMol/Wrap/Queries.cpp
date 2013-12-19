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
  
  /*
    NOTE: it looks like there is a typo in the below code
    ATOM_GREATER_QUERY is intentionally being used for the LessQueryAtom
    and ATOM_LESS_QUERY for GreaterQueryAtom in the python API.
    The C++ API is internally consistent and logical, but having
    AtomNumLessQueryAtom(6) return atoms where 6 is less than their atomic
    number feels backwards in Python.
    
   */
#define QAFUNC1(_funcname_,_func_,_typ_) QueryAtom *_funcname_##EqualsQueryAtom(_typ_ val,bool negate) { \
    QueryAtom *res=new QueryAtom();\
    res->setQuery(_func_(val));\
    if(negate) res->getQuery()->setNegation(true);\
    return res;\
  }\
  QueryAtom *_funcname_##LessQueryAtom(_typ_ val,bool negate) { \
    QueryAtom *res=new QueryAtom();\
    res->setQuery(_func_<ATOM_GREATER_QUERY>(val,std::string(#_funcname_ "Less"))); \
    if(negate) res->getQuery()->setNegation(true);\
    return res;\
  }\
  QueryAtom *_funcname_##GreaterQueryAtom(_typ_ val,bool negate) { \
    QueryAtom *res=new QueryAtom();\
    res->setQuery(_func_<ATOM_LESS_QUERY>(val,std::string(#_funcname_ "Greater"))); \
    if(negate) res->getQuery()->setNegation(true);\
    return res;\
  }

#define QAFUNC2(_funcname_,_func_,_typ_) QueryAtom *_funcname_(bool negate) { \
    QueryAtom *res=new QueryAtom();\
    res->setQuery(_func_());\
    if(negate) res->getQuery()->setNegation(true);\
    return res;\
  }
  QAFUNC1(AtomNum,makeAtomNumQuery,int);
  QAFUNC1(ExplicitValence,makeAtomExplicitValenceQuery,int);
  QAFUNC1(TotalValence,makeAtomTotalValenceQuery,int);
  QAFUNC1(ExplicitDegree,makeAtomExplicitDegreeQuery,int);
  QAFUNC1(TotalDegree,makeAtomTotalDegreeQuery,int);

  QAFUNC1(HCount,makeAtomHCountQuery,int);
  QAFUNC1(Mass,makeAtomMassQuery,int);
  QAFUNC1(Isotope,makeAtomIsotopeQuery,int);
  QAFUNC1(FormalCharge,makeAtomFormalChargeQuery,int);

  QAFUNC1(Hybridization,makeAtomHybridizationQuery,int);
  QAFUNC1(InNRings,makeAtomInNRingsQuery,int);
  QAFUNC1(MinRingSize,makeAtomMinRingSizeQuery,int);
  QAFUNC1(RingBondCount,makeAtomRingBondCountQuery,int);

  QAFUNC2(IsUnsaturatedQueryAtom,makeAtomUnsaturatedQuery,int);
  QAFUNC2(IsAromaticQueryAtom,makeAtomAromaticQuery,int);
  QAFUNC2(IsAliphaticQueryAtom,makeAtomAliphaticQuery,int);
  QAFUNC2(IsInRingQueryAtom,makeAtomInRingQuery,int);




struct queries_wrapper {
  static void wrap(){
#define QADEF1(_funcname_)     python::def(# _funcname_ "EqualsQueryAtom",_funcname_##EqualsQueryAtom,\
                                          (python::arg("val"),python::arg("negate")=false), \
                                           "Returns a QueryAtom that matches atoms where " #_funcname_ " is equal to the target value.", \
                                         python::return_value_policy<python::manage_new_object>());\
                               python::def(# _funcname_ "LessQueryAtom",_funcname_##LessQueryAtom,\
                                          (python::arg("val"),python::arg("negate")=false), \
                                           "Returns a QueryAtom that matches atoms where " #_funcname_ " is less than the target value.\n" \
                                           "NOTE: the direction of comparison is reversed relative to the C++ API", \
                                           python::return_value_policy<python::manage_new_object>()); \
                               python::def(# _funcname_ "GreaterQueryAtom",_funcname_##GreaterQueryAtom, \
                                          (python::arg("val"),python::arg("negate")=false), \
                                           "Returns a QueryAtom that matches atoms where " #_funcname_ " is equal to the target value.\n" \
                                           "NOTE: the direction of comparison is reversed relative to the C++ API", \
                                           python::return_value_policy<python::manage_new_object>());
#define QADEF2(_funcname_)     python::def(# _funcname_ "QueryAtom",_funcname_##QueryAtom,\
                                          (python::arg("negate")=false), \
                                           "Returns a QueryAtom that matches atoms when " #_funcname_ " is True.", \
                                           python::return_value_policy<python::manage_new_object>());

    QADEF1(AtomNum);
    QADEF1(ExplicitValence);
    QADEF1(TotalValence);
    QADEF1(ExplicitDegree);
    QADEF1(TotalDegree);
    QADEF1(HCount);
    QADEF1(Mass);
    QADEF1(Isotope);
    QADEF1(FormalCharge);
    QADEF1(Hybridization);
    QADEF1(InNRings);
    QADEF1(MinRingSize);
    QADEF1(RingBondCount);

    QADEF2(IsUnsaturated);
    QADEF2(IsAromatic);
    QADEF2(IsAliphatic);
    QADEF2(IsInRing);

    
  };
};
}// end of namespace


void wrap_queries() {
  RDKit::queries_wrapper::wrap();
}


