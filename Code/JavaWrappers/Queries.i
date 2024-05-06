//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

%{
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/types.h>
%}

%include <GraphMol/RDKitBase.h>;
%include <GraphMol/RDKitQueries.h>;
%include <RDGeneral/types.h>;

/*
  Adapted from the Python Wrappers

  NOTE: it looks like there is a typo in the below code
  ATOM_GREATER_QUERY is intentionally being used for the LessQueryAtom
  and ATOM_LESS_QUERY for GreaterQueryAtom in the python API.
  The C++ API is internally consistent and logical, but having
  AtomNumLessQueryAtom(6) return atoms where 6 is less than their atomic
  number feels backwards in Python.

 */
 %define QAFUNC1(funcname, func, type)
 %inline %{

   RDKit::QueryAtom* funcname ## EqualsQueryAtom(type val, bool negate=false) {         
    auto *res = new RDKit::QueryAtom();                                      
    res->setQuery(RDKit:: ## func(val));                                            
    if (negate) res->getQuery()->setNegation(true);                        
    return res;                                                            
  }
  
  RDKit::QueryAtom* funcname ## LessQueryAtom(type val, bool negate=false) {           
    auto *res = new RDKit::QueryAtom();                                      
    res->setQuery(                                                         
        RDKit:: ## func <RDKit::ATOM_GREATER_QUERY>(val, std::string( # funcname "Less"))); 
    if (negate) res->getQuery()->setNegation(true);                        
    return res;                                                            
  }                                                                        

  RDKit::QueryAtom* funcname ## GreaterQueryAtom(type val, bool negate=false) {        
    auto *res = new RDKit::QueryAtom();                                      
    res->setQuery(                                                         
        RDKit:: ## func <RDKit::ATOM_LESS_QUERY>(val, std::string(# funcname "Greater"))); 
    if (negate) res->getQuery()->setNegation(true);                        
    return res;                                                            
  }
%}
%enddef

%define QAFUNC2(funcname, func, type)          
%inline %{
  RDKit::QueryAtom* funcname(bool negate) {              
    auto *res = new RDKit::QueryAtom();               
    res->setQuery(RDKit:: ## func());                        
    if (negate) res->getQuery()->setNegation(true); 
    return res;                                     
  }
%}
%enddef

%newobject AtomNumEqualsQueryAtom;
%newobject AtomNumLessQueryAtom;
%newobject AtomNumGreaterQueryAtom;
QAFUNC1(AtomNum, makeAtomNumQuery, int);

QAFUNC1(ExplicitValence, makeAtomExplicitValenceQuery, int);
QAFUNC1(TotalValence, makeAtomTotalValenceQuery, int);
QAFUNC1(ExplicitDegree, makeAtomExplicitDegreeQuery, int);
QAFUNC1(TotalDegree, makeAtomTotalDegreeQuery, int);
QAFUNC1(NonHydrogenDegree, makeAtomNonHydrogenDegreeQuery, int);

QAFUNC1(HCount, makeAtomHCountQuery, int);
QAFUNC1(Mass, makeAtomMassQuery, int);
QAFUNC1(Isotope, makeAtomIsotopeQuery, int);
QAFUNC1(FormalCharge, makeAtomFormalChargeQuery, int);

QAFUNC1(Hybridization, makeAtomHybridizationQuery, int);
QAFUNC1(InNRings, makeAtomInNRingsQuery, int);
QAFUNC1(MinRingSize, makeAtomMinRingSizeQuery, int);
QAFUNC1(RingBondCount, makeAtomRingBondCountQuery, int);
QAFUNC1(NumRadicalElectrons, makeAtomNumRadicalElectronsQuery, int);
QAFUNC1(NumHeteroatomNeighbors, makeAtomNumHeteroatomNbrsQuery, int);
QAFUNC1(NumAliphaticHeteroatomNeighbors,
        makeAtomNumAliphaticHeteroatomNbrsQuery, int);

%newobject IsUnsaturatedQueryAtom;
QAFUNC2(IsUnsaturatedQueryAtom, makeAtomUnsaturatedQuery, int);

QAFUNC2(IsAromaticQueryAtom, makeAtomAromaticQuery, int);
QAFUNC2(IsAliphaticQueryAtom, makeAtomAliphaticQuery, int);
QAFUNC2(IsInRingQueryAtom, makeAtomInRingQuery, int);
QAFUNC2(HasChiralTagQueryAtom, makeAtomHasChiralTagQuery, int);
QAFUNC2(MissingChiralTagQueryAtom, makeAtomMissingChiralTagQuery, int);
QAFUNC2(IsBridgeheadQueryAtom, makeAtomIsBridgeheadQuery, int);

QAFUNC2(AAtomQueryAtom, makeAAtomQuery, int);
QAFUNC2(AHAtomQueryAtom, makeAHAtomQuery, int);
QAFUNC2(QAtomQueryAtom, makeQAtomQuery, int);
QAFUNC2(QHAtomQueryAtom, makeQHAtomQuery, int);
QAFUNC2(XAtomQueryAtom, makeXAtomQuery, int);
QAFUNC2(XHAtomQueryAtom, makeXHAtomQuery, int);
QAFUNC2(MAtomQueryAtom, makeMAtomQuery, int);
QAFUNC2(MHAtomQueryAtom, makeMHAtomQuery, int);

%extend RDKit::QueryAtom {
  void ExpandQuery(const RDKit::QueryAtom *other, Queries::CompositeQueryType how=Queries::COMPOSITE_AND, bool maintainOrder=true) {
  if (other->hasQuery()) {
	  const RDKit::QueryAtom::QUERYATOM_QUERY *qry = other->getQuery();
	  ($self)->expandQuery(qry->copy(), how, maintainOrder);
	}
  }
}

%extend RDKit::Atom {
  void ExpandQuery(const RDKit::QueryAtom *other, Queries::CompositeQueryType how=Queries::COMPOSITE_AND, bool maintainOrder=true) {
  if (other->hasQuery()) {
	  const RDKit::QueryAtom::QUERYATOM_QUERY *qry = other->getQuery();
	  ($self)->expandQuery(qry->copy(), how, maintainOrder);
	}
  }
}


