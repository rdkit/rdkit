//
//  Copyright (C) 2024 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// Adapted from the Python Wrappers

%{
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/types.h>
%}

%include <GraphMol/RDKitBase.h>;
%include <GraphMol/RDKitQueries.h>;
%include <RDGeneral/types.h>;

/*
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
    std::unique_ptr<RDKit::QueryAtom> res(new RDKit::QueryAtom());
    res->setQuery(RDKit:: ## func(val));                                            
    if (negate) {
      res->getQuery()->setNegation(true);
    }
    return res.release();
  }
  
  RDKit::QueryAtom* funcname ## LessQueryAtom(type val, bool negate=false) {           
    std::unique_ptr<RDKit::QueryAtom> res(new RDKit::QueryAtom());
    res->setQuery(                                                         
        RDKit:: ## func <RDKit::ATOM_GREATER_QUERY>(val, std::string( # funcname "Less"))); 
    if (negate) {
      res->getQuery()->setNegation(true);                        
    }
    return res.release();
  }                                                                        

  RDKit::QueryAtom* funcname ## GreaterQueryAtom(type val, bool negate=false) {        
    std::unique_ptr<RDKit::QueryAtom> res(new RDKit::QueryAtom());
    res->setQuery(                                                         
        RDKit:: ## func <RDKit::ATOM_LESS_QUERY>(val, std::string(# funcname "Greater"))); 
    if (negate) {
      res->getQuery()->setNegation(true);
    }
    return res.release();                      
  }
%}
%enddef

%define QAFUNC2(funcname, func, type)          
%inline %{
  RDKit::QueryAtom* funcname(bool negate=false) {              
    std::unique_ptr<RDKit::QueryAtom> res(new RDKit::QueryAtom());
    res->setQuery(RDKit:: ## func());                        
    if (negate) {
      res->getQuery()->setNegation(true); 
    }
    return res.release();
  }
%}
%enddef

%newobject AtomNumEqualsQueryAtom;
%newobject AtomNumLessQueryAtom;
%newobject AtomNumGreaterQueryAtom;
QAFUNC1(AtomNum, makeAtomNumQuery, int);

%newobject ExplicitValenceEqualsQueryAtom;
%newobject ExplicitValenceLessQueryAtom;
%newobject AtomExplicitValenceGreaterQueryAtom;
QAFUNC1(ExplicitValence, makeAtomExplicitValenceQuery, int);

%newobject TotalValenceEqualsQueryAtom;
%newobject TotalValenceLessQueryAtom;
%newobject TotalValanceGreaterQueryAtom;
QAFUNC1(TotalValence, makeAtomTotalValenceQuery, int);

%newobject ExplicitDegreeEqualsQueryAtom;
%newobject ExplicitDegreeLessQueryAtom;
%newobject ExplicitDegreeGreaterQueryAtom;
QAFUNC1(ExplicitDegree, makeAtomExplicitDegreeQuery, int);

%newobject TotalDegreeEqualsQueryAtom;
%newobject TotalDegreeLessQueryAtom;
%newobject TotalDegreeGreaterQueryAtom;
QAFUNC1(TotalDegree, makeAtomTotalDegreeQuery, int);

%newobject NonHydrogenDegreeEqualsQueryAtom;
%newobject NonHydrogenDegreeLessQueryAtom;
%newobject NonHydrogenDegreeGreaterQueryAtom;
QAFUNC1(NonHydrogenDegree, makeAtomNonHydrogenDegreeQuery, int);

%newobject HCountEqualsQueryAtom;
%newobject HCountLessQueryAtom;
%newobject HCountGreaterQueryAtom;
QAFUNC1(HCount, makeAtomHCountQuery, int);

%newobject MassEqualsQueryAtom;
%newobject MassLessQueryAtom;
%newobject MassGreaterQueryAtom;
QAFUNC1(Mass, makeAtomMassQuery, int);

%newobject IsotopeEqualsQueryAtom;
%newobject IsotopeLessQueryAtom;
%newobject IsotopeGreaterQueryAtom;
QAFUNC1(Isotope, makeAtomIsotopeQuery, int);

%newobject FormalChargeEqualsQueryAtom;
%newobject FormalChargeLessQueryAtom;
%newobject FormalChargeGreaterQueryAtom;
QAFUNC1(FormalCharge, makeAtomFormalChargeQuery, int);

%newobject HybridizationEqualsQueryAtom;
%newobject HybridizationLessQueryAtom;
%newobject HybridizationGreaterQueryAtom;
QAFUNC1(Hybridization, makeAtomHybridizationQuery, int);

%newobject InNRingsEqualsQueryAtom;
%newobject InNRingsLessQueryAtom;
%newobject InNRingsGreaterQueryAtom;
QAFUNC1(InNRings, makeAtomInNRingsQuery, int);

%newobject MinRingSizeEqualsQueryAtom;
%newobject MinRingSizeLessQueryAtom;
%newobject MinRingSizeGreaterQueryAtom;
QAFUNC1(MinRingSize, makeAtomMinRingSizeQuery, int);

%newobject RingBondCountEqualsQueryAtom;
%newobject RingBondCountLessQueryAtom;
%newobject RingBondCountGreaterQueryAtom;
QAFUNC1(RingBondCount, makeAtomRingBondCountQuery, int);

%newobject NumRadicalElectronsEqualsQueryAtom;
%newobject NumRadicalElectronsLessQueryAtom;
%newobject NumRadicalElectronsGreaterQueryAtom;
QAFUNC1(NumRadicalElectrons, makeAtomNumRadicalElectronsQuery, int);

%newobject NumHeteroatomNeighborsEqualsQueryAtom;
%newobject NumHeteroatomNeighborsLessQueryAtom;
%newobject NumHeteroatomNeighborsGreaterQueryAtom;
QAFUNC1(NumHeteroatomNeighbors, makeAtomNumHeteroatomNbrsQuery, int);

%newobject NumAliphaticHeteroatomNeighborsEqualsQueryAtom;
%newobject NumAliphaticHeteroatomNeighborsLessQueryAtom;
%newobject NumAliphaticHeteroatomNeighborsGreaterQueryAtom;
QAFUNC1(NumAliphaticHeteroatomNeighbors,
        makeAtomNumAliphaticHeteroatomNbrsQuery, int);

%newobject IsUnsaturatedQueryAtom;
QAFUNC2(IsUnsaturatedQueryAtom, makeAtomUnsaturatedQuery, int);

%newobject IsAromaticQueryAtom;
QAFUNC2(IsAromaticQueryAtom, makeAtomAromaticQuery, int);

%newobject IsAliphaticQueryAtom;
QAFUNC2(IsAliphaticQueryAtom, makeAtomAliphaticQuery, int);

%newobject IsInRingQueryAtom;
QAFUNC2(IsInRingQueryAtom, makeAtomInRingQuery, int);

%newobject HasChiralTagQueryAtom;
QAFUNC2(HasChiralTagQueryAtom, makeAtomHasChiralTagQuery, int);

%newobject MissingChiralTagQueryAtom;
QAFUNC2(MissingChiralTagQueryAtom, makeAtomMissingChiralTagQuery, int);

%newobject IsBridgeheadQueryAtom;
QAFUNC2(IsBridgeheadQueryAtom, makeAtomIsBridgeheadQuery, int);

%newobject AAtomQueryAtom;
QAFUNC2(AAtomQueryAtom, makeAAtomQuery, int);

%newobject AHAtomQueryAtom;
QAFUNC2(AHAtomQueryAtom, makeAHAtomQuery, int);

%newobject QAtomQueryAtom;
QAFUNC2(QAtomQueryAtom, makeQAtomQuery, int);

%newobject QHAtomQueryAtom;
QAFUNC2(QHAtomQueryAtom, makeQHAtomQuery, int);

%newobject XAtomQueryAtom;
QAFUNC2(XAtomQueryAtom, makeXAtomQuery, int);

%newobject XHAtomQueryAtom;
QAFUNC2(XHAtomQueryAtom, makeXHAtomQuery, int);

%newobject MAtomQueryAtom;
QAFUNC2(MAtomQueryAtom, makeMAtomQuery, int);

%newobject MHAtomQueryAtom;
QAFUNC2(MHAtomQueryAtom, makeMHAtomQuery, int);

%{

  template <class Ob, class Ret, class T>
  Ret *PropQuery(const std::string &propname, const T &v, bool negate) {
	std::unique_ptr<Ret> res(new Ret());
	res->setQuery(RDKit::makePropQuery<Ob, T>(propname, v));
	if (negate) {
	  res->getQuery()->setNegation(true);
	}
	return res.release();
  }

  template <class Ob, class Ret, class T>
  Ret *PropQueryWithTol(const std::string &propname, const T &v, bool negate,
						const T &tol = T()) {
	std::unique_ptr<Ret> res(new Ret());
	res->setQuery(RDKit::makePropQuery<Ob, T>(propname, v, tol));
	if (negate) {
	  res->getQuery()->setNegation(true);
	}
	return res.release();
  }

  template <class Ob, class Ret>
  Ret *PropQueryWithTol(const std::string &propname, const ExplicitBitVect &v,
						bool negate=false, float tol = 0.0) {
	std::unique_ptr<Ret> res(new Ret());
	res->setQuery(RDKit::makePropQuery<Ob>(propname, v, tol));
	if (negate) {
	  res->getQuery()->setNegation(true);
	}
	return res.release();
  }

  RDKit::QueryAtom *HasPropQueryAtom(const std::string &propname, bool negate=false) {
	std::unique_ptr<Ret> res(new Ret());
	res->setQuery(RDKit::makeHasPropQuery<RDKit::Atom>(propname));
	if (negate) {
	  res->getQuery()->setNegation(true);
	}
	return res.release();
  }

  RDKit::QueryAtom *HasIntPropWithValueQueryAtom(const std::string &propname, int val, bool negate=false) {
 	return PropQuery<RDKit::Atom, RDKit::QueryAtom, int>(propname, val, negate); 
  }

  RDKit::QueryAtom *HasBoolPropWithValueQueryAtom(const std::string &propname, bool val, bool negate=false) {
	return PropQuery<RDKit::Atom, RDKit::QueryAtom, bool>(propname, val, negate);
  }

  RDKit::QueryAtom *HasStringPropWithValueQueryAtom(const std::string &propname, const std::string &val, bool negate=false) {
	return PropQuery<RDKit::Atom, RDKit::QueryAtom, std::string>(propname, val, negate);
  }

  RDKit::QueryAtom *HasDoublePropWithValueQueryAtom(const std::string &propname, double val, bool negate=false, double tol=0) {
	return PropQueryWithTol<RDKit::Atom, RDKit::QueryAtom, double>(propname, val, negate, tol);
  }

  RDKit::QueryAtom *HasBitVectPropWithValueQueryAtom(const std::string  &propname, const ExplicitBitVect &val, bool negate=false, float tol=0) {
	return PropQueryWithTol<RDKit::Atom, RDKit::QueryAtom>(propname, val, negate, tol);
  }

  RDKit::QueryBond *HasPropQueryBond(const std::string &propname, bool negate=false) {
	std::unique_ptr<QueryBond> res(new RDKit::QueryBond());
	res->setQuery(RDKit::makeHasPropQuery<RDKit::Bond>(propname));
	if (negate) {
	  res->getQuery()->setNegation(true);
	}
	return res.release();
  }

  RDKit::QueryBond *HasIntPropWithValueQueryBond(const std::string &propname, int val, bool negate=false) {
    return PropQuery<RDKit::Bond, RDKit::QueryBond, int>(propname, val, negate);
  }

  RDKit::QueryBond *HasBoolPropWithValueQueryBond(const std::string &propname, bool val, bool negate=false) {
    return PropQuery<RDKit::Bond, RDKit::QueryBond, bool>(propname, val, negate);
  }

  RDKit::QueryBond *HasStringPropWithValueQueryBond(const std::string &propname, const std::string &val, bool negate=false) {
    return PropQuery<RDKit::Bond, RDKit::QueryBond, std::string>(propname, val, negate);
  }

  RDKit::QueryBond *HasDoublePropWithValueQueryBond(const std::string &propname, double val, bool negate=false, double tol=0) {
    return PropQueryWithTol<RDKit::Bond, RDKit::QueryBond, double>(propname, val, negate, tol);
  }

%}

%newobject HasPropQueryAtom;
RDKit::QueryAtom *HasPropQueryAtom(const std::string &propname, bool negate=false);

%newobject HasPropQueryBond;
RDKit::QueryBond *HasPropQueryBond(const std::string &propname, bool negate=false);

%newobject HasIntPropWithValueQueryAtom;
RDKit::QueryAtom *HasIntPropWithValueQueryAtom(const std::string &propname, int val, bool negate=false);

%newobject HasBoolPropWithValueQueryAtom;
RDKit::QueryBond *HasBoolPropWithValueQueryAtom(const std::string &propname, bool val, bool negate=false);

%newobject HasStringPropWithValueQueryAtom;
RDKit::QueryAtom *HasStringPropWithValueQueryAtom(const std::string &propname, const std::string &val, bool negate=false);

%newobject HasDoublePropWithValueQueryAtom;
RDKit::QueryAtom *HasDoublePropWithValueQueryAtom(const std::string &propname, double val, bool negate=false, double tol=0);

%newobject HasBitVectPropWithValueQueryAtom;
RDKit::QueryAtom *HasBitVectPropWithValueQueryAtom(const std::string  &propname, const ExplicitBitVect &val, bool negate=false, float tol=0);

%newobject HasPropQueryBond;
RDKit::QueryBond *HasPropQueryBond(const std::string &propname, bool negate=false);

%newobject HasIntPropWithValueQueryBond;
RDKit::QueryBond *HasIntPropWithValueQueryBond(const std::string &propname, int val, bool negate=false);

%newobject HasBoolPropWithValueQueryBond;
RDKit::QueryBond *HasBoolPropWithValueQueryBond(const std::string &propname, bool val, bool negate=false);

%newobject HasStringPropWithValueQueryBond;
RDKit::QueryBond *HasStringPropWithValueQueryBond(const std::string &propname, const std::string &val, bool negate=false);

%newobject HasDoublePropWithValueQueryBond;
RDKit::QueryBond *HasDoublePropWithValueQueryBond(const std::string &propname, double val, bool negate=false, double tol=0);

%extend RDKit::QueryAtom {
  void ExpandQuery(const RDKit::QueryAtom *other, Queries::CompositeQueryType how=Queries::COMPOSITE_AND, bool maintainOrder=true) {
	PRECONDITION(other, "bad atoms");
	if (other->hasQuery()) {
	  const RDKit::QueryAtom::QUERYATOM_QUERY *qry = other->getQuery();
	  ($self)->expandQuery(qry->copy(), how, maintainOrder);
	}
  }

  void setQuery(const RDKit::QueryAtom *other) {
	PRECONDITION(other, "bad atoms");
	if (other->hasQuery()) {
	  ($self)->setQuery(other->getQuery()->copy());
	}
  }
}

%extend RDKit::Atom {
  void ExpandQuery(const RDKit::QueryAtom *other, Queries::CompositeQueryType how=Queries::COMPOSITE_AND, bool maintainOrder=true) {
  PRECONDITION(other, "bad atoms");
  if (other->hasQuery()) {
	  const RDKit::QueryAtom::QUERYATOM_QUERY *qry = other->getQuery();
	  ($self)->expandQuery(qry->copy(), how, maintainOrder);
	}
  }

  void setQuery(const RDKit::QueryAtom *other) {
	PRECONDITION(other, "bad atoms");
	if (other->hasQuery()) {
	  ($self)->setQuery(other->getQuery()->copy());
	}
  }
}

%extend RDKit::QueryBond {
  void ExpandQuery(const RDKit::QueryBond *other, Queries::CompositeQueryType how=Queries::COMPOSITE_AND, bool maintainOrder=true) {
	PRECONDITION(other, "bad bonds");
	if (other->hasQuery()) {
	  const RDKit::QueryBond::QUERYBOND_QUERY *qry = other->getQuery();
	  ($self)->expandQuery(qry->copy(), how, maintainOrder);
	}
  }

  void SetQuery(const RDKit::QueryBond *other) {
	PRECONDITION(other, "bad bonds");
	if (other->hasQuery()) {
	  ($self)->setQuery(other->getQuery()->copy());
	}
  }
}

%extend RDKit::Bond {
  void ExpandQuery(const RDKit::QueryBond *other, Queries::CompositeQueryType how=Queries::COMPOSITE_AND, bool maintainOrder=true) {
	PRECONDITION(other, "bad bonds");
	if (other->hasQuery()) {
	  const RDKit::QueryBond::QUERYBOND_QUERY *qry = other->getQuery();
	  ($self)->expandQuery(qry->copy(), how, maintainOrder);
	}
  }

  void SetQuery(const RDKit::QueryBond *other) {
	PRECONDITION(other, "bad bonds");
	if (other->hasQuery()) {
	  ($self)->setQuery(other->getQuery()->copy());
	}
  }
}
