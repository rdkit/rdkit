//
//  Copyright (C) 2003-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <string>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/types.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {

/*
  NOTE: it looks like there is a typo in the below code
  ATOM_GREATER_QUERY is intentionally being used for the LessQueryAtom
  and ATOM_LESS_QUERY for GreaterQueryAtom in the python API.
  The C++ API is internally consistent and logical, but having
  AtomNumLessQueryAtom(6) return atoms where 6 is less than their atomic
  number feels backwards in Python.

 */
#define QAFUNC1(_funcname_, _func_, _typ_)                                 \
  QueryAtom *_funcname_##EqualsQueryAtom(_typ_ val, bool negate) {         \
    QueryAtom *res = new QueryAtom();                                      \
    res->setQuery(_func_(val));                                            \
    if (negate) {                                                          \
      res->getQuery()->setNegation(true);                                  \
    }                                                                      \
    return res;                                                            \
  }                                                                        \
  QueryAtom *_funcname_##LessQueryAtom(_typ_ val, bool negate) {           \
    QueryAtom *res = new QueryAtom();                                      \
    res->setQuery(                                                         \
        _func_<ATOM_GREATER_QUERY>(val, std::string(#_funcname_ "Less"))); \
    if (negate) {                                                          \
      res->getQuery()->setNegation(true);                                  \
    }                                                                      \
    return res;                                                            \
  }                                                                        \
  QueryAtom *_funcname_##GreaterQueryAtom(_typ_ val, bool negate) {        \
    QueryAtom *res = new QueryAtom();                                      \
    res->setQuery(                                                         \
        _func_<ATOM_LESS_QUERY>(val, std::string(#_funcname_ "Greater"))); \
    if (negate) {                                                          \
      res->getQuery()->setNegation(true);                                  \
    }                                                                      \
    return res;                                                            \
  }

#define QAFUNC2(_funcname_, _func_, _typ_) \
  QueryAtom *_funcname_(bool negate) {     \
    QueryAtom *res = new QueryAtom();      \
    res->setQuery(_func_());               \
    if (negate) {                          \
      res->getQuery()->setNegation(true);  \
    }                                      \
    return res;                            \
  }
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

QueryAtom *HasPropQueryAtom(const std::string &propname, bool negate) {
  auto *res = new QueryAtom();
  res->setQuery(makeHasPropQuery<Atom>(propname));
  if (negate) {
    res->getQuery()->setNegation(true);
  }
  return res;
}

QueryBond *HasPropQueryBond(const std::string &propname, bool negate) {
  auto *res = new QueryBond();
  res->setQuery(makeHasPropQuery<Bond>(propname));
  if (negate) {
    res->getQuery()->setNegation(true);
  }
  return res;
}

template <class Ob, class Ret, class T>
Ret *PropQuery(const std::string &propname, const T &v, bool negate) {
  auto *res = new Ret();
  res->setQuery(makePropQuery<Ob, T>(propname, v));
  if (negate) {
    res->getQuery()->setNegation(true);
  }
  return res;
}

template <class Ob, class Ret, class T>
Ret *PropQueryWithTol(const std::string &propname, const T &v, bool negate,
                      const T &tol = T()) {
  auto *res = new Ret();
  res->setQuery(makePropQuery<Ob, T>(propname, v, tol));
  if (negate) {
    res->getQuery()->setNegation(true);
  }
  return res;
}

template <class Ob, class Ret>
Ret *PropQueryWithTol(const std::string &propname, const ExplicitBitVect &v,
                      bool negate, float tol = 0.0) {
  auto *res = new Ret();
  res->setQuery(makePropQuery<Ob>(propname, v, tol));
  if (negate) {
    res->getQuery()->setNegation(true);
  }
  return res;
}

namespace {
Atom *replaceAtomWithQueryAtomHelper(ROMol &mol, Atom &atom) {
  return QueryOps::replaceAtomWithQueryAtom(static_cast<RWMol *>(&mol), &atom);
}
}  // namespace

struct queries_wrapper {
  static void wrap(nb::module_ &m) {
#define QADEF1(_funcname_)                                                     \
  m.def(#_funcname_ "EqualsQueryAtom", _funcname_##EqualsQueryAtom, "val"_a,   \
        "negate"_a = false,                                                    \
        R"DOC(Returns a QueryAtom that matches atoms where )DOC" #_funcname_   \
        R"DOC( is equal to the target value.)DOC",                             \
        nb::rv_policy::take_ownership);                                        \
  m.def(#_funcname_ "LessQueryAtom", _funcname_##LessQueryAtom, "val"_a,       \
        "negate"_a = false,                                                    \
        R"DOC(Returns a QueryAtom that matches atoms where )DOC" #_funcname_   \
        R"DOC( is less than the target value.                      \
NOTE: the direction of comparison is reversed relative to the C++ API)DOC",    \
        nb::rv_policy::take_ownership);                                        \
  m.def(#_funcname_ "GreaterQueryAtom", _funcname_##GreaterQueryAtom, "val"_a, \
        "negate"_a = false,                                                    \
        R"DOC(Returns a QueryAtom that matches atoms where )DOC" #_funcname_   \
        R"DOC( is equal to the target value.                       \
NOTE: the direction of comparison is reversed relative to the C++ API)DOC",    \
        nb::rv_policy::take_ownership);
#define QADEF2(_funcname_)                                                  \
  m.def(#_funcname_ "QueryAtom", _funcname_##QueryAtom, "negate"_a = false, \
        R"DOC(Returns a QueryAtom that matches atoms when )DOC" #_funcname_ \
        R"DOC( is True.)DOC",                                               \
        nb::rv_policy::take_ownership);

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
    QADEF1(NumRadicalElectrons)
    QADEF1(NumHeteroatomNeighbors)
    QADEF1(NumAliphaticHeteroatomNeighbors)
    QADEF1(NonHydrogenDegree)

    QADEF2(IsUnsaturated);
    QADEF2(IsAromatic);
    QADEF2(IsAliphatic);
    QADEF2(IsInRing);
    QADEF2(HasChiralTag);
    QADEF2(MissingChiralTag);
    QADEF2(IsBridgehead);

    QADEF2(AAtom);
    QADEF2(AHAtom);
    QADEF2(XAtom);
    QADEF2(XHAtom);
    QADEF2(QAtom);
    QADEF2(QHAtom);
    QADEF2(MAtom);
    QADEF2(MHAtom);

    m.def("HasPropQueryAtom", HasPropQueryAtom, "propname"_a,
          "negate"_a = false,
          R"DOC(Returns a QueryAtom that matches when the property 'propname'
exists in the atom.)DOC",
          nb::rv_policy::take_ownership);

    m.def("HasPropQueryBond", HasPropQueryBond, "propname"_a,
          "negate"_a = false,
          R"DOC(Returns a QueryBond that matches when the property 'propname'
exists in the bond.)DOC",
          nb::rv_policy::take_ownership);

    m.def("HasIntPropWithValueQueryAtom",
          PropQueryWithTol<Atom, QueryAtom, int>, "propname"_a, "val"_a,
          "negate"_a = false, "tolerance"_a = 0,
          R"DOC(Returns a QueryAtom that matches when the property 'propname'
has the specified int value.)DOC",
          nb::rv_policy::take_ownership);

    m.def("HasBoolPropWithValueQueryAtom", PropQuery<Atom, QueryAtom, bool>,
          "propname"_a, "val"_a, "negate"_a = false,
          R"DOC(Returns a QueryAtom that matches when the property 'propname'
has the specified boolean value.)DOC",
          nb::rv_policy::take_ownership);

    m.def("HasStringPropWithValueQueryAtom",
          PropQuery<Atom, QueryAtom, std::string>, "propname"_a, "val"_a,
          "negate"_a = false,
          R"DOC(Returns a QueryAtom that matches when the property 'propname'
has the specified string value.)DOC",
          nb::rv_policy::take_ownership);

    m.def("HasDoublePropWithValueQueryAtom",
          PropQueryWithTol<Atom, QueryAtom, double>, "propname"_a, "val"_a,
          "negate"_a = false, "tolerance"_a = 0.0,
          R"DOC(Returns a QueryAtom that matches when the property 'propname'
has the specified value +- tolerance)DOC",
          nb::rv_policy::take_ownership);

    m.def("HasBitVectPropWithValueQueryAtom", PropQueryWithTol<Atom, QueryAtom>,
          "propname"_a, "val"_a, "negate"_a = false, "tolerance"_a = 0,
          R"DOC(Returns a QueryAtom that matches when the property 'propname'
has the specified explicit bit vector value.  The Tolerance is the allowed
Tanimoto difference)DOC",
          nb::rv_policy::take_ownership);

    // Bond Queries
    m.def("HasPropQueryBond", HasPropQueryBond, "propname"_a,
          "negate"_a = false,
          R"DOC(Returns a QueryBond that matches when the property 'propname'
exists in the bond.)DOC",
          nb::rv_policy::take_ownership);

    m.def("HasIntPropWithValueQueryBond",
          PropQueryWithTol<Bond, QueryBond, int>, "propname"_a, "val"_a,
          "negate"_a = false, "tolerance"_a = 0,
          R"DOC(Returns a QueryBond that matches when the property 'propname'
has the specified int value.)DOC",
          nb::rv_policy::take_ownership);

    m.def("HasBoolPropWithValueQueryBond", PropQuery<Bond, QueryBond, bool>,
          "propname"_a, "val"_a, "negate"_a = false,
          R"DOC(Returns a QueryBond that matches when the property 'propname'
has the specified boolean value.)DOC",
          nb::rv_policy::take_ownership);

    m.def("HasStringPropWithValueQueryBond",
          PropQuery<Bond, QueryBond, std::string>, "propname"_a, "val"_a,
          "negate"_a = false,
          R"DOC(Returns a QueryBond that matches when the property 'propname'
has the specified string value.)DOC",
          nb::rv_policy::take_ownership);

    m.def("HasDoublePropWithValueQueryBond",
          PropQueryWithTol<Bond, QueryBond, double>, "propname"_a, "val"_a,
          "negate"_a = false, "tolerance"_a = 0.0,
          R"DOC(Returns a QueryBond that matches when the property 'propname'
has the specified value +- tolerance)DOC",
          nb::rv_policy::take_ownership);

    m.def(
        "ReplaceAtomWithQueryAtom", replaceAtomWithQueryAtomHelper, "mol"_a,
        "atom"_a,
        R"DOC(Changes the given atom in the molecule to a query atom and returns
the atom which can then be modified, for example with additional query
constraints added.  The new atom is otherwise a copy of the old.
If the atom already has a query, nothing will be changed.)DOC",
        nb::rv_policy::reference_internal);
  };
};
void wrap_queries(nb::module_ &m) { RDKit::queries_wrapper::wrap(m); }
}  // namespace RDKit
