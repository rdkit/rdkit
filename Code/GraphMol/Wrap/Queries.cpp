//
//  Copyright (C) 2003-2017 Rational Discovery LLC and Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define NO_IMPORT_ARRAY
#include <RDBoost/python.h>
#include <string>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/types.h>

namespace python = boost::python;
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
    if (negate) res->getQuery()->setNegation(true);                        \
    return res;                                                            \
  }                                                                        \
  QueryAtom *_funcname_##LessQueryAtom(_typ_ val, bool negate) {           \
    QueryAtom *res = new QueryAtom();                                      \
    res->setQuery(                                                         \
        _func_<ATOM_GREATER_QUERY>(val, std::string(#_funcname_ "Less"))); \
    if (negate) res->getQuery()->setNegation(true);                        \
    return res;                                                            \
  }                                                                        \
  QueryAtom *_funcname_##GreaterQueryAtom(_typ_ val, bool negate) {        \
    QueryAtom *res = new QueryAtom();                                      \
    res->setQuery(                                                         \
        _func_<ATOM_LESS_QUERY>(val, std::string(#_funcname_ "Greater"))); \
    if (negate) res->getQuery()->setNegation(true);                        \
    return res;                                                            \
  }

#define QAFUNC2(_funcname_, _func_, _typ_)          \
  QueryAtom *_funcname_(bool negate) {              \
    QueryAtom *res = new QueryAtom();               \
    res->setQuery(_func_());                        \
    if (negate) res->getQuery()->setNegation(true); \
    return res;                                     \
  }
QAFUNC1(AtomNum, makeAtomNumQuery, int);
QAFUNC1(ExplicitValence, makeAtomExplicitValenceQuery, int);
QAFUNC1(TotalValence, makeAtomTotalValenceQuery, int);
QAFUNC1(ExplicitDegree, makeAtomExplicitDegreeQuery, int);
QAFUNC1(TotalDegree, makeAtomTotalDegreeQuery, int);

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

QueryAtom *HasPropQueryAtom(const std::string &propname, bool negate) {
  auto *res = new QueryAtom();
  res->setQuery(makeHasPropQuery<Atom>(propname));
  if (negate) res->getQuery()->setNegation(true);
  return res;
}

QueryBond *HasPropQueryBond(const std::string &propname, bool negate) {
  auto *res = new QueryBond();
  res->setQuery(makeHasPropQuery<Bond>(propname));
  if (negate) res->getQuery()->setNegation(true);
  return res;
}

template <class Ob, class Ret, class T>
Ret *PropQuery(const std::string &propname, const T &v, bool negate) {
  auto *res = new Ret();
  res->setQuery(makePropQuery<Ob, T>(propname, v));
  if (negate) res->getQuery()->setNegation(true);
  return res;
}

template <class Ob, class Ret, class T>
Ret *PropQueryWithTol(const std::string &propname, const T &v, bool negate,
                      const T &tol = T()) {
  auto *res = new Ret();
  res->setQuery(makePropQuery<Ob, T>(propname, v, tol));
  if (negate) res->getQuery()->setNegation(true);
  return res;
}

struct queries_wrapper {
  static void wrap() {
#define QADEF1(_funcname_)                                                     \
  python::def(#_funcname_ "EqualsQueryAtom", _funcname_##EqualsQueryAtom,      \
              (python::arg("val"), python::arg("negate") = false),             \
              "Returns a QueryAtom that matches atoms where " #_funcname_      \
              " is equal to the target value.",                                \
              python::return_value_policy<python::manage_new_object>());       \
  python::def(                                                                 \
      #_funcname_ "LessQueryAtom", _funcname_##LessQueryAtom,                  \
      (python::arg("val"), python::arg("negate") = false),                     \
      "Returns a QueryAtom that matches atoms where " #_funcname_              \
      " is less than the target value.\n"                                      \
      "NOTE: the direction of comparison is reversed relative to the C++ API", \
      python::return_value_policy<python::manage_new_object>());               \
  python::def(                                                                 \
      #_funcname_ "GreaterQueryAtom", _funcname_##GreaterQueryAtom,            \
      (python::arg("val"), python::arg("negate") = false),                     \
      "Returns a QueryAtom that matches atoms where " #_funcname_              \
      " is equal to the target value.\n"                                       \
      "NOTE: the direction of comparison is reversed relative to the C++ API", \
      python::return_value_policy<python::manage_new_object>());
#define QADEF2(_funcname_)                                               \
  python::def(#_funcname_ "QueryAtom", _funcname_##QueryAtom,            \
              (python::arg("negate") = false),                           \
              "Returns a QueryAtom that matches atoms when " #_funcname_ \
              " is True.",                                               \
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
    QADEF1(NumRadicalElectrons)
    QADEF1(NumHeteroatomNeighbors)
    QADEF1(NumAliphaticHeteroatomNeighbors)

    QADEF2(IsUnsaturated);
    QADEF2(IsAromatic);
    QADEF2(IsAliphatic);
    QADEF2(IsInRing);
    QADEF2(HasChiralTag);
    QADEF2(MissingChiralTag);

    python::def("HasPropQueryAtom", HasPropQueryAtom,
                (python::arg("propname"), python::arg("negate") = false),
                "Returns a QueryAtom that matches when the propery 'propname' "
                "exists in the atom.",
                python::return_value_policy<python::manage_new_object>());

    python::def("HasPropQueryBond", HasPropQueryBond,
                (python::arg("propname"), python::arg("negate") = false),
                "Returns a QueryBond that matches when the propery 'propname' "
                "exists in the bond.",
                python::return_value_policy<python::manage_new_object>());

    python::def("HasIntPropWithValueQueryAtom",
                PropQueryWithTol<Atom, QueryAtom, int>,
                (python::arg("propname"), python::arg("val"),
                 python::arg("negate") = false, python::arg("tolerance") = 0),
                "Returns a QueryAtom that matches when the propery 'propname' "
                "has the specified int value.",
                python::return_value_policy<python::manage_new_object>());

    python::def("HasBoolPropWithValueQueryAtom",
                PropQuery<Atom, QueryAtom, bool>,
                (python::arg("propname"), python::arg("val"),
                 python::arg("negate") = false),
                "Returns a QueryAtom that matches when the propery 'propname' "
                "has the specified boolean"
                " value.",
                python::return_value_policy<python::manage_new_object>());

    python::def("HasStringPropWithValueQueryAtom",
                PropQuery<Atom, QueryAtom, std::string>,
                (python::arg("propname"), python::arg("val"),
                 python::arg("negate") = false),
                "Returns a QueryAtom that matches when the propery 'propname' "
                "has the specified string "
                "value.",
                python::return_value_policy<python::manage_new_object>());

    python::def("HasDoublePropWithValueQueryAtom",
                PropQueryWithTol<Atom, QueryAtom, double>,
                (python::arg("propname"), python::arg("val"),
                 python::arg("negate") = false, python::arg("tolerance") = 0.0),
                "Returns a QueryAtom that matches when the propery 'propname' "
                "has the specified "
                "value +- tolerance",
                python::return_value_policy<python::manage_new_object>());

    /////////////////////////////////////////////////////////////////////////////////////
    //  Bond Queries
    python::def("HasPropQueryBond", HasPropQueryBond,
                (python::arg("propname"), python::arg("negate") = false),
                "Returns a QueryBond that matches when the propery 'propname' "
                "exists in the bond.",
                python::return_value_policy<python::manage_new_object>());

    python::def("HasPropQueryBond", HasPropQueryBond,
                (python::arg("propname"), python::arg("negate") = false),
                "Returns a QueryBond that matches when the propery 'propname' "
                "exists in the bond.",
                python::return_value_policy<python::manage_new_object>());

    python::def("HasIntPropWithValueQueryBond",
                PropQueryWithTol<Bond, QueryBond, int>,
                (python::arg("propname"), python::arg("val"),
                 python::arg("negate") = false, python::arg("tolerance") = 0),
                "Returns a QueryBond that matches when the propery 'propname' "
                "has the specified int value.",
                python::return_value_policy<python::manage_new_object>());

    python::def("HasBoolPropWithValueQueryBond",
                PropQuery<Bond, QueryBond, bool>,
                (python::arg("propname"), python::arg("val"),
                 python::arg("negate") = false),
                "Returns a QueryBond that matches when the propery 'propname' "
                "has the specified boolean"
                " value.",
                python::return_value_policy<python::manage_new_object>());

    python::def("HasStringPropWithValueQueryBond",
                PropQuery<Bond, QueryBond, std::string>,
                (python::arg("propname"), python::arg("val"),
                 python::arg("negate") = false),
                "Returns a QueryBond that matches when the propery 'propname' "
                "has the specified string "
                "value.",
                python::return_value_policy<python::manage_new_object>());

    python::def("HasDoublePropWithValueQueryBond",
                PropQueryWithTol<Bond, QueryBond, double>,
                (python::arg("propname"), python::arg("val"),
                 python::arg("negate") = false, python::arg("tolerance") = 0.0),
                "Returns a QueryBond that matches when the propery 'propname' "
                "has the specified "
                "value +- tolerance",
                python::return_value_policy<python::manage_new_object>());
  };
};
}  // end of namespace

void wrap_queries() { RDKit::queries_wrapper::wrap(); }
