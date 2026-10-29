
//  Copyright (C) 2001-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryOps.h>
#include <Query/NullQueryAlgebra.h>

namespace RDKit {

QueryAtom::~QueryAtom() {
  delete dp_query;
  dp_query = nullptr;
};

Atom *QueryAtom::copy() const {
  auto *res = new QueryAtom(*this);
  return static_cast<Atom *>(res);
}

void QueryAtom::expandQuery(QUERYATOM_QUERY *what,
                            Queries::CompositeQueryType how,
                            bool maintainOrder) {
  PRECONDITION(dp_query, "Can't expand empty query");
  bool thisIsNullQuery = dp_query->getDescription() == "AtomNull";
  bool otherIsNullQuery = what->getDescription() == "AtomNull";

  if (thisIsNullQuery || otherIsNullQuery) {
    mergeNullQueries(dp_query, thisIsNullQuery, what, otherIsNullQuery, how);
    delete what;
    return;
  }

  QUERYATOM_QUERY *origQ = dp_query;
  std::string descrip;
  switch (how) {
    case Queries::COMPOSITE_AND:
      dp_query = new ATOM_AND_QUERY;
      descrip = "AtomAnd";
      break;
    case Queries::COMPOSITE_OR:
      dp_query = new ATOM_OR_QUERY;
      descrip = "AtomOr";
      break;
    case Queries::COMPOSITE_XOR:
      dp_query = new ATOM_XOR_QUERY;
      descrip = "AtomXor";
      break;
    default:
      UNDER_CONSTRUCTION("unrecognized combination query");
  }
  dp_query->setDescription(descrip);
  if (maintainOrder) {
    dp_query->addChild(QUERYATOM_QUERY::CHILD_TYPE(origQ));
    dp_query->addChild(QUERYATOM_QUERY::CHILD_TYPE(what));
  } else {
    dp_query->addChild(QUERYATOM_QUERY::CHILD_TYPE(what));
    dp_query->addChild(QUERYATOM_QUERY::CHILD_TYPE(origQ));
  }
}

namespace {
struct AtomQueryTraits {
  using BaseQuery = QueryAtom::QUERYATOM_QUERY;
  using RangeQuery = ATOM_RANGE_QUERY;
  using EqQuery = ATOM_EQUALS_QUERY;
  static constexpr std::string_view NullDesc = "AtomNull";
  static constexpr std::string_view OrDesc = "AtomOr";
  static constexpr std::string_view AndDesc = "AtomAnd";
  static bool isEqualityQuery(std::string_view s) {
    static const std::string queries[] = {"AtomType",
                                          "AtomRingBondCount",
                                          "AtomRingSize",
                                          "AtomMinRingSize",
                                          "AtomImplicitValence",
                                          "AtomExplicitValence",
                                          "AtomTotalValence",
                                          "AtomAtomicNum",
                                          "AtomExplicitDegree",
                                          "AtomTotalDegree",
                                          "AtomHCount",
                                          "AtomIsAromatic",
                                          "AtomIsAliphatic",
                                          "AtomUnsaturated",
                                          "AtomMass",
                                          "AtomFormalCharge",
                                          "AtomNegativeFormalCharge",
                                          "AtomHybridization",
                                          "AtomInRing",
                                          "AtomInNRings"};
    return std::find(std::begin(queries), std::end(queries), s) !=
           std::end(queries);
  }
};
}  // namespace

bool QueryAtom::Match(Atom const *what) const {
  PRECONDITION(what, "bad query atom");
  PRECONDITION(dp_query, "no query set");
  std::cerr << "match: " << dp_query->getDescription() << " vs "
            << what->hasQuery() << std::endl;
  return dp_query->Match(what);
}
bool QueryAtom::QueryMatch(QueryAtom const *what) const {
  PRECONDITION(what, "bad query atom");
  PRECONDITION(dp_query, "no query set");
  if (!what->hasQuery()) {
    return dp_query->Match(what);
  } else {
    return QueryOps::queriesMatch<AtomQueryTraits>(dp_query, what->getQuery());
  }
}

namespace detail {
bool hasRecursiveQuery(const QueryAtom::QUERYATOM_QUERY *q,
                       bool checkInitialized) {
  if (!q) {
    return false;
  }
  if (q->getDescription() == "RecursiveStructure" &&
      (!checkInitialized ||
       !static_cast<RecursiveStructureQuery const *>(q)->getInitialized())) {
    return true;
  }
  for (auto iter = q->beginChildren(); iter != q->endChildren(); ++iter) {
    if (hasRecursiveQuery(iter->get(), checkInitialized)) {
      return true;
    }
  }
  return false;
}
}  // namespace detail

bool hasUninitializedRecursiveQuery(const Atom &atom) {
  if (!atom.hasQuery()) {
    return false;
  }
  return detail::hasRecursiveQuery(atom.getQuery(), true);
}
bool hasRecursiveQuery(const Atom &atom) {
  if (!atom.hasQuery()) {
    return false;
  }
  return detail::hasRecursiveQuery(atom.getQuery(), false);
}

}  // namespace RDKit
