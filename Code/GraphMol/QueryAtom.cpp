
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
#include <boost/algorithm/string/predicate.hpp>
#include <limits>

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
bool localMatch(ATOM_EQUALS_QUERY const *q1, ATOM_EQUALS_QUERY const *q2) {
  if (q1->getNegation() == q2->getNegation()) {
    return q1->getVal() == q2->getVal();
  } else {
    return q1->getVal() != q2->getVal();
  }
}

bool getAtomInterval(const QueryAtom::QUERYATOM_QUERY *q, int &lo, int &hi) {
  auto addSafe = [](int a, int b) {
    if (b > 0 && a > std::numeric_limits<int>::max() - b)
      return std::numeric_limits<int>::max();
    return a + b;
  };
  auto subSafe = [](int a, int b) {
    if (b > 0 && a < std::numeric_limits<int>::min() + b)
      return std::numeric_limits<int>::min();
    return a - b;
  };

  auto *r = dynamic_cast<const ATOM_RANGE_QUERY *>(q);
  if (r) {
    int tol = r->getTol();
    lo = subSafe(r->getLower(), tol);
    hi = addSafe(r->getUpper(), tol);
    auto ends = r->getEndsOpen();
    if (ends.first) {
      if (lo >= std::numeric_limits<int>::max()) return false;
      lo++;
    }
    if (ends.second) {
      if (hi <= std::numeric_limits<int>::min()) return false;
      hi--;
    }
    if (lo > hi) return false;
    return true;
  }
  auto *le = dynamic_cast<const ATOM_LESSEQUAL_QUERY *>(q);
  if (le) {
    int tol = le->getTol();
    lo = subSafe(le->getVal(), tol);
    hi = std::numeric_limits<int>::max();
    return true;
  }
  auto *ge = dynamic_cast<const ATOM_GREATEREQUAL_QUERY *>(q);
  if (ge) {
    int tol = ge->getTol();
    lo = std::numeric_limits<int>::min();
    hi = addSafe(ge->getVal(), tol);
    return true;
  }
  auto *gt = dynamic_cast<const ATOM_GREATER_QUERY *>(q);
  if (gt) {
    int tol = gt->getTol();
    if (gt->getVal() <= std::numeric_limits<int>::min() + tol) return false;
    lo = std::numeric_limits<int>::min();
    hi = subSafe(gt->getVal(), tol) - 1;
    return true;
  }
  auto *lt = dynamic_cast<const ATOM_LESS_QUERY *>(q);
  if (lt) {
    int tol = lt->getTol();
    if (lt->getVal() >= std::numeric_limits<int>::max() - tol) return false;
    lo = addSafe(lt->getVal(), tol) + 1;
    hi = std::numeric_limits<int>::max();
    return true;
  }
  auto *e = dynamic_cast<const ATOM_EQUALS_QUERY *>(q);
  if (e) {
    int tol = e->getTol();
    lo = subSafe(e->getVal(), tol);
    hi = addSafe(e->getVal(), tol);
    return true;
  }
  return false;
}

bool queriesMatch(QueryAtom::QUERYATOM_QUERY const *q1,
                  QueryAtom::QUERYATOM_QUERY const *q2) {
  PRECONDITION(q1, "no q1");
  PRECONDITION(q2, "no q2");

  static const unsigned int nQueries = 20;
  static std::string equalityQueries[nQueries] = {"AtomType",
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

  bool res = false;
  std::string d1 = q1->getDescription();
  std::string d2 = q2->getDescription();
  if (d1 == "AtomNull" || d2 == "AtomNull") {
    res = true;
  } else if (d1 == "AtomOr") {
    // FIX: handle negation on AtomOr and AtomAnd
    for (auto iter1 = q1->beginChildren(); iter1 != q1->endChildren();
         ++iter1) {
      if (d2 == "AtomOr") {
        for (auto iter2 = q2->beginChildren(); iter2 != q2->endChildren();
             ++iter2) {
          if (queriesMatch(iter1->get(), iter2->get())) {
            res = true;
            break;
          }
        }
      } else {
        if (queriesMatch(iter1->get(), q2)) {
          res = true;
        }
      }
      if (res) {
        break;
      }
    }
  } else if (d1 == "AtomAnd") {
    res = true;
    for (auto iter1 = q1->beginChildren(); iter1 != q1->endChildren();
         ++iter1) {
      bool matched = false;
      if (d2 == "AtomAnd") {
        for (auto iter2 = q2->beginChildren(); iter2 != q2->endChildren();
             ++iter2) {
          if (queriesMatch(iter1->get(), iter2->get())) {
            matched = true;
            break;
          }
        }
      } else {
        matched = queriesMatch(iter1->get(), q2);
      }
      if (!matched) {
        res = false;
        break;
      }
    }
    // FIX : handle AtomXOr
  } else if (d2 == "AtomOr") {
    // FIX: handle negation on AtomOr and AtomAnd
    for (auto iter2 = q2->beginChildren(); iter2 != q2->endChildren();
         ++iter2) {
      if (queriesMatch(q1, iter2->get())) {
        res = true;
        break;
      }
    }
  } else if (d2 == "AtomAnd") {
    res = true;
    for (auto iter2 = q2->beginChildren(); iter2 != q2->endChildren();
         ++iter2) {
      if (!queriesMatch(q1, iter2->get())) {
        res = false;
        break;
      }
    }
  } else {
    // Strip prefix to get base description
    auto stripPrefix = [](const std::string &s) -> std::string {
      if (boost::starts_with(s, "range_")) return s.substr(6);
      if (boost::starts_with(s, "less_")) return s.substr(5);
      if (boost::starts_with(s, "greater_")) return s.substr(8);
      return s;
    };
    std::string bd1 = stripPrefix(d1);
    std::string bd2 = stripPrefix(d2);
    if (bd1 == bd2 && std::find(&equalityQueries[0], &equalityQueries[nQueries],
                                bd1) != &equalityQueries[nQueries]) {
      int lo1 = 0, hi1 = 0, lo2 = 0, hi2 = 0;
      bool hasI1 = getAtomInterval(q1, lo1, hi1);
      bool hasI2 = getAtomInterval(q2, lo2, hi2);
      if (hasI1 && hasI2) {
        if (q1->getNegation() == q2->getNegation()) {
          if (!q1->getNegation()) {
            // Both positive: pattern interval must be subset of target interval
            res = lo2 <= lo1 && hi1 <= hi2;
          } else {
            // Both negated: accepted sets are complements; pattern accepted
            // set ⊆ target accepted set iff target excluded range ⊆ pattern
            // excluded range
            res = lo1 <= lo2 && hi2 <= hi1;
          }
        } else {
          if (q1->getNegation()) {
            // Pattern negated, target positive: complement of pattern's
            // excluded range must fit inside target. Complement has up to
            // two pieces: [MIN, lo1-1] and [hi1+1, MAX].
            {
              bool resLeft = true, resRight = true;
              if (lo1 > std::numeric_limits<int>::min()) {
                resLeft =
                    lo2 <= std::numeric_limits<int>::min() && lo1 - 1 <= hi2;
              }
              if (hi1 < std::numeric_limits<int>::max()) {
                resRight =
                    lo2 <= hi1 + 1 && std::numeric_limits<int>::max() <= hi2;
              }
              res = resLeft && resRight;
            }
          } else {
            // Pattern positive, target negated: pattern interval must lie
            // entirely outside target's excluded range
            res = hi1 < lo2 || lo1 > hi2;
          }
        }
      } else {
        // One or both queries have empty interval (e.g. GreaterQuery(INT_MIN)
        // or LessQuery(INT_MAX)). Empty interval means the query matches no
        // values; its negation matches all values.
        if (!hasI1 && !hasI2) {
          // Both empty: pattern accepted set ⊆ target accepted set
          // Empty ⊆ empty → true; All ⊆ All → true; Empty ⊆ All → true;
          // All ⊆ Empty → false
          if (q1->getNegation() && !q2->getNegation()) {
            res = false;
          } else {
            res = true;
          }
        } else if (!hasI1) {
          // Pattern is NONE or ALL; target is a proper subset → only NONE
          // matches.
          res = !q1->getNegation();
        } else {
          // Target is NONE or ALL; pattern is nonempty → only ALL matches.
          res = q2->getNegation();
        }
      }
    }
  }
  return res;
}
}  // namespace

bool QueryAtom::Match(Atom const *what) const {
  PRECONDITION(what, "bad query atom");
  PRECONDITION(dp_query, "no query set");
  return dp_query->Match(what);
}
bool QueryAtom::QueryMatch(QueryAtom const *what) const {
  PRECONDITION(what, "bad query atom");
  PRECONDITION(dp_query, "no query set");
  if (!what->hasQuery()) {
    return dp_query->Match(what);
  } else {
    return queriesMatch(dp_query, what->getQuery());
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
