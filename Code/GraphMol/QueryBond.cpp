//
//  Copyright (C) 2001-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/QueryBond.h>
#include <Query/NullQueryAlgebra.h>
#include <boost/algorithm/string/predicate.hpp>
#include <limits>

namespace RDKit {

QueryBond::QueryBond(BondType bT) : Bond(bT) {
  if (bT != Bond::UNSPECIFIED) {
    dp_query = makeBondOrderEqualsQuery(bT);
  } else {
    dp_query = makeBondNullQuery();
  }
};

QueryBond::~QueryBond() {
  delete dp_query;
  dp_query = nullptr;
};

QueryBond &QueryBond::operator=(const QueryBond &other) {
  // FIX: should we copy molecule ownership?  I don't think so.
  // FIX: how to deal with atom indices?
  dp_mol = nullptr;
  d_bondType = other.d_bondType;
  if (other.dp_query) {
    dp_query = other.dp_query->copy();
  } else {
    dp_query = nullptr;
  }
  d_props = other.d_props;
  return *this;
}

Bond *QueryBond::copy() const {
  auto *res = new QueryBond(*this);
  return res;
}

void QueryBond::setBondType(BondType bT) {
  // NOTE: calling this blows out any existing query
  d_bondType = bT;
  delete dp_query;
  dp_query = nullptr;

  dp_query = makeBondOrderEqualsQuery(bT);
}

void QueryBond::setBondDir(BondDir bD) {
  // NOTE: calling this blows out any existing query
  //
  //   Ignoring bond orders (which this implicitly does by blowing out
  //   any bond order query) is ok for organic molecules, where the
  //   only bonds assigned directions are single.  It'll fail in other
  //   situations, whatever those may be.
  //
  d_dirTag = bD;
}

void QueryBond::expandQuery(QUERYBOND_QUERY *what,
                            Queries::CompositeQueryType how,
                            bool maintainOrder) {
  bool thisIsNullQuery = dp_query->getDescription() == "BondNull";
  bool otherIsNullQuery = what->getDescription() == "BondNull";

  if (thisIsNullQuery || otherIsNullQuery) {
    mergeNullQueries(dp_query, thisIsNullQuery, what, otherIsNullQuery, how);
    delete what;
    return;
  }

  QUERYBOND_QUERY *origQ = dp_query;
  std::string descrip;
  switch (how) {
    case Queries::COMPOSITE_AND:
      dp_query = new BOND_AND_QUERY;
      descrip = "BondAnd";
      break;
    case Queries::COMPOSITE_OR:
      dp_query = new BOND_OR_QUERY;
      descrip = "BondOr";
      break;
    case Queries::COMPOSITE_XOR:
      dp_query = new BOND_XOR_QUERY;
      descrip = "BondXor";
      break;
    default:
      UNDER_CONSTRUCTION("unrecognized combination query");
  }
  dp_query->setDescription(descrip);
  if (maintainOrder) {
    dp_query->addChild(QUERYBOND_QUERY::CHILD_TYPE(origQ));
    dp_query->addChild(QUERYBOND_QUERY::CHILD_TYPE(what));
  } else {
    dp_query->addChild(QUERYBOND_QUERY::CHILD_TYPE(what));
    dp_query->addChild(QUERYBOND_QUERY::CHILD_TYPE(origQ));
  }
}

namespace {
bool localMatch(BOND_EQUALS_QUERY const *q1, BOND_EQUALS_QUERY const *q2) {
  if (q1->getNegation() == q2->getNegation()) {
    return q1->getVal() == q2->getVal();
  } else {
    return q1->getVal() != q2->getVal();
  }
}

// Extract interval [lo, hi] from a bond query. Returns false if the query
// type doesn't support interval extraction. GreaterEqualQuery(N) means
// N >= bond → [MIN,N]. LessEqualQuery(N) means N <= bond → [N,MAX].
bool getBondInterval(const QueryBond::QUERYBOND_QUERY *q, int &lo, int &hi) {
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

  auto *r = dynamic_cast<const BOND_RANGE_QUERY *>(q);
  if (r) {
    long long tol = r->getTol();
    auto ends = r->getEndsOpen();
    long long l = (long long)r->getLower() + (ends.first ? tol + 1 : -tol);
    long long h = (long long)r->getUpper() + (ends.second ? -tol - 1 : tol);
    if (l > h || l > (long long)std::numeric_limits<int>::max() ||
        h < (long long)std::numeric_limits<int>::min())
      return false;
    lo = (int)std::max(l, (long long)std::numeric_limits<int>::min());
    hi = (int)std::min(h, (long long)std::numeric_limits<int>::max());
    return true;
  }
  auto *le = dynamic_cast<const BOND_LESSEQUAL_QUERY *>(q);
  if (le) {
    int tol = le->getTol();
    lo = subSafe(le->getVal(), tol);
    hi = std::numeric_limits<int>::max();
    return true;
  }
  auto *ge = dynamic_cast<const BOND_GREATEREQUAL_QUERY *>(q);
  if (ge) {
    int tol = ge->getTol();
    lo = std::numeric_limits<int>::min();
    hi = addSafe(ge->getVal(), tol);
    return true;
  }
  auto *gt = dynamic_cast<const BOND_GREATER_QUERY *>(q);
  if (gt) {
    int tol = gt->getTol();
    if (gt->getVal() <= std::numeric_limits<int>::min() + tol) return false;
    lo = std::numeric_limits<int>::min();
    hi = subSafe(gt->getVal(), tol) - 1;
    return true;
  }
  auto *lt = dynamic_cast<const BOND_LESS_QUERY *>(q);
  if (lt) {
    int tol = lt->getTol();
    if (lt->getVal() >= std::numeric_limits<int>::max() - tol) return false;
    lo = addSafe(lt->getVal(), tol) + 1;
    hi = std::numeric_limits<int>::max();
    return true;
  }
  auto *e = dynamic_cast<const BOND_EQUALS_QUERY *>(q);
  if (e) {
    int tol = e->getTol();
    lo = subSafe(e->getVal(), tol);
    hi = addSafe(e->getVal(), tol);
    return true;
  }
  return false;
}

bool queriesMatch(QueryBond::QUERYBOND_QUERY const *q1,
                  QueryBond::QUERYBOND_QUERY const *q2) {
  PRECONDITION(q1, "no q1");
  PRECONDITION(q2, "no q2");

  static const unsigned int nQueries = 6;
  static std::string equalityQueries[nQueries] = {
      "BondRingSize", "BondMinRingSize", "BondOrder",
      "BondDir",      "BondInRing",      "BondInNRings"};

  bool res = false;
  std::string d1 = q1->getDescription();
  std::string d2 = q2->getDescription();
  if (d1 == "BondNull" || d2 == "BondNull") {
    res = true;
  } else if (d1 == "BondOr") {
    // FIX: handle negation on BondOr and BondAnd
    for (auto iter1 = q1->beginChildren(); iter1 != q1->endChildren();
         ++iter1) {
      if (d2 == "BondOr") {
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
  } else if (d1 == "BondAnd") {
    res = true;
    for (auto iter1 = q1->beginChildren(); iter1 != q1->endChildren();
         ++iter1) {
      bool matched = false;
      if (d2 == "BondAnd") {
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
    // FIX : handle BondXOr
  } else if (d2 == "BondOr") {
    // FIX: handle negation on BondOr and BondAnd
    for (auto iter2 = q2->beginChildren(); iter2 != q2->endChildren();
         ++iter2) {
      if (queriesMatch(q1, iter2->get())) {
        res = true;
        break;
      }
    }
  } else if (d2 == "BondAnd") {
    res = true;
    for (auto iter2 = q2->beginChildren(); iter2 != q2->endChildren();
         ++iter2) {
      if (queriesMatch(q1, iter2->get())) {
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
      bool hasI1 = getBondInterval(q1, lo1, hi1);
      bool hasI2 = getBondInterval(q2, lo2, hi2);
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

bool QueryBond::Match(Bond const *what) const {
  PRECONDITION(what, "bad query bond");
  PRECONDITION(dp_query, "no query set");
  return dp_query->Match(what);
}
bool QueryBond::QueryMatch(QueryBond const *what) const {
  PRECONDITION(what, "bad query bond");
  PRECONDITION(dp_query, "no query set");
  if (!what->hasQuery()) {
    return dp_query->Match(what);
  } else {
    return queriesMatch(dp_query, what->getQuery());
  }
}

double QueryBond::getValenceContrib(const Atom *atom) const {
  if (!hasQuery() || !QueryOps::hasComplexBondTypeQuery(*getQuery())) {
    return Bond::getValenceContrib(atom);
  }
  return 0;
}
}  // namespace RDKit
