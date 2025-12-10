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

namespace RDKit {

QueryBond::QueryBond(BondType bT) : Bond(bT) {
  if (bT != Bond::UNSPECIFIED) {
    dp_query.reset(makeBondOrderEqualsQuery(bT));
  } else {
    dp_query.reset(makeBondNullQuery());
  }
};

QueryBond &QueryBond::operator=(const QueryBond &other) {
  if (this == &other) {
    return *this;
  }

  // FIX: should we copy molecule ownership?  I don't think so.
  // FIX: how to deal with atom indices?
  // FIXME: The previous implementation only set dp_mol to null,
  // copied the bond type, copied the query, and copied the props.
  // Should this implementation also skip copying the other bond members?
  *(static_cast<Bond *>(this)) = other;

  setQuery(other.dp_query ? other.dp_query->copy() : nullptr);

  return *this;
}

Bond *QueryBond::copy() const {
  auto *res = new QueryBond(*this);
  return res;
}

void QueryBond::setBondType(BondType bT) {
  // NOTE: calling this blows out any existing query

  Bond::setBondType(bT);

  dp_query.reset(makeBondOrderEqualsQuery(bT));
}

void QueryBond::setBondDir(BondDir bD) {
  // NOTE: calling this blows out any existing query
  //
  //   Ignoring bond orders (which this implicitly does by blowing out
  //   any bond order query) is ok for organic molecules, where the
  //   only bonds assigned directions are single.  It'll fail in other
  //   situations, whatever those may be.
  //
  Bond::setBondDir(bD);
#if 0
  dp_query.reset(makeBondDirEqualsQuery(bD));
#endif
}

void QueryBond::setQuery(QUERYBOND_QUERY *what) {
  dp_query.reset(what);
  getDataRDMol().setBondQueryCompat(getIdx(), what);
  PRECONDITION(getDataRDMol().getBondQuery(getIdx()) != nullptr,
               "Missing new query");
}

void QueryBond::expandQuery(QUERYBOND_QUERY *what,
                            Queries::CompositeQueryType how,
                            bool maintainOrder) {
  std::unique_ptr<QUERYBOND_QUERY> whatPtr(what);
  QueryOps::expandQuery(dp_query, std::move(whatPtr), how, maintainOrder);
  getDataRDMol().setBondQueryCompat(getIdx(), dp_query.get());
}

namespace {
bool localMatch(BOND_EQUALS_QUERY const *q1, BOND_EQUALS_QUERY const *q2) {
  if (q1->getNegation() == q2->getNegation()) {
    return q1->getVal() == q2->getVal();
  } else {
    return q1->getVal() != q2->getVal();
  }
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
  } else if (std::find(&equalityQueries[0], &equalityQueries[nQueries], d1) !=
             &equalityQueries[nQueries]) {
    res = localMatch(static_cast<BOND_EQUALS_QUERY const *>(q1),
                     static_cast<BOND_EQUALS_QUERY const *>(q2));
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
    return queriesMatch(dp_query.get(), what->getQuery());
  }
}

double QueryBond::getValenceContrib(const Atom *atom) const {
  if (!hasQuery() || !QueryOps::hasComplexBondTypeQuery(*getQuery())) {
    return Bond::getValenceContrib(atom);
  }
  return 0;
}
}  // namespace RDKit
