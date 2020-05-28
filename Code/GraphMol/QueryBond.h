//
//  Copyright (C) 2001-2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_QUERYBOND_H
#define _RD_QUERYBOND_H

#include <Query/QueryObjects.h>
#include "Bond.h"
#include "QueryOps.h"

namespace RDKit {

//! Class for storing Bond queries
/*!
  QueryBond objects are derived from Bond objects, so they can be
  added to molecules and the like, but they have much fancier
  querying capabilities.

 */

class RDKIT_GRAPHMOL_EXPORT QueryBond : public Bond {
 public:
  typedef Queries::Query<int, Bond const *, true> QUERYBOND_QUERY;

  QueryBond() : Bond(){};
  //! initialize with a particular bond order
  explicit QueryBond(BondType bT);
  //! initialize from a bond
  explicit QueryBond(const Bond &other)
      : Bond(other), dp_query(makeBondOrderEqualsQuery(other.getBondType())){};
  QueryBond(const QueryBond &other)
      : Bond(other), dp_query(other.dp_query->copy()){};

  ~QueryBond();

  //! returns a copy of this query, owned by the caller
  virtual Bond *copy() const;

  QueryBond &operator=(const QueryBond &other);

  //! sets the BondType of this query:
  void setBondType(BondType bT);
  //! sets the BondDir of this query:
  void setBondDir(BondDir bD);

  //! returns true if we match Bond \c what
  bool Match(Bond const *what) const;

  //! returns true if our query details match those of QueryBond \c what
  bool QueryMatch(QueryBond const *what) const;

  // This method can be used to distinguish query bonds from standard bonds
  bool hasQuery() const { return dp_query != nullptr; };

  //! returns our current query
  QUERYBOND_QUERY *getQuery() const { return dp_query; };
  //! replaces our current query with the value passed in
  void setQuery(QUERYBOND_QUERY *what) {
    // free up any existing query (Issue255):
    delete dp_query;
    dp_query = what;
  };

  //! expands our current query
  /*!
    \param what          the Queries::Query to be added. The ownership of
                         the query is passed to the current object, where it
                         might be deleted, so that the pointer should not be
                         used again in the calling code.
    \param how           the operator to be used in the expansion
    \param maintainOrder (optional) flags whether the relative order of
                         the queries needs to be maintained, if this is
                         false, the order is reversed

    <b>Notes:</b>
      - \c what should probably be constructed using one of the functions
         defined in QueryOps.h
      - the \c maintainOrder option can be useful because the combination
        operators short circuit when possible.

  */
  void expandQuery(QUERYBOND_QUERY *what,
                   Queries::CompositeQueryType how = Queries::COMPOSITE_AND,
                   bool maintainOrder = true);

 protected:
  QUERYBOND_QUERY *dp_query{nullptr};
};

namespace detail {
inline std::string qhelper(Bond::QUERYBOND_QUERY *q, unsigned int depth) {
  std::string res;
  if (q) {
    for (unsigned int i = 0; i < depth; ++i) res += "  ";
    res += q->getFullDescription() + "\n";
    for (Bond::QUERYBOND_QUERY::CHILD_VECT_CI ci = q->beginChildren();
         ci != q->endChildren(); ++ci) {
      res += qhelper((*ci).get(), depth + 1);
    }
  }
  return res;
}
}  // namespace detail
inline std::string describeQuery(const Bond *bond) {
  PRECONDITION(bond, "bad bond");
  std::string res = "";
  if (bond->hasQuery()) {
    res = detail::qhelper(bond->getQuery(), 0);
  }
  return res;
}
};  // namespace RDKit

#endif
