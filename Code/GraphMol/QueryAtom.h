//
//  Copyright (C) 2001-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_QUERYATOM_H
#define RD_QUERYATOM_H

#include <utility>
#include "Atom.h"
#include <Query/QueryObjects.h>
#include <GraphMol/QueryOps.h>

namespace RDKit {

//! Class for storing atomic queries
/*!
  QueryAtom objects are derived from Atom objects, so they can be
  added to molecules and the like, but they have much fancier
  querying capabilities.

 */
class RDKIT_GRAPHMOL_EXPORT QueryAtom : public Atom {
 public:
  typedef Queries::Query<int, Atom const *, true> QUERYATOM_QUERY;

  QueryAtom() : Atom() {}
  explicit QueryAtom(int num) : Atom(num), dp_query(makeAtomNumQuery(num)) {}
  explicit QueryAtom(const Atom &other)
      : Atom(other), dp_query(makeAtomNumQuery(other.getAtomicNum())) {
    if (other.getIsotope()) {
      this->expandQuery(makeAtomIsotopeQuery(other.getIsotope()),
                        Queries::CompositeQueryType::COMPOSITE_AND);
    }
    if (other.getFormalCharge()) {
      this->expandQuery(makeAtomFormalChargeQuery(other.getFormalCharge()),
                        Queries::CompositeQueryType::COMPOSITE_AND);
    }
    if (other.getNumRadicalElectrons()) {
      this->expandQuery(
          makeAtomNumRadicalElectronsQuery(other.getNumRadicalElectrons()),
          Queries::CompositeQueryType::COMPOSITE_AND);
    }
  }
  QueryAtom(const QueryAtom &other) : Atom(other) {
    if (other.dp_query) {
      dp_query = other.dp_query->copy();
    } else {
      dp_query = nullptr;
    }
  }
  QueryAtom &operator=(const QueryAtom &other) {
    if (this == &other) {
      return *this;
    }
    Atom::operator=(other);
    delete dp_query;
    if (other.dp_query) {
      dp_query = other.dp_query->copy();
    } else {
      dp_query = nullptr;
    }
    return *this;
  }

  QueryAtom(QueryAtom &&other) noexcept : Atom(std::move(other)) {
    dp_query = std::exchange(other.dp_query, nullptr);
  }
  QueryAtom &operator=(QueryAtom &&other) noexcept {
    if (this == &other) {
      return *this;
    }
    QueryAtom::operator=(std::move(other));
    dp_query = std::exchange(other.dp_query, nullptr);
    return *this;
  }

  ~QueryAtom() override;

  //! returns a copy of this query, owned by the caller
  Atom *copy() const override;

  // This method can be used to distinguish query atoms from standard atoms:
  bool hasQuery() const override { return dp_query != nullptr; }

  //! returns the label associated to this query
  std::string getQueryType() const override { return dp_query->getTypeLabel(); }

  //! replaces our current query with the value passed in
  void setQuery(QUERYATOM_QUERY *what) override {
    delete dp_query;
    dp_query = what;
  }
  //! returns our current query
  QUERYATOM_QUERY *getQuery() const override { return dp_query; }

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
  void expandQuery(QUERYATOM_QUERY *what,
                   Queries::CompositeQueryType how = Queries::COMPOSITE_AND,
                   bool maintainOrder = true) override;

  //! returns true if we match Atom \c what
  bool Match(Atom const *what) const override;

  //! returns true if our query details match those of QueryAtom \c what
  bool QueryMatch(QueryAtom const *what) const;

 private:
  QUERYATOM_QUERY *dp_query{nullptr};

};  // end o' class

namespace detail {
inline std::string qhelper(const Atom::QUERYATOM_QUERY *q, unsigned int depth) {
  std::string res = "";
  if (q) {
    for (unsigned int i = 0; i < depth; ++i) {
      res += "  ";
    }
    res += q->getFullDescription() + "\n";
    for (const auto &child :
         boost::make_iterator_range(q->beginChildren(), q->endChildren())) {
      res += qhelper(child.get(), depth + 1);
    }
  }
  return res;
}
}  // namespace detail
inline std::string describeQuery(const Atom *atom) {
  PRECONDITION(atom, "bad atom");
  std::string res = "";
  if (atom->hasQuery()) {
    res = detail::qhelper(atom->getQuery(), 0);
  }
  return res;
}

};  // namespace RDKit

#endif
