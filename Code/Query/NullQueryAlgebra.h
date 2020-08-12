//
//
//  Copyright (C) 2019 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_NULLQUERYALGEBRA_H
#define RD_NULLQUERYALGEBRA_H

#include <GraphMol/QueryOps.h>

namespace RDKit {
namespace {
template <class T>
void mergeBothNullQ(T *&returnQuery, T *&otherNullQ,
                    Queries::CompositeQueryType how) {
  bool negatedQ = returnQuery->getNegation();
  bool negatedOtherQ = otherNullQ->getNegation();

  if (how == Queries::COMPOSITE_AND) {
    // This is the only case in which we need to do anything
    if (!negatedQ && negatedOtherQ) {
      returnQuery->setNegation(true);
    }
  } else if (how == Queries::COMPOSITE_OR) {
    // This is the only case in which we need to do anything
    if (negatedQ && !negatedOtherQ) {
      returnQuery->setNegation(false);
    }
  } else if (how == Queries::COMPOSITE_XOR) {
    if (!negatedQ && !negatedOtherQ) {
      returnQuery->setNegation(true);
    } else if (negatedQ + negatedOtherQ == 1) {
      returnQuery->setNegation(false);
    }
  }
}

template <class T>
void mergeNullQFirst(T *&returnQuery, T *&otherQ,
                     Queries::CompositeQueryType how) {
  bool negatedQ = returnQuery->getNegation();

  if (how == Queries::COMPOSITE_AND) {
    if (!negatedQ) {
      std::swap(returnQuery, otherQ);
    }
  } else if (how == Queries::COMPOSITE_OR) {
    if (negatedQ) {
      std::swap(returnQuery, otherQ);
    }
  } else if (how == Queries::COMPOSITE_XOR) {
    std::swap(returnQuery, otherQ);
    if (!negatedQ) {
      returnQuery->setNegation(!returnQuery->getNegation());
    }
  }
}

}  // namespace

template <class T>
void mergeNullQueries(T *&returnQuery, bool isQueryNull, T *&otherQuery,
                      bool isOtherQNull, Queries::CompositeQueryType how) {
  PRECONDITION(returnQuery, "bad query");
  PRECONDITION(otherQuery, "bad query");
  PRECONDITION(how == Queries::COMPOSITE_AND || how == Queries::COMPOSITE_OR ||
                   how == Queries::COMPOSITE_XOR,
               "bad combination op");

  if (isQueryNull && isOtherQNull) {
    mergeBothNullQ(returnQuery, otherQuery, how);
  } else if (isQueryNull) {
    mergeNullQFirst(returnQuery, otherQuery, how);
  } else if (isOtherQNull) {
    std::swap(returnQuery, otherQuery);
    mergeNullQFirst(returnQuery, otherQuery, how);
  }
}
}  // namespace RDKit

#endif