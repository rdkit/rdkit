//
// Copyright (c) 2003-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_ANDQUERY_H__
#define __RD_ANDQUERY_H__

#include "Query.h"

namespace Queries {

//! a Query implementing AND: requires all children to be \c true
template <class MatchFuncArgType, class DataFuncArgType = MatchFuncArgType,
          bool needsConversion = false>
class AndQuery
    : public Query<MatchFuncArgType, DataFuncArgType, needsConversion> {
 public:
  typedef Query<MatchFuncArgType, DataFuncArgType, needsConversion> BASE;
  AndQuery() { this->df_negate = false; };

  bool Match(const DataFuncArgType what) const {
    bool res = true;
    typename BASE::CHILD_VECT_CI it1;
    for (it1 = this->beginChildren(); it1 != this->endChildren(); ++it1) {
      bool tmp = (*it1)->Match(what);
      if (!tmp) {
        res = false;
        break;
      }
    }
    if (this->getNegation()) res = !res;
    return res;
  };
  Query<MatchFuncArgType, DataFuncArgType, needsConversion> *copy() const {
    AndQuery<MatchFuncArgType, DataFuncArgType, needsConversion> *res =
        new AndQuery<MatchFuncArgType, DataFuncArgType, needsConversion>();
    typename BASE::CHILD_VECT_CI i;
    for (i = this->beginChildren(); i != this->endChildren(); ++i) {
      res->addChild(typename BASE::CHILD_TYPE(i->get()->copy()));
    }
    res->setNegation(this->getNegation());
    res->d_description = this->d_description;
    return res;
  };
};
}
#endif
