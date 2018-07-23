//
// Copyright (c) 2003-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef __RD_GREATERQUERY_H__
#define __RD_GREATERQUERY_H__
#include "Query.h"
#include "EqualityQuery.h"

namespace Queries {
//! \brief a Query implementing > using a particular
//!  value (and an optional tolerance)
template <class MatchFuncArgType, class DataFuncArgType = MatchFuncArgType,
          bool needsConversion = false>
class GreaterQuery
    : public EqualityQuery<MatchFuncArgType, DataFuncArgType, needsConversion> {
 public:
  GreaterQuery() { this->d_tol = 0; };
  //! constructs with our target value
  explicit GreaterQuery(DataFuncArgType v) {
    this->d_val = v;
    this->d_tol = 0;
    this->df_negate = false;
  };
  //! constructs with our target value and a tolerance
  GreaterQuery(DataFuncArgType v, DataFuncArgType t) {
    this->d_val = v;
    this->d_tol = t;
    this->df_negate = false;
  };

  bool Match(const DataFuncArgType what) const {
    MatchFuncArgType mfArg =
        this->TypeConvert(what, Int2Type<needsConversion>());
    if (queryCmp(this->d_val, mfArg, this->d_tol) > 0) {
      if (this->getNegation())
        return false;
      else
        return true;
    } else {
      if (this->getNegation())
        return true;
      else
        return false;
    }
  };

  Query<MatchFuncArgType, DataFuncArgType, needsConversion> *copy() const {
    GreaterQuery<MatchFuncArgType, DataFuncArgType, needsConversion> *res =
        new GreaterQuery<MatchFuncArgType, DataFuncArgType, needsConversion>();
    res->setVal(this->d_val);
    res->setTol(this->d_tol);
    res->setNegation(this->getNegation());
    res->setDataFunc(this->d_dataFunc);
    res->d_description = this->d_description;
    return res;
  };

  std::string getFullDescription() const {
    std::ostringstream res;
    res << this->getDescription();
    res << " " << this->d_val;
    if (this->getNegation())
      res << " ! > ";
    else
      res << " > ";
    return res.str();
  };
};
}
#endif
