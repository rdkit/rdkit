//
// Copyright (c) 2003-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
/*!
  \file QueryObjects.h

  \brief Pulls in all the query types.
*/  
  
#ifndef __RD_QUERYOBJECTS_H__
#define __RD_QUERYOBJECTS_H__

#include "Query.h"

#include "EqualityQuery.h"
#include "LessQuery.h"
#include "LessEqualQuery.h"
#include "GreaterQuery.h"
#include "GreaterEqualQuery.h"

#include "RangeQuery.h"
#include "SetQuery.h"

#include "AndQuery.h"
#include "OrQuery.h"
#include "XOrQuery.h"

namespace Queries {
  typedef enum {
    COMPOSITE_AND,
    COMPOSITE_OR,
    COMPOSITE_XOR
  } CompositeQueryType;
}
#endif
