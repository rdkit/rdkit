//
// Copyright (c) 2003-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
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
