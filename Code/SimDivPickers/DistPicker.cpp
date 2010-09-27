// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "DistPicker.h"

namespace RDPickers {
  double getDistFromLTM(const double *distMat, unsigned int i, unsigned int j){
    CHECK_INVARIANT(distMat, "");
    if (i == j) {
      return 0.0;
    } else if (i > j) {
      return distMat[i*(i-1)/2 + j];
    } else {
      return distMat[j*(j-1)/2 + i];
    }
  }
}
