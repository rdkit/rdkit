// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
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
