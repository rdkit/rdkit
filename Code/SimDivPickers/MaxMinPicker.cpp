// $Id$
//
//  Copyright (C) 2003-2007 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include "MaxMinPicker.h"
#include <RDGeneral/Invariant.h>
#include <RDGeneral/types.h>
#include <RDBoost/Exceptions.h>
#include <cstdlib>

namespace RDPickers {
  RDKit::INT_VECT MaxMinPicker::pick(const double *distMat, 
                                     unsigned int poolSize, unsigned int pickSize) const {
    CHECK_INVARIANT(distMat, "Invalid Distance Matrix");
    if(poolSize<pickSize)
      throw ValueErrorException("pickSize cannot be larger than the poolSize");

    RDKit::INT_VECT picks;
    picks.reserve(pickSize);

    // enter the pool into a list so that we can pick out of it easily
    RDKit::INT_LIST pool;
    for (unsigned int i = 0; i < poolSize; i++) {
      pool.push_back(i);
    }

    // pick a random entry
    unsigned int pick = rand()%poolSize;
    // add the pick to the picks
    picks.push_back(pick);
    // and remove it from the pool
    pool.remove(pick);

    // now pick 1 compound at a time
    while (picks.size() < pickSize) {
      RDKit::INT_LIST_I plri;
      double maxOFmin = 0.0;
      for (RDKit::INT_LIST_I pli = pool.begin(); pli != pool.end(); ++pli) {
        unsigned int i = (*pli);
        double minTOi = RDKit::MAX_DOUBLE;
        for (RDKit::INT_VECT_CI pi = picks.begin(); pi != picks.end(); ++pi) {
          unsigned int j = (*pi);
          CHECK_INVARIANT(i!=j,"");
          double dist = getDistFromLTM(distMat, i, j);
          if (dist < minTOi) {
            minTOi = dist;
          }
        }
        if (minTOi > maxOFmin) {
          maxOFmin = minTOi;
          pick = i;
          plri = pli;
        }
      }
      // add the new pick to picks and remove it from the pool
      picks.push_back(pick);
      pool.erase(plri);
    }
    return picks;
  }
}
