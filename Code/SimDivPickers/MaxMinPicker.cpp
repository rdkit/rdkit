// $Id: MaxMinPicker.cpp 4959 2006-02-17 23:53:31Z glandrum $
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include "MaxMinPicker.h"
#include <RDGeneral/Invariant.h>
#include <RDGeneral/types.h>
#include <cstdlib>

namespace RDPickers {
  RDKit::INT_VECT MaxMinPicker::pick(const double *distMat, 
                                     unsigned int poolSize, unsigned int pickSize) {
    CHECK_INVARIANT(distMat, "Invalid Distance Matrix");
    CHECK_INVARIANT((poolSize >= pickSize), "pickSize cannot be larger than the poolSize");

    RDKit::INT_LIST pool;
    RDKit::INT_LIST_I plri, pli;
    

    RDKit::INT_VECT picks;
    RDKit::INT_VECT_CI pi;
    picks.reserve(pickSize);
    unsigned int i, j;
    unsigned int pick;

    // enter the pool into a list so that we can pick out of it easily
    for (i = 0; i < poolSize; i++) {
      pool.push_back(i);
    }

    // pick the first entry
    pick = rand()%poolSize;
    // add the pick to the picks
    picks.push_back(pick);
    // and remove it from the pool
    pool.remove(pick);
    double minTOi, maxOFmin, dist;
    // now pick 1 compound at a time
    while (picks.size() < pickSize) {
      maxOFmin = 0.0;
      for (pli = pool.begin(); pli != pool.end(); pli++) {
        i = (*pli);
        minTOi = RDKit::MAX_DOUBLE;
        for (pi = picks.begin(); pi != picks.end(); pi++) {
          j = (*pi);
	  CHECK_INVARIANT(i!=j,"");
          dist = getDistFromLTM(distMat, i, j);
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
      
      // now add teh new pick to  picks and remove it from the pool
      picks.push_back(pick);
      pool.erase(plri);
    }
    return picks;
  }
}
