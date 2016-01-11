//
// Copyright (C) 2015 Greg Landrum
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDTHREADS_H_2015
#define RDTHREADS_H_2015

#ifdef RDK_THREADSAFE_SSS
#include "RDGeneral/Invariant.h"
#include <RDGeneral/BoostStartInclude.h>
#include <boost/thread.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {
inline unsigned int getNumThreadsToUse(int target) {
  if (target >= 1) {
    return static_cast<unsigned int>(target);
  }
  unsigned int res = boost::thread::hardware_concurrency();
  if (res > rdcast<unsigned int>(-target)) {
    return res + target;
  } else {
    return 1;
  }
}
}

#else

namespace RDKit {
inline unsigned int getNumThreadsToUse(int target) {
  RDUNUSED_PARAM(target);
  return 1;
}
}
#endif

#endif
