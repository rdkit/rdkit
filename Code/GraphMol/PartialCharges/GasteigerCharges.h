//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_GASTEIGERCHARGES_H
#define _RD_GASTEIGERCHARGES_H

namespace RDKit {
  class ROMol;
  void computeGasteigerCharges(const ROMol *mol, int nIter=12,
			       bool throwOnParamFailure=false);
}

#endif

    
