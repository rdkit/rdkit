//
//  Copyright (C) 2003-2011 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_GASTEIGERCHARGES_H
#define _RD_GASTEIGERCHARGES_H

#include <vector>

namespace RDKit {
  class ROMol;
  void computeGasteigerCharges(const ROMol *mol, int nIter=12,
			       bool throwOnParamFailure=false);
  void computeGasteigerCharges(const ROMol &mol, int nIter=12,
			       bool throwOnParamFailure=false);
  void computeGasteigerCharges(const ROMol &mol,std::vector<double> &charges,
                               int nIter=12,bool throwOnParamFailure=false);
}

#endif

    
