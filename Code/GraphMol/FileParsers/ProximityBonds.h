//
//  Copyright (C) 2013 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_PROXIMITYBONDS_H_
#define _RD_PROXIMITYBONDS_H_
#include <GraphMol/RWMol.h>

namespace RDKit {
  void ConnectTheDots(RWMol *mol);
  void StandardPDBResidueBondOrders(RWMol *mol);
}

#endif  // _RD_PROXIMITYBONDS_H_

