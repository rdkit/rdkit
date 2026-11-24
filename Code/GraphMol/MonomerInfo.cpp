//
//  Copyright (C) 2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <ostream>
#include "MonomerInfo.h"

using namespace RDKit;

AtomMonomerInfo *AtomMonomerInfo::copy() const {
  if (d_monomerType == PDBRESIDUE) {
    return new AtomPDBResidueInfo(d_name, 0, "", d_residueName, d_residueNumber,
                                  d_chainId, "", 1.0, 0.0, false, 0, 0,
                                  d_monomerClass);
  }
  return new AtomMonomerInfo(*this);
}

//! allows AtomPDBResidueInfo objects to be dumped to streams
std::ostream &operator<<(std::ostream &target, const AtomPDBResidueInfo &apri) {
  target << apri.getSerialNumber() << " " << apri.getName() << " "
         << apri.getResidueName() << " " << apri.getChainId() << " "
         << apri.getResidueNumber();
  return target;
}
