// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_MOL_FILE_STEREOCHEM_H_
#define _RD_MOL_FILE_STEREOCHEM_H_

#include <GraphMol/RDKitBase.h>

namespace RDKit {
  void DetectAtomStereoChemistry(RWMol &mol, const Conformer *conf);
  void DetectBondStereoChemistry(ROMol &mol, const Conformer *conf);
  void WedgeMolBonds(ROMol &mol, const Conformer *conf);
  INT_MAP_INT pickBondsToWedge(const ROMol &mol);
  void ClearSingleBondDirFlags(ROMol &mol);
  Bond::BondDir DetermineBondWedgeState(const Bond *bond,
                                        const INT_MAP_INT &wedgeBonds,
                                        const Conformer *conf);

}
#endif
