/*==============================================*/
/* Copyright (C)  2019       NextMove Software  */
/* All rights reserved.                         */
/*                                              */
/* This file is part of molhash.                */
/*                                              */
/* The contents are covered by the terms of the */
/* BSD license, which is included in the file   */
/* license.txt.                                 */
/*==============================================*/

#include "toolkit.h"
#include "molhash.h"

void Strip(NMS_pMOL mol, unsigned int striptype)
{
  // The order of these operations is significant to some degree
  // - Hydrogens should be at the end as that is the common
  //   use case

  if (striptype & StripType::AtomStereo) {
    NMS_FOR_ATOM_IN_MOL(atom, mol) {
      NMS_pATOM aptr = NMS_ITER_MOL_ATOM(atom, mol);
      if (NMS_ATOM_HAS_STEREO(aptr))
        NMS_ATOM_REMOVE_STEREO(aptr);
    }
  }
  if (striptype & StripType::BondStereo) {
    NMS_FOR_BOND_IN_MOL(bond, mol) {
      NMS_pBOND bptr = NMS_ITER_MOL_BOND(bond, mol);
      if (NMS_BOND_HAS_STEREO(bptr))
        NMS_BOND_REMOVE_STEREO(bptr);
    }
  }
  if (striptype & StripType::Isotope) {
    NMS_FOR_ATOM_IN_MOL(atom, mol) {
      NMS_pATOM aptr = NMS_ITER_MOL_ATOM(atom, mol);
      NMS_ATOM_SET_ISOTOPE(aptr, 0);
    }
  }
  if (striptype & StripType::AtomMap) {
    NMS_FOR_ATOM_IN_MOL(atom, mol) {
      NMS_pATOM aptr = NMS_ITER_MOL_ATOM(atom, mol);
      NMS_ATOM_SET_MAPIDX(aptr, 0);
    }
  }
  if (striptype & StripType::Hydrogen) {
    NMS_MOL_SUPPRESS_HYDROGENS(mol);
  }
}

void SplitMolecule(NMS_pMOL mol, std::vector<NMS_MOL*> &molv)
{
  NMS_MOL_SPLIT_INTO_FRAGMENTS(mol, molv);
}
