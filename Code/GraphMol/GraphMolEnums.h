//
//  Copyright (C) 2025 NVIDIA Corporation & Affiliates and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include <RDGeneral/export.h>

namespace RDKit {

namespace AtomEnums {
//! store hybridization
enum HybridizationType {
  UNSPECIFIED = 0,  //!< hybridization that hasn't been specified
  S,
  SP,
  SP2,
  SP3,
  SP2D,
  SP3D,
  SP3D2,
  OTHER  //!< unrecognized hybridization
};

//! store type of chirality
enum ChiralType {
  CHI_UNSPECIFIED = 0,      //!< chirality that hasn't been specified
  CHI_TETRAHEDRAL_CW,       //!< tetrahedral: clockwise rotation (SMILES \@\@)
  CHI_TETRAHEDRAL_CCW,      //!< tetrahedral: counter-clockwise rotation (SMILES \@)
  CHI_OTHER,                //!< some unrecognized type of chirality
  CHI_TETRAHEDRAL,          //!< tetrahedral, use permutation flag
  CHI_ALLENE,               //!< allene, use permutation flag
  CHI_SQUAREPLANAR,         //!< square planar, use permutation flag
  CHI_TRIGONALBIPYRAMIDAL,  //!< trigonal bipyramidal, use permutation flag
  CHI_OCTAHEDRAL            //!< octahedral, use permutation flag
};
}  // namespace AtomEnums

namespace BondEnums {
//! the type of Bond
enum BondType {
  UNSPECIFIED = 0,
  SINGLE,
  DOUBLE,
  TRIPLE,
  QUADRUPLE,
  QUINTUPLE,
  HEXTUPLE,
  ONEANDAHALF,
  TWOANDAHALF,
  THREEANDAHALF,
  FOURANDAHALF,
  FIVEANDAHALF,
  AROMATIC,
  IONIC,
  HYDROGEN,
  THREECENTER,
  DATIVEONE,  //!< one-electron dative (e.g. from a C in a Cp ring to a metal)
  DATIVE,     //!< standard two-electron dative
  DATIVEL,    //!< standard two-electron dative
  DATIVER,    //!< standard two-electron dative
  OTHER,
  ZERO  //!< Zero-order bond (from
  // http://pubs.acs.org/doi/abs/10.1021/ci200488k)
};

//! the bond's direction (for chirality)
enum BondDir {
  NONE = 0,    //!< no special style
  BEGINWEDGE,  //!< wedged: narrow at begin
  BEGINDASH,   //!< dashed: narrow at begin
  // FIX: this may not really be adequate
  ENDDOWNRIGHT,  //!< for cis/trans
  ENDUPRIGHT,    //!<  ditto
  EITHERDOUBLE,  //!< a "crossed" double bond
  UNKNOWN,       //!< intentionally unspecified stereochemistry
};

//! the nature of the bond's stereochem (for cis/trans)
enum BondStereo {  // stereochemistry of double bonds
  STEREONONE = 0,  // no special style
  STEREOANY,       // intentionally unspecified
  // -- Put any true specifications about this point so
  // that we can do comparisons like if(bond->getStereo()>BondStereo::STEREOANY)
  STEREOZ,         // Z double bond
  STEREOE,         // E double bond
  STEREOCIS,       // cis double bond
  STEREOTRANS,     // trans double bond
  STEREOATROPCW,   //  atropisomer clockwise rotation
  STEREOATROPCCW,  //  atropisomer counter clockwise rotation
};
}  // namespace BondEnums

//! returns true if the atom is to the left of C
RDKIT_GRAPHMOL_EXPORT bool isEarlyAtom(int atomicNum);

}  // namespace RDKit
