//
//  Copyright (C) 2020 Brian P. Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDKIT_BCUT_H
#define RDKIT_BCUT_H
#ifdef RDK_HAS_EIGEN3
#include <RDGeneral/export.h>
#include <vector>
#include <string>
namespace RDKit {
class ROMol;
namespace Descriptors {
const std::string BCUT2DVersion = "1.0.0";

//! Implements BCUT descriptors From J. Chem. Inf. Comput. Sci., Vol. 39, No. 1, 1999
/*!
   \param mol         the molecule of interest
   \param atom_props  A vector of double for the atom properties to use for the diagonal elements
                      of the BCUT matrix.  Must be equal to he number of atoms in mol.
   \return            std::pair<double,double> pair.first is the high eval, pair.second the lowest
*/  
RDKIT_DESCRIPTORS_EXPORT
std::pair<double,double> BCUT2D(const ROMol &mol, const std::vector<double> &atom_props);

//! Implements BCUT descriptors From J. Chem. Inf. Comput. Sci., Vol. 39, No. 1, 1999
/*!
   \param mol         the molecule of interest
   \param atom_props  An atom property holding a double value on the atom
                      If atom_propname does not exist on the atom, raises
                      If atom_propname cannot be coerced into a double, raises
   \return            std::pair<double,double> pair.first is the high eval, pair.second the lowest
*/  
RDKIT_DESCRIPTORS_EXPORT
std::pair<double,double> BCUT2D(const ROMol &mol, const std::string &atom_propname);

//! Implements BCUT descriptors From J. Chem. Inf. Comput. Sci., Vol. 39, No. 1, 1999
//!  Diagonal elements are (currently) atomic mass, gasteiger charge,
//!  crippen logP and crippen MR
/*!
  \param mol           the molecule of interest

  \return              The bcut vector (mass high eval, mass low eval,
                                        gasteiger charge high eval, gasteiger charge low eval,
                                        LogP high eval, LogP low eval,
                                        MR high eval, MR low eval)
*/
RDKIT_DESCRIPTORS_EXPORT    
std::vector<double> BCUT2D(const ROMol &m);
}
}
 
#endif
#endif
