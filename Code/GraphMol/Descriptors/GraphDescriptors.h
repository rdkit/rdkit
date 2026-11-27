//
//  Copyright (C) 2026 Glysade LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file GraphDescriptors.h

  \brief C++ ports of the topological descriptors from
  rdkit.Chem.GraphDescriptors 

*/
#include <RDGeneral/export.h>
#ifndef RD_GRAPHDESCRIPTORS_H
#define RD_GRAPHDESCRIPTORS_H

#include <string>

namespace RDKit {
class ROMol;
namespace Descriptors {

//! Chi0: From equations (1),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422,
//! (1991)
/*!
  Uses the heavy-atom degree of each atom; atoms with degree 0 are skipped.
*/
RDKIT_DESCRIPTORS_EXPORT double calcChi0(const ROMol &mol);
const std::string chi0Version = "1.0.0";

//! Chi1: From equations (1),(11) and (12) of Rev. Comp. Chem. vol 2, 367-422,
//! (1991)
RDKIT_DESCRIPTORS_EXPORT double calcChi1(const ROMol &mol);
const std::string chi1Version = "1.0.0";

//! Balaban's J value for a molecule
/*!
  Chem. Phys. Lett. vol 89, 399-404, (1982)

  \param mol     the molecule of interest
  \param distMat (optional) a row-major nAtoms x nAtoms distance matrix to
                 use instead of the bond-order weighted distance matrix of
                 the molecule
*/
RDKIT_DESCRIPTORS_EXPORT double calcBalabanJ(const ROMol &mol,
                                             const double *distMat = nullptr);
const std::string balabanJVersion = "1.0.0";

//! Bertz's complexity index
/*!
  From S. H. Bertz, J. Am. Chem. Soc., vol 103, 3599-3601 (1981)

  \param mol    the molecule of interest
  \param cutoff vertices are considered topologically identical if their
                sorted distance vectors are equal out to the cutoff-th
                nearest neighbor
  \param distMat (optional) a row-major nAtoms x nAtoms distance matrix to
                 use, instead of the bond-order weighted distance matrix of
                 the molecule, when assigning the symmetry classes
*/
RDKIT_DESCRIPTORS_EXPORT double calcBertzCT(const ROMol &mol,
                                            unsigned int cutoff = 100,
                                            const double *distMat = nullptr);
const std::string bertzCTVersion = "2.0.0";

//! Information content of the coefficients of the characteristic polynomial
//! of the adjacency matrix of a hydrogen-suppressed graph of a molecule.
/*!
  From Eq 6 of D. Bonchev & N. Trinajstic, J. Chem. Phys. vol 67, 4517-4533
  (1977)

  \param mol  the molecule of interest
  \param avg  if set, the information content divided by the total
              population (AvgIpc) is returned
  \param distMat (optional) a row-major nAtoms x nAtoms distance matrix;
                 the adjacency matrix is taken from the entries of its lower
                 triangle, diagonal included, that equal 1. By default the
                 bonds of the molecule are used.

  The coefficients of the characteristic polynomial are integers. They are
  computed exactly and each rounded to double once, so they do not depend on
  the platform or on an eigensolver; only the final entropy is evaluated in
  floating point. NaN is returned for graphs so large (about 3,000 atoms of a
  molecular graph) that the coefficients exceed the built-in table of primes;
  the double range is exceeded long before that.
*/
RDKIT_DESCRIPTORS_EXPORT double calcIpc(const ROMol &mol, bool avg = false,
                                        const double *distMat = nullptr);
const std::string ipcVersion = "2.0.0";

//! Average information content (Eq 7 of Bonchev & Trinajstic), see calcIpc()
RDKIT_DESCRIPTORS_EXPORT double calcAvgIpc(const ROMol &mol,
                                           const double *distMat = nullptr);
const std::string avgIpcVersion = "2.0.0";

}  // namespace Descriptors
}  // namespace RDKit
#endif
