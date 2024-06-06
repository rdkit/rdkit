/*=================================================================*/
/* Copyright (C)  2024  Greg Landrum and other RDKit contributors  */
/* Contributed by NextMove Software, Cambridge, UK.                */
/*                                                                 */
/*                                                                 */
/* @@ All Rights Reserved @@                                       */
/* The contents are covered by the terms of the                    */
/* BSD license, which is included in the file                      */
/* license.txt.                                                    */
/*=================================================================*/

#ifndef RDKIT_DCLV_H
#define RDKIT_DCLV_H

#include <iostream>
#include <string>
#include <list>
#include <cmath>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <Geometry/point.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/export.h>

namespace RDKit {
namespace Descriptors {

class RDKIT_DESCRIPTORS_EXPORT DoubleCubicLatticeVolume {
 public:
  /*!

    \param mol: input molecule or protein
    \param isProtein: flag to calculate burried surface area of a protein ligand 
    complex [default=false, free ligand]
    \param includeLigand: flag to trigger
    inclusion of bound ligand in surface area and volume calculations where
    molecule is a protein [default=true]
    \param probeRadius: radius of the
    sphere representing the probe solvent atom
    \param depth: controls the number
    of dots per atom
    \param dotDensity: controls density of dots per atom

  */
  DoubleCubicLatticeVolume(const ROMol& mol, bool isProtein = false,
                           bool includeLigand = true, double probeRadius = 1.2,
                           int depth = 4, int dotDensity = 0);
  //! Class for calculation of the Shrake and Rupley surface area and volume
  //! using the Double Cubic Lattice Method.
  //!
  //! Frank Eisenhaber, Philip Lijnzaad, Patrick Argos, Chris Sander and
  //! Michael Scharf, "The Double Cubic Lattice Method: Efficient Approaches
  //! to Numerical Integration of Surface Area and Volume and to Dot Surface
  //! Contouring of Molecular Assemblies", Journal of Computational Chemistry,
  //! Vol. 16, No. 3, pp. 273-284, 1995.

  // value returns
  double getSurfaceArea() {
    /*! \return Solvent Accessible Surface Area */
    return surfaceArea;
  }

  double getVolume() {
    /*! \return Volume bound by probe sphere */
    return totalVolume;
  }

  double getVDWVolume() { /*! \return van der Waals Volume */
    return vdwVolume;
  }

  double getCompactness() {
    /*! \return Compactness of the protein */
    return compactness;
  }

  double getPackingDensity() {
    /*! \return Packing Density of the protein */
    return packingDensity;
  }

 private:
  // outputs
  double surfaceArea;
  double totalVolume;
  double vdwVolume;
  double compactness;
  double packingDensity;
};
}  // namespace Descriptors
}  // namespace RDKit
#endif