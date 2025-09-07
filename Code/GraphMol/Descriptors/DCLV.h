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
    \param confId: conformer ID to consider [default=-1]

  */

  // default params assume a small molecule and default conformer
  const ROMol &mol;
  bool isProtein = false;
  bool includeLigand = true;
  double probeRadius = 1.4;
  int confId = -1;

  DoubleCubicLatticeVolume(const ROMol &mol, bool isProtein = false,
                           bool includeLigand = true, double probeRadius = 1.2,
                           int confId = -1);
  //! Class for calculation of the Shrake and Rupley surface area and volume
  //! using the Double Cubic Lattice Method.
  //!
  //! Frank Eisenhaber, Philip Lijnzaad, Patrick Argos, Chris Sander and
  //! Michael Scharf, "The Double Cubic Lattice Method: Efficient Approaches
  //! to Numerical Integration of Surface Area and Volume and to Dot Surface
  //! Contouring of Molecular Assemblies", Journal of Computational Chemistry,
  //! Vol. 16, No. 3, pp. 273-284, 1995.

  // value returns

  /*! \return Solvent Accessible Surface Area */
  double getSurfaceArea();

  /*! \return Polar Surface Area */
  double getPolarSurfaceArea(const bool &includeSandP);

  /*! \return Solvent Accessible Surface Area for specified atom */
  double getAtomSurfaceArea(const unsigned int &atom_idx);
  
  /*! \return Volume for specified atom */
  double getAtomVolume(const unsigned int &atom_idx, const double &solventRadius);
  
  /*! \return Volume bound by probe sphere */
  double getVolume();
  
  /*! \return van der Waals Volume */
  double getVDWVolume(); 
  
  /*! \return Compactness of the protein */
  double getCompactness();

  /*! \return Packing Density of the protein */
  double getPackingDensity();

  /*! \return Set of Points representing the surface */
  
  private:
  
  // used by methods
  unsigned int numAtoms = 0;
  std::vector<RDGeom::Point3D> positions;
  std::vector<double> radii;
  std::vector<std::vector<unsigned int>> neighbours;
  RDGeom::Point3D centreOfGravity;

  // outputs
  double surfaceArea = 0.0;
  double polarSurfaceArea = 0.0;
  double totalVolume = 0.0;
  double vdwVolume = 0.0;

  // helpers
  bool testPoint(double *vect, const double &solvrad, const std::vector<unsigned int> &nbrs);

};
}  // namespace Descriptors
}  // namespace RDKit
#endif