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

#include <string>
#include <list>
#include <cmath>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <Geometry/point.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/export.h>
#include <boost/dynamic_bitset.hpp>

namespace RDKit {
namespace Descriptors {

//! Class for calculation of the Shrake and Rupley surface area and volume
//! using the Double Cubic Lattice Method.
//!
//! Frank Eisenhaber, Philip Lijnzaad, Patrick Argos, Chris Sander and
//! Michael Scharf, "The Double Cubic Lattice Method: Efficient Approaches
//! to Numerical Integration of Surface Area and Volume and to Dot Surface
//! Contouring of Molecular Assemblies", Journal of Computational Chemistry,
//! Vol. 16, No. 3, pp. 273-284, 1995.
class RDKIT_DESCRIPTORS_EXPORT DoubleCubicLatticeVolume {
 public:
  // default params assume a small molecule and default conformer
  const ROMol &mol;
  std::vector<double> radii_;
  bool isProtein = false;
  bool includeLigand = true;
  double probeRadius = 1.4;
  int confId = -1;
  double maxRadius = 1.7;  // treat default max radius as Carbon

  /*!

    \param mol: input molecule or protein
    \param radii: radii for atoms of input mol, empty for default values
    \param isProtein: flag to calculate burried surface area of a protein ligand
    complex [default=false, free ligand]
    \param includeLigand: flag to trigger
    inclusion of bound ligand in surface area and volume calculations where
    molecule is a protein [default=true]
    \param probeRadius: radius of the
    sphere representing the probe solvent atom  [default=1.4]
    \param confId: conformer ID to consider [default=-1]

  */
  DoubleCubicLatticeVolume(const ROMol &mol, std::vector<double> radii,
                           bool isProtein = false, bool includeLigand = true,
                           double probeRadius = 1.4, int confId = -1);

  //! \overload uses default vdw radii
  DoubleCubicLatticeVolume(const ROMol &mol, bool isProtein = false,
                           bool includeLigand = true, double probeRadius = 1.4,
                           int confId = -1)
      : DoubleCubicLatticeVolume(mol, std::vector<double>(), isProtein,
                                 includeLigand, probeRadius, confId) {};

  // value returns

  /*! \return Solvent Accessible Surface Area */
  double getSurfaceArea();

  /*! \return Polar Surface Area */
  double getPolarSurfaceArea(bool includeSandP, bool includeHs);

  /*! \return Surface Area from specified atoms */
  double getPartialSurfaceArea(const boost::dynamic_bitset<> &incAtoms);

  /*! \return Solvent Accessible Surface Area for specified atom */
  double getAtomSurfaceArea(unsigned int atomIdx);

  /*! \return Set of Points representing the surface */
  std::map<unsigned int, std::vector<RDGeom::Point3D>> &getSurfacePoints();

  /*! \return Volume bound by probe sphere */
  double getVolume();

  /*! \return van der Waals Volume */
  double getVDWVolume();

  /*! \return Polar Volume */
  double getPolarVolume(bool includeSandP, bool includeHs);

  /*! \return Volume from specified atoms */
  double getPartialVolume(const boost::dynamic_bitset<> &incAtoms);

  /*! \return Volume for specified atom */
  double getAtomVolume(unsigned int atomIdx, double solventRadius);

  /*! \return Compactness of the protein */
  double getCompactness();

  /*! \return Packing Density of the protein */
  double getPackingDensity();

 private:
  // used by methods
  unsigned int numAtoms = 0;
  std::vector<RDGeom::Point3D> positions;
  std::vector<std::vector<unsigned int>> neighbours;
  RDGeom::Point3D centreOfGravity;

  // outputs
  double surfaceArea = 0.0;
  double totalVolume = 0.0;
  double vdwVolume = 0.0;
  std::map<unsigned int, std::vector<RDGeom::Point3D>> surfacePoints;

  // helpers
  bool testPoint(const RDGeom::Point3D &vect, double solvrad,
                 const std::vector<unsigned int> &nbrs);
};
}  // namespace Descriptors
}  // namespace RDKit
#endif