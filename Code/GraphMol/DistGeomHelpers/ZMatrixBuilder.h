//
//  Copyright (C) 2026 ETH Zurich
//  Created by: Katharina Buchthal
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_Z_MATRIX_BUILDER_H
#define RD_Z_MATRIX_BUILDER_H

#include <DistGeom/ZMatrix.h>
#include <boost/unordered/unordered_flat_map.hpp>

namespace RDKit {
class ROMol;
namespace DGeomHelpers {

using TorsionInfo =
    boost::unordered_flat_map<std::size_t, DistGeom::TorsionCandidates>;
using BondAngleInfo = boost::unordered_flat_map<std::size_t, double>;
using BondLengthInfo = std::vector<double>;

//! A struct to store internal coordinates by unified ids through the ids from
//! the invovled bonds
struct InternalCoordinates {
  TorsionInfo torsionRange;
  BondAngleInfo angles;
  BondLengthInfo lengths;
  InternalCoordinates(unsigned int numBonds)
      : torsionRange{}, angles{}, lengths(numBonds) {}
};

//! Generates a Z-Matrix via a DFS walk through the molecule, as starting point,
//! a non-ring bond (if existing) with minimal degree is selected
/*!
  \param mol  Molecule for which the matrix should be constructes
  \param zmat Pointer to an initialized (!) but empty Z-Matrix that should be
  filled
  \param internalCoords Internal coordinates to be used for the ZMatrix (NOTE:
  the molecule that was given here must be the same used to construct the
  internal coordinates, since they are accessed by their bond ids)
*/
RDKIT_DISTGEOMHELPERS_EXPORT void setMoleculeDFS(
    const ROMol &mol, std::shared_ptr<DistGeom::ZMatrix> zmat,
    const InternalCoordinates &internalCoords);

//! Generates a Z-Matrix via a DFS walk through the molecule, as starting point,
//! the given two atoms are used
/*!
  \param mol  Molecule for which the matrix should be constructes
  \param zmat Pointer to an initialized (!) but empty Z-Matrix that should be
  filled
  \param internalCoords Internal coordinates to be used for the ZMatrix (NOTE:
  the molecule that was given here must be the same used to construct the
  internal coordinates, since they are accessed by their bond ids)
  \param firstAtomIdx First row in zmatrix
  \param secondAtomIdx Second row in zmatrix
*/
RDKIT_DISTGEOMHELPERS_EXPORT void setMoleculeDFS(
    const ROMol &mol, std::shared_ptr<DistGeom::ZMatrix> zmat,
    const InternalCoordinates &internalCoords, unsigned int firstAtomIdx,
    unsigned int secondAtomIdx);

//! Corrects sign for offsets of improper torsions to account for chiral
//! chemistry
/*!
  \param mol  Corresponding molecule
  \param zmat ZMatrix
*/
RDKIT_DISTGEOMHELPERS_EXPORT void correctChiralCenters(
    const ROMol &mol, std::shared_ptr<DistGeom::ZMatrix> zmat);

}  // namespace DGeomHelpers
}  // namespace RDKit
#endif
