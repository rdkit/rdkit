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

// TODO as details    // typedef struct InternalCoordinates InternalCoordinates;
typedef struct InternalCoordinates {
  TorsionInfo torsionRange;
  BondAngleInfo angles;
  BondLengthInfo lengths;
  InternalCoordinates(unsigned int numBonds)
      : torsionRange{}, angles{}, lengths(numBonds) {}
} InternalCoordinates;

/*!
  \param TODO
*/
RDKIT_DISTGEOMHELPERS_EXPORT void setMoleculeDFS(
    const ROMol &mol, std::shared_ptr<DistGeom::ZMatrix> zmat,
    const InternalCoordinates &internalCoords);

//! TODO
/*!
  \param TODO
*/
RDKIT_DISTGEOMHELPERS_EXPORT void setMoleculeDFS(
    const ROMol &mol, std::shared_ptr<DistGeom::ZMatrix> zmat,
    const InternalCoordinates &internalCoords, const unsigned int startAtomIdx,
    const unsigned int atomIdx2, const unsigned int atomIdx3);

//! TODO
/*!
  \param TODO
*/
RDKIT_DISTGEOMHELPERS_EXPORT void correctChiralCenters(
    const ROMol &mol, std::shared_ptr<DistGeom::ZMatrix> zmat);

}  // namespace DGeomHelpers
}  // namespace RDKit
#endif
