//
//   Copyright (C) 2005-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <Geometry/Transform3D.h>
#include <Geometry/UniformGrid3D.h>

#include "ShapeEncoder.h"
#include <Geometry/Transform3D.h>
#include <Geometry/point.h>
#include <Geometry/UniformGrid3D.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/RDKitBase.h>
namespace RDKit {
namespace MolShapes {
void EncodeShape(const ROMol &mol, RDGeom::UniformGrid3D &grid, int confId,
                 const RDGeom::Transform3D *trans, double vdwScale,
                 double stepSize, int maxLayers, bool ignoreHs) {
  const auto &conf = mol.getConformer(confId);
  EncodeShape(conf, grid, trans, vdwScale, stepSize, maxLayers, ignoreHs);
}

void EncodeShape(const Conformer &conf, RDGeom::UniformGrid3D &grid,
                 const RDGeom::Transform3D *trans, double vdwScale,
                 double stepSize, int maxLayers, bool ignoreHs) {
  const auto &mol = conf.getOwningMol();
  for (const auto &atom : mol.atoms()) {
    auto anum = atom->getAtomicNum();
    if (anum == 1 && ignoreHs) {  // ignore hydrigens
      continue;
    }
    auto rad = PeriodicTable::getTable()->getRvdw(anum);
    RDGeom::Point3D loc = conf.getAtomPos(atom->getIdx());
    if (trans) {
      trans->TransformPoint(loc);
    }
    grid.setSphereOccupancy(loc, vdwScale * rad, stepSize, maxLayers);
  }
}
}  // namespace MolShapes
}  // namespace RDKit
