// $Id$
// 
//   Copyright (C) 2005-2006 Rational Discovery LLC
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
//#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
//#include <GraphMol/PeriodicTable.h>
#include <GraphMol/RDKitBase.h>
namespace RDKit {
  namespace MolShapes {
    void EncodeShape(const ROMol &mol, RDGeom::UniformGrid3D &grid, int confId, 
                     const RDGeom::Transform3D *trans,
                     double vdwScale, double stepSize, int maxLayers, bool ignoreHs) {
      const Conformer &conf = mol.getConformer(confId);
      EncodeShape(conf, grid, trans, vdwScale, stepSize, maxLayers, ignoreHs);
    }

    void EncodeShape(const Conformer &conf, RDGeom::UniformGrid3D &grid, 
                     const RDGeom::Transform3D *trans, double vdwScale, 
                     double stepSize, int maxLayers, bool ignoreHs) {
      const ROMol &mol = conf.getOwningMol();
      ROMol::ConstAtomIterator ai;
      double rad;
      unsigned int aid, anum;
      for (ai = mol.beginAtoms(); ai != mol.endAtoms(); ai++) {
        anum = (*ai)->getAtomicNum();
        if ((anum == 1) && (ignoreHs)) { //ignore hydrigens
          continue;
        }
        aid = (*ai)->getIdx();
        RDGeom::Point3D loc = conf.getAtomPos(aid);
        rad = PeriodicTable::getTable()->getRvdw(anum);
        if (trans) {
          trans->TransformPoint(loc);
        }
        grid.setSphereOccupancy(loc, vdwScale*rad, stepSize, maxLayers);
      }
    }
  }
}        
        
      
      
