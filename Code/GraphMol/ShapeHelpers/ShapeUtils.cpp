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
#include "ShapeUtils.h"
#include "ShapeEncoder.h"
#include <Geometry/UniformGrid3D.h>
#include <GraphMol/RDKitBase.h>
#include <Geometry/Transform3D.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <Geometry/GridUtils.h>

namespace RDKit {
  namespace MolShapes {

    void computeConfBox(const Conformer &conf, RDGeom::Point3D &leftBottom, 
                         RDGeom::Point3D &rightTop, const RDGeom::Transform3D *trans,
                         double padding) {
      double xmin, xmax, ymin, ymax, zmin, zmax;
      xmin = ymin = zmin = 1.e8;
      xmax = ymax = zmax = -1.e8;
      unsigned int i, nAtms = conf.getNumAtoms();
      for (i = 0; i < nAtms; ++i) {
        RDGeom::Point3D loc = conf.getAtomPos(i);
        if (trans) {
          trans->TransformPoint(loc);
        }
        xmax = std::max(xmax, loc.x);
        xmin = std::min(xmin, loc.x);
        ymax = std::max(ymax, loc.y);
        ymin = std::min(ymin, loc.y);
        zmax = std::max(zmax, loc.z);
        zmin = std::min(zmin, loc.z);
        
      }
      RDGeom::Point3D padPt(padding, padding, padding);
      leftBottom.x = xmin; leftBottom.y = ymin; leftBottom.z = zmin;
      rightTop.x = xmax; rightTop.y = ymax; rightTop.z = zmax;
      leftBottom -= padPt;
      rightTop += padPt;
    }
    
    void computeConfDimsAndOffset(const Conformer &conf, RDGeom::Point3D &dims, 
                                  RDGeom::Point3D &offSet, const RDGeom::Transform3D *trans, 
                                  double padding) {
      //RDGeom::Point3D lb, rb;
      computeConfBox(conf, offSet, dims, trans, padding);
      dims -= offSet;
    }

    std::vector<double> getConfDimensions(const Conformer &conf, double padding, 
                                          const RDGeom::Point3D *center, bool ignoreHs) {
      RDGeom::Point3D lb, rb;
      computeConfBox(conf, lb, rb, 0, padding);
      
      if (!center) {
        RDGeom::Point3D cpt = MolTransforms::computeCentroid(conf, ignoreHs);  
        rb -= cpt;
        lb -= cpt;
      } else {
        rb -= (*center);
        lb -= (*center);
      }
      lb *= -1.0;
      double dimX = 2.0*std::max(rb.x, lb.x);
      double dimY = 2.0*std::max(rb.y, lb.y);
      double dimZ = 2.0*std::max(rb.z, lb.z);
      std::vector<double> res;
      res.reserve(3);
      res.push_back(dimX);
      res.push_back(dimY);
      res.push_back(dimZ);
      return res;
    }

    void computeUnionBox(const RDGeom::Point3D &leftBottom1, const RDGeom::Point3D &rightTop1, 
                         const RDGeom::Point3D &leftBottom2, const RDGeom::Point3D &rightTop2, 
                         RDGeom::Point3D &uLeftBottom, RDGeom::Point3D &uRightTop) {
      uLeftBottom.x = std::min(leftBottom1.x, leftBottom2.x);
      uLeftBottom.y = std::min(leftBottom1.y, leftBottom2.y);
      uLeftBottom.z = std::min(leftBottom1.z, leftBottom2.z);

      uRightTop.x = std::max(rightTop1.x, rightTop2.x);
      uRightTop.y = std::max(rightTop1.y, rightTop2.y);
      uRightTop.z = std::max(rightTop1.z, rightTop2.z);
    }

    double tanimotoDistance(const ROMol &mol1, const ROMol &mol2, int confId1, int confId2,
                            double gridSpacing, DiscreteValueVect::DiscreteValueType bitsPerPoint,
                            double vdwScale, double stepSize, int maxLayers, bool ignoreHs) {
      const Conformer &conf1 = mol1.getConformer(confId1);
      const Conformer &conf2 = mol2.getConformer(confId2);
      return tanimotoDistance(conf1, conf2, gridSpacing=0.5, bitsPerPoint, vdwScale,
                       stepSize, maxLayers, ignoreHs);
    }

    double tanimotoDistance(const Conformer &conf1, const Conformer &conf2, double gridSpacing, 
                            DiscreteValueVect::DiscreteValueType bitsPerPoint, double vdwScale,
                            double stepSize, int maxLayers, bool ignoreHs) {
      RDGeom::Transform3D *trans = MolTransforms::computeCanonicalTransform(conf1);
      
      // now use this transform and figure out what size grid we will need
      // find the lower-left and upper-right corners for each of the conformers
      // and take a union of these boxes - we will use this fo grid dimensions
      RDGeom::Point3D leftBottom1, rightTop1, leftBottom2, rightTop2, uLeftBottom, uRightTop;
      computeConfBox(conf1, leftBottom1, rightTop1, trans);
      computeConfBox(conf2, leftBottom2, rightTop2, trans);
      
      computeUnionBox(leftBottom1, rightTop1, leftBottom2, rightTop2, uLeftBottom, uRightTop);
      
      // make the grid object to store the encoding
      uRightTop -= uLeftBottom; // uRightTop now has grid dimensions
      
      RDGeom::UniformGrid3D grd1(uRightTop.x, uRightTop.y, uRightTop.z, gridSpacing, bitsPerPoint,
                                 &uLeftBottom);
      RDGeom::UniformGrid3D grd2(uRightTop.x, uRightTop.y, uRightTop.z, gridSpacing, bitsPerPoint,
                                 &uLeftBottom);

      EncodeShape(conf1, grd1, trans, vdwScale, stepSize, maxLayers, ignoreHs);
      EncodeShape(conf2, grd2, trans, vdwScale, stepSize, maxLayers, ignoreHs);
      return RDGeom::tanimotoDistance(grd1, grd2);
    }


    double protrudeDistance(const ROMol &mol1, const ROMol &mol2, int confId1, int confId2,
                            double gridSpacing, DiscreteValueVect::DiscreteValueType bitsPerPoint,
                            double vdwScale, double stepSize, int maxLayers, bool ignoreHs,
                            bool allowReordering) {
      const Conformer &conf1 = mol1.getConformer(confId1);
      const Conformer &conf2 = mol2.getConformer(confId2);
      return protrudeDistance(conf1, conf2, gridSpacing=0.5, bitsPerPoint, vdwScale,
                       stepSize, maxLayers, ignoreHs, allowReordering);
    }
        
    double protrudeDistance(const Conformer &conf1, const Conformer &conf2, double gridSpacing, 
                            DiscreteValueVect::DiscreteValueType bitsPerPoint, double vdwScale,
                            double stepSize, int maxLayers, bool ignoreHs,
                            bool allowReordering) {
      //
      // FIX: all this duplicated code needs to be refactored out.
      //
      RDGeom::Transform3D *trans = MolTransforms::computeCanonicalTransform(conf1);
      
      // now use this transform and figure out what size grid we will need
      // find the lower-left and upper-right corners for each of the conformers
      // and take a union of these boxes - we will use this fo grid dimensions
      RDGeom::Point3D leftBottom1, rightTop1, leftBottom2, rightTop2, uLeftBottom, uRightTop;
      computeConfBox(conf1, leftBottom1, rightTop1, trans);
      computeConfBox(conf2, leftBottom2, rightTop2, trans);
      
      computeUnionBox(leftBottom1, rightTop1, leftBottom2, rightTop2, uLeftBottom, uRightTop);
      
      // make the grid object to store the encoding
      uRightTop -= uLeftBottom; // uRightTop now has grid dimensions
      
      RDGeom::UniformGrid3D grd1(uRightTop.x, uRightTop.y, uRightTop.z, gridSpacing, bitsPerPoint,
                                 &uLeftBottom);
      RDGeom::UniformGrid3D grd2(uRightTop.x, uRightTop.y, uRightTop.z, gridSpacing, bitsPerPoint,
                                 &uLeftBottom);

      EncodeShape(conf1, grd1, trans, vdwScale, stepSize, maxLayers, ignoreHs);
      EncodeShape(conf2, grd2, trans, vdwScale, stepSize, maxLayers, ignoreHs);
      double res;
      if(allowReordering &&
         ( grd2.getOccupancyVect()->getTotalVal() < grd1.getOccupancyVect()->getTotalVal() ) ){
        res=RDGeom::protrudeDistance(grd2, grd1);
      } else {
        res=RDGeom::protrudeDistance(grd1, grd2);
      }
      return res;
    }
  }
}
      
      

      
      
