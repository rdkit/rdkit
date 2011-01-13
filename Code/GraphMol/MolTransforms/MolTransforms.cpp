//  $Id$
// 
//   Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MolTransforms.h"
#include <GraphMol/RDKitBase.h>
#include <Numerics/EigenSolvers/PowerEigenSolver.h>
#include <Numerics/SymmMatrix.h>
#include <Numerics/Matrix.h>
#include <Geometry/Transform3D.h>

#define EIGEN_TOLERANCE 1.0e-2
namespace MolTransforms {
  
  using namespace RDKit;
  void transformAtom(Atom *atom,RDGeom::Transform3D &tform){
    PRECONDITION(atom,"no atom");
    ROMol &mol = atom->getOwningMol();
    for (ROMol::ConstConformerIterator ci = mol.beginConformers();
	 ci != mol.endConformers(); ci++) {
      RDGeom::Point3D &pos = (*ci)->getAtomPos(atom->getIdx());
      tform.TransformPoint(pos);
    }
    //atom->setPos(pos);
  }
  void transformMolsAtoms(ROMol *mol,RDGeom::Transform3D &tform){
    PRECONDITION(mol,"no molecule");

    ROMol::AtomIterator atomIt;
    for(atomIt=mol->beginAtoms();atomIt!=mol->endAtoms();atomIt++){
      transformAtom(*atomIt,tform);
    }
  }

  RDGeom::Point3D computeCentroid(const Conformer &conf, bool ignoreHs) {
    RDGeom::Point3D res(0.0, 0.0, 0.0);
    const ROMol &mol = conf.getOwningMol();
    ROMol::ConstAtomIterator cai;
    unsigned int nAtms = 0;
    
    for (cai = mol.beginAtoms(); cai != mol.endAtoms(); cai++) {
      if (((*cai)->getAtomicNum() == 1) && (ignoreHs)) {
        continue;
      }
      res += conf.getAtomPos((*cai)->getIdx());
      nAtms++;
    }
    res /= nAtms;
    return res;
  }

  RDNumeric::DoubleSymmMatrix *computeCovarianceMatrix(const Conformer &conf, 
                                                       const RDGeom::Point3D &center,
                                                       bool normalize, bool ignoreHs) {
    double xx, xy, xz, yy, yz, zz;
    xx = xy = xz = yy = yz = zz = 0.0;
    const ROMol &mol = conf.getOwningMol();
    ROMol::ConstAtomIterator cai;
    unsigned int nAtms = 0;
    for (cai = mol.beginAtoms(); cai != mol.endAtoms(); cai++) {
      if (((*cai)->getAtomicNum() == 1) && (ignoreHs) ) {
        continue;
      }
      RDGeom::Point3D loc = conf.getAtomPos((*cai)->getIdx());
      loc -= center;
      xx += loc.x*loc.x;
      xy += loc.x*loc.y;
      xz += loc.x*loc.z;
      yy += loc.y*loc.y;
      yz += loc.y*loc.z;
      zz += loc.z*loc.z;
      nAtms++;
    }
    if (normalize) {
      xx /= nAtms;
      xy /= nAtms;
      xz /= nAtms;
      yy /= nAtms;
      yz /= nAtms;
      zz /= nAtms;
    }
    RDNumeric::DoubleSymmMatrix *res = new RDNumeric::DoubleSymmMatrix(3,3);
    res->setVal(0,0, xx);
    res->setVal(0,1, xy);
    res->setVal(0,2, xz);
    res->setVal(1,1, yy);
    res->setVal(1,2, yz);
    res->setVal(2,2, zz);
    return res;
  }

  RDGeom::Transform3D *computeCanonicalTransform(const Conformer &conf,
                                                 const RDGeom::Point3D *center,
                                                 bool normalizeCovar,
                                                 bool ignoreHs) {
    RDGeom::Point3D origin;
    if (!center) {
      origin = computeCentroid(conf, ignoreHs);
    } else {
      origin = (*center);
    }
    RDNumeric::DoubleSymmMatrix *covMat = computeCovarianceMatrix(conf, origin, 
                                                                  normalizeCovar, ignoreHs);
    // find the eigen values and eigen vectors for the covMat
    RDNumeric::DoubleMatrix eigVecs(3,3);
    RDNumeric::DoubleVector eigVals(3);
    // if we have a single atom system we don't need to do anyhting other than setting translation
    // translation
    unsigned int nAtms = conf.getNumAtoms();
    RDGeom::Transform3D *trans = new RDGeom::Transform3D;
    
    // set the translation
    origin *= -1.0;
    //trans->SetTranslation(origin);
    // if we have a single atom system we don't need to do anyhting setting translation is sufficient
    if (nAtms > 1) {
      RDNumeric::EigenSolvers::powerEigenSolver(3, *covMat, eigVals, eigVecs,
                                                 conf.getNumAtoms());
      // deal with zero eigen value systems
      unsigned int i, j, dim = 3;
      for (i = 0; i < 3; ++i) {
        if (fabs(eigVals.getVal(i)) < EIGEN_TOLERANCE) {
          dim--;
        }
      }
      CHECK_INVARIANT(dim >= 1, "");
      if (dim < 3) {
        RDGeom::Point3D first(eigVecs.getVal(0,0), eigVecs.getVal(0,1), eigVecs.getVal(0,2));
        if (dim == 1) {
          // pick an arbitrary eigen vector perpendicular to the first vector
          RDGeom::Point3D  second(first.getPerpendicular());
          eigVecs.setVal(1,0, second.x);
          eigVecs.setVal(1,1, second.y);
          eigVecs.setVal(1,2, second.z);
          if (eigVals.getVal(0) > 1.0) {
            eigVals.setVal(1, 1.0);
          } else {
            eigVals.setVal(1, eigVals.getVal(0)/2.0);
          }
        }
        RDGeom::Point3D second(eigVecs.getVal(1,0), eigVecs.getVal(1,1), eigVecs.getVal(1,2));
        // pick the third eigen vector perpendicular to the first two
        RDGeom::Point3D third = first.crossProduct(second);
        eigVecs.setVal(2,0, third.x);
        eigVecs.setVal(2,1, third.y);
        eigVecs.setVal(2,2, third.z);
        if (eigVals.getVal(1) > 1.0) {
          eigVals.setVal(2, 1.0);
        } else {
          eigVals.setVal(2, eigVals.getVal(1)/2.0);
        }
      }
      // now set the transformation
      for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
          trans->setVal(i, j, eigVecs.getVal(i,j));
        }
      }
    }// end of multiple atom system
    trans->TransformPoint(origin);
    trans->SetTranslation(origin);
    delete covMat;

    return trans;
  }
  
  void transformConformer(Conformer &conf, const RDGeom::Transform3D &trans) {
    RDGeom::POINT3D_VECT &positions = conf.getPositions();
    RDGeom::POINT3D_VECT_I pi;
    for (pi = positions.begin(); pi != positions.end(); ++pi) {
      trans.TransformPoint(*pi);  
    }
  }

  void canonicalizeConformer(Conformer &conf, const RDGeom::Point3D *center,
                             bool normalizeCovar, bool ignoreHs) {
    RDGeom::Transform3D *trans = computeCanonicalTransform(conf, center,
                                                           normalizeCovar, ignoreHs);
    transformConformer(conf, *trans);
    delete trans;
  }

  void canonicalizeMol(RDKit::ROMol &mol, bool normalizeCovar, bool ignoreHs) {
    ROMol::ConformerIterator ci;
    for (ci = mol.beginConformers(); ci != mol.endConformers(); ci++) {
      canonicalizeConformer(*(*ci), 0, normalizeCovar, ignoreHs);
    }
  }
}

  
                               
                               
