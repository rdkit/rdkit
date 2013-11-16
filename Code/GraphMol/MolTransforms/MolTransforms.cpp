//  $Id$
// 
//   Copyright (C) 2003-2013 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MolTransforms.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <Numerics/Alignment/AlignPoints.h>
#include <Numerics/SymmMatrix.h>
#include <Numerics/Matrix.h>
#include <Geometry/Transform3D.h>
#include <stack>
#include <boost/dynamic_bitset.hpp>
#include <RDBoost/Exceptions.h>

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
    const ROMol &mol = conf.getOwningMol();
    std::vector<RDGeom::Point3D const *> pts;
    pts.reserve(mol.getNumAtoms());
    unsigned int nAtms = 0;
    for (ROMol::ConstAtomIterator cai = mol.beginAtoms(); cai != mol.endAtoms();
         ++cai) {
      if (((*cai)->getAtomicNum() == 1) && (ignoreHs) ) {
        continue;
      }
      pts.push_back(&(conf.getAtomPos((*cai)->getIdx())));
    }
    return RDNumeric::computeCentroid(pts);
  }

  RDNumeric::DoubleSymmMatrix *computeCovarianceMatrix(const Conformer &conf, 
                                                       const RDGeom::Point3D &center,
                                                       bool normalize, bool ignoreHs) {
    const ROMol &mol = conf.getOwningMol();
    std::vector<RDGeom::Point3D const *> pts;
    pts.reserve(mol.getNumAtoms());
    unsigned int nAtms = 0;
    for (ROMol::ConstAtomIterator cai = mol.beginAtoms(); cai != mol.endAtoms();
         ++cai) {
      if (((*cai)->getAtomicNum() == 1) && (ignoreHs) ) {
        continue;
      }
      pts.push_back(&(conf.getAtomPos((*cai)->getIdx())));
    }
    
    return RDNumeric::computeCovarianceMatrix(pts,center,normalize);
  }

  RDGeom::Transform3D *computeCanonicalTransform(const Conformer &conf,
                                                 const RDGeom::Point3D *center,
                                                 bool normalizeCovar,
                                                 bool ignoreHs) {
    const ROMol &mol = conf.getOwningMol();
    std::vector<RDGeom::Point3D const *> pts;
    pts.reserve(mol.getNumAtoms());
    unsigned int nAtms = 0;
    for (ROMol::ConstAtomIterator cai = mol.beginAtoms(); cai != mol.endAtoms();
         ++cai) {
      if (((*cai)->getAtomicNum() == 1) && (ignoreHs) ) {
        continue;
      }
      pts.push_back(&(conf.getAtomPos((*cai)->getIdx())));
    }
    return RDNumeric::computeCanonicalTransform(pts,center,normalizeCovar);
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

  namespace {
    void _toBeMovedIdxList(const ROMol &mol, unsigned int iAtomId, unsigned int jAtomId,
                           std::list<unsigned int> &alist) {
      const Atom *iAtom = mol.getAtomWithIdx(iAtomId);
      const Atom *jAtom = mol.getAtomWithIdx(jAtomId);
      unsigned int nAtoms = mol.getNumAtoms();
      boost::dynamic_bitset<> visitedIdx(nAtoms);
      std::stack<unsigned int> stack;
      stack.push(jAtomId);
      visitedIdx[iAtomId] = 1;
      visitedIdx[jAtomId] = 1;
      unsigned int tIdx;
      unsigned int wIdx;
      ROMol::ADJ_ITER nbrIdx;
      ROMol::ADJ_ITER endNbrs;
      bool doMainLoop;
      while (stack.size()) {
        doMainLoop = false;
        tIdx = stack.top();
        const Atom *tAtom = mol.getAtomWithIdx(tIdx);
        boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(tAtom);
        unsigned int eIdx;
        for (eIdx = 0; nbrIdx != endNbrs; ++nbrIdx, ++eIdx) {
          wIdx = (mol[*nbrIdx].get())->getIdx();
          if (!visitedIdx[wIdx]) {
            visitedIdx[wIdx] = 1;
            stack.push(wIdx);
            doMainLoop = true;
            break;
          }
        }
        if (doMainLoop) {
          continue;
        }
        visitedIdx[tIdx] = 1;
        stack.pop();
      }
      alist.clear();
      for (unsigned int i = 0; i < nAtoms; ++i) {
        if (visitedIdx[i] && (i != iAtomId)) {
          alist.push_back(i);
        }
      }
    }
  }
  double getBondLength(Conformer &conf,
    unsigned int iAtomId, unsigned int jAtomId) {
    RDGeom::POINT3D_VECT &pos = conf.getPositions();
    RANGE_CHECK(0, iAtomId, pos.size() - 1);
    RANGE_CHECK(0, jAtomId, pos.size() - 1);
    
    return (pos[iAtomId] - pos[jAtomId]).length();
  }

  void setBondLength(Conformer &conf,
    unsigned int iAtomId, unsigned int jAtomId, double value) {
    RDGeom::POINT3D_VECT &pos = conf.getPositions();
    RANGE_CHECK(0, iAtomId, pos.size() - 1);
    RANGE_CHECK(0, jAtomId, pos.size() - 1);
    ROMol &mol = conf.getOwningMol();
    Bond *bond = mol.getBondBetweenAtoms(iAtomId, jAtomId);
    if(!bond) throw ValueErrorException("atoms i and j must be bonded");
    if(queryIsBondInRing(bond)) throw ValueErrorException("bond (i,j) must not belong to a ring");
    RDGeom::Point3D v = pos[iAtomId] - pos[jAtomId];
    double origValue = v.length();
    if(origValue <= 1.e-8) throw ValueErrorException("atoms i and j have identical 3D coordinates");
    
    // get all atoms bonded to j
    std::list<unsigned int> alist;
    _toBeMovedIdxList(mol, iAtomId, jAtomId, alist);
    v *= (value / origValue - 1.);
    for (std::list<unsigned int>::iterator it = alist.begin(); it != alist.end(); ++it) {
      pos[*it] -= v;
    }
  }
  
  double getAngleRad(Conformer &conf,
    unsigned int iAtomId, unsigned int jAtomId, unsigned int kAtomId) {
    RDGeom::POINT3D_VECT &pos = conf.getPositions();
    RANGE_CHECK(0, iAtomId, pos.size() - 1);
    RANGE_CHECK(0, jAtomId, pos.size() - 1);
    RANGE_CHECK(0, kAtomId, pos.size() - 1);
    RDGeom::Point3D rJI = pos[iAtomId] - pos[jAtomId];
    double rJISqLength = rJI.lengthSq();
    if(rJISqLength <= 1.e-16) throw ValueErrorException("atoms i and j have identical 3D coordinates");
    RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
    double rJKSqLength = rJK.lengthSq();
    if(rJKSqLength <= 1.e-16) throw ValueErrorException("atoms j and k have identical 3D coordinates");
    return rJI.angleTo(rJK);
  }

  void setAngleRad(Conformer &conf, unsigned int iAtomId,
    unsigned int jAtomId, unsigned int kAtomId, double value) {
    RDGeom::POINT3D_VECT &pos = conf.getPositions();
    RANGE_CHECK(0, iAtomId, pos.size() - 1);
    RANGE_CHECK(0, jAtomId, pos.size() - 1);
    RANGE_CHECK(0, kAtomId, pos.size() - 1);
    ROMol &mol = conf.getOwningMol();
    Bond *bondJI = mol.getBondBetweenAtoms(jAtomId, iAtomId);
    if(!bondJI) throw ValueErrorException("atoms i and j must be bonded");
    Bond *bondJK = mol.getBondBetweenAtoms(jAtomId, kAtomId);
    if(!bondJK) throw ValueErrorException("atoms j and k must be bonded");
    if(queryIsBondInRing(bondJI) && queryIsBondInRing(bondJK))
       throw ValueErrorException("bonds (i,j) and (j,k) must not both belong to a ring");

    RDGeom::Point3D rJI = pos[iAtomId] - pos[jAtomId];
    double rJISqLength = rJI.lengthSq();
    if(rJISqLength <= 1.e-16) throw ValueErrorException("atoms i and j have identical 3D coordinates");
    RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
    double rJKSqLength = rJK.lengthSq();
    if(rJKSqLength <= 1.e-16) throw ValueErrorException("atoms j and k have identical 3D coordinates");
    
    // we only need to rotate by delta with respect to the current angle value
    value -= rJI.angleTo(rJK);
    RDGeom::Point3D &rotAxisBegin = pos[jAtomId];
    // our rotation axis is the normal to the plane of atoms i, j, k
    RDGeom::Point3D rotAxisEnd = rJI.crossProduct(rJK) + pos[jAtomId];
    RDGeom::Point3D rotAxis = rotAxisEnd - rotAxisBegin;
    rotAxis.normalize();
    // get all atoms bonded to j and loop through them
    std::list<unsigned int> alist;
    _toBeMovedIdxList(mol, jAtomId, kAtomId,alist);
    for (std::list<unsigned int>::iterator it = alist.begin(); it != alist.end(); ++it) {
      // translate atom so that it coincides with the origin of rotation
      pos[*it] -= rotAxisBegin;
      // rotate around our rotation axis
      RDGeom::Transform3D rotByAngle;
      rotByAngle.SetRotation(value, rotAxis);
      rotByAngle.TransformPoint(pos[*it]);
      // translate atom back
      pos[*it] += rotAxisBegin;
    }
  }
  
  double getDihedralRad(Conformer &conf, unsigned int iAtomId,
    unsigned int jAtomId, unsigned int kAtomId, unsigned int lAtomId) {
    RDGeom::POINT3D_VECT &pos = conf.getPositions();
    RANGE_CHECK(0, iAtomId, pos.size() - 1);
    RANGE_CHECK(0, jAtomId, pos.size() - 1);
    RANGE_CHECK(0, kAtomId, pos.size() - 1);
    RANGE_CHECK(0, lAtomId, pos.size() - 1);
    RDGeom::Point3D rIJ = pos[jAtomId] - pos[iAtomId];
    double rIJSqLength = rIJ.lengthSq();
    if(rIJSqLength <= 1.e-16) throw ValueErrorException("atoms i and j have identical 3D coordinates");
    RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
    double rJKSqLength = rJK.lengthSq();
    if(rJKSqLength <= 1.e-16) throw ValueErrorException("atoms j and k have identical 3D coordinates");
    RDGeom::Point3D rKL = pos[lAtomId] - pos[kAtomId];
    double rKLSqLength = rKL.lengthSq();
    if(rKLSqLength <= 1.e-16) throw ValueErrorException("atoms k and l have identical 3D coordinates");
    
    RDGeom::Point3D nIJK = rIJ.crossProduct(rJK);
    double nIJKSqLength = nIJK.lengthSq();
    RDGeom::Point3D nJKL = rJK.crossProduct(rKL);
    double nJKLSqLength = nJKL.lengthSq();
    RDGeom::Point3D m = nIJK.crossProduct(rJK);
    // we want a signed dihedral, that's why we use atan2 instead of acos
    return -atan2(m.dotProduct(nJKL) / sqrt(nJKLSqLength * m.lengthSq()),
      nIJK.dotProduct(nJKL) / sqrt(nIJKSqLength * nJKLSqLength));
  }

  void setDihedralRad(Conformer &conf, unsigned int iAtomId,
    unsigned int jAtomId, unsigned int kAtomId, unsigned int lAtomId,
    double value) {
    RDGeom::POINT3D_VECT &pos = conf.getPositions();
    RANGE_CHECK(0, iAtomId, pos.size() - 1);
    RANGE_CHECK(0, jAtomId, pos.size() - 1);
    RANGE_CHECK(0, kAtomId, pos.size() - 1);
    RANGE_CHECK(0, lAtomId, pos.size() - 1);
    ROMol &mol = conf.getOwningMol();
    Bond *bondIJ = mol.getBondBetweenAtoms(iAtomId, jAtomId);
    if(!bondIJ) throw ValueErrorException("atoms i and j must be bonded");
    Bond *bondJK = mol.getBondBetweenAtoms(jAtomId, kAtomId);
    if(!bondJK) throw ValueErrorException("atoms j and k must be bonded");
    Bond *bondKL = mol.getBondBetweenAtoms(kAtomId, lAtomId);
    if(!bondKL) throw ValueErrorException("atoms k and l must be bonded");

    if(queryIsBondInRing(bondJK)) throw ValueErrorException("bond (j,k) must not belong to a ring");
    RDGeom::Point3D rIJ = pos[jAtomId] - pos[iAtomId];
    double rIJSqLength = rIJ.lengthSq();
    if(rIJSqLength <= 1.e-16) throw ValueErrorException("atoms i and j have identical 3D coordinates");
    RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
    double rJKSqLength = rJK.lengthSq();
    if(rJKSqLength <= 1.e-16) throw ValueErrorException("atoms j and k have identical 3D coordinates");
    RDGeom::Point3D rKL = pos[lAtomId] - pos[kAtomId];
    double rKLSqLength = rKL.lengthSq();
    if(rKLSqLength <= 1.e-16) throw ValueErrorException("atoms k and l have identical 3D coordinates");

    RDGeom::Point3D nIJK = rIJ.crossProduct(rJK);
    double nIJKSqLength = nIJK.lengthSq();
    RDGeom::Point3D nJKL = rJK.crossProduct(rKL);
    double nJKLSqLength = nJKL.lengthSq();
    RDGeom::Point3D m = nIJK.crossProduct(rJK);
    // we only need to rotate by delta with respect to the current dihedral value
    value -= -atan2(m.dotProduct(nJKL) / sqrt(nJKLSqLength * m.lengthSq()),
      nIJK.dotProduct(nJKL) / sqrt(nIJKSqLength * nJKLSqLength));
    // our rotation axis is the (j,k) bond
    RDGeom::Point3D &rotAxisBegin = pos[jAtomId];
    RDGeom::Point3D &rotAxisEnd = pos[kAtomId];
    RDGeom::Point3D rotAxis = rotAxisEnd - rotAxisBegin;
    rotAxis.normalize();
    // get all atoms bonded to k and loop through them
    std::list<unsigned int> alist;
    _toBeMovedIdxList(mol, jAtomId, kAtomId,alist);
    for (std::list<unsigned int>::iterator it = alist.begin(); it != alist.end(); ++it) {
      // translate atom so that it coincides with the origin of rotation
      pos[*it] -= rotAxisBegin;
      // rotate around our rotation axis
      RDGeom::Transform3D rotByAngle;
      rotByAngle.SetRotation(value, rotAxis);
      rotByAngle.TransformPoint(pos[*it]);
      // translate atom back
      pos[*it] += rotAxisBegin;
    }
  }
}
