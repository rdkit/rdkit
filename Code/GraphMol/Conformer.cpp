// $Id$
//
//  Copyright (C) 2001-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Conformer.h"
#include "ROMol.h"
#include "QueryOps.h"
#include <stack>
#include <boost/dynamic_bitset.hpp>

namespace RDKit {

  Conformer::Conformer(const Conformer &conf) {
    dp_mol = 0;
    int i, nat = conf.getNumAtoms();
    d_positions.reserve(nat);
    
    for (i = 0; i < nat; i++) {
      d_positions.push_back(conf.getAtomPos(i));
    }
    d_id = conf.getId();
    df_is3D = conf.is3D();
  }

  void Conformer::setOwningMol(ROMol *mol) {
    PRECONDITION(mol, "");
    dp_mol = mol;
  }

  void Conformer::setOwningMol(ROMol &mol) {setOwningMol(&mol);}

  const RDGeom::POINT3D_VECT &Conformer::getPositions() const {
    if (dp_mol) {
      PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
    }
    return d_positions;
  }

  RDGeom::POINT3D_VECT &Conformer::getPositions() {
    if (dp_mol) {
      PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
    }
    return d_positions;
  }

  const RDGeom::Point3D &Conformer::getAtomPos(unsigned int atomId) const {
    if (dp_mol) {
      PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
    }
    RANGE_CHECK(0,atomId,d_positions.size()-1);
    return d_positions[atomId];
  } 
  
  RDGeom::Point3D &Conformer::getAtomPos(unsigned int atomId) {
    if (dp_mol) {
      PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
    }
    RANGE_CHECK(0,atomId,d_positions.size()-1);
    return d_positions[atomId];
  }


  std::list<unsigned int> _toBeMovedIdxList(ROMol &mol, unsigned int iAtomId, unsigned int jAtomId) {
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
    std::list<unsigned int> list;
    for (unsigned int i = 0; i < nAtoms; ++i) {
      if (visitedIdx[i] && (i != iAtomId)) {
        list.push_back(i);
      }
    }
    
    return list;
  }


  void _rotateAroundAxis(RDGeom::Point3D &p,
    RDGeom::Point3D &rotAxis, double dist, double value)
  {
    // if rotAxis does not coincide with zAxis we need to
    // perform two rotations around xAxis and yAxis to
    // make it coincide
    RDGeom::Point3D t = p;
    if (dist > 1.e-8) {
      t.y = (p.y * rotAxis.z - p.z * rotAxis.y) / dist;
      t.z = (p.y * rotAxis.y + p.z * rotAxis.z) / dist;
    }
    p.x = t.x * dist - t.z * rotAxis.x;
    p.y = t.y;
    p.z = t.x * rotAxis.x + t.z * dist;
    // now we perform the desired rotation about zAxis
    t.x = p.x * cos(value) - p.y * sin(value);
    t.y = p.x * sin(value) + p.y * cos(value);
    t.z = p.z;
    // and now we rotate things back as they where
    p.x = t.x * dist + t.z * rotAxis.x;
    p.y = t.y;
    p.z = -t.x * rotAxis.x + t.z * dist;
    if (dist > 1.e-8) {
      t.y = (p.y * rotAxis.z + p.z * rotAxis.y) / dist;
      t.z = (-p.y * rotAxis.y + p.z * rotAxis.z) / dist;
      p.y = t.y;
      p.z = t.z;
    }
  }


  const double Conformer::getBondLength(unsigned int iAtomId, unsigned int jAtomId) const {
    RANGE_CHECK(0, iAtomId, d_positions.size() - 1);
    RANGE_CHECK(0, jAtomId, d_positions.size() - 1);
    
    return (d_positions[iAtomId] - d_positions[jAtomId]).length();
  }


  void Conformer::setBondLength(unsigned int iAtomId,
    unsigned int jAtomId, double value) {
    RANGE_CHECK(0, iAtomId, d_positions.size() - 1);
    RANGE_CHECK(0, jAtomId, d_positions.size() - 1);
    Bond *bond = dp_mol->getBondBetweenAtoms(iAtomId, jAtomId);
    PRECONDITION(bond, "atoms i and j must be bonded");
    PRECONDITION(!queryIsBondInRing(bond), "bond (i,j) must not belong to a ring");
    RDGeom::Point3D v = d_positions[iAtomId] - d_positions[jAtomId];
    double origValue = v.length();
    PRECONDITION(origValue > 1.e-8, "atoms i and j have identical 3D coordinates");
    
    // get all atoms bonded to j
    std::list<unsigned int> list = _toBeMovedIdxList(*dp_mol, iAtomId, jAtomId);
    v *= (value / origValue - 1.);
    for (std::list<unsigned int>::iterator it = list.begin(); it != list.end(); ++it) {
      d_positions[*it] -= v;
    }
  }

  
  const double Conformer::getAngleRad(unsigned int iAtomId,
    unsigned int jAtomId, unsigned int kAtomId) const {
    RANGE_CHECK(0, iAtomId, d_positions.size() - 1);
    RANGE_CHECK(0, jAtomId, d_positions.size() - 1);
    RANGE_CHECK(0, kAtomId, d_positions.size() - 1);
    RDGeom::Point3D rJI = d_positions[iAtomId] - d_positions[jAtomId];
    double rJISqLength = rJI.lengthSq();
    PRECONDITION(rJISqLength > 1.e-16, "atoms i and j have identical 3D coordinates");
    RDGeom::Point3D rJK = d_positions[kAtomId] - d_positions[jAtomId];
    double rJKSqLength = rJK.lengthSq();
    PRECONDITION(rJKSqLength > 1.e-16, "atoms j and k have identical 3D coordinates");
    
    return acos(rJI.dotProduct(rJK) / sqrt(rJISqLength * rJKSqLength));
  }


  void Conformer::setAngleRad(unsigned int iAtomId,
    unsigned int jAtomId, unsigned int kAtomId, double value) {
    RANGE_CHECK(0, iAtomId, d_positions.size() - 1);
    RANGE_CHECK(0, jAtomId, d_positions.size() - 1);
    RANGE_CHECK(0, kAtomId, d_positions.size() - 1);
    Bond *bondJI = dp_mol->getBondBetweenAtoms(jAtomId, iAtomId);
    PRECONDITION(bondJI, "atoms i and j must be bonded");
    Bond *bondJK = dp_mol->getBondBetweenAtoms(jAtomId, kAtomId);
    PRECONDITION(bondJK, "atoms j and k must be bonded");
    PRECONDITION(!(queryIsBondInRing(bondJI) && queryIsBondInRing(bondJK)),
      "bonds (i,j) and (j,k) must not both belong to a ring");
    RDGeom::Point3D rJI = d_positions[iAtomId] - d_positions[jAtomId];
    double rJISqLength = rJI.lengthSq();
    PRECONDITION(rJISqLength > 1.e-16, "atoms i and j have identical 3D coordinates");
    RDGeom::Point3D rJK = d_positions[kAtomId] - d_positions[jAtomId];
    double rJKSqLength = rJK.lengthSq();
    PRECONDITION(rJKSqLength > 1.e-16, "atoms j and k have identical 3D coordinates");
    
    // we only need to rotate by delta with respect to the current angle value
    value -= acos(rJI.dotProduct(rJK) / sqrt(rJISqLength * rJKSqLength));
    RDGeom::Point3D &rotAxisBegin = d_positions[jAtomId];
    // our rotation axis is the normal to the plane of atoms i, j, k
    RDGeom::Point3D rotAxisEnd = rJI.crossProduct(rJK) + d_positions[jAtomId];
    RDGeom::Point3D rotAxis = rotAxisEnd - rotAxisBegin;
    rotAxis.normalize();
    double dist = sqrt(rotAxis.y * rotAxis.y + rotAxis.z * rotAxis.z);
    // get all atoms bonded to j and loop through them
    std::list<unsigned int> list = _toBeMovedIdxList(*dp_mol, jAtomId, kAtomId);
    for (std::list<unsigned int>::iterator it = list.begin(); it != list.end(); ++it) {
      // translate atom so that it coincides with the origin of rotation
      RDGeom::Point3D p = d_positions[*it] - rotAxisBegin;
      _rotateAroundAxis(p, rotAxis, dist, value);
      // translate atom back
      d_positions[*it] = p + rotAxisBegin;
    }
  }

  
  const double Conformer::getDihedralRad(unsigned int iAtomId,
    unsigned int jAtomId, unsigned int kAtomId, unsigned int lAtomId) const {
    RANGE_CHECK(0, iAtomId, d_positions.size() - 1);
    RANGE_CHECK(0, jAtomId, d_positions.size() - 1);
    RANGE_CHECK(0, kAtomId, d_positions.size() - 1);
    RANGE_CHECK(0, lAtomId, d_positions.size() - 1);
    RDGeom::Point3D rIJ = d_positions[jAtomId] - d_positions[iAtomId];
    double rIJSqLength = rIJ.lengthSq();
    PRECONDITION(rIJSqLength > 1.e-16, "atoms i and j have identical 3D coordinates");
    RDGeom::Point3D rJK = d_positions[kAtomId] - d_positions[jAtomId];
    double rJKSqLength = rJK.lengthSq();
    PRECONDITION(rJKSqLength > 1.e-16, "atoms j and k have identical 3D coordinates");
    RDGeom::Point3D rKL = d_positions[lAtomId] - d_positions[kAtomId];
    double rKLSqLength = rKL.lengthSq();
    PRECONDITION(rKLSqLength > 1.e-16, "atoms k and l have identical 3D coordinates");
    
    RDGeom::Point3D nIJK = rIJ.crossProduct(rJK);
    double nIJKSqLength = nIJK.lengthSq();
    RDGeom::Point3D nJKL = rJK.crossProduct(rKL);
    double nJKLSqLength = nJKL.lengthSq();
    RDGeom::Point3D m = nIJK.crossProduct(rJK);
    // we want a signed dihedral, that's why we use atan2 instead of acos
    return -atan2(m.dotProduct(nJKL) / sqrt(nJKLSqLength * m.lengthSq()),
      nIJK.dotProduct(nJKL) / sqrt(nIJKSqLength * nJKLSqLength));
  }


  void Conformer::setDihedralRad(unsigned int iAtomId,
    unsigned int jAtomId, unsigned int kAtomId, unsigned int lAtomId,
    double value) {
    RANGE_CHECK(0, iAtomId, d_positions.size() - 1);
    RANGE_CHECK(0, jAtomId, d_positions.size() - 1);
    RANGE_CHECK(0, kAtomId, d_positions.size() - 1);
    RANGE_CHECK(0, lAtomId, d_positions.size() - 1);
    Bond *bondIJ = dp_mol->getBondBetweenAtoms(iAtomId, jAtomId);
    PRECONDITION(bondIJ, "atoms i and j must be bonded");
    Bond *bondJK = dp_mol->getBondBetweenAtoms(jAtomId, kAtomId);
    PRECONDITION(bondJK, "atoms j and k must be bonded");
    Bond *bondKL = dp_mol->getBondBetweenAtoms(kAtomId, lAtomId);
    PRECONDITION(bondJK, "atoms k and l must be bonded");
    PRECONDITION(!queryIsBondInRing(bondJK), "bond (j,k) must not belong to a ring");
    RDGeom::Point3D rIJ = d_positions[jAtomId] - d_positions[iAtomId];
    double rIJSqLength = rIJ.lengthSq();
    PRECONDITION(rIJSqLength > 1.e-16, "atoms i and j have identical 3D coordinates");
    RDGeom::Point3D rJK = d_positions[kAtomId] - d_positions[jAtomId];
    double rJKSqLength = rJK.lengthSq();
    PRECONDITION(rJKSqLength > 1.e-16, "atoms j and k have identical 3D coordinates");
    RDGeom::Point3D rKL = d_positions[lAtomId] - d_positions[kAtomId];
    double rKLSqLength = rKL.lengthSq();
    PRECONDITION(rKLSqLength > 1.e-16, "atoms k and l have identical 3D coordinates");
    
    RDGeom::Point3D nIJK = rIJ.crossProduct(rJK);
    double nIJKSqLength = nIJK.lengthSq();
    RDGeom::Point3D nJKL = rJK.crossProduct(rKL);
    double nJKLSqLength = nJKL.lengthSq();
    RDGeom::Point3D m = nIJK.crossProduct(rJK);
    // we only need to rotate by delta with respect to the current dihedral value
    value -= -atan2(m.dotProduct(nJKL) / sqrt(nJKLSqLength * m.lengthSq()),
      nIJK.dotProduct(nJKL) / sqrt(nIJKSqLength * nJKLSqLength));
    // our rotation axis is the (j,k) bond
    RDGeom::Point3D &rotAxisBegin = d_positions[jAtomId];
    RDGeom::Point3D &rotAxisEnd = d_positions[kAtomId];
    RDGeom::Point3D rotAxis = rotAxisEnd - rotAxisBegin;
    rotAxis.normalize();
    double dist = sqrt(rotAxis.y * rotAxis.y + rotAxis.z * rotAxis.z);
    // get all atoms bonded to k and loop through them
    std::list<unsigned int> list = _toBeMovedIdxList(*dp_mol, jAtomId, kAtomId);
    for (std::list<unsigned int>::iterator it = list.begin(); it != list.end(); ++it) {
      // translate atom so that it coincides with the origin of rotation
      RDGeom::Point3D p = d_positions[*it] - rotAxisBegin;
      _rotateAroundAxis(p, rotAxis, dist, value);
      // translate atom back
      d_positions[*it] = p + rotAxisBegin;
    }
  }
}
