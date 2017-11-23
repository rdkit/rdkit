//
//   Copyright (C) 2003-2016 Greg Landrum and Rational Discovery LLC
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
#include <Numerics/EigenSolvers/PowerEigenSolver.h>
#include <Numerics/SymmMatrix.h>
#include <Numerics/Matrix.h>
#include <Geometry/Transform3D.h>
#include <stack>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/Exceptions.h>

#define EIGEN_TOLERANCE 1.0e-2
namespace MolTransforms {

using namespace RDKit;
void transformAtom(Atom *atom, RDGeom::Transform3D &tform) {
  PRECONDITION(atom, "no atom");
  ROMol &mol = atom->getOwningMol();
  for (ROMol::ConstConformerIterator ci = mol.beginConformers();
       ci != mol.endConformers(); ci++) {
    RDGeom::Point3D &pos = (*ci)->getAtomPos(atom->getIdx());
    tform.TransformPoint(pos);
  }
  // atom->setPos(pos);
}
void transformMolsAtoms(ROMol *mol, RDGeom::Transform3D &tform) {
  PRECONDITION(mol, "no molecule");

  ROMol::AtomIterator atomIt;
  for (atomIt = mol->beginAtoms(); atomIt != mol->endAtoms(); atomIt++) {
    transformAtom(*atomIt, tform);
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
namespace {
void computeCovarianceTerms(const Conformer &conf,
                            const RDGeom::Point3D &center, double &xx,
                            double &xy, double &xz, double &yy, double &yz,
                            double &zz, bool normalize, bool ignoreHs,
                            const std::vector<double> *weights) {
  PRECONDITION(!weights || weights->size() >= conf.getNumAtoms(),
               "bad weights vector");

  xx = xy = xz = yy = yz = zz = 0.0;
  const ROMol &mol = conf.getOwningMol();
  double wSum = 0.0;
  for (ROMol::ConstAtomIterator cai = mol.beginAtoms(); cai != mol.endAtoms();
       cai++) {
    if (((*cai)->getAtomicNum() == 1) && (ignoreHs)) {
      continue;
    }
    RDGeom::Point3D loc = conf.getAtomPos((*cai)->getIdx());
    loc -= center;
    double w = 1.0;
    if (weights) {
      w = (*weights)[(*cai)->getIdx()];
    }
    wSum += w;
    xx += w * loc.x * loc.x;
    xy += w * loc.x * loc.y;
    xz += w * loc.x * loc.z;
    yy += w * loc.y * loc.y;
    yz += w * loc.y * loc.z;
    zz += w * loc.z * loc.z;
  }
  if (normalize) {
    xx /= wSum;
    xy /= wSum;
    xz /= wSum;
    yy /= wSum;
    yz /= wSum;
    zz /= wSum;
  }
}

RDNumeric::DoubleSymmMatrix *computeCovarianceMatrix(
    const Conformer &conf, const RDGeom::Point3D &center, bool normalize,
    bool ignoreHs) {
  double xx, xy, xz, yy, yz, zz;
  computeCovarianceTerms(conf, center, xx, xy, xz, yy, yz, zz, normalize,
                         ignoreHs, NULL);
  RDNumeric::DoubleSymmMatrix *res = new RDNumeric::DoubleSymmMatrix(3, 3);
  res->setVal(0, 0, xx);
  res->setVal(0, 1, xy);
  res->setVal(0, 2, xz);
  res->setVal(1, 1, yy);
  res->setVal(1, 2, yz);
  res->setVal(2, 2, zz);
  return res;
}

void computeInertiaTerms(const Conformer &conf, const RDGeom::Point3D &center,
                         double &xx, double &xy, double &xz, double &yy,
                         double &yz, double &zz, bool ignoreHs,
                         const std::vector<double> *weights) {
  PRECONDITION(!weights || weights->size() >= conf.getNumAtoms(),
               "bad weights vector");

  xx = xy = xz = yy = yz = zz = 0.0;
  const ROMol &mol = conf.getOwningMol();
  for (ROMol::ConstAtomIterator cai = mol.beginAtoms(); cai != mol.endAtoms();
       cai++) {
    if (((*cai)->getAtomicNum() == 1) && (ignoreHs)) {
      continue;
    }
    RDGeom::Point3D loc = conf.getAtomPos((*cai)->getIdx());
    loc -= center;
    double w = 1.0;
    if (weights) {
      w = (*weights)[(*cai)->getIdx()];
    }
    xx += w * (loc.y * loc.y + loc.z * loc.z);
    yy += w * (loc.x * loc.x + loc.z * loc.z);
    zz += w * (loc.y * loc.y + loc.x * loc.x);
    xy -= w * loc.x * loc.y;
    xz -= w * loc.x * loc.z;
    yz -= w * loc.z * loc.y;
  }
}
}
#ifdef RDK_HAS_EIGEN3
#include <Eigen/Dense>

bool computePrincipalAxesAndMoments(const RDKit::Conformer &conf,
                                    Eigen::Matrix3d &axes,
                                    Eigen::Vector3d &moments, bool ignoreHs,
                                    bool force,
                                    const std::vector<double> *weights) {
  PRECONDITION((!weights || weights->size() >= conf.getNumAtoms()),
               "bad weights vector");
  const char *axesPropName = ignoreHs ? "_principalAxes_noH" : "_principalAxes";
  const char *momentsPropName =
      ignoreHs ? "_principalMoments_noH" : "_principalMoments";
  if (!weights && !force && conf.getOwningMol().hasProp(axesPropName) &&
      conf.getOwningMol().hasProp(momentsPropName)) {
    conf.getOwningMol().getProp(axesPropName, axes);
    conf.getOwningMol().getProp(momentsPropName, moments);
    return true;
  }
  const ROMol &mol = conf.getOwningMol();
  RDGeom::Point3D origin(0, 0, 0);
  double wSum = 0.0;
  for (unsigned int i = 0; i < conf.getNumAtoms(); ++i) {
    if (ignoreHs && mol.getAtomWithIdx(i)->getAtomicNum() == 1) continue;
    double w = 1.0;
    if (weights) {
      w = (*weights)[i];
    }
    wSum += w;
    origin += conf.getAtomPos(i) * w;
  }
  // std::cerr<<"  origin: "<<origin<<" "<<wSum<<std::endl;
  origin /= wSum;

  double sumXX, sumXY, sumXZ, sumYY, sumYZ, sumZZ;
  computeInertiaTerms(conf, origin, sumXX, sumXY, sumXZ, sumYY, sumYZ, sumZZ,
                      ignoreHs, weights);

  Eigen::Matrix3d mat;
  mat << sumXX, sumXY, sumXZ, sumXY, sumYY, sumYZ, sumXZ, sumYZ, sumZZ;
  // std::cerr<<"  matrix: "<<mat<<std::endl;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(mat);
  if (eigensolver.info() != Eigen::Success) {
    BOOST_LOG(rdErrorLog) << "eigenvalue calculation did not converge"
                          << std::endl;
    return false;
  }

  axes = eigensolver.eigenvectors();
  moments = eigensolver.eigenvalues();
  if (!weights) {
    conf.getOwningMol().setProp(axesPropName, axes, true);
    conf.getOwningMol().setProp(momentsPropName, moments, true);
  }
  return true;
}
bool computePrincipalAxesAndMomentsFromGyrationMatrix(
    const RDKit::Conformer &conf, Eigen::Matrix3d &axes,
    Eigen::Vector3d &moments, bool ignoreHs, bool force,
    const std::vector<double> *weights) {
  PRECONDITION((!weights || weights->size() >= conf.getNumAtoms()),
               "bad weights vector");
  const char *axesPropName =
      ignoreHs ? "_principalAxes_noH_cov" : "_principalAxes_cov";
  const char *momentsPropName =
      ignoreHs ? "_principalMoments_noH_cov" : "_principalMoments_cov";
  if (!weights && !force && conf.getOwningMol().hasProp(axesPropName) &&
      conf.getOwningMol().hasProp(momentsPropName)) {
    conf.getOwningMol().getProp(axesPropName, axes);
    conf.getOwningMol().getProp(momentsPropName, moments);
    return true;
  }
  const ROMol &mol = conf.getOwningMol();
  RDGeom::Point3D origin(0, 0, 0);
  double wSum = 0.0;
  for (unsigned int i = 0; i < conf.getNumAtoms(); ++i) {
    if (ignoreHs && mol.getAtomWithIdx(i)->getAtomicNum() == 1) continue;
    double w = 1.0;
    if (weights) {
      w = (*weights)[i];
    }
    wSum += w;
    origin += conf.getAtomPos(i) * w;
  }
  // std::cerr<<"  origin: "<<origin<<" "<<wSum<<std::endl;
  origin /= wSum;

  double sumXX, sumXY, sumXZ, sumYY, sumYZ, sumZZ;
  computeCovarianceTerms(conf, origin, sumXX, sumXY, sumXZ, sumYY, sumYZ, sumZZ,
                         true, ignoreHs, weights);

  Eigen::Matrix3d mat;
  mat << sumXX, sumXY, sumXZ, sumXY, sumYY, sumYZ, sumXZ, sumYZ, sumZZ;
  // std::cerr<<"  matrix: "<<mat<<std::endl;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(mat);
  if (eigensolver.info() != Eigen::Success) {
    BOOST_LOG(rdErrorLog) << "eigenvalue calculation did not converge"
                          << std::endl;
    return false;
  }

  axes = eigensolver.eigenvectors();
  moments = eigensolver.eigenvalues();
  if (!weights) {
    conf.getOwningMol().setProp(axesPropName, axes, true);
    conf.getOwningMol().setProp(momentsPropName, moments, true);
  }
  return true;
}
#endif

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
  RDNumeric::DoubleSymmMatrix *covMat =
      computeCovarianceMatrix(conf, origin, normalizeCovar, ignoreHs);
  // find the eigen values and eigen vectors for the covMat
  RDNumeric::DoubleMatrix eigVecs(3, 3);
  RDNumeric::DoubleVector eigVals(3);
  // if we have a single atom system we don't need to do anyhting other than
  // setting translation
  // translation
  unsigned int nAtms = conf.getNumAtoms();
  RDGeom::Transform3D *trans = new RDGeom::Transform3D;

  // set the translation
  origin *= -1.0;
  // trans->SetTranslation(origin);
  // if we have a single atom system we don't need to do anyhting setting
  // translation is sufficient
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
      RDGeom::Point3D first(eigVecs.getVal(0, 0), eigVecs.getVal(0, 1),
                            eigVecs.getVal(0, 2));
      if (dim == 1) {
        // pick an arbitrary eigen vector perpendicular to the first vector
        RDGeom::Point3D second(first.getPerpendicular());
        eigVecs.setVal(1, 0, second.x);
        eigVecs.setVal(1, 1, second.y);
        eigVecs.setVal(1, 2, second.z);
        if (eigVals.getVal(0) > 1.0) {
          eigVals.setVal(1, 1.0);
        } else {
          eigVals.setVal(1, eigVals.getVal(0) / 2.0);
        }
      }
      RDGeom::Point3D second(eigVecs.getVal(1, 0), eigVecs.getVal(1, 1),
                             eigVecs.getVal(1, 2));
      // pick the third eigen vector perpendicular to the first two
      RDGeom::Point3D third = first.crossProduct(second);
      eigVecs.setVal(2, 0, third.x);
      eigVecs.setVal(2, 1, third.y);
      eigVecs.setVal(2, 2, third.z);
      if (eigVals.getVal(1) > 1.0) {
        eigVals.setVal(2, 1.0);
      } else {
        eigVals.setVal(2, eigVals.getVal(1) / 2.0);
      }
    }
    // now set the transformation
    for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) {
        trans->setVal(i, j, eigVecs.getVal(i, j));
      }
    }
  }  // end of multiple atom system
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
  RDGeom::Transform3D *trans =
      computeCanonicalTransform(conf, center, normalizeCovar, ignoreHs);
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
void _toBeMovedIdxList(const ROMol &mol, unsigned int iAtomId,
                       unsigned int jAtomId, std::list<unsigned int> &alist) {
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
double getBondLength(const Conformer &conf, unsigned int iAtomId,
                     unsigned int jAtomId) {
  const RDGeom::POINT3D_VECT &pos = conf.getPositions();
  URANGE_CHECK(iAtomId, pos.size());
  URANGE_CHECK(jAtomId, pos.size());

  return (pos[iAtomId] - pos[jAtomId]).length();
}

void setBondLength(Conformer &conf, unsigned int iAtomId, unsigned int jAtomId,
                   double value) {
  RDGeom::POINT3D_VECT &pos = conf.getPositions();
  URANGE_CHECK(iAtomId, pos.size());
  URANGE_CHECK(jAtomId, pos.size());
  ROMol &mol = conf.getOwningMol();
  Bond *bond = mol.getBondBetweenAtoms(iAtomId, jAtomId);
  if (!bond) throw ValueErrorException("atoms i and j must be bonded");
  if (queryIsBondInRing(bond))
    throw ValueErrorException("bond (i,j) must not belong to a ring");
  RDGeom::Point3D v = pos[iAtomId] - pos[jAtomId];
  double origValue = v.length();
  if (origValue <= 1.e-8)
    throw ValueErrorException("atoms i and j have identical 3D coordinates");

  // get all atoms bonded to j
  std::list<unsigned int> alist;
  _toBeMovedIdxList(mol, iAtomId, jAtomId, alist);
  v *= (value / origValue - 1.);
  for (std::list<unsigned int>::iterator it = alist.begin(); it != alist.end();
       ++it) {
    pos[*it] -= v;
  }
}

double getAngleRad(const Conformer &conf, unsigned int iAtomId,
                   unsigned int jAtomId, unsigned int kAtomId) {
  const RDGeom::POINT3D_VECT &pos = conf.getPositions();
  URANGE_CHECK(iAtomId, pos.size());
  URANGE_CHECK(jAtomId, pos.size());
  URANGE_CHECK(kAtomId, pos.size());
  RDGeom::Point3D rJI = pos[iAtomId] - pos[jAtomId];
  double rJISqLength = rJI.lengthSq();
  if (rJISqLength <= 1.e-16)
    throw ValueErrorException("atoms i and j have identical 3D coordinates");
  RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
  double rJKSqLength = rJK.lengthSq();
  if (rJKSqLength <= 1.e-16)
    throw ValueErrorException("atoms j and k have identical 3D coordinates");
  return rJI.angleTo(rJK);
}

void setAngleRad(Conformer &conf, unsigned int iAtomId, unsigned int jAtomId,
                 unsigned int kAtomId, double value) {
  RDGeom::POINT3D_VECT &pos = conf.getPositions();
  URANGE_CHECK(iAtomId, pos.size());
  URANGE_CHECK(jAtomId, pos.size());
  URANGE_CHECK(kAtomId, pos.size());
  ROMol &mol = conf.getOwningMol();
  Bond *bondJI = mol.getBondBetweenAtoms(jAtomId, iAtomId);
  if (!bondJI) throw ValueErrorException("atoms i and j must be bonded");
  Bond *bondJK = mol.getBondBetweenAtoms(jAtomId, kAtomId);
  if (!bondJK) throw ValueErrorException("atoms j and k must be bonded");
  if (queryIsBondInRing(bondJI) && queryIsBondInRing(bondJK))
    throw ValueErrorException(
        "bonds (i,j) and (j,k) must not both belong to a ring");

  RDGeom::Point3D rJI = pos[iAtomId] - pos[jAtomId];
  double rJISqLength = rJI.lengthSq();
  if (rJISqLength <= 1.e-16)
    throw ValueErrorException("atoms i and j have identical 3D coordinates");
  RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
  double rJKSqLength = rJK.lengthSq();
  if (rJKSqLength <= 1.e-16)
    throw ValueErrorException("atoms j and k have identical 3D coordinates");

  // we only need to rotate by delta with respect to the current angle value
  value -= rJI.angleTo(rJK);
  RDGeom::Point3D &rotAxisBegin = pos[jAtomId];
  // our rotation axis is the normal to the plane of atoms i, j, k
  RDGeom::Point3D rotAxisEnd = rJI.crossProduct(rJK) + pos[jAtomId];
  RDGeom::Point3D rotAxis = rotAxisEnd - rotAxisBegin;
  rotAxis.normalize();
  // get all atoms bonded to j and loop through them
  std::list<unsigned int> alist;
  _toBeMovedIdxList(mol, jAtomId, kAtomId, alist);
  for (std::list<unsigned int>::iterator it = alist.begin(); it != alist.end();
       ++it) {
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

double getDihedralRad(const Conformer &conf, unsigned int iAtomId,
                      unsigned int jAtomId, unsigned int kAtomId,
                      unsigned int lAtomId) {
  const RDGeom::POINT3D_VECT &pos = conf.getPositions();
  URANGE_CHECK(iAtomId, pos.size());
  URANGE_CHECK(jAtomId, pos.size());
  URANGE_CHECK(kAtomId, pos.size());
  URANGE_CHECK(lAtomId, pos.size());
  RDGeom::Point3D rIJ = pos[jAtomId] - pos[iAtomId];
  double rIJSqLength = rIJ.lengthSq();
  if (rIJSqLength <= 1.e-16)
    throw ValueErrorException("atoms i and j have identical 3D coordinates");
  RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
  double rJKSqLength = rJK.lengthSq();
  if (rJKSqLength <= 1.e-16)
    throw ValueErrorException("atoms j and k have identical 3D coordinates");
  RDGeom::Point3D rKL = pos[lAtomId] - pos[kAtomId];
  double rKLSqLength = rKL.lengthSq();
  if (rKLSqLength <= 1.e-16)
    throw ValueErrorException("atoms k and l have identical 3D coordinates");

  RDGeom::Point3D nIJK = rIJ.crossProduct(rJK);
  double nIJKSqLength = nIJK.lengthSq();
  RDGeom::Point3D nJKL = rJK.crossProduct(rKL);
  double nJKLSqLength = nJKL.lengthSq();
  RDGeom::Point3D m = nIJK.crossProduct(rJK);
  // we want a signed dihedral, that's why we use atan2 instead of acos
  return -atan2(m.dotProduct(nJKL) / sqrt(nJKLSqLength * m.lengthSq()),
                nIJK.dotProduct(nJKL) / sqrt(nIJKSqLength * nJKLSqLength));
}

void setDihedralRad(Conformer &conf, unsigned int iAtomId, unsigned int jAtomId,
                    unsigned int kAtomId, unsigned int lAtomId, double value) {
  RDGeom::POINT3D_VECT &pos = conf.getPositions();
  URANGE_CHECK(iAtomId, pos.size());
  URANGE_CHECK(jAtomId, pos.size());
  URANGE_CHECK(kAtomId, pos.size());
  URANGE_CHECK(lAtomId, pos.size());
  ROMol &mol = conf.getOwningMol();
  Bond *bondJK = mol.getBondBetweenAtoms(jAtomId, kAtomId);
  if (!bondJK) throw ValueErrorException("atoms j and k must be bonded");

  if (queryIsBondInRing(bondJK))
    throw ValueErrorException("bond (j,k) must not belong to a ring");
  RDGeom::Point3D rIJ = pos[jAtomId] - pos[iAtomId];
  double rIJSqLength = rIJ.lengthSq();
  if (rIJSqLength <= 1.e-16)
    throw ValueErrorException("atoms i and j have identical 3D coordinates");
  RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
  double rJKSqLength = rJK.lengthSq();
  if (rJKSqLength <= 1.e-16)
    throw ValueErrorException("atoms j and k have identical 3D coordinates");
  RDGeom::Point3D rKL = pos[lAtomId] - pos[kAtomId];
  double rKLSqLength = rKL.lengthSq();
  if (rKLSqLength <= 1.e-16)
    throw ValueErrorException("atoms k and l have identical 3D coordinates");

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
  _toBeMovedIdxList(mol, jAtomId, kAtomId, alist);
  for (std::list<unsigned int>::iterator it = alist.begin(); it != alist.end();
       ++it) {
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
