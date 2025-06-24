//
//   Copyright (C) 2003-2024 Greg Landrum and other RDKit contributors
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

#ifndef RDK_HAS_EIGEN3
constexpr double EIGEN_TOLERANCE = 5.0e-2;
#endif
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

RDGeom::Point3D computeCentroid(const Conformer &conf, bool ignoreHs,
                                const std::vector<double> *weights) {
  PRECONDITION(!weights || weights->size() >= conf.getNumAtoms(),
               "bad weights vector");
  RDGeom::Point3D res(0.0, 0.0, 0.0);
  const ROMol &mol = conf.getOwningMol();
  double wSum = 0.0;
  for (unsigned int i = 0; i < conf.getNumAtoms(); ++i) {
    if (ignoreHs && mol.getAtomWithIdx(i)->getAtomicNum() == 1) {
      continue;
    }
    double w = (weights ? weights->at(i) : 1.0);
    wSum += w;
    res += conf.getAtomPos(i) * w;
  }
  res /= wSum;
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

#ifndef RDK_HAS_EIGEN3
RDNumeric::DoubleSymmMatrix *computeCovarianceMatrix(
    const Conformer &conf, const RDGeom::Point3D &center, bool normalize,
    bool ignoreHs) {
  double xx, xy, xz, yy, yz, zz;
  computeCovarianceTerms(conf, center, xx, xy, xz, yy, yz, zz, normalize,
                         ignoreHs, nullptr);
  auto *res = new RDNumeric::DoubleSymmMatrix(3, 3);
  res->setVal(0, 0, xx);
  res->setVal(0, 1, xy);
  res->setVal(0, 2, xz);
  res->setVal(1, 1, yy);
  res->setVal(1, 2, yz);
  res->setVal(2, 2, zz);
  return res;
}
#endif

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
}  // namespace

#ifdef RDK_HAS_EIGEN3
#include <Eigen/Dense>

namespace {
bool getEigenValEigenVectHelper(Eigen::Matrix3d &eigVecs,
                                Eigen::Vector3d &eigVals, double sumXX,
                                double sumXY, double sumXZ, double sumYY,
                                double sumYZ, double sumZZ) {
  Eigen::Matrix3d mat;
  mat << sumXX, sumXY, sumXZ, sumXY, sumYY, sumYZ, sumXZ, sumYZ, sumZZ;
  // std::cerr<<"getEigenValEigenVectHelper  matrix: "<<mat<<std::endl;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(mat);
  if (eigensolver.info() != Eigen::Success) {
    BOOST_LOG(rdErrorLog) << "eigenvalue calculation did not converge"
                          << std::endl;
    return false;
  }
  eigVecs = eigensolver.eigenvectors();
  // std::cerr<<"getEigenValEigenVectHelper  eigVecs: "<<eigVecs<<std::endl;
  eigVals = eigensolver.eigenvalues();
  // std::cerr<<"getEigenValEigenVectHelper  eigVals: "<<eigVals<<std::endl;
  return true;
}

bool getEigenValEigenVectFromCovMat(const RDKit::Conformer &conf,
                                    Eigen::Matrix3d &eigVecs,
                                    Eigen::Vector3d &eigVals,
                                    const RDGeom::Point3D &origin,
                                    bool ignoreHs, bool normalizeCovar,
                                    const std::vector<double> *weights) {
  PRECONDITION((!weights || weights->size() >= conf.getNumAtoms()),
               "bad weights vector");

  // std::cerr << "getEigenValEigenVectFromCovMat ignoreHs " << ignoreHs << "
  // normalizeCovar " << normalizeCovar << " weights " << weights << " origin "
  // << origin.x << "," << origin.y << ","  << origin.z << std::endl;
  double sumXX, sumXY, sumXZ, sumYY, sumYZ, sumZZ;
  computeCovarianceTerms(conf, origin, sumXX, sumXY, sumXZ, sumYY, sumYZ, sumZZ,
                         normalizeCovar, ignoreHs, weights);

  return getEigenValEigenVectHelper(eigVecs, eigVals, sumXX, sumXY, sumXZ,
                                    sumYY, sumYZ, sumZZ);
}
}  // namespace

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
  const ROMol &mol = conf.getOwningMol();
  if (!weights && !force && mol.hasProp(axesPropName) &&
      mol.hasProp(momentsPropName)) {
    mol.getProp(axesPropName, axes);
    mol.getProp(momentsPropName, moments);
    return true;
  }
  auto origin = computeCentroid(conf, ignoreHs, weights);

  double sumXX, sumXY, sumXZ, sumYY, sumYZ, sumZZ;
  computeInertiaTerms(conf, origin, sumXX, sumXY, sumXZ, sumYY, sumYZ, sumZZ,
                      ignoreHs, weights);

  if (!getEigenValEigenVectHelper(axes, moments, sumXX, sumXY, sumXZ, sumYY,
                                  sumYZ, sumZZ)) {
    return false;
  }

  if (!weights) {
    mol.setProp(axesPropName, axes, true);
    mol.setProp(momentsPropName, moments, true);
  }
  return true;
}

bool computePrincipalAxesAndMomentsFromGyrationMatrix(
    const RDKit::Conformer &conf, Eigen::Matrix3d &axes,
    Eigen::Vector3d &moments, bool ignoreHs, bool force,
    const std::vector<double> *weights) {
  const char *axesPropName =
      ignoreHs ? "_principalAxes_noH_cov" : "_principalAxes_cov";
  const char *momentsPropName =
      ignoreHs ? "_principalMoments_noH_cov" : "_principalMoments_cov";
  const ROMol &mol = conf.getOwningMol();
  if (!weights && !force && mol.hasProp(axesPropName) &&
      mol.hasProp(momentsPropName)) {
    mol.getProp(axesPropName, axes);
    mol.getProp(momentsPropName, moments);
    return true;
  }
  auto origin = computeCentroid(conf, ignoreHs, weights);
  bool res = getEigenValEigenVectFromCovMat(conf, axes, moments, origin,
                                            ignoreHs, true, weights);
  if (res && !weights) {
    conf.getOwningMol().setProp(axesPropName, axes, true);
    conf.getOwningMol().setProp(momentsPropName, moments, true);
  }
  return res;
}

RDGeom::Transform3D *computeCanonicalTransform(const Conformer &conf,
                                               const RDGeom::Point3D *center,
                                               bool normalizeCovar,
                                               bool ignoreHs) {
  constexpr unsigned int DIM = 3;
  RDGeom::Point3D origin;
  if (!center) {
    origin = computeCentroid(conf, ignoreHs);
  } else {
    origin = (*center);
  }
  unsigned int nAtms = conf.getNumAtoms();
  auto *trans = new RDGeom::Transform3D;
  trans->setToIdentity();

  // if we have a single atom system we don't need to do anyhting setting
  // translation is sufficient
  if (nAtms > 1) {
    Eigen::Matrix3d eigVecs;
    Eigen::Vector3d eigVals;
    if (getEigenValEigenVectFromCovMat(conf, eigVecs, eigVals, origin, ignoreHs,
                                       normalizeCovar, nullptr)) {
      std::vector<std::pair<unsigned int, double>> eigValsSorted;
      CHECK_INVARIANT(eigVals.size() == DIM, "less eigenvalues than expected");
      eigValsSorted.reserve(DIM);
      for (unsigned int i = 0; i < DIM; ++i) {
        eigValsSorted.emplace_back(i, eigVals(i));
      }
      std::sort(eigValsSorted.begin(), eigValsSorted.end(),
                [](const std::pair<unsigned int, double> &a,
                   const std::pair<unsigned int, double> &b) {
                  return (a.second > b.second);
                });
      for (unsigned int col = 0; col < DIM; ++col) {
        unsigned int colSorted = eigValsSorted.at(col).first;
        for (unsigned int row = 0; row < DIM; ++row) {
          trans->setVal(col, row, eigVecs(row, colSorted));
        }
      }
    }
  }
  origin *= -1.0;

  // In some situations we can end up with one or more negative values on the
  //  diagonal. An odd number of these will result in an inversion of the
  //  structure, so we need to check for that and, if necessary, correct by
  //  negating one row
  if (trans->getVal(0, 0) * trans->getVal(1, 1) * trans->getVal(2, 2) < 0) {
    for (auto i = 0; i < 3; ++i) {
      trans->setVal(2, i, trans->getVal(2, i) * -1);
    }
  }

  trans->TransformPoint(origin);
  trans->SetTranslation(origin);

  return trans;
}
#else
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
  auto *trans = new RDGeom::Transform3D;

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

  // In some situations we can end up with one or more negative values on the
  //  diagonal. An odd number of these will result in an inversion of the
  //  structure, so we need to check for that and, if necessary, correct by
  //  negating one row
  if (trans->getVal(0, 0) * trans->getVal(1, 1) * trans->getVal(2, 2) < 0) {
    for (auto i = 0; i < 3; ++i) {
      trans->setVal(2, i, trans->getVal(2, i) * -1);
    }
  }
  trans->TransformPoint(origin);

  trans->SetTranslation(origin);
  delete covMat;

  return trans;
}
#endif

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
    canonicalizeConformer(*(*ci), nullptr, normalizeCovar, ignoreHs);
  }
}

void transformConformer(Conformer &conf, const RDGeom::Transform3D &trans) {
  RDGeom::POINT3D_VECT &positions = conf.getPositions();
  RDGeom::POINT3D_VECT_I pi;
  for (pi = positions.begin(); pi != positions.end(); ++pi) {
    trans.TransformPoint(*pi);
  }
}

void transformMolSubstanceGroups(ROMol &mol, const RDGeom::Transform3D &trans) {
  auto &sgs = getSubstanceGroups(mol);
  for (auto &sg : sgs) {
    for (auto &brk : sg.getBrackets()) {
      trans.TransformPoint(brk[0]);
      trans.TransformPoint(brk[1]);
      trans.TransformPoint(brk[2]);
    }
    for (auto &cs : sg.getCStates()) {
      trans.TransformPoint(cs.vector);
    }
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
      wIdx = (mol[*nbrIdx])->getIdx();
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
}  // namespace
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
  if (!bond) {
    throw ValueErrorException("atoms i and j must be bonded");
  }
  if (queryIsBondInRing(bond)) {
    throw ValueErrorException("bond (i,j) must not belong to a ring");
  }
  RDGeom::Point3D v = pos[iAtomId] - pos[jAtomId];
  double origValue = v.length();
  if (origValue <= 1.e-8) {
    throw ValueErrorException("atoms i and j have identical 3D coordinates");
  }

  // get all atoms bonded to j
  std::list<unsigned int> alist;
  _toBeMovedIdxList(mol, iAtomId, jAtomId, alist);
  v *= (value / origValue - 1.);
  for (unsigned int &it : alist) {
    pos[it] -= v;
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
  if (rJISqLength <= 1.e-16) {
    throw ValueErrorException("atoms i and j have identical 3D coordinates");
  }
  RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
  double rJKSqLength = rJK.lengthSq();
  if (rJKSqLength <= 1.e-16) {
    throw ValueErrorException("atoms j and k have identical 3D coordinates");
  }
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
  if (!bondJI) {
    throw ValueErrorException("atoms i and j must be bonded");
  }
  Bond *bondJK = mol.getBondBetweenAtoms(jAtomId, kAtomId);
  if (!bondJK) {
    throw ValueErrorException("atoms j and k must be bonded");
  }
  if (queryIsBondInRing(bondJI) && queryIsBondInRing(bondJK)) {
    throw ValueErrorException(
        "bonds (i,j) and (j,k) must not both belong to a ring");
  }

  RDGeom::Point3D rJI = pos[iAtomId] - pos[jAtomId];
  double rJISqLength = rJI.lengthSq();
  if (rJISqLength <= 1.e-16) {
    throw ValueErrorException("atoms i and j have identical 3D coordinates");
  }
  RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
  double rJKSqLength = rJK.lengthSq();
  if (rJKSqLength <= 1.e-16) {
    throw ValueErrorException("atoms j and k have identical 3D coordinates");
  }

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
  for (unsigned int &it : alist) {
    // translate atom so that it coincides with the origin of rotation
    pos[it] -= rotAxisBegin;
    // rotate around our rotation axis
    RDGeom::Transform3D rotByAngle;
    rotByAngle.SetRotation(value, rotAxis);
    rotByAngle.TransformPoint(pos[it]);
    // translate atom back
    pos[it] += rotAxisBegin;
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
  if (rIJSqLength <= 1.e-16) {
    throw ValueErrorException("atoms i and j have identical 3D coordinates");
  }
  RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
  double rJKSqLength = rJK.lengthSq();
  if (rJKSqLength <= 1.e-16) {
    throw ValueErrorException("atoms j and k have identical 3D coordinates");
  }
  RDGeom::Point3D rKL = pos[lAtomId] - pos[kAtomId];
  double rKLSqLength = rKL.lengthSq();
  if (rKLSqLength <= 1.e-16) {
    throw ValueErrorException("atoms k and l have identical 3D coordinates");
  }

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
  if (!bondJK) {
    throw ValueErrorException("atoms j and k must be bonded");
  }

  if (queryIsBondInRing(bondJK)) {
    throw ValueErrorException("bond (j,k) must not belong to a ring");
  }
  RDGeom::Point3D rIJ = pos[jAtomId] - pos[iAtomId];
  double rIJSqLength = rIJ.lengthSq();
  if (rIJSqLength <= 1.e-16) {
    throw ValueErrorException("atoms i and j have identical 3D coordinates");
  }
  RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
  double rJKSqLength = rJK.lengthSq();
  if (rJKSqLength <= 1.e-16) {
    throw ValueErrorException("atoms j and k have identical 3D coordinates");
  }
  RDGeom::Point3D rKL = pos[lAtomId] - pos[kAtomId];
  double rKLSqLength = rKL.lengthSq();
  if (rKLSqLength <= 1.e-16) {
    throw ValueErrorException("atoms k and l have identical 3D coordinates");
  }

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
  for (unsigned int &it : alist) {
    // translate atom so that it coincides with the origin of rotation
    pos[it] -= rotAxisBegin;
    // rotate around our rotation axis
    RDGeom::Transform3D rotByAngle;
    rotByAngle.SetRotation(value, rotAxis);
    rotByAngle.TransformPoint(pos[it]);
    // translate atom back
    pos[it] += rotAxisBegin;
  }
}
}  // namespace MolTransforms
