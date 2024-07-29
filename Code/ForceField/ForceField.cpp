// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "ForceField.h"
#include "Contrib.h"

#include <RDGeneral/Invariant.h>
#include <Numerics/Optimizer/BFGSOpt.h>

namespace RDKit {
namespace ForceFieldsHelper {
void normalizeAngleDeg(double &angleDeg) {
  angleDeg = fmod(angleDeg, 360.0);
  if (angleDeg < -180.0) {
    angleDeg += 360.0;
  } else if (angleDeg > 180.0) {
    angleDeg -= 360.0;
  }
}

void computeDihedral(const RDGeom::PointPtrVect &pos, unsigned int idx1,
                     unsigned int idx2, unsigned int idx3, unsigned int idx4,
                     double *dihedral, double *cosPhi, RDGeom::Point3D r[4],
                     RDGeom::Point3D t[2], double d[2]) {
  computeDihedral(static_cast<RDGeom::Point3D *>(pos[idx1]),
                  static_cast<RDGeom::Point3D *>(pos[idx2]),
                  static_cast<RDGeom::Point3D *>(pos[idx3]),
                  static_cast<RDGeom::Point3D *>(pos[idx4]), dihedral, cosPhi,
                  r, t, d);
}

void computeDihedral(const double *pos, unsigned int idx1, unsigned int idx2,
                     unsigned int idx3, unsigned int idx4, double *dihedral,
                     double *cosPhi, RDGeom::Point3D r[4], RDGeom::Point3D t[2],
                     double d[2]) {
  RDGeom::Point3D p1(pos[3 * idx1], pos[3 * idx1 + 1], pos[3 * idx1 + 2]);
  RDGeom::Point3D p2(pos[3 * idx2], pos[3 * idx2 + 1], pos[3 * idx2 + 2]);
  RDGeom::Point3D p3(pos[3 * idx3], pos[3 * idx3 + 1], pos[3 * idx3 + 2]);
  RDGeom::Point3D p4(pos[3 * idx4], pos[3 * idx4 + 1], pos[3 * idx4 + 2]);
  computeDihedral(&p1, &p2, &p3, &p4, dihedral, cosPhi, r, t, d);
}

void computeDihedral(const RDGeom::Point3D *p1, const RDGeom::Point3D *p2,
                     const RDGeom::Point3D *p3, const RDGeom::Point3D *p4,
                     double *dihedral, double *cosPhi, RDGeom::Point3D r[4],
                     RDGeom::Point3D t[2], double d[2]) {
  PRECONDITION(p1, "p1 must not be null");
  PRECONDITION(p2, "p2 must not be null");
  PRECONDITION(p3, "p3 must not be null");
  PRECONDITION(p4, "p4 must not be null");
  RDGeom::Point3D rLocal[4];
  RDGeom::Point3D tLocal[2];
  double dLocal[2];
  if (!r) {
    r = rLocal;
  }
  if (!t) {
    t = tLocal;
  }
  if (!d) {
    d = dLocal;
  }
  r[0] = *p1 - *p2;
  r[1] = *p3 - *p2;
  r[2] = -r[1];
  r[3] = *p4 - *p3;

  t[0] = r[0].crossProduct(r[1]);
  d[0] = (std::max)(t[0].length(), 1.0e-5);
  t[0] /= d[0];
  t[1] = r[2].crossProduct(r[3]);
  d[1] = (std::max)(t[1].length(), 1.0e-5);
  t[1] /= d[1];
  double cosPhiLocal;
  if (!cosPhi) {
    cosPhi = &cosPhiLocal;
  }
  *cosPhi = (std::max)(-1.0, (std::min)(t[0].dotProduct(t[1]), 1.0));
  // we want a signed dihedral, that's why we use atan2 instead of acos
  if (dihedral) {
    RDGeom::Point3D m = t[0].crossProduct(r[1]);
    double mLength = (std::max)(m.length(), 1.0e-5);
    *dihedral = -atan2(m.dotProduct(t[1]) / mLength, *cosPhi);
  }
}
}  // namespace ForceFieldsHelper
}  // namespace RDKit

namespace ForceFieldsHelper {
class calcEnergy {
 public:
  calcEnergy(ForceFields::ForceField *ffHolder) : mp_ffHolder(ffHolder){};
  double operator()(double *pos) const { return mp_ffHolder->calcEnergy(pos); }

 private:
  ForceFields::ForceField *mp_ffHolder;
};

class calcGradient {
 public:
  calcGradient(ForceFields::ForceField *ffHolder) : mp_ffHolder(ffHolder){};
  double operator()(double *pos, double *grad) const {
    double res = 1.0;
    // the contribs to the gradient function use +=, so we need
    // to zero the grad out before moving on:
    for (unsigned int i = 0;
         i < mp_ffHolder->numPoints() * mp_ffHolder->dimension(); i++) {
      grad[i] = 0.0;
    }
    mp_ffHolder->calcGrad(pos, grad);

#if 1
    // FIX: this hack reduces the gradients so that the
    // minimizer is more efficient.
    double maxGrad = -1e8;
    double gradScale = 0.1;
    for (unsigned int i = 0;
         i < mp_ffHolder->numPoints() * mp_ffHolder->dimension(); i++) {
      grad[i] *= gradScale;
      if (grad[i] > maxGrad) {
        maxGrad = grad[i];
      }
    }
    // this is a continuation of the same hack to avoid
    // some potential numeric instabilities:
    if (maxGrad > 10.0) {
      while (maxGrad * gradScale > 10.0) {
        gradScale *= .5;
      }
      for (unsigned int i = 0;
           i < mp_ffHolder->numPoints() * mp_ffHolder->dimension(); i++) {
        grad[i] *= gradScale;
      }
    }
    res = gradScale;
#endif
    return res;
  }

 private:
  ForceFields::ForceField *mp_ffHolder;
};
}  // namespace ForceFieldsHelper

namespace ForceFields {
ForceField::~ForceField() {
  d_numPoints = 0;
  d_positions.clear();
  d_contribs.clear();
  delete[] dp_distMat;
  dp_distMat = nullptr;
}

ForceField::ForceField(const ForceField &other)
    : d_dimension(other.d_dimension),
      df_init(false),
      d_numPoints(other.d_numPoints),
      dp_distMat(nullptr) {
  d_contribs.clear();
  for (const auto &contrib : other.d_contribs) {
    ForceFieldContrib *ncontrib = contrib->copy();
    ncontrib->dp_forceField = this;
    d_contribs.push_back(ContribPtr(ncontrib));
  }
};

double ForceField::distance(unsigned int i, unsigned int j, double *pos) {
  PRECONDITION(df_init, "not initialized");
  URANGE_CHECK(i, d_numPoints);
  URANGE_CHECK(j, d_numPoints);
  if (j < i) {
    int tmp = j;
    j = i;
    i = tmp;
  }
  unsigned int idx = i + j * (j + 1) / 2;
  CHECK_INVARIANT(idx < d_matSize, "Bad index");
  double &res = dp_distMat[idx];
  if (res < 0.0) {
    // we need to calculate this distance:
    if (!pos) {
      res = 0.0;
      for (unsigned int idx = 0; idx < d_dimension; ++idx) {
        double tmp =
            (*(this->positions()[i]))[idx] - (*(this->positions()[j]))[idx];
        res += tmp * tmp;
      }
    } else {
      res = 0.0;
#if 0
      for (unsigned int idx = 0; idx < d_dimension; idx++) {
        double tmp = pos[d_dimension * i + idx] - pos[d_dimension * j + idx];
        res += tmp * tmp;
      }
#else
      double *pi = &(pos[d_dimension * i]), *pj = &(pos[d_dimension * j]);
      for (unsigned int idx = 0; idx < d_dimension; ++idx, ++pi, ++pj) {
        double tmp = *pi - *pj;
        res += tmp * tmp;
      }
#endif
    }
    res = sqrt(res);
  }
  return res;
}

double ForceField::distance2(unsigned int i, unsigned int j,
                             double *pos) const {
  PRECONDITION(df_init, "not initialized");
  URANGE_CHECK(i, d_numPoints);
  URANGE_CHECK(j, d_numPoints);
  if (j < i) {
    int tmp = j;
    j = i;
    i = tmp;
  }
  double res;
  if (!pos) {
    res = 0.0;
    for (unsigned int idx = 0; idx < d_dimension; ++idx) {
      double tmp =
          (*(this->positions()[i]))[idx] - (*(this->positions()[j]))[idx];
      res += tmp * tmp;
    }
  } else {
    res = 0.0;
#if 0
    for (unsigned int idx = 0; idx < d_dimension; idx++) {
      double tmp = pos[d_dimension * i + idx] - pos[d_dimension * j + idx];
      res += tmp * tmp;
    }
#else
    double *pi = &(pos[d_dimension * i]), *pj = &(pos[d_dimension * j]);
    for (unsigned int idx = 0; idx < d_dimension; ++idx, ++pi, ++pj) {
      double tmp = *pi - *pj;
      res += tmp * tmp;
    }
#endif
  }
  return res;
}
double ForceField::distance(unsigned int i, unsigned int j, double *pos) const {
  auto res = sqrt(distance2(i, j, pos));
  return res;
}

void ForceField::initialize() {
  // clean up if we have used this already:
  df_init = false;
  delete[] dp_distMat;
  dp_distMat = nullptr;

  d_numPoints = d_positions.size();
  d_matSize = d_numPoints * (d_numPoints + 1) / 2;
  dp_distMat = new double[d_matSize];
  this->initDistanceMatrix();
  df_init = true;
}

int ForceField::minimize(unsigned int maxIts, double forceTol,
                         double energyTol) {
  return minimize(0, nullptr, maxIts, forceTol, energyTol);
}

int ForceField::minimize(unsigned int snapshotFreq,
                         RDKit::SnapshotVect *snapshotVect, unsigned int maxIts,
                         double forceTol, double energyTol) {
  PRECONDITION(df_init, "not initialized");
  PRECONDITION(static_cast<unsigned int>(d_numPoints) == d_positions.size(),
               "size mismatch");
  if (d_contribs.empty()) {
    return 0;
  }

  unsigned int numIters = 0;
  unsigned int dim = this->d_numPoints * d_dimension;
  double finalForce = 0.0;
  auto *points = new double[dim];

  this->scatter(points);
  ForceFieldsHelper::calcEnergy eCalc(this);
  ForceFieldsHelper::calcGradient gCalc(this);

  int res =
      BFGSOpt::minimize(dim, points, forceTol, numIters, finalForce, eCalc,
                        gCalc, snapshotFreq, snapshotVect, energyTol, maxIts);
  this->gather(points);

  delete[] points;
  return res;
}

double ForceField::calcEnergy(std::vector<double> *contribs) const {
  PRECONDITION(df_init, "not initialized");
  double res = 0.0;
  if (d_contribs.empty()) {
    return res;
  }
  if (contribs) {
    contribs->clear();
    contribs->reserve(d_contribs.size());
  }

  unsigned int N = d_positions.size();
  auto *pos = new double[d_dimension * N];
  this->scatter(pos);
  // now loop over the contribs
  for (const auto &d_contrib : d_contribs) {
    double e = d_contrib->getEnergy(pos);
    res += e;
    if (contribs) {
      contribs->push_back(e);
    }
  }
  delete[] pos;
  return res;
}

double ForceField::calcEnergy(double *pos) {
  PRECONDITION(df_init, "not initialized");
  PRECONDITION(pos, "bad position vector");
  double res = 0.0;

  this->initDistanceMatrix();
  if (d_contribs.empty()) {
    return res;
  }

  // now loop over the contribs
  for (ContribPtrVect::const_iterator contrib = d_contribs.begin();
       contrib != d_contribs.end(); contrib++) {
    double E = (*contrib)->getEnergy(pos);
    res += E;
  }
  return res;
}

void ForceField::calcGrad(double *grad) const {
  PRECONDITION(df_init, "not initialized");
  PRECONDITION(grad, "bad gradient vector");
  if (d_contribs.empty()) {
    return;
  }

  unsigned int N = d_positions.size();
  auto *pos = new double[d_dimension * N];
  this->scatter(pos);
  for (const auto &d_contrib : d_contribs) {
    d_contrib->getGrad(pos, grad);
  }
  // zero out gradient values for any fixed points:
  for (int d_fixedPoint : d_fixedPoints) {
    CHECK_INVARIANT(static_cast<unsigned int>(d_fixedPoint) < d_numPoints,
                    "bad fixed point index");
    unsigned int idx = d_dimension * d_fixedPoint;
    for (unsigned int di = 0; di < this->dimension(); ++di) {
      grad[idx + di] = 0.0;
    }
  }
  delete[] pos;
}
void ForceField::calcGrad(double *pos, double *grad) {
  PRECONDITION(df_init, "not initialized");
  PRECONDITION(pos, "bad position vector");
  PRECONDITION(grad, "bad gradient vector");
  if (d_contribs.empty()) {
    return;
  }

  for (ContribPtrVect::const_iterator contrib = d_contribs.begin();
       contrib != d_contribs.end(); contrib++) {
    (*contrib)->getGrad(pos, grad);
  }

  for (INT_VECT::const_iterator it = d_fixedPoints.begin();
       it != d_fixedPoints.end(); it++) {
    CHECK_INVARIANT(static_cast<unsigned int>(*it) < d_numPoints,
                    "bad fixed point index");
    unsigned int idx = d_dimension * (*it);
    for (unsigned int di = 0; di < this->dimension(); ++di) {
      grad[idx + di] = 0.0;
    }
  }
}

void ForceField::scatter(double *pos) const {
  PRECONDITION(df_init, "not initialized");
  PRECONDITION(pos, "bad position vector");

  unsigned int tab = 0;
  for (auto d_position : d_positions) {
    for (unsigned int di = 0; di < this->dimension(); ++di) {
      pos[tab + di] = (*d_position)[di];  //->x;
    }
    tab += this->dimension();
  }
  POSTCONDITION(tab == this->dimension() * d_positions.size(), "bad index");
}

void ForceField::gather(double *pos) {
  PRECONDITION(df_init, "not initialized");
  PRECONDITION(pos, "bad position vector");

  unsigned int tab = 0;
  for (auto &d_position : d_positions) {
    for (unsigned int di = 0; di < this->dimension(); ++di) {
      (*d_position)[di] = pos[tab + di];
    }
    tab += this->dimension();
  }
}

void ForceField::initDistanceMatrix() {
  PRECONDITION(d_numPoints, "no points");
  PRECONDITION(dp_distMat, "no distance matrix");
  PRECONDITION(static_cast<unsigned int>(d_numPoints * (d_numPoints + 1) / 2) <=
                   d_matSize,
               "matrix size mismatch");
  for (unsigned int i = 0; i < d_numPoints * (d_numPoints + 1) / 2; i++) {
    dp_distMat[i] = -1.0;
  }
}
}  // namespace ForceFields
