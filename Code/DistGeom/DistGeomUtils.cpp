//
//  Copyright (C) 2004-2019 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "BoundsMatrix.h"
#include "DistGeomUtils.h"
#include "DistViolationContrib.h"
#include "ChiralViolationContrib.h"
#include "FourthDimContrib.h"
#include <Numerics/Matrix.h>
#include <Numerics/SymmMatrix.h>
#include <Numerics/Vector.h>
#include <RDGeneral/Invariant.h>
#include <Numerics/EigenSolvers/PowerEigenSolver.h>
#include <RDGeneral/utils.h>
#include <ForceField/ForceField.h>
#include <ForceField/UFF/DistanceConstraint.h>
#include <ForceField/UFF/AngleConstraint.h>
#include <ForceField/UFF/Inversion.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionAngleM6.h>
#include <boost/dynamic_bitset.hpp>
#include <ForceField/MMFF/Nonbonded.h>

namespace DistGeom {
const double EIGVAL_TOL = 0.001;

double pickRandomDistMat(const BoundsMatrix &mmat,
                         RDNumeric::SymmMatrix<double> &distMat, int seed) {
  if (seed > 0) {
    RDKit::getRandomGenerator(seed);
  }
  return pickRandomDistMat(mmat, distMat, RDKit::getDoubleRandomSource());
}

double pickRandomDistMat(const BoundsMatrix &mmat,
                         RDNumeric::SymmMatrix<double> &distMat,
                         RDKit::double_source_type &rng) {
  // make sure the sizes match up
  unsigned int npt = mmat.numRows();
  CHECK_INVARIANT(npt == distMat.numRows(), "Size mismatch");

  double largestVal = -1.0;
  double *ddata = distMat.getData();
  for (unsigned int i = 1; i < npt; i++) {
    unsigned int id = i * (i + 1) / 2;
    for (unsigned int j = 0; j < i; j++) {
      double ub = mmat.getUpperBound(i, j);
      double lb = mmat.getLowerBound(i, j);
      CHECK_INVARIANT(ub >= lb, "");
      double rval = rng();
      // std::cerr<<i<<"-"<<j<<": "<<rval<<std::endl;
      double d = lb + (rval) * (ub - lb);
      ddata[id + j] = d;
      if (d > largestVal) {
        largestVal = d;
      }
    }
  }
  return largestVal;
}

bool computeInitialCoords(const RDNumeric::SymmMatrix<double> &distMat,
                          RDGeom::PointPtrVect &positions, bool randNegEig,
                          unsigned int numZeroFail, int seed) {
  if (seed > 0) {
    RDKit::getRandomGenerator(seed);
  }
  return computeInitialCoords(distMat, positions,
                              RDKit::getDoubleRandomSource(), randNegEig,
                              numZeroFail);
}
bool computeInitialCoords(const RDNumeric::SymmMatrix<double> &distMat,
                          RDGeom::PointPtrVect &positions,
                          RDKit::double_source_type &rng, bool randNegEig,
                          unsigned int numZeroFail) {
  unsigned int N = distMat.numRows();
  unsigned int nPt = positions.size();
  CHECK_INVARIANT(nPt == N, "Size mismatch");

  unsigned int dim = positions.front()->dimension();

  const double *data = distMat.getData();
  RDNumeric::SymmMatrix<double> sqMat(N), T(N, 0.0);
  RDNumeric::DoubleMatrix eigVecs(dim, N);
  RDNumeric::DoubleVector eigVals(dim);

  double *sqDat = sqMat.getData();

  unsigned int dSize = distMat.getDataSize();
  double sumSqD2 = 0.0;
  for (unsigned int i = 0; i < dSize; i++) {
    sqDat[i] = data[i] * data[i];
    sumSqD2 += sqDat[i];
  }
  sumSqD2 /= (N * N);

  RDNumeric::DoubleVector sqD0i(N, 0.0);
  double *sqD0iData = sqD0i.getData();
  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = 0; j < N; j++) {
      sqD0iData[i] += sqMat.getVal(i, j);
    }
    sqD0iData[i] /= N;
    sqD0iData[i] -= sumSqD2;

    if ((sqD0iData[i] < EIGVAL_TOL) && (N > 3)) {
      return false;
    }
  }

  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = 0; j <= i; j++) {
      double val = 0.5 * (sqD0iData[i] + sqD0iData[j] - sqMat.getVal(i, j));
      T.setVal(i, j, val);
    }
  }
  unsigned int nEigs = (dim < N) ? dim : N;
  RDNumeric::EigenSolvers::powerEigenSolver(nEigs, T, eigVals, eigVecs,
                                            (int)(sumSqD2 * N));

  double *eigData = eigVals.getData();
  bool foundNeg = false;
  unsigned int zeroEigs = 0;
  for (unsigned int i = 0; i < dim; i++) {
    if (eigData[i] > EIGVAL_TOL) {
      eigData[i] = sqrt(eigData[i]);
    } else if (fabs(eigData[i]) < EIGVAL_TOL) {
      eigData[i] = 0.0;
      zeroEigs++;
    } else {
      foundNeg = true;
    }
  }
  if ((foundNeg) && (!randNegEig)) {
    return false;
  }

  if ((zeroEigs >= numZeroFail) && (N > 3)) {
    return false;
  }

  for (unsigned int i = 0; i < N; i++) {
    RDGeom::Point *pt = positions[i];
    for (unsigned int j = 0; j < dim; ++j) {
      if (eigData[j] >= 0.0) {
        (*pt)[j] = eigData[j] * eigVecs.getVal(j, i);
      } else {
        // std::cerr<<"!!! "<<i<<"-"<<j<<": "<<eigData[j]<<std::endl;
        (*pt)[j] = 1.0 - 2.0 * rng();
      }
    }
  }
  return true;
}

bool computeRandomCoords(RDGeom::PointPtrVect &positions, double boxSize,
                         int seed) {
  if (seed > 0) {
    RDKit::getRandomGenerator(seed);
  }
  return computeRandomCoords(positions, boxSize,
                             RDKit::getDoubleRandomSource());
}
bool computeRandomCoords(RDGeom::PointPtrVect &positions, double boxSize,
                         RDKit::double_source_type &rng) {
  CHECK_INVARIANT(boxSize > 0.0, "bad boxSize");

  for (auto pt : positions) {
    for (unsigned int i = 0; i < pt->dimension(); ++i) {
      (*pt)[i] = boxSize * (rng() - 0.5);
    }
  }
  return true;
}

ForceFields::ForceField *constructForceField(
    const BoundsMatrix &mmat, RDGeom::PointPtrVect &positions,
    const VECT_CHIRALSET &csets, double weightChiral, double weightFourthDim,
    std::map<std::pair<int, int>, double> *extraWeights, double basinSizeTol,
    boost::dynamic_bitset<> *fixedPts) {
  unsigned int N = mmat.numRows();
  CHECK_INVARIANT(N == positions.size(), "");
  auto *field = new ForceFields::ForceField(positions[0]->dimension());
  for (unsigned int i = 0; i < N; i++) {
    field->positions().push_back(positions[i]);
  }

  for (unsigned int i = 1; i < N; i++) {
    for (unsigned int j = 0; j < i; j++) {
      if (fixedPts != nullptr && (*fixedPts)[i] && (*fixedPts)[j]) {
        continue;
      }
      double w = 1.0;
      double l = mmat.getLowerBound(i, j);
      double u = mmat.getUpperBound(i, j);
      bool includeIt = false;
      if (extraWeights) {
        std::map<std::pair<int, int>, double>::const_iterator mapIt;
        mapIt = extraWeights->find(std::make_pair(i, j));
        if (mapIt != extraWeights->end()) {
          w = mapIt->second;
          includeIt = true;
        }
      }
      if (u - l <= basinSizeTol) {
        includeIt = true;
      }
      if (includeIt) {
        auto *contrib = new DistViolationContrib(field, i, j, u, l, w);
        field->contribs().push_back(ForceFields::ContribPtr(contrib));
      }
    }
  }

  // now add chiral constraints
  if (weightChiral > 1.e-8) {
    for (const auto &cset : csets) {
      auto *contrib =
          new ChiralViolationContrib(field, cset.get(), weightChiral);
      field->contribs().push_back(ForceFields::ContribPtr(contrib));
    }
  }

  // finally the contribution from the fourth dimension if we need to
  if ((field->dimension() == 4) && (weightFourthDim > 1.e-8)) {
    for (unsigned int i = 0; i < N; i++) {
      auto *contrib = new FourthDimContrib(field, i, weightFourthDim);
      field->contribs().push_back(ForceFields::ContribPtr(contrib));
    }
  }
  return field;
}  // constructForceField

ForceFields::ForceField *construct3DForceField(
    const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails) {
  unsigned int N = mmat.numRows();
  CHECK_INVARIANT(N == positions.size(), "");
  CHECK_INVARIANT(etkdgDetails.expTorsionAtoms.size() ==
                      etkdgDetails.expTorsionAngles.size(),
                  "");
  auto *field = new ForceFields::ForceField(positions[0]->dimension());
  for (unsigned int i = 0; i < N; ++i) {
    field->positions().push_back(positions[i]);
  }

  // keep track which atoms are 1,2- or 1,3-restrained
  boost::dynamic_bitset<> atomPairs(N * N);

  // torsion constraints
  for (unsigned int t = 0; t < etkdgDetails.expTorsionAtoms.size(); ++t) {
    int i = etkdgDetails.expTorsionAtoms[t][0];
    int j = etkdgDetails.expTorsionAtoms[t][1];
    int k = etkdgDetails.expTorsionAtoms[t][2];
    int l = etkdgDetails.expTorsionAtoms[t][3];
    if (i < j) {
      atomPairs[i * N + j] = 1;
    } else {
      atomPairs[j * N + i] = 1;
    }
    // etkdgDetails.expTorsionAngles[t][0] = (signs, V's)
    auto *contrib = new ForceFields::CrystalFF::TorsionAngleContribM6(
        field, i, j, k, l, etkdgDetails.expTorsionAngles[t].second,
        etkdgDetails.expTorsionAngles[t].first);
    field->contribs().push_back(ForceFields::ContribPtr(contrib));
  }  // torsion constraints

  // improper torsions / out-of-plane bend / inversion
  double oobForceScalingFactor = 10.0;
  for (const auto &improperAtom : etkdgDetails.improperAtoms) {
    std::vector<int> n(4);
    for (unsigned int i = 0; i < 3; ++i) {
      n[1] = 1;
      switch (i) {
        case 0:
          n[0] = 0;
          n[2] = 2;
          n[3] = 3;
          break;

        case 1:
          n[0] = 0;
          n[2] = 3;
          n[3] = 2;
          break;

        case 2:
          n[0] = 2;
          n[2] = 3;
          n[3] = 0;
          break;
      }
      auto *contrib = new ForceFields::UFF::InversionContrib(
          field, improperAtom[n[0]], improperAtom[n[1]], improperAtom[n[2]],
          improperAtom[n[3]], improperAtom[4],
          static_cast<bool>(improperAtom[5]), oobForceScalingFactor);
      field->contribs().push_back(ForceFields::ContribPtr(contrib));
    }
  }

  double fdist = 100.0;  // force constant
  // 1,2 distance constraints
  for (const auto &bnd : etkdgDetails.bonds) {
    unsigned int i = bnd.first;
    unsigned int j = bnd.second;
    if (i < j) {
      atomPairs[i * N + j] = 1;
    } else {
      atomPairs[j * N + i] = 1;
    }
    double d = ((*positions[i]) - (*positions[j])).length();
    double l = d - 0.01;
    double u = d + 0.01;
    auto *contrib = new ForceFields::UFF::DistanceConstraintContrib(
        field, i, j, l, u, fdist);
    field->contribs().push_back(ForceFields::ContribPtr(contrib));
  }

  // 1,3 distance constraints
  for (const auto &angle : etkdgDetails.angles) {
    unsigned int i = angle[0];
    unsigned int j = angle[1];
    unsigned int k = angle[2];
    if (i < j) {
      atomPairs[i * N + j] = 1;
    } else {
      atomPairs[j * N + i] = 1;
    }
    // check for triple bonds
    if (angle[3]) {
      auto *contrib = new ForceFields::UFF::AngleConstraintContrib(
          field, i, j, k, 179.0, 180.0, fdist);
      field->contribs().push_back(ForceFields::ContribPtr(contrib));
    } else {
      double d = ((*positions[i]) - (*positions[k])).length();
      double l = d - 0.01;
      double u = d + 0.01;
      auto *contrib = new ForceFields::UFF::DistanceConstraintContrib(
          field, i, k, l, u, fdist);
      field->contribs().push_back(ForceFields::ContribPtr(contrib));
    }
  }

  // minimum distance for all other atom pairs
  fdist = etkdgDetails.boundsMatForceScaling * 10.0;
  for (unsigned int i = 1; i < N; ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      if (!atomPairs[j * N + i]) {
        double l = mmat.getLowerBound(i, j);
        double u = mmat.getUpperBound(i, j);
        auto *contrib = new ForceFields::UFF::DistanceConstraintContrib(
            field, i, j, l, u, fdist);
        field->contribs().push_back(ForceFields::ContribPtr(contrib));
      }
    }
  }

  return field;
}  // construct3DForceField

ForceFields::ForceField *construct3DForceField(
    const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
    const std::map<std::pair<unsigned int, unsigned int>, double> &CPCI) {
  auto *field = construct3DForceField(mmat, positions, etkdgDetails);

  bool is1_4 = false;
  // double dielConst = 1.0;
  boost::uint8_t dielModel = 1;
  for (const auto &charge : CPCI) {
    auto *contrib = new ForceFields::MMFF::EleContrib(
        field, charge.first.first, charge.first.second, charge.second,
        dielModel, is1_4);
    field->contribs().push_back(ForceFields::ContribPtr(contrib));
  }

  return field;
}

ForceFields::ForceField *constructPlain3DForceField(
    const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails) {
  unsigned int N = mmat.numRows();
  CHECK_INVARIANT(N == positions.size(), "");
  CHECK_INVARIANT(etkdgDetails.expTorsionAtoms.size() ==
                      etkdgDetails.expTorsionAngles.size(),
                  "");
  auto *field = new ForceFields::ForceField(positions[0]->dimension());
  for (unsigned int i = 0; i < N; ++i) {
    field->positions().push_back(positions[i]);
  }

  // keep track which atoms are 1,2- or 1,3-restrained
  boost::dynamic_bitset<> atomPairs(N * N);

  // torsion constraints
  for (unsigned int t = 0; t < etkdgDetails.expTorsionAtoms.size(); ++t) {
    int i = etkdgDetails.expTorsionAtoms[t][0];
    int j = etkdgDetails.expTorsionAtoms[t][1];
    int k = etkdgDetails.expTorsionAtoms[t][2];
    int l = etkdgDetails.expTorsionAtoms[t][3];
    if (i < j) {
      atomPairs[i * N + j] = 1;
    } else {
      atomPairs[j * N + i] = 1;
    }
    // etkdgDetails.expTorsionAngles[t][0] = (signs, V's)
    auto *contrib = new ForceFields::CrystalFF::TorsionAngleContribM6(
        field, i, j, k, l, etkdgDetails.expTorsionAngles[t].second,
        etkdgDetails.expTorsionAngles[t].first);
    field->contribs().push_back(ForceFields::ContribPtr(contrib));
  }  // torsion constraints

  double fdist = 100.0;  // force constant
  // 1,2 distance constraints
  for (const auto &bnd : etkdgDetails.bonds) {
    unsigned int i = bnd.first;
    unsigned int j = bnd.second;
    if (i < j) {
      atomPairs[i * N + j] = 1;
    } else {
      atomPairs[j * N + i] = 1;
    }
    double d = ((*positions[i]) - (*positions[j])).length();
    double l = d - 0.01;
    double u = d + 0.01;
    auto *contrib = new ForceFields::UFF::DistanceConstraintContrib(
        field, i, j, l, u, fdist);
    field->contribs().push_back(ForceFields::ContribPtr(contrib));
  }

  // 1,3 distance constraints
  for (unsigned int a = 1; a < etkdgDetails.angles.size(); ++a) {
    unsigned int i = etkdgDetails.angles[a][0];
    unsigned int j = etkdgDetails.angles[a][2];
    if (i < j) {
      atomPairs[i * N + j] = 1;
    } else {
      atomPairs[j * N + i] = 1;
    }
    double d = ((*positions[i]) - (*positions[j])).length();
    double l = d - 0.01;
    double u = d + 0.01;
    auto *contrib = new ForceFields::UFF::DistanceConstraintContrib(
        field, i, j, l, u, fdist);
    field->contribs().push_back(ForceFields::ContribPtr(contrib));
  }

  // minimum distance for all other atom pairs
  fdist = 10.0;
  for (unsigned int i = 1; i < N; ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      if (!atomPairs[j * N + i]) {
        double l = mmat.getLowerBound(i, j);
        double u = mmat.getUpperBound(i, j);
        auto *contrib = new ForceFields::UFF::DistanceConstraintContrib(
            field, i, j, l, u, fdist);
        field->contribs().push_back(ForceFields::ContribPtr(contrib));
      }
    }
  }

  return field;
}  // constructPlain3DForceField

ForceFields::ForceField *construct3DImproperForceField(
    const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    const std::vector<std::vector<int>> &improperAtoms,
    const std::vector<int> &atomNums) {
  (void)atomNums;
  unsigned int N = mmat.numRows();
  CHECK_INVARIANT(N == positions.size(), "");
  auto *field = new ForceFields::ForceField(positions[0]->dimension());
  for (unsigned int i = 0; i < N; ++i) {
    field->positions().push_back(positions[i]);
  }

  // improper torsions / out-of-plane bend / inversion
  double oobForceScalingFactor = 10.0;
  for (const auto &improperAtom : improperAtoms) {
    std::vector<int> n(4);
    for (unsigned int i = 0; i < 3; ++i) {
      n[1] = 1;
      switch (i) {
        case 0:
          n[0] = 0;
          n[2] = 2;
          n[3] = 3;
          break;

        case 1:
          n[0] = 0;
          n[2] = 3;
          n[3] = 2;
          break;

        case 2:
          n[0] = 2;
          n[2] = 3;
          n[3] = 0;
          break;
      }
      auto *contrib = new ForceFields::UFF::InversionContrib(
          field, improperAtom[n[0]], improperAtom[n[1]], improperAtom[n[2]],
          improperAtom[n[3]], improperAtom[4],
          static_cast<bool>(improperAtom[5]), oobForceScalingFactor);
      field->contribs().push_back(ForceFields::ContribPtr(contrib));
    }
  }

  return field;

}  // construct3DImproperForceField
}  // namespace DistGeom
