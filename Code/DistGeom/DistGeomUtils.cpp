//
//  Copyright (C) 2004-2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "BoundsMatrix.h"
#include "DistGeomUtils.h"
#include "DistViolationContribs.h"
#include "ChiralViolationContribs.h"
#include "FourthDimContribs.h"
#include <Numerics/Matrix.h>
#include <Numerics/SymmMatrix.h>
#include <Numerics/Vector.h>
#include <RDGeneral/Invariant.h>
#include <Numerics/EigenSolvers/PowerEigenSolver.h>
#include <RDGeneral/utils.h>
#include <ForceField/ForceField.h>
#include <ForceField/UFF/DistanceConstraints.h>
#include <ForceField/UFF/AngleConstraints.h>
#include <ForceField/UFF/Inversions.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionAngleContribs.h>
#include <boost/dynamic_bitset.hpp>
#include <ForceField/MMFF/Nonbonded.h>

namespace DistGeom {
const double EIGVAL_TOL = 0.001;
const double KNOWN_DIST_TOL = 0.01;

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
  field->positions().insert(field->positions().begin(), positions.begin(),
                            positions.end());

  auto contrib = new DistViolationContribs(field);
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
        contrib->addContrib(i, j, u, l, w);
      }
    }
  }
  if (!contrib->empty()) {
    field->contribs().emplace_back(contrib);
  } else {
    delete contrib;
  }
  // now add chiral constraints
  if (weightChiral > 1.e-8) {
    auto contrib = new ChiralViolationContribs(field);

    for (const auto &cset : csets) {
      contrib->addContrib(cset.get(), weightChiral);
    }
    if (!contrib->empty()) {
      field->contribs().emplace_back(contrib);
    } else {
      delete contrib;
    }
  }

  // finally the contribution from the fourth dimension if we need to
  if ((field->dimension() == 4) && (weightFourthDim > 1.e-8)) {
    auto contrib = new FourthDimContribs(field);
    for (unsigned int i = 0; i < N; i++) {
      contrib->addContrib(i, weightFourthDim);
    }
    if (!contrib->empty()) {
      field->contribs().emplace_back(contrib);
    } else {
      delete contrib;
    }
  }
  return field;
}  // constructForceField

//! Add improper torsion contributions to a force field
/*!

  \param ff                 Force field to add contributions to
  \param forceScalingFactor Force constant to use for contrib
  \param improperAtoms      Indices of atoms to be used in improper torsion
  terms.
  \param is13Constrained    The bit at every position where the atom with
  the same index in the molecule is the center of an improper torsion is set to
  one.

*/
void addImproperTorsionTerms(ForceFields::ForceField *ff,
                             double forceScalingFactor,
                             const std::vector<std::vector<int>> &improperAtoms,
                             boost::dynamic_bitset<> &is13Constrained) {
  auto contrib = new ForceFields::UFF::InversionContribs(ff);
  for (const auto &improperAtom : improperAtoms) {
    std::vector<int> n(4);
    for (unsigned int atomIdx = 0; atomIdx < 3; ++atomIdx) {
      n[1] = 1;
      switch (atomIdx) {
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
      contrib->addContrib(improperAtom[n[0]], improperAtom[n[1]],
                          improperAtom[n[2]], improperAtom[n[3]],
                          improperAtom[4], static_cast<bool>(improperAtom[5]),
                          forceScalingFactor);
      is13Constrained[improperAtom[n[1]]] = 1;
    }
  }
  if (contrib->empty()) {
    delete contrib;
  } else {
    ff->contribs().emplace_back(contrib);
  }
}

//! Add experimental torsion angle contributions to a force field
/*!

  \param ff           Force field to add contributions to
  \param etkdgDetails Contains information about the ETKDG force field
  \param atomPairs    bitset for which atom pairs will be set to one for all
  the 1-4 atom pairs that are part of an experimental torsion contribution
  \param numAtoms     number of atoms in the molecule

 */
void addExperimentalTorsionTerms(
    ForceFields::ForceField *ff,
    const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
    boost::dynamic_bitset<> &atomPairs, unsigned int numAtoms) {
  auto contrib = new ForceFields::CrystalFF::TorsionAngleContribs(ff);
  for (unsigned int torsIdx = 0; torsIdx < etkdgDetails.expTorsionAtoms.size();
       ++torsIdx) {
    int i = etkdgDetails.expTorsionAtoms[torsIdx][0];
    int j = etkdgDetails.expTorsionAtoms[torsIdx][1];
    int k = etkdgDetails.expTorsionAtoms[torsIdx][2];
    int l = etkdgDetails.expTorsionAtoms[torsIdx][3];
    if (i < l) {
      atomPairs[i * numAtoms + l] = 1;
    } else {
      atomPairs[l * numAtoms + i] = 1;
    }
    contrib->addContrib(i, j, k, l,
                        etkdgDetails.expTorsionAngles[torsIdx].second,
                        etkdgDetails.expTorsionAngles[torsIdx].first);
  }
  if (contrib->empty()) {
    delete contrib;
  } else {
    ff->contribs().emplace_back(contrib);
  }
}

//! Add bond constraints with padding at current positions to force field
/*!

  \param ff Force field to add contributions to
  \param etkdgDetails Contains information about the ETKDG force field
  \param atomPairs bitset for which atom pairs will be set to one for all the
  1-4 atom pairs that are part of an experimental torsion contribution
  \param positions       A vector of pointers to 3D Points to write out the
  resulting coordinates
  \param forceConstant force constant with which to constrain bond distances
  \param numAtoms number of atoms in molecule

*/
void add12Terms(ForceFields::ForceField *ff,
                const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
                boost::dynamic_bitset<> &atomPairs,
                RDGeom::Point3DPtrVect &positions, double forceConstant,
                unsigned int numAtoms) {
  auto contrib = new ForceFields::UFF::DistanceConstraintContribs(ff);
  for (const auto &bond : etkdgDetails.bonds) {
    unsigned int i = bond.first;
    unsigned int j = bond.second;
    if (i < j) {
      atomPairs[i * numAtoms + j] = 1;
    } else {
      atomPairs[j * numAtoms + i] = 1;
    }
    double d = ((*positions[i]) - (*positions[j])).length();
    contrib->addContrib(i, j, d - KNOWN_DIST_TOL, d + KNOWN_DIST_TOL,
                        forceConstant);
  }
  if (contrib->empty()) {
    delete contrib;
  } else {
    ff->contribs().emplace_back(contrib);
  }
}
//! Add 1-3 distance constraints with padding at current positions to force
/// field
/*!

  \param ff Force field to add contributions to
  \param etkdgDetails Contains information about the ETKDG force field
  \param atomPairs bitset for which atom pairs will be set to one for all the
  1-4 atom pairs that are part of an experimental torsion contribution
  \param positions       A vector of pointers to 3D Points to write out the
  resulting coordinates
  \param forceConstant force constant with which to constrain bond distances
  \param is13Constrained  The bit at every position where the atom with the
  same index in the molecule is the center of an improper torsion is set to
  one. \param useBasicKnowledge whether to use basic knowledge terms \param
  numAtoms number of atoms in molecule

*/
void add13Terms(ForceFields::ForceField *ff,
                const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
                boost::dynamic_bitset<> &atomPairs,
                RDGeom::Point3DPtrVect &positions, double forceConstant,
                boost::dynamic_bitset<> &is13Constrained,
                bool useBasicKnowledge, unsigned int numAtoms) {
  auto distanceContribs = new ForceFields::UFF::DistanceConstraintContribs(ff);
  auto angleContribs = new ForceFields::UFF::AngleConstraintContribs(ff);
  for (const auto &angle : etkdgDetails.angles) {
    unsigned int i = angle[0];
    unsigned int j = angle[1];
    unsigned int k = angle[2];

    if (i < k) {
      atomPairs[i * numAtoms + k] = 1;
    } else {
      atomPairs[k * numAtoms + i] = 1;
    }
    // check for triple bonds
    if (useBasicKnowledge && angle[3]) {
      angleContribs->addContrib(i, j, k, 179, 180, 1);
    } else if (!is13Constrained[j]) {
      double d = ((*positions[i]) - (*positions[k])).length();
      distanceContribs->addContrib(i, k, d - KNOWN_DIST_TOL, d + KNOWN_DIST_TOL,
                                   forceConstant);
    }
  }
  if (distanceContribs->empty()) {
    delete distanceContribs;
  } else {
    ff->contribs().emplace_back(distanceContribs);
  }
  if (angleContribs->empty()) {
    delete angleContribs;
  } else {
    ff->contribs().emplace_back(angleContribs);
  }
}

//! Add long distance constraints to bounds matrix borders or constrained
//! atoms
/// when provideds
/*!

  \param ff Force field to add contributions to
  \param etkdgDetails Contains information about the ETKDG force field
  \param atomPairs bitset for which atom pairs will be set to one for all the
  1-4 atom pairs that are part of an experimental torsion contribution
  \param positions       A vector of pointers to 3D Points to write out the
  resulting coordinates
  \param knownDistanceForceConstant force constant with which to constrain
  bond distances \param mmat  Bounds matrix to use bounds from for constraints
  \param numAtoms number of atoms in molecule

*/
void addLongRangeDistanceConstraints(
    ForceFields::ForceField *ff,
    const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
    boost::dynamic_bitset<> &atomPairs, RDGeom::Point3DPtrVect &positions,
    double knownDistanceForceConstant, const BoundsMatrix &mmat,
    unsigned int numAtoms) {
  auto contrib = new ForceFields::UFF::DistanceConstraintContribs(ff);
  double fdist = knownDistanceForceConstant;
  for (unsigned int i = 1; i < numAtoms; ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      if (!atomPairs[j * numAtoms + i]) {
        fdist = etkdgDetails.boundsMatForceScaling * 10.0;
        double l = mmat.getLowerBound(i, j);
        double u = mmat.getUpperBound(i, j);
        if (!etkdgDetails.constrainedAtoms.empty() &&
            etkdgDetails.constrainedAtoms[i] &&
            etkdgDetails.constrainedAtoms[j]) {
          // we're constrained, so use very tight bounds
          l = u = ((*positions[i]) - (*positions[j])).length();
          l -= KNOWN_DIST_TOL;
          u += KNOWN_DIST_TOL;
          fdist = knownDistanceForceConstant;
        }
        contrib->addContrib(i, j, l, u, fdist);
      }
    }
  }
  if (contrib->empty()) {
    delete contrib;
  } else {
    ff->contribs().emplace_back(contrib);
  }
}

ForceFields::ForceField *construct3DForceField(
    const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails) {
  unsigned int N = mmat.numRows();
  CHECK_INVARIANT(N == positions.size(), "");
  CHECK_INVARIANT(etkdgDetails.expTorsionAtoms.size() ==
                      etkdgDetails.expTorsionAngles.size(),
                  "");
  auto *field = new ForceFields::ForceField(positions[0]->dimension());
  field->positions().insert(field->positions().begin(), positions.begin(),
                            positions.end());

  // keep track which atoms are 1,2- or 1,3-restrained
  // don't add 1-3 Distances constraints for angles where the
  // central atom of the angle is the central atom of an improper torsion.
  boost::dynamic_bitset<> atomPairs(N * N);
  boost::dynamic_bitset<> is13Constrained(N);

  // ET/K torsion constraints
  addExperimentalTorsionTerms(field, etkdgDetails, atomPairs, N);
  addImproperTorsionTerms(field, 10.0, etkdgDetails.improperAtoms,
                          is13Constrained);

  constexpr double knownDistanceConstraintForce = 100.0;
  add12Terms(field, etkdgDetails, atomPairs, positions,
             knownDistanceConstraintForce, N);
  add13Terms(field, etkdgDetails, atomPairs, positions,
             knownDistanceConstraintForce, is13Constrained, true, N);

  // minimum distance for all other atom pairs that aren't constrained
  addLongRangeDistanceConstraints(field, etkdgDetails, atomPairs, positions,
                                  knownDistanceConstraintForce, mmat, N);
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
    field->contribs().emplace_back(contrib);
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
  field->positions().insert(field->positions().begin(), positions.begin(),
                            positions.end());

  // keep track which atoms are 1,2- or 1,3-restrained
  boost::dynamic_bitset<> atomPairs(N * N);
  boost::dynamic_bitset<> is13Constrained(N);

  // ET torsion constraints
  addExperimentalTorsionTerms(field, etkdgDetails, atomPairs, N);

  constexpr double knownDistanceConstraintForce = 100.0;
  add12Terms(field, etkdgDetails, atomPairs, positions,
             knownDistanceConstraintForce, N);
  add13Terms(field, etkdgDetails, atomPairs, positions,
             knownDistanceConstraintForce, is13Constrained, false, N);
  // minimum distance for all other atom pairs that aren't constrained
  addLongRangeDistanceConstraints(field, etkdgDetails, atomPairs, positions,
                                  knownDistanceConstraintForce, mmat, N);

  return field;
}  // constructPlain3DForceField

ForceFields::ForceField *construct3DImproperForceField(
    const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    const std::vector<std::vector<int>> &improperAtoms,
    const std::vector<std::vector<int>> &angles,
    const std::vector<int> &atomNums) {
  RDUNUSED_PARAM(atomNums);
  unsigned int N = mmat.numRows();
  CHECK_INVARIANT(N == positions.size(), "");
  auto *field = new ForceFields::ForceField(positions[0]->dimension());
  field->positions().insert(field->positions().begin(), positions.begin(),
                            positions.end());

  // improper torsions / out-of-plane bend / inversion
  double oobForceScalingFactor = 10.0;
  boost::dynamic_bitset<> is13Constrained(N);
  addImproperTorsionTerms(field, oobForceScalingFactor, improperAtoms,
                          is13Constrained);

  // Check that SP Centers have an angle of 180 degrees.
  auto contrib = new ForceFields::UFF::AngleConstraintContribs(field);
  for (const auto &angle : angles) {
    if (angle[3]) {
      contrib->addContrib(angle[0], angle[1], angle[2], 179.0, 180.0,
                          oobForceScalingFactor);
    }
  }
  if (contrib->empty()) {
    delete contrib;
  } else {
    field->contribs().emplace_back(contrib);
  }
  return field;
}  // construct3DImproperForceField
}  // namespace DistGeom
