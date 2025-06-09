//
//  Copyright (C) 2004-2025 Greg Landrum and other RDKit contributors
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
#include <ForceField/DistanceConstraints.h>
#include <ForceField/AngleConstraints.h>
#include <ForceField/UFF/Inversions.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionAngleContribs.h>
#include <boost/dynamic_bitset.hpp>
#include <ForceField/MMFF/Nonbonded.h>

namespace DistGeom {
constexpr double EIGVAL_TOL = 0.001;
constexpr double KNOWN_DIST_TOL = 0.01;
constexpr double KNOWN_DIST_FORCE_CONSTANT = 100;

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
    field->contribs().push_back(ForceFields::ContribPtr(contrib));
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
      field->contribs().push_back(ForceFields::ContribPtr(contrib));
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
      field->contribs().push_back(ForceFields::ContribPtr(contrib));
    } else {
      delete contrib;
    }
  }
  return field;
}  // constructForceField

//! Add basic knowledge improper torsion contributions to a force field
/*!

  \param ff Force field to add contributions to
  \param forceScalingFactor Force scaling factor to use in inversion contrib
  \param improperAtoms Indices of atoms to be used in improper torsion
  terms.
  \param isImproperConstrained bit vector with length of total num atoms of
  the molecule where index of every central atom of improper torsion is set to
  one

*/
void addImproperTorsionTerms(ForceFields::ForceField *ff,
                             double forceScalingFactor,
                             const std::vector<std::vector<int>> &improperAtoms,
                             boost::dynamic_bitset<> &isImproperConstrained) {
  PRECONDITION(ff, "bad force field");
  auto inversionContribs =
      std::make_unique<ForceFields::UFF::InversionContribs>(ff);
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

      inversionContribs->addContrib(
          improperAtom[n[0]], improperAtom[n[1]], improperAtom[n[2]],
          improperAtom[n[3]], improperAtom[4],
          static_cast<bool>(improperAtom[5]), forceScalingFactor);
      isImproperConstrained[improperAtom[n[1]]] = 1;
    }
  }
  if (!inversionContribs->empty()) {
    ff->contribs().push_back(std::move(inversionContribs));
  }
}

//! Add experimental torsion angle contributions to a force field
/*!

  \param ff Force field to add contributions to
  \param etkdgDetails Contains information about the ETKDG force field
  \param atomPairs bit set for every atom pair in the molecule where
  a bit is set to one when the atom pair are the end atoms of a torsion
  angle contribution
  \param numAtoms number of atoms in the molecule

 */
void addExperimentalTorsionTerms(
    ForceFields::ForceField *ff,
    const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
    boost::dynamic_bitset<> &atomPairs, unsigned int numAtoms) {
  PRECONDITION(ff, "bad force field");
  auto torsionContribs =
      std::make_unique<ForceFields::CrystalFF::TorsionAngleContribs>(ff);
  for (unsigned int t = 0; t < etkdgDetails.expTorsionAtoms.size(); ++t) {
    int i = etkdgDetails.expTorsionAtoms[t][0];
    int j = etkdgDetails.expTorsionAtoms[t][1];
    int k = etkdgDetails.expTorsionAtoms[t][2];
    int l = etkdgDetails.expTorsionAtoms[t][3];
    if (i < l) {
      atomPairs[i * numAtoms + l] = 1;
    } else {
      atomPairs[l * numAtoms + i] = 1;
    }
    torsionContribs->addContrib(i, j, k, l,
                                etkdgDetails.expTorsionAngles[t].second,
                                etkdgDetails.expTorsionAngles[t].first);
  }
  if (!torsionContribs->empty()) {
    ff->contribs().push_back(std::move(torsionContribs));
  }
}

//! Add bond constraints with padding at current positions to force field
/*!

  \param ff Force field to add contributions to
  \param etkdgDetails Contains information about the ETKDG force field
  \param atomPairs bit set for every atom pair in the molecule where
  a bit is set to one when the atom pair is a bond that is constrained here
  \param positions A vector of pointers to 3D Points to write out the
  resulting coordinates
  \param forceConstant force constant with which to constrain bond distances
  \param numAtoms number of atoms in molecule

*/
void add12Terms(ForceFields::ForceField *ff,
                const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
                boost::dynamic_bitset<> &atomPairs,
                RDGeom::Point3DPtrVect &positions, double forceConstant,
                unsigned int numAtoms) {
  PRECONDITION(ff, "bad force field");
  auto distContribs =
      std::make_unique<ForceFields::DistanceConstraintContribs>(ff);
  for (const auto &bond : etkdgDetails.bonds) {
    unsigned int i = bond.first;
    unsigned int j = bond.second;
    if (i < j) {
      atomPairs[i * numAtoms + j] = 1;
    } else {
      atomPairs[j * numAtoms + i] = 1;
    }
    double d = ((*positions[i]) - (*positions[j])).length();
    distContribs->addContrib(i, j, d - KNOWN_DIST_TOL, d + KNOWN_DIST_TOL,
                             forceConstant);
  }
  if (!distContribs->empty()) {
    ff->contribs().push_back(std::move(distContribs));
  }
}
//! Add 1-3 distance constraints with padding at current positions to force
/// field
/*!

  \param ff Force field to add contributions to
  \param etkdgDetails Contains information about the ETKDG force field
  \param atomPairs bit set for every atom pair in the molecule where
  a bit is set to one when the atom pair is the both end atoms of a 13
  contribution that is constrained here
  \param positions A vector of pointers to 3D Points to write out the resulting
  coordinates \param forceConstant force constant with which to constrain bond
  distances \param isImproperConstrained bit vector with length of total num
  atoms of the molecule where index of every central atom of improper torsion is
  set to one \param useBasicKnowledge whether to use basic knowledge terms
  \param mmat Bounds matrix from which 13 distances are used in case an angle
  is part of an improper torsion
  \param numAtoms number of atoms in molecule

*/
void add13Terms(ForceFields::ForceField *ff,
                const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
                boost::dynamic_bitset<> &atomPairs,
                RDGeom::Point3DPtrVect &positions, double forceConstant,
                const boost::dynamic_bitset<> &isImproperConstrained,
                bool useBasicKnowledge, const BoundsMatrix &mmat,
                unsigned int numAtoms) {
  PRECONDITION(ff, "bad force field");
  auto distContribs =
      std::make_unique<ForceFields::DistanceConstraintContribs>(ff);
  auto angleContribs =
      std::make_unique<ForceFields::AngleConstraintContribs>(ff);
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
      angleContribs->addContrib(i, j, k, 179.0, 180.0, 1);
    } else if (isImproperConstrained[j]) {
      distContribs->addContrib(i, k, mmat.getLowerBound(i, k),
                               mmat.getUpperBound(i, k), forceConstant);
    } else {
      double d = ((*positions[i]) - (*positions[k])).length();
      distContribs->addContrib(i, k, d - KNOWN_DIST_TOL, d + KNOWN_DIST_TOL,
                               forceConstant);
    }
  }
  if (!angleContribs->empty()) {
    ff->contribs().push_back(std::move(angleContribs));
  }
  if (!distContribs->empty()) {
    ff->contribs().push_back(std::move(distContribs));
  }
}

//! Add long distance constraints to bounds matrix borders or constrained atoms
/// when provideds
/*!

  \param ff Force field to add contributions to
  \param etkdgDetails Contains information about the ETKDG force field
  \param atomPairs bit set for every atom pair in the molecule where
  a bit is set to one when the two atoms in the pair are distance constrained
  with respect to each other
  \param positions A vector of pointers to 3D Points to write out the
  resulting coordinates
  \param knownDistanceForceConstant force constant with which to constrain bond
  distances
  \param mmat  Bounds matrix to use bounds from for constraints
  \param numAtoms number of atoms in molecule

*/
void addLongRangeDistanceConstraints(
    ForceFields::ForceField *ff,
    const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
    const boost::dynamic_bitset<> &atomPairs, RDGeom::Point3DPtrVect &positions,
    double knownDistanceForceConstant, const BoundsMatrix &mmat,
    unsigned int numAtoms) {
  PRECONDITION(ff, "bad force field");
  auto distContribs =
      std::make_unique<ForceFields::DistanceConstraintContribs>(ff);
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
        distContribs->addContrib(i, j, l, u, fdist);
      }
    }
  }
  if (!distContribs->empty()) {
    ff->contribs().push_back(std::move(distContribs));
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

  // keep track which atoms are 1,2-, 1,3- or 1,4-restrained
  boost::dynamic_bitset<> atomPairs(N * N);
  // don't add 1-3 Distances constraints for angles where the
  // central atom of the angle is the central atom of an improper torsion.
  boost::dynamic_bitset<> isImproperConstrained(N);

  addExperimentalTorsionTerms(field, etkdgDetails, atomPairs, N);
  addImproperTorsionTerms(field, 10.0, etkdgDetails.improperAtoms,
                          isImproperConstrained);
  add12Terms(field, etkdgDetails, atomPairs, positions,
             KNOWN_DIST_FORCE_CONSTANT, N);
  add13Terms(field, etkdgDetails, atomPairs, positions,
             KNOWN_DIST_FORCE_CONSTANT, isImproperConstrained, true, mmat, N);
  // minimum distance for all other atom pairs that aren't constrained
  addLongRangeDistanceConstraints(field, etkdgDetails, atomPairs, positions,
                                  KNOWN_DIST_FORCE_CONSTANT, mmat, N);
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
  auto *contrib = new ForceFields::MMFF::EleContrib(field);
  field->contribs().emplace_back(contrib);
  for (const auto &charge : CPCI) {
    contrib->addTerm(charge.first.first, charge.first.second, charge.second,
                     dielModel, is1_4);
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

  // keep track which atoms are 1,2-, 1,3- or 1,4-restrained
  boost::dynamic_bitset<> atomPairs(N * N);
  // don't add 1-3 Distances constraints for angles where the
  // central atom of the angle is the central atom of an improper torsion.
  boost::dynamic_bitset<> isImproperConstrained(N);

  addExperimentalTorsionTerms(field, etkdgDetails, atomPairs, N);
  add12Terms(field, etkdgDetails, atomPairs, positions,
             KNOWN_DIST_FORCE_CONSTANT, N);
  add13Terms(field, etkdgDetails, atomPairs, positions,
             KNOWN_DIST_FORCE_CONSTANT, isImproperConstrained, false, mmat, N);
  // minimum distance for all other atom pairs that aren't constrained
  addLongRangeDistanceConstraints(field, etkdgDetails, atomPairs, positions,
                                  KNOWN_DIST_FORCE_CONSTANT, mmat, N);

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
  boost::dynamic_bitset<> isImproperConstrained(N);
  addImproperTorsionTerms(field, oobForceScalingFactor, improperAtoms,
                          isImproperConstrained);

  // Check that SP Centers have an angle of 180 degrees.
  auto angleContribs =
      std::make_unique<ForceFields::AngleConstraintContribs>(field);
  for (const auto &angle : angles) {
    if (angle[3]) {
      angleContribs->addContrib(angle[0], angle[1], angle[2], 179.0, 180.0,
                                oobForceScalingFactor);
    }
  }
  if (!angleContribs->empty()) {
    field->contribs().push_back(std::move(angleContribs));
  }
  return field;
}  // construct3DImproperForceField
}  // namespace DistGeom
