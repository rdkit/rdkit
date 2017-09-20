//
//  Copyright (C) 2004-2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// define DEBUG_EMBEDDING 1
#include "Embedder.h"
#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/DistGeomUtils.h>
#include <DistGeom/TriangleSmooth.h>
#include <DistGeom/ChiralViolationContrib.h>
#include "BoundsMatrixBuilder.h"
#include <ForceField/ForceField.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/RingInfo.h>

#include <GraphMol/Conformer.h>
#include <RDGeneral/types.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Exceptions.h>

#include <Geometry/Transform3D.h>
#include <Numerics/Alignment/AlignPoints.h>
#include <DistGeom/ChiralSet.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>
#include <boost/dynamic_bitset.hpp>
#include <iomanip>
#include <RDGeneral/RDThreads.h>

//#define DEBUG_EMBEDDING 1

#define ERROR_TOL 0.00001
// these tolerances, all to detect and filter out bogus conformations, are a
// delicate balance between sensitive enough to detect obviously bad
// conformations but not so sensitive that a bunch of ok conformations get
// filtered out, which slows down the whole conformation generation process
#define MAX_MINIMIZED_E_PER_ATOM 0.05
#define MAX_MINIMIZED_E_CONTRIB 0.20
#define MIN_TETRAHEDRAL_CHIRAL_VOL 0.50
#define TETRAHEDRAL_CENTERINVOLUME_TOL 0.30

namespace RDKit {
namespace DGeomHelpers {
typedef std::pair<int, int> INT_PAIR;
typedef std::vector<INT_PAIR> INT_PAIR_VECT;

//! Parameters corresponding to Sereina Riniker's KDG approach
const EmbedParameters KDG(0,      // maxIterations
                          1,      // numThreads
                          -1,     // randomSeed
                          true,   // clearConfs
                          false,  // useRandomCoords
                          2.0,    // boxSizeMult
                          true,   // randNegEig
                          1,      // numZeroFail
                          NULL,   // coordMap
                          1e-3,   // optimizerForceTol
                          false,  // ignoreSmoothingFailures
                          true,   // enforceChirality
                          false,  // useExpTorsionAnglePrefs
                          true,   // useBasicKnowledge
                          false,  // verbose
                          5.0,    // basinThresh
                          -1.0,   // pruneRmsThresh
                          true    // onlyHeavyAtomsForRMS
                          );

//! Parameters corresponding to Sereina Riniker's ETDG approach
const EmbedParameters ETDG(0,      // maxIterations
                           1,      // numThreads
                           -1,     // randomSeed
                           true,   // clearConfs
                           false,  // useRandomCoords
                           2.0,    // boxSizeMult
                           true,   // randNegEig
                           1,      // numZeroFail
                           NULL,   // coordMap
                           1e-3,   // optimizerForceTol
                           false,  // ignoreSmoothingFailures
                           false,  // enforceChirality
                           true,   // useExpTorsionAnglePrefs
                           false,  // useBasicKnowledge
                           false,  // verbose
                           5.0,    // basinThresh
                           -1.0,   // pruneRmsThresh
                           true    // onlyHeavyAtomsForRMS
                           );
//! Parameters corresponding to Sereina Riniker's ETKDG approach
const EmbedParameters ETKDG(0,      // maxIterations
                            1,      // numThreads
                            -1,     // randomSeed
                            true,   // clearConfs
                            false,  // useRandomCoords
                            2.0,    // boxSizeMult
                            true,   // randNegEig
                            1,      // numZeroFail
                            NULL,   // coordMap
                            1e-3,   // optimizerForceTol
                            false,  // ignoreSmoothingFailures
                            true,   // enforceChirality
                            true,   // useExpTorsionAnglePrefs
                            true,   // useBasicKnowledge
                            false,  // verbose
                            5.0,    // basinThresh
                            -1.0,   // pruneRmsThresh
                            true    // onlyHeavyAtomsForRMS
                            );

bool _volumeTest(const DistGeom::ChiralSetPtr &chiralSet,
                 const RDGeom::PointPtrVect &positions, bool verbose = false) {
  RDGeom::Point3D p0((*positions[chiralSet->d_idx0])[0],
                     (*positions[chiralSet->d_idx0])[1],
                     (*positions[chiralSet->d_idx0])[2]);
  RDGeom::Point3D p1((*positions[chiralSet->d_idx1])[0],
                     (*positions[chiralSet->d_idx1])[1],
                     (*positions[chiralSet->d_idx1])[2]);
  RDGeom::Point3D p2((*positions[chiralSet->d_idx2])[0],
                     (*positions[chiralSet->d_idx2])[1],
                     (*positions[chiralSet->d_idx2])[2]);
  RDGeom::Point3D p3((*positions[chiralSet->d_idx3])[0],
                     (*positions[chiralSet->d_idx3])[1],
                     (*positions[chiralSet->d_idx3])[2]);
  RDGeom::Point3D p4((*positions[chiralSet->d_idx4])[0],
                     (*positions[chiralSet->d_idx4])[1],
                     (*positions[chiralSet->d_idx4])[2]);

  // even if we are minimizing in higher dimension the chiral volume is
  // calculated using only the first 3 dimensions
  RDGeom::Point3D v1 = p0 - p1;
  v1.normalize();
  RDGeom::Point3D v2 = p0 - p2;
  v2.normalize();
  RDGeom::Point3D v3 = p0 - p3;
  v3.normalize();
  RDGeom::Point3D v4 = p0 - p4;
  v4.normalize();

  RDGeom::Point3D crossp = v1.crossProduct(v2);
  double vol = crossp.dotProduct(v3);
  if (verbose) std::cerr << "   " << fabs(vol) << std::endl;
  if (fabs(vol) < MIN_TETRAHEDRAL_CHIRAL_VOL) return false;
  crossp = v1.crossProduct(v2);
  vol = crossp.dotProduct(v4);
  if (verbose) std::cerr << "   " << fabs(vol) << std::endl;
  if (fabs(vol) < MIN_TETRAHEDRAL_CHIRAL_VOL) return false;
  crossp = v1.crossProduct(v3);
  vol = crossp.dotProduct(v4);
  if (verbose) std::cerr << "   " << fabs(vol) << std::endl;
  if (fabs(vol) < MIN_TETRAHEDRAL_CHIRAL_VOL) return false;
  crossp = v2.crossProduct(v3);
  vol = crossp.dotProduct(v4);
  if (verbose) std::cerr << "   " << fabs(vol) << std::endl;
  if (fabs(vol) < MIN_TETRAHEDRAL_CHIRAL_VOL) return false;

  return true;
}

bool _sameSide(const RDGeom::Point3D &v1, const RDGeom::Point3D &v2,
               const RDGeom::Point3D &v3, const RDGeom::Point3D &v4,
               const RDGeom::Point3D &p0, double tol = 0.1) {
  RDGeom::Point3D normal = (v2 - v1).crossProduct(v3 - v1);
  double d1 = normal.dotProduct(v4 - v1);
  double d2 = normal.dotProduct(p0 - v1);
  // std::cerr << "     " << d1 << " - " << d2 << std::endl;
  if (fabs(d1) < tol || fabs(d2) < tol) return false;
  return !((d1 < 0.) ^ (d2 < 0.));
}
bool _centerInVolume(unsigned int idx0, unsigned int idx1, unsigned int idx2,
                     unsigned int idx3, unsigned int idx4,
                     const RDGeom::PointPtrVect &positions, double tol,
                     bool verbose = false) {
  RDGeom::Point3D p0((*positions[idx0])[0], (*positions[idx0])[1],
                     (*positions[idx0])[2]);
  RDGeom::Point3D p1((*positions[idx1])[0], (*positions[idx1])[1],
                     (*positions[idx1])[2]);
  RDGeom::Point3D p2((*positions[idx2])[0], (*positions[idx2])[1],
                     (*positions[idx2])[2]);
  RDGeom::Point3D p3((*positions[idx3])[0], (*positions[idx3])[1],
                     (*positions[idx3])[2]);
  RDGeom::Point3D p4((*positions[idx4])[0], (*positions[idx4])[1],
                     (*positions[idx4])[2]);
  // RDGeom::Point3D centroid = (p1+p2+p3+p4)/4.;
  if (verbose) {
    std::cerr << _sameSide(p1, p2, p3, p4, p0, tol) << " "
              << _sameSide(p2, p3, p4, p1, p0, tol) << " "
              << _sameSide(p3, p4, p1, p2, p0, tol) << " "
              << _sameSide(p4, p1, p2, p3, p0, tol) << std::endl;
  }
  bool res = _sameSide(p1, p2, p3, p4, p0, tol) &&
             _sameSide(p2, p3, p4, p1, p0, tol) &&
             _sameSide(p3, p4, p1, p2, p0, tol) &&
             _sameSide(p4, p1, p2, p3, p0, tol);
  return res;
}

bool _centerInVolume(const DistGeom::ChiralSetPtr &chiralSet,
                     const RDGeom::PointPtrVect &positions, double tol = 0.1,
                     bool verbose = false) {
  if (chiralSet->d_idx0 ==
      chiralSet->d_idx4) {  // this happens for three-coordinate centers
    return true;
  }
  return _centerInVolume(chiralSet->d_idx0, chiralSet->d_idx1,
                         chiralSet->d_idx2, chiralSet->d_idx3,
                         chiralSet->d_idx4, positions, tol, verbose);
}
bool _boundsFulfilled(const std::vector<int> &atoms,
                      const DistGeom::BoundsMatrix &mmat,
                      const RDGeom::PointPtrVect &positions) {
  // unsigned int N = mmat.numRows();
  // std::cerr << N << " " << atoms.size() << std::endl;
  // loop over all pair of atoms
  for (unsigned int i = 0; i < atoms.size() - 1; ++i) {
    for (unsigned int j = i + 1; j < atoms.size(); ++j) {
      int a1 = atoms[i];
      int a2 = atoms[j];
      RDGeom::Point3D p0((*positions[a1])[0], (*positions[a1])[1],
                         (*positions[a1])[2]);
      RDGeom::Point3D p1((*positions[a2])[0], (*positions[a2])[1],
                         (*positions[a2])[2]);
      double d2 = (p0 - p1).length();  // distance
      double lb = mmat.getLowerBound(a1, a2);
      double ub = mmat.getUpperBound(a1, a2);  // bounds
      if (((d2 < lb) && (fabs(d2 - lb) > 0.1 * ub)) ||
          ((d2 > ub) && (fabs(d2 - ub) > 0.1 * ub))) {
#ifdef DEBUG_EMBEDDING
        std::cerr << a1 << " " << a2 << ":" << d2 << " " << lb << " " << ub
                  << " " << fabs(d2 - lb) << " " << fabs(d2 - ub) << std::endl;
#endif
        return false;
      }
    }
  }
  return true;
}

// the minimization using experimental torsion angle preferences
bool _minimizeWithExpTorsions(
    RDGeom::PointPtrVect &positions, DistGeom::BoundsMatPtr mmat,
    double optimizerForceTol, double basinThresh,
    const std::vector<std::pair<int, int> > &bonds,
    const std::vector<std::vector<int> > &angles,
    const std::vector<std::vector<int> > &expTorsionAtoms,
    const std::vector<std::pair<std::vector<int>, std::vector<double> > >
        &expTorsionAngles,
    const std::vector<std::vector<int> > &improperAtoms,
    const std::vector<int> &atomNums, bool useBasicKnowledge) {
  RDUNUSED_PARAM(basinThresh);

  bool planar = true;

  // convert to 3D positions and create coordMap
  RDGeom::Point3DPtrVect positions3D;
  for (unsigned int p = 0; p < positions.size(); ++p) {
    positions3D.push_back(new RDGeom::Point3D(
        (*positions[p])[0], (*positions[p])[1], (*positions[p])[2]));
  }

  // create the force field
  ForceFields::ForceField *field;
  if (useBasicKnowledge) {  // ETKDG or KDG
    field = DistGeom::construct3DForceField(*mmat, positions3D, bonds, angles,
                                            expTorsionAtoms, expTorsionAngles,
                                            improperAtoms, atomNums);
  } else {  // plain ETDG
    field = DistGeom::constructPlain3DForceField(*mmat, positions3D, bonds,
                                                 angles, expTorsionAtoms,
                                                 expTorsionAngles, atomNums);
  }

  // minimize!
  field->initialize();
  // std::cout << "Field with torsion constraints: " << field->calcEnergy() << "
  // " << ERROR_TOL << std::endl;
  if (field->calcEnergy() > ERROR_TOL) {
    // while (needMore) {
    field->minimize(300, optimizerForceTol);
    //      ++nPasses;
    //}
  }
  // std::cout << field->calcEnergy() << std::endl;
  delete field;

  // check for planarity if ETKDG or KDG
  if (useBasicKnowledge) {
    // create a force field with only the impropers
    ForceFields::ForceField *field2;
    field2 = DistGeom::construct3DImproperForceField(*mmat, positions3D,
                                                     improperAtoms, atomNums);
    field2->initialize();
    // check if the energy is low enough
    double planarityTolerance = 0.7;
    if (field2->calcEnergy() > improperAtoms.size() * planarityTolerance) {
#ifdef DEBUG_EMBEDDING
      std::cerr << "   planar fail: " << field2->calcEnergy() << " "
                << improperAtoms.size() * planarityTolerance << std::endl;
#endif
      planar = false;
    }
    delete field2;
  }

  // overwrite positions and delete the 3D ones
  for (unsigned int i = 0; i < positions3D.size(); ++i) {
    (*positions[i])[0] = (*positions3D[i])[0];
    (*positions[i])[1] = (*positions3D[i])[1];
    (*positions[i])[2] = (*positions3D[i])[2];
    delete positions3D[i];
  }

  return planar;
}

bool _embedPoints(
    RDGeom::PointPtrVect *positions, const DistGeom::BoundsMatPtr mmat,
    bool useRandomCoords, double boxSizeMult, bool randNegEig,
    unsigned int numZeroFail, double optimizerForceTol, double basinThresh,
    int seed, unsigned int maxIterations,
    const DistGeom::VECT_CHIRALSET *chiralCenters,
    const DistGeom::VECT_CHIRALSET *tetrahedralCarbons, bool enforceChirality,
    bool useExpTorsionAnglePrefs, bool useBasicKnowledge,
    const std::vector<std::pair<int, int> > &bonds,
    const std::vector<std::vector<int> > &angles,
    const std::vector<std::vector<int> > &expTorsionAtoms,
    const std::vector<std::pair<std::vector<int>, std::vector<double> > >
        &expTorsionAngles,
    const std::vector<std::vector<int> > &improperAtoms,
    const std::vector<int> &atomNums) {
  unsigned int nat = positions->size();
  if (maxIterations == 0) {
    maxIterations = 10 * nat;
  }
  RDNumeric::DoubleSymmMatrix distMat(nat, 0.0);

  // The basin threshold just gets us into trouble when we're using
  // random coordinates since it ends up ignoring 1-4 (and higher)
  // interactions. This causes us to get folded-up (and self-penetrating)
  // conformations for large flexible molecules
  if (useRandomCoords) basinThresh = 1e8;

  RDKit::double_source_type *rng = 0;
  RDKit::rng_type *generator;
  RDKit::uniform_double *distrib;
  if (seed > 0) {
    generator = new RDKit::rng_type(42u);
    generator->seed(seed);
    distrib = new RDKit::uniform_double(0.0, 1.0);
    rng = new RDKit::double_source_type(*generator, *distrib);
  } else {
    rng = &RDKit::getDoubleRandomSource();
  }

  bool gotCoords = false;
  unsigned int iter = 0;
  double largestDistance = -1.0;
  RDUNUSED_PARAM(largestDistance);
  while ((gotCoords == false) && (iter < maxIterations)) {
    ++iter;
    if (!useRandomCoords) {
      largestDistance = DistGeom::pickRandomDistMat(*mmat, distMat, *rng);
      gotCoords = DistGeom::computeInitialCoords(distMat, *positions, *rng,
                                                 randNegEig, numZeroFail);
    } else {
      double boxSize;
      if (boxSizeMult > 0) {
        boxSize = 5. * boxSizeMult;
      } else {
        boxSize = -1 * boxSizeMult;
      }
      gotCoords = DistGeom::computeRandomCoords(*positions, boxSize, *rng);
    }
#ifdef DEBUG_EMBEDDING
    if (!gotCoords) {
      std::cerr << "Initial embedding failed!, Iter: " << iter << std::endl;
    }
#endif
    // std::cerr << " ITER: " << iter << " gotCoords: " << gotCoords <<
    // std::endl;
    if (gotCoords) {
      boost::scoped_ptr<ForceFields::ForceField> field(
          DistGeom::constructForceField(*mmat, *positions, *chiralCenters, 1.0,
                                        0.1, 0, basinThresh));
      unsigned int nPasses = 0;
      field->initialize();
      // std::cerr << "FIELD E: " << field->calcEnergy() << std::endl;
      if (field->calcEnergy() > ERROR_TOL) {
        int needMore = 1;
        while (needMore) {
          needMore = field->minimize(400, optimizerForceTol);
          ++nPasses;
        }
      }
      std::vector<double> e_contribs;
      double local_e = field->calcEnergy(&e_contribs);
// if (e_contribs.size()) {
//   std::cerr << "        check: " << local_e / nat << " "
//             << *(std::max_element(e_contribs.begin(),
//             e_contribs.end()))
//             << std::endl;
// }

#ifdef DEBUG_EMBEDDING
      std::cerr << " Energy : " << local_e / nat << " "
                << *(std::max_element(e_contribs.begin(), e_contribs.end()))
                << std::endl;
// std::copy(e_contribs.begin(), e_contribs.end(),
//           std::ostream_iterator<double>(std::cerr, " "));
// std::cerr << std::endl;
#endif

      // check that neither the energy nor any of the contributions to it are
      // too high (this is part of github #971)
      if (local_e / nat >= MAX_MINIMIZED_E_PER_ATOM ||
          (e_contribs.size() &&
           *(std::max_element(e_contribs.begin(), e_contribs.end())) >
               MAX_MINIMIZED_E_CONTRIB)) {
#ifdef DEBUG_EMBEDDING
        std::cerr << " Energy fail: " << local_e / nat << " "
                  << *(std::max_element(e_contribs.begin(), e_contribs.end()))
                  << std::endl;
#endif
        gotCoords = false;
        continue;
      }
      // for each of the atoms in the "tetrahedralCarbons" list, make sure
      // that there is a minimum volume around them and that they are inside
      // that volume. (this is part of github #971)
      BOOST_FOREACH (DistGeom::ChiralSetPtr tetSet, *tetrahedralCarbons) {
        // it could happen that the centroid is outside the volume defined
        // by the other
        // four points. That is also a fail.
        if (!_volumeTest(tetSet, *positions) ||
            !_centerInVolume(tetSet, *positions,
                             TETRAHEDRAL_CENTERINVOLUME_TOL)) {
#ifdef DEBUG_EMBEDDING
          std::cerr << " fail2! (" << tetSet->d_idx0 << ") iter: " << iter
                    << " vol: " << _volumeTest(tetSet, *positions, true)
                    << " center: "
                    << _centerInVolume(tetSet, *positions,
                                       TETRAHEDRAL_CENTERINVOLUME_TOL, true)
                    << std::endl;
#endif
          gotCoords = false;
          continue;
        }
      }

      // Check if any of our chiral centers are badly out of whack. If so, try
      // again
      if (gotCoords && enforceChirality && chiralCenters->size() > 0) {
        // check the chiral volume:
        BOOST_FOREACH (DistGeom::ChiralSetPtr chiralSet, *chiralCenters) {
          double vol = DistGeom::ChiralViolationContrib::calcChiralVolume(
              chiralSet->d_idx1, chiralSet->d_idx2, chiralSet->d_idx3,
              chiralSet->d_idx4, *positions);
          double lb = chiralSet->getLowerVolumeBound();
          double ub = chiralSet->getUpperVolumeBound();
          if ((lb > 0 && vol < lb && (lb - vol) / lb > .2) ||
              (ub < 0 && vol > ub && (vol - ub) / ub > .2)) {
#ifdef DEBUG_EMBEDDING
            std::cerr << " fail! (" << chiralSet->d_idx0 << ") iter: " << iter
                      << " " << vol << " " << lb << "-" << ub << std::endl;
#endif
            gotCoords = false;
            break;
          }
        }
      }
      // now redo the minimization if we have a chiral center
      // or have started from random coords. This
      // time removing the chiral constraints and
      // increasing the weight on the fourth dimension
      if (gotCoords && (chiralCenters->size() > 0 || useRandomCoords)) {
        boost::scoped_ptr<ForceFields::ForceField> field2(
            DistGeom::constructForceField(*mmat, *positions, *chiralCenters,
                                          0.2, 1.0, 0, basinThresh));
        field2->initialize();
        // std::cerr<<"FIELD2 E: "<<field2->calcEnergy()<<std::endl;
        if (field2->calcEnergy() > ERROR_TOL) {
          int needMore = 1;
          int nPasses2 = 0;
          while (needMore) {
            needMore = field2->minimize(200, optimizerForceTol);
            ++nPasses2;
          }
          // std::cerr<<"   "<<field2->calcEnergy()<<" after npasses2:
          // "<<nPasses2<<std::endl;
        }
      }

      // (ET)(K)DG
      if (gotCoords && (useExpTorsionAnglePrefs || useBasicKnowledge)) {
        gotCoords = _minimizeWithExpTorsions(
            *positions, mmat, optimizerForceTol, basinThresh, bonds, angles,
            expTorsionAtoms, expTorsionAngles, improperAtoms, atomNums,
            useBasicKnowledge);
      }
      // test if chirality is correct
      if (enforceChirality && gotCoords && (chiralCenters->size() > 0)) {
        // "distance matrix" chirality test
        std::set<int> atoms;
        BOOST_FOREACH (DistGeom::ChiralSetPtr chiralSet, *chiralCenters) {
          if (chiralSet->d_idx0 != chiralSet->d_idx4) {
            atoms.insert(chiralSet->d_idx0);
            atoms.insert(chiralSet->d_idx1);
            atoms.insert(chiralSet->d_idx2);
            atoms.insert(chiralSet->d_idx3);
            atoms.insert(chiralSet->d_idx4);
          }
        }
        std::vector<int> atomsToCheck(atoms.begin(), atoms.end());
        if (atomsToCheck.size() > 0) {
          if (!_boundsFulfilled(atomsToCheck, *mmat, *positions)) {
            gotCoords = false;

#ifdef DEBUG_EMBEDDING
            std::cerr << " fail3a! (" << atomsToCheck[0] << ") iter: " << iter
                      << std::endl;
#endif
          }
        }

        // "center in volume" chirality test
        if (gotCoords) {
          BOOST_FOREACH (DistGeom::ChiralSetPtr chiralSet, *chiralCenters) {
            // it could happen that the centroid is outside the volume defined
            // by the other
            // four points. That is also a fail.
            if (!_centerInVolume(chiralSet, *positions)) {
#ifdef DEBUG_EMBEDDING
              std::cerr << " fail3b! (" << chiralSet->d_idx0
                        << ") iter: " << iter << std::endl;
#endif
              gotCoords = false;
              break;
            }
          }
        }
      }
    }  // if(gotCoords)
  }    // while
  if (seed > 0 && rng) {
    delete rng;
    delete generator;
    delete distrib;
  }
  return gotCoords;
}

void _findChiralSets(const ROMol &mol, DistGeom::VECT_CHIRALSET &chiralCenters,
                     DistGeom::VECT_CHIRALSET &tetrahedralCenters,
                     const std::map<int, RDGeom::Point3D> *coordMap) {
  ROMol::ConstAtomIterator ati;
  INT_VECT nbrs;
  ROMol::OEDGE_ITER beg, end;
  // Atom *oatom;
  for (ati = mol.beginAtoms(); ati != mol.endAtoms(); ati++) {
    if ((*ati)->getAtomicNum() != 1) {  // skip hydrogens
      Atom::ChiralType chiralType = (*ati)->getChiralTag();
      if ((chiralType == Atom::CHI_TETRAHEDRAL_CW ||
           chiralType == Atom::CHI_TETRAHEDRAL_CCW) ||
          (((*ati)->getAtomicNum() == 6 || (*ati)->getAtomicNum() == 7) &&
           (*ati)->getDegree() == 4)) {
        // make a chiral set from the neighbors
        nbrs.clear();
        nbrs.reserve(4);
        // find the neighbors of this atom and enter them into the
        // nbr list
        boost::tie(beg, end) = mol.getAtomBonds(*ati);
        while (beg != end) {
          nbrs.push_back(mol[*beg]->getOtherAtom(*ati)->getIdx());
          ++beg;
        }
        // if we have less than 4 heavy atoms as neighbors,
        // we need to include the chiral center into the mix
        // we should at least have 3 though
        bool includeSelf = false;
        RDUNUSED_PARAM(includeSelf);
        CHECK_INVARIANT(nbrs.size() >= 3, "Cannot be a chiral center");

        if (nbrs.size() < 4) {
          nbrs.insert(nbrs.end(), (*ati)->getIdx());
          includeSelf = true;
        }

        // now create a chiral set and set the upper and lower bound on the
        // volume
        if (chiralType == Atom::CHI_TETRAHEDRAL_CCW) {
          // postive chiral volume
          DistGeom::ChiralSet *cset = new DistGeom::ChiralSet(
              (*ati)->getIdx(), nbrs[0], nbrs[1], nbrs[2], nbrs[3], 5.0, 100.0);
          DistGeom::ChiralSetPtr cptr(cset);
          chiralCenters.push_back(cptr);
        } else if (chiralType == Atom::CHI_TETRAHEDRAL_CW) {
          DistGeom::ChiralSet *cset =
              new DistGeom::ChiralSet((*ati)->getIdx(), nbrs[0], nbrs[1],
                                      nbrs[2], nbrs[3], -100.0, -5.0);
          DistGeom::ChiralSetPtr cptr(cset);
          chiralCenters.push_back(cptr);
        } else {
          if ((coordMap &&
               coordMap->find((*ati)->getIdx()) != coordMap->end()) ||
              (mol.getRingInfo()->isInitialized() &&
               (mol.getRingInfo()->numAtomRings((*ati)->getIdx()) < 2 ||
                mol.getRingInfo()->isAtomInRingOfSize((*ati)->getIdx(), 3)))) {
            // we only want to these tests for ring atoms that are not part of
            // the coordMap
            // there's no sense doing 3-rings because those are a nightmare
          } else {
            DistGeom::ChiralSet *cset = new DistGeom::ChiralSet(
                (*ati)->getIdx(), nbrs[0], nbrs[1], nbrs[2], nbrs[3], 0.0, 0.0);
            DistGeom::ChiralSetPtr cptr(cset);
            tetrahedralCenters.push_back(cptr);
          }
        }
      }  // if block -chirality check
    }    // if block - heavy atom check
  }      // for loop over atoms
}  // end of _findChiralSets

void _fillAtomPositions(RDGeom::Point3DConstPtrVect &pts, const Conformer &conf,
                        const ROMol &mol, bool onlyHeavyAtomsForRMS) {
  unsigned int na = conf.getNumAtoms();
  pts.clear();
  unsigned int ai;
  pts.reserve(na);
  for (ai = 0; ai < na; ++ai) {
    // FIX: should we include D and T here?
    if (onlyHeavyAtomsForRMS && mol.getAtomWithIdx(ai)->getAtomicNum() == 1) {
      continue;
    }
    pts.push_back(&conf.getAtomPos(ai));
  }
}

bool _isConfFarFromRest(const ROMol &mol, const Conformer &conf,
                        double threshold, bool onlyHeavyAtomsForRMS) {
  // NOTE: it is tempting to use some triangle inequality to prune
  // conformations here but some basic testing has shown very
  // little advantage and given that the time for pruning fades in
  // comparison to embedding - we will use a simple for loop below
  // over all conformation until we find a match
  ROMol::ConstConformerIterator confi;

  RDGeom::Point3DConstPtrVect refPoints, prbPoints;
  _fillAtomPositions(refPoints, conf, mol, onlyHeavyAtomsForRMS);

  bool res = true;
  unsigned int na = conf.getNumAtoms();
  double ssrThres = na * threshold * threshold;

  RDGeom::Transform3D trans;
  double ssr;
  for (confi = mol.beginConformers(); confi != mol.endConformers(); confi++) {
    _fillAtomPositions(prbPoints, *(*confi), mol, onlyHeavyAtomsForRMS);
    ssr = RDNumeric::Alignments::AlignPoints(refPoints, prbPoints, trans);
    if (ssr < ssrThres) {
      res = false;
      break;
    }
  }
  return res;
}

int EmbedMolecule(ROMol &mol, unsigned int maxIterations, int seed,
                  bool clearConfs, bool useRandomCoords, double boxSizeMult,
                  bool randNegEig, unsigned int numZeroFail,
                  const std::map<int, RDGeom::Point3D> *coordMap,
                  double optimizerForceTol, bool ignoreSmoothingFailures,
                  bool enforceChirality, bool useExpTorsionAnglePrefs,
                  bool useBasicKnowledge, bool verbose, double basinThresh,
                  bool onlyHeavyAtomsForRMS) {
  INT_VECT confIds;
  EmbedMultipleConfs(
      mol, confIds, 1, 1, maxIterations, seed, clearConfs, useRandomCoords,
      boxSizeMult, randNegEig, numZeroFail, -1.0, coordMap, optimizerForceTol,
      ignoreSmoothingFailures, enforceChirality, useExpTorsionAnglePrefs,
      useBasicKnowledge, verbose, basinThresh, onlyHeavyAtomsForRMS);

  int res;
  if (confIds.size()) {
    res = confIds[0];
  } else {
    res = -1;
  }
  return res;
}

void adjustBoundsMatFromCoordMap(
    DistGeom::BoundsMatPtr mmat, unsigned int nAtoms,
    const std::map<int, RDGeom::Point3D> *coordMap) {
  RDUNUSED_PARAM(nAtoms);
  // std::cerr<<std::endl;
  // for(unsigned int i=0;i<nAtoms;++i){
  //   for(unsigned int j=0;j<nAtoms;++j){
  //     std::cerr<<"  "<<std::setprecision(3)<<mmat->getVal(i,j);
  //   }
  //   std::cerr<<std::endl;
  // }
  // std::cerr<<std::endl;
  for (std::map<int, RDGeom::Point3D>::const_iterator iIt = coordMap->begin();
       iIt != coordMap->end(); ++iIt) {
    int iIdx = iIt->first;
    const RDGeom::Point3D &iPoint = iIt->second;

    std::map<int, RDGeom::Point3D>::const_iterator jIt = iIt;
    while (++jIt != coordMap->end()) {
      int jIdx = jIt->first;
      const RDGeom::Point3D &jPoint = jIt->second;
      double dist = (iPoint - jPoint).length();
      mmat->setUpperBound(iIdx, jIdx, dist);
      mmat->setLowerBound(iIdx, jIdx, dist);
    }
  }
  // std::cerr<<std::endl;
  // for(unsigned int i=0;i<nAtoms;++i){
  //   for(unsigned int j=0;j<nAtoms;++j){
  //     std::cerr<<"  "<<std::setprecision(3)<<mmat->getVal(i,j);
  //   }
  //   std::cerr<<std::endl;
  // }
  // std::cerr<<std::endl;
}

namespace detail {
typedef struct {
  boost::dynamic_bitset<> *confsOk;
  bool fourD;
  INT_VECT *fragMapping;
  std::vector<Conformer *> *confs;
  unsigned int fragIdx;
  DistGeom::BoundsMatPtr mmat;
  bool useRandomCoords;
  double boxSizeMult;
  bool randNegEig;
  unsigned int numZeroFail;
  double optimizerForceTol;
  double basinThresh;
  int seed;
  unsigned int maxIterations;
  DistGeom::VECT_CHIRALSET const *chiralCenters;
  DistGeom::VECT_CHIRALSET const *tetrahedralCarbons;
  bool enforceChirality;
  bool useExpTorsionAnglePrefs;
  bool useBasicKnowledge;
  std::vector<std::pair<int, int> > *bonds;
  std::vector<std::vector<int> > *angles;
  std::vector<std::vector<int> > *expTorsionAtoms;
  std::vector<std::pair<std::vector<int>, std::vector<double> > >
      *expTorsionAngles;
  std::vector<std::vector<int> > *improperAtoms;
  std::vector<int> *atomNums;
} EmbedArgs;
void embedHelper_(int threadId, int numThreads, EmbedArgs *eargs) {
  unsigned int nAtoms = eargs->mmat->numRows();
  RDGeom::PointPtrVect positions;
  for (unsigned int i = 0; i < nAtoms; ++i) {
    if (eargs->fourD) {
      positions.push_back(new RDGeom::PointND(4));
    } else {
      positions.push_back(new RDGeom::Point3D());
    }
  }
  for (size_t ci = 0; ci < eargs->confs->size(); ci++) {
    if (rdcast<int>(ci % numThreads) != threadId) continue;
    if (!(*eargs->confsOk)[ci]) {
      // if one of the fragments here has already failed, there's no
      // sense in embedding this one
      continue;
    }
    bool gotCoords = _embedPoints(
        &positions, eargs->mmat, eargs->useRandomCoords, eargs->boxSizeMult,
        eargs->randNegEig, eargs->numZeroFail, eargs->optimizerForceTol,
        eargs->basinThresh, (ci + 1) * eargs->seed, eargs->maxIterations,
        eargs->chiralCenters, eargs->tetrahedralCarbons,
        eargs->enforceChirality, eargs->useExpTorsionAnglePrefs,
        eargs->useBasicKnowledge, *eargs->bonds, *eargs->angles,
        *eargs->expTorsionAtoms, *eargs->expTorsionAngles,
        *eargs->improperAtoms, *eargs->atomNums);

    if (gotCoords) {
      Conformer *conf = (*eargs->confs)[ci];
      unsigned int fragAtomIdx = 0;
      for (unsigned int i = 0; i < (*eargs->confs)[0]->getNumAtoms(); ++i) {
        if ((*eargs->fragMapping)[i] == static_cast<int>(eargs->fragIdx)) {
          conf->setAtomPos(i, RDGeom::Point3D((*positions[fragAtomIdx])[0],
                                              (*positions[fragAtomIdx])[1],
                                              (*positions[fragAtomIdx])[2]));
          ++fragAtomIdx;
        }
      }
    } else {
      (*eargs->confsOk)[ci] = 0;
    }
  }
  for (unsigned int i = 0; i < nAtoms; ++i) {
    delete positions[i];
  }
}
}  // end of namespace detail

void EmbedMultipleConfs(ROMol &mol, INT_VECT &res, unsigned int numConfs,
                        int numThreads, unsigned int maxIterations, int seed,
                        bool clearConfs, bool useRandomCoords,
                        double boxSizeMult, bool randNegEig,
                        unsigned int numZeroFail, double pruneRmsThresh,
                        const std::map<int, RDGeom::Point3D> *coordMap,
                        double optimizerForceTol, bool ignoreSmoothingFailures,
                        bool enforceChirality, bool useExpTorsionAnglePrefs,
                        bool useBasicKnowledge, bool verbose,
                        double basinThresh, bool onlyHeavyAtomsForRMS) {
  if (!mol.getNumAtoms()) {
    throw ValueErrorException("molecule has no atoms");
  }

  INT_VECT fragMapping;
  std::vector<ROMOL_SPTR> molFrags =
      MolOps::getMolFrags(mol, true, &fragMapping);
  if (molFrags.size() > 1 && coordMap) {
    BOOST_LOG(rdWarningLog)
        << "Constrained conformer generation (via the coordMap argument) does "
           "not work with molecules that have multiple fragments."
        << std::endl;
    coordMap = 0;
  }
  std::vector<Conformer *> confs;
  confs.reserve(numConfs);
  for (unsigned int i = 0; i < numConfs; ++i) {
    confs.push_back(new Conformer(mol.getNumAtoms()));
  }
  boost::dynamic_bitset<> confsOk(numConfs);
  confsOk.set();

  if (clearConfs) {
    res.clear();
    mol.clearConformers();
  }

  for (unsigned int fragIdx = 0; fragIdx < molFrags.size(); ++fragIdx) {
    ROMOL_SPTR piece = molFrags[fragIdx];
    unsigned int nAtoms = piece->getNumAtoms();
    DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(nAtoms);
    DistGeom::BoundsMatPtr mmat(mat);
    initBoundsMat(mmat);

    double tol = 0.0;
    std::vector<std::vector<int> > expTorsionAtoms;
    std::vector<std::pair<std::vector<int>, std::vector<double> > >
        expTorsionAngles;
    std::vector<std::vector<int> > improperAtoms;
    std::vector<std::pair<int, int> > bonds;
    std::vector<std::vector<int> > angles;
    std::vector<int> atomNums(nAtoms);
    if (useExpTorsionAnglePrefs || useBasicKnowledge) {
      ForceFields::CrystalFF::getExperimentalTorsions(
          *piece, expTorsionAtoms, expTorsionAngles, improperAtoms,
          useExpTorsionAnglePrefs, useBasicKnowledge, verbose);
      setTopolBounds(*piece, mmat, bonds, angles, true, false);
      for (unsigned int i = 0; i < nAtoms; ++i) {
        atomNums[i] = (*piece).getAtomWithIdx(i)->getAtomicNum();
      }
    } else {
      setTopolBounds(*piece, mmat, true, false);
    }
    if (coordMap) {
      adjustBoundsMatFromCoordMap(mmat, nAtoms, coordMap);
      tol = 0.05;
    }
    if (!DistGeom::triangleSmoothBounds(mmat, tol)) {
      // ok this bound matrix failed to triangle smooth - re-compute the bounds
      // matrix
      // without 15 bounds and with VDW scaling
      initBoundsMat(mmat);
      setTopolBounds(*piece, mmat, false, true);

      if (coordMap) {
        adjustBoundsMatFromCoordMap(mmat, nAtoms, coordMap);
      }

      // try triangle smoothing again
      if (!DistGeom::triangleSmoothBounds(mmat, tol)) {
        // ok, we're not going to be able to smooth this,
        if (ignoreSmoothingFailures) {
          // proceed anyway with the more relaxed bounds matrix
          initBoundsMat(mmat);
          setTopolBounds(*piece, mmat, false, true);

          if (coordMap) {
            adjustBoundsMatFromCoordMap(mmat, nAtoms, coordMap);
          }
        } else {
          BOOST_LOG(rdWarningLog)
              << "Could not triangle bounds smooth molecule." << std::endl;
          return;
        }
      }
    }
#if 0
        for(unsigned int li=0;li<piece->getNumAtoms();++li){
          for(unsigned int lj=li+1;lj<piece->getNumAtoms();++lj){
            std::cerr<<" ("<<li<<","<<lj<<"): "<<mat->getLowerBound(li,lj)<<" -> "<<mat->getUpperBound(li,lj)<<std::endl;
          }
        }
#endif
    // find all the chiral centers in the molecule
    DistGeom::VECT_CHIRALSET chiralCenters;
    DistGeom::VECT_CHIRALSET tetrahedralCarbons;

    MolOps::assignStereochemistry(*piece);
    _findChiralSets(*piece, chiralCenters, tetrahedralCarbons, coordMap);

    // if we have any chiral centers or are using random coordinates, we will
    // first embed the molecule in four dimensions, otherwise we will use 3D
    bool fourD = false;
    if (useRandomCoords || chiralCenters.size() > 0) {
      fourD = true;
    }
#ifdef RDK_THREADSAFE_SSS
    boost::thread_group tg;
#endif
    numThreads = getNumThreadsToUse(numThreads);

    detail::EmbedArgs eargs = {&confsOk,
                               fourD,
                               &fragMapping,
                               &confs,
                               fragIdx,
                               mmat,
                               useRandomCoords,
                               boxSizeMult,
                               randNegEig,
                               numZeroFail,
                               optimizerForceTol,
                               basinThresh,
                               seed,
                               maxIterations,
                               &chiralCenters,
                               &tetrahedralCarbons,
                               enforceChirality,
                               useExpTorsionAnglePrefs,
                               useBasicKnowledge,
                               &bonds,
                               &angles,
                               &expTorsionAtoms,
                               &expTorsionAngles,
                               &improperAtoms,
                               &atomNums};
    if (numThreads == 1) {
      detail::embedHelper_(0, 1, &eargs);
    }
#ifdef RDK_THREADSAFE_SSS
    else {
      for (int tid = 0; tid < numThreads; ++tid) {
        tg.add_thread(
            new boost::thread(detail::embedHelper_, tid, numThreads, &eargs));
      }
      tg.join_all();
    }
#endif
  }
  for (unsigned int ci = 0; ci < confs.size(); ++ci) {
    Conformer *conf = confs[ci];
    if (confsOk[ci]) {
      // check if we are pruning away conformations and
      // a closeby conformation has already been chosen :
      if (pruneRmsThresh > 0.0 &&
          !_isConfFarFromRest(mol, *conf, pruneRmsThresh,
                              onlyHeavyAtomsForRMS)) {
        delete conf;
      } else {
        int confId = (int)mol.addConformer(conf, true);
        res.push_back(confId);
      }
    } else {
      delete conf;
    }
  }
}

INT_VECT EmbedMultipleConfs(
    ROMol &mol, unsigned int numConfs, unsigned int maxIterations, int seed,
    bool clearConfs, bool useRandomCoords, double boxSizeMult, bool randNegEig,
    unsigned int numZeroFail, double pruneRmsThresh,
    const std::map<int, RDGeom::Point3D> *coordMap, double optimizerForceTol,
    bool ignoreSmoothingFailures, bool enforceChirality,
    bool useExpTorsionAnglePrefs, bool useBasicKnowledge, bool verbose,
    double basinThresh, bool onlyHeavyAtomsForRMS) {
  INT_VECT res;
  EmbedMultipleConfs(mol, res, numConfs, 1, maxIterations, seed, clearConfs,
                     useRandomCoords, boxSizeMult, randNegEig, numZeroFail,
                     pruneRmsThresh, coordMap, optimizerForceTol,
                     ignoreSmoothingFailures, enforceChirality,
                     useExpTorsionAnglePrefs, useBasicKnowledge, verbose,
                     basinThresh, onlyHeavyAtomsForRMS);
  return res;
}
}  // end of namespace DGeomHelpers
}  // end of namespace RDKit
