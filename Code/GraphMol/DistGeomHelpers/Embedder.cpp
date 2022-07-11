//
//  Copyright (C) 2004-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//#define DEBUG_EMBEDDING 1
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
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/dynamic_bitset.hpp>
#include <iomanip>
#include <RDGeneral/RDThreads.h>

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <future>
#endif

//#define DEBUG_EMBEDDING 1

namespace {
const double ERROR_TOL = 0.00001;
// these tolerances, all to detect and filter out bogus conformations, are a
// delicate balance between sensitive enough to detect obviously bad
// conformations but not so sensitive that a bunch of ok conformations get
// filtered out, which slows down the whole conformation generation process
const double MAX_MINIMIZED_E_PER_ATOM = 0.05;
const double MAX_MINIMIZED_E_CONTRIB = 0.20;
const double MIN_TETRAHEDRAL_CHIRAL_VOL = 0.50;
const double TETRAHEDRAL_CENTERINVOLUME_TOL = 0.30;
}  // namespace

namespace RDKit {
namespace DGeomHelpers {

//! Parameters corresponding to Sereina Riniker's KDG approach
const EmbedParameters KDG(0,        // maxIterations
                          1,        // numThreads
                          -1,       // randomSeed
                          true,     // clearConfs
                          false,    // useRandomCoords
                          2.0,      // boxSizeMult
                          true,     // randNegEig
                          1,        // numZeroFail
                          nullptr,  // coordMap
                          1e-3,     // optimizerForceTol
                          false,    // ignoreSmoothingFailures
                          true,     // enforceChirality
                          false,    // useExpTorsionAnglePrefs
                          true,     // useBasicKnowledge
                          false,    // verbose
                          5.0,      // basinThresh
                          -1.0,     // pruneRmsThresh
                          true,     // onlyHeavyAtomsForRMS
                          1,        // ETversion
                          nullptr,  // boundsMat
                          true,     // embedFragmentsSeparately
                          false,    // useSmallRingTorsions
                          false,    // useMacrocycleTorsions
                          false,    // useMacrocycle14config
                          nullptr,  // CPCI
                          nullptr   // callback
);

//! Parameters corresponding to Sereina Riniker's ETDG approach
const EmbedParameters ETDG(0,        // maxIterations
                           1,        // numThreads
                           -1,       // randomSeed
                           true,     // clearConfs
                           false,    // useRandomCoords
                           2.0,      // boxSizeMult
                           true,     // randNegEig
                           1,        // numZeroFail
                           nullptr,  // coordMap
                           1e-3,     // optimizerForceTol
                           false,    // ignoreSmoothingFailures
                           false,    // enforceChirality
                           true,     // useExpTorsionAnglePrefs
                           false,    // useBasicKnowledge
                           false,    // verbose
                           5.0,      // basinThresh
                           -1.0,     // pruneRmsThresh
                           true,     // onlyHeavyAtomsForRMS
                           1,        // ETversion
                           nullptr,  // boundsMat
                           true,     // embedFragmentsSeparately
                           false,    // useSmallRingTorsions
                           false,    // useMacrocycleTorsions
                           false,    // useMacrocycle14config
                           nullptr,  // CPCI
                           nullptr   // callback
);
//! Parameters corresponding to Sereina Riniker's ETKDG approach
const EmbedParameters ETKDG(0,        // maxIterations
                            1,        // numThreads
                            -1,       // randomSeed
                            true,     // clearConfs
                            false,    // useRandomCoords
                            2.0,      // boxSizeMult
                            true,     // randNegEig
                            1,        // numZeroFail
                            nullptr,  // coordMap
                            1e-3,     // optimizerForceTol
                            false,    // ignoreSmoothingFailures
                            true,     // enforceChirality
                            true,     // useExpTorsionAnglePrefs
                            true,     // useBasicKnowledge
                            false,    // verbose
                            5.0,      // basinThresh
                            -1.0,     // pruneRmsThresh
                            true,     // onlyHeavyAtomsForRMS
                            1,        // ETversion
                            nullptr,  // boundsMat
                            true,     // embedFragmentsSeparately
                            false,    // useSmallRingTorsions
                            false,    // useMacrocycleTorsions
                            false,    // useMacrocycle14config
                            nullptr,  // CPCI
                            nullptr   // callback
);

//! Parameters corresponding to Sereina Riniker's ETKDG approach - version 2
const EmbedParameters ETKDGv2(0,        // maxIterations
                              1,        // numThreads
                              -1,       // randomSeed
                              true,     // clearConfs
                              false,    // useRandomCoords
                              2.0,      // boxSizeMult
                              true,     // randNegEig
                              1,        // numZeroFail
                              nullptr,  // coordMap
                              1e-3,     // optimizerForceTol
                              false,    // ignoreSmoothingFailures
                              true,     // enforceChirality
                              true,     // useExpTorsionAnglePrefs
                              true,     // useBasicKnowledge
                              false,    // verbose
                              5.0,      // basinThresh
                              -1.0,     // pruneRmsThresh
                              true,     // onlyHeavyAtomsForRMS
                              2,        // ETversion
                              nullptr,  // boundsMat
                              true,     // embedFragmentsSeparately
                              false,    // useSmallRingTorsions
                              false,    // useMacrocycleTorsions
                              false,    // useMacrocycle14config
                              nullptr,  // CPCI
                              nullptr   // callback
);

//! Parameters corresponding improved ETKDG by Wang, Witek, Landrum and Riniker
//! (10.1021/acs.jcim.0c00025) - the macrocycle part
const EmbedParameters ETKDGv3(0,        // maxIterations
                              1,        // numThreads
                              -1,       // randomSeed
                              true,     // clearConfs
                              false,    // useRandomCoords
                              2.0,      // boxSizeMult
                              true,     // randNegEig
                              1,        // numZeroFail
                              nullptr,  // coordMap
                              1e-3,     // optimizerForceTol
                              false,    // ignoreSmoothingFailures
                              true,     // enforceChirality
                              true,     // useExpTorsionAnglePrefs
                              true,     // useBasicKnowledge
                              false,    // verbose
                              5.0,      // basinThresh
                              -1.0,     // pruneRmsThresh
                              true,     // onlyHeavyAtomsForRMS
                              2,        // ETversion
                              nullptr,  // boundsMat
                              true,     // embedFragmentsSeparately
                              false,    // useSmallRingTorsions
                              true,     // useMacrocycleTorsions
                              true,     // useMacrocycle14config
                              nullptr,  // CPCI
                              nullptr   // callback
);

//! Parameters corresponding improved ETKDG by Wang, Witek, Landrum and Riniker
//! (10.1021/acs.jcim.0c00025) - the small ring part
const EmbedParameters srETKDGv3(0,        // maxIterations
                                1,        // numThreads
                                -1,       // randomSeed
                                true,     // clearConfs
                                false,    // useRandomCoords
                                2.0,      // boxSizeMult
                                true,     // randNegEig
                                1,        // numZeroFail
                                nullptr,  // coordMap
                                1e-3,     // optimizerForceTol
                                false,    // ignoreSmoothingFailures
                                true,     // enforceChirality
                                true,     // useExpTorsionAnglePrefs
                                true,     // useBasicKnowledge
                                false,    // verbose
                                5.0,      // basinThresh
                                -1.0,     // pruneRmsThresh
                                true,     // onlyHeavyAtomsForRMS
                                2,        // ETversion
                                nullptr,  // boundsMat
                                true,     // embedFragmentsSeparately
                                true,     // useSmallRingTorsions
                                false,    // useMacrocycleTorsions
                                false,    // useMacrocycle14config
                                nullptr,  // CPCI
                                nullptr   // callback
);

namespace detail {
struct EmbedArgs {
  boost::dynamic_bitset<> *confsOk;
  bool fourD;
  INT_VECT *fragMapping;
  std::vector<std::unique_ptr<Conformer>> *confs;
  unsigned int fragIdx;
  DistGeom::BoundsMatPtr mmat;
  DistGeom::VECT_CHIRALSET const *chiralCenters;
  DistGeom::VECT_CHIRALSET const *tetrahedralCarbons;
  ForceFields::CrystalFF::CrystalFFDetails *etkdgDetails;
};
}  // namespace detail

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
  if (verbose) {
    std::cerr << "   " << fabs(vol) << std::endl;
  }
  if (fabs(vol) < MIN_TETRAHEDRAL_CHIRAL_VOL) {
    return false;
  }
  crossp = v1.crossProduct(v2);
  vol = crossp.dotProduct(v4);
  if (verbose) {
    std::cerr << "   " << fabs(vol) << std::endl;
  }
  if (fabs(vol) < MIN_TETRAHEDRAL_CHIRAL_VOL) {
    return false;
  }
  crossp = v1.crossProduct(v3);
  vol = crossp.dotProduct(v4);
  if (verbose) {
    std::cerr << "   " << fabs(vol) << std::endl;
  }
  if (fabs(vol) < MIN_TETRAHEDRAL_CHIRAL_VOL) {
    return false;
  }
  crossp = v2.crossProduct(v3);
  vol = crossp.dotProduct(v4);
  if (verbose) {
    std::cerr << "   " << fabs(vol) << std::endl;
  }
  return fabs(vol) >= MIN_TETRAHEDRAL_CHIRAL_VOL;
}

bool _sameSide(const RDGeom::Point3D &v1, const RDGeom::Point3D &v2,
               const RDGeom::Point3D &v3, const RDGeom::Point3D &v4,
               const RDGeom::Point3D &p0, double tol = 0.1) {
  RDGeom::Point3D normal = (v2 - v1).crossProduct(v3 - v1);
  double d1 = normal.dotProduct(v4 - v1);
  double d2 = normal.dotProduct(p0 - v1);
  // std::cerr << "     " << d1 << " - " << d2 << std::endl;
  if (fabs(d1) < tol || fabs(d2) < tol) {
    return false;
  }
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

namespace EmbeddingOps {
bool generateInitialCoords(RDGeom::PointPtrVect *positions,
                           const detail::EmbedArgs &eargs,
                           const EmbedParameters &embedParams,
                           RDNumeric::DoubleSymmMatrix &distMat,
                           RDKit::double_source_type *rng) {
  bool gotCoords = false;
  if (!embedParams.useRandomCoords) {
    double largestDistance =
        DistGeom::pickRandomDistMat(*eargs.mmat, distMat, *rng);
    RDUNUSED_PARAM(largestDistance);
    gotCoords = DistGeom::computeInitialCoords(distMat, *positions, *rng,
                                               embedParams.randNegEig,
                                               embedParams.numZeroFail);
  } else {
    double boxSize;
    if (embedParams.boxSizeMult > 0) {
      boxSize = 5. * embedParams.boxSizeMult;
    } else {
      boxSize = -1 * embedParams.boxSizeMult;
    }
    gotCoords = DistGeom::computeRandomCoords(*positions, boxSize, *rng);
    if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
      for (const auto &v : *embedParams.coordMap) {
        auto p = positions->at(v.first);
        for (unsigned int ci = 0; ci < v.second.dimension(); ++ci) {
          (*p)[ci] = v.second[ci];
        }
        // zero out any higher dimensional components:
        for (unsigned int ci = v.second.dimension(); ci < p->dimension();
             ++ci) {
          (*p)[ci] = 0.0;
        }
      }
    }
  }
  return gotCoords;
}
bool firstMinimization(RDGeom::PointPtrVect *positions,
                       const detail::EmbedArgs &eargs,
                       const EmbedParameters &embedParams) {
  bool gotCoords = true;
  boost::dynamic_bitset<> fixedPts(positions->size());
  if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
    for (const auto &v : *embedParams.coordMap) {
      fixedPts.set(v.first);
    }
  }
  std::unique_ptr<ForceFields::ForceField> field(DistGeom::constructForceField(
      *eargs.mmat, *positions, *eargs.chiralCenters, 1.0, 0.1, nullptr,
      embedParams.basinThresh, &fixedPts));
  if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
    for (const auto &v : *embedParams.coordMap) {
      field->fixedPoints().push_back(v.first);
    }
  }
  unsigned int nPasses = 0;
  field->initialize();
  if (field->calcEnergy() > ERROR_TOL) {
    int needMore = 1;
    while (needMore) {
      needMore = field->minimize(400, embedParams.optimizerForceTol);
      ++nPasses;
    }
  }
  std::vector<double> e_contribs;
  double local_e = field->calcEnergy(&e_contribs);

#ifdef DEBUG_EMBEDDING
  std::cerr << " Energy : " << local_e / positions->size() << " "
            << *(std::max_element(e_contribs.begin(), e_contribs.end()))
            << std::endl;
#endif

  // check that neither the energy nor any of the contributions to it are
  // too high (this is part of github #971)
  if (local_e / positions->size() >= MAX_MINIMIZED_E_PER_ATOM ||
      (e_contribs.size() &&
       *(std::max_element(e_contribs.begin(), e_contribs.end())) >
           MAX_MINIMIZED_E_CONTRIB)) {
#ifdef DEBUG_EMBEDDING
    std::cerr << " Energy fail: " << local_e / positions->size() << " "
              << *(std::max_element(e_contribs.begin(), e_contribs.end()))
              << std::endl;
#endif
    gotCoords = false;
  }
  return gotCoords;
}

bool checkTetrahedralCenters(const RDGeom::PointPtrVect *positions,
                             const detail::EmbedArgs &eargs,
                             const EmbedParameters &) {
  // for each of the atoms in the "tetrahedralCarbons" list, make sure
  // that there is a minimum volume around them and that they are inside
  // that volume. (this is part of github #971)
  for (const auto &tetSet : *eargs.tetrahedralCarbons) {
    // it could happen that the centroid is outside the volume defined
    // by the other
    // four points. That is also a fail.
    if (!_volumeTest(tetSet, *positions) ||
        !_centerInVolume(tetSet, *positions, TETRAHEDRAL_CENTERINVOLUME_TOL)) {
#ifdef DEBUG_EMBEDDING
      std::cerr << " fail2! (" << tetSet->d_idx0 << ") iter: "  //<< iter
                << " vol: " << _volumeTest(tetSet, *positions, true)
                << " center: "
                << _centerInVolume(tetSet, *positions,
                                   TETRAHEDRAL_CENTERINVOLUME_TOL, true)
                << std::endl;
#endif
      return false;
    }
  }
  return true;
}
bool checkChiralCenters(const RDGeom::PointPtrVect *positions,
                        const detail::EmbedArgs &eargs,
                        const EmbedParameters &) {
  // check the chiral volume:
  for (const auto &chiralSet : *eargs.chiralCenters) {
    double vol = DistGeom::ChiralViolationContrib::calcChiralVolume(
        chiralSet->d_idx1, chiralSet->d_idx2, chiralSet->d_idx3,
        chiralSet->d_idx4, *positions);
    double lb = chiralSet->getLowerVolumeBound();
    double ub = chiralSet->getUpperVolumeBound();
    if ((lb > 0 && vol < lb && (lb - vol) / lb > .2) ||
        (ub < 0 && vol > ub && (vol - ub) / ub > .2)) {
#ifdef DEBUG_EMBEDDING
      std::cerr << " fail! (" << chiralSet->d_idx0 << ") iter: "
                << " " << vol << " " << lb << "-" << ub << std::endl;
#endif
      return false;
    }
  }
  return true;
}
bool minimizeFourthDimension(RDGeom::PointPtrVect *positions,
                             const detail::EmbedArgs &eargs,
                             const EmbedParameters &embedParams) {
  // now redo the minimization if we have a chiral center
  // or have started from random coords. This
  // time removing the chiral constraints and
  // increasing the weight on the fourth dimension

  std::unique_ptr<ForceFields::ForceField> field2(DistGeom::constructForceField(
      *eargs.mmat, *positions, *eargs.chiralCenters, 0.2, 1.0, nullptr,
      embedParams.basinThresh));
  if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
    for (const auto &v : *embedParams.coordMap) {
      field2->fixedPoints().push_back(v.first);
    }
  }

  field2->initialize();
  // std::cerr<<"FIELD2 E: "<<field2->calcEnergy()<<std::endl;
  if (field2->calcEnergy() > ERROR_TOL) {
    int needMore = 1;
    int nPasses2 = 0;
    while (needMore) {
      needMore = field2->minimize(200, embedParams.optimizerForceTol);
      ++nPasses2;
    }
    // std::cerr<<"   "<<field2->calcEnergy()<<" after npasses2:
    // "<<nPasses2<<std::endl;
  }
  return true;
}
// the minimization using experimental torsion angle preferences
bool minimizeWithExpTorsions(RDGeom::PointPtrVect &positions,
                             const detail::EmbedArgs &eargs,
                             const EmbedParameters &embedParams) {
  PRECONDITION(eargs.etkdgDetails, "bogus etkdgDetails pointer");
  bool planar = true;

  // convert to 3D positions and create coordMap
  RDGeom::Point3DPtrVect positions3D;
  for (auto &position : positions) {
    positions3D.push_back(
        new RDGeom::Point3D((*position)[0], (*position)[1], (*position)[2]));
  }

  // create the force field
  std::unique_ptr<ForceFields::ForceField> field;
  if (embedParams.useBasicKnowledge) {  // ETKDG or KDG
    if (embedParams.CPCI != nullptr) {
      field.reset(DistGeom::construct3DForceField(
          *eargs.mmat, positions3D, *eargs.etkdgDetails, *embedParams.CPCI));
    } else {
      field.reset(DistGeom::construct3DForceField(*eargs.mmat, positions3D,
                                                  *eargs.etkdgDetails));
    }
  } else {  // plain ETDG
    field.reset(DistGeom::constructPlain3DForceField(*eargs.mmat, positions3D,
                                                     *eargs.etkdgDetails));
  }
  if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
    for (const auto &v : *embedParams.coordMap) {
      field->fixedPoints().push_back(v.first);
    }
  }

  // minimize!
  field->initialize();
  if (field->calcEnergy() > ERROR_TOL) {
    // while (needMore) {
    field->minimize(300, embedParams.optimizerForceTol);
    //      ++nPasses;
    //}
  }
  // std::cout << field->calcEnergy() << std::endl;

  // check for planarity if ETKDG or KDG
  if (embedParams.useBasicKnowledge) {
    // create a force field with only the impropers
    std::unique_ptr<ForceFields::ForceField> field2(
        DistGeom::construct3DImproperForceField(
            *eargs.mmat, positions3D, eargs.etkdgDetails->improperAtoms,
            eargs.etkdgDetails->atomNums));
    if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
      for (const auto &v : *embedParams.coordMap) {
        field2->fixedPoints().push_back(v.first);
      }
    }

    field2->initialize();
    // check if the energy is low enough
    double planarityTolerance = 0.7;
    if (field2->calcEnergy() >
        eargs.etkdgDetails->improperAtoms.size() * planarityTolerance) {
#ifdef DEBUG_EMBEDDING
      std::cerr << "   planar fail: " << field2->calcEnergy() << " "
                << eargs.etkdgDetails->improperAtoms.size() * planarityTolerance
                << std::endl;
#endif
      planar = false;
    }
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

bool finalChiralChecks(RDGeom::PointPtrVect *positions,
                       const detail::EmbedArgs &eargs,
                       const EmbedParameters &) {
  // "distance matrix" chirality test
  std::set<int> atoms;
  for (const auto &chiralSet : *eargs.chiralCenters) {
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
    if (!_boundsFulfilled(atomsToCheck, *eargs.mmat, *positions)) {
#ifdef DEBUG_EMBEDDING
      std::cerr << " fail3a! (" << atomsToCheck[0] << ") iter: "  //<< iter
                << std::endl;
#endif
      return false;
    }
  }

  // "center in volume" chirality test
  for (const auto &chiralSet : *eargs.chiralCenters) {
    // it could happen that the centroid is outside the volume defined
    // by the other four points. That is also a fail.
    if (!_centerInVolume(chiralSet, *positions)) {
#ifdef DEBUG_EMBEDDING
      std::cerr << " fail3b! (" << chiralSet->d_idx0 << ") iter: "  //<< iter
                << std::endl;
#endif
      return false;
    }
  }

  return true;
}

bool embedPoints(RDGeom::PointPtrVect *positions, detail::EmbedArgs eargs,
                 EmbedParameters embedParams, int seed) {
  PRECONDITION(positions, "bogus positions");
  if (embedParams.maxIterations == 0) {
    embedParams.maxIterations = 10 * positions->size();
  }
  RDNumeric::DoubleSymmMatrix distMat(positions->size(), 0.0);

  // The basin threshold just gets us into trouble when we're using
  // random coordinates since it ends up ignoring 1-4 (and higher)
  // interactions. This causes us to get folded-up (and self-penetrating)
  // conformations for large flexible molecules
  if (embedParams.useRandomCoords) {
    embedParams.basinThresh = 1e8;
  }

  RDKit::double_source_type *rng = nullptr;
  RDKit::rng_type *generator;
  RDKit::uniform_double *distrib;
  CHECK_INVARIANT(seed >= -1,
                  "random seed must either be positive, zero, or negative one");
  if (seed > -1) {
    generator = new RDKit::rng_type(42u);
    generator->seed(seed);
    distrib = new RDKit::uniform_double(0.0, 1.0);
    rng = new RDKit::double_source_type(*generator, *distrib);
  } else {
    rng = &RDKit::getDoubleRandomSource();
  }

  bool gotCoords = false;
  unsigned int iter = 0;
  while (!gotCoords && iter < embedParams.maxIterations) {
    ++iter;
    if (embedParams.callback != nullptr) {
      embedParams.callback(iter);
    }
    gotCoords = EmbeddingOps::generateInitialCoords(positions, eargs,
                                                    embedParams, distMat, rng);
#ifdef DEBUG_EMBEDDING
    if (!gotCoords) {
      std::cerr << "Initial embedding failed!, Iter: " << iter << std::endl;
    }
#endif
    if (gotCoords) {
      gotCoords =
          EmbeddingOps::firstMinimization(positions, eargs, embedParams);
      if (gotCoords) {
        gotCoords = EmbeddingOps::checkTetrahedralCenters(positions, eargs,
                                                          embedParams);
      }

      // Check if any of our chiral centers are badly out of whack.
      if (gotCoords && embedParams.enforceChirality &&
          eargs.chiralCenters->size() > 0) {
        gotCoords =
            EmbeddingOps::checkChiralCenters(positions, eargs, embedParams);
      }
      // redo the minimization if we have a chiral center
      // or have started from random coords.
      if (gotCoords &&
          (eargs.chiralCenters->size() > 0 || embedParams.useRandomCoords)) {
        gotCoords = EmbeddingOps::minimizeFourthDimension(positions, eargs,
                                                          embedParams);
      }

      // (ET)(K)DG
      if (gotCoords && (embedParams.useExpTorsionAnglePrefs ||
                        embedParams.useBasicKnowledge)) {
        gotCoords = EmbeddingOps::minimizeWithExpTorsions(*positions, eargs,
                                                          embedParams);
      }
      // test if chirality is correct
      if (embedParams.enforceChirality && gotCoords &&
          (eargs.chiralCenters->size() > 0)) {
        gotCoords =
            EmbeddingOps::finalChiralChecks(positions, eargs, embedParams);
      }
    }  // if(gotCoords)
  }    // while
  if (seed > -1) {
    delete rng;
    delete generator;
    delete distrib;
  }
  return gotCoords;
}

void findChiralSets(const ROMol &mol, DistGeom::VECT_CHIRALSET &chiralCenters,
                    DistGeom::VECT_CHIRALSET &tetrahedralCenters,
                    const std::map<int, RDGeom::Point3D> *coordMap) {
  for (const auto &atom : mol.atoms()) {
    if (atom->getAtomicNum() != 1) {  // skip hydrogens
      Atom::ChiralType chiralType = atom->getChiralTag();
      if ((chiralType == Atom::CHI_TETRAHEDRAL_CW ||
           chiralType == Atom::CHI_TETRAHEDRAL_CCW) ||
          ((atom->getAtomicNum() == 6 || atom->getAtomicNum() == 7) &&
           atom->getDegree() == 4)) {
        // make a chiral set from the neighbors
        INT_VECT nbrs;
        nbrs.reserve(4);
        // find the neighbors of this atom and enter them into the
        // nbr list
        ROMol::OEDGE_ITER beg, end;
        boost::tie(beg, end) = mol.getAtomBonds(atom);
        while (beg != end) {
          nbrs.push_back(mol[*beg]->getOtherAtom(atom)->getIdx());
          ++beg;
        }
        // if we have less than 4 heavy atoms as neighbors,
        // we need to include the chiral center into the mix
        // we should at least have 3 though
        CHECK_INVARIANT(nbrs.size() >= 3, "Cannot be a chiral center");

        if (nbrs.size() < 4) {
          nbrs.insert(nbrs.end(), atom->getIdx());
        }

        // now create a chiral set and set the upper and lower bound on the
        // volume
        if (chiralType == Atom::CHI_TETRAHEDRAL_CCW) {
          // positive chiral volume
          auto *cset = new DistGeom::ChiralSet(atom->getIdx(), nbrs[0], nbrs[1],
                                               nbrs[2], nbrs[3], 5.0, 100.0);
          DistGeom::ChiralSetPtr cptr(cset);
          chiralCenters.push_back(cptr);
        } else if (chiralType == Atom::CHI_TETRAHEDRAL_CW) {
          auto *cset = new DistGeom::ChiralSet(atom->getIdx(), nbrs[0], nbrs[1],
                                               nbrs[2], nbrs[3], -100.0, -5.0);
          DistGeom::ChiralSetPtr cptr(cset);
          chiralCenters.push_back(cptr);
        } else {
          if ((coordMap && coordMap->find(atom->getIdx()) != coordMap->end()) ||
              (mol.getRingInfo()->isInitialized() &&
               (mol.getRingInfo()->numAtomRings(atom->getIdx()) < 2 ||
                mol.getRingInfo()->isAtomInRingOfSize(atom->getIdx(), 3)))) {
            // we only want to these tests for ring atoms that are not part of
            // the coordMap
            // there's no sense doing 3-rings because those are a nightmare
          } else {
            auto *cset = new DistGeom::ChiralSet(
                atom->getIdx(), nbrs[0], nbrs[1], nbrs[2], nbrs[3], 0.0, 0.0);
            DistGeom::ChiralSetPtr cptr(cset);
            tetrahedralCenters.push_back(cptr);
          }
        }
      }  // if block -chirality check
    }    // if block - heavy atom check
  }      // for loop over atoms
}  // end of _findChiralSets

void adjustBoundsMatFromCoordMap(
    DistGeom::BoundsMatPtr mmat, unsigned int,
    const std::map<int, RDGeom::Point3D> *coordMap) {
  for (auto iIt = coordMap->begin(); iIt != coordMap->end(); ++iIt) {
    unsigned int iIdx = iIt->first;
    const RDGeom::Point3D &iPoint = iIt->second;
    auto jIt = iIt;
    while (++jIt != coordMap->end()) {
      unsigned int jIdx = jIt->first;
      const RDGeom::Point3D &jPoint = jIt->second;
      double dist = (iPoint - jPoint).length();
      mmat->setUpperBound(iIdx, jIdx, dist);
      mmat->setLowerBound(iIdx, jIdx, dist);
    }
  }
}

void initETKDG(ROMol *mol, const EmbedParameters &params,
               ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails) {
  PRECONDITION(mol, "bad molecule");
  unsigned int nAtoms = mol->getNumAtoms();
  if (params.useExpTorsionAnglePrefs || params.useBasicKnowledge) {
    ForceFields::CrystalFF::getExperimentalTorsions(
        *mol, etkdgDetails, params.useExpTorsionAnglePrefs,
        params.useSmallRingTorsions, params.useMacrocycleTorsions,
        params.useBasicKnowledge, params.ETversion, params.verbose);
    etkdgDetails.atomNums.resize(nAtoms);
    for (unsigned int i = 0; i < nAtoms; ++i) {
      etkdgDetails.atomNums[i] = mol->getAtomWithIdx(i)->getAtomicNum();
    }
  }
  etkdgDetails.boundsMatForceScaling = params.boundsMatForceScaling;
}

bool setupInitialBoundsMatrix(
    ROMol *mol, DistGeom::BoundsMatPtr mmat,
    const std::map<int, RDGeom::Point3D> *coordMap,
    const EmbedParameters &params,
    ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails) {
  PRECONDITION(mol, "bad molecule");
  unsigned int nAtoms = mol->getNumAtoms();
  if (params.useExpTorsionAnglePrefs || params.useBasicKnowledge) {
    setTopolBounds(*mol, mmat, etkdgDetails.bonds, etkdgDetails.angles, true,
                   false, params.useMacrocycle14config,
                   params.forceTransAmides);
  } else {
    setTopolBounds(*mol, mmat, true, false, params.useMacrocycle14config,
                   params.forceTransAmides);
  }
  double tol = 0.0;
  if (coordMap) {
    adjustBoundsMatFromCoordMap(mmat, nAtoms, coordMap);
    tol = 0.05;
  }
  if (!DistGeom::triangleSmoothBounds(mmat, tol)) {
    // ok this bound matrix failed to triangle smooth - re-compute the
    // bounds matrix without 15 bounds and with VDW scaling
    initBoundsMat(mmat);
    setTopolBounds(*mol, mmat, false, true, params.useMacrocycle14config,
                   params.forceTransAmides);

    if (coordMap) {
      adjustBoundsMatFromCoordMap(mmat, nAtoms, coordMap);
    }

    // try triangle smoothing again
    if (!DistGeom::triangleSmoothBounds(mmat, tol)) {
      // ok, we're not going to be able to smooth this,
      if (params.ignoreSmoothingFailures) {
        // proceed anyway with the more relaxed bounds matrix
        initBoundsMat(mmat);
        setTopolBounds(*mol, mmat, false, true, params.useMacrocycle14config,
                       params.forceTransAmides);

        if (coordMap) {
          adjustBoundsMatFromCoordMap(mmat, nAtoms, coordMap);
        }
      } else {
        BOOST_LOG(rdWarningLog)
            << "Could not triangle bounds smooth molecule." << std::endl;
        return false;
      }
    }
  }
  return true;
}
}  // namespace EmbeddingOps

void _fillAtomPositions(RDGeom::Point3DConstPtrVect &pts, const Conformer &conf,
                        const ROMol &, const std::vector<unsigned int> &match) {
  PRECONDITION(pts.size() == match.size(), "bad pts size");
  for (unsigned int i = 0; i < match.size(); i++) {
    pts[i] = &conf.getAtomPos(match[i]);
  }
}

bool _isConfFarFromRest(
    const ROMol &mol, const Conformer &conf, double threshold,
    const std::vector<std::vector<unsigned int>> &selfMatches) {
  // NOTE: it is tempting to use some triangle inequality to prune
  // conformations here but some basic testing has shown very
  // little advantage and given that the time for pruning fades in
  // comparison to embedding - we will use a simple for loop below
  // over all conformation until we find a match
  RDGeom::Point3DConstPtrVect refPoints(selfMatches[0].size());
  RDGeom::Point3DConstPtrVect prbPoints(selfMatches[0].size());
  _fillAtomPositions(refPoints, conf, mol, selfMatches[0]);

  double ssrThres = conf.getNumAtoms() * threshold * threshold;
  for (const auto &match : selfMatches) {
    for (auto confi = mol.beginConformers(); confi != mol.endConformers();
         ++confi) {
      _fillAtomPositions(prbPoints, *(*confi), mol, match);
      RDGeom::Transform3D trans;
      auto ssr =
          RDNumeric::Alignments::AlignPoints(refPoints, prbPoints, trans);
      if (ssr < ssrThres) {
        return false;
      }
    }
  }
  return true;
}

namespace detail {

template <class T>
bool multiplication_overflows_(T a, T b) {
  // a * b > c if and only if a > c / b
  if (a == 0 || b == 0) {
    return false;
  }
  return a > std::numeric_limits<T>::max() / b;
}

void embedHelper_(int threadId, int numThreads, EmbedArgs *eargs,
                  const EmbedParameters *params) {
  PRECONDITION(eargs, "bogus eargs");
  PRECONDITION(params, "bogus params");
  unsigned int nAtoms = eargs->mmat->numRows();
  RDGeom::PointPtrVect positions(nAtoms);
  // we might thrown an exception in a callback
  // in order to avoid leaking the points we're working with
  // allocate them with unique_ptrs and then work with the naked
  // pointers from those
  std::vector<std::unique_ptr<RDGeom::Point>> positionsStore;
  positionsStore.reserve(nAtoms);
  for (unsigned int i = 0; i < nAtoms; ++i) {
    if (eargs->fourD) {
      positionsStore.emplace_back(new RDGeom::PointND(4));
    } else {
      positionsStore.emplace_back(new RDGeom::Point3D());
    }
    positions[i] = positionsStore[i].get();
  }
  for (size_t ci = 0; ci < eargs->confs->size(); ci++) {
    if (rdcast<int>(ci % numThreads) != threadId) {
      continue;
    }
    if (!(*eargs->confsOk)[ci]) {
      // we call this function for each fragment in a molecule,
      // if one of the fragments has already failed, there's no
      // sense in embedding this one
      continue;
    }

    CHECK_INVARIANT(
        params->randomSeed >= -1,
        "random seed must either be positive, zero, or negative one");
    int new_seed = params->randomSeed;
    if (new_seed > -1) {
      if (!multiplication_overflows_(rdcast<int>(ci + 1), params->randomSeed)) {
        // old method of computing a new seed
        new_seed = (ci + 1) * params->randomSeed;
      } else {
        // If the above simple multiplication will overflow, use a
        // cheap and easy way to hash the conformer index and seed
        // together: for N'ary numerical system, where N is the
        // maximum possible value of the pair of numbers. The
        // following will generate unique integers:
        // hash(a, b) = a + b * N
        auto big_seed = rdcast<size_t>(params->randomSeed);
        size_t max_val = std::max(ci + 1, big_seed);
        size_t big_num = big_seed + max_val * (ci + 1);
        // only grab the first 31 bits xor'd with the next 31 bits to
        // make sure its positive, careful, the 'ULL' is important
        // here, 0x7fffffff is the 'int' type because of C default
        // number semantics and that we definitely don't want!
        const size_t positive_int_mask = 0x7fffffffULL;
        size_t folded_num = (big_num & positive_int_mask) ^ (big_num >> 31ULL);
        new_seed = rdcast<int>(folded_num & positive_int_mask);
      }
    }
    CHECK_INVARIANT(new_seed >= -1,
                    "Something went wrong calculating a new seed");

    bool gotCoords =
        EmbeddingOps::embedPoints(&positions, *eargs, *params, new_seed);

    // copy the coordinates into the correct conformer
    if (gotCoords) {
      auto &conf = (*eargs->confs)[ci];
      unsigned int fragAtomIdx = 0;
      for (unsigned int i = 0; i < conf->getNumAtoms(); ++i) {
        if (!eargs->fragMapping ||
            (*eargs->fragMapping)[i] == static_cast<int>(eargs->fragIdx)) {
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
}

std::vector<std::vector<unsigned int>> getMolSelfMatches(
    const ROMol &mol, const EmbedParameters &params) {
  std::vector<std::vector<unsigned int>> res;
  if (params.pruneRmsThresh && params.useSymmetryForPruning) {
    RWMol tmol(mol);
    MolOps::RemoveHsParameters ps;
    bool sanitize = false;
    MolOps::removeHs(tmol, ps, sanitize);
    SubstructMatchParameters sssps;
    sssps.maxMatches = 1;
    // provides the atom indices in the molecule corresponding
    // to the indices in the H-stripped version
    auto strippedMatch = SubstructMatch(mol, tmol, sssps);
    CHECK_INVARIANT(strippedMatch.size() == 1, "expected match not found");

    sssps.maxMatches = 1000;
    sssps.uniquify = false;
    auto heavyAtomMatches = SubstructMatch(tmol, tmol, sssps);
    for (const auto &match : heavyAtomMatches) {
      res.emplace_back(0);
      res.back().reserve(match.size());
      for (auto midx : match) {
        res.back().push_back(strippedMatch[0][midx.second].second);
      }
    }
  } else if (params.onlyHeavyAtomsForRMS) {
    res.emplace_back(0);
    for (const auto &at : mol.atoms()) {
      if (at->getAtomicNum() != 1) {
        res.back().push_back(at->getIdx());
      }
    }
  } else {
    res.emplace_back(0);
    res.back().reserve(mol.getNumAtoms());
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      res.back().push_back(i);
    }
  }
  return res;
}

}  // end of namespace detail

void EmbedMultipleConfs(ROMol &mol, INT_VECT &res, unsigned int numConfs,
                        const EmbedParameters &params) {
  if (!mol.getNumAtoms()) {
    throw ValueErrorException("molecule has no atoms");
  }
  if (params.ETversion < 1 || params.ETversion > 2) {
    throw ValueErrorException(
        "Only version 1 and 2 of the experimental "
        "torsion-angle preferences (ETversion) supported");
  }

  if (MolOps::needsHs(mol)) {
    BOOST_LOG(rdWarningLog)
        << "Molecule does not have explicit Hs. Consider calling AddHs()"
        << std::endl;
  }

  // initialize the conformers we're going to be creating:
  if (params.clearConfs) {
    res.clear();
    mol.clearConformers();
  }
  std::vector<std::unique_ptr<Conformer>> confs;
  confs.reserve(numConfs);
  for (unsigned int i = 0; i < numConfs; ++i) {
    confs.emplace_back(new Conformer(mol.getNumAtoms()));
  }

  boost::dynamic_bitset<> confsOk(numConfs);
  confsOk.set();

  INT_VECT fragMapping;
  std::vector<ROMOL_SPTR> molFrags;
  if (params.embedFragmentsSeparately) {
    molFrags = MolOps::getMolFrags(mol, true, &fragMapping);
  } else {
    molFrags.push_back(ROMOL_SPTR(new ROMol(mol)));
    fragMapping.resize(mol.getNumAtoms());
    std::fill(fragMapping.begin(), fragMapping.end(), 0);
  }
  const std::map<int, RDGeom::Point3D> *coordMap = params.coordMap;
  if (molFrags.size() > 1 && coordMap) {
    BOOST_LOG(rdWarningLog)
        << "Constrained conformer generation (via the coordMap argument) "
           "does not work with molecules that have multiple fragments."
        << std::endl;
    coordMap = nullptr;
  }

  if (molFrags.size() > 1 && params.boundsMat != nullptr) {
    BOOST_LOG(rdWarningLog)
        << "Conformer generation using a user-provided boundsMat "
           "does not work with molecules that have multiple fragments. The "
           "boundsMat will be ignored."
        << std::endl;
    coordMap = nullptr;  // FIXME not directly related to ETKDG, but here I
                         // think it should be params.boundsMat = nullptr
  }

  // we will generate conformations for each fragment in the molecule
  // separately, so loop over them:
  for (unsigned int fragIdx = 0; fragIdx < molFrags.size(); ++fragIdx) {
    ROMOL_SPTR piece = molFrags[fragIdx];
    unsigned int nAtoms = piece->getNumAtoms();

    ForceFields::CrystalFF::CrystalFFDetails etkdgDetails;
    EmbeddingOps::initETKDG(piece.get(), params, etkdgDetails);

    DistGeom::BoundsMatPtr mmat;
    if (params.boundsMat == nullptr || molFrags.size() > 1) {
      // The user didn't provide one, so create and initialize the distance
      // bounds matrix
      mmat.reset(new DistGeom::BoundsMatrix(nAtoms));
      initBoundsMat(mmat);
      if (!EmbeddingOps::setupInitialBoundsMatrix(piece.get(), mmat, coordMap,
                                                  params, etkdgDetails)) {
        // return if we couldn't setup the bounds matrix
        // possible causes include a triangle smoothing failure
        return;
      }
    } else {
      // just use what they gave us
      // first make sure it's the right size though:
      if (params.boundsMat->numRows() != nAtoms) {
        throw ValueErrorException(
            "size of boundsMat provided does not match the number of atoms in "
            "the molecule.");
      }
      collectBondsAndAngles((*piece.get()), etkdgDetails.bonds,
                            etkdgDetails.angles);
      mmat.reset(new DistGeom::BoundsMatrix(*params.boundsMat));
    }

    // find all the chiral centers in the molecule
    MolOps::assignStereochemistry(*piece);
    DistGeom::VECT_CHIRALSET chiralCenters;
    DistGeom::VECT_CHIRALSET tetrahedralCarbons;
    EmbeddingOps::findChiralSets(*piece, chiralCenters, tetrahedralCarbons,
                                 coordMap);

    // if we have any chiral centers or are using random coordinates, we will
    // first embed the molecule in four dimensions, otherwise we will use 3D
    bool fourD = false;
    if (params.useRandomCoords || chiralCenters.size() > 0) {
      fourD = true;
    }
    int numThreads = getNumThreadsToUse(params.numThreads);

    // do the embedding, using multiple threads if requested
    detail::EmbedArgs eargs = {
        &confsOk, fourD,          &fragMapping,        &confs,       fragIdx,
        mmat,     &chiralCenters, &tetrahedralCarbons, &etkdgDetails};
    if (numThreads == 1) {
      detail::embedHelper_(0, 1, &eargs, &params);
    }
#ifdef RDK_BUILD_THREADSAFE_SSS
    else {
      std::vector<std::future<void>> tg;
      for (int tid = 0; tid < numThreads; ++tid) {
        tg.emplace_back(std::async(std::launch::async, detail::embedHelper_,
                                   tid, numThreads, &eargs, &params));
      }
      for (auto &fut : tg) {
        fut.get();
      }
    }
#endif
  }
  auto selfMatches = detail::getMolSelfMatches(mol, params);
  for (unsigned int ci = 0; ci < confs.size(); ++ci) {
    auto &conf = confs[ci];
    if (confsOk[ci]) {
      // check if we are pruning away conformations and
      // a close-by conformation has already been chosen :
      if (params.pruneRmsThresh <= 0.0 ||
          _isConfFarFromRest(mol, *conf, params.pruneRmsThresh, selfMatches)) {
        int confId = (int)mol.addConformer(conf.release(), true);
        res.push_back(confId);
      }
    }
  }
}

}  // end of namespace DGeomHelpers
}  // end of namespace RDKit
