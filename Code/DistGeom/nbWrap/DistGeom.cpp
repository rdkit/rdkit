//
//  Copyright (C) 2004-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <Geometry/point.h>
#include <Numerics/Matrix.h>
#include <Numerics/SymmMatrix.h>

#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/TriangleSmooth.h>
#include <DistGeom/DistGeomUtils.h>
#include <DistGeom/ChiralSet.h>

#include <ForceField/ForceField.h>

#include <vector>
#include <map>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {

bool doTriangleSmoothing(
    nb::ndarray<nb::numpy, double, nb::ndim<2>, nb::c_contig> boundsMat,
    double tol) {
  if (boundsMat.shape(0) != boundsMat.shape(1)) {
    throw nb::value_error("The array has to be square");
  }
  if (boundsMat.shape(0) == 0) {
    throw nb::value_error("The array has to have a nonzero size");
  }

  unsigned int nrows = (unsigned int)boundsMat.shape(0);
  unsigned int dSize = nrows * nrows;
  auto *cData = new double[dSize];
  auto *inData = boundsMat.data();
  memcpy(static_cast<void *>(cData), static_cast<const void *>(inData),
         dSize * sizeof(double));
  DistGeom::BoundsMatrix::DATA_SPTR sdata(cData);
  DistGeom::BoundsMatrix bm(nrows, sdata);

  bool res = DistGeom::triangleSmoothBounds(&bm, tol);
  memcpy(static_cast<void *>(inData), static_cast<const void *>(cData),
         dSize * sizeof(double));
  return res;
}

nb::ndarray<nb::numpy, double, nb::ndim<2>> embedBoundsMatrix(
    nb::ndarray<nb::numpy, double, nb::ndim<2>, nb::c_contig> boundsMat,
    int maxIters = 10, int randomizeOnFailure = 0, int numZeroFail = 2,
    nb::list weights = nb::list(), int randomSeed = -1) {
  if (boundsMat.shape(0) != boundsMat.shape(1)) {
    throw nb::value_error("The array has to be square");
  }
  if (boundsMat.shape(0) == 0) {
    throw nb::value_error("The array has to have a nonzero size");
  }

  unsigned int nrows = (unsigned int)boundsMat.shape(0);
  unsigned int dSize = nrows * nrows;
  auto *cData = new double[dSize];
  auto *inData = boundsMat.data();
  memcpy(static_cast<void *>(cData), static_cast<const void *>(inData),
         dSize * sizeof(double));

  DistGeom::BoundsMatrix::DATA_SPTR sdata(cData);
  DistGeom::BoundsMatrix bm(nrows, sdata);

  std::unique_ptr<RDGeom::Point3D[]> positions(new RDGeom::Point3D[nrows]);
  std::vector<RDGeom::Point *> posPtrs;
  for (unsigned int i = 0; i < nrows; i++) {
    posPtrs.push_back(&positions[i]);
  }

  RDNumeric::DoubleSymmMatrix distMat(nrows, 0.0);

  // ---- ---- ---- ---- ---- ---- ---- ---- ----
  // start the embedding:
  bool gotCoords = false;
  for (int iter = 0; iter < maxIters && !gotCoords; iter++) {
    // pick a random distance matrix
    DistGeom::pickRandomDistMat(bm, distMat, randomSeed);

    // and embed it:
    gotCoords = DistGeom::computeInitialCoords(
        distMat, posPtrs, (bool)randomizeOnFailure, numZeroFail, randomSeed);

    // update the seed:
    if (randomSeed >= 0) {
      randomSeed += iter * 999;
    }
  }

  if (gotCoords) {
    std::map<std::pair<int, int>, double> weightMap;
    unsigned int nElems = nb::len(weights);
    for (unsigned int entryIdx = 0; entryIdx < nElems; entryIdx++) {
      nb::object entry = weights[entryIdx];
      nb::sequence seq = nb::cast<nb::sequence>(entry);
      if (nb::len(seq) != 3) {
        throw nb::value_error(
            "weights argument must be a sequence of 3-sequences");
      }
      int idx1 = nb::cast<int>(seq[0]);
      int idx2 = nb::cast<int>(seq[1]);
      double w = nb::cast<double>(seq[2]);
      weightMap[std::make_pair(idx1, idx2)] = w;
    }
    DistGeom::VECT_CHIRALSET csets;
    ForceFields::ForceField *field =
        DistGeom::constructForceField(bm, posPtrs, csets, 0.0, 0.0, &weightMap);
    CHECK_INVARIANT(field, "could not build dgeom force field");
    field->initialize();
    if (field->calcEnergy() > 1e-5) {
      int needMore = 1;
      while (needMore) {
        needMore = field->minimize();
      }
    }
    delete field;
  } else {
    throw nb::value_error("could not embed matrix");
  }

  // ---- ---- ---- ---- ---- ---- ---- ---- ----
  // construct the results matrix:
  auto *resData = new double[nrows * 3];
  for (unsigned int i = 0; i < nrows; i++) {
    unsigned int iTab = i * 3;
    for (unsigned int j = 0; j < 3; ++j) {
      resData[iTab + j] = positions[i][j];
    }
  }
  nb::capsule owner(resData, [](void *f) noexcept {
    delete[] reinterpret_cast<double *>(f);
  });
  return nb::ndarray<nb::numpy, double, nb::ndim<2>>(resData, {nrows, 3},
                                                     owner);
}

}  // namespace RDKit

NB_MODULE(DistGeom, m) {
  m.doc() = "Module containing functions for basic distance geometry operations";

  m.def(
      "DoTriangleSmoothing", &RDKit::doTriangleSmoothing,
      "boundsMatrix"_a, "tol"_a = 0.,
      R"DOC(Do triangle smoothing on a bounds matrix

ARGUMENTS:

   - mat: a square Numeric array of doubles containing the bounds matrix, this matrix
          *is* modified by the smoothing

RETURNS:

   a boolean indicating whether or not the smoothing worked.
)DOC");

  m.def(
      "EmbedBoundsMatrix", &RDKit::embedBoundsMatrix,
      "boundsMatrix"_a, "maxIters"_a = 10, "randomizeOnFailure"_a = 0,
      "numZeroFail"_a = 2, "weights"_a = nb::list(), "randomSeed"_a = -1,
      R"DOC(Embed a bounds matrix and return the coordinates

ARGUMENTS:

   - boundsMatrix: a square Numeric array of doubles containing the bounds matrix, this matrix
          should already be smoothed
   - maxIters: (optional) the maximum number of random distance matrices to try
   - randomizeOnFailure: (optional) toggles using random coords if a matrix fails to embed
   - numZeroFail: (optional) sets the number of zero eigenvalues to be considered a failure
   - weights: (optional) a sequence of 3 sequences (i,j,weight) indicating elements of
      the bounds matrix whose weights should be adjusted
   - randomSeed: (optional) sets the random number seed used for embedding

RETURNS:

   a Numeric array of doubles with the coordinates
)DOC");
}
