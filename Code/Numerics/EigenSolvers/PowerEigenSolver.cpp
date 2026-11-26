//
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "PowerEigenSolver.h"
#include <Numerics/Vector.h>
#include <Numerics/Matrix.h>
#include <Numerics/SymmMatrix.h>
#include <RDGeneral/Invariant.h>
#include <algorithm>
#include <ctime>
#include <vector>

namespace RDNumeric {
namespace EigenSolvers {
namespace {
// copy the packed lower triangle of a symmetric matrix into full storage
void fillDense(const DoubleSymmMatrix &mat, std::vector<double> &dense) {
  const unsigned int N = mat.numRows();
  const double *data = mat.getData();
  for (unsigned int i = 0; i < N; i++) {
    const unsigned int id = i * (i + 1) / 2;
    for (unsigned int j = 0; j <= i; j++) {
      dense[static_cast<std::size_t>(i) * N + j] = data[id + j];
      dense[static_cast<std::size_t>(j) * N + i] = data[id + j];
    }
  }
}

// y = A*x for a symmetric matrix A in full storage.
// This produces exactly the same result as RDNumeric::multiply() on the packed
// matrix: because A is symmetric, y[i] = sum_j A[j][i] * x[j], and every y[i]
// receives its terms in the same order (j = 0, 1, ...). Sweeping the matrix
// row by row means the inner loop reads memory contiguously and has no serial
// dependency, so it is much faster.
void denseSymmMultiply(const std::vector<double> &dense, const DoubleVector &x,
                       DoubleVector &y) {
  const unsigned int N = x.size();
  const double *xData = x.getData();
  double *yData = y.getData();
  std::fill(yData, yData + N, 0.0);
  for (unsigned int j = 0; j < N; j++) {
    const double *row = dense.data() + static_cast<std::size_t>(j) * N;
    const double xj = xData[j];
    for (unsigned int i = 0; i < N; i++) {
      yData[i] += (row[i] * xj);
    }
  }
}
}  // namespace

bool powerEigenSolver(unsigned int numEig, DoubleSymmMatrix &mat,
                      DoubleVector &eigenValues, DoubleMatrix *eigenVectors,
                      int seed) {
  const unsigned int MAX_ITERATIONS = 1000;
  const double TOLERANCE = 0.001;
  const double HUGE_EIGVAL = 1.0e10;
  const double TINY_EIGVAL = 1.0e-10;

  // first check all the sizes
  unsigned int N = mat.numRows();
  CHECK_INVARIANT(eigenValues.size() >= numEig, "");
  CHECK_INVARIANT(numEig <= N, "");
  if (eigenVectors) {
    unsigned int evRows, evCols;
    evRows = eigenVectors->numRows();
    evCols = eigenVectors->numCols();
    CHECK_INVARIANT(evCols >= N, "");
    CHECK_INVARIANT(evRows >= numEig, "");
  }

  unsigned int ei;
  double eigVal, prevVal;
  bool converged = false;
  unsigned int i, j, id, iter, evalId;

  DoubleVector v(N), z(N);
  // A dense (full storage) copy of mat. Multiplying with it is much faster
  // than multiplying with the packed symmetric storage, see denseSymmMultiply.
  std::vector<double> dense(static_cast<std::size_t>(N) * N);
  if (seed <= 0) {
    seed = clock();
  }
  for (ei = 0; ei < numEig; ei++) {
    eigVal = -HUGE_EIGVAL;
    seed += ei;
    v.setToRandom(seed);
    fillDense(mat, dense);

    converged = false;
    for (iter = 0; iter < MAX_ITERATIONS; iter++) {
      // z = mat*v
      denseSymmMultiply(dense, v, z);
      prevVal = eigVal;
      evalId = z.largestAbsValId();
      eigVal = z.getVal(evalId);

      if (fabs(eigVal) < TINY_EIGVAL) {
        break;
      }

      // compute the next estimate for the eigen vector
      v.assign(z);
      v /= eigVal;
      if (fabs(eigVal - prevVal) < TOLERANCE) {
        converged = true;
        break;
      }
    }
    if (!converged) {
      break;
    }
    v.normalize();

    // save this is a eigen vector and value
    // directly access the data instead of setVal so that we save time
    double *vdata = v.getData();
    if (eigenVectors) {
      id = ei * eigenVectors->numCols();
      double *eigVecData = eigenVectors->getData();
      for (i = 0; i < N; i++) {
        eigVecData[id + i] = vdata[i];
      }
    }
    eigenValues[ei] = eigVal;

    // now remove this eigen vector space out of the matrix
    double *matData = mat.getData();
    for (i = 0; i < N; i++) {
      id = i * (i + 1) / 2;
      for (j = 0; j < i + 1; j++) {
        matData[id + j] -= (eigVal * vdata[i] * vdata[j]);
      }
    }
  }
  return converged;
}
}  // namespace EigenSolvers
}  // namespace RDNumeric
