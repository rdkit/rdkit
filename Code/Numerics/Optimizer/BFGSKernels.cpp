//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Scalar implementation of the BFGS kernels. Builds with RDK_USE_HIGHWAY use
// BFGSKernelsHwy.cpp instead.
//
#include "BFGSKernels.h"

namespace BFGSOpt::detail {

double initDirection(unsigned int dim, const double *grad, double *xi,
                     const double *pos) {
  double sum = 0.0;
  for (unsigned int i = 0; i < dim; i++) {
    xi[i] = -grad[i];
    sum += pos[i] * pos[i];
  }
  return sum;
}

UpdateDots hessianTimesDGrad(unsigned int dim, const double *invHessian,
                             const double *dGrad, const double *xi,
                             double *hessDGrad) {
  UpdateDots dots;
  // fused matrix-vector multiply and dot-product accumulation.
  // Pointer arithmetic (++ivh, ++dgj) avoids repeated index computations
  // and helps the compiler generate efficient load sequences.
  for (unsigned int i = 0; i < dim; i++) {
    const double *ivh = &(invHessian[i * dim]);
    double &hdgradi = hessDGrad[i];
    const double *dgj = dGrad;
    hdgradi = 0.0;
    for (unsigned int j = 0; j < dim; ++j, ++ivh, ++dgj) {
      hdgradi += *ivh * *dgj;
    }
    dots.fac += dGrad[i] * xi[i];
    dots.fae += dGrad[i] * hessDGrad[i];
    dots.sumDGrad += dGrad[i] * dGrad[i];
    dots.sumXi += xi[i] * xi[i];
  }
  return dots;
}

void updateInverseHessian(unsigned int dim, double *invHessian,
                          const double *xi, const double *hessDGrad,
                          const double *dGrad, double fac, double fad,
                          double fae) {
  // upper-triangle-only update (j >= i) followed by explicit
  // symmetrisation. This halves the number of Hessian writes at the cost of
  // one additional pass over a row to mirror elements.
  for (unsigned int i = 0; i < dim; i++) {
    unsigned int itab = i * dim;
    double pxi = fac * xi[i], hdgi = fad * hessDGrad[i], dgi = fae * dGrad[i];
    const double *pxj = &(xi[i]), *hdgj = &(hessDGrad[i]), *dgj = &(dGrad[i]);
    for (unsigned int j = i; j < dim; ++j, ++pxj, ++hdgj, ++dgj) {
      invHessian[itab + j] += pxi * *pxj - hdgi * *hdgj + dgi * *dgj;
      invHessian[j * dim + i] = invHessian[itab + j];
    }
  }
}

void newDirection(unsigned int dim, const double *invHessian,
                  const double *grad, double *xi) {
  for (unsigned int i = 0; i < dim; i++) {
    unsigned int itab = i * dim;
    xi[i] = 0.0;
    double &pxi = xi[i];
    const double *ivh = &(invHessian[itab]);
    const double *gj = grad;
    for (unsigned int j = 0; j < dim; ++j, ++ivh, ++gj) {
      pxi -= *ivh * *gj;
    }
  }
}

}  // namespace BFGSOpt::detail
