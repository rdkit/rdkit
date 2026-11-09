//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_BFGSKERNELS_H
#define RD_BFGSKERNELS_H

#include <RDGeneral/export.h>

//! The O(dim) and O(dim^2) array kernels used by BFGSOpt::minimize().
/*!
  These are compiled into the Optimizer library rather than living in the
  header so that the implementation (scalar or Highway, depending on
  the build) is selected in one place. All matrices are dense, row-major and
  dim x dim.
*/
namespace BFGSOpt::detail {

//! The dot products needed to decide on and apply the inverse Hessian update.
struct UpdateDots {
  double fac = 0.0;       //!< dGrad . xi
  double fae = 0.0;       //!< dGrad . hessDGrad
  double sumDGrad = 0.0;  //!< dGrad . dGrad
  double sumXi = 0.0;     //!< xi . xi
};

//! Sets xi = -grad and returns pos . pos
RDKIT_OPTIMIZER_EXPORT double initDirection(unsigned int dim,
                                            const double *grad, double *xi,
                                            const double *pos);

//! Sets hessDGrad = invHessian * dGrad and returns the update dot products
RDKIT_OPTIMIZER_EXPORT UpdateDots hessianTimesDGrad(unsigned int dim,
                                                    const double *invHessian,
                                                    const double *dGrad,
                                                    const double *xi,
                                                    double *hessDGrad);

//! Applies the BFGS update to the (symmetric) inverse Hessian
/*!
  invHessian += fac * xi xi^T - fad * hessDGrad hessDGrad^T
              + fae * dGrad dGrad^T
*/
RDKIT_OPTIMIZER_EXPORT void updateInverseHessian(
    unsigned int dim, double *invHessian, const double *xi,
    const double *hessDGrad, const double *dGrad, double fac, double fad,
    double fae);

//! Sets xi = -(invHessian * grad)
RDKIT_OPTIMIZER_EXPORT void newDirection(unsigned int dim,
                                         const double *invHessian,
                                         const double *grad, double *xi);

}  // namespace BFGSOpt::detail

#endif
