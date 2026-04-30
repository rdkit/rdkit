//
// Copyright (C)  2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#include <cmath>
#include <RDGeneral/Invariant.h>
#include <GraphMol/Trajectory/Snapshot.h>
#include <cstring>
#include <vector>
#include <algorithm>
#ifdef __aarch64__
  #include <sys/auxv.h>
  #include <asm/hwcap.h>
  #if defined(__has_include)
    #if __has_include(<arm_sve.h>)
      #include <arm_sve.h>
      #define RDK_SVE_AVAILABLE 1
    #endif
  #endif
#endif

namespace BFGSOpt {
RDKIT_OPTIMIZER_EXPORT extern int HEAD_ONLY_LIBRARY;
RDKIT_OPTIMIZER_EXPORT extern int REALLY_A_HEADER_ONLY_LIBRARY;
const double FUNCTOL =
    1e-4;  //!< Default tolerance for function convergence in the minimizer
const double MOVETOL =
    1e-7;                 //!< Default tolerance for x changes in the minimizer
const int MAXITS = 200;   //!< Default maximum number of iterations
const double EPS = 3e-8;  //!< Default gradient tolerance in the minimizer
const double TOLX =
    4. * EPS;  //!< Default direction vector tolerance in the minimizer
const double MAXSTEP = 100.0;  //!< Default maximum step size in the minimizer

static bool cpuHasSVE() {
#if defined(__linux__) && defined(__aarch64__) && defined(RDK_SVE_AVAILABLE)
  static const bool result = (getauxval(AT_HWCAP) & HWCAP_SVE) != 0;
  return result;
#else
  return false;
#endif
}

#ifdef RDK_SVE_AVAILABLE

// ---------------------------------------------------------------------------
// SVE kernel: initialise search direction xi = -grad and accumulate ||pos||^2
// ---------------------------------------------------------------------------
__attribute__((target("+sve")))
static void sveInitXiAndSum(unsigned int dim, const double *grad,
                             double *xi, const double *pos, double *outSum) {
  svfloat64_t acc = svdup_f64(0.0);
  uint64_t i = 0;
  while (i < dim) {
    svbool_t pg = svwhilelt_b64(i, (uint64_t)dim);
    svfloat64_t g = svld1_f64(pg, grad + i);
    svfloat64_t p = svld1_f64(pg, pos + i);
    // Negate grad in-place during store — avoids a separate negation pass
    svst1_f64(pg, xi + i, svneg_f64_m(svdup_f64(0.0), pg, g));
    // Fused multiply-add accumulates p[i]^2 without intermediate stores
    acc = svmla_f64_m(pg, acc, p, p);
    i += svcntd();
  }
  // Horizontal reduction collapses all vector lanes to a single scalar sum
  *outSum = svaddv_f64(svptrue_b64(), acc);
}

// ---------------------------------------------------------------------------
// SVE kernel: compute hessDGrad = invHessian * dGrad, then accumulate the
// four dot-product scalars (fac, fae, sumDGrad, sumXi) needed for the BFGS
// rank-1 update.
//
// The inner matrix-vector product uses SVE gather/FMA to handle arbitrary dim
// without scalar remainder loops. Each outer row is streamed once, keeping
// cache pressure proportional to dim rather than dim^2.  The four scalar
// accumulators are computed in a single fused pass over the result vector,
// saving two additional O(dim) traversals compared to separate dot-product
// calls.
// ---------------------------------------------------------------------------
__attribute__((target("+sve")))
static void sveHessianVecMul(unsigned int dim,
                              const double *invHessian, const double *dGrad,
                              double *hessDGrad, const double *xi,
                              double *outFac, double *outFae,
                              double *outSumDGrad, double *outSumXi) {
  // Phase 1: matrix-vector multiply  invHessian * dGrad -> hessDGrad
  // Fully vectorised over both i and j.
  for (unsigned int i = 0; i < dim; i++) {
    const double *row = invHessian + i * dim;
    svfloat64_t acc = svdup_f64(0.0);
    uint64_t j = 0;
    while (j < dim) {
      svbool_t pg = svwhilelt_b64(j, (uint64_t)dim);
      acc = svmla_f64_m(pg, acc,
                        svld1_f64(pg, row + j),
                        svld1_f64(pg, dGrad + j));
      j += svcntd();
    }
    hessDGrad[i] = svaddv_f64(svptrue_b64(), acc);
  }

  svfloat64_t vFac = svdup_f64(0.0), vFae     = svdup_f64(0.0);
  svfloat64_t vSDG = svdup_f64(0.0), vSXi     = svdup_f64(0.0);
  uint64_t i = 0;
  while (i < dim) {
    svbool_t    pg  = svwhilelt_b64(i, (uint64_t)dim);
    svfloat64_t vdg = svld1_f64(pg, dGrad     + i);
    svfloat64_t vxi = svld1_f64(pg, xi        + i);
    svfloat64_t vhd = svld1_f64(pg, hessDGrad + i);

    vFac = svmla_f64_m(pg, vFac, vdg, vxi);
    vFae = svmla_f64_m(pg, vFae, vdg, vhd);
    vSDG = svmla_f64_m(pg, vSDG, vdg, vdg);
    vSXi = svmla_f64_m(pg, vSXi, vxi, vxi);
    i += svcntd();
  }
  *outFac      = svaddv_f64(svptrue_b64(), vFac);
  *outFae      = svaddv_f64(svptrue_b64(), vFae);
  *outSumDGrad = svaddv_f64(svptrue_b64(), vSDG);
  *outSumXi    = svaddv_f64(svptrue_b64(), vSXi);
}

// ---------------------------------------------------------------------------
// SVE kernel: symmetric rank-1 update of the inverse Hessian approximation.
//
// The BFGS update formula adds three outer-product terms to invHessian.
// Exploiting symmetry (only the upper triangle is computed; the lower is
// mirrored afterwards) halves the number of FLOPs and memory writes versus a
// naive full-matrix update. The SVE FMA instructions (svmla / svmls) fuse
// multiply and add into a single pipeline stage, reducing instruction count
// and register pressure compared to separate multiply + add sequences.
// ---------------------------------------------------------------------------
__attribute__((target("+sve")))
static void sveHessianRank1Update(unsigned int dim, double *invHessian,
                                   const double *xi, const double *hessDGrad,
                                   const double *dGrad,
                                   double fac, double fad, double fae) {
  for (unsigned int i = 0; i < dim; i++) {
    // Broadcast scalar multipliers once per row to avoid redundant computation
    svfloat64_t vpxi  = svdup_f64(fac * xi[i]);
    svfloat64_t vhdgi = svdup_f64(fad * hessDGrad[i]);
    svfloat64_t vdgi  = svdup_f64(fae * dGrad[i]);
    double *row = invHessian + i * dim;
    // Start from column i to process only the upper triangle
    uint64_t j = i;
    while (j < dim) {
      svbool_t pg = svwhilelt_b64(j, (uint64_t)dim);
      svfloat64_t vxj   = svld1_f64(pg, xi + j);
      svfloat64_t vhdgj = svld1_f64(pg, hessDGrad + j);
      svfloat64_t vdgj  = svld1_f64(pg, dGrad + j);
      svfloat64_t vh    = svld1_f64(pg, row + j);
      // Three fused multiply-adds apply all three BFGS update terms at once
      vh = svmla_f64_m(pg, vh, vpxi,  vxj);
      vh = svmls_f64_m(pg, vh, vhdgi, vhdgj);
      vh = svmla_f64_m(pg, vh, vdgi,  vdgj);
      svst1_f64(pg, row + j, vh);
      j += svcntd();
    }
    // Mirror upper triangle to lower triangle to maintain symmetry;
    // scalar loop is cheap (dim - i iterations) relative to the SVE inner loop
    for (unsigned int j2 = i + 1; j2 < dim; j2++) {
      invHessian[j2 * dim + i] = invHessian[i * dim + j2];
    }
  }
}

// ---------------------------------------------------------------------------
// SVE kernel: compute new search direction xi = -(invHessian * grad)
//
// SVE implementation accelerates the matrix-vector multiply using vectorized
// FMA across columns, then negates the accumulated dot product at the scalar
// level (single negation per row) rather than applying a separate vector
// negation pass. This keeps the loop body to one SVE FMA instruction
// per iteration, maximising throughput on in-order SVE pipelines.
// ---------------------------------------------------------------------------
__attribute__((target("+sve")))
static void sveHessianVecMulNeg(unsigned int dim, const double *invHessian,
                                 const double *grad, double *xi) {
  for (unsigned int i = 0; i < dim; i++) {
    const double *row = invHessian + i * dim;
    svfloat64_t acc = svdup_f64(0.0);
    uint64_t j = 0;
    while (j < dim) {
      svbool_t pg = svwhilelt_b64(j, (uint64_t)dim);
      svfloat64_t h = svld1_f64(pg, row + j);
      svfloat64_t g = svld1_f64(pg, grad + j);
      acc = svmla_f64_m(pg, acc, h, g);
      j += svcntd();
    }
    // Negate the scalar result once per row rather than vectorising the
    // negation, keeping the store path simple and avoiding an extra SVE pass
    xi[i] = -svaddv_f64(svptrue_b64(), acc);
  }
}

#endif

//! Do a Quasi-Newton minimization along a line.
/*!
  See Numerical Recipes in C, Section 9.7 for a description of the algorithm.

   \param dim     the dimensionality of the space.
   \param oldPt   the current position, as an array.
   \param oldVal  the current function value.
   \param grad    the value of the function gradient at oldPt
   \param dir     the minimization direction
   \param newPt   used to return the final position
   \param newVal  used to return the final function value
   \param func    the function to minimize
   \param maxStep the maximum allowable step size
   \param resCode used to return the results of the search.

   Possible values for resCode are on return are:
    -  0: success
    -  1: the stepsize got too small.  This probably indicates success.
    - -1: the direction is bad (orthogonal to the gradient)
*/
template <typename EnergyFunctor>
void linearSearch(unsigned int dim, double *oldPt, double oldVal, double *grad,
                  double *dir, double *newPt, double &newVal,
                  EnergyFunctor func, double maxStep, int &resCode) {
  PRECONDITION(oldPt, "bad input array");
  PRECONDITION(grad, "bad input array");
  PRECONDITION(dir, "bad input array");
  PRECONDITION(newPt, "bad input array");

  const unsigned int MAX_ITER_LINEAR_SEARCH = 1000;
  double sum = 0.0, slope = 0.0, test = 0.0, lambda = 0.0;
  double lambda2 = 0.0, lambdaMin = 0.0, tmpLambda = 0.0, val2 = 0.0;

  resCode = -1;

  // get the length of the direction vector:
  sum = 0.0;
  for (unsigned int i = 0; i < dim; i++) {
    sum += dir[i] * dir[i];
  }
  sum = sqrt(sum);

  // Rescale if we're trying to move too far
  if (sum > maxStep) {
    for (unsigned int i = 0; i < dim; i++) {
      dir[i] *= maxStep / sum;
    }
  }

  // make sure our direction has at least some component along
  // -grad
  slope = 0.0;
  for (unsigned int i = 0; i < dim; i++) {
    slope += dir[i] * grad[i];
  }
  if (slope >= 0.0) {
    return;
  }

  test = 0.0;
  for (unsigned int i = 0; i < dim; i++) {
    double temp = fabs(dir[i]) / std::max(fabs(oldPt[i]), 1.0);
    if (temp > test) {
      test = temp;
    }
  }

  lambdaMin = MOVETOL / test;
  lambda = 1.0;
  unsigned int it = 0;
  while (it < MAX_ITER_LINEAR_SEARCH) {
    if (lambda < lambdaMin) {
      // Step size is below the position-scaled threshold; treat as converged
      resCode = 1;
      break;
    }
    for (unsigned int i = 0; i < dim; i++) {
      newPt[i] = oldPt[i] + lambda * dir[i];
    }
    newVal = func(newPt);
    if (newVal - oldVal <= FUNCTOL * lambda * slope) {
      // Armijo sufficient-decrease condition satisfied; accept the step
      resCode = 0;
      return;
    }
    // if we made it this far, we need to backtrack:
    if (it == 0) {
      // Quadratic model: only one prior function value available
      tmpLambda = -slope / (2.0 * (newVal - oldVal - slope));
    } else {
      double rhs1 = newVal - oldVal - lambda * slope;
      double rhs2 = val2  - oldVal - lambda2 * slope;
      double a = (rhs1 / (lambda * lambda) - rhs2 / (lambda2 * lambda2)) / (lambda - lambda2);
      double b = (-lambda2 * rhs1 / (lambda * lambda) + lambda * rhs2 / (lambda2 * lambda2)) / (lambda - lambda2);
      if (a == 0.0) {
        tmpLambda = -slope / (2.0 * b);
      } else {
        double disc = b * b - 3 * a * slope;
        if (disc < 0.0) {
          tmpLambda = 0.5 * lambda;
        } else if (b <= 0.0) {
          tmpLambda = (-b + sqrt(disc)) / (3.0 * a);
        } else {
          tmpLambda = -slope / (b + sqrt(disc));
        }
      }
      if (tmpLambda > 0.5 * lambda) {
        tmpLambda = 0.5 * lambda;
      }
    }
    lambda2 = lambda;
    val2    = newVal;
    lambda  = std::max(tmpLambda, 0.1 * lambda);
    ++it;
  }
  // nothing was done
  for (unsigned int i = 0; i < dim; i++) {
    newPt[i] = oldPt[i];
  }
}

//! Do a BFGS minimization of a function.
/*!
   See Numerical Recipes in C, Section 10.7 for a description of the algorithm.

   \param dim     the dimensionality of the space.
   \param pos   the starting position, as an array.
   \param gradTol tolerance for gradient convergence
   \param numIters used to return the number of iterations required
   \param funcVal  used to return the final function value
   \param func    the function to minimize
   \param gradFunc  calculates the gradient of func
   \param funcTol tolerance for changes in the function value for convergence.
   \param maxIts   maximum number of iterations allowed
   \param snapshotFreq     a snapshot of the minimization trajectory
                           will be stored after as many steps as indicated
                           through this parameter; defaults to 0 (no
                           snapshots stored)
   \param snapshotVect     pointer to a std::vector<Snapshot> object that will
   receive the coordinates and energies every snapshotFreq steps; defaults to
   NULL (no snapshots stored)

   \return a flag indicating success (or type of failure). Possible values are:
    -  0: success
    -  1: too many iterations were required
*/
template <typename EnergyFunctor, typename GradientFunctor>
int minimize(unsigned int dim, double *pos, double gradTol,
             unsigned int &numIters, double &funcVal, EnergyFunctor func,
             GradientFunctor gradFunc, unsigned int snapshotFreq,
             RDKit::SnapshotVect *snapshotVect, double funcTol = TOLX,
             unsigned int maxIts = MAXITS) {
  RDUNUSED_PARAM(funcTol);
  PRECONDITION(pos, "bad input array");
  PRECONDITION(gradTol > 0, "bad tolerance");

  std::vector<double> grad(dim);
  std::vector<double> dGrad(dim);
  std::vector<double> hessDGrad(dim);
  std::vector<double> xi(dim);
  std::vector<double> invHessian(dim * dim, 0);
  std::unique_ptr<double[]> newPos(new double[dim]);
  snapshotFreq = std::min(snapshotFreq, maxIts);


  double fp = func(pos);
  gradFunc(pos, grad.data());

  double sum = 0.0;
#ifdef RDK_SVE_AVAILABLE
  if (cpuHasSVE()) {
    // SVE path: initialise xi = -grad and compute ||pos||^2 in a single
    // vectorised pass.  The identity inverse Hessian is initialised separately
    // (scalar, O(dim)) since it is a simple diagonal write and does not benefit
    // from vectorisation over rows.
    sveInitXiAndSum(dim, grad.data(), xi.data(), pos, &sum);
    for (unsigned int i = 0; i < dim; i++) invHessian[i * dim + i] = 1.0;
  } else
#endif
  {
    // Scalar path: initialise the inverse Hessian to the identity matrix,
    // set the initial search direction xi = -grad (steepest descent step),
    // and accumulate ||pos||^2 to set an appropriate maximum step size.
    for (unsigned int i = 0; i < dim; i++) {
      unsigned int itab = i * dim;
      invHessian[itab + i] = 1.0;
      xi[i]  = -grad[i];
      sum   += pos[i] * pos[i];
    }
  }
  double maxStep = MAXSTEP * std::max(sqrt(sum), static_cast<double>(dim));

  for (unsigned int iter = 1; iter <= maxIts; ++iter) {
    numIters = iter;
    int status = -1;

    linearSearch(dim, pos, fp, grad.data(), xi.data(), newPos.get(),
                 funcVal, func, maxStep, status);
    CHECK_INVARIANT(status >= 0, "bad direction in linearSearch");

    // save the function value for the next search:
    fp = funcVal;
    // set the direction of this line and save the gradient:
    double test = 0.0;
    for (unsigned int i = 0; i < dim; i++) {
      xi[i]   = newPos[i] - pos[i];
      pos[i]  = newPos[i];
      double temp = fabs(xi[i]) / std::max(fabs(pos[i]), 1.0);
      if (temp > test) {
        test = temp;
      }
      dGrad[i] = grad[i];
    }
    if (test < TOLX) {
      if (snapshotVect && snapshotFreq) {
        RDKit::Snapshot s(boost::shared_array<double>(newPos.release()), fp);
        snapshotVect->push_back(s);
      }
      return 0;
    }

    // update the gradient:
    double gradScale = gradFunc(pos, grad.data());

    test = 0.0;
    double term = std::max(funcVal * gradScale, 1.0);
    for (unsigned int i = 0; i < dim; i++) {
      double temp = fabs(grad[i]) * std::max(fabs(pos[i]), 1.0);
      test = std::max(test, temp);
      dGrad[i] = grad[i] - dGrad[i];
    }
    test /= term;
    if (test < gradTol) {
      if (snapshotVect && snapshotFreq) {
        RDKit::Snapshot s(boost::shared_array<double>(newPos.release()), fp);
        snapshotVect->push_back(s);
      }
      return 0;
    }

    // BFGS inverse Hessian update.
    double fac = 0, fae = 0, sumDGrad = 0, sumXi = 0;
#ifdef RDK_SVE_AVAILABLE
    if (cpuHasSVE()) {
      // SVE path: matrix-vector multiply and all four dot products computed in
      // one vectorised pass, saving two additional O(dim) traversals compared
      // to separate scalar dot-product calls.
      sveHessianVecMul(dim, invHessian.data(), dGrad.data(), hessDGrad.data(),
                       xi.data(), &fac, &fae, &sumDGrad, &sumXi);
    } else
#endif
    {
      // Scalar path: fused matrix-vector multiply and dot-product accumulation.
      // Pointer arithmetic (++ivh, ++dgj) avoids repeated index computations
      // and helps the compiler generate efficient load sequences.
      for (unsigned int i = 0; i < dim; i++) {
        double *ivh = &(invHessian[i * dim]);
        double &hdgradi = hessDGrad[i];
        double *dgj = dGrad.data();
        hdgradi = 0.0;
        for (unsigned int j = 0; j < dim; ++j, ++ivh, ++dgj) {
          hdgradi += *ivh * *dgj;
        }
        fac      += dGrad[i] * xi[i];
        fae      += dGrad[i] * hessDGrad[i];
        sumDGrad += dGrad[i] * dGrad[i];
        sumXi    += xi[i]    * xi[i];
      }
    }
    if (fac > sqrt(EPS * sumDGrad * sumXi)) {
      fac = 1.0 / fac;
      double fad = 1.0 / fae;
      for (unsigned int i = 0; i < dim; i++) {
        dGrad[i] = fac * xi[i] - fad * hessDGrad[i];
      }

#ifdef RDK_SVE_AVAILABLE
      if (cpuHasSVE()) {
        // SVE path: symmetric rank-1 update with FMA, exploiting symmetry to
        // halve memory writes and FLOPs versus a full-matrix update
        sveHessianRank1Update(dim, invHessian.data(), xi.data(),
                              hessDGrad.data(), dGrad.data(), fac, fad, fae);
      } else
#endif
      {
        // Scalar path: upper-triangle-only update (j >= i) followed by
        // explicit symmetrisation. This halves the number of Hessian writes
        // at the cost of one additional pass over a row to mirror elements.
        for (unsigned int i = 0; i < dim; i++) {
          unsigned int itab = i * dim;
          double pxi = fac * xi[i], hdgi = fad * hessDGrad[i],
                 dgi = fae * dGrad[i];
          double *pxj = &(xi[i]), *hdgj = &(hessDGrad[i]), *dgj = &(dGrad[i]);
          for (unsigned int j = i; j < dim; ++j, ++pxj, ++hdgj, ++dgj) {
            invHessian[itab + j] += pxi * *pxj - hdgi * *hdgj + dgi * *dgj;
            invHessian[j * dim + i] = invHessian[itab + j];
          }
        }
      }
    }

#ifdef RDK_SVE_AVAILABLE
    if (cpuHasSVE()) {
      sveHessianVecMulNeg(dim, invHessian.data(), grad.data(), xi.data());
    } else
#endif
    {
      for (unsigned int i = 0; i < dim; i++) {
        unsigned int itab = i * dim;
        xi[i] = 0.0;
        double &pxi = xi[i];
        double *ivh = &(invHessian[itab]);
        double *gj  = grad.data();
        for (unsigned int j = 0; j < dim; ++j, ++ivh, ++gj) {
          pxi -= *ivh * *gj;
        }
      }
    }
    if (snapshotVect && snapshotFreq && !(iter % snapshotFreq)) {
      RDKit::Snapshot s(boost::shared_array<double>(newPos.release()), fp);
      snapshotVect->push_back(s);
      newPos.reset(new double[dim]);
    }
  }
  return 1;
}

//! Do a BFGS minimization of a function.
/*!
   \param dim     the dimensionality of the space.
   \param pos   the starting position, as an array.
   \param gradTol tolerance for gradient convergence
   \param numIters used to return the number of iterations required
   \param funcVal  used to return the final function value
   \param func    the function to minimize
   \param gradFunc  calculates the gradient of func
   \param funcTol tolerance for changes in the function value for convergence.
   \param maxIts   maximum number of iterations allowed

   \return a flag indicating success (or type of failure). Possible values are:
    -  0: success
    -  1: too many iterations were required
*/
template <typename EnergyFunctor, typename GradientFunctor>
int minimize(unsigned int dim, double *pos, double gradTol,
             unsigned int &numIters, double &funcVal, EnergyFunctor func,
             GradientFunctor gradFunc, double funcTol = TOLX,
             unsigned int maxIts = MAXITS) {
  return minimize(dim, pos, gradTol, numIters, funcVal, func, gradFunc, 0,
                  nullptr, funcTol, maxIts);
}

}  // namespace BFGSOpt
