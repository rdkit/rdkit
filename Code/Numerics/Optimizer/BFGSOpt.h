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
#include <math.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/Trajectory/Snapshot.h>
#include <cstring>
#include <vector>
#include <algorithm>

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
const double MAXSTEP = 100.0;  //!< Default maximim step size in the minimizer

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
  for (unsigned int i = 0; i < dim; i++) sum += dir[i] * dir[i];
  sum = sqrt(sum);

  // rescale if we're trying to move too far:
  if (sum > maxStep) {
    for (unsigned int i = 0; i < dim; i++) dir[i] *= maxStep / sum;
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
    if (temp > test) test = temp;
  }

  lambdaMin = MOVETOL / test;
  lambda = 1.0;
  unsigned int it = 0;
  while (it < MAX_ITER_LINEAR_SEARCH) {
    // std::cerr << "\t" << it<<" : "<<lambda << " " << lambdaMin << std::endl;
    if (lambda < lambdaMin) {
      // the position change is too small
      resCode = 1;
      break;
    }
    for (unsigned int i = 0; i < dim; i++) {
      newPt[i] = oldPt[i] + lambda * dir[i];
    }
    newVal = func(newPt);

    if (newVal - oldVal <= FUNCTOL * lambda * slope) {
      // we're converged on the function:
      resCode = 0;
      return;
    }
    // if we made it this far, we need to backtrack:
    if (it == 0) {
      // it's the first step:
      tmpLambda = -slope / (2.0 * (newVal - oldVal - slope));
    } else {
      double rhs1 = newVal - oldVal - lambda * slope;
      double rhs2 = val2 - oldVal - lambda2 * slope;
      double a = (rhs1 / (lambda * lambda) - rhs2 / (lambda2 * lambda2)) /
                 (lambda - lambda2);
      double b = (-lambda2 * rhs1 / (lambda * lambda) +
                  lambda * rhs2 / (lambda2 * lambda2)) /
                 (lambda - lambda2);
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
    val2 = newVal;
    lambda = std::max(tmpLambda, 0.1 * lambda);
    ++it;
  }
  // nothing was done
  // std::cerr<<"  RETURN AT END: "<<it<<" "<<resCode<<std::endl;
  for (unsigned int i = 0; i < dim; i++) {
    newPt[i] = oldPt[i];
  }
}

#define CLEANUP()        \
  {                      \
    delete[] grad;       \
    delete[] dGrad;      \
    delete[] hessDGrad;  \
    delete[] newPos;     \
    delete[] xi;         \
    delete[] invHessian; \
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
   \param snapshotVect     pointer to a std::vector<Snapshot> object that will receive
                           the coordinates and energies every snapshotFreq steps;
                           defaults to NULL (no snapshots stored)

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

  double sum, maxStep, fp;

  double *grad, *dGrad, *hessDGrad;
  double *newPos, *xi;
  double *invHessian;

  grad = new double[dim];
  dGrad = new double[dim];
  hessDGrad = new double[dim];
  newPos = new double[dim];
  xi = new double[dim];
  invHessian = new double[dim * dim];
  snapshotFreq = std::min(snapshotFreq, maxIts);

  // evaluate the function and gradient in our current position:
  fp = func(pos);
  gradFunc(pos, grad);

  sum = 0.0;
  memset(invHessian, 0, dim * dim * sizeof(double));
  for (unsigned int i = 0; i < dim; i++) {
    unsigned int itab = i * dim;
    // initialize the inverse hessian to be identity:
    invHessian[itab + i] = 1.0;
    // the first line dir is -grad:
    xi[i] = -grad[i];
    sum += pos[i] * pos[i];
  }
  // pick a max step size:
  maxStep = MAXSTEP * std::max(sqrt(sum), static_cast<double>(dim));

  for (unsigned int iter = 1; iter <= maxIts; iter++) {
    numIters = iter;
    int status;

    // do the line search:
    linearSearch(dim, pos, fp, grad, xi, newPos, funcVal, func, maxStep,
                 status);
    CHECK_INVARIANT(status >= 0, "bad direction in linearSearch");

    // save the function value for the next search:
    fp = funcVal;

    // set the direction of this line and save the gradient:
    double test = 0.0;
    for (unsigned int i = 0; i < dim; i++) {
      xi[i] = newPos[i] - pos[i];
      pos[i] = newPos[i];
      double temp = fabs(xi[i]) / std::max(fabs(pos[i]), 1.0);
      if (temp > test) test = temp;
      dGrad[i] = grad[i];
    }
    // std::cerr<<"      iter: "<<iter<<" "<<fp<<" "<<test<<"
    // "<<TOLX<<std::endl;
    if (test < TOLX) {
      if (snapshotVect && snapshotFreq) {
        RDKit::Snapshot s(boost::shared_array<double>(newPos), fp);
        snapshotVect->push_back(s);
        newPos = NULL;
      }
      CLEANUP();
      return 0;
    }

    // update the gradient:
    double gradScale = gradFunc(pos, grad);

    // is the gradient converged?
    test = 0.0;
    double term = std::max(funcVal * gradScale, 1.0);
    for (unsigned int i = 0; i < dim; i++) {
      double temp = fabs(grad[i]) * std::max(fabs(pos[i]), 1.0);
      test = std::max(test, temp);
      dGrad[i] = grad[i] - dGrad[i];
    }
    test /= term;
    // std::cerr<<"              "<<gradScale<<" "<<test<<"
    // "<<gradTol<<std::endl;
    if (test < gradTol) {
      if (snapshotVect && snapshotFreq) {
        RDKit::Snapshot s(boost::shared_array<double>(newPos), fp);
        snapshotVect->push_back(s);
        newPos = NULL;
      }
      CLEANUP();
      return 0;
    }

    // for(unsigned int i=0;i<dim;i++){
    // figure out how much the gradient changed:
    //}

    // compute hessian*dGrad:
    double fac = 0, fae = 0, sumDGrad = 0, sumXi = 0;
    for (unsigned int i = 0; i < dim; i++) {
#if 0
      unsigned int itab = i * dim;
      hessDGrad[i] = 0.0;
      for (unsigned int j = 0; j < dim; j++) {
        hessDGrad[i] += invHessian[itab + j] * dGrad[j];
      }

#else
      double *ivh = &(invHessian[i * dim]);
      double &hdgradi = hessDGrad[i];
      double *dgj = dGrad;
      hdgradi = 0.0;
      for (unsigned int j = 0; j < dim; ++j, ++ivh, ++dgj) {
        hdgradi += *ivh * *dgj;
      }
#endif
      fac += dGrad[i] * xi[i];
      fae += dGrad[i] * hessDGrad[i];
      sumDGrad += dGrad[i] * dGrad[i];
      sumXi += xi[i] * xi[i];
    }
    if (fac > sqrt(EPS * sumDGrad * sumXi)) {
      fac = 1.0 / fac;
      double fad = 1.0 / fae;
      for (unsigned int i = 0; i < dim; i++) {
        dGrad[i] = fac * xi[i] - fad * hessDGrad[i];
      }

      for (unsigned int i = 0; i < dim; i++) {
        unsigned int itab = i * dim;
#if 0
        for (unsigned int j = i; j < dim; j++) {
          invHessian[itab + j] += fac * xi[i] * xi[j] -
                                  fad * hessDGrad[i] * hessDGrad[j] +
                                  fae * dGrad[i] * dGrad[j];
          invHessian[j * dim + i] = invHessian[itab + j];
        }
#else
        double pxi = fac * xi[i], hdgi = fad * hessDGrad[i],
               dgi = fae * dGrad[i];
        double *pxj = &(xi[i]), *hdgj = &(hessDGrad[i]), *dgj = &(dGrad[i]);
        for (unsigned int j = i; j < dim; ++j, ++pxj, ++hdgj, ++dgj) {
          invHessian[itab + j] += pxi * *pxj - hdgi * *hdgj + dgi * *dgj;
          invHessian[j * dim + i] = invHessian[itab + j];
        }
#endif
      }
    }
    // generate the next direction to move:
    for (unsigned int i = 0; i < dim; i++) {
      unsigned int itab = i * dim;
      xi[i] = 0.0;
#if 0
      for (unsigned int j = 0; j < dim; j++) {
        xi[i] -= invHessian[itab + j] * grad[j];
      }
#else
      double &pxi = xi[i], *ivh = &(invHessian[itab]), *gj = grad;
      for (unsigned int j = 0; j < dim; ++j, ++ivh, ++gj) {
        pxi -= *ivh * *gj;
      }
#endif
    }
    if (snapshotVect && snapshotFreq && !(iter % snapshotFreq)) {
      RDKit::Snapshot s(boost::shared_array<double>(newPos), fp);
      snapshotVect->push_back(s);
      newPos = new double[dim];
    }
  }
  CLEANUP();
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
  return minimize(dim, pos, gradTol, numIters, funcVal, func,
             gradFunc, 0, NULL, funcTol, maxIts);
}

}
