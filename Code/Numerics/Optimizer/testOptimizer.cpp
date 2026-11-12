//
// Copyright (C)  2004-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <RDGeneral/Invariant.h>
#include <catch2/catch_all.hpp>

#include "BFGSOpt.h"

#ifdef RDK_USE_HIGHWAY
#include <hwy/targets.h>
#endif

double circ_0_0(double *v) {
  double dx = v[0];
  double dy = v[1];

  return dx * dx + dy * dy;
}

double circ_0_0_grad(double *v, double *grad) {
  double dx = v[0];
  double dy = v[1];
  grad[0] = 2 * dx;
  grad[1] = 2 * dy;
  return 1.0;
}

double circ_1_0(double *v) {
  double dx = v[0] - 1;
  double dy = v[1];

  return dx * dx + dy * dy;
}

double circ_1_0_grad(double *v, double *grad) {
  double dx = v[0] - 1;
  double dy = v[1];
  grad[0] = 2 * dx;
  grad[1] = 2 * dy;
  return 1.0;
}

double func2(double *v) {
  double weight = .5;
  double dx = v[0] - 1;
  double dy = v[1];
  double term1 = dx * dx - dy * dy;
  double term2 = dx * dx + dy * dy;
  return term1 * term1 + weight * term2;
}

double grad2(double *v, double *grad) {
  double weight = .5;
  double dx = v[0] - 1;
  double dy = v[1];
  double term1 = dx * dx - dy * dy;
  grad[0] = 4 * dx * term1 + 2 * weight * dx;
  grad[1] = -4 * dy * term1 + 2 * weight * dy;
  return 1.0;
}

// A quadratic whose value is always negative on the convergence path:
//   f(x, y) = (x - 3)^2 + (y + 1)^2 - 100
// The minimum is f = -100 at (3, -1). Useful for exercising the
// gradient-convergence denominator when the energy is negative.
double circ_neg(double *v) {
  double dx = v[0] - 3.0;
  double dy = v[1] + 1.0;
  return dx * dx + dy * dy - 100.0;
}

double circ_neg_grad(double *v, double *grad) {
  double dx = v[0] - 3.0;
  double dy = v[1] + 1.0;
  grad[0] = 2 * dx;
  grad[1] = 2 * dy;
  return 1.0;
}

TEST_CASE("testLinearSearch") {
  int dim = 2;
  double oLoc[2], oVal;
  double grad[2], dir[2];
  double nLoc[2], nVal;
  int resCode;
  double (*func)(double *);
  double (*gradFunc)(double *, double *);

  func = circ_0_0;
  gradFunc = circ_0_0_grad;
  oLoc[0] = 0;
  oLoc[1] = 1.0;
  oVal = func(oLoc);
  REQUIRE_THAT(oVal, Catch::Matchers::WithinAbs(1.0, 1e-4));
  gradFunc(oLoc, grad);
  dir[0] = 0;
  dir[1] = -.5;

  BFGSOpt::linearSearch(dim, oLoc, oVal, grad, dir, nLoc, nVal, func, 0.5,
                        resCode);
  REQUIRE(resCode == 0);
  REQUIRE_THAT(nVal, Catch::Matchers::WithinAbs(0.25, 1e-4));
  REQUIRE_THAT(nLoc[0], Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(nLoc[1], Catch::Matchers::WithinAbs(0.5, 1e-4));

  oLoc[0] = 1.0;
  oLoc[1] = 1.0;
  oVal = func(oLoc);
  REQUIRE_THAT(oVal, Catch::Matchers::WithinAbs(2.0, 1e-4));
  gradFunc(oLoc, grad);
  dir[0] = -.5;
  dir[1] = -.5;

  BFGSOpt::linearSearch(dim, oLoc, oVal, grad, dir, nLoc, nVal, func, 1.0,
                        resCode);
  REQUIRE(resCode == 0);
  REQUIRE_THAT(nVal, Catch::Matchers::WithinAbs(0.5, 1e-4));
  REQUIRE_THAT(nLoc[0], Catch::Matchers::WithinAbs(0.5, 1e-4));
  REQUIRE_THAT(nLoc[1], Catch::Matchers::WithinAbs(0.5, 1e-4));

  func = circ_0_0;
  oLoc[0] = 0;
  oLoc[1] = 1.0;
  oVal = func(oLoc);
  REQUIRE_THAT(oVal, Catch::Matchers::WithinAbs(1.0, 1e-4));
  gradFunc(oLoc, grad);
  dir[0] = 0;
  dir[1] = -2;

  BFGSOpt::linearSearch(dim, oLoc, oVal, grad, dir, nLoc, nVal, func, 2.0,
                        resCode);
  REQUIRE(resCode == 0);
  REQUIRE_THAT(nVal, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(nLoc[0], Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(nLoc[1], Catch::Matchers::WithinAbs(0.0, 1e-4));
}

TEST_CASE("testBFGSOptimization") {
  unsigned int dim = 2;
  double oLoc[2], oVal;
  double nVal;
  unsigned int nIters;
  double (*func)(double *);
  double (*gradFunc)(double *, double *);

  func = circ_0_0;
  gradFunc = circ_0_0_grad;
  oLoc[0] = 0;
  oLoc[1] = 1.0;
  oVal = func(oLoc);
  REQUIRE_THAT(oVal, Catch::Matchers::WithinAbs(1.0, 1e-4));

  BFGSOpt::minimize(dim, oLoc, 1e-4, nIters, nVal, func, gradFunc);
  REQUIRE(nIters == 1);
  REQUIRE_THAT(nVal, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(oLoc[0], Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(oLoc[1], Catch::Matchers::WithinAbs(0.0, 1e-4));

  func = func2;
  gradFunc = grad2;
  oLoc[0] = 2.0;
  oLoc[1] = 0.5;
  BFGSOpt::minimize(dim, oLoc, 1e-4, nIters, nVal, func, gradFunc);
  REQUIRE_THAT(nVal, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(oLoc[0], Catch::Matchers::WithinAbs(1.0, 1e-3));
  REQUIRE_THAT(oLoc[1], Catch::Matchers::WithinAbs(0.0, 1e-3));

  oLoc[0] = 2.0;
  oLoc[1] = 0.5;
  BFGSOpt::minimize(dim, oLoc, 1e-4, nIters, nVal, func, gradFunc, 1e-8);
  REQUIRE_THAT(nVal, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(oLoc[0], Catch::Matchers::WithinAbs(1.0, 1e-4));
  REQUIRE_THAT(oLoc[1], Catch::Matchers::WithinAbs(0.0, 1e-4));
}

// Regression test for the gradient-convergence denominator: when the
// energy is negative the denominator used to be
// std::max(funcVal * gradScale, 1.0), which clamps to 1.0 for any
// funcVal < 0 and so over-tightens the test. With the fix
// (std::max(fabs(funcVal) * gradScale, 1.0)) this case must still
// converge to the analytic minimum.
TEST_CASE("testBFGSOptimizationNegativeEnergy") {
  unsigned int dim = 2;
  double oLoc[2] = {0.0, 0.0};
  double nVal = 0.0;
  unsigned int nIters = 0;
  BFGSOpt::minimize(dim, oLoc, 1e-4, nIters, nVal, circ_neg, circ_neg_grad);
  REQUIRE_THAT(nVal, Catch::Matchers::WithinAbs(-100.0, 1e-4));
  REQUIRE_THAT(oLoc[0], Catch::Matchers::WithinAbs(3.0, 1e-3));
  REQUIRE_THAT(oLoc[1], Catch::Matchers::WithinAbs(-1.0, 1e-3));
}

namespace {
// Compares the BFGS kernels against straightforward reference
// implementations. The dimensions cover the tail handling of every vector
// width the kernels might be compiled for.
void checkBFGSKernels() {
  for (unsigned int dim : {1u, 2u, 3u, 5u, 8u, 13u, 64u, 131u}) {
    INFO("dim = " << dim);
    std::vector<double> a(dim), b(dim), c(dim), sym(dim * dim);
    for (unsigned int i = 0; i < dim; ++i) {
      a[i] = std::sin(1.0 + i);
      b[i] = std::cos(0.5 * i) + 0.1;
      c[i] = 0.3 * i - 1.0;
      for (unsigned int j = 0; j <= i; ++j) {
        sym[i * dim + j] = sym[j * dim + i] = 1.0 / (1.0 + i + j);
      }
    }
    auto dot = [dim](const double *x, const double *y) {
      double res = 0.0;
      for (unsigned int i = 0; i < dim; ++i) {
        res += x[i] * y[i];
      }
      return res;
    };
    // summation order differs between implementations
    auto near = [](double expected) {
      return Catch::Matchers::WithinAbs(
          expected, 1e-12 * std::max(1.0, std::fabs(expected)));
    };

    std::vector<double> xi(dim);
    double sum =
        BFGSOpt::detail::initDirection(dim, a.data(), xi.data(), b.data());
    CHECK_THAT(sum, near(dot(b.data(), b.data())));
    for (unsigned int i = 0; i < dim; ++i) {
      CHECK(xi[i] == -a[i]);
    }

    std::vector<double> hd(dim);
    auto dots = BFGSOpt::detail::hessianTimesDGrad(dim, sym.data(), a.data(),
                                                   b.data(), hd.data());
    for (unsigned int i = 0; i < dim; ++i) {
      CHECK_THAT(hd[i], near(dot(&sym[i * dim], a.data())));
    }
    CHECK_THAT(dots.fac, near(dot(a.data(), b.data())));
    CHECK_THAT(dots.fae, near(dot(a.data(), hd.data())));
    CHECK_THAT(dots.sumDGrad, near(dot(a.data(), a.data())));
    CHECK_THAT(dots.sumXi, near(dot(b.data(), b.data())));

    const double fac = 0.7, fad = 0.2, fae = 1.3;
    auto updated = sym;
    BFGSOpt::detail::updateInverseHessian(dim, updated.data(), a.data(),
                                          b.data(), c.data(), fac, fad, fae);
    for (unsigned int i = 0; i < dim; ++i) {
      for (unsigned int j = 0; j < dim; ++j) {
        double expected = sym[i * dim + j] + fac * a[i] * a[j] -
                          fad * b[i] * b[j] + fae * c[i] * c[j];
        CHECK_THAT(updated[i * dim + j], near(expected));
        CHECK(updated[i * dim + j] == updated[j * dim + i]);
      }
    }

    BFGSOpt::detail::newDirection(dim, sym.data(), c.data(), xi.data());
    for (unsigned int i = 0; i < dim; ++i) {
      CHECK_THAT(xi[i], near(-dot(&sym[i * dim], c.data())));
    }
  }
}
}  // namespace

TEST_CASE("BFGS kernels match reference implementations") {
#ifdef RDK_USE_HIGHWAY
  for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
    INFO("Highway target " << hwy::TargetName(target));
    hwy::SetSupportedTargetsForTest(target);
    checkBFGSKernels();
  }
  hwy::SetSupportedTargetsForTest(0);
#else
  checkBFGSKernels();
#endif
}
