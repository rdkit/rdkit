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
#include <cmath>
#include <RDGeneral/Invariant.h>
#include <catch2/catch_all.hpp>

#include "BFGSOpt.h"

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
