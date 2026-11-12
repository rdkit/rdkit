//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Highway implementation of the BFGS kernels, used in place of
// BFGSKernels.cpp when RDK_USE_HIGHWAY is enabled. foreach_target.h
// re-includes this file once per SIMD target and the best one supported by
// the CPU is chosen at runtime.
//
// Note that results are not bitwise identical across targets: the vector
// width changes the summation order of the dot products, and FMA is only
// fused on targets that have it.
//
#include "BFGSKernels.h"

#include <algorithm>
#include <cstddef>

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "Numerics/Optimizer/BFGSKernelsHwy.cpp"
#include <hwy/foreach_target.h>  // must come before highway.h
#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();
namespace BFGSOpt::detail::HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

using D = hn::ScalableTag<double>;
using V = hn::Vec<D>;

// Calls body(j, count) for each vector-sized chunk of [begin, end); count is
// less than the number of lanes only for the final chunk, which should be
// loaded and stored with LoadN/StoreN.
template <typename Body>
HWY_INLINE void forEachChunk(D d, size_t begin, size_t end, Body &&body) {
  const size_t N = hn::Lanes(d);
  size_t j = begin;
  for (; j + N <= end; j += N) {
    body(j, N);
  }
  if (j < end) {
    body(j, end - j);
  }
}

// LoadN zero-fills the lanes past count, so partial chunks contribute
// nothing to sums.
HWY_INLINE V load(D d, const double *HWY_RESTRICT p, size_t count) {
  return count == hn::Lanes(d) ? hn::LoadU(d, p) : hn::LoadN(d, p, count);
}

HWY_INLINE void store(D d, V v, double *HWY_RESTRICT p, size_t count) {
  if (count == hn::Lanes(d)) {
    hn::StoreU(v, d, p);
  } else {
    hn::StoreN(v, d, p, count);
  }
}

HWY_INLINE double dot(D d, const double *HWY_RESTRICT a,
                      const double *HWY_RESTRICT b, size_t n) {
  V acc = hn::Zero(d);
  forEachChunk(d, 0, n, [&](size_t j, size_t count) HWY_ATTR {
    acc = hn::MulAdd(load(d, a + j, count), load(d, b + j, count), acc);
  });
  return hn::ReduceSum(d, acc);
}

double InitDirection(unsigned int dim, const double *HWY_RESTRICT grad,
                     double *HWY_RESTRICT xi, const double *HWY_RESTRICT pos) {
  const D d;
  V acc = hn::Zero(d);
  forEachChunk(d, 0, dim, [&](size_t i, size_t count) HWY_ATTR {
    store(d, hn::Neg(load(d, grad + i, count)), xi + i, count);
    const V p = load(d, pos + i, count);
    acc = hn::MulAdd(p, p, acc);
  });
  return hn::ReduceSum(d, acc);
}

UpdateDots HessianTimesDGrad(unsigned int dim,
                             const double *HWY_RESTRICT invHessian,
                             const double *HWY_RESTRICT dGrad,
                             const double *HWY_RESTRICT xi,
                             double *HWY_RESTRICT hessDGrad) {
  const D d;
  for (size_t i = 0; i < dim; ++i) {
    hessDGrad[i] = dot(d, invHessian + i * dim, dGrad, dim);
  }

  // all four dot products in a single pass
  V fac = hn::Zero(d), fae = hn::Zero(d);
  V sumDGrad = hn::Zero(d), sumXi = hn::Zero(d);
  forEachChunk(d, 0, dim, [&](size_t i, size_t count) HWY_ATTR {
    const V dg = load(d, dGrad + i, count);
    const V x = load(d, xi + i, count);
    fac = hn::MulAdd(dg, x, fac);
    fae = hn::MulAdd(dg, load(d, hessDGrad + i, count), fae);
    sumDGrad = hn::MulAdd(dg, dg, sumDGrad);
    sumXi = hn::MulAdd(x, x, sumXi);
  });
  return {hn::ReduceSum(d, fac), hn::ReduceSum(d, fae),
          hn::ReduceSum(d, sumDGrad), hn::ReduceSum(d, sumXi)};
}

void UpdateInverseHessian(unsigned int dim, double *HWY_RESTRICT invHessian,
                          const double *HWY_RESTRICT xi,
                          const double *HWY_RESTRICT hessDGrad,
                          const double *HWY_RESTRICT dGrad, double fac,
                          double fad, double fae) {
  const D d;
  // update the upper triangle (j >= i) and mirror it, which keeps the matrix
  // exactly symmetric
  for (size_t i = 0; i < dim; ++i) {
    const V pxi = hn::Set(d, fac * xi[i]);
    const V hdgi = hn::Set(d, fad * hessDGrad[i]);
    const V dgi = hn::Set(d, fae * dGrad[i]);
    double *row = invHessian + i * dim;
    forEachChunk(d, i, dim, [&](size_t j, size_t count) HWY_ATTR {
      V h = load(d, row + j, count);
      h = hn::MulAdd(pxi, load(d, xi + j, count), h);
      h = hn::NegMulAdd(hdgi, load(d, hessDGrad + j, count), h);
      h = hn::MulAdd(dgi, load(d, dGrad + j, count), h);
      store(d, h, row + j, count);
    });
    for (size_t j = i + 1; j < dim; ++j) {
      invHessian[j * dim + i] = row[j];
    }
  }
}

void NewDirection(unsigned int dim, const double *HWY_RESTRICT invHessian,
                  const double *HWY_RESTRICT grad, double *HWY_RESTRICT xi) {
  const D d;
  for (size_t i = 0; i < dim; ++i) {
    xi[i] = -dot(d, invHessian + i * dim, grad, dim);
  }
}

}  // namespace BFGSOpt::detail::HWY_NAMESPACE
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace BFGSOpt::detail {
HWY_EXPORT(InitDirection);
HWY_EXPORT(HessianTimesDGrad);
HWY_EXPORT(UpdateInverseHessian);
HWY_EXPORT(NewDirection);

double initDirection(unsigned int dim, const double *grad, double *xi,
                     const double *pos) {
  return HWY_DYNAMIC_DISPATCH(InitDirection)(dim, grad, xi, pos);
}

UpdateDots hessianTimesDGrad(unsigned int dim, const double *invHessian,
                             const double *dGrad, const double *xi,
                             double *hessDGrad) {
  return HWY_DYNAMIC_DISPATCH(HessianTimesDGrad)(dim, invHessian, dGrad, xi,
                                                 hessDGrad);
}

void updateInverseHessian(unsigned int dim, double *invHessian,
                          const double *xi, const double *hessDGrad,
                          const double *dGrad, double fac, double fad,
                          double fae) {
  HWY_DYNAMIC_DISPATCH(UpdateInverseHessian)
  (dim, invHessian, xi, hessDGrad, dGrad, fac, fad, fae);
}

void newDirection(unsigned int dim, const double *invHessian,
                  const double *grad, double *xi) {
  HWY_DYNAMIC_DISPATCH(NewDirection)(dim, invHessian, grad, xi);
}
}  // namespace BFGSOpt::detail
#endif  // HWY_ONCE
