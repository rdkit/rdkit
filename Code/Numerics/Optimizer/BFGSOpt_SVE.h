#ifndef RDKIT_NUMERICS_OPTIMIZER_BFGSOPT_SVE_H
#define RDKIT_NUMERICS_OPTIMIZER_BFGSOPT_SVE_H

#if defined(__linux__) && defined(__aarch64__)
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

}  // namespace BFGSOpt

#endif  // RDKIT_NUMERICS_OPTIMIZER_BFGSOPT_SVE_H
