//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <cstdlib>

#include <RDGeneral/Invariant.h>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

typedef double real;

namespace nb = nanobind;
using namespace nb::literals;

extern "C" void distdriver_(long *n, long *len, real *dists, long *toggle,
                            long *ia, long *ib, real *crit);

//
// Rather than deal with any nonsense like trying to get
// the distance matrix built properly on the f2c side of things
// (thus drowning in the waves of f2c hate), we'll generate
// the distance matrix on our own here and then call distdriver_
//
static void clusterit(const real *dataP, long n, long m, long iopt, long *ia,
                      long *ib, real *crit) {
  real *dists;
  long len;
  long pos = 0;
  long i, j, k, iTab, jTab;
  double tmp;
  len = (n * (n - 1)) / 2;
  dists = static_cast<real *>(calloc(static_cast<size_t>(len), sizeof(real)));
  CHECK_INVARIANT(dists, "failed to allocate memory");
  for (i = 1; i < n; i++) {
    iTab = i * m;
    for (j = 0; j < i; j++) {
      jTab = j * m;
      for (k = 0; k < m; k++) {
        tmp = dataP[iTab + k] - dataP[jTab + k];
        dists[pos] += tmp * tmp;
      }
      pos++;
    }
  }
  distdriver_(&n, &len, dists, &iopt, ia, ib, crit);
  free(dists);
};

static nb::tuple Clustering_MurtaghCluster(
    nb::ndarray<nb::numpy, const double, nb::ndim<2>, nb::c_contig> data,
    int nPts, int sz, int option) {
  auto dataRows = static_cast<long>(data.shape(0));
  auto dataCols = static_cast<long>(data.shape(1));
  CHECK_INVARIANT(dataRows == nPts, "data row count does not match nPts");
  CHECK_INVARIANT(dataCols == sz, "data column count does not match sz");

  auto *ia =
      static_cast<long *>(calloc(static_cast<size_t>(nPts), sizeof(long)));
  CHECK_INVARIANT(ia, "failed to allocate memory");
  auto *ib =
      static_cast<long *>(calloc(static_cast<size_t>(nPts), sizeof(long)));
  CHECK_INVARIANT(ib, "failed to allocate memory");
  auto *crit =
      static_cast<real *>(calloc(static_cast<size_t>(nPts), sizeof(real)));
  CHECK_INVARIANT(crit, "failed to allocate memory");

  clusterit(data.data(), static_cast<long>(nPts), static_cast<long>(sz),
            static_cast<long>(option), ia, ib, crit);

  nb::list res;

  nb::capsule iaOwner(ia, [](void *p) noexcept { free(p); });
  res.append(nb::ndarray<nb::numpy, long, nb::ndim<1>>(
      ia, {static_cast<size_t>(nPts)}, iaOwner));

  nb::capsule ibOwner(ib, [](void *p) noexcept { free(p); });
  res.append(nb::ndarray<nb::numpy, long, nb::ndim<1>>(
      ib, {static_cast<size_t>(nPts)}, ibOwner));

  nb::capsule critOwner(crit, [](void *p) noexcept { free(p); });
  res.append(nb::ndarray<nb::numpy, double, nb::ndim<1>>(
      crit, {static_cast<size_t>(nPts)}, critOwner));

  return nb::tuple(res);
}

void distclusterit(real *dists, long n, long iopt, long *ia, long *ib,
                   real *crit) {
  long len;

  len = (n * (n - 1)) / 2;
  distdriver_(&n, &len, dists, &iopt, ia, ib, crit);
};

static nb::tuple Clustering_MurtaghDistCluster(
    nb::ndarray<nb::numpy, const double, nb::ndim<1>, nb::c_contig> data,
    int nPts, int option) {
  auto dataLen = static_cast<long>(data.shape(0));
  CHECK_INVARIANT(dataLen > 0, "distance matrix must not be empty");

  auto *ia =
      static_cast<long *>(calloc(static_cast<size_t>(nPts), sizeof(long)));
  CHECK_INVARIANT(ia, "failed to allocate memory");
  auto *ib =
      static_cast<long *>(calloc(static_cast<size_t>(nPts), sizeof(long)));
  CHECK_INVARIANT(ib, "failed to allocate memory");
  auto *crit =
      static_cast<real *>(calloc(static_cast<size_t>(nPts), sizeof(real)));
  CHECK_INVARIANT(crit, "failed to allocate memory");

  distclusterit(const_cast<real *>(data.data()), static_cast<long>(nPts),
                static_cast<long>(option), ia, ib, crit);

  nb::list res;

  nb::capsule iaOwner(ia, [](void *p) noexcept { free(p); });
  res.append(nb::ndarray<nb::numpy, long, nb::ndim<1>>(
      ia, {static_cast<size_t>(nPts)}, iaOwner));

  nb::capsule ibOwner(ib, [](void *p) noexcept { free(p); });
  res.append(nb::ndarray<nb::numpy, long, nb::ndim<1>>(
      ib, {static_cast<size_t>(nPts)}, ibOwner));

  nb::capsule critOwner(crit, [](void *p) noexcept { free(p); });
  res.append(nb::ndarray<nb::numpy, double, nb::ndim<1>>(
      crit, {static_cast<size_t>(nPts)}, critOwner));

  return nb::tuple(res);
}

NB_MODULE(Clustering, m) {
  m.def("MurtaghCluster", Clustering_MurtaghCluster, "data"_a, "nPts"_a, "sz"_a,
        "option"_a,
        R"DOC(Cluster points using Murtagh's hierarchical clustering algorithm.

ARGUMENTS:
  - data: 2D NumPy array of point coordinates
  - nPts: number of points in the array
  - sz: number of coordinate values per point
  - option: clustering option passed to the underlying driver

RETURNS:
  - tuple of three 1D NumPy arrays containing the cluster results
)DOC");
  m.def("MurtaghDistCluster", Clustering_MurtaghDistCluster, "data"_a, "nPts"_a,
        "option"_a,
        R"DOC(Cluster using a precomputed condensed distance matrix.

ARGUMENTS:
  - data: 1D NumPy array containing the lower triangle of the distance matrix
  - nPts: number of points in the array
  - option: clustering option passed to the underlying driver

RETURNS:
  - tuple of three 1D NumPy arrays containing the cluster results
)DOC");
}
