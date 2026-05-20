//
//  Copyright (C) 2003-2026 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/tuple.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <RDBoost/import_array.h>

#include "PickerHelpers.h"

#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>
#include <SimDivPickers/DistPicker.h>
#include <SimDivPickers/MaxMinPicker.h>
#include <SimDivPickers/HierarchicalClusterPicker.h>
#include <utility>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDPickers {

// REVIEW: the poolSize can be pulled from the numeric array
RDKit::INT_VECT MaxMinPicks(MaxMinPicker *picker, nb::object distMat,
                            int poolSize, int pickSize,
                            nb::object firstPicks, int seed) {
  if (pickSize >= poolSize) {
    throw nb::value_error("pickSize must be less than poolSize");
  }

  if (!PyArray_Check(distMat.ptr())) {
    throw nb::value_error("distance mat argument must be a numpy matrix");
  }

  PyArrayObject *copy;
  copy = (PyArrayObject *)PyArray_ContiguousFromObject(distMat.ptr(),
                                                       NPY_DOUBLE, 1, 1);
  auto *dMat = (double *)PyArray_DATA(copy);

  RDKit::INT_VECT firstPickVect;
  auto len = nb::len(firstPicks);
  for (size_t i = 0; i < len; ++i) {
    firstPickVect.push_back(nb::cast<int>(firstPicks[i]));
  }
  RDKit::INT_VECT res;
  try {
    res = picker->pick(dMat, poolSize, pickSize, firstPickVect, seed);
  } catch (...) {
    Py_DECREF(copy);
    throw;
  }
  Py_DECREF(copy);
  return res;
}

namespace {
template <typename T>
void LazyMaxMinHelper(MaxMinPicker *picker, T functor, unsigned int poolSize,
                      unsigned int pickSize, nb::object firstPicks,
                      int seed, RDKit::INT_VECT &res, double &threshold) {
  RDKit::INT_VECT firstPickVect;
  auto len = nb::len(firstPicks);
  for (size_t i = 0; i < len; ++i) {
    firstPickVect.push_back(nb::cast<int>(firstPicks[i]));
  }
  res = picker->lazyPick(functor, poolSize, pickSize, firstPickVect, seed,
                         threshold);
}
}  // end of anonymous namespace

RDKit::INT_VECT LazyMaxMinPicks(MaxMinPicker *picker, nb::object distFunc,
                                int poolSize, int pickSize,
                                nb::object firstPicks, int seed,
                                nb::object useCache) {
  if (!useCache.is_none()) {
    BOOST_LOG(rdWarningLog)
        << "the useCache argument is deprecated and ignored" << std::endl;
  }
  pyobjFunctor functor(distFunc);
  RDKit::INT_VECT res;
  double threshold = -1.;
  LazyMaxMinHelper(picker, functor, poolSize, pickSize, firstPicks, seed, res,
                   threshold);
  return res;
}

std::tuple<RDKit::INT_VECT, double> LazyMaxMinPicksWithThreshold(
    MaxMinPicker *picker, nb::object distFunc, int poolSize, int pickSize,
    double threshold, nb::object firstPicks, int seed) {
  pyobjFunctor functor(distFunc);
  RDKit::INT_VECT res;
  LazyMaxMinHelper(picker, functor, poolSize, pickSize, firstPicks, seed, res,
                   threshold);
  return std::make_tuple(res, threshold);
}

RDKit::INT_VECT LazyVectorMaxMinPicks(MaxMinPicker *picker, nb::object objs,
                                      int poolSize, int pickSize,
                                      nb::object firstPicks, int seed,
                                      nb::object useCache) {
  if (!useCache.is_none()) {
    BOOST_LOG(rdWarningLog)
        << "the useCache argument is deprecated and ignored" << std::endl;
  }
  std::vector<const ExplicitBitVect *> bvs(poolSize);
  for (int i = 0; i < poolSize; ++i) {
    bvs[i] = nb::cast<const ExplicitBitVect *>(objs[i]);
  }
  pyBVFunctor<ExplicitBitVect> functor(bvs, TANIMOTO);

  RDKit::INT_VECT res;
  double threshold = -1.;
  LazyMaxMinHelper(picker, functor, poolSize, pickSize, firstPicks, seed, res,
                   threshold);
  return res;
}

std::tuple<RDKit::INT_VECT, double> LazyVectorMaxMinPicksWithThreshold(
    MaxMinPicker *picker, nb::object objs, int poolSize, int pickSize,
    double threshold, nb::object firstPicks, int seed) {
  std::vector<const ExplicitBitVect *> bvs(poolSize);
  for (int i = 0; i < poolSize; ++i) {
    bvs[i] = nb::cast<const ExplicitBitVect *>(objs[i]);
  }
  pyBVFunctor<ExplicitBitVect> functor(bvs, TANIMOTO);

  RDKit::INT_VECT res;
  LazyMaxMinHelper(picker, functor, poolSize, pickSize, firstPicks, seed, res,
                   threshold);
  return std::make_tuple(res, threshold);
}

}  // end of namespace RDPickers

void wrap_maxminpick(nb::module_ &m) {
  nb::class_<RDPickers::MaxMinPicker>(
      m, "MaxMinPicker",
      "A class for diversity picking of items using the MaxMin Algorithm\n")
      .def("Pick", RDPickers::MaxMinPicks,
           "distMat"_a, "poolSize"_a, "pickSize"_a,
           "firstPicks"_a = nb::tuple(), "seed"_a = -1,
           R"DOC(Pick a subset of items from a pool of items using the MaxMin Algorithm
Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604

ARGUMENTS:
  - distMat: 1D distance matrix (only the lower triangle elements)
  - poolSize: number of items in the pool
  - pickSize: number of items to pick from the pool
  - firstPicks: (optional) the first items to be picked (seeds the list)
  - seed: (optional) seed for the random number generator
)DOC")

      .def("LazyPick", RDPickers::LazyMaxMinPicks,
           "distFunc"_a, "poolSize"_a, "pickSize"_a,
           "firstPicks"_a = nb::tuple(), "seed"_a = -1,
           "useCache"_a = nb::none(),
           R"DOC(Pick a subset of items from a pool of items using the MaxMin Algorithm
Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
ARGUMENTS:

  - distFunc: a function that should take two indices and return the
              distance between those two points.
              NOTE: the implementation caches distance values, so the
              client code does not need to do so; indeed, it should not.
  - poolSize: number of items in the pool
  - pickSize: number of items to pick from the pool
  - firstPicks: (optional) the first items to be picked (seeds the list)
  - seed: (optional) seed for the random number generator
  - useCache: IGNORED
)DOC")
      .def("LazyBitVectorPick", RDPickers::LazyVectorMaxMinPicks,
           "objects"_a, "poolSize"_a, "pickSize"_a,
           "firstPicks"_a = nb::tuple(), "seed"_a = -1,
           "useCache"_a = nb::none(),
           R"DOC(Pick a subset of items from a pool of bit vectors using the MaxMin Algorithm
Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
ARGUMENTS:

  - vectors: a sequence of the bit vectors that should be picked from.
  - poolSize: number of items in the pool
  - pickSize: number of items to pick from the pool
  - firstPicks: (optional) the first items to be picked (seeds the list)
  - seed: (optional) seed for the random number generator
  - useCache: IGNORED.
)DOC")

      .def("LazyPickWithThreshold", RDPickers::LazyMaxMinPicksWithThreshold,
           "distFunc"_a, "poolSize"_a, "pickSize"_a,
           "threshold"_a,
           "firstPicks"_a = nb::tuple(), "seed"_a = -1,
           R"DOC(Pick a subset of items from a pool of items using the MaxMin Algorithm
Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
ARGUMENTS:

  - distFunc: a function that should take two indices and return the
              distance between those two points.
              NOTE: the implementation caches distance values, so the
              client code does not need to do so; indeed, it should not.
  - poolSize: number of items in the pool
  - pickSize: number of items to pick from the pool
  - threshold: stop picking when the distance goes below this value
  - firstPicks: (optional) the first items to be picked (seeds the list)
  - seed: (optional) seed for the random number generator
)DOC")
      .def("LazyBitVectorPickWithThreshold",
           RDPickers::LazyVectorMaxMinPicksWithThreshold,
           "objects"_a, "poolSize"_a, "pickSize"_a,
           "threshold"_a,
           "firstPicks"_a = nb::tuple(), "seed"_a = -1,
           R"DOC(Pick a subset of items from a pool of bit vectors using the MaxMin Algorithm
Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
ARGUMENTS:

  - vectors: a sequence of the bit vectors that should be picked from.
  - poolSize: number of items in the pool
  - pickSize: number of items to pick from the pool
  - threshold: stop picking when the distance goes below this value
  - firstPicks: (optional) the first items to be picked (seeds the list)
  - seed: (optional) seed for the random number generator
)DOC");
}
