//
//  Copyright (C) 2003-2017  Greg Landrum and Rational Discovery LLC
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define NO_IMPORT_ARRAY

#define PY_ARRAY_UNIQUE_SYMBOL rdpicker_array_API
#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>
#include <RDBoost/boost_numpy.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <map>

#include "PickerHelpers.h"

#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>
#include <SimDivPickers/DistPicker.h>
#include <SimDivPickers/MaxMinPicker.h>
#include <SimDivPickers/HierarchicalClusterPicker.h>
#include <iostream>
#include <utility>

namespace python = boost::python;
namespace RDPickers {

// REVIEW: the poolSize can be pulled from the numeric array
RDKit::INT_VECT MaxMinPicks(MaxMinPicker *picker, python::object distMat,
                            int poolSize, int pickSize,
                            python::object firstPicks, int seed) {
  if (pickSize >= poolSize) {
    throw ValueErrorException("pickSize must be less than poolSize");
  }

  if (!PyArray_Check(distMat.ptr())) {
    throw ValueErrorException("distance mat argument must be a numpy matrix");
  }

  PyArrayObject *copy;
  copy = (PyArrayObject *)PyArray_ContiguousFromObject(distMat.ptr(),
                                                       NPY_DOUBLE, 1, 1);
  auto *dMat = (double *)PyArray_DATA(copy);

  RDKit::INT_VECT firstPickVect;
  for (unsigned int i = 0;
       i < python::extract<unsigned int>(firstPicks.attr("__len__")()); ++i) {
    firstPickVect.push_back(python::extract<int>(firstPicks[i]));
  }
  RDKit::INT_VECT res =
      picker->pick(dMat, poolSize, pickSize, firstPickVect, seed);
  Py_DECREF(copy);
  return res;
}

namespace {
template <typename T>
void LazyMaxMinHelper(MaxMinPicker *picker, T functor, unsigned int poolSize,
                      unsigned int pickSize, python::object firstPicks,
                      int seed, RDKit::INT_VECT &res, double &threshold) {
  RDKit::INT_VECT firstPickVect;
  for (unsigned int i = 0;
       i < python::extract<unsigned int>(firstPicks.attr("__len__")()); ++i) {
    firstPickVect.push_back(python::extract<int>(firstPicks[i]));
  }
  res = picker->lazyPick(functor, poolSize, pickSize, firstPickVect, seed,
                         threshold);
}
}  // end of anonymous namespace

RDKit::INT_VECT LazyMaxMinPicks(MaxMinPicker *picker, python::object distFunc,
                                int poolSize, int pickSize,
                                python::object firstPicks, int seed,
                                python::object useCache) {
  if (useCache != python::object()) {
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
python::tuple LazyMaxMinPicksWithThreshold(
    MaxMinPicker *picker, python::object distFunc, int poolSize, int pickSize,
    double threshold, python::object firstPicks, int seed) {
  pyobjFunctor functor(distFunc);
  RDKit::INT_VECT res;
  LazyMaxMinHelper(picker, functor, poolSize, pickSize, firstPicks, seed, res,
                   threshold);
  return python::make_tuple(res, threshold);
}

RDKit::INT_VECT LazyVectorMaxMinPicks(MaxMinPicker *picker, python::object objs,
                                      int poolSize, int pickSize,
                                      python::object firstPicks, int seed,
                                      python::object useCache) {
  if (useCache != python::object()) {
    BOOST_LOG(rdWarningLog)
        << "the useCache argument is deprecated and ignored" << std::endl;
  }
  std::vector<const ExplicitBitVect *> bvs(poolSize);
  for (int i = 0; i < poolSize; ++i) {
    bvs[i] = python::extract<const ExplicitBitVect *>(objs[i]);
  }
  pyBVFunctor<ExplicitBitVect> functor(bvs, TANIMOTO);

  RDKit::INT_VECT res;
  double threshold = -1.;
  LazyMaxMinHelper(picker, functor, poolSize, pickSize, firstPicks, seed, res,
                   threshold);
  return res;
}

python::tuple LazyVectorMaxMinPicksWithThreshold(
    MaxMinPicker *picker, python::object objs, int poolSize, int pickSize,
    double threshold, python::object firstPicks, int seed) {
  std::vector<const ExplicitBitVect *> bvs(poolSize);
  for (int i = 0; i < poolSize; ++i) {
    bvs[i] = python::extract<const ExplicitBitVect *>(objs[i]);
  }
  pyBVFunctor<ExplicitBitVect> functor(bvs, TANIMOTO);

  RDKit::INT_VECT res;
  LazyMaxMinHelper(picker, functor, poolSize, pickSize, firstPicks, seed, res,
                   threshold);
  return python::make_tuple(res, threshold);
}

}  // end of namespace RDPickers

struct MaxMin_wrap {
  static void wrap() {
    python::class_<RDPickers::MaxMinPicker>(
        "MaxMinPicker",
        "A class for diversity picking of items using the MaxMin Algorithm\n")
        .def("Pick", RDPickers::MaxMinPicks,
             (python::arg("self"), python::arg("distMat"),
              python::arg("poolSize"), python::arg("pickSize"),
              python::arg("firstPicks") = python::tuple(),
              python::arg("seed") = -1),
             "Pick a subset of items from a pool of items using the MaxMin "
             "Algorithm\n"
             "Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), "
             "598-604 \n\n"
             "ARGUMENTS:\n"
             "  - distMat: 1D distance matrix (only the lower triangle "
             "elements)\n"
             "  - poolSize: number of items in the pool\n"
             "  - pickSize: number of items to pick from the pool\n"
             "  - firstPicks: (optional) the first items to be picked (seeds "
             "the list)\n"
             "  - seed: (optional) seed for the random number generator\n")

        .def("LazyPick", RDPickers::LazyMaxMinPicks,
             (python::arg("self"), python::arg("distFunc"),
              python::arg("poolSize"), python::arg("pickSize"),
              python::arg("firstPicks") = python::tuple(),
              python::arg("seed") = -1,
              python::arg("useCache") = python::object()),
             "Pick a subset of items from a pool of items using the MaxMin "
             "Algorithm\n"
             "Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), "
             "598-604 \n"
             "ARGUMENTS:\n\n"
             "  - distFunc: a function that should take two indices and return "
             "the\n"
             "              distance between those two points.\n"
             "              NOTE: the implementation caches distance values, "
             "so the\n"
             "              client code does not need to do so; indeed, it "
             "should not.\n"
             "  - poolSize: number of items in the pool\n"
             "  - pickSize: number of items to pick from the pool\n"
             "  - firstPicks: (optional) the first items to be picked (seeds "
             "the list)\n"
             "  - seed: (optional) seed for the random number generator\n"
             "  - useCache: IGNORED\n")
        .def("LazyBitVectorPick", RDPickers::LazyVectorMaxMinPicks,
             (python::arg("self"), python::arg("objects"),
              python::arg("poolSize"), python::arg("pickSize"),
              python::arg("firstPicks") = python::tuple(),
              python::arg("seed") = -1,
              python::arg("useCache") = python::object()),
             "Pick a subset of items from a pool of bit vectors using the "
             "MaxMin Algorithm\n"
             "Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), "
             "598-604 \n"
             "ARGUMENTS:\n\n"
             "  - vectors: a sequence of the bit vectors that should be picked "
             "from.\n"
             "  - poolSize: number of items in the pool\n"
             "  - pickSize: number of items to pick from the pool\n"
             "  - firstPicks: (optional) the first items to be picked (seeds "
             "the list)\n"
             "  - seed: (optional) seed for the random number generator\n"
             "  - useCache: IGNORED.\n")

        .def("LazyPickWithThreshold", RDPickers::LazyMaxMinPicksWithThreshold,
             (python::arg("self"), python::arg("distFunc"),
              python::arg("poolSize"), python::arg("pickSize"),
              python::arg("threshold"),
              python::arg("firstPicks") = python::tuple(),
              python::arg("seed") = -1),
             "Pick a subset of items from a pool of items using the MaxMin "
             "Algorithm\n"
             "Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), "
             "598-604 \n"
             "ARGUMENTS:\n\n"
             "  - distFunc: a function that should take two indices and return "
             "the\n"
             "              distance between those two points.\n"
             "              NOTE: the implementation caches distance values, "
             "so the\n"
             "              client code does not need to do so; indeed, it "
             "should not.\n"
             "  - poolSize: number of items in the pool\n"
             "  - pickSize: number of items to pick from the pool\n"
             "  - threshold: stop picking when the distance goes below this "
             "value\n"
             "  - firstPicks: (optional) the first items to be picked (seeds "
             "the list)\n"
             "  - seed: (optional) seed for the random number generator\n")
        .def("LazyBitVectorPickWithThreshold",
             RDPickers::LazyVectorMaxMinPicksWithThreshold,
             (python::arg("self"), python::arg("objects"),
              python::arg("poolSize"), python::arg("pickSize"),
              python::arg("threshold"),
              python::arg("firstPicks") = python::tuple(),
              python::arg("seed") = -1),
             "Pick a subset of items from a pool of bit vectors using the "
             "MaxMin Algorithm\n"
             "Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), "
             "598-604 \n"
             "ARGUMENTS:\n\n"
             "  - vectors: a sequence of the bit vectors that should be picked "
             "from.\n"
             "  - poolSize: number of items in the pool\n"
             "  - pickSize: number of items to pick from the pool\n"
             "  - threshold: stop picking when the distance goes below this "
             "value\n"
             "  - firstPicks: (optional) the first items to be picked (seeds "
             "the list)\n"
             "  - seed: (optional) seed for the random number generator\n");
  };
};

void wrap_maxminpick() { MaxMin_wrap::wrap(); }
