//
//  Copyright (C) 2019  Greg Landrum
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
#include <SimDivPickers/LeaderPicker.h>
#include <iostream>
#include <utility>

namespace python = boost::python;
namespace RDPickers {
namespace {
template <typename T>
void LazyLeaderHelper(LeaderPicker *picker, T functor, unsigned int poolSize,
                      double &threshold, unsigned int pickSize,
                      python::object firstPicks, RDKit::INT_VECT &res,
                      int nThreads) {
  RDKit::INT_VECT firstPickVect;
  for (unsigned int i = 0;
       i < python::extract<unsigned int>(firstPicks.attr("__len__")()); ++i) {
    firstPickVect.push_back(python::extract<int>(firstPicks[i]));
  }
  res = picker->lazyPick(functor, poolSize, pickSize, firstPickVect, threshold,
                         nThreads);
}
}  // end of anonymous namespace

RDKit::INT_VECT LazyVectorLeaderPicks(LeaderPicker *picker, python::object objs,
                                      int poolSize, double threshold,
                                      int pickSize, python::object firstPicks,
                                      int numThreads) {
  std::vector<const ExplicitBitVect *> bvs(poolSize);
  for (int i = 0; i < poolSize; ++i) {
    bvs[i] = python::extract<const ExplicitBitVect *>(objs[i]);
  }
  pyBVFunctor<ExplicitBitVect> functor(bvs, TANIMOTO);
  RDKit::INT_VECT res;
  LazyLeaderHelper(picker, functor, poolSize, threshold, pickSize, firstPicks,
                   res, numThreads);
  return res;
}

RDKit::INT_VECT LazyLeaderPicks(LeaderPicker *picker, python::object distFunc,
                                int poolSize, double threshold, int pickSize,
                                python::object firstPicks, int numThreads) {
  pyobjFunctor functor(distFunc);
  RDKit::INT_VECT res;
  LazyLeaderHelper(picker, functor, poolSize, threshold, pickSize, firstPicks,
                   res, numThreads);
  return res;
}

}  // end of namespace RDPickers

struct LeaderPicker_wrap {
  static void wrap() {
    python::class_<RDPickers::LeaderPicker>(
        "LeaderPicker",
        "A class for diversity picking of items using Roger Sayle's Leader "
        "algorithm (analogous to sphere exclusion). The algorithm is "
        "currently unpublished, but a description is available in this "
        "presentation from the 2019 RDKit UGM: "
        "https://github.com/rdkit/UGM_2019/raw/master/Presentations/"
        "Sayle_Clustering.pdf\n")
        .def("LazyBitVectorPick", RDPickers::LazyVectorLeaderPicks,
             (python::arg("self"), python::arg("objects"),
              python::arg("poolSize"), python::arg("threshold"),
              python::arg("pickSize") = 0,
              python::arg("firstPicks") = python::tuple(),
              python::arg("numThreads") = 1),
             "Pick a subset of items from a collection of bit vectors using "
             "Tanimoto distance. The threshold value is a "
             "*distance* (i.e. 1-similarity). Note that the numThreads "
             "argument is currently ignored.")
        .def("LazyPick", RDPickers::LazyLeaderPicks,
             (python::arg("self"), python::arg("distFunc"),
              python::arg("poolSize"), python::arg("threshold"),
              python::arg("pickSize") = 0,
              python::arg("firstPicks") = python::tuple(),
              python::arg("numThreads") = 1),
             "Pick a subset of items from a pool of items using the "
             "user-provided function to determine distances. Note that the "
             "numThreads argument is currently ignored.");
  };
};

void wrap_leaderpick() { LeaderPicker_wrap::wrap(); }
