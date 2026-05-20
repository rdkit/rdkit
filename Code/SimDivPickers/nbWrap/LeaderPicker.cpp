//
//  Copyright (C) 2019-2026 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

#include "PickerHelpers.h"

#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>
#include <SimDivPickers/DistPicker.h>
#include <SimDivPickers/LeaderPicker.h>
#include <utility>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDPickers {
namespace {
template <typename T>
void LazyLeaderHelper(LeaderPicker *picker, T functor, unsigned int poolSize,
                      double &threshold, unsigned int pickSize,
                      nb::object firstPicks, RDKit::INT_VECT &res,
                      int nThreads) {
  RDKit::INT_VECT firstPickVect;
  auto len = nb::len(firstPicks);
  for (size_t i = 0; i < len; ++i) {
    firstPickVect.push_back(nb::cast<int>(firstPicks[i]));
  }
  res = picker->lazyPick(functor, poolSize, pickSize, firstPickVect, threshold,
                         nThreads);
}
}  // end of anonymous namespace

RDKit::INT_VECT LazyVectorLeaderPicks(LeaderPicker *picker, nb::object objs,
                                      int poolSize, double threshold,
                                      int pickSize, nb::object firstPicks,
                                      int numThreads) {
  std::vector<const ExplicitBitVect *> bvs(poolSize);
  for (int i = 0; i < poolSize; ++i) {
    bvs[i] = nb::cast<const ExplicitBitVect *>(objs[i]);
  }
  pyBVFunctor<ExplicitBitVect> functor(bvs, TANIMOTO);
  RDKit::INT_VECT res;
  LazyLeaderHelper(picker, functor, poolSize, threshold, pickSize, firstPicks,
                   res, numThreads);
  return res;
}

RDKit::INT_VECT LazyLeaderPicks(LeaderPicker *picker, nb::object distFunc,
                                int poolSize, double threshold, int pickSize,
                                nb::object firstPicks, int numThreads) {
  pyobjFunctor functor(distFunc);
  RDKit::INT_VECT res;
  LazyLeaderHelper(picker, functor, poolSize, threshold, pickSize, firstPicks,
                   res, numThreads);
  return res;
}

}  // end of namespace RDPickers

void wrap_leaderpick(nb::module_ &m) {
  nb::class_<RDPickers::LeaderPicker>(
      m, "LeaderPicker",
      R"DOC(A class for diversity picking of items using Roger Sayle's Leader
algorithm (analogous to sphere exclusion). The algorithm is
currently unpublished, but a description is available in this
presentation from the 2019 RDKit UGM:
https://github.com/rdkit/UGM_2019/raw/master/Presentations/Sayle_Clustering.pdf
)DOC")
      .def("LazyBitVectorPick", RDPickers::LazyVectorLeaderPicks,
           "objects"_a, "poolSize"_a, "threshold"_a,
           "pickSize"_a = 0, "firstPicks"_a = nb::tuple(),
           "numThreads"_a = 1,
           "Pick a subset of items from a collection of bit vectors using "
           "Tanimoto distance. The threshold value is a "
           "*distance* (i.e. 1-similarity). Note that the numThreads "
           "argument is currently ignored.")
      .def("LazyPick", RDPickers::LazyLeaderPicks,
           "distFunc"_a, "poolSize"_a, "threshold"_a,
           "pickSize"_a = 0, "firstPicks"_a = nb::tuple(),
           "numThreads"_a = 1,
           "Pick a subset of items from a pool of items using the "
           "user-provided function to determine distances. Note that the "
           "numThreads argument is currently ignored.");
}
