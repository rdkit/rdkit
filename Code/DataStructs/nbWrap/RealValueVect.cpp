//
//  Copyright (C) 2026 greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <DataStructs/RealValueVect.h>
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <RDBoost/Wrap_nb.h>

using namespace RDKit;
namespace nb = nanobind;
using namespace nb::literals;

std::string realValVectDoc =
    "A container class for storing real\n\
values.\n\
\n\
The length of the vector is set at construction time.\n\
\n\
As you would expect, _RealValueVects_ support a set of binary operations\n\
so you can do things like:\n\
  rvv3 = rvv1 & rvv2  the result contains the smallest value in each entry\n\
  rvv3 = rvv1 | rvv2  the result contains the largest value in each entry\n\
  rvv1 += rvv2     \n\
  rvv3 = rvv1 + rvv2    \n\
  rvv1 -= rvv3    \n\
  rvv3 = rvv1 - rvv2    \n\
\n\
Elements can be set and read using indexing (i.e. bv[i] = 4 or val=bv[i])\n\
\n";

namespace {
std::string rvvToString(const RealValueVect &rvv) { return rvv.toString(); }
}  // namespace

struct realValVec_wrapper {
  static void wrap(nb::module_ &m) {
    nb::class_<RealValueVect>(m, "RealValueVect", realValVectDoc.c_str())
        .def(nb::new_([]() { return new RealValueVect(); }))
        .def(nb::new_(
                 [](unsigned int length) { return new RealValueVect(length); }),
             "length"_a)
        .def(nb::new_([](nb::bytes b) {
               return new RealValueVect(
                   std::string(static_cast<const char *>(b.data()),
                               static_cast<size_t>(b.size())));
             }),
             "pkl"_a)
        .def("__len__", &RealValueVect::getLength,
             "Get the number of entries in the vector")
        .def("__setitem__", &RealValueVect::setVal,
             "Set the value at a specified location")
        .def("__getitem__", &RealValueVect::getVal,
             "Get the value at a specified location")
        .def(nb::self & nb::self)
        .def(nb::self | nb::self)
        .def(nb::self - nb::self)
        .def(nb::self -= nb::self)
        .def(nb::self + nb::self)
        .def(nb::self += nb::self)
        .def("GetTotalVal", &RealValueVect::getTotalVal,
             "Get the sum of the values in the vector, basically L1 norm")
        .def("__getstate__", getObjectState<RealValueVect, rvvToString>)
        .def("__setstate__", setObjectState<RealValueVect>)
        .doc() = realValVectDoc.c_str();
    m.def("ComputeL1Norm", computeL1Norm,
          "Compute the distance between two real vector values\n");
  }
};

void wrap_realValVect(nb::module_ &m) { realValVec_wrapper::wrap(m); }
