//
//  Copyright (C) 2026 Greg Landrum
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <DataStructs/DiscreteValueVect.h>
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>

using namespace RDKit;
namespace nb = nanobind;
using namespace nb::literals;

std::string disValVectDoc =
    "A container class for storing unsigned integer\n\
values within a particular range.\n\
\n\
The length of the vector and type of its elements (determines the maximum value\n\
that can be stored) are both set at construction time.\n\
\n\
As you would expect, _DiscreteValueVects_ support a set of binary operations\n\
so you can do things like:\n\
  dvv3 = dvv1 & dvv2  the result contains the smallest value in each entry\n\
  dvv3 = dvv1 | dvv2  the result contains the largest value in each entry\n\
  dvv1 += dvv2     values are truncated when necessary\n\
  dvv3 = dvv1 + dvv2    values are truncated when necessary\n\
  dvv1 -= dvv3    would-be negative values are set to zero\n\
  dvv3 = dvv1 - dvv2    would-be negative values are set to zero\n\
\n\
Elements can be set and read using indexing (i.e. bv[i] = 4 or val=bv[i])\n\
\n";

struct discreteValVec_wrapper {
  static void wrap(nb::module_ &m) {
    nb::enum_<DiscreteValueVect::DiscreteValueType>(m, "DiscreteValueType")
        .value("ONEBITVALUE", DiscreteValueVect::ONEBITVALUE)
        .value("TWOBITVALUE", DiscreteValueVect::TWOBITVALUE)
        .value("FOURBITVALUE", DiscreteValueVect::FOURBITVALUE)
        .value("EIGHTBITVALUE", DiscreteValueVect::EIGHTBITVALUE)
        .value("SIXTEENBITVALUE", DiscreteValueVect::SIXTEENBITVALUE)
        .export_values();

    nb::class_<DiscreteValueVect>(m, "DiscreteValueVect")
        .def(nb::init<DiscreteValueVect::DiscreteValueType, unsigned int>(),
             "valType"_a, "length"_a, "Constructor")
        .def("__init__",
             [](DiscreteValueVect *t, nb::bytes b) {
               new (t) DiscreteValueVect(
                   std::string(static_cast<const char *>(b.data()),
                               static_cast<size_t>(b.size())));
             })
        .def("__len__", &DiscreteValueVect::getLength,
             "Get the number of entries in the vector")
        .def("__setitem__", &DiscreteValueVect::setVal, "i"_a, "val"_a,
             "Set the value at a specified location")
        .def("__getitem__", &DiscreteValueVect::getVal, "i"_a,
             "Get the value at a specified location")
        .def(nb::self & nb::self)
        .def(nb::self | nb::self)
        .def(nb::self - nb::self)
#ifdef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wself-assign-overloaded"
#endif
        .def(nb::self -= nb::self)  // clang warns incorrectly on this construct
#ifdef __clang__
#pragma GCC diagnostic pop
#endif
        .def(nb::self + nb::self)
        .def(nb::self += nb::self)

        .def("GetValueType", &DiscreteValueVect::getValueType,
             "Get the type of value stored in the vector")
        .def("GetTotalVal", &DiscreteValueVect::getTotalVal,
             "Get the sum of the values in the vector, basically L1 norm")

        // FIX: probably want to include helper functionality for
        // working with ctors that expect binary strings. nanobind 
        // works with bytes.
        .def("__getstate__",
             [](const DiscreteValueVect &dvv) {
               const auto pkl = dvv.toString();
               return std::make_tuple(nb::bytes(pkl.c_str(), pkl.size()));
             })
        .def("__setstate__",
             [](DiscreteValueVect &dvv, const std::tuple<nb::bytes> &state) {
               std::string pkl = std::string(
                   static_cast<const char *>(std::get<0>(state).data()),
                   static_cast<size_t>(std::get<0>(state).size()));
               new (&dvv) DiscreteValueVect(pkl);
             })
        .doc() = disValVectDoc.c_str()
        //
        ;

    m.def("ComputeL1Norm", computeL1Norm, "v1"_a, "v2"_a,
          "Compute the distance between two discrete vector values\n");
  }
};

void wrap_discreteValVect(nb::module_ &m) { discreteValVec_wrapper::wrap(m); }
