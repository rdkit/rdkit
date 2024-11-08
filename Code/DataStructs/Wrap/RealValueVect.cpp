//
//  Copyright (c) 2014-2024, Novartis Institutes for BioMedical Research and
//  other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <boost/python.hpp>

#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <RDBoost/PySequenceHolder.h>
#include <DataStructs/RealValueVect.h>

using namespace RDKit;

struct rvv_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(const RealValueVect &self) {
    std::string res = self.toString();
    python::object retval = python::object(
        python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
    return python::make_tuple(retval);
  };
};

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

struct realValVec_wrapper {
  static void wrap() {
    python::class_<RealValueVect>("RealValueVect", realValVectDoc.c_str(),
                                  python::init<unsigned int>("Constructor"))
        .def(python::init<std::string>())
        .def("__len__", &RealValueVect::getLength,
             "Get the number of entries in the vector")
        .def("__setitem__", &RealValueVect::setVal,
             "Set the value at a specified location")
        .def("__getitem__", &RealValueVect::getVal,
             "Get the value at a specified location")
        .def(python::self & python::self)
        .def(python::self | python::self)
        .def(python::self - python::self)
        .def(python::self -= python::self)
        .def(python::self + python::self)
        .def(python::self += python::self)
        .def("GetTotalVal", &RealValueVect::getTotalVal,
             "Get the sum of the values in the vector, basically L1 norm")
        .def_pickle(rvv_pickle_suite());
    python::def("ComputeL1Norm", computeL1Norm,
                "Compute the distance between two real vector values\n");
  }
};

void wrap_realValVect() { realValVec_wrapper::wrap(); }
