//  RealValueVect.cpp
//  Created on: Apr 14, 2014
//  Author: hahnda6
//
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <boost/python.hpp>

#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <RDBoost/PySequenceHolder.h>
#include <DataStructs/RealValueVect.h>

using namespace RDKit;

struct rvv_pickle_suite : python::pickle_suite
{
  static python::tuple
  getinitargs(const RealValueVect& self)
  {
    return python::make_tuple(self.toString());
  };
};

std::string realValVectDoc="A container class for storing real\n\
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

    python::class_<RealValueVect>("RealValueVect",
                                  realValVectDoc.c_str(),
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
              .def_pickle(rvv_pickle_suite())
                                  ;
    python::def("ComputeL1Norm", computeL1Norm,
                "Compute the distance between two real vector values\n");
  }
};

void wrap_realValVect() {
  realValVec_wrapper::wrap();
}




