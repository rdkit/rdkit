// $Id$
//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <boost/python.hpp>

#include <RDGeneral/types.h>
#include<RDGeneral/Invariant.h>
#include <RDBoost/PySequenceHolder.h>
#include <DataStructs/DiscreteValueVect.h>

using namespace RDKit;

struct dvv_pickle_suite : python::pickle_suite
{
  static python::tuple
  getinitargs(const DiscreteValueVect& self)
  {
    std::string res=self.toString();
    python::object retval = python::object(python::handle<>(PyBytes_FromStringAndSize(res.c_str(),res.length())));
    return python::make_tuple(retval);
  };
};

std::string disValVectDoc="A container class for storing unsigned integer\n\
values within a particular range.\n\
\n\
The length of the vector and type of its elements (determines the maxium value\n\
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
  static void wrap() {
    python::enum_<DiscreteValueVect::DiscreteValueType>("DiscreteValueType")
      .value("ONEBITVALUE", DiscreteValueVect::ONEBITVALUE)
      .value("TWOBITVALUE", DiscreteValueVect::TWOBITVALUE)
      .value("FOURBITVALUE", DiscreteValueVect::FOURBITVALUE)
      .value("EIGHTBITVALUE", DiscreteValueVect::EIGHTBITVALUE)
      .value("SIXTEENBITVALUE", DiscreteValueVect::SIXTEENBITVALUE)
      ;

    python::class_<DiscreteValueVect>("DiscreteValueVect", 
                                      disValVectDoc.c_str(),
                                      python::init<DiscreteValueVect::DiscreteValueType, unsigned int>("Constructor"))
      .def(python::init<std::string>())
      .def("__len__", &DiscreteValueVect::getLength,
           "Get the number of entries in the vector")
      .def("__setitem__", &DiscreteValueVect::setVal,
           "Set the value at a specified location")
      .def("__getitem__", &DiscreteValueVect::getVal,
           "Get the value at a specified location")
      .def(python::self & python::self)
      .def(python::self | python::self)
      .def(python::self - python::self)
      .def(python::self -= python::self)
      .def(python::self + python::self)
      .def(python::self += python::self)

      .def("GetValueType", &DiscreteValueVect::getValueType,
           "Get the type of value stored in the vector")
      .def("GetTotalVal", &DiscreteValueVect::getTotalVal,
           "Get the sum of the values in the vector, basically L1 norm")


      .def_pickle(dvv_pickle_suite())
      ;

    python::def("ComputeL1Norm", computeL1Norm,
                "Compute the distance between two discrete vector values\n");
  }
};

void wrap_discreteValVect() {
  discreteValVec_wrapper::wrap();
}

