// $Id$
//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
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
    return python::make_tuple(self.toString());
  };
};

std::string disValVectDoc = "A container class for storing discrete integer values\n";

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

