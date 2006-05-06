// $Id: FragFPGenerator.cpp 4968 2006-02-18 00:27:15Z glandrum $
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <boost/python.hpp>
#include <DataStructs/BitVects.h>

#include <GraphMol/FragCatalog/FragFPGenerator.h>

namespace python = boost::python;
namespace RDKit{
  struct fragFPgen_wrapper {
    static void wrap() {
      python::class_<FragFPGenerator>("FragFPGenerator", python::init<>())
        .def("GetFPForMol", &FragFPGenerator::getFPForMol,
             python::return_value_policy<python::manage_new_object>())
        ;
    };
  };
}

void wrap_fragFPgen() {
  RDKit::fragFPgen_wrapper::wrap();
}

