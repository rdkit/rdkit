// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
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

