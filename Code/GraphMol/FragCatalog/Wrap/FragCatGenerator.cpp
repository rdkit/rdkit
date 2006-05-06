// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <boost/python.hpp>
#include <DataStructs/BitVects.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FragCatalog/FragCatGenerator.h>

namespace python = boost::python;
namespace RDKit{
  struct fragcatgen_wrapper {
    static void wrap() {
      python::class_<FragCatGenerator>("FragCatGenerator", python::init<>())
	.def("AddFragsFromMol", &FragCatGenerator::addFragsFromMol)
      ;
    };
  }; // end of struct
} // end of namespace

void wrap_fragcatgen() {
  RDKit::fragcatgen_wrapper::wrap();
}
