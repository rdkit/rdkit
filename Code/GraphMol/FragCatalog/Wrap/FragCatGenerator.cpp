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
