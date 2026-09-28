//
//  Copyright (C) 2003-2026 Rational Discovery LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/FragCatalog/FragCatGenerator.h>

namespace nb = nanobind;
using namespace nb::literals;

void wrap_fragcatgen(nb::module_ &m) {
  nb::class_<RDKit::FragCatGenerator>(m, "FragCatGenerator")
      .def(nb::init<>())
      .def("AddFragsFromMol", &RDKit::FragCatGenerator::addFragsFromMol,
           "mol"_a, "fcat"_a);
}
