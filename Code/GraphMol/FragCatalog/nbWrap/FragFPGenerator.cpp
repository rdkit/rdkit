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

#include <DataStructs/BitVects.h>
#include <GraphMol/FragCatalog/FragFPGenerator.h>

namespace nb = nanobind;
using namespace nb::literals;

void wrap_fragFPgen(nb::module_ &m) {
  nb::class_<RDKit::FragFPGenerator>(m, "FragFPGenerator")
      .def(nb::init<>())
      .def("GetFPForMol", &RDKit::FragFPGenerator::getFPForMol,
           nb::rv_policy::take_ownership, "mol"_a, "fcat"_a);
}
