//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// #include <RDBoost/Wrap.h>
// #include <RDBoost/python.h>

#include <Geometry/point.h>
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/tuple.h>

namespace nb = nanobind;
void wrap_point(nb::module_ &m);
void wrap_uniformGrid(nb::module_ &m);
void wrap_uniformrealvalueGrid(nb::module_ &m);

NB_MODULE(rdGeometry, m) {
  m.doc() = "Module containing geometry objects like points, grids, etc.";
  wrap_point(m);
  wrap_uniformGrid(m);
  wrap_uniformrealvalueGrid(m);
}
