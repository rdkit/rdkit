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
#include <nanobind/stl/function.h>

using namespace RDGeom;
namespace nb = nanobind;
using namespace nb::literals;

namespace {
Point2D add(const Point2D &p1, const Point2D &p2) {
    auto res = p1 + p2;
    return res;
}
// for some reason if I use a reference for the argument here, the
// object is not modified when passed to python. Maybe it's copied?
void foo(std::function<void(Point2D *)> func, Point2D *pt){
    func(pt);
}
}  // namespace

NB_MODULE(rdGeometryTest, m) {
  m.doc() = "Testing cross-module objects";
  m.def("add", &add);
  m.def("foo", &foo);
}
