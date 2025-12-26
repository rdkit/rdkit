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

using namespace RDGeom;
namespace nb = nanobind;
using namespace nb::literals;

namespace {
int add(int a, int b) { return a + b; }
double point2dGetItem(const Point2D &self, int idx) {
  switch (idx) {
    case 0:
    case -2:
      return self.x;
      break;
    case 1:
    case -1:
      return self.y;
      break;
    default:
      throw IndexErrorException(idx);
  }
}

}  // namespace

NB_MODULE(rdGeometry, m) {
  m.doc() = "Module containing geometry objects like points, grids, etc.";
  m.def("add", &add);

  nb::class_<Point2D>(m, "Point2D")
      .def(nb::init<>())
      .def(nb::init<double, double>(), "x"_a, "y"_a)
      .def(nb::init<const Point3D &>())
      .def_rw("x", &Point2D::x)
      .def_rw("y", &Point2D::y)
      .def("__getitem__", point2dGetItem)
      .def("__len__", &Point2D::dimension)
      .def(nb::self - nb::self)
#ifdef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wself-assign-overloaded"
#endif
      .def(nb::self -= nb::self)  // clang warns incorrectly on this construct
#ifdef __clang__
#pragma GCC diagnostic pop
#endif
      .def(nb::self + nb::self)
      .def(nb::self += nb::self)
      .def(nb::self * double())
      .def(nb::self / double())
      .def(nb::self *= double())
      .def(nb::self /= double())
      .def("Normalize", &Point2D::normalize,
           "Normalize the vector (using L2 norm)")
      .def("Length", &Point2D::length, "Length of the vector")
      .def("LengthSq", &Point2D::lengthSq, "Square of the length")
      .def("DotProduct", &Point2D::dotProduct, "Dot product with another point")
      .def("AngleTo", &Point2D::angleTo,
           "determines the angle between a vector to this point (between 0 "
           "and PI)")
      .def("SignedAngleTo", &Point2D::signedAngleTo,
           "determines the signed angle between a vector to this point "
           "(between 0 and 2*PI)")
      .def("DirectionVector", &Point2D::directionVector,
           "return a normalized direction vector from this point to another")
      .def("__getstate__",
           [](const Point2D &pt) { return std::make_tuple(pt.x, pt.y); })
      .def("__setstate__",
           [](Point2D &pt, const std::tuple<double, double> &state) {
             new (&pt) Point2D(std::get<0>(state), std::get<1>(state));
           });
}

// namespace python = boost::python;
// void wrap_point();
// void wrap_uniformGrid();
// void wrap_uniformrealvalueGrid();

// BOOST_PYTHON_MODULE(rdGeometry) {
//   python::scope().attr("__doc__") =
//       "Module containing geometry objects like points, grids, etc\n";

//   wrap_point();
//   wrap_uniformGrid();
//   wrap_uniformrealvalueGrid();
// }
