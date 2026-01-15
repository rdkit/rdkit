//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <Geometry/point.h>
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

using namespace RDGeom;
namespace nb = nanobind;
using namespace nb::literals;

namespace {
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
double point3dGetItem(const Point3D &self, int idx) {
  switch (idx) {
    case 0:
    case -3:
      return self.x;
      break;
    case 1:
    case -2:
      return self.y;
      break;
    case 2:
    case -1:
      return self.z;
      break;
    default:
      throw IndexErrorException(idx);
  }
}
double pointNdGetItem(const PointND &self, int idx) {
  if (idx >= static_cast<int>(self.dimension()) ||
      idx < -1 * static_cast<int>(self.dimension())) {
    throw IndexErrorException(idx);
  }
  if (idx < 0) {
    idx = self.dimension() + idx;
  }
  return self[idx];
}
double pointNdSetItem(PointND &self, int idx, double val) {
  if (idx >= static_cast<int>(self.dimension()) ||
      idx < -1 * static_cast<int>(self.dimension())) {
    throw IndexErrorException(idx);
  }
  if (idx < 0) {
    idx = self.dimension() + idx;
  }
  self[idx] = val;
  return val;
}

}  // namespace

void wrap_point(nb::module_ &m) {
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

  nb::class_<Point3D>(m, "Point3D")
      .def(nb::init<>())
      .def(nb::init<double, double, double>(), "x"_a, "y"_a, "z"_a)
      .def(nb::init<const Point3D &>())
      .def_rw("x", &Point3D::x)
      .def_rw("y", &Point3D::y)
      .def_rw("z", &Point3D::z)
      .def("__getitem__", point3dGetItem)
      .def("__len__", &Point3D::dimension)
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
      .def("Normalize", &Point3D::normalize,
           "Normalize the vector (using L2 norm)")
      .def("Length", &Point3D::length, "Length of the vector")
      .def("LengthSq", &Point3D::lengthSq, "Square of the length")
      .def("DotProduct", &Point3D::dotProduct, "Dot product with another point")
      .def("CrossProduct", &Point3D::crossProduct,
           "Get the cross product between two points")
      .def("AngleTo", &Point3D::angleTo,
           "determines the angle between a vector to this point (between 0 "
           "and PI)")
      .def("SignedAngleTo", &Point3D::signedAngleTo,
           "determines the signed angle between a vector to this point "
           "(between 0 and 2*PI)")
      .def("DirectionVector", &Point3D::directionVector,
           "return a normalized direction vector from this point to another")
      .def("__getstate__",
           [](const Point3D &pt) { return std::make_tuple(pt.x, pt.y, pt.z); })
      .def("__setstate__",
           [](Point3D &pt, const std::tuple<double, double, double> &state) {
             new (&pt) Point3D(std::get<0>(state), std::get<1>(state),
                               std::get<2>(state));
           });

  nb::class_<PointND>(m, "PointND")
      .def(nb::init<unsigned int>(), "dim"_a)
      // .def(nb::init<double, double>(), "x"_a, "y"_a)
      // .def(nb::init<const Point3D &>())
      .def("__getitem__", pointNdGetItem)
      .def("__setitem__", pointNdSetItem)
      .def("__len__", &PointND::dimension)
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
      .def("Normalize", &PointND::normalize,
           "Normalize the vector (using L2 norm)")
      .def("Length", &PointND::length, "Length of the vector")
      .def("LengthSq", &PointND::lengthSq, "Square of the length")
      .def("DotProduct", &PointND::dotProduct, "Dot product with another point")
      .def("AngleTo", &PointND::angleTo,
           "determines the angle between a vector to this point (between 0 "
           "and PI)")
      // .def("SignedAngleTo", &PointND::signedAngleTo,
      //      "determines the signed angle between a vector to this point "
      //      "(between 0 and 2*PI)")
      .def("DirectionVector", &PointND::directionVector,
           "return a normalized direction vector from this point to another")
      .def("__getstate__",
           [](const PointND &pt) {
             std::vector<double> state;
             for (unsigned int i = 0; i < pt.dimension(); ++i) {
               state.push_back(pt[i]);
             }
             return state;
           })
      .def("__setstate__", [](PointND &pt, const std::vector<double> &state) {
        new (&pt) PointND(state);
      });
  ;

  m.def("ComputeDihedralAngle", computeDihedralAngle, "pt1"_a, "pt2"_a, "pt3"_a,
        "pt4"_a,
        "calculates the dihedral angle determined by four Point3D objects");
  m.def("ComputeSignedDihedralAngle", computeSignedDihedralAngle, "pt1"_a,
        "pt2"_a, "pt3"_a, "pt4"_a,
        "calculates the signed dihedral angle determined by four "
        "Point3D objects");
}
