// $Id$
//
//  Copyright (C) 2005 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>

#include <RDBoost/Wrap.h>
#include <Geometry/point.h>
namespace python = boost::python;

namespace {
struct Point3D_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(RDGeom::Point3D const &pt) {
    return python::make_tuple(pt.x, pt.y, pt.z);
  }
};
struct Point2D_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(RDGeom::Point2D const &pt) {
    return python::make_tuple(pt.x, pt.y);
  }
};
struct PointND_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(RDGeom::PointND const &pt) {
    return python::make_tuple(pt.dimension());
  }
  static python::tuple getstate(RDGeom::PointND const &pt) {
    python::list res;
    for (unsigned int i = 0; i < pt.dimension(); ++i) {
      res.append(pt[i]);
    }
    return python::tuple(res);
  }
  static void setstate(RDGeom::PointND &pt, python::tuple state) {
    unsigned int sz = python::extract<unsigned int>(state.attr("__len__")());
    for (unsigned int i = 0; i < sz; ++i) {
      pt[i] = python::extract<double>(state[i]);
    }
  }
};
}  // namespace

namespace RDGeom {

std::string Point3Ddoc =
    "A class to represent a three-dimensional point\n\
The x, y, and z coordinates can be read and written using either attributes\n\
(i.e. pt.x = 4) or indexing (i.e. pt[0] = 4).\n";
std::string Point2Ddoc = "A class to represent a two-dimensional point";
std::string PointNDdoc = "A class to represent an N-dimensional point";

double point3Ddist(const Point3D &pt1, const Point3D &pt2) {
  Point3D tpt(pt1);
  tpt -= pt2;
  return tpt.length();
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

struct Point_wrapper {
  static void wrap() {
    python::class_<Point3D>(
        "Point3D", Point3Ddoc.c_str(),
        python::init<>(python::args("self"), "Default Constructor"))
        .def(python::init<double, double, double>(
            python::args("self", "xv", "yv", "zv")))
        .def_readwrite("x", &Point3D::x)
        .def_readwrite("y", &Point3D::y)
        .def_readwrite("z", &Point3D::z)
        .def("__getitem__", point3dGetItem, python::args("self", "idx"))
        .def("__len__", &Point3D::dimension, python::args("self"))
        .def("__iadd__", &Point3D::operator+=,
             python::return_value_policy<python::copy_non_const_reference>(),
             python::args("self", "other"), "Addition to another point")
        .def("__isub__", &Point3D::operator-=,
             python::return_value_policy<python::copy_non_const_reference>(),
             python::args("self", "other"),
             "Vector difference")
        .def(python::self - python::self)
#ifdef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wself-assign-overloaded"
#endif
        .def(python::self -=
             python::self)  // clang warns incorrectly on this construct
#ifdef __clang__
#pragma GCC diagnostic pop
#endif
        .def(python::self + python::self)
        .def(python::self += python::self)
        .def(python::self * double())
        .def(python::self / double())
        .def("__imul__", &Point3D::operator*=,
             python::return_value_policy<python::copy_non_const_reference>(),
             python::args("self", "scale"), "Scalar multiplication")
        .def("__idiv__", &Point3D::operator/=,
             python::return_value_policy<python::copy_non_const_reference>(),
             python::args("self", "scale"), "Scalar division")
        .def("Normalize", &Point3D::normalize, python::args("self"),
             "Normalize the vector (using L2 norm)")
        .def("Length", &Point3D::length, python::args("self"),
             "Length of the vector")
        .def("Distance", point3Ddist, python::args("self", "pt2"),
             "Distance from this point to another point")
        .def("LengthSq", &Point3D::lengthSq, python::args("self"),
             "Square of the length")
        .def("DotProduct", &Point3D::dotProduct, python::args("self", "other"),
             "Dot product with another point")
        .def("AngleTo", &Point3D::angleTo, python::args("self", "other"),
             "determines the angle between a vector to this point (between 0 "
             "and PI)")
        .def("SignedAngleTo", &Point3D::signedAngleTo,
             python::args("self", "other"),
             "determines the signed angle between a vector to this point "
             "(between 0 and 2*PI)")
        .def("DirectionVector", &Point3D::directionVector,
             python::args("self", "other"),
             "return a normalized direction vector from this point to another")
        .def("CrossProduct", &Point3D::crossProduct,
             python::args("self", "other"),
             "Get the cross product between two points")

        .def_pickle(Point3D_pickle_suite());

    python::class_<Point2D>(
        "Point2D", Point2Ddoc.c_str(),
        python::init<>(python::args("self"), "Default Constructor"))
        .def(python::init<double, double>(python::args("self", "xv", "yv")))
        .def(python::init<const Point3D &>(
            (python::args("self"), python::args("other")),
            "construct from a Point3D (ignoring the z component)"))
        .def_readwrite("x", &Point2D::x)
        .def_readwrite("y", &Point2D::y)
        .def("__getitem__", point2dGetItem, python::args("self", "idx"))
        .def("__len__", &Point2D::dimension, python::args("self"))
        .def(python::self - python::self)
#ifdef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wself-assign-overloaded"
#endif
        .def(python::self -=
             python::self)  // clang warns incorrectly on this construct
#ifdef __clang__
#pragma GCC diagnostic pop
#endif
        .def(python::self + python::self)
        .def(python::self += python::self)
        .def(python::self * double())
        .def(python::self / double())
        .def("__imul__", &Point2D::operator*=,
             python::return_value_policy<python::copy_non_const_reference>(),
             python::args("self", "scale"), "Scalar multiplication")
        .def("__idiv__", &Point2D::operator/=,
             python::return_value_policy<python::copy_non_const_reference>(),
             python::args("self", "scale"), "Scalar division")
        .def("Normalize", &Point2D::normalize, python::args("self"),
             "Normalize the vector (using L2 norm)")
        .def("Length", &Point2D::length, python::args("self"),
             "Length of the vector")
        .def("LengthSq", &Point2D::lengthSq, python::args("self"),
             "Square of the length")
        .def("DotProduct", &Point2D::dotProduct, python::args("self", "other"),
             "Dot product with another point")
        .def("AngleTo", &Point2D::angleTo, python::args("self", "other"),
             "determines the angle between a vector to this point (between 0 "
             "and PI)")
        .def("SignedAngleTo", &Point2D::signedAngleTo,
             python::args("self", "other"),
             "determines the signed angle between a vector to this point "
             "(between 0 and 2*PI)")
        .def("DirectionVector", &Point2D::directionVector,
             python::args("self", "other"),
             "return a normalized direction vector from this point to another")

        .def_pickle(Point2D_pickle_suite());

    python::class_<PointND>(
        "PointND", PointNDdoc.c_str(),
        python::init<unsigned int>(python::args("self", "dim")))
        .def("__getitem__", pointNdGetItem, python::args("self", "idx"))
        .def("__setitem__", pointNdSetItem, python::args("self", "idx", "val"))
        .def("__len__", &PointND::dimension, python::args("self"))
        .def("__iadd__", &PointND::operator+=,
             python::return_value_policy<python::copy_non_const_reference>(),
             python::args("self", "other"), "Addition to another point")
        .def("__isub__", &PointND::operator-=,
             python::return_value_policy<python::copy_non_const_reference>(),
             python::args("self", "other"),
             "Vector difference")
        .def(python::self - python::self)
#ifdef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wself-assign-overloaded"
#endif
        .def(python::self -=
             python::self)  // clang warns incorrectly on this construct
#ifdef __clang__
#pragma GCC diagnostic pop
#endif
        .def(python::self + python::self)
        .def(python::self += python::self)
        .def(python::self * double())
        .def(python::self / double())
        .def("__imul__", &PointND::operator*=,
             python::return_value_policy<python::copy_non_const_reference>(),
             python::args("self", "scale"), "Scalar multiplication")
        .def("__idiv__", &PointND::operator/=,
             python::return_value_policy<python::copy_non_const_reference>(),
             python::args("self", "scale"), "Scalar division")
        .def("Normalize", &PointND::normalize, python::args("self"),
             "Normalize the vector (using L2 norm)")
        .def("Length", &PointND::length, python::args("self"),
             "Length of the vector")
        .def("Distance", point3Ddist, python::args("self", "pt2"),
             "Distance from this point to another point")
        .def("LengthSq", &PointND::lengthSq, python::args("self"),
             "Square of the length")
        .def("DotProduct", &PointND::dotProduct, python::args("self", "other"),
             "Dot product with another point")
        .def("AngleTo", &PointND::angleTo, python::args("self", "other"),
             "determines the angle between a vector to this point (between 0 "
             "and PI)")
        .def("DirectionVector", &PointND::directionVector,
             python::args("self", "other"),
             "return a normalized direction vector from this point to another")
        //.def("SignedAngleTo", &PointND::signedAngleTo,
        //     "determines the signed angle between a vector to this point
        //     (between 0 and 2*PI)")
        //.def("CrossProduct", &PointND::crossProduct,
        //     "Get the cross product between two points")
        .def_pickle(PointND_pickle_suite());

    python::def(
        "ComputeDihedralAngle", computeDihedralAngle,
        python::args("pt1", "pt2", "pt3", "pt4"),
        "calculates the dihedral angle determined by four Point3D objects");
    python::def("ComputeSignedDihedralAngle", computeSignedDihedralAngle,
                python::args("pt1", "pt2", "pt3", "pt4"),
                "calculates the signed dihedral angle determined by four "
                "Point3D objects");
  }
};
}  // namespace RDGeom

void wrap_point() { RDGeom::Point_wrapper::wrap(); }
