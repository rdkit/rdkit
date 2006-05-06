// $Id: Point.cpp 5116 2006-04-06 18:18:37Z glandrum $
//
//  Copyright (C) 2005 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <Geometry/point.h>
#include <boost/python.hpp>
#include <RDBoost/Wrap.h>
namespace python = boost::python;


namespace RDGeom {
  
  std::string Point3Ddoc = "A class to represent a three-dimensional point";
  std::string Point2Ddoc = "A class to represent a two-dimensional point";
  
  double point3Ddist(const Point3D &pt1, const Point3D &pt2) {
    Point3D tpt(pt1);
    tpt -= pt2;
    return tpt.length();
  }

  double point3dGetLen(const Point3D &self){
    return 3;
  }
  double point3dGetItem(const Point3D &self,int idx){
    switch(idx){
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
  
  double point2dGetLen(const Point2D &self){
    return 2;
  }
  double point2dGetItem(const Point2D &self,int idx){
    switch(idx){
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
  
  
  struct Point3D_wrapper {
    static void wrap() {
      python::class_<Point3D>("Point3D", Point3Ddoc.c_str(),
                              python::init<>("Default Constructor"))
        .def(python::init<double, double, double>())
        .def_readwrite("x", &Point3D::x)
        .def_readwrite("y", &Point3D::y)
        .def_readwrite("z", &Point3D::z)
        .def("__getitem__", point3dGetItem)
        .def("__len__", point3dGetLen)
        .def("__iadd__", &Point3D::operator+=,
             python::return_value_policy<python::copy_non_const_reference>(),
             "Addition to another point")
        .def("__isub__", &Point3D::operator-=,
             python::return_value_policy<python::copy_non_const_reference>(),
             "Vector difference")
        .def(python::self - python::self)
        .def(python::self -= python::self)
        .def(python::self + python::self)
        .def(python::self += python::self)
        .def("__imul__", &Point3D::operator*=,
             python::return_value_policy<python::copy_non_const_reference>(),
             "Scalar multiplication")
        .def("__idiv__", &Point3D::operator/=,
             python::return_value_policy<python::copy_non_const_reference>(),
             "Scalar division")
        .def("Normalize", &Point3D::normalize,
             "Normalize the vector (using L2 norm)")
        .def("Length", &Point3D::length,
             "Length of the vector")
        .def("Distance", point3Ddist,
             "Distance from this point to another point")
        .def("LengthSq", &Point3D::lengthSq,
             "Square of the length")
        .def("DotProduct", &Point3D::dotProduct,
             "Dot product with another point")
        .def("AngleTo", &Point3D::angleTo,
             "determines the angle between a vector to this point (between 0 and PI)")
        .def("SignedAngleTo", &Point3D::signedAngleTo,
             "determines the signed angle between a vector to this point (between 0 and 2*PI)")
        .def("DirectionVector", &Point3D::directionVector,
             "return a normalized direction vector from this point to another")
        .def("CrossProduct", &Point3D::crossProduct,
             "Get the cross product between two points")
        ;

      
      python::class_<Point2D>("Point2D", Point2Ddoc.c_str(),
                              python::init<>("Default Constructor"))
        .def(python::init<double, double>())
        .def_readwrite("x", &Point2D::x)
        .def_readwrite("y", &Point2D::y)
        .def("__getitem__", point2dGetItem)
        .def("__len__", point2dGetLen)
        .def(python::self - python::self)
        .def(python::self -= python::self)
        .def(python::self + python::self)
        .def(python::self += python::self)
        .def("__imul__", &Point2D::operator*=,
             python::return_value_policy<python::copy_non_const_reference>(),
             "Scalar multiplication")
        .def("__idiv__", &Point2D::operator/=,
             python::return_value_policy<python::copy_non_const_reference>(),
             "Scalar division")
        .def("Normalize", &Point2D::normalize,
             "Normalize the vector (using L2 norm)")
        .def("Length", &Point2D::length,
             "Length of the vector")
        .def("LengthSq", &Point2D::lengthSq,
             "Square of the length")
        .def("DotProduct", &Point2D::dotProduct,
             "Dot product with another point")
        .def("AngleTo", &Point2D::angleTo,
             "determines the angle between a vector to this point (between 0 and PI)")
        .def("SignedAngleTo", &Point2D::signedAngleTo,
             "determines the signed angle between a vector to this point (between 0 and 2*PI)")
        .def("DirectionVector", &Point2D::directionVector,
             "return a normalized direction vector from this point to another")
        ;
    }
  };
}

void wrap_point3D() {
  RDGeom::Point3D_wrapper::wrap();
}
