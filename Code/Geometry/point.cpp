// $Id: point.cpp 4955 2006-02-17 23:37:53Z glandrum $
//
// Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "point.h"

std::ostream & operator<<(std::ostream& target, const RDGeom::Point3D &pt){
  target << pt.x << " " << pt.y << " " << pt.z;
  return target;
}

std::ostream & operator<<(std::ostream& target, const RDGeom::Point2D &pt){
  target << pt.x << " " << pt.y;
  return target;
}

RDGeom::Point3D operator+ (const RDGeom::Point3D& p1, const RDGeom::Point3D& p2) {
  RDGeom::Point3D res;
  res.x = p1.x + p2.x;
  res.y = p1.y + p2.y;
  res.z = p1.z + p2.z;
  return res;
}

RDGeom::Point3D operator- (const RDGeom::Point3D& p1, const RDGeom::Point3D& p2) {
  RDGeom::Point3D res;
  res.x = p1.x - p2.x;
  res.y = p1.y - p2.y;
  res.z = p1.z - p2.z;
  return res;
}

RDGeom::Point3D operator* (const RDGeom::Point3D& p1, const double v) {
  RDGeom::Point3D res;
  res.x = p1.x * v;
  res.y = p1.y * v;
  res.z = p1.z * v;
  return res;
}

RDGeom::Point3D operator/ (const RDGeom::Point3D& p1, const double v) {
  RDGeom::Point3D res;
  res.x = p1.x / v;
  res.y = p1.y / v;
  res.z = p1.z / v;
  return res;
}

RDGeom::Point2D operator+ (const RDGeom::Point2D& p1, const RDGeom::Point2D& p2) {
  RDGeom::Point2D res;
  res.x = p1.x + p2.x;
  res.y = p1.y + p2.y;
  return res;
}

RDGeom::Point2D operator- (const RDGeom::Point2D& p1, const RDGeom::Point2D& p2) {
  RDGeom::Point2D res;
  res.x = p1.x - p2.x;
  res.y = p1.y - p2.y;
  return res;
}

RDGeom::Point2D operator* (const RDGeom::Point2D& p1, const double v) {
  RDGeom::Point2D res;
  res.x = p1.x * v;
  res.y = p1.y * v;
  return res;
}

RDGeom::Point2D operator/ (const RDGeom::Point2D& p1, const double v) {
  RDGeom::Point2D res;
  res.x = p1.x / v;
  res.y = p1.y / v;
  return res;
}

namespace RDGeom {
  double computeDihedralAngle(Point3D pt1, Point3D pt2, Point3D pt3, Point3D pt4) {
    Point3D begEndVec = pt3 - pt2;
    Point3D begNbrVec = pt1 - pt2;
    Point3D crs1 = begNbrVec.crossProduct(begEndVec);
    
    begEndVec *= -1.0;
    Point3D endNbrVec = pt4 - pt3;
    Point3D crs2 = begEndVec.crossProduct(endNbrVec);

    double ang = crs1.angleTo(crs2);
    return ang;
  }
}

