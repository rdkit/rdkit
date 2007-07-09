// $Id$
//
// Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "point.h"
//#include <Numerics/Vector.h>

std::ostream & operator<<(std::ostream& target, const RDGeom::Point &pt){
  for (unsigned int di = 0; di < pt.dimension(); ++di) {
    target << pt[di] << " ";
  }
  
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


RDGeom::PointND operator+ (const RDGeom::PointND& p1, const RDGeom::PointND& p2) {
  unsigned int dim;
  if(p1.dimension()<p2.dimension()) {
    dim=p1.dimension();
  } else {
    dim=p1.dimension();
  }
  RDGeom::PointND res(dim);
  for(unsigned int i=0;i<dim;++i){
    res[i] = p1[i]+p2[i];
  }
  return res;
}
RDGeom::PointND operator- (const RDGeom::PointND& p1, const RDGeom::PointND& p2) {
  unsigned int dim;
  if(p1.dimension()<p2.dimension()) {
    dim=p1.dimension();
  } else {
    dim=p1.dimension();
  }
  RDGeom::PointND res(dim);
  for(unsigned int i=0;i<dim;++i){
    res[i] = p1[i]-p2[i];
  }
  return res;
}

RDGeom::PointND operator* (const RDGeom::PointND& p1, const double v) {
  RDGeom::PointND res(p1.dimension());
  for(unsigned int i=0;i<p1.dimension();++i){
    res[i] = p1[i]*v;
  }
  return res;
}

RDGeom::PointND operator/ (const RDGeom::PointND& p1, const double v) {
  RDGeom::PointND res(p1.dimension());
  for(unsigned int i=0;i<p1.dimension();++i){
    res[i] = p1[i]/v;
  }
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

