// $Id$
//
// Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "point.h"
//#include <Numerics/Vector.h>

namespace RDGeom {
  double computeSignedDihedralAngle(const Point3D &pt1, const Point3D &pt2,
                              const Point3D &pt3, const Point3D &pt4) {
    Point3D begEndVec = pt3 - pt2;
    Point3D begNbrVec = pt1 - pt2;
    Point3D crs1 = begNbrVec.crossProduct(begEndVec);
    
    Point3D endNbrVec = pt4 - pt3;
    Point3D crs2 = endNbrVec.crossProduct(begEndVec);

    double ang = crs1.angleTo(crs2);

    // now calculate the sign:
    Point3D crs3 = crs1.crossProduct(crs2);
    double dot = crs3.dotProduct(begEndVec);
    if(dot<0.0) ang*=-1;
    
    return ang;
  }
  double computeDihedralAngle(const Point3D &pt1, const Point3D &pt2,
                              const Point3D &pt3, const Point3D &pt4) {
    Point3D begEndVec = pt3 - pt2;
    Point3D begNbrVec = pt1 - pt2;
    Point3D crs1 = begNbrVec.crossProduct(begEndVec);
    
    Point3D endNbrVec = pt4 - pt3;
    Point3D crs2 = endNbrVec.crossProduct(begEndVec);

    double ang = crs1.angleTo(crs2);
    return ang;
  }


std::ostream & operator<<(std::ostream& target, const Point &pt){
  for (unsigned int di = 0; di < pt.dimension(); ++di) {
    target << pt[di] << " ";
  }
  
  return target;
}

Point3D operator+ (const Point3D& p1, const Point3D& p2) {
  Point3D res;
  res.x = p1.x + p2.x;
  res.y = p1.y + p2.y;
  res.z = p1.z + p2.z;
  return res;
}

Point3D operator- (const Point3D& p1, const Point3D& p2) {
  Point3D res;
  res.x = p1.x - p2.x;
  res.y = p1.y - p2.y;
  res.z = p1.z - p2.z;
  return res;
}

Point3D operator* (const Point3D& p1, double v) {
  Point3D res;
  res.x = p1.x * v;
  res.y = p1.y * v;
  res.z = p1.z * v;
  return res;
}

Point3D operator/ (const Point3D& p1, double v) {
  Point3D res;
  res.x = p1.x / v;
  res.y = p1.y / v;
  res.z = p1.z / v;
  return res;
}

Point2D operator+ (const Point2D& p1, const Point2D& p2) {
  Point2D res;
  res.x = p1.x + p2.x;
  res.y = p1.y + p2.y;
  return res;
}

Point2D operator- (const Point2D& p1, const Point2D& p2) {
  Point2D res;
  res.x = p1.x - p2.x;
  res.y = p1.y - p2.y;
  return res;
}

Point2D operator* (const Point2D& p1, double v) {
  Point2D res;
  res.x = p1.x * v;
  res.y = p1.y * v;
  return res;
}

Point2D operator/ (const Point2D& p1, double v) {
  Point2D res;
  res.x = p1.x / v;
  res.y = p1.y / v;
  return res;
}


PointND operator+ (const PointND& p1, const PointND& p2) {
  unsigned int dim;
  if(p1.dimension()<p2.dimension()) {
    dim=p1.dimension();
  } else {
    dim=p1.dimension();
  }
  PointND res(dim);
  for(unsigned int i=0;i<dim;++i){
    res[i] = p1[i]+p2[i];
  }
  return res;
}
PointND operator- (const PointND& p1, const PointND& p2) {
  unsigned int dim;
  if(p1.dimension()<p2.dimension()) {
    dim=p1.dimension();
  } else {
    dim=p1.dimension();
  }
  PointND res(dim);
  for(unsigned int i=0;i<dim;++i){
    res[i] = p1[i]-p2[i];
  }
  return res;
}

PointND operator* (const PointND& p1, double v) {
  PointND res(p1.dimension());
  for(unsigned int i=0;i<p1.dimension();++i){
    res[i] = p1[i]*v;
  }
  return res;
}

PointND operator/ (const PointND& p1, double v) {
  PointND res(p1.dimension());
  for(unsigned int i=0;i<p1.dimension();++i){
    res[i] = p1[i]/v;
  }
  return res;
}



}

