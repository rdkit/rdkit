//  $Id$
// 
//   Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/Invariant.h>
#include <cstring>
#include "Transform.h"
#include "Transform3D.h"
#include <math.h>
#include "point.h"

namespace RDGeom {
  void Transform3D::setToIdentity() {
    unsigned int i,id;
    double *data = d_data.get();
    memset(static_cast<void *>(data),0,d_dataSize*sizeof(double));
    for (i = 0; i < DIM_3D; i++) {
      id = i*(DIM_3D+1);
      data[id] = 1.0;
    }
  }

  void Transform3D::TransformPoint(Point3D &pt) const {
    double x, y, z;
    double *data = d_data.get();
    x = data[0]*pt.x + data[1]*pt.y + data[2]*pt.z + data[3];
    y = data[4]*pt.x + data[5]*pt.y + data[6]*pt.z + data[7];
    z = data[8]*pt.x + data[9]*pt.y + data[10]*pt.z + data[11];
    pt.x = x; pt.y = y; pt.z = z;
  }
    
  void Transform3D::SetTranslation(const Point3D &move) {
    unsigned int i = DIM_3D -1;
    double *data = d_data.get();
    data[i] = move.x; i += DIM_3D;
    data[i] = move.y; i += DIM_3D;
    data[i] = move.z; i += DIM_3D;
    data[i] = 1.0;
  }

  void Transform3D::SetRotation(double angle, AxisType axis) {
    double cosT = cos(angle);
    double sinT = sin(angle);
    this->setToIdentity();
    double *data = d_data.get();
    if(axis==Z_Axis){
      data[0] = cosT;
      data[1] = -sinT;
      data[4] = sinT;
      data[5] = cosT;
    } else if(axis==X_Axis){
      data[5] = cosT;
      data[6] = -sinT;
      data[9] = sinT;
      data[10] = cosT;
    } else if(axis==Y_Axis){
      data[0] = cosT;
      data[2] = sinT;
      data[8] = -sinT;
      data[10] = cosT;
    }
  }

  void Transform3D::SetRotation(double cosT,double sinT, const Point3D &axis) {
    double t = 1-cosT;
    double X = axis.x,Y=axis.y,Z=axis.z;
    double *data = d_data.get();
    data[0] = t*X*X + cosT;
    data[1] = t*X*Y - sinT*Z;
    data[2] = t*X*Z + sinT*Y;
    
    data[4] = t*X*Y + sinT*Z;
    data[5] = t*Y*Y + cosT;
    data[6] = t*Y*Z - sinT*X;

    data[8] = t*X*Z - sinT*Y;
    data[9] = t*Y*Z + sinT*X;
    data[10] = t*Z*Z + cosT;
  }

  void Transform3D::SetRotation(double angle, const Point3D &axis) {
    this->setToIdentity();
    double c = cos(angle);
    double s = sin(angle);
    this->SetRotation(c,s,axis);
  }

  void Transform3D::SetRotationFromQuaternion(double quaternion[4]) {
    double q00 = quaternion[0]*quaternion[0];
    double q11 = quaternion[1]*quaternion[1];
    double q22 = quaternion[2]*quaternion[2];
    double q33 = quaternion[3]*quaternion[3];
    double sumSq = q00 + q11 + q22 + q33;
    
    double q01 = 2*quaternion[0]*quaternion[1];
    double q02 = 2*quaternion[0]*quaternion[2];
    double q03 = 2*quaternion[0]*quaternion[3];
    double q12 = 2*quaternion[1]*quaternion[2];
    double q13 = 2*quaternion[1]*quaternion[3];
    double q23 = 2*quaternion[2]*quaternion[3];
    double *data = d_data.get();
    data[0] = (q00 + q11 - q22 - q33)/sumSq;
    data[1] = (q12 + q03)/sumSq; 
    data[2] = (q13 - q02)/sumSq; 
    
    data[4] = (q12 - q03)/sumSq; 
    data[5] = (q00 - q11 + q22 - q33)/sumSq;
    data[6] = (q23 + q01)/sumSq; 
    
    data[8] = (q13 + q02)/sumSq; 
    data[9] = (q23 - q01)/sumSq; 
    data[10] = (q00 - q11 - q22 + q33)/sumSq;
  }

  void Transform3D::Reflect() {
    unsigned int i, j, id;
    double *data = d_data.get();
    for (i = 0; i < DIM_3D-1; i++) {
      id = i*DIM_3D;
      for (j = 0; j < DIM_3D-1; j++) {
        data[id+j] *= -1.0;
      }
    }
  }
}

RDGeom::Transform3D operator* (const RDGeom::Transform3D &t1, const RDGeom::Transform3D &t2) {
  RDGeom::Transform3D res;
  RDNumeric::multiply(t1, t2, res);
  return res;
};

RDGeom::Point3D operator* (const RDGeom::Transform3D &t, const RDGeom::Point3D &pt) {
  RDGeom::Point3D res=pt;
  t.TransformPoint(res);
  return res;
};

