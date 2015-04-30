// $Id$
//
//  Copyright (C) 2013 Paolo Tosco
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "AngleBend.h"
#include "AngleConstraint.h"
#include "Params.h"
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>

namespace ForceFields {
  namespace UFF {
    void _pretreatAngles(double &minAngleDeg, double &maxAngleDeg)
    {
      minAngleDeg = fmod(minAngleDeg, 360.0);
      maxAngleDeg = fmod(maxAngleDeg, 360.0);
      if (minAngleDeg > 180.0) minAngleDeg -= 360.0;
      if (maxAngleDeg > 180.0) maxAngleDeg -= 360.0;
      if ((minAngleDeg < 0.0) && (!(maxAngleDeg < 0.0))) {
        maxAngleDeg = std::max(fabs(maxAngleDeg), fabs(minAngleDeg));
        minAngleDeg = 0.0;
      }
      minAngleDeg = fabs(minAngleDeg);
      maxAngleDeg = fabs(maxAngleDeg);
      if (minAngleDeg > maxAngleDeg) {
        double t = minAngleDeg;
        minAngleDeg = maxAngleDeg;
        maxAngleDeg = t;
      }
    }
    
    AngleConstraintContrib::AngleConstraintContrib(ForceField *owner,
      unsigned int idx1, unsigned int idx2, unsigned int idx3,
      double minAngleDeg, double maxAngleDeg, double forceConst)
    {
      PRECONDITION(owner,"bad owner");
      RANGE_CHECK(0, idx1, owner->positions().size() - 1);
      RANGE_CHECK(0, idx2, owner->positions().size() - 1);
      RANGE_CHECK(0, idx3, owner->positions().size() - 1);
      PRECONDITION(maxAngleDeg >= minAngleDeg, "allowedDeltaDeg must be >= 0.0");
      _pretreatAngles(minAngleDeg, maxAngleDeg);

      dp_forceField = owner;
      d_at1Idx = idx1;
      d_at2Idx = idx2;
      d_at3Idx = idx3;
      d_minAngleDeg = minAngleDeg;
      d_maxAngleDeg = maxAngleDeg;
      d_forceConstant = forceConst;
    }

    AngleConstraintContrib::AngleConstraintContrib(ForceField *owner,
      unsigned int idx1, unsigned int idx2, unsigned int idx3,
      bool relative, double minAngleDeg, double maxAngleDeg,
      double forceConst)
    {
      PRECONDITION(owner,"bad owner");
      const RDGeom::PointPtrVect &pos = owner->positions();
      RANGE_CHECK(0, idx1, pos.size() - 1);
      RANGE_CHECK(0, idx2, pos.size() - 1);
      RANGE_CHECK(0, idx3, pos.size() - 1);
      PRECONDITION(maxAngleDeg >= minAngleDeg, "allowedDeltaDeg must be >= 0.0");

      double angle = 0.0;
      if (relative) {
        RDGeom::Point3D p1 = *((RDGeom::Point3D *)pos[idx1]);
        RDGeom::Point3D p2 = *((RDGeom::Point3D *)pos[idx2]);
        RDGeom::Point3D p3 = *((RDGeom::Point3D *)pos[idx3]);
        double dist1 = (p1 - p2).length();
        double dist2 = (p3 - p2).length();
        RDGeom::Point3D p12 = (p1 - p2) / dist1;
        RDGeom::Point3D p32 = (p3 - p2) / dist2;
        double cosTheta = p12.dotProduct(p32);
        clipToOne(cosTheta);
        angle = RAD2DEG * acos(cosTheta);
      }
      dp_forceField = owner;
      d_at1Idx = idx1;
      d_at2Idx = idx2;
      d_at3Idx = idx3;
      minAngleDeg += angle;
      maxAngleDeg += angle;
      _pretreatAngles(minAngleDeg, maxAngleDeg);
      d_minAngleDeg = minAngleDeg;
      d_maxAngleDeg = maxAngleDeg;
      d_forceConstant = forceConst;
    }

    double AngleConstraintContrib::getEnergy(double *pos) const
    {
      PRECONDITION(dp_forceField, "no owner");
      PRECONDITION(pos, "bad vector");
 
      double dist1 = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
      double dist2 = dp_forceField->distance(d_at2Idx, d_at3Idx, pos);
      
      RDGeom::Point3D p1(pos[3 * d_at1Idx],
			  pos[3 * d_at1Idx + 1], pos[3 * d_at1Idx + 2]);
      RDGeom::Point3D p2(pos[3 * d_at2Idx],
			  pos[3 * d_at2Idx + 1], pos[3 * d_at2Idx + 2]);
      RDGeom::Point3D p3(pos[3 * d_at3Idx],
			  pos[3 * d_at3Idx + 1], pos[3 * d_at3Idx + 2]);
      
      RDGeom::Point3D p12 = (p1 - p2) / dist1;
      RDGeom::Point3D p32 = (p3 - p2) / dist2;
      double cosTheta = p12.dotProduct(p32);
      clipToOne(cosTheta);
      double angle = RAD2DEG * acos(cosTheta);
      double angleTerm = 0.0;
      if (angle < d_minAngleDeg) {
        angleTerm = angle - d_minAngleDeg;
      }
      else if (angle > d_maxAngleDeg) {
        angleTerm = angle - d_maxAngleDeg;
      }
      double const c = 0.5 * DEG2RAD * DEG2RAD;
      double res = c * d_forceConstant * angleTerm * angleTerm;

      return res;
    }
    
    void AngleConstraintContrib::getGrad(double *pos, double *grad) const
    {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");
      PRECONDITION(grad,"bad vector");

      double dist[2] = {
        dp_forceField->distance(d_at1Idx, d_at2Idx, pos),
        dp_forceField->distance(d_at2Idx, d_at3Idx, pos)
      };

      RDGeom::Point3D p1(pos[3 * d_at1Idx],
        pos[3 * d_at1Idx + 1], pos[3 * d_at1Idx + 2]);
      RDGeom::Point3D p2(pos[3 * d_at2Idx],
        pos[3 * d_at2Idx + 1], pos[3 * d_at2Idx + 2]);
      RDGeom::Point3D p3(pos[3 * d_at3Idx],
        pos[3 * d_at3Idx + 1], pos[3 * d_at3Idx + 2]);
      double *g[3] = {
        &(grad[3 * d_at1Idx]),
        &(grad[3 * d_at2Idx]),
        &(grad[3 * d_at3Idx])
      };
      RDGeom::Point3D r[2] = {
        (p1 - p2) / dist[0],
        (p3 - p2) / dist[1]
      };
      double cosTheta = r[0].dotProduct(r[1]);
      clipToOne(cosTheta);
      double sinThetaSq = 1.0 - cosTheta * cosTheta;
      double sinTheta = std::max(((sinThetaSq > 0.0) ? sqrt(sinThetaSq) : 0.0), 1.0e-8);

      // use the chain rule:
      // dE/dx = dE/dTheta * dTheta/dx

      // dE/dTheta is independent of cartesians:
      double angle = RAD2DEG * acos(cosTheta);
      double angleTerm = 0.0;
      if (angle < d_minAngleDeg) {
        angleTerm = angle - d_minAngleDeg;
      }
      else if (angle > d_maxAngleDeg) {
        angleTerm = angle - d_maxAngleDeg;
      }
        
      double dE_dTheta = DEG2RAD * d_forceConstant * angleTerm;
    
      Utils::calcAngleBendGrad(r, dist, g, dE_dTheta, cosTheta, sinTheta);
    }
  }
}  
