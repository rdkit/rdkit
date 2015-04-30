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
#include "TorsionAngle.h"
#include "TorsionConstraint.h"
#include "Params.h"
#include <cmath>
#include <boost/math/special_functions/round.hpp>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>

namespace ForceFields {
  namespace MMFF {
    void _pretreatDihedrals(double &minDihedralDeg, double &maxDihedralDeg)
    {
      if (minDihedralDeg < 0.0) minDihedralDeg += 360.0;
      if (maxDihedralDeg < 0.0) maxDihedralDeg += 360.0;
      minDihedralDeg = fmod(minDihedralDeg, 360.0);
      maxDihedralDeg = fmod(maxDihedralDeg, 360.0);
      if (maxDihedralDeg < minDihedralDeg) maxDihedralDeg += 360.0;
    }
    
    TorsionConstraintContrib::TorsionConstraintContrib(ForceField *owner,
      unsigned int idx1, unsigned int idx2, unsigned int idx3,
      unsigned int idx4, double minDihedralDeg, double maxDihedralDeg,
      double forceConst)
    {
      PRECONDITION(owner,"bad owner");
      RANGE_CHECK(0, idx1, owner->positions().size() - 1);
      RANGE_CHECK(0, idx2, owner->positions().size() - 1);
      RANGE_CHECK(0, idx3, owner->positions().size() - 1);
      RANGE_CHECK(0, idx4, owner->positions().size() - 1);
      PRECONDITION((!(maxDihedralDeg < minDihedralDeg))
        && ((maxDihedralDeg - minDihedralDeg) < 360.0), "bad bounds");
      _pretreatDihedrals(minDihedralDeg, maxDihedralDeg);

      dp_forceField = owner;
      d_at1Idx = idx1;
      d_at2Idx = idx2;
      d_at3Idx = idx3;
      d_at4Idx = idx4;
      d_minDihedralDeg = minDihedralDeg;
      d_maxDihedralDeg = maxDihedralDeg;
      d_forceConstant = forceConst;
    }

    TorsionConstraintContrib::TorsionConstraintContrib(ForceField *owner,
      unsigned int idx1, unsigned int idx2, unsigned int idx3,
      unsigned int idx4, bool relative, double minDihedralDeg,
      double maxDihedralDeg, double forceConst)
    {
      PRECONDITION(owner,"bad owner");
      const RDGeom::PointPtrVect &pos = owner->positions();
      RANGE_CHECK(0, idx1, pos.size() - 1);
      RANGE_CHECK(0, idx2, pos.size() - 1);
      RANGE_CHECK(0, idx3, pos.size() - 1);
      RANGE_CHECK(0, idx4, pos.size() - 1);
      PRECONDITION((!(maxDihedralDeg < minDihedralDeg))
        && ((maxDihedralDeg - minDihedralDeg) < 360.0), "bad bounds");

      double dihedral = 0.0;
      if (relative) {
        RDGeom::Point3D p1 = *((RDGeom::Point3D *)pos[idx1]);
        RDGeom::Point3D p2 = *((RDGeom::Point3D *)pos[idx2]);
        RDGeom::Point3D p3 = *((RDGeom::Point3D *)pos[idx3]);
        RDGeom::Point3D p4 = *((RDGeom::Point3D *)pos[idx4]);
        RDGeom::Point3D r12 = p2 - p1;
        RDGeom::Point3D r23 = p3 - p2;
        RDGeom::Point3D r34 = p4 - p3;
        
        RDGeom::Point3D n123 = r12.crossProduct(r23);
        double nIJKSqLength = n123.lengthSq();
        RDGeom::Point3D n234 = r23.crossProduct(r34);
        double nJKLSqLength = n234.lengthSq();
        RDGeom::Point3D m = n123.crossProduct(r23);
        // we want a signed dihedral, that's why we use atan2 instead of acos
        dihedral = RAD2DEG
          * (-atan2(m.dotProduct(n234) / sqrt(nJKLSqLength * m.lengthSq()),
          n123.dotProduct(n234) / sqrt(nIJKSqLength * nJKLSqLength)));
      }
      dp_forceField = owner;
      d_at1Idx = idx1;
      d_at2Idx = idx2;
      d_at3Idx = idx3;
      d_at4Idx = idx4;
      minDihedralDeg += dihedral;
      maxDihedralDeg += dihedral;
      _pretreatDihedrals(minDihedralDeg, maxDihedralDeg);
      d_minDihedralDeg = minDihedralDeg;
      d_maxDihedralDeg = maxDihedralDeg;
      d_forceConstant = forceConst;
    }

    double TorsionConstraintContrib::getEnergy(double *pos) const
    {
      PRECONDITION(dp_forceField, "no owner");
      PRECONDITION(pos, "bad vector");
 
      RDGeom::Point3D p1(pos[3 * d_at1Idx],
			  pos[3 * d_at1Idx + 1], pos[3 * d_at1Idx + 2]);
      RDGeom::Point3D p2(pos[3 * d_at2Idx],
			  pos[3 * d_at2Idx + 1], pos[3 * d_at2Idx + 2]);
      RDGeom::Point3D p3(pos[3 * d_at3Idx],
			  pos[3 * d_at3Idx + 1], pos[3 * d_at3Idx + 2]);
      RDGeom::Point3D p4(pos[3 * d_at4Idx],
			  pos[3 * d_at4Idx + 1], pos[3 * d_at4Idx + 2]);
      
      RDGeom::Point3D r1 = p1 - p2;
      RDGeom::Point3D r2 = p3 - p2;
      RDGeom::Point3D r3 = p2 - p3;
      RDGeom::Point3D r4 = p4 - p3;
      RDGeom::Point3D t1 = r1.crossProduct(r2);
      RDGeom::Point3D t2 = r3.crossProduct(r4);
      double d1 = std::max(t1.length(), 0.0);
      double d2 = std::max(t2.length(), 0.0);
      t1 /= d1;
      t2 /= d2;
      
      RDGeom::Point3D n123 = (-r1).crossProduct(r2);
      double n123SqLength = n123.lengthSq();
      RDGeom::Point3D n234 = r2.crossProduct(r4);
      double n234SqLength = n234.lengthSq();
      RDGeom::Point3D m = n123.crossProduct(r2);
      // we want a signed dihedral, that's why we use atan2 instead of acos
      double dihedral = RAD2DEG * (-atan2(m.dotProduct(n234) / sqrt(n234SqLength * m.lengthSq()),
        n123.dotProduct(n234) / sqrt(n123SqLength * n234SqLength)));
      if (dihedral < 0.0) dihedral += 360.0;
      double dihedralTerm = 0.0;
      if (dihedral < d_minDihedralDeg) {
        dihedralTerm = dihedral - d_minDihedralDeg;
      }
      else if (dihedral > d_maxDihedralDeg) {
        dihedralTerm = dihedral - d_maxDihedralDeg;
      }
      double const c = 0.5 * DEG2RAD * DEG2RAD;
      double res = c * d_forceConstant * dihedralTerm * dihedralTerm;

      return res;
    }
    
    void TorsionConstraintContrib::getGrad(double *pos, double *grad) const
    {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");
      PRECONDITION(grad,"bad vector");

      RDGeom::Point3D p1(pos[3 * d_at1Idx],
        pos[3 * d_at1Idx + 1], pos[3 * d_at1Idx + 2]);
      RDGeom::Point3D p2(pos[3 * d_at2Idx],
        pos[3 * d_at2Idx + 1], pos[3 * d_at2Idx + 2]);
      RDGeom::Point3D p3(pos[3 * d_at3Idx],
        pos[3 * d_at3Idx + 1], pos[3 * d_at3Idx + 2]);
      RDGeom::Point3D p4(pos[3 * d_at4Idx],
        pos[3 * d_at4Idx + 1], pos[3 * d_at4Idx + 2]);
      double *g[4] = {
        &(grad[3 * d_at1Idx]),
        &(grad[3 * d_at2Idx]),
        &(grad[3 * d_at3Idx]),
        &(grad[3 * d_at4Idx])
      };
      
      RDGeom::Point3D r[4] = {
        p1 - p2,
        p3 - p2,
        p2 - p3,
        p4 - p3
      };
      RDGeom::Point3D t[2] = {
        r[0].crossProduct(r[1]),
        r[2].crossProduct(r[3])
      };
      double d[2] = {
        t[0].length(),
        t[1].length()
      };
      if (isDoubleZero(d[0]) || isDoubleZero(d[1])) {
        return;
      }
      t[0] /= d[0];
      t[1] /= d[1];
      double cosPhi = t[0].dotProduct(t[1]);
      clipToOne(cosPhi);
      double sinPhiSq = 1.0 - cosPhi * cosPhi;
      double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);
      // dE/dPhi is independent of cartesians:
      
      RDGeom::Point3D n123 = (-r[0]).crossProduct(r[1]);
      double n123SqLength = n123.lengthSq();
      RDGeom::Point3D n234 = r[1].crossProduct(r[3]);
      double n234SqLength = n234.lengthSq();
      RDGeom::Point3D m = n123.crossProduct(r[1]);
      // we want a signed dihedral, that's why we use atan2 instead of acos
      double dihedral = RAD2DEG * (-atan2(m.dotProduct(n234) / sqrt(n234SqLength * m.lengthSq()),
        n123.dotProduct(n234) / sqrt(n123SqLength * n234SqLength)));
      if (dihedral < 0.0) dihedral += 360.0;
      //double dihedral = RAD2DEG * acos(cosPhi);
      double dihedralTerm = 0.0;
      if (dihedral < d_minDihedralDeg) {
        dihedralTerm = dihedral - d_minDihedralDeg;
      }
      else if (dihedral > d_maxDihedralDeg) {
        dihedralTerm = dihedral - d_maxDihedralDeg;
      }
      if (dihedral > 180.0) dihedralTerm = -dihedralTerm;
      double dE_dPhi = DEG2RAD * d_forceConstant * dihedralTerm;
      
      // FIX: use a tolerance here
      // this is hacky, but it's per the
      // recommendation from Niketic and Rasmussen:
      double sinTerm = -dE_dPhi * (isDoubleZero(sinPhi)
        ? (1.0 / cosPhi) : (1.0 / sinPhi));
      Utils::calcTorsionGrad(r, t, d, g, sinTerm, cosPhi);
    }
  }
}  
