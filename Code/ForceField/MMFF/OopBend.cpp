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
#include "OopBend.h"
#include "Params.h"
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>

namespace ForceFields {
  namespace MMFF {
    namespace Utils {
      double calcOopChi(const RDGeom::Point3D &iPoint, const RDGeom::Point3D &jPoint,
        const RDGeom::Point3D &kPoint, const RDGeom::Point3D &lPoint)
      {
        RDGeom::Point3D rJI = iPoint - jPoint;
        RDGeom::Point3D rJK = kPoint - jPoint;
        RDGeom::Point3D rJL = lPoint - jPoint;
        rJI /= rJI.length();
        rJK /= rJK.length();
        rJL /= rJL.length();

        RDGeom::Point3D n = rJI.crossProduct(rJK);
        n /= n.length();
        double sinChi = n.dotProduct(rJL);
        clipToOne(sinChi);

        return RAD2DEG * asin(sinChi);
      }
      
      double calcOopBendForceConstant(const MMFFOop *mmffOopParams)
      {
        PRECONDITION(mmffOopParams, "no OOP parameters");

        return mmffOopParams->koop;
      }

      double calcOopBendEnergy(const double chi, const double koop)
      {
        double const c2 = 0.043844;
        return (0.5 * c2 * koop * chi * chi);
      }
    } // end of namespace Utils


    OopBendContrib::OopBendContrib(ForceField *owner,
      unsigned int idx1, unsigned int idx2, unsigned int idx3, unsigned int idx4,
      const MMFFOop *mmffOopParams)
    {
      PRECONDITION(owner, "bad owner");
      PRECONDITION(mmffOopParams, "no OOP parameters");
      PRECONDITION((idx1 != idx2) && (idx1 != idx3) && (idx1 != idx4)
        && (idx2 != idx3) && (idx2 != idx4) && (idx3 != idx4), "degenerate points");
      RANGE_CHECK(0, idx1, owner->positions().size() - 1);
      RANGE_CHECK(0, idx2, owner->positions().size() - 1);
      RANGE_CHECK(0, idx3, owner->positions().size() - 1);
      RANGE_CHECK(0, idx4, owner->positions().size() - 1);

      dp_forceField = owner;
      d_at1Idx = idx1;
      d_at2Idx = idx2;
      d_at3Idx = idx3;
      d_at4Idx = idx4;
      d_koop = mmffOopParams->koop;
    }

  
    double OopBendContrib::getEnergy(double *pos) const
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

      return Utils::calcOopBendEnergy
        (Utils::calcOopChi(p1, p2, p3, p4), d_koop);
    }

    void OopBendContrib::getGrad(double *pos, double *grad) const
    {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");
      PRECONDITION(grad,"bad vector");

      RDGeom::Point3D iPoint(pos[3 * d_at1Idx],
        pos[3 * d_at1Idx + 1], pos[3 * d_at1Idx + 2]);
      RDGeom::Point3D jPoint(pos[3 * d_at2Idx],
        pos[3 * d_at2Idx + 1], pos[3 * d_at2Idx + 2]);
      RDGeom::Point3D kPoint(pos[3 * d_at3Idx],
        pos[3 * d_at3Idx + 1], pos[3 * d_at3Idx + 2]);
      RDGeom::Point3D lPoint(pos[3 * d_at4Idx],
        pos[3 * d_at4Idx + 1], pos[3 * d_at4Idx + 2]);
      double *g1 = &(grad[3 * d_at1Idx]);
      double *g2 = &(grad[3 * d_at2Idx]);
      double *g3 = &(grad[3 * d_at3Idx]);
      double *g4 = &(grad[3 * d_at4Idx]);

      RDGeom::Point3D rJI = iPoint - jPoint;
      RDGeom::Point3D rJK = kPoint - jPoint;
      RDGeom::Point3D rJL = lPoint - jPoint;
      double dJI = rJI.length();
      double dJK = rJK.length();
      double dJL = rJL.length();
      if (isDoubleZero(dJI) || isDoubleZero(dJK) || isDoubleZero(dJL)) {
        return;
      }
      rJI /= dJI;
      rJK /= dJK;
      rJL /= dJL;
      
      RDGeom::Point3D n = (-rJI).crossProduct(rJK);
      n /= n.length();
      double const c2 = 0.043844;
      double sinChi = rJL.dotProduct(n);
      clipToOne(sinChi);
      double cosChiSq = 1.0 - sinChi * sinChi;
      double cosChi = std::max(((cosChiSq > 0.0) ? sqrt(cosChiSq) : 0.0), 1.0e-8);
      double chi = RAD2DEG * asin(sinChi);
      double cosTheta = rJI.dotProduct(rJK);
      clipToOne(cosTheta);
      double sinThetaSq = std::max(1.0 - cosTheta * cosTheta, 1.0e-8);
      double sinTheta = std::max(((sinThetaSq > 0.0) ? sqrt(sinThetaSq) : 0.0), 1.0e-8);
      
      double dE_dChi = RAD2DEG * c2 * d_koop * chi;
      RDGeom::Point3D t1 = rJL.crossProduct(rJK);
      RDGeom::Point3D t2 = rJI.crossProduct(rJL);
      RDGeom::Point3D t3 = rJK.crossProduct(rJI);
      double term1 = cosChi * sinTheta;
      double term2 = sinChi / (cosChi * sinThetaSq);
      double tg1[3] = { (t1.x / term1 - (rJI.x - rJK.x * cosTheta) * term2) / dJI,
        (t1.y / term1 - (rJI.y - rJK.y * cosTheta) * term2) / dJI,
        (t1.z / term1 - (rJI.z - rJK.z * cosTheta) * term2) / dJI };
      double tg3[3] = { (t2.x / term1 - (rJK.x - rJI.x * cosTheta) * term2) / dJK,
        (t2.y / term1 - (rJK.y - rJI.y * cosTheta) * term2) / dJK,
        (t2.z / term1 - (rJK.z - rJI.z * cosTheta) * term2) / dJK };
      double tg4[3] = { (t3.x / term1 - rJL.x * sinChi / cosChi) / dJL,
        (t3.y / term1 - rJL.y * sinChi / cosChi) / dJL,
        (t3.z / term1 - rJL.z * sinChi / cosChi) / dJL };
      for (unsigned int i = 0; i < 3; ++i) {
        g1[i] += dE_dChi * tg1[i];
        g2[i] += -dE_dChi * (tg1[i] + tg3[i] + tg4[i]);
        g3[i] += dE_dChi * tg3[i];
        g4[i] += dE_dChi * tg4[i];
      }
    }
  }
}
