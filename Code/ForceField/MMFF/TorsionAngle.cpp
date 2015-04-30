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
#include "Params.h"
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>

namespace ForceFields {
  namespace MMFF {
    namespace Utils {
      double calcTorsionCosPhi(const RDGeom::Point3D &iPoint,
        const RDGeom::Point3D &jPoint, const RDGeom::Point3D &kPoint,
        const RDGeom::Point3D &lPoint)
      {
        RDGeom::Point3D r1 = iPoint - jPoint;
        RDGeom::Point3D r2 = kPoint - jPoint;
        RDGeom::Point3D r3 = jPoint - kPoint;
        RDGeom::Point3D r4 = lPoint - kPoint;
        RDGeom::Point3D t1 = r1.crossProduct(r2);
        RDGeom::Point3D t2 = r3.crossProduct(r4);
        double cosPhi = t1.dotProduct(t2) / (t1.length() * t2.length());
        clipToOne(cosPhi);

        return cosPhi;
      }
      
      boost::tuple<double, double, double>
        calcTorsionForceConstant(const MMFFTor *mmffTorParams)
      {
        return boost::make_tuple(mmffTorParams->V1,
          mmffTorParams->V2, mmffTorParams->V3);
      }

      double calcTorsionEnergy(const double V1,
        const double V2, const double V3, const double cosPhi)
      {
        double cos2Phi = 2.0 * cosPhi * cosPhi - 1.0;
        double cos3Phi = cosPhi * (2.0 * cos2Phi - 1.0);

        return (0.5 * (V1 * (1.0 + cosPhi)
          + V2 * (1.0 - cos2Phi) + V3 * (1.0 + cos3Phi)));
      }
      
      void calcTorsionGrad(RDGeom::Point3D *r, RDGeom::Point3D *t,
        double *d, double **g, double &sinTerm, double &cosPhi)
      {
        // -------
        // dTheta/dx is trickier:
        double dCos_dT[6] = {
          1.0 / d[0] * (t[1].x - cosPhi * t[0].x),
          1.0 / d[0] * (t[1].y - cosPhi * t[0].y),
          1.0 / d[0] * (t[1].z - cosPhi * t[0].z),
          1.0 / d[1] * (t[0].x - cosPhi * t[1].x),
          1.0 / d[1] * (t[0].y - cosPhi * t[1].y),
          1.0 / d[1] * (t[0].z - cosPhi * t[1].z)
        };
      
        g[0][0] += sinTerm * (dCos_dT[2] * r[1].y - dCos_dT[1] * r[1].z);
        g[0][1] += sinTerm * (dCos_dT[0] * r[1].z - dCos_dT[2] * r[1].x);
        g[0][2] += sinTerm * (dCos_dT[1] * r[1].x - dCos_dT[0] * r[1].y);

        g[1][0] += sinTerm * (dCos_dT[1] * (r[1].z - r[0].z)
          + dCos_dT[2] * (r[0].y - r[1].y) + dCos_dT[4] * (-r[3].z) + dCos_dT[5] * (r[3].y));
        g[1][1] += sinTerm * (dCos_dT[0] * (r[0].z - r[1].z)
          + dCos_dT[2] * (r[1].x - r[0].x) + dCos_dT[3] * (r[3].z) + dCos_dT[5] * (-r[3].x));
        g[1][2] += sinTerm * (dCos_dT[0] * (r[1].y - r[0].y)
          + dCos_dT[1] * (r[0].x - r[1].x) + dCos_dT[3] * (-r[3].y) + dCos_dT[4] * (r[3].x));

        g[2][0] += sinTerm * (dCos_dT[1] * (r[0].z) + dCos_dT[2] * (-r[0].y) +
          dCos_dT[4] * (r[3].z - r[2].z) + dCos_dT[5] * (r[2].y - r[3].y));
        g[2][1] += sinTerm * (dCos_dT[0] * (-r[0].z) + dCos_dT[2] * (r[0].x) +
          dCos_dT[3] * (r[2].z - r[3].z) + dCos_dT[5] * (r[3].x - r[2].x));
        g[2][2] += sinTerm * (dCos_dT[0] * (r[0].y) + dCos_dT[1] * (-r[0].x) +
          dCos_dT[3] * (r[3].y - r[2].y) + dCos_dT[4] * (r[2].x - r[3].x));

        g[3][0] += sinTerm * (dCos_dT[4] * r[2].z - dCos_dT[5] * r[2].y);
        g[3][1] += sinTerm * (dCos_dT[5] * r[2].x - dCos_dT[3] * r[2].z);
        g[3][2] += sinTerm * (dCos_dT[3] * r[2].y - dCos_dT[4] * r[2].x);
      }
    }


    TorsionAngleContrib::TorsionAngleContrib(ForceField *owner,
      unsigned int idx1, unsigned int idx2, unsigned int idx3, unsigned int idx4,
      const MMFFTor *mmffTorParams)
    {
      PRECONDITION(owner, "bad owner");
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
      d_V1 = mmffTorParams->V1;
      d_V2 = mmffTorParams->V2;
      d_V3 = mmffTorParams->V3;
    }

  
    double TorsionAngleContrib::getEnergy(double *pos) const
    {
      PRECONDITION(dp_forceField, "no owner");
      PRECONDITION(pos, "bad vector");

      RDGeom::Point3D iPoint(pos[3 * d_at1Idx],
        pos[3 * d_at1Idx + 1], pos[3 * d_at1Idx + 2]);
      RDGeom::Point3D jPoint(pos[3 * d_at2Idx],
        pos[3 * d_at2Idx + 1], pos[3 * d_at2Idx + 2]);
      RDGeom::Point3D kPoint(pos[3 * d_at3Idx],
        pos[3 * d_at3Idx + 1], pos[3 * d_at3Idx + 2]);
      RDGeom::Point3D lPoint(pos[3 * d_at4Idx],
        pos[3 * d_at4Idx + 1], pos[3 * d_at4Idx + 2]);

      return Utils::calcTorsionEnergy(d_V1, d_V2, d_V3,
        Utils::calcTorsionCosPhi(iPoint, jPoint, kPoint, lPoint));
    }

    void TorsionAngleContrib::getGrad(double *pos, double *grad) const
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
      double *g[4] = {
        &(grad[3 * d_at1Idx]),
        &(grad[3 * d_at2Idx]),
        &(grad[3 * d_at3Idx]),
        &(grad[3 * d_at4Idx])
      };

      RDGeom::Point3D r[4] = {
        iPoint - jPoint,
        kPoint - jPoint,
        jPoint - kPoint,
        lPoint - kPoint
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
      double sin2Phi = 2.0 * sinPhi * cosPhi;
      double sin3Phi = 3.0 * sinPhi - 4.0 * sinPhi * sinPhiSq;
      // dE/dPhi is independent of cartesians:
      double dE_dPhi = 0.5 * (-(d_V1) * sinPhi + 2.0
        * d_V2 * sin2Phi - 3.0 * d_V3 * sin3Phi);
#if 0
      if(dE_dPhi!=dE_dPhi){
        std::cout << "\tNaN in Torsion("<<d_at1Idx<<","<<d_at2Idx<<","<<d_at3Idx<<","<<d_at4Idx<<")"<< std::endl;
        std::cout << "sin: " << sinPhi << std::endl;
        std::cout << "cos: " << cosPhi << std::endl;
      } 
      
#endif
      // FIX: use a tolerance here
      // this is hacky, but it's per the
      // recommendation from Niketic and Rasmussen:
      double sinTerm = -dE_dPhi * (isDoubleZero(sinPhi)
        ? (1.0 / cosPhi) : (1.0 / sinPhi));

      Utils::calcTorsionGrad(r, t, d, g, sinTerm, cosPhi);
    }
  }
}
