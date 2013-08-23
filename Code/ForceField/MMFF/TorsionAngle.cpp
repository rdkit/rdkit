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

        return (t1.dotProduct(t2) / (t1.length() * t2.length()));
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

      RDGeom::Point3D iPoint(pos[3 * this->d_at1Idx],
        pos[3 * this->d_at1Idx + 1], pos[3 * this->d_at1Idx + 2]);
      RDGeom::Point3D jPoint(pos[3 * this->d_at2Idx],
        pos[3 * this->d_at2Idx + 1], pos[3 * this->d_at2Idx + 2]);
      RDGeom::Point3D kPoint(pos[3 * this->d_at3Idx],
        pos[3 * this->d_at3Idx + 1], pos[3 * this->d_at3Idx + 2]);
      RDGeom::Point3D lPoint(pos[3 * this->d_at4Idx],
        pos[3 * this->d_at4Idx + 1], pos[3 * this->d_at4Idx + 2]);

      return Utils::calcTorsionEnergy(this->d_V1, this->d_V2, this->d_V3,
        Utils::calcTorsionCosPhi(iPoint, jPoint, kPoint, lPoint));
    }

    void TorsionAngleContrib::getGrad(double *pos, double *grad) const
    {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");
      PRECONDITION(grad,"bad vector");

      RDGeom::Point3D iPoint(pos[3 * this->d_at1Idx],
        pos[3 * this->d_at1Idx + 1], pos[3 * this->d_at1Idx + 2]);
      RDGeom::Point3D jPoint(pos[3 * this->d_at2Idx],
        pos[3 * this->d_at2Idx + 1], pos[3 * this->d_at2Idx + 2]);
      RDGeom::Point3D kPoint(pos[3 * this->d_at3Idx],
        pos[3 * this->d_at3Idx + 1], pos[3 * this->d_at3Idx + 2]);
      RDGeom::Point3D lPoint(pos[3 * this->d_at4Idx],
        pos[3 * this->d_at4Idx + 1], pos[3 * this->d_at4Idx + 2]);
      double *g1 = &(grad[3 * this->d_at1Idx]);
      double *g2 = &(grad[3 * this->d_at2Idx]);
      double *g3 = &(grad[3 * this->d_at3Idx]);
      double *g4 = &(grad[3 * this->d_at4Idx]);

      RDGeom::Point3D r1 = iPoint - jPoint;
      RDGeom::Point3D r2 = kPoint - jPoint;
      RDGeom::Point3D r3 = jPoint - kPoint;
      RDGeom::Point3D r4 = lPoint - kPoint;
      RDGeom::Point3D t1 = r1.crossProduct(r2);
      RDGeom::Point3D t2 = r3.crossProduct(r4);
      double d1 = t1.length();
      t1 /= d1;
      double d2 = t2.length();
      t2 /= d2;
      if (isDoubleZero(d1) || isDoubleZero(d2)) {
        return;
      }
      double cosPhi = t1.dotProduct(t2);
      double sinPhiSq = 1.0 - cosPhi * cosPhi;
      double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);
      double sin2Phi = 2.0 * sinPhi * cosPhi;
      double sin3Phi = 3.0 * sinPhi - 4.0 * sinPhi * sinPhiSq;
      // dE/dPhi is independent of cartesians:
      double dE_dPhi = 0.5 * (-(this->d_V1) * sinPhi + 2.0
        * this->d_V2 * sin2Phi - 3.0 * this->d_V3 * sin3Phi);
#if 0
      if(dE_dPhi!=dE_dPhi){
        std::cout << "\tNaN in Torsion("<<this->d_at1Idx<<","<<this->d_at2Idx<<","<<this->d_at3Idx<<","<<this->d_at4Idx<<")"<< std::endl;
        std::cout << "sin: " << sinPhi << std::endl;
        std::cout << "cos: " << cosPhi << std::endl;
      } 
      
#endif
      
      // -------
      // dTheta/dx is trickier:
      double dCos_dT1 = 1.0 / d1 * (t2.x - cosPhi * t1.x);
      double dCos_dT2 = 1.0 / d1 * (t2.y - cosPhi * t1.y);
      double dCos_dT3 = 1.0 / d1 * (t2.z - cosPhi * t1.z);
                                                    
      double dCos_dT4 = 1.0 / d2 * (t1.x - cosPhi * t2.x);
      double dCos_dT5 = 1.0 / d2 * (t1.y - cosPhi * t2.y);
      double dCos_dT6 = 1.0 / d2 * (t1.z - cosPhi * t2.z);
    
      // FIX: use a tolerance here
      // this is hacky, but it's per the
      // recommendation from Niketic and Rasmussen:
      double sinTerm = -dE_dPhi * (isDoubleZero(sinPhi)
        ? (1.0 / cosPhi) : (1.0 / sinPhi));

      g1[0] += sinTerm * (dCos_dT3 * r2.y - dCos_dT2 * r2.z);
      g1[1] += sinTerm * (dCos_dT1 * r2.z - dCos_dT3 * r2.x);
      g1[2] += sinTerm * (dCos_dT2 * r2.x - dCos_dT1 * r2.y);

      g2[0] += sinTerm * (dCos_dT2 * (r2.z - r1.z)
        + dCos_dT3 * (r1.y - r2.y) + dCos_dT5 * (-1.0 * r4.z) + dCos_dT6 * (r4.y));
      g2[1] += sinTerm * (dCos_dT1 * (r1.z - r2.z)
        + dCos_dT3 * (r2.x - r1.x) + dCos_dT4 * (r4.z) + dCos_dT6 * (-1.0 * r4.x));
      g2[2] += sinTerm * (dCos_dT1 * (r2.y - r1.y)
        + dCos_dT2 * (r1.x - r2.x) + dCos_dT4 * (-1.0 * r4.y) + dCos_dT5 * (r4.x));

      g3[0] += sinTerm * (dCos_dT2 * (r1.z) + dCos_dT3 * (-1.0 * r1.y) +
        dCos_dT5 * (r4.z - r3.z) + dCos_dT6 * (r3.y - r4.y));
      g3[1] += sinTerm * (dCos_dT1 * (-1.0 * r1.z) + dCos_dT3 * (r1.x) +
        dCos_dT4 * (r3.z - r4.z) + dCos_dT6 * (r4.x - r3.x));
      g3[2] += sinTerm * (dCos_dT1 * (r1.y) + dCos_dT2 * (-1.0 * r1.x) +
        dCos_dT4 * (r4.y - r3.y) + dCos_dT5 * (r3.x - r4.x));

      g4[0] += sinTerm * (dCos_dT5 * r3.z - dCos_dT6 * r3.y);
      g4[1] += sinTerm * (dCos_dT6 * r3.x - dCos_dT4 * r3.z);
      g4[2] += sinTerm * (dCos_dT4 * r3.y - dCos_dT5 * r3.x);
    }
  }
}
