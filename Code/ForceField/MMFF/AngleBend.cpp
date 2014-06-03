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
#include "BondStretch.h"
#include "Params.h"
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>


namespace ForceFields {
  namespace MMFF {
    namespace Utils {

      double calcAngleRestValue(const MMFFAngle *mmffAngleParams)
      {
        PRECONDITION(mmffAngleParams, "angle parameters not found");
        
        return mmffAngleParams->theta0;
      }
      
      double calcCosTheta(RDGeom::Point3D p1, RDGeom::Point3D p2,
        RDGeom::Point3D p3, double dist1, double dist2)
      {
        RDGeom::Point3D p12 = p1 - p2;
        RDGeom::Point3D p32 = p3 - p2;
        double cosTheta = p12.dotProduct(p32) / (dist1 * dist2);
        clipToOne(cosTheta);

        return cosTheta;
      }
      
      double calcAngleForceConstant(const MMFFAngle *mmffAngleParams)
      {
        PRECONDITION(mmffAngleParams, "angle parameters not found");
        
        return mmffAngleParams->ka;
      }
      
      double calcAngleBendEnergy(const double theta0,
        const double ka, bool isLinear, const double cosTheta)
      {
        double angle = RAD2DEG * acos(cosTheta) - theta0;
        double const cb = -0.006981317;
        double const c1 = 143.9325;
        double const c2 = 0.043844;
        double res = 0.0;
        
        if (isLinear) {
          res = c1 * ka * (1.0 + cosTheta);
        }
        else {
          res = 0.5 * c2 * ka * angle * angle * (1.0 + cb * angle);
        }

        return res;
      }
      
      void calcAngleBendGrad(RDGeom::Point3D *r, double *dist,
        double **g, double &dE_dTheta, double &cosTheta, double &sinTheta)
      {
        // -------
        // dTheta/dx is trickier:
        double dCos_dS[6] = {
          1.0 / dist[0] * (r[1].x - cosTheta * r[0].x),
          1.0 / dist[0] * (r[1].y - cosTheta * r[0].y),
          1.0 / dist[0] * (r[1].z - cosTheta * r[0].z),
          1.0 / dist[1] * (r[0].x - cosTheta * r[1].x),
          1.0 / dist[1] * (r[0].y - cosTheta * r[1].y),
          1.0 / dist[1] * (r[0].z - cosTheta * r[1].z)
        };

        g[0][0] += dE_dTheta * dCos_dS[0] / (-sinTheta);
        g[0][1] += dE_dTheta * dCos_dS[1] / (-sinTheta);
        g[0][2] += dE_dTheta * dCos_dS[2] / (-sinTheta);

        g[1][0] += dE_dTheta * (-dCos_dS[0] - dCos_dS[3]) / (-sinTheta);
        g[1][1] += dE_dTheta * (-dCos_dS[1] - dCos_dS[4]) / (-sinTheta);
        g[1][2] += dE_dTheta * (-dCos_dS[2] - dCos_dS[5]) / (-sinTheta);
      
        g[2][0] += dE_dTheta * dCos_dS[3] / (-sinTheta);
        g[2][1] += dE_dTheta * dCos_dS[4] / (-sinTheta);
        g[2][2] += dE_dTheta * dCos_dS[5] / (-sinTheta);
      }
    } // end of namespace Utils
    
    AngleBendContrib::AngleBendContrib(ForceField *owner,
      unsigned int idx1, unsigned int idx2, unsigned int idx3,
      const MMFFAngle *mmffAngleParams, const MMFFProp *mmffPropParamsCentralAtom)
    {
      PRECONDITION(owner, "bad owner");
      PRECONDITION(((idx1 != idx2) && (idx2 != idx3) && (idx1 != idx3)),"degenerate points");
      RANGE_CHECK(0, idx1, owner->positions().size() - 1);
      RANGE_CHECK(0, idx2, owner->positions().size() - 1);
      RANGE_CHECK(0, idx3, owner->positions().size() - 1);
        
      dp_forceField = owner;
      d_at1Idx = idx1;
      d_at2Idx = idx2;
      d_at3Idx = idx3;
      d_isLinear = (mmffPropParamsCentralAtom->linh ? true : false);

      d_theta0 = mmffAngleParams->theta0;
      d_ka = mmffAngleParams->ka;
    }

    double AngleBendContrib::getEnergy(double *pos) const {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");

      double dist1 = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
      double dist2 = dp_forceField->distance(d_at2Idx, d_at3Idx, pos);

      RDGeom::Point3D p1(pos[3 * d_at1Idx],
			  pos[3 * d_at1Idx + 1],
			  pos[3 * d_at1Idx + 2]);
      RDGeom::Point3D p2(pos[3 * d_at2Idx],
			  pos[3 * d_at2Idx + 1],
			  pos[3 * d_at2Idx + 2]);
      RDGeom::Point3D p3(pos[3 * d_at3Idx],
			  pos[3 * d_at3Idx + 1],
			  pos[3 * d_at3Idx + 2]);
      
      return Utils::calcAngleBendEnergy(d_theta0,
        d_ka, d_isLinear, Utils::calcCosTheta(p1, p2, p3, dist1, dist2));
    }

    void AngleBendContrib::getGrad(double *pos, double *grad) const
    {
      PRECONDITION(dp_forceField, "no owner");
      PRECONDITION(pos, "bad vector");
      PRECONDITION(grad, "bad vector");
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
      double angleTerm = RAD2DEG * acos(cosTheta) - d_theta0;
      double const cb = -0.006981317;
      double const c1 = 143.9325;
      double const c2 = 0.043844;
        
      double dE_dTheta = (d_isLinear ? -c1 * d_ka * sinTheta
        : RAD2DEG * c2 * d_ka * angleTerm * (1.0 + 1.5 * cb * angleTerm));
    
      Utils::calcAngleBendGrad(r, dist, g, dE_dTheta, cosTheta, sinTheta);
    }
  }
}
