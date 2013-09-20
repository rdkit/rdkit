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
        
        return p12.dotProduct(p32) / (dist1 * dist2);
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
        double res = 0.0;
        
        if (isLinear) {
          res = 143.9325 * ka * (1.0 + cosTheta);
        }
        else {
          res = 0.043844 * ka / 2.0 * angle * angle * (1.0 + cb * angle);
        }

        return res;
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

      this->d_theta0 = mmffAngleParams->theta0;
      this->d_ka = mmffAngleParams->ka;
    }

    double AngleBendContrib::getEnergy(double *pos) const {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");

      double dist1 = this->dp_forceField->distance(this->d_at1Idx, this->d_at2Idx, pos);
      double dist2 = this->dp_forceField->distance(this->d_at2Idx, this->d_at3Idx, pos);
      double res = 0.0;
      
      RDGeom::Point3D p1(pos[3 * this->d_at1Idx],
			  pos[3 * this->d_at1Idx + 1],
			  pos[3 * this->d_at1Idx + 2]);
      RDGeom::Point3D p2(pos[3 * this->d_at2Idx],
			  pos[3 * this->d_at2Idx + 1],
			  pos[3 * this->d_at2Idx + 2]);
      RDGeom::Point3D p3(pos[3 * this->d_at3Idx],
			  pos[3 * this->d_at3Idx + 1],
			  pos[3 * this->d_at3Idx + 2]);
      
      return Utils::calcAngleBendEnergy(this->d_theta0,
        this->d_ka, this->d_isLinear, Utils::calcCosTheta(p1, p2, p3, dist1, dist2));
    }

    void AngleBendContrib::getGrad(double *pos, double *grad) const
    {
      PRECONDITION(dp_forceField, "no owner");
      PRECONDITION(pos, "bad vector");
      PRECONDITION(grad, "bad vector");
      double dist1 = this->dp_forceField->distance(this->d_at1Idx, this->d_at2Idx, pos);
      double dist2 = this->dp_forceField->distance(this->d_at2Idx, this->d_at3Idx, pos);

      RDGeom::Point3D p1(pos[3 * this->d_at1Idx],
        pos[3 * this->d_at1Idx + 1], pos[3 * this->d_at1Idx + 2]);
      RDGeom::Point3D p2(pos[3 * this->d_at2Idx],
        pos[3 * this->d_at2Idx + 1], pos[3 * this->d_at2Idx + 2]);
      RDGeom::Point3D p3(pos[3 * this->d_at3Idx],
        pos[3 * this->d_at3Idx + 1], pos[3 * this->d_at3Idx + 2]);
      double *g1 = &(grad[3 * this->d_at1Idx]);
      double *g2 = &(grad[3 * this->d_at2Idx]);
      double *g3 = &(grad[3 * this->d_at3Idx]);
      RDGeom::Point3D p12 = (p1 - p2) / dist1;
      RDGeom::Point3D p32 = (p3 - p2) / dist2;
      double cosTheta = p12.dotProduct(p32);
      double sinTheta = std::max(sqrt(1.0 - cosTheta * cosTheta), 1.0e-8);

      // use the chain rule:
      // dE/dx = dE/dTheta * dTheta/dx

      // dE/dTheta is independent of cartesians:
      double angleTerm = RAD2DEG * acos(cosTheta) - this->d_theta0;
      double const cb = -0.006981317;
        
      double dE_dTheta = (this->d_isLinear ? -143.9325 * this->d_ka * sinTheta
        : RAD2DEG * 0.043844 * this->d_ka * angleTerm * (1.0 + 1.5 * cb * angleTerm));
    
      // -------
      // dTheta/dx is trickier:
      double dCos_dS1 = 1.0 / dist1 * (p32.x - cosTheta * p12.x);
      double dCos_dS2 = 1.0 / dist1 * (p32.y - cosTheta * p12.y);
      double dCos_dS3 = 1.0 / dist1 * (p32.z - cosTheta * p12.z);

      double dCos_dS4 = 1.0 / dist2 * (p12.x - cosTheta * p32.x);
      double dCos_dS5 = 1.0 / dist2 * (p12.y - cosTheta * p32.y);
      double dCos_dS6 = 1.0 / dist2 * (p12.z - cosTheta * p32.z);


      g1[0] += dE_dTheta * dCos_dS1 / (-sinTheta);
      g1[1] += dE_dTheta * dCos_dS2 / (-sinTheta);
      g1[2] += dE_dTheta * dCos_dS3 / (-sinTheta);

      g2[0] += dE_dTheta * (-dCos_dS1 - dCos_dS4) / (-sinTheta);
      g2[1] += dE_dTheta * (-dCos_dS2 - dCos_dS5) / (-sinTheta);
      g2[2] += dE_dTheta * (-dCos_dS3 - dCos_dS6) / (-sinTheta);
    
      g3[0] += dE_dTheta * dCos_dS4 / (-sinTheta);
      g3[1] += dE_dTheta * dCos_dS5 / (-sinTheta);
      g3[2] += dE_dTheta * dCos_dS6 / (-sinTheta);
    }
  }
}
