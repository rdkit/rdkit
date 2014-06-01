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

#include "BondStretch.h"
#include "Params.h"
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>

namespace ForceFields {
  namespace MMFF {
    namespace Utils {

      double calcBondRestLength(const MMFFBond *mmffBondParams)
      {
        PRECONDITION(mmffBondParams, "bond parameters not found");
        
        return mmffBondParams->r0;
      }
  
      double calcBondForceConstant(const MMFFBond *mmffBondParams)
      {
        PRECONDITION(mmffBondParams, "bond parameters not found");

        return mmffBondParams->kb;
      }

      double calcBondStretchEnergy(const double r0, const double kb, const double distance)
      {
        double distTerm = distance - r0;
        double distTerm2 = distTerm * distTerm;
        double const c1 = 143.9325;
        double const cs = -2.0;
        double const c3 = 7.0 / 12.0;
        
        return (0.5 * c1 * kb * distTerm2
          * (1.0 + cs * distTerm + c3 * cs * cs * distTerm2));
      }
    } // end of namespace Utils
  
    BondStretchContrib::BondStretchContrib(ForceField *owner,
      const unsigned int idx1, const unsigned int idx2, 
      const MMFFBond *mmffBondParams)
    {
      PRECONDITION(owner,"bad owner");
      RANGE_CHECK(0, idx1, owner->positions().size() - 1);
      RANGE_CHECK(0, idx2, owner->positions().size() - 1);

      dp_forceField = owner;
      d_at1Idx = idx1;
      d_at2Idx = idx2;
      d_r0 = mmffBondParams->r0;
      d_kb = mmffBondParams->kb;
    }

    double BondStretchContrib::getEnergy(double *pos) const
    {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");

      return Utils::calcBondStretchEnergy(d_r0, d_kb,
        dp_forceField->distance(d_at1Idx, d_at2Idx, pos));
    }

    void BondStretchContrib::getGrad(double *pos, double *grad) const
    {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos, "bad vector");
      PRECONDITION(grad, "bad vector");

      double dist = dp_forceField->distance
        (d_at1Idx, d_at2Idx, pos);

      double *at1Coords = &(pos[3 * d_at1Idx]);
      double *at2Coords = &(pos[3 * d_at2Idx]);
      double *g1 = &(grad[3 * d_at1Idx]);
      double *g2 = &(grad[3 * d_at2Idx]);
      double const cs = -2.0;
      double const c1 = 143.9325;
      double const c3 = 7.0 / 12.0;
      double distTerm = dist - d_r0;
      double dE_dr = c1 * d_kb * distTerm
        * (1.0 + 1.5 * cs * distTerm + 2.0 * c3 * cs * cs * distTerm * distTerm);
      double dGrad;
      for (unsigned int i = 0; i < 3; ++i) {
        dGrad = ((dist > 0.0)
          ? (dE_dr * (at1Coords[i] - at2Coords[i]) / dist) : d_kb * 0.01);
        g1[i] += dGrad;
        g2[i] -= dGrad;
      }    
    }
  }
}  
