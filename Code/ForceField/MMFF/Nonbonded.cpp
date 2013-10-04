// $Id$
//
//  Copyright (C) 2013 Paolo Tosco
//
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Nonbonded.h"
#include "Params.h"
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>

namespace ForceFields {
  namespace MMFF {
    namespace Utils {
      double calcUnscaledVdWMinimum(MMFFVdWCollection *mmffVdW,
        const MMFFVdW *mmffVdWParamsIAtom, const MMFFVdW *mmffVdWParamsJAtom)
      {
        double gamma_ij = (mmffVdWParamsIAtom->R_star - mmffVdWParamsJAtom->R_star)
          / (mmffVdWParamsIAtom->R_star + mmffVdWParamsJAtom->R_star);
        
        return (0.5 * (mmffVdWParamsIAtom->R_star + mmffVdWParamsJAtom->R_star)
          * (1.0 + (((mmffVdWParamsIAtom->DA == 'D') || (mmffVdWParamsJAtom->DA == 'D'))
            ? 0.0 : mmffVdW->B * (1.0 - exp(-(mmffVdW->Beta) * gamma_ij * gamma_ij)))));
      }

      double calcUnscaledVdWWellDepth(double R_star_ij,
        const MMFFVdW *mmffVdWParamsIAtom, const MMFFVdW *mmffVdWParamsJAtom)
      {
        double R_star_ij2 = R_star_ij * R_star_ij;
        
        return (181.16 * mmffVdWParamsIAtom->G_i * mmffVdWParamsJAtom->G_i
          * mmffVdWParamsIAtom->alpha_i * mmffVdWParamsJAtom->alpha_i
          / ((sqrt(mmffVdWParamsIAtom->alpha_i / mmffVdWParamsIAtom->N_i)
          + sqrt(mmffVdWParamsJAtom->alpha_i / mmffVdWParamsJAtom->N_i))
          * R_star_ij2 * R_star_ij2 * R_star_ij2));
      }

      double calcVdWEnergy(const double dist,
        const double R_star_ij, const double wellDepth)
      {
        double dist2 = dist * dist;
        double dist7 = dist2 * dist2 * dist2 * dist;
        double aTerm = 1.07 * R_star_ij / (dist + 0.07 * R_star_ij);
        double aTerm2 = aTerm * aTerm;
        double aTerm7 = aTerm2 * aTerm2 * aTerm2 * aTerm;
        double R_star_ij2 = R_star_ij * R_star_ij;
        double R_star_ij7 = R_star_ij2 * R_star_ij2 * R_star_ij2 * R_star_ij;
        double bTerm = 1.12 * R_star_ij7 / (dist7 + 0.12 * R_star_ij7) - 2.0;
        double res = wellDepth * aTerm7 * bTerm;
        
        return res;
      }

      void scaleVdWParams(double &R_star_ij, double &wellDepth,
        MMFFVdWCollection *mmffVdW, const MMFFVdW *mmffVdWParamsIAtom,
        const MMFFVdW *mmffVdWParamsJAtom)
      {
        if (((mmffVdWParamsIAtom->DA == 'D') && (mmffVdWParamsJAtom->DA == 'A'))
          || ((mmffVdWParamsIAtom->DA == 'A') && (mmffVdWParamsJAtom->DA == 'D'))) {
          R_star_ij *= mmffVdW->DARAD;
          wellDepth *= mmffVdW->DAEPS;
        }
      }
      
      double calcEleEnergy(unsigned int idx1, unsigned int idx2, double dist,
        double chargeTerm, boost::uint8_t dielModel, bool is1_4)
      {
        double corr_dist = dist + 0.05;
        if (dielModel == RDKit::MMFF::DISTANCE) {
          corr_dist *= corr_dist;
        }
        return (332.0716 * chargeTerm / corr_dist * (is1_4 ? 0.75 : 1.0));
      }
    } // end of namespace utils
    
    VdWContrib::VdWContrib(ForceField *owner, unsigned int idx1, unsigned int idx2,
      MMFFVdWCollection *mmffVdW, const MMFFVdW *mmffVdWParamsIAtom,
      const MMFFVdW *mmffVdWParamsJAtom)
    {
      PRECONDITION(owner, "bad owner");
      PRECONDITION(mmffVdW, "bad MMFFVdWCollection");
      PRECONDITION(mmffVdWParamsIAtom, "bad MMFFVdW parameters for atom " + idx1);
      PRECONDITION(mmffVdWParamsJAtom, "bad MMFFVdW parameters for atom " + idx2);
      RANGE_CHECK(0, idx1, owner->positions().size() - 1);
      RANGE_CHECK(0, idx2, owner->positions().size() - 1);

      dp_forceField = owner;
      d_at1Idx = idx1;
      d_at2Idx = idx2;

      this->d_R_star_ij = Utils::calcUnscaledVdWMinimum
          (mmffVdW, mmffVdWParamsIAtom, mmffVdWParamsJAtom);
      this->d_wellDepth = Utils::calcUnscaledVdWWellDepth
          (this->d_R_star_ij, mmffVdWParamsIAtom, mmffVdWParamsJAtom);
      Utils::scaleVdWParams(this->d_R_star_ij, this->d_wellDepth,
        mmffVdW, mmffVdWParamsIAtom, mmffVdWParamsJAtom);
    }

    double VdWContrib::getEnergy(double *pos) const
    {
      PRECONDITION(dp_forceField, "no owner");
      PRECONDITION(pos, "bad vector");

      double dist = this->dp_forceField->distance
        (this->d_at1Idx, this->d_at2Idx, pos);
      return Utils::calcVdWEnergy(dist, this->d_R_star_ij, this->d_wellDepth);
    }


    void VdWContrib::getGrad(double *pos, double *grad) const {
      PRECONDITION(dp_forceField, "no owner");
      PRECONDITION(pos, "bad vector");
      PRECONDITION(grad, "bad vector");

      double dist = this->dp_forceField->distance
        (this->d_at1Idx, this->d_at2Idx, pos);
      double *at1Coords = &(pos[3 * this->d_at1Idx]);
      double *at2Coords = &(pos[3 * this->d_at2Idx]);
      double *g1 = &(grad[3 * this->d_at1Idx]);
      double *g2 = &(grad[3 * this->d_at2Idx]);
      double q = dist / this->d_R_star_ij;
      double q2 = q * q;
      double q6 = q2 * q2 * q2;
      double q7 = q6 * q;
      double t = 1.07 / (q + 0.07);
      double t2 = t * t;
      double t7 = t2 * t2 * t2 * t;
      double dE_dr = this->d_wellDepth / this->d_R_star_ij
        * t7 * (-7.84 * q6 / ((q7 + 0.12) * (q7 + 0.12))
        + ((-7.84 / (q7 + 0.12) + 14.0) / (q + 0.07)));
      for (unsigned int i = 0; i < 3; ++i) {
        double dGrad;
        dGrad = ((dist > 0.0)
          ? (dE_dr * (at1Coords[i] - at2Coords[i]) / dist) : this->d_R_star_ij * 0.01);
        g1[i] += dGrad;
        g2[i] -= dGrad;
      }    
    }
  
    EleContrib::EleContrib(ForceField *owner, unsigned int idx1, unsigned int idx2,
      double chargeTerm, boost::uint8_t dielModel, bool is1_4)
    {
      PRECONDITION(owner, "bad owner");
      RANGE_CHECK(0, idx1, owner->positions().size() - 1);
      RANGE_CHECK(0, idx2, owner->positions().size() - 1);

      dp_forceField = owner;
      this->d_at1Idx = idx1;
      this->d_at2Idx = idx2;
      this->d_chargeTerm = chargeTerm;
      this->d_dielModel = dielModel;
      this->d_is1_4 = is1_4;
    }

    double EleContrib::getEnergy(double *pos) const
    {
      PRECONDITION(dp_forceField, "no owner");
      PRECONDITION(pos, "bad vector");

      return Utils::calcEleEnergy(this->d_at1Idx, this->d_at2Idx,
        this->dp_forceField->distance(this->d_at1Idx, this->d_at2Idx, pos),
        this->d_chargeTerm, this->d_dielModel, this->d_is1_4);
    }

    void EleContrib::getGrad(double *pos, double *grad) const
    {
      PRECONDITION(dp_forceField, "no owner");
      PRECONDITION(pos, "bad vector");
      PRECONDITION(grad, "bad vector");

      double dist = this->dp_forceField->distance
        (this->d_at1Idx, this->d_at2Idx, pos);
      double *at1Coords = &(pos[3 * this->d_at1Idx]);
      double *at2Coords = &(pos[3 * this->d_at2Idx]);
      double *g1 = &(grad[3 * this->d_at1Idx]);
      double *g2 = &(grad[3 * this->d_at2Idx]);
      double corr_dist = dist + 0.05;
      corr_dist *= ((this->d_dielModel == RDKit::MMFF::DISTANCE)
        ? corr_dist * corr_dist : corr_dist);
      double dE_dr = -332.0716 * (double)(this->d_dielModel)
        * this->d_chargeTerm / corr_dist  * (this->d_is1_4 ? 0.75 : 1.0);
      for (unsigned int i = 0; i < 3; ++i) {
        double dGrad;
        dGrad = ((dist > 0.0)
          ? (dE_dr * (at1Coords[i] - at2Coords[i]) / dist) : 0.02);
        g1[i] += dGrad;
        g2[i] -= dGrad;
      }    
    }
  }
}
