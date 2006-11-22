//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef __RD_FOURTHDIMCONTRIB_H__
#define __RD_FOURTHDIMCONTRIB_H__

#include <ForceField/Contrib.h>
#include <ForceField/ForceField.h>

namespace DistGeom {
  class FourthDimContrib : public ForceFields::ForceFieldContrib {
  public:
    FourthDimContrib() : d_idx(0), d_weight(0.0) {};

    FourthDimContrib(unsigned int idx, double weight) {
      d_idx = idx;
      d_weight = weight;
    }

    double getEnergy(double *pos) const {
      unsigned int pid = d_idx*dp_forceField->dimension() + 3;
      return d_weight*pos[pid]*pos[pid];
    }
    
    void getGrad(double *pos, double *grad) const {
      unsigned int pid = d_idx*dp_forceField->dimension() + 3;
      grad[pid] += d_weight*pos[pid];
    }
  private:
    unsigned int d_idx;
    double d_weight;
  };
}

#endif
    
