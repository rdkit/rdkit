//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef __RD_CHIRALVIOLATIONCONTRIB_H__
#define __RD_CHIRALVIOLATIONCONTRIB_H__

#include <ForceField/Contrib.h>
#include "EmbedObject.h"
namespace DistGeom {

  
  //! A term to capture the violation of the chirality at an atom center
  //!
  class ChiralViolationContrib : public ForceFields::ForceFieldContrib {
  public:
    ChiralViolationContrib() : d_idx1(0), d_idx2(0), d_idx3(0), d_idx4(0), 
      d_volLower(0.0), d_volUpper(0.0), d_weight(0.0){};
    
    //! Constructor
    /*!
      \param owner      pointer to the owning forcefield
      \param cset       a chiral set containing the four chiral atom ids (insequence)
                        and the upper and lower limits on the signed chiral volume
    */
    ChiralViolationContrib(ForceFields::ForceField *owner, const ChiralSet *cset, double weight=1.0);
    
    double getEnergy(double *pos) const;
    
    void getGrad(double *pos, double *grad) const;

  private:
    unsigned int d_idx1, d_idx2, d_idx3, d_idx4;
    double d_volLower;
    double d_volUpper;
    double d_weight;
  };
}

#endif
