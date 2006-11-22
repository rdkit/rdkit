//
// Created by Santosh Putta, Nov 2006
//
#ifndef __RD_CHIRALVIOLATIONCONTRIB_H__
#define __RD_CHIRALVIOLATIONCONTRIB_H__

#include <ForceField/Contrib.h>

namespace DistGeom {
  class ChiralSet;
  
  //! A term to capture the violation of chirality at an atom center
  //!
  class ChiralViolationContrib : public ForceFields::ForceFieldContrib {
  public:
    ChiralViolationContrib() : d_idx1(0), d_idx2(0), d_idx3(0), d_idx4(0), 
      d_volLower(0.0), d_volUpper(0.0), d_weight(0.0){};
    
    //! Constructor
    /*!
      \param owner      pointer to the owning forcefield
      \param cset       a chiral set containing the four chiral atom ids (in sequence)
                        and the upper and lower limits on the signed chiral volume
      \param weight     (optional) the weight to be used for this contrib
                        
    */
    ChiralViolationContrib(ForceFields::ForceField *owner, const ChiralSet *cset, double weight=1.0);
    
    //! return the contribution of this contrib to the energy of a given state 
    double getEnergy(double *pos) const;
    
    //! calculate the contribution of this contrib to the gradient at a given state
    void getGrad(double *pos, double *grad) const;

  private:
    unsigned int d_idx1, d_idx2, d_idx3, d_idx4;
    double d_volLower;
    double d_volUpper;
    double d_weight;
  };
}

#endif
