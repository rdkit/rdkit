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
#ifndef __RD_MMFFBONDSTRETCH_H__
#define __RD_MMFFBONDSTRETCH_H__
#include <ForceField/Contrib.h>

namespace ForceFields {
  namespace MMFF {
    class MMFFBond;
    class MMFFBondStretchEmpiricalRule;

    //! The bond-stretch term for MMFF
    class BondStretchContrib : public ForceFieldContrib {
    public:
      BondStretchContrib() : d_at1Idx(-1), d_at2Idx(-1) {};
      //! Constructor
      /*!
	\param owner       pointer to the owning ForceField
	\param idx1        index of end1 in the ForceField's positions
	\param idx2        index of end2 in the ForceField's positions
	\param bondType    MMFF94 type of the bond (as an unsigned int)
	\param end1Params  pointer to the parameters for end1
	\param end2Params  pointer to the parameters for end2
	
      */
      BondStretchContrib(ForceField *owner,
        const unsigned int idx1, const unsigned int idx2,
        const MMFFBond *mmffBondParams);

      double getEnergy(double *pos) const;

      void getGrad(double *pos,double *grad) const;
    
    private:
      int d_at1Idx, d_at2Idx; //!< indices of end points
      double d_r0;        //!< rest length of the bond
      double d_kb;  //!< force constant of the bond

    };
  
    namespace Utils {
      //! returns the MMFF rest length for a bond 
      double calcBondRestLength(const MMFFBond *mmffBondParams);
      //! returns the MMFF force constant for a bond 
      double calcBondForceConstant(const MMFFBond *mmffBondParams);
      //! calculates and returns the bond stretching MMFF energy
      double calcBondStretchEnergy(const double r0, const double kb, const double distance);
    }  
  }
}
#endif
