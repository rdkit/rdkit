//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_BONDSTRETCH_H__
#define __RD_BONDSTRETCH_H__
#include <ForceField/Contrib.h>

namespace ForceFields {
  namespace UFF {
    class AtomicParams;

    //! The bond-stretch term for the Universal Force Field
    class BondStretchContrib : public ForceFieldContrib {
    public:
      BondStretchContrib() : d_end1Idx(-1), d_end2Idx(-1) {};
      //! Constructor
      /*!
	\param owner       pointer to the owning ForceField
	\param idx1        index of end1 in the ForceField's positions
	\param idx2        index of end2 in the ForceField's positions
	\param bondOrder   order of the bond (as a double)
	\param end1Params  pointer to the parameters for end1
	\param end2Params  pointer to the parameters for end2
	
      */
      BondStretchContrib(ForceField *owner,unsigned int idx1,unsigned int idx2,
			 double bondOrder,
			 const AtomicParams *end1Params,
			 const AtomicParams *end2Params);

      double getEnergy(double *pos) const;

      void getGrad(double *pos,double *grad) const;
    
    private:
      int d_end1Idx,d_end2Idx; //!< indices of end points
      double d_restLen;        //!< rest length of the bond
      double d_forceConstant;  //!< force constant of the bond

    };
  
    namespace Utils {
      //! calculates and returns the UFF rest length for a bond 
      /*!

	\param bondOrder the order of the bond (as a double)
	\param end1Params  pointer to the parameters for end1
	\param end2Params  pointer to the parameters for end2

	\return the rest length

      */
      double calcBondRestLength(double bondOrder,
				const AtomicParams *end1Params,
				const AtomicParams *end2Params);

      //! calculates and returns the UFF force constant for a bond 
      /*!

	\param restLength  the rest length of the bond
	\param end1Params  pointer to the parameters for end1
	\param end2Params  pointer to the parameters for end2

	\return the force constant
	
      */
      double calcBondForceConstant(double restLength,
				   const AtomicParams *end1Params,
				   const AtomicParams *end2Params);

    }  
  }
}
#endif
