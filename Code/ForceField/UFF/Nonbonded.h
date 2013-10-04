//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_NONBONDED_H__
#define __RD_NONBONDED_H__
#include <ForceField/Contrib.h>

namespace ForceFields {
  namespace UFF {
    class AtomicParams;

    //! the van der Waals term for the Universal Force Field
    /*!
      <b>The Distance Threshold</b>
       For the sake of efficiency, each vdwContrib maintains a threshold
       distance.  When the distance between the two atoms exceeds this
       threshold, the vdwContrib makes no contribution to either the
       energy or the gradient.

       The threshold is set to a multiple of the vdW distance's preferred
       length. This multiplier can be supplied to the constructor.
       
     */
    class vdWContrib : public ForceFieldContrib {
    public:
      vdWContrib() : d_at1Idx(-1), d_at2Idx(-1) {};

      //! Constructor
      /*!
	\param owner       pointer to the owning ForceField
	\param idx1        index of end1 in the ForceField's positions
	\param idx2        index of end2 in the ForceField's positions
	\param at1Params   pointer to the parameters for end1
	\param at2Params   pointer to the parameters for end2
	\param threshMultiplier (optional) multiplier for the threshold
	       calculation. See class documentation for details.

      */
      vdWContrib(ForceField *owner,unsigned int idx1,unsigned int idx2,
		 const AtomicParams *at1Params,
		 const AtomicParams *at2Params,
		 double threshMultiplier=10.0);
      double getEnergy(double *pos) const;
      void getGrad(double *pos,double *grad) const;
    
    private:
      int d_at1Idx,d_at2Idx;
      double d_xij;       //!< the preferred length of the contact
      double d_wellDepth; //!< the vdW well depth (strength of the interaction)
      double d_thresh;    //!< the distance threshold

    };
    namespace Utils {
      //! calculates and returns the UFF minimum position for a vdW contact
      /*!

	\param at1Params  pointer to the parameters for end1
	\param at2Params  pointer to the parameters for end2

	\return the position of the minimum

      */
      double calcNonbondedMinimum(const AtomicParams *at1Params,
				  const AtomicParams *at2Params);

      //! calculates and returns the UFF well depth for a vdW contact
      /*!

	\param at1Params  pointer to the parameters for end1
	\param at2Params  pointer to the parameters for end2

	\return the depth of the well

      */
      double calcNonbondedDepth(const AtomicParams *at1Params,
				const AtomicParams *at2Params);


    }

  }
}
#endif
