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
#ifndef __RD_MMFFNONBONDED_H__
#define __RD_MMFFNONBONDED_H__
#include <ForceField/Contrib.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>

namespace ForceFields {
  namespace MMFF {
    class MMFFVdWCollection;
    class MMFFVdW;
    //! the van der Waals term for MMFF
    class VdWContrib : public ForceFieldContrib {
    public:
      VdWContrib() : d_at1Idx(-1), d_at2Idx(-1) {};

      //! Constructor
      /*!
	\param owner       pointer to the owning ForceField
	\param idx1        index of end1 in the ForceField's positions
	\param idx2        index of end2 in the ForceField's positions

      */
      VdWContrib(ForceField *owner, unsigned int idx1, unsigned int idx2,
        MMFFVdWCollection *mmffVdW, const MMFFVdW *mmffVdWParamsAtom1,
        const MMFFVdW *mmffVdWParamsAtom2);
      double getEnergy(double *pos) const;
      void getGrad(double *pos, double *grad) const;
    
    private:
      int d_at1Idx, d_at2Idx;
      double d_R_star_ij;       //!< the preferred length of the contact
      double d_wellDepth; //!< the vdW well depth (strength of the interaction)

    };

    //! the electrostatic term for MMFF
    class EleContrib : public ForceFieldContrib {
    public:
      EleContrib() : d_at1Idx(-1), d_at2Idx(-1) {};

      //! Constructor
      /*!
	\param owner       pointer to the owning ForceField
	\param idx1        index of end1 in the ForceField's positions
	\param idx2        index of end2 in the ForceField's positions

      */
      EleContrib(ForceField *owner, unsigned int idx1, unsigned int idx2,
        double chargeTerm, boost::uint8_t dielModel, bool is1_4);
      double getEnergy(double *pos) const;
      void getGrad(double *pos, double *grad) const;
    
    private:
      int d_at1Idx, d_at2Idx;
      double d_chargeTerm;    //!< q1 * q2 / D
      boost::uint8_t d_dielModel;    //!< dielectric model (1: constant; 2: distance-dependent)
      bool d_is1_4;    //!< flag set for atoms in a 1,4 relationship

    };

    namespace Utils {
      //! calculates and returns the unscaled minimum distance (R*ij) for a MMFF VdW contact
      double calcUnscaledVdWMinimum(MMFFVdWCollection *mmffVdW,
        const MMFFVdW *mmffVdWParamsAtom1, const MMFFVdW *mmffVdWParamsAtom2);
      //! calculates and returns the unscaled well depth (epsilon) for a MMFF VdW contact
      double calcUnscaledVdWWellDepth(double R_star_ij,
        const MMFFVdW *mmffVdWParamsIAtom, const MMFFVdW *mmffVdWParamsJAtom);
      //! scales the VdW parameters
      void scaleVdWParams(double &R_star_ij, double &wellDepth,
        MMFFVdWCollection *mmffVdW, const MMFFVdW *mmffVdWParamsIAtom,
        const MMFFVdW *mmffVdWParamsJAtom);
      //! calculates and returns the Van der Waals MMFF energy
      double calcVdWEnergy(const double dist,
        const double R_star_ij, const double wellDepth);
      //! calculates and returns the electrostatic MMFF energy
      double calcEleEnergy(unsigned int idx1, unsigned int idx2, double dist,
        double chargeTerm, boost::uint8_t dielModel, bool is1_4);
    }
  }
}
#endif
