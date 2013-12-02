//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_TORSIONANGLE_H__
#define __RD_TORSIONANGLE_H__

#include <ForceField/Contrib.h>
#include <Geometry/point.h>

// we need this so that we get the hybridizations:
#include <GraphMol/Atom.h>

namespace RDGeom {
  class Point3D;
}

namespace ForceFields {
  namespace UFF {
    class AtomicParams;

    //! the torsion term for the Universal Force Field
    class TorsionAngleContrib : public ForceFieldContrib {
    public:
      TorsionAngleContrib() : d_at1Idx(-1), d_at2Idx(-1), d_at3Idx(-1), d_at4Idx(-1), d_order(0) {};
      //! Constructor
      /*!
	The torsion is between atom1 - atom2 - atom3 - atom4
	(i.e the angle between bond atom1-atom2 and bond atom3-atom4
	while looking down bond atom2-atom3)
	
	\param owner       pointer to the owning ForceField
	\param idx1        index of atom1 in the ForceField's positions
	\param idx2        index of atom2 in the ForceField's positions
	\param idx3        index of atom3 in the ForceField's positions
	\param idx4        index of atom4 in the ForceField's positions
	\param bondOrder23 order of the torsional bond between atoms 2 and 3 (as a double)
	\param atNum2      atomic number of atom2
	\param atNum3      atomic number of atom3
	\param hyb2        hybridization of atom2
	\param hyb3        hybridization of atom3
	\param at2Params   pointer to the parameters for atom 2
	\param at3Params   pointer to the parameters for atom 3
	\param endAtomIsSP2 (optional)
	  This boolean is used to signal whether either atom1 or atom4 are
	  RDKit::Atom::SP2 hybridized.
	  This triggers a special case when either of these cases holds:
	   - atom1 is RDKit::Atom::SP2, atom2 is RDKit::Atom::SP2
	     and atom3 is RDKit::Atom::SP3
	   - atom4 is RDKit::Atom::SP2, atom3 is RDKit::Atom::SP2
	     and atom2 is RDKit::Atom::SP3
      */
      TorsionAngleContrib(ForceField *owner,
			  unsigned int idx1,unsigned int idx2,
			  unsigned int idx3,unsigned int idx4,
			  double bondOrder23,
			  int atNum2,int atNum3,
			  RDKit::Atom::HybridizationType hyb2,
			  RDKit::Atom::HybridizationType hyb3,
			  const AtomicParams *at2Params,
			  const AtomicParams *at3Params,
			  bool endAtomIsSP2=false);
      double getEnergy(double *pos) const;
      void getGrad(double *pos,double *grad) const;
      void scaleForceConstant(unsigned int count) { this->d_forceConstant /= static_cast<double>(count); };
    private:
      int d_at1Idx,d_at2Idx,d_at3Idx,d_at4Idx;
      unsigned int d_order;
      double d_forceConstant,d_cosTerm;

      //! returns dE/dTheta
      double getThetaDeriv(double cosTheta,double sinTheta) const;

      //! calculate default values of the torsion parameters.
      /*!
	 see the constructor for an explanation of the arguments
      */
      void calcTorsionParams(double bondOrder23,
			     int atNum2,int atNum3,
			     RDKit::Atom::HybridizationType hyb2,
			     RDKit::Atom::HybridizationType hyb3,
			     const AtomicParams *at2Params,
			     const AtomicParams *at3Params,
			     bool endAtomIsSP2);

    };

    namespace Utils {
      //! calculates and returns the cosine of a torsion angle
      double calculateCosTorsion(const RDGeom::Point3D &p1,const RDGeom::Point3D &p2,
				 const RDGeom::Point3D &p3,const RDGeom::Point3D &p4);
      void calcTorsionGrad(RDGeom::Point3D *r, RDGeom::Point3D *t,
        double *d, double **g, double &sinTerm, double &cosPhi);
    }
  }
}
#endif
