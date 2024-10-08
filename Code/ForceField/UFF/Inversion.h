//
//  Copyright (C) 2013 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_UFFINVERSION_H
#define RD_UFFINVERSION_H
#include <ForceField/Contrib.h>
#include <tuple>
#include <Geometry/point.h>

namespace ForceFields {
namespace UFF {
class AtomicParams;

//! The inversion term for the Universal Force Field
class RDKIT_FORCEFIELD_EXPORT InversionContrib : public ForceFieldContrib {
 public:
  InversionContrib() {}
  //! Constructor
  /*!
    \param owner          pointer to the owning ForceField
    \param idx1           index of atom1 in the ForceField's positions
    \param idx2           index of atom2 in the ForceField's positions
    \param idx3           index of atom3 in the ForceField's positions
    \param idx4           index of atom4 in the ForceField's positions
    \param at2AtomicNum   atomic number for atom 2
    \param isCBoundToO    boolean flag; true if atom 2 is sp2 carbon bound to
    sp2 oxygen

  */
  InversionContrib(ForceField *owner, unsigned int idx1, unsigned int idx2,
                   unsigned int idx3, unsigned int idx4, int at2AtomicNum,
                   bool isCBoundToO, double oobForceScalingFactor = 1.0);

  double getEnergy(double *pos) const override;

  void getGrad(double *pos, double *grad) const override;
  InversionContrib *copy() const override {
    return new InversionContrib(*this);
  }

 private:
  int d_at1Idx{-1};
  int d_at2Idx{-1};
  int d_at3Idx{-1};
  int d_at4Idx{-1};
  double d_forceConstant, d_C0, d_C1, d_C2;
};
}  // namespace UFF
}  // namespace ForceFields
#endif
