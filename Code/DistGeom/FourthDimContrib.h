//
// Created by Santosh Putta, Nov 2006
//
#include <RDGeneral/export.h>
#ifndef __RD_FOURTHDIMCONTRIB_H__
#define __RD_FOURTHDIMCONTRIB_H__

#include <RDGeneral/Invariant.h>
#include <ForceField/Contrib.h>
#include <ForceField/ForceField.h>

namespace DistGeom {
//! DEPRECATED: use FourthDimContribs instead
//! A term used in penalizing chirality violations
//!
class RDKIT_DISTGEOMETRY_EXPORT FourthDimContrib
    : public ForceFields::ForceFieldContrib {
 public:
  FourthDimContrib() {}

  //! Constructor
  /*!
    \param owner       pointer to the owning ForceField
    \param idx        the index of the atom to be considered
    \param weight     (optional) the weight to be used for this contrib

  */
  FourthDimContrib(ForceFields::ForceField *owner, unsigned int idx,
                   double weight)
      : d_idx(idx), d_weight(weight) {
    PRECONDITION(owner, "bad force field");
    PRECONDITION(owner->dimension() == 4, "force field has wrong dimension");
    dp_forceField = owner;
  }

  //! return the contribution of this contrib to the energy of a given state
  double getEnergy(double *pos) const override {
    PRECONDITION(dp_forceField, "no owner");
    PRECONDITION(dp_forceField->dimension() == 4,
                 "force field has wrong dimension");
    PRECONDITION(pos, "bad vector");
    unsigned int pid = d_idx * dp_forceField->dimension() + 3;
    return d_weight * pos[pid] * pos[pid];
  }

  //! calculate the contribution of this contrib to the gradient at a given
  /// state
  void getGrad(double *pos, double *grad) const override {
    PRECONDITION(dp_forceField, "no owner");
    PRECONDITION(dp_forceField->dimension() == 4,
                 "force field has wrong dimension");
    PRECONDITION(pos, "bad vector");
    unsigned int pid = d_idx * dp_forceField->dimension() + 3;
    grad[pid] += d_weight * pos[pid];
  }
  FourthDimContrib *copy() const override {
    return new FourthDimContrib(*this);
  }

 private:
  unsigned int d_idx{0};
  double d_weight{0.0};
};
}  // namespace DistGeom

#endif
