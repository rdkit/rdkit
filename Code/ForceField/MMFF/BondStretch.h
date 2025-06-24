//  Copyright (C) 2013-2025 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef __RD_MMFFBONDSTRETCH_H__
#define __RD_MMFFBONDSTRETCH_H__
#include <ForceField/Contrib.h>

#include <vector>

namespace ForceFields {
namespace MMFF {
class MMFFBond;
class MMFFBondStretchEmpiricalRule;

//! The bond-stretch term for MMFF
class RDKIT_FORCEFIELD_EXPORT BondStretchContrib : public ForceFieldContrib {
 public:
  BondStretchContrib() {}
  //! Constructor
  BondStretchContrib(ForceField *owner);
  /*! Adds a bond stretch to the contrib
    \param idx1        index of end1 in the ForceField's positions
    \param idx2        index of end2 in the ForceField's positions
    \param bondType    MMFF94 type of the bond (as an unsigned int)
    \param end1Params  pointer to the parameters for end1
    \param end2Params  pointer to the parameters for end2

  */
  void addTerm(const unsigned int idx1,
               const unsigned int idx2,
               const MMFFBond *mmffBondParams);

  double getEnergy(double *pos) const override;

  void getGrad(double *pos, double *grad) const override;

  BondStretchContrib *copy() const override {
    return new BondStretchContrib(*this);
  }

 private:
  std::vector<int> d_at1Idxs, d_at2Idxs;  //!< indices of end points
  std::vector<double> d_r0;               //!< rest length of the bond
  std::vector<double> d_kb;               //!< force constant of the bond
};

namespace Utils {
//! returns the MMFF rest length for a bond
RDKIT_FORCEFIELD_EXPORT double calcBondRestLength(
    const MMFFBond *mmffBondParams);
//! returns the MMFF force constant for a bond
RDKIT_FORCEFIELD_EXPORT double calcBondForceConstant(
    const MMFFBond *mmffBondParams);
//! calculates and returns the bond stretching MMFF energy
RDKIT_FORCEFIELD_EXPORT double calcBondStretchEnergy(const double r0,
                                                     const double kb,
                                                     const double distance);
}  // namespace Utils
}  // namespace MMFF
}  // namespace ForceFields
#endif
