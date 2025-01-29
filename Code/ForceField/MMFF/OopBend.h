//  Copyright (C) 2013-2025 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef __RD_MMFFOopBend_H__
#define __RD_MMFFOopBend_H__

#include <ForceField/Contrib.h>
#include <Geometry/point.h>

namespace ForceFields {
namespace MMFF {
class MMFFOop;

//! the out-of-plane term for MMFF
class RDKIT_FORCEFIELD_EXPORT OopBendContrib : public ForceFieldContrib {
 public:
  OopBendContrib() {}
  //! Constructor
  OopBendContrib(ForceField *owner);
  /*! Adds an out-of-plane term to the force field contrib.

  The Wilson angle is between the vector formed by atom2-atom4
  and the angle formed by atom1-atom2-atom3

    \param idx1        index of atom1 in the ForceField's positions
    \param idx2        index of atom2 in the ForceField's positions
    \param idx3        index of atom3 in the ForceField's positions
    \param idx4        index of atom4 in the ForceField's positions
  */
  void addTerm(unsigned int idx1, unsigned int idx2,
               unsigned int idx3, unsigned int idx4,
               const MMFFOop *mmffOopParams);

  double getEnergy(double *pos) const override;
  void getGrad(double *pos, double *grad) const override;
  OopBendContrib *copy() const override { return new OopBendContrib(*this); }

  void getSingleGrad(double* pos, double* grad, unsigned int termIdx) const;

 private:
  std::vector<int> d_at1Idxs, d_at2Idxs, d_at3Idxs, d_at4Idxs;
  std::vector<double> d_koop;
};

namespace Utils {
//! calculates and returns the Wilson angle (in degrees)
RDKIT_FORCEFIELD_EXPORT double calcOopChi(const RDGeom::Point3D &iPoint,
                                          const RDGeom::Point3D &jPoint,
                                          const RDGeom::Point3D &kPoint,
                                          const RDGeom::Point3D &lPoint);
//! returns the out-of-plane force constant koop
RDKIT_FORCEFIELD_EXPORT double calcOopBendForceConstant(
    const MMFFOop *mmffOopParams);
//! calculates and returns the out-of-plane MMFF energy
RDKIT_FORCEFIELD_EXPORT double calcOopBendEnergy(const double chi,
                                                 const double koop);
}  // namespace Utils
}  // namespace MMFF
}  // namespace ForceFields
#endif
