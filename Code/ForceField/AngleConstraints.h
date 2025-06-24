//
//  Copyright (C) 2024 Niels Maeder and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RD_ANGLECONSTRAINTS_H
#define RD_ANGLECONSTRAINTS_H
#include <vector>
#include <iostream>
#include "ForceField.h"
#include "Contrib.h"

namespace ForceFields {

struct AngleConstraintContribsParams {
  unsigned int idx1{0};       //!< index of atom1 of the angle constraint
  unsigned int idx2{0};       //!< index of atom2 of the angle constraint
  unsigned int idx3{0};       //!< index of atom3 of the angle constraint
  double minAngle{0.0};       //!< lower bound of the flat bottom potential
  double maxAngle{0.0};       //!< upper bound of the flat bottom potential
  double forceConstant{1.0};  //!< force constant for angle constraint
  AngleConstraintContribsParams(unsigned int idx1, unsigned int idx2,
                                unsigned int idx3, double minAngle,
                                double maxAngle, double forceConstant = 1.0)
      : idx1(idx1),
        idx2(idx2),
        idx3(idx3),
        minAngle(minAngle),
        maxAngle(maxAngle),
        forceConstant(forceConstant) {};
};

//! A term to capture all flat bottom angle constraint potentials.
class RDKIT_FORCEFIELD_EXPORT AngleConstraintContribs
    : public ForceFieldContrib {
 public:
  AngleConstraintContribs() = default;

  //! Constructor
  /*!
    \param owner  pointer to the owning ForceField
  */
  AngleConstraintContribs(ForceField *owner);
  ~AngleConstraintContribs() override = default;
  //! Add a contribution to this contrib collection
  /*!
  \param idx1        index of atom1 in the ForceField's positions
  \param idx2        index of atom2 in the ForceField's positions
  \param idx3        index of atom3 in the ForceField's positions
  \param minAngle    minimum angle
  \param maxAngle    maximum angle
  \param forceConst  force Constant

  */
  void addContrib(unsigned int idx1, unsigned int idx2, unsigned int idx3,
                  double minAngleDeg, double maxAngleDeg, double forceConst);
  //! Add a contribution to this contrib collection
  /*!
  \param idx1        index of atom1 in the ForceField's positions
  \param idx2        index of atom2 in the ForceField's positions
  \param idx3        index of atom3 in the ForceField's positions
  \param relative    whether to add the provided angle to the current angle
  \param minAngle    minimum angle
  \param maxAngle    maximum angle
  \param forceConst  force Constant

  */
  void addContrib(unsigned int idx1, unsigned int idx2, unsigned int idx3,
                  bool relative, double minAngleDeg, double maxAngleDeg,
                  double forceConst);
  //! return the contribution of this contrib to the energy of a given state
  /*!
    \param pos  positions of the atoms in the current state
  */
  double getEnergy(double *pos) const override;
  //! calculate the contribution of this contrib to the gradient at a given
  /// state
  /*!
    \param pos  positions of the atoms in the current state
    \param grad gradients to be adapted
  */
  void getGrad(double *pos, double *grad) const override;
  //! Copy constructor

  AngleConstraintContribs *copy() const override {
    return new AngleConstraintContribs(*this);
  }

  //! Return true if there are no contributions in this contrib
  bool empty() const { return d_contribs.empty(); }

  //! Get number of contributions in this contrib
  unsigned int size() const { return d_contribs.size(); }

 private:
  std::vector<AngleConstraintContribsParams> d_contribs;
  double computeAngleTerm(const double &angle,
                          const AngleConstraintContribsParams &contrib) const;
};
}  // namespace ForceFields
#endif
