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
#ifndef RD_DISTANCECONSTRAINTS_H
#define RD_DISTANCECONSTRAINTS_H
#include <vector>
#include <iostream>
#include "Contrib.h"

namespace ForceFields {

struct DistanceConstraintContribsParams {
  unsigned int idx1{0};       //!< index of atom1 of the distance constraint
  unsigned int idx2{0};       //!< index of atom2 of the distance constraint
  double minLen{0.0};         //!< lower bound of the flat bottom potential
  double maxLen{0.0};         //!< upper bound of the flat bottom potential
  double forceConstant{1.0};  //!< force constant for distance constraint
  DistanceConstraintContribsParams(unsigned int idx1, unsigned int idx2,
                                   double minLen, double maxLen,
                                   double forceConstant = 1.0)
      : idx1(idx1),
        idx2(idx2),
        minLen(minLen),
        maxLen(maxLen),
        forceConstant(forceConstant) {};
};

//! A term to capture all flat bottom distance constraint potentials.
class RDKIT_FORCEFIELD_EXPORT DistanceConstraintContribs
    : public ForceFieldContrib {
 public:
  DistanceConstraintContribs() = default;

  //! Constructor
  /*!
    \param owner  pointer to the owning ForceField
  */
  DistanceConstraintContribs(ForceField *owner);
  ~DistanceConstraintContribs() override = default;
  //! Add contribution to this contrib.
  /*!
    \param idx1        index of atom1 in the ForceField's positions
    \param idx2        index of atom2 in the ForceField's positions
    \param minLen      minimum distance
    \param maxLen      maximum distance
    \param forceConst  force Constant

  */
  void addContrib(unsigned int idx1, unsigned int idx2, double minLen,
                  double maxLen, double forceConstant);
  //! Add contribution to this contrib.
  /*!
    \param idx1        index of atom1 in the ForceField's positions
    \param idx2        index of atom2 in the ForceField's positions
    \param relative    whether to add the provided distance to the
    current distance
    \param minLen      minimum distance
    \param maxLen      maximum distance
    \param forceConst  force Constant

  */
  void addContrib(unsigned int idx1, unsigned int idx2, bool relative,
                  double minLen, double maxLen, double forceConstant);

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
  DistanceConstraintContribs *copy() const override {
    return new DistanceConstraintContribs(*this);
  }

  //! Return true if there are no contributions in this contrib
  bool empty() const { return d_contribs.empty(); }

  //! Get number of contributions in this contrib
  unsigned int size() const { return d_contribs.size(); }

 private:
  std::vector<DistanceConstraintContribsParams> d_contribs;
};

}  // namespace ForceFields

#endif