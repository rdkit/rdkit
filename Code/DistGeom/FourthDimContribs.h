//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_FOURTHDIMCONTRIBS_H
#define RD_FOURTHDIMCONTRIBS_H

#include <vector>
#include <RDGeneral/Invariant.h>
#include <ForceField/Contrib.h>
#include <ForceField/ForceField.h>

namespace DistGeom {

struct FourthDimContribsParams {
  unsigned int idx{0};
  double weight{0.0};
  FourthDimContribsParams(unsigned int idx, double w) : idx(idx), weight(w) {};
};

//! A term used in penalizing the 4th dimension in order to move from 4D->3D
//!
class RDKIT_DISTGEOMETRY_EXPORT FourthDimContribs
    : public ForceFields::ForceFieldContrib {
 public:
  FourthDimContribs() = default;

  //! Constructor
  /*!
    \param owner       pointer to the owning ForceField
    \param idx        the index of the atom to be considered
    \param weight     (optional) the weight to be used for this contrib

  */
  FourthDimContribs(ForceFields::ForceField *owner) {
    PRECONDITION(owner, "bad force field");
    PRECONDITION(owner->dimension() == 4, "force field has wrong dimension");
    dp_forceField = owner;
  }

  void addContrib(unsigned int idx, double weight) {
    d_contribs.emplace_back(idx, weight);
  }

  //! return the contribution of this contrib to the energy of a given state
  double getEnergy(double *pos) const override {
    PRECONDITION(pos, "bad vector");
    constexpr unsigned int ffdim = 4;
    double res = 0.0;
    for (const auto &contrib : d_contribs) {
      unsigned int pid = contrib.idx * ffdim + 3;
      res += contrib.weight * pos[pid] * pos[pid];
    }
    return res;
  }

  //! calculate the contribution of this contrib to the gradient at a given
  /// state
  void getGrad(double *pos, double *grad) const override {
    PRECONDITION(pos, "bad vector");
    constexpr unsigned int ffdim = 4;
    for (const auto &contrib : d_contribs) {
      unsigned int pid = contrib.idx * ffdim + 3;
      grad[pid] += contrib.weight * pos[pid];
    }
  }
  FourthDimContribs *copy() const override {
    return new FourthDimContribs(*this);
  }
  bool empty() const { return d_contribs.empty(); }
  unsigned int size() const { return d_contribs.size(); }

 private:
  std::vector<FourthDimContribsParams> d_contribs;
};
}  // namespace DistGeom

#endif
