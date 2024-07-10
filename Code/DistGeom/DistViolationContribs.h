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
#ifndef RD_DISTVIOLATIONCONTRIBS_H
#define RD_DISTVIOLATIONCONTRIBS_H

#include <vector>
#include <ForceField/Contrib.h>

namespace DistGeom {

struct DistViolationContribsParams {
  unsigned int idx1{0};  //!< index of end1 in the ForceField's positions
  unsigned int idx2{0};  //!< index of end2 in the ForceField's positions
  double ub{1000.0};     //!< upper bound on the distance
  double lb{0.0};        //!< lower bound on the distance
  double weight{1.0};    //!< used to adjust relative contribution weights
  DistViolationContribsParams(unsigned int i1, unsigned int i2, double u,
                              double l, double w = 1.0)
      : idx1(i1), idx2(i2), ub(u), lb(l), weight(w) {};
};
//! A term to capture all violations of the upper and lower bounds by
//! distance between two points
class RDKIT_DISTGEOMETRY_EXPORT DistViolationContribs
    : public ForceFields::ForceFieldContrib {
 public:
  DistViolationContribs() = default;

  //! Constructor
  /*!
    \param owner       pointer to the owning ForceField
  */
  DistViolationContribs(ForceFields::ForceField *owner);

  double getEnergy(double *pos) const override;

  void getGrad(double *pos, double *grad) const override;
  DistViolationContribs *copy() const override {
    return new DistViolationContribs(*this);
  }

  void addContrib(unsigned int idx1, unsigned int idx2, double ub, double lb,
                  double weight = 1.0) {
    d_contribs.emplace_back(idx1, idx2, ub, lb, weight);
  }
  bool empty() const { return d_contribs.empty(); }
  unsigned int size() const { return d_contribs.size(); }

 private:
  std::vector<DistViolationContribsParams> d_contribs;
};
}  // namespace DistGeom

#endif
