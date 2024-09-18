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
#ifndef RD_CHIRALVIOLATIONCONTRIBS_H
#define RD_CHIRALVIOLATIONCONTRIBS_H

#include <vector>
#include <ForceField/Contrib.h>
#include <Geometry/point.h>

namespace DistGeom {
class ChiralSet;

struct ChiralViolationContribsParams {
  unsigned int idx1{0}, idx2{0}, idx3{0}, idx4{0};
  double volUpper{0.0};
  double volLower{0.0};
  double weight{1.0};
  ChiralViolationContribsParams(unsigned int i1, unsigned int i2,
                                unsigned int i3, unsigned int i4, double u,
                                double l, double w = 1.0)
      : idx1(i1),
        idx2(i2),
        idx3(i3),
        idx4(i4),
        volUpper(u),
        volLower(l),
        weight(w) {};
};
//! A term to capture the violation of chirality at atom centers
//!
class RDKIT_DISTGEOMETRY_EXPORT ChiralViolationContribs
    : public ForceFields::ForceFieldContrib {
 public:
  ChiralViolationContribs() = default;

  //! Constructor
  /*!
    \param owner      pointer to the owning forcefield
    \param cset       a chiral set containing the four chiral atom ids (in
    sequence)
                      and the upper and lower limits on the signed chiral
    volume \param weight     (optional) the weight to be used for this contrib

  */
  ChiralViolationContribs(ForceFields::ForceField *owner);

  //! adds a new chiral constraint
  /*!
    \param cset       a chiral set containing the four chiral atom ids (in
    sequence)
                      and the upper and lower limits on the signed chiral
    volume \param weight     (optional) the weight to be used for this contrib

  */
  void addContrib(const ChiralSet *cset, double weight = 1.0);

  //! return the contribution of this contrib to the energy of a given state
  double getEnergy(double *pos) const override;

  //! calculate the contribution of this contrib to the gradient at a given
  /// state
  void getGrad(double *pos, double *grad) const override;
  ChiralViolationContribs *copy() const override {
    return new ChiralViolationContribs(*this);
  }
  bool empty() const { return d_contribs.empty(); }
  unsigned int size() const { return d_contribs.size(); }

 private:
  std::vector<ChiralViolationContribsParams> d_contribs;
};
}  // namespace DistGeom

#endif
