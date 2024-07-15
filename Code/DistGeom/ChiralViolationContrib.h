//
// Created by Santosh Putta, Nov 2006
//
#include <RDGeneral/export.h>
#ifndef __RD_CHIRALVIOLATIONCONTRIB_H__
#define __RD_CHIRALVIOLATIONCONTRIB_H__

#include <ForceField/Contrib.h>
#include <Geometry/point.h>

namespace DistGeom {
class ChiralSet;

//! DEPRECATED: use ChiralViolationContribs instead
//! A term to capture the violation of chirality at an atom center
//!
class RDKIT_DISTGEOMETRY_EXPORT ChiralViolationContrib
    : public ForceFields::ForceFieldContrib {
 public:
  ChiralViolationContrib() {}

  //! Constructor
  /*!
    \param owner      pointer to the owning forcefield
    \param cset       a chiral set containing the four chiral atom ids (in
    sequence)
                      and the upper and lower limits on the signed chiral volume
    \param weight     (optional) the weight to be used for this contrib

  */
  ChiralViolationContrib(ForceFields::ForceField *owner, const ChiralSet *cset,
                         double weight = 1.0);

  //! return the contribution of this contrib to the energy of a given state
  double getEnergy(double *pos) const override;

  //! calculate the contribution of this contrib to the gradient at a given
  /// state
  void getGrad(double *pos, double *grad) const override;
  ChiralViolationContrib *copy() const override {
    return new ChiralViolationContrib(*this);
  }

  static double calcChiralVolume(unsigned int idx1, unsigned int idx2,
                                 unsigned int idx3, unsigned int idx4,
                                 const double *pos, unsigned int dim) {
    // even if we are minimizing in higher dimension the chiral volume is
    // calculated using only the first 3 dimensions
    RDGeom::Point3D v1(pos[idx1 * dim] - pos[idx4 * dim],
                       pos[idx1 * dim + 1] - pos[idx4 * dim + 1],
                       pos[idx1 * dim + 2] - pos[idx4 * dim + 2]);

    RDGeom::Point3D v2(pos[idx2 * dim] - pos[idx4 * dim],
                       pos[idx2 * dim + 1] - pos[idx4 * dim + 1],
                       pos[idx2 * dim + 2] - pos[idx4 * dim + 2]);

    RDGeom::Point3D v3(pos[idx3 * dim] - pos[idx4 * dim],
                       pos[idx3 * dim + 1] - pos[idx4 * dim + 1],
                       pos[idx3 * dim + 2] - pos[idx4 * dim + 2]);

    RDGeom::Point3D v2xv3 = v2.crossProduct(v3);

    double vol = v1.dotProduct(v2xv3);
    return vol;
  }
  static double calcChiralVolume(unsigned int idx1, unsigned int idx2,
                                 unsigned int idx3, unsigned int idx4,
                                 const RDGeom::PointPtrVect &pts) {
    // even if we are minimizing in higher dimension the chiral volume is
    // calculated using only the first 3 dimensions
    RDGeom::Point3D v1((*pts[idx1])[0] - (*pts[idx4])[0],
                       (*pts[idx1])[1] - (*pts[idx4])[1],
                       (*pts[idx1])[2] - (*pts[idx4])[2]);

    RDGeom::Point3D v2((*pts[idx2])[0] - (*pts[idx4])[0],
                       (*pts[idx2])[1] - (*pts[idx4])[1],
                       (*pts[idx2])[2] - (*pts[idx4])[2]);

    RDGeom::Point3D v3((*pts[idx3])[0] - (*pts[idx4])[0],
                       (*pts[idx3])[1] - (*pts[idx4])[1],
                       (*pts[idx3])[2] - (*pts[idx4])[2]);

    RDGeom::Point3D v2xv3 = v2.crossProduct(v3);

    double vol = v1.dotProduct(v2xv3);
    return vol;
  }

 private:
  unsigned int d_idx1{0}, d_idx2{0}, d_idx3{0}, d_idx4{0};
  double d_volLower{0.0};
  double d_volUpper{0.0};
  double d_weight{0.0};
};
}  // namespace DistGeom

#endif
