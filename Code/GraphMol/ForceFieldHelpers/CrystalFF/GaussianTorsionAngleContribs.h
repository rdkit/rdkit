//
//  Copyright (C) 2026 Niels Maeder and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RD_GAUSSIANTORSIONANGLECONTRIBS_H
#define RD_GAUSSIANTORSIONANGLECONTRIBS_H
#include <ForceField/Contrib.h>
#include <vector>

namespace RDGeom {
class Point3D;
}

namespace ForceFields {
class ForceField;
class ForceFieldContrib;
}  // namespace ForceFields

namespace ForceFields {
namespace CrystalFF {

//! A term to capture all torsion constraint potentials.
//!
struct RDKIT_FORCEFIELDHELPERS_EXPORT GaussianTorsionAngleContribsParams {
  std::size_t idx1{0};
  std::size_t idx2{0};
  std::size_t idx3{0};
  std::size_t idx4{0};
  std::vector<double> energies{};
  std::vector<double> gradients{};
  double scaling{1.0};
  GaussianTorsionAngleContribsParams(std::size_t idx1, std::size_t idx2,
                                     std::size_t idx3, std::size_t idx4,
                                     std::vector<double> energies,
                                     std::vector<double> gradients,
                                     double scaling = 1.0)
      : idx1(idx1),
        idx2(idx2),
        idx3(idx3),
        idx4(idx4),
        energies(energies),
        gradients(gradients),
        scaling(scaling) {}
};

class RDKIT_FORCEFIELDHELPERS_EXPORT GaussianTorsionAngleContribs
    : public ForceFieldContrib {
 public:
  GaussianTorsionAngleContribs() = default;

  //! Constructor
  /*!
    \param owner  pointer to the owning ForceField
  */
  GaussianTorsionAngleContribs(ForceField *owner);
  ~GaussianTorsionAngleContribs() = default;
  //! Add contribution to this collection.
  /*!
    \param idx1           index of atom1 in the ForceField's positions
    \param idx2           index of atom2 in the ForceField's positions
    \param idx3           index of atom3 in the ForceField's positions
    \param idx4           index of atom4 in the ForceField's positions
    \param heights        heights of the gaussians
    \param positions      positions of the gaussians
    \param widths         widths of the gaussians
    \param scaling        Scaling factor for energy and gradient (for K terms)
  */
  void addContrib(std::size_t idx1, std::size_t idx2, std::size_t idx3,
                  std::size_t idx4, std::vector<double> energies,
                  std::vector<double> gradients, double scaling = 1.0);
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
  GaussianTorsionAngleContribs *copy() const override {
    return new GaussianTorsionAngleContribs(*this);
  }

  //! Return true if there are no contributions in this contrib
  bool empty() const { return d_contribs.empty(); }

  //! Get number of contributions in this contrib
  unsigned int size() const { return d_contribs.size(); }

 private:
  std::vector<GaussianTorsionAngleContribsParams> d_contribs;
};

//! Calculate the torsion energy as described in 10.1021/acs.jcim.5b00654, this
//! can be used with any i > 0.
/*!
 \param heights        heights of the gaussians
 \param positions      positions of the gaussians
 \param widths         widths of the gaussians
 \param cosPhi         cosine of the torsion angle phi
*/
RDKIT_FORCEFIELDHELPERS_EXPORT double getEnergy(
    const std::vector<double> &heights, const std::vector<double> &positions,
    const std::vector<double> &widths, const double cosPhi);

RDKIT_FORCEFIELDHELPERS_EXPORT double getdEdPhi(
    const std::vector<double> &heights, const std::vector<double> &positions,
    const std::vector<double> &widths, const double cosPhi);
}  // namespace CrystalFF
}  // namespace ForceFields

#endif
