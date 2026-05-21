//
//  Copyright (C) 2026 Niels Maeder and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <RDGeneral/export.h>
#ifndef RD_CRYSTALFF_OOP_H
#define RD_CRYSTALFF_OOP_H
#include <ForceField/Contrib.h>
#include <ForceField/ForceField.h>
#include <vector>

namespace ForceFields {
namespace CrystalFF {

struct RDKIT_FORCEFIELDHELPERS_EXPORT PlanarityContribsParams {
  std::size_t idx1{0};        //!< index of atom1 in the ForceField's positions
  std::size_t idx2{0};        //!< index of atom2 in the ForceField's positions
  std::size_t idx3{0};        //!< index of atom3 in the ForceField's positions
  std::size_t idx4{0};        //!< index of atom4 in the ForceField's positions
  double forceConstant{1.0};  //!< force constant
  PlanarityContribsParams(std::size_t idx1, std::size_t idx2, std::size_t idx3,
                          std::size_t idx4, double forceConstant = 1.0)
      : idx1(idx1),
        idx2(idx2),
        idx3(idx3),
        idx4(idx4),
        forceConstant(forceConstant){};
};
//! A term to capture all Inversion Contributionss.
class RDKIT_FORCEFIELDHELPERS_EXPORT PlanarityContribs
    : public ForceFieldContrib {
 public:
  PlanarityContribs() = default;
  //! Constructor
  /*!
    \param owner  pointer to the owning ForceField
  */
  PlanarityContribs(ForceField *owner);

  ~PlanarityContribs() override = default;
  //! Add contribution to this contrib.
  /*!
    \param idx1        index of atom1 in the ForceField's positions
    \param idx2        index of atom2 in the ForceField's positions
    \param idx3        index of atom3 in the ForceField's positions
    \param idx4        index of atom4 in the ForceField's positions
    \param forceConst  scaling factor for force constant

  */
  void addContrib(std::size_t idx1, std::size_t idx2, std::size_t idx3,
                  std::size_t idx4, double forceConstant = 1.0);
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
  PlanarityContribs *copy() const override {
    return new PlanarityContribs(*this);
  }

  //! Return true if there are no contributions in this contrib
  bool empty() const { return d_contribs.empty(); }

  //! Get number of contributions in this contrib
  unsigned int size() const { return d_contribs.size(); }

 private:
  std::vector<PlanarityContribsParams> d_contribs;
};

}  // namespace CrystalFF
}  // namespace ForceFields

#endif
