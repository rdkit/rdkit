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
#ifndef RD_TORSIONANGLECONTRIBS_H
#define RD_TORSIONANGLECONTRIBS_H
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
struct RDKIT_FORCEFIELDHELPERS_EXPORT TorsionAngleContribsParams {
  unsigned int idx1{0};
  unsigned int idx2{0};
  unsigned int idx3{0};
  unsigned int idx4{0};
  std::vector<double> forceConstants{6, 1.0};
  std::vector<int> signs{6, 1};
  TorsionAngleContribsParams(unsigned int idx1, unsigned int idx2,
                             unsigned int idx3, unsigned int idx4,
                             std::vector<double> forceConstants,
                             std::vector<int> signs)
      : idx1(idx1),
        idx2(idx2),
        idx3(idx3),
        idx4(idx4),
        forceConstants(forceConstants),
        signs(signs) {}
};

class RDKIT_FORCEFIELDHELPERS_EXPORT TorsionAngleContribs
    : public ForceFieldContrib {
 public:
  TorsionAngleContribs() = default;

  //! Constructor
  /*!
    \param owner  pointer to the owning ForceField
  */
  TorsionAngleContribs(ForceField *owner);
  ~TorsionAngleContribs() = default;
  //! Add contribution to this collection.
  /*!
    \param idx1           index of atom1 in the ForceField's positions
    \param idx2           index of atom2 in the ForceField's positions
    \param idx3           index of atom3 in the ForceField's positions
    \param idx4           index of atom4 in the ForceField's positions
    \param forceConstants force constants for the torsion potentials
    \param signs          phase for the torsion potentials (-1 or 1)

  */
  void addContrib(unsigned int idx1, unsigned int idx2, unsigned int idx3,
                  unsigned int idx4, std::vector<double> forceConstants,
                  std::vector<int> signs);
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
  TorsionAngleContribs *copy() const override {
    return new TorsionAngleContribs(*this);
  }

  //! Return true if there are no contributions in this contrib
  bool empty() const { return d_contribs.empty(); }

  //! Get number of contributions in this contrib
  unsigned int size() const { return d_contribs.size(); }

 private:
  std::vector<TorsionAngleContribsParams> d_contribs;
};

//! Calculate the torsion energy as described in 10.1021/acs.jcim.5b00654, this
//! can be used with any i > 0.
/*!
 \param forceConstants Force constants for the different cosine fits
 \param signs          Phases of the cosine fits
 \param cosPhi         cosine of the torsion angle phi
*/
RDKIT_FORCEFIELDHELPERS_EXPORT double calcTorsionEnergy(
    const std::vector<double> &forceConstants, const std::vector<int> &signs,
    const double cosPhi);

}  // namespace CrystalFF
}  // namespace ForceFields

#endif