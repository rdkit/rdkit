//  Copyright (C) 2013-2025 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef __RD_MMFFNONBONDED_H__
#define __RD_MMFFNONBONDED_H__
#include <ForceField/Contrib.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>

namespace ForceFields {
namespace MMFF {
class MMFFVdWCollection;
class MMFFVdW;
//! the van der Waals term for MMFF
class RDKIT_FORCEFIELD_EXPORT VdWContrib : public ForceFieldContrib {
 public:
  VdWContrib() {}
  VdWContrib(ForceField *owner);
  //! Track a new VdW pair
  void addTerm(unsigned int idx1, unsigned int idx2, const MMFFVdWRijstarEps *mmffVdWConstants);
  double getEnergy(double *pos) const override;
  void getGrad(double *pos, double *grad) const override;
  VdWContrib *copy() const override { return new VdWContrib(*this); }

 private:
  std::vector<int16_t> d_at1Idxs;
  std::vector<int16_t> d_at2Idxs;
  std::vector<double> d_R_ij_stars; //!< the preferred length of the contact
  std::vector<double> d_wellDepths; //!< the vdW well depth (strength of the interaction)
};

//! the electrostatic term for MMFF
class RDKIT_FORCEFIELD_EXPORT EleContrib : public ForceFieldContrib {
 public:
  EleContrib() {}

  //! Constructor
  /*!
    \param owner       pointer to the owning ForceField
    \param idx1        index of end1 in the ForceField's positions
    \param idx2        index of end2 in the ForceField's positions

  */
  EleContrib(ForceField *owner);
  void addTerm(unsigned int idx1, unsigned int idx2,
               double chargeTerm, std::uint8_t dielModel, bool is1_4);
  double getEnergy(double *pos) const override;
  void getGrad(double *pos, double *grad) const override;

  EleContrib *copy() const override { return new EleContrib(*this); }

 private:
  std::vector<int16_t> d_at1Idxs;
  std::vector<int16_t> d_at2Idxs;
  std::vector<double> d_chargeTerms;
  std::vector<std::uint8_t> d_is_1_4s;
  std::vector<std::uint8_t>
      d_dielModels;  //!< dielectric model (1: constant; 2: distance-dependent)
};

namespace Utils {
//! calculates and returns the unscaled minimum distance (R*ij) for a MMFF VdW
/// contact
RDKIT_FORCEFIELD_EXPORT double calcUnscaledVdWMinimum(
    const MMFFVdWCollection *mmffVdW, const MMFFVdW *mmffVdWParamsAtom1,
    const MMFFVdW *mmffVdWParamsAtom2);
//! calculates and returns the unscaled well depth (epsilon) for a MMFF VdW
/// contact
RDKIT_FORCEFIELD_EXPORT double calcUnscaledVdWWellDepth(
    double R_star_ij, const MMFFVdW *mmffVdWParamsIAtom,
    const MMFFVdW *mmffVdWParamsJAtom);
//! scales the VdW parameters
RDKIT_FORCEFIELD_EXPORT void scaleVdWParams(double &R_star_ij,
                                            double &wellDepth,
                                            const MMFFVdWCollection *mmffVdW,
                                            const MMFFVdW *mmffVdWParamsIAtom,
                                            const MMFFVdW *mmffVdWParamsJAtom);
//! calculates and returns the Van der Waals MMFF energy
RDKIT_FORCEFIELD_EXPORT double calcVdWEnergy(const double dist,
                                             const double R_star_ij,
                                             const double wellDepth);
//! calculates and returns the electrostatic MMFF energy
// FIX: idx1 and idx2 are not used
RDKIT_FORCEFIELD_EXPORT double calcEleEnergy(unsigned int idx1,
                                             unsigned int idx2, double dist,
                                             double chargeTerm,
                                             std::uint8_t dielModel,
                                             bool is1_4);
}  // namespace Utils
}  // namespace MMFF
}  // namespace ForceFields
#endif
