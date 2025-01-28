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
#ifndef RD_UFFINVERSIONS_H
#define RD_UFFINVERSIONS_H
#include <ForceField/Contrib.h>
#include <vector>

namespace ForceFields {
namespace UFF {
class AtomicParams;

struct RDKIT_FORCEFIELD_EXPORT InversionContribsParams {
  unsigned int idx1{0};       //!< index of atom1 in the ForceField's positions
  unsigned int idx2{0};       //!< index of atom2 in the ForceField's positions
  unsigned int idx3{0};       //!< index of atom3 in the ForceField's positions
  unsigned int idx4{0};       //!< index of atom4 in the ForceField's positions
  int at2AtomicNum{0};        //!< atomic number for atom 2
  bool isCBoundToO{false};    //!< boolean flag; true if atom 2 is sp2 carbon
                              //!< bound to sp2 oxygen
  double C0{0.0};             //!< inversion coefficient 0
  double C1{0.0};             //!< inversion coefficient 1
  double C2{0.0};             //!< inversion coefficient 2
  double forceConstant{1.0};  //!< force constant
  InversionContribsParams(unsigned int idx1, unsigned int idx2,
                          unsigned int idx3, unsigned int idx4,
                          int at2AtomicNum, bool isCBoundToO, double C0,
                          double C1, double C2, double forceConstant = 1.0)
      : idx1(idx1),
        idx2(idx2),
        idx3(idx3),
        idx4(idx4),
        at2AtomicNum(at2AtomicNum),
        isCBoundToO(isCBoundToO),
        C0(C0),
        C1(C1),
        C2(C2),
        forceConstant(forceConstant) {};
};
//! A term to capture all Inversion Contributionss.
class RDKIT_FORCEFIELD_EXPORT InversionContribs : public ForceFieldContrib {
 public:
  InversionContribs() = default;
  //! Constructor
  /*!
    \param owner  pointer to the owning ForceField
  */
  InversionContribs(ForceField *owner);

  ~InversionContribs() override = default;
  //! Add contribution to this contrib.
  /*!
    \param idx1        index of atom1 in the ForceField's positions
    \param idx2        index of atom2 in the ForceField's positions
    \param idx3        index of atom3 in the ForceField's positions
    \param idx4        index of atom4 in the ForceField's positions
    \param at2AtomicNum     atomic number for atom 2
    \param isCBoundToO      boolean flag; true if atom 2 is sp2 C bound to
    sp2 O
    \param oobForceScalingFactor  scaling factor for force constant

  */
  void addContrib(unsigned int idx1, unsigned int idx2, unsigned int idx3,
                  unsigned int idx4, int at2AtomicNum, bool isCBoundToO,
                  double oobForceScalingFactor = 1.0);
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
  InversionContribs *copy() const override {
    return new InversionContribs(*this);
  }

  //! Return true if there are no contributions in this contrib
  bool empty() const { return d_contribs.empty(); }

  //! Get number of contributions in this contrib
  unsigned int size() const { return d_contribs.size(); }

 private:
  std::vector<InversionContribsParams> d_contribs;
};

}  // namespace UFF
}  // namespace ForceFields

#endif