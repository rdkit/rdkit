//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDKIT_MOLENUMERATOR_H
#define RDKIT_MOLENUMERATOR_H

#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolBundle.h>

#include <vector>
#include <map>
#include <string>
#include <memory>

namespace RDKit {
class ChemicalReaction;
namespace MolEnumerator {

//! abstract base class for the a molecule enumeration operation
class RDKIT_MOLENUMERATOR_EXPORT MolEnumeratorOp {
 public:
  MolEnumeratorOp(){};
  virtual ~MolEnumeratorOp(){};
  //! returns a vector of the number of possible variations at variability point
  //! covered by this operation
  virtual std::vector<size_t> getVariationCounts() const = 0;
  //! returns a the molecule corresponding to a particular variation
  /*!  which.size() should be equal to the number of variation counts.
   */
  virtual std::unique_ptr<ROMol> operator()(
      const std::vector<size_t> &which) const = 0;
  //! initializes this operation to work on a particular molecule
  virtual void initFromMol(const ROMol &mol) = 0;
  //! polymorphic copy
  virtual std::unique_ptr<MolEnumeratorOp> copy() const = 0;
};

//! Molecule enumeration operation corresponding to position variation bonds
/*! This uses ATTACH and ENDPTS properties on bonds and requires that the bond
 * has one dummy atom (which will be discarded). The other atom of the bond will
 * be connected to the atoms listed in the ENDPTS property
 */
class RDKIT_MOLENUMERATOR_EXPORT PositionVariationOp : public MolEnumeratorOp {
 public:
  PositionVariationOp(){};
  PositionVariationOp(const std::shared_ptr<ROMol> mol) : dp_mol(mol) {
    PRECONDITION(mol, "bad molecule");
    initFromMol();
  };
  PositionVariationOp(const ROMol &mol) : dp_mol(new ROMol(mol)) {
    initFromMol();
  };
  PositionVariationOp(const PositionVariationOp &other)
      : dp_mol(other.dp_mol), d_variationPoints(other.d_variationPoints){};
  PositionVariationOp &operator=(const PositionVariationOp &other) {
    if (&other == this) {
      return *this;
    }
    dp_mol = other.dp_mol;
    d_variationPoints = other.d_variationPoints;
    return *this;
  };
  //! \override
  std::vector<size_t> getVariationCounts() const override;

  //! \override
  std::unique_ptr<ROMol> operator()(
      const std::vector<size_t> &which) const override;

  //! \override
  void initFromMol(const ROMol &mol) override;

  //! \override
  std::unique_ptr<MolEnumeratorOp> copy() const override {
    return std::unique_ptr<MolEnumeratorOp>(new PositionVariationOp(*this));
  }

 private:
  std::shared_ptr<ROMol> dp_mol{nullptr};
  std::vector<std::pair<unsigned int, std::vector<unsigned int>>>
      d_variationPoints{};
  std::vector<size_t> d_dummiesAtEachPoint{};
  void initFromMol();
};

//! Molecule enumeration operation corresponding to LINKNODES
/*!
 */
class RDKIT_MOLENUMERATOR_EXPORT LinkNodeOp : public MolEnumeratorOp {
 public:
  LinkNodeOp(){};
  LinkNodeOp(const std::shared_ptr<ROMol> mol) : dp_mol(mol) {
    PRECONDITION(mol, "bad molecule");
    initFromMol();
  };
  LinkNodeOp(const ROMol &mol) : dp_mol(new ROMol(mol)) { initFromMol(); };
  LinkNodeOp(const LinkNodeOp &other)
      : dp_mol(other.dp_mol),
        dp_frame(other.dp_frame),
        d_countAtEachPoint(other.d_countAtEachPoint),
        d_variations(other.d_variations),
        d_pointRanges(other.d_pointRanges),
        d_isotopeMap(other.d_isotopeMap){};
  LinkNodeOp &operator=(const LinkNodeOp &other) {
    if (&other == this) {
      return *this;
    }
    dp_mol = other.dp_mol;
    dp_frame = other.dp_frame;
    d_countAtEachPoint = other.d_countAtEachPoint;
    d_variations = other.d_variations;
    d_pointRanges = other.d_pointRanges;
    d_isotopeMap = other.d_isotopeMap;
    return *this;
  };
  //! \override
  std::vector<size_t> getVariationCounts() const override;

  //! \override
  std::unique_ptr<ROMol> operator()(
      const std::vector<size_t> &which) const override;

  //! \override
  void initFromMol(const ROMol &mol) override;

  //! \override
  std::unique_ptr<MolEnumeratorOp> copy() const override {
    return std::unique_ptr<MolEnumeratorOp>(new LinkNodeOp(*this));
  }

 private:
  std::shared_ptr<ROMol> dp_mol{nullptr};
  std::shared_ptr<RWMol> dp_frame{nullptr};
  std::vector<size_t> d_countAtEachPoint{};
  std::vector<std::tuple<unsigned, unsigned, unsigned>> d_variations;
  std::vector<std::pair<unsigned, unsigned>> d_pointRanges;
  std::map<unsigned, unsigned> d_isotopeMap;

  void initFromMol();
};

//! Parameters used to control the molecule enumeration
struct RDKIT_MOLENUMERATOR_EXPORT MolEnumeratorParams {
  bool sanitize = false;
  size_t maxToEnumerate = 1000;
  bool doRandom = false;  //< not yet implemented
  int randomSeed = -1;    //< not yet implemented
  std::shared_ptr<MolEnumeratorOp> dp_operation;
};

//! Returns a MolBundle containing the molecules resulting from applying the
//! operator contained in \c params to \c mol.
RDKIT_MOLENUMERATOR_EXPORT MolBundle
enumerate(const ROMol &mol, const MolEnumeratorParams &params);
}  // namespace MolEnumerator
}  // namespace RDKit

#endif
