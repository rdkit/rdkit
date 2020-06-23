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
namespace MolEnumerator {

//! abstract base class
class RDKIT_MOLENUMERATOR_EXPORT MolEnumeratorOp {
 public:
  MolEnumeratorOp(){};
  ~MolEnumeratorOp(){};
  virtual std::vector<size_t> getVariationCounts() const = 0;
  // caller owns the memory
  virtual ROMol *operator()(const std::vector<size_t> &which) const = 0;
};

class RDKIT_MOLENUMERATOR_EXPORT PositionVariationOp : public MolEnumeratorOp {
 public:
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
  };
  std::vector<size_t> getVariationCounts() const override;

  ROMol *operator()(const std::vector<size_t> &which) const override;

 private:
  std::shared_ptr<ROMol> dp_mol;
  std::vector<std::pair<unsigned int, std::vector<unsigned int>>>
      d_variationPoints;
  void initFromMol();
};

struct RDKIT_MOLENUMERATOR_EXPORT MolEnumeratorParams {
  bool sanitize = false;
  size_t maxToEnumerate = 0;
  bool doRandom = false;
  int randomSeed = -1;
  std::vector<MolEnumeratorOp> operations;
};

}  // namespace MolEnumerator
}  // namespace RDKit

#endif