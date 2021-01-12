//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MolEnumerator.h"
#include <RDGeneral/Exceptions.h>

namespace RDKit {
namespace MolEnumerator {

namespace {

//! recursively builds the variations
void getVariations(size_t level, std::vector<size_t> base,
                   std::vector<std::vector<size_t>> &variations,
                   const std::vector<size_t> &variationCounts,
                   size_t maxToEnumerate, bool doRandom) {
  PRECONDITION(level < variationCounts.size(), "bad recursion");
  for (size_t i = 0; i < variationCounts[level]; ++i) {
    base[level] = i;
    if (level + 1 >= variationCounts.size()) {
      // at the bottom of the recursion
      variations.push_back(base);
    } else {
      getVariations(level + 1, base, variations, variationCounts,
                    maxToEnumerate, doRandom);
    }
    if (variations.size() >= maxToEnumerate) return;
  }
}
void enumerateVariations(std::vector<std::vector<size_t>> &variations,
                         const std::vector<size_t> &variationCounts,
                         const MolEnumeratorParams &params) {
  if (params.doRandom) {
    UNDER_CONSTRUCTION("random enumeration not yet supported");
  }
  variations.clear();
  std::vector<size_t> base(variationCounts.size(), 0);
  getVariations(0, base, variations, variationCounts, params.maxToEnumerate,
                params.doRandom);
}
}  // namespace

MolBundle enumerate(const ROMol &mol,
                    const std::vector<MolEnumeratorParams> &paramLists) {
  if (paramLists.empty()) {
    return MolBundle();
  }
  std::unique_ptr<MolBundle> accum{new MolBundle()};
  boost::shared_ptr<ROMol> molCpy{new ROMol(mol)};
  accum->addMol(molCpy);
  bool variationsFound = false;
  for (const auto &params : paramLists) {
    if (!params.dp_operation) {
      throw ValueErrorException("MolEnumeratorParams has no operation set");
    }
    std::unique_ptr<MolBundle> thisRound{new MolBundle()};
    for (const auto &tmol : accum->getMols()) {
      // copy the op since we will modify it:
      auto op = params.dp_operation->copy();
      op->initFromMol(*tmol);
      auto variationCounts = op->getVariationCounts();
      if (variationCounts.empty()) {
        thisRound->addMol(tmol);
        continue;
      }
      std::vector<std::vector<size_t>> variations;
      enumerateVariations(variations, variationCounts, params);
      for (const auto &variation : variations) {
        variationsFound = true;
        auto newMol = (*op)(variation);
        newMol->updatePropertyCache(false);
        thisRound->addMol(ROMOL_SPTR(newMol.release()));
      }
    }
    accum.swap(thisRound);
  }
  if (!variationsFound) {
    return MolBundle();
  }
  return *accum;
}

MolBundle enumerate(const ROMol &mol) {
  std::vector<MolEnumeratorParams> paramsList;
  MolEnumerator::MolEnumeratorParams posVariationParams;
  posVariationParams.dp_operation =
      std::shared_ptr<MolEnumerator::MolEnumeratorOp>(
          new MolEnumerator::PositionVariationOp());
  paramsList.push_back(posVariationParams);
  MolEnumerator::MolEnumeratorParams linkParams;
  linkParams.dp_operation = std::shared_ptr<MolEnumerator::MolEnumeratorOp>(
      new MolEnumerator::LinkNodeOp());
  paramsList.push_back(linkParams);
  return enumerate(mol, paramsList);
}

}  // namespace MolEnumerator

}  // namespace RDKit
