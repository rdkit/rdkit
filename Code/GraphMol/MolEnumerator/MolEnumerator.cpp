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
    if (variations.size() >= maxToEnumerate) {
      return;
    }
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

namespace detail {
const std::string idxPropName = "_enumeratorOrigIdx";
void preserveOrigIndices(ROMol &mol) {
  mol.setProp(idxPropName, 1);
  for (auto atom : mol.atoms()) {
    atom->setProp(idxPropName, atom->getIdx());
  }
}
void removeOrigIndices(ROMol &mol) {
  for (auto atom : mol.atoms()) {
    atom->clearProp(idxPropName);
  }
  mol.clearProp(idxPropName);
}
}  // namespace detail

namespace {
void clearReactionProps(ROMol &mol) {
  bool includeRings = false;
  mol.clearComputedProps(includeRings);
  for (auto atom : mol.atoms()) {
    atom->clearProp(common_properties::reactantAtomIdx);
    atom->clearProp(common_properties::reactantIdx);
    atom->clearProp("was_dummy");
    atom->clearProp(common_properties::reactionMapNum);
  }
}
}  // namespace

MolBundle enumerate(const ROMol &mol,
                    const std::vector<MolEnumeratorParams> &paramLists) {
  std::unique_ptr<MolBundle> accum{new MolBundle()};
  boost::shared_ptr<ROMol> molCpy{new ROMol(mol)};
  detail::preserveOrigIndices(*molCpy);
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
        RDKit::ROMOL_SPTR mcp(new ROMol(*tmol));
        clearReactionProps(*mcp);
        thisRound->addMol(mcp);
        continue;
      }
      std::vector<std::vector<size_t>> variations;
      enumerateVariations(variations, variationCounts, params);
      for (const auto &variation : variations) {
        variationsFound = true;
        auto newMol = (*op)(variation);
        newMol->updatePropertyCache(false);
        clearReactionProps(*newMol);
        thisRound->addMol(ROMOL_SPTR(newMol.release()));
      }
    }
    accum.swap(thisRound);
  }
  if (!variationsFound) {
    return MolBundle();
  }
  for (auto rmol : accum->getMols()) {
    detail::removeOrigIndices(*rmol);
  }
  return *accum;
}

MolBundle enumerate(const ROMol &mol, size_t maxPerOperation) {
  std::vector<MolEnumeratorParams> paramsList;

  // position variation first
  MolEnumerator::MolEnumeratorParams posVariationParams;
  auto pvOp = new MolEnumerator::PositionVariationOp();
  posVariationParams.dp_operation =
      std::shared_ptr<MolEnumerator::MolEnumeratorOp>(pvOp);
  if (maxPerOperation > 0) {
    posVariationParams.maxToEnumerate = maxPerOperation;
  }
  paramsList.push_back(posVariationParams);

  // now repeat units
  MolEnumerator::MolEnumeratorParams repeatUnitParams;
  auto sruOp = new MolEnumerator::RepeatUnitOp();
  repeatUnitParams.dp_operation =
      std::shared_ptr<MolEnumerator::MolEnumeratorOp>(sruOp);
  if (maxPerOperation > 0) {
    sruOp->d_maxNumRounds = maxPerOperation;
    repeatUnitParams.maxToEnumerate = maxPerOperation;
  }
  paramsList.push_back(repeatUnitParams);

  // linknodes last because we can only enumerate mols with a single
  // fragment
  MolEnumerator::MolEnumeratorParams linkParams;
  auto lpOp = new MolEnumerator::LinkNodeOp();
  linkParams.dp_operation =
      std::shared_ptr<MolEnumerator::MolEnumeratorOp>(lpOp);
  if (maxPerOperation > 0) {
    linkParams.maxToEnumerate = maxPerOperation;
  }
  paramsList.push_back(linkParams);

  return enumerate(mol, paramsList);
}

}  // namespace MolEnumerator

}  // namespace RDKit
