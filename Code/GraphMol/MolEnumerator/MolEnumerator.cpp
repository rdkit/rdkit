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

MolBundle enumerate(const ROMol &mol, const MolEnumeratorParams &params) {
  MolBundle res;
  PRECONDITION(params.dp_operation, "no operation set");
  // copy the op since we will modify it:
  auto op = params.dp_operation->copy();
  op->initFromMol(mol);
  auto variationCounts = op->getVariationCounts();
  if (variationCounts.empty()) {
    return res;
  }
  std::vector<std::vector<size_t>> variations;
  enumerateVariations(variations, variationCounts, params);
  for (const auto &variation : variations) {
    auto newMol = (*op)(variation);
    newMol->updatePropertyCache(false);
    res.addMol(ROMOL_SPTR(newMol.release()));
  }
  return res;
}

}  // namespace MolEnumerator

}  // namespace RDKit
