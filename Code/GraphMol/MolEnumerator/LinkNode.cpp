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

void LinkNodeOp::initFromMol(const ROMol &mol) {
  dp_mol.reset(new ROMol(mol));
  initFromMol();
}
void LinkNodeOp::initFromMol() {
  // d_variationPoints.clear();
  if (!dp_mol) {
    return;
  }
}

std::vector<size_t> LinkNodeOp::getVariationCounts() const {
  // std::vector<size_t> res(d_variationPoints.size());
  // std::transform(
  //     d_variationPoints.begin(), d_variationPoints.end(), res.begin(),
  //     [](std::pair<unsigned int, std::vector<unsigned int>> pr) -> size_t {
  //       return pr.second.size();
  //     });
  // return res;
  std::vector<size_t> res;
  std::string pval;
  if (!dp_mol->getPropIfPresent(common_properties::molFileLinkNodes, pval)) {
    return res;
  }

  return res;
}

std::unique_ptr<ROMol> LinkNodeOp::operator()(
    const std::vector<size_t> &which) const {
  PRECONDITION(dp_mol, "no molecule");
  // if (which.size() != d_variationPoints.size()) {
  //   throw ValueErrorException("bad element choice in enumeration");
  // }
  RWMol *res = new RWMol(*dp_mol);
  return std::unique_ptr<ROMol>(static_cast<ROMol *>(res));
}

}  // namespace MolEnumerator

}  // namespace RDKit
