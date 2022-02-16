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
#include <GraphMol/FileParsers/MolSGroupParsing.h>

namespace RDKit {
namespace MolEnumerator {

void PositionVariationOp::initFromMol(const ROMol &mol) {
  dp_mol.reset(new ROMol(mol));
  initFromMol();
}
void PositionVariationOp::initFromMol() {
  d_variationPoints.clear();
  if (!dp_mol) {
    return;
  }
  if (!dp_mol->hasProp(detail::idxPropName)) {
    detail::preserveOrigIndices(*dp_mol);
  }
  for (const auto bond : dp_mol->bonds()) {
    std::string endpts;
    std::string attach;
    if (bond->getPropIfPresent(common_properties::_MolFileBondEndPts, endpts) &&
        bond->getPropIfPresent(common_properties::_MolFileBondAttach, attach) &&
        attach == "ANY") {
      const Atom *atom = bond->getBeginAtom();
      if (atom->getAtomicNum() == 0) {
        atom = bond->getEndAtom();
        if (atom->getAtomicNum() == 0) {
          // marvin sketch seems to place the position-variation dummy at the
          // beginning of the bond, so we're going to favor taking the end atom.
          // In case other tools construct this differently, we have an
          // exception to that if the end atom is an AtomNull query and the
          // beginning atom is not one.
          if (atom->hasQuery() &&
              atom->getQuery()->getDescription() == "AtomNull" &&
              bond->getBeginAtom()->hasQuery() &&
              bond->getBeginAtom()->getQuery()->getDescription() !=
                  "AtomNull") {
            atom = bond->getBeginAtom();
          }
        }
      }
      d_dummiesAtEachPoint.push_back(bond->getOtherAtomIdx(atom->getIdx()));
      std::vector<unsigned int> oats =
          RDKit::SGroupParsing::ParseV3000Array<unsigned int>(
              endpts, dp_mol->getNumAtoms(), false);
      // decrement the indices and do error checking and whatever additional
      // cleanup is required:
      for (auto &oat : oats) {
        if (oat == 0 || oat > dp_mol->getNumAtoms()) {
          throw ValueErrorException("Bad variation point index");
        }
        --oat;
        // github #4381: if we're connecting to an aromatic heteroatom which
        // has implicit Hs, we should remove those
        auto attachAtom = dp_mol->getAtomWithIdx(oat);
        if (attachAtom->getIsAromatic() && attachAtom->getAtomicNum() != 6) {
          attachAtom->setNumExplicitHs(0);
        }
      }
      d_variationPoints.push_back(std::make_pair(atom->getIdx(), oats));
    }
  }
}

std::vector<size_t> PositionVariationOp::getVariationCounts() const {
  std::vector<size_t> res(d_variationPoints.size());
  std::transform(
      d_variationPoints.begin(), d_variationPoints.end(), res.begin(),
      [](std::pair<unsigned int, std::vector<unsigned int>> pr) -> size_t {
        return pr.second.size();
      });
  return res;
}

std::unique_ptr<ROMol> PositionVariationOp::operator()(
    const std::vector<size_t> &which) const {
  PRECONDITION(dp_mol, "no molecule");
  if (which.size() != d_variationPoints.size()) {
    throw ValueErrorException("bad element choice in enumeration");
  }
  // a bit of quick error checking before starting the real work
  for (unsigned int i = 0; i < d_variationPoints.size(); ++i) {
    if (which[i] >= d_variationPoints[i].second.size()) {
      throw ValueErrorException("bad element value in enumeration");
    }
  }
  RWMol *res = new RWMol(*dp_mol);
  for (unsigned int i = 0; i < d_variationPoints.size(); ++i) {
    const auto tpl = d_variationPoints[i];
    auto begAtomIdx = tpl.first;
    auto endAtomIdx = tpl.second[which[i]];
    // do we already have a bond?
    if (res->getBondBetweenAtoms(begAtomIdx, endAtomIdx)) {
      continue;
    }
    res->addBond(begAtomIdx, endAtomIdx, Bond::BondType::SINGLE);
  }
  // now remove the dummies:
  res->beginBatchEdit();
  for (auto idx : d_dummiesAtEachPoint) {
    res->removeAtom(idx);
  }
  res->commitBatchEdit();
  return std::unique_ptr<ROMol>(static_cast<ROMol *>(res));
}

}  // namespace MolEnumerator

}  // namespace RDKit
