//
//  Copyright (C) 2019 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "ScaffoldNetwork.h"
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MolStandardize/Fragment.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <RDGeneral/Exceptions.h>

namespace RDKit {
namespace ScaffoldNetwork {

ScaffoldNetworkParams::ScaffoldNetworkParams(
    const std::vector<std::string> &bondBreakersSmarts) {
  bondBreakersRxns.clear();
  bondBreakersRxns.reserve(bondBreakersSmarts.size());
  for (auto sma : bondBreakersSmarts) {
    std::shared_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(sma));
    if (rxn.get() == nullptr) {
      throw ValueErrorException("could not parse reaction smarts: " + sma);
    }
    if (rxn->getNumReactantTemplates() != 1) {
      throw ValueErrorException(
          "bond breaker reactions must have exactly one reactant");
    }
    if (rxn->getNumProductTemplates() != 2) {
      throw ValueErrorException(
          "bond breaker reactions must have exactly two products");
    }
    bondBreakersRxns.push_back(rxn);
  }
};

namespace detail {
ROMol *pruneMol(const ROMol &mol, const ScaffoldNetworkParams &params){
  ROMol *res = MurckoDecompose(mol);
  res->updatePropertyCache();
  MolOps::fastFindRings(*res);
  return res;
}

ROMol *flattenMol(const ROMol &mol, const ScaffoldNetworkParams &params) {
  RWMol *res;
  if (params.flattenKeepLargest) {
    MolStandardize::LargestFragmentChooser fragmenter;
    res = static_cast<RWMol *>(fragmenter.choose(mol));
  } else {
    res = new RWMol(mol);
  }
  for (auto atom : res->atoms()) {
    if (params.flattenIsotopes) atom->setIsotope(0);
    if (params.flattenChirality) {
      if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
        atom->setChiralTag(Atom::CHI_UNSPECIFIED);
        if (atom->getNoImplicit() &&
            (atom->getAtomicNum() == 6 || atom->getAtomicNum() == 7 ||
             atom->getAtomicNum() == 15 || atom->getAtomicNum() == 16)) {
          atom->setNoImplicit(false);
          atom->setNumExplicitHs(0);
        }
      }
    }
  }
  if (params.flattenChirality) {
    for (auto bond : res->bonds()) {
      bond->setBondDir(Bond::NONE);
      bond->setStereo(Bond::STEREONONE);
    }
  }
  return static_cast<ROMol *>(res);
}
void addMolToNetwork(const ROMol &mol, ScaffoldNetwork &network,
                     const ScaffoldNetworkParams &params) {
  auto ismi = MolToSmiles(mol);
  boost::shared_ptr<ROMol> fmol(flattenMol(mol, params));
}
}  // namespace detail

template <typename T>
void updateScaffoldNetwork(const T &mols, ScaffoldNetwork &network,
                           const ScaffoldNetworkParams &params) {
  for (const auto &mol : mols) {
    addMolToNetwork(*mol, network, params);
  }
}

}  // namespace ScaffoldNetwork
}  // namespace RDKit