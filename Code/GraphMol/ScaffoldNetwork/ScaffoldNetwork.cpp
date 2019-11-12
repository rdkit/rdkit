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
#include <boost/shared_ptr.hpp>
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
    rxn->initReactantMatchers();
    bondBreakersRxns.push_back(rxn);
  }
};

namespace detail {

void updateMolProps(RWMol &mol, const ScaffoldNetworkParams &params) {
  RDUNUSED_PARAM(params);
  MolOps::sanitizeMol(mol);
}

std::vector<std::pair<std::string, ROMOL_SPTR>> getMolFragments(
    const ROMol &mol, const ScaffoldNetworkParams &params) {
  std::vector<std::pair<std::string, ROMOL_SPTR>> res;
  std::deque<ROMOL_SPTR> stack;
  stack.push_back(ROMOL_SPTR(new ROMol(mol)));
  while (!stack.empty()) {
    auto wmol = stack.front();
    stack.pop_front();
    auto parentSmi = MolToSmiles(*wmol);
    for (auto rxn : params.bondBreakersRxns) {
      auto ps = rxn->runReactant(wmol, 0);
      for (auto p : ps) {
        updateMolProps(*static_cast<RWMol *>(p[0].get()), params);
        stack.push_back(p[0]);
        res.push_back(std::make_pair(parentSmi, p[0]));
        if (!params.keepOnlyFirstFragment) {
          updateMolProps(*static_cast<RWMol *>(p[1].get()), params);
          stack.push_back(p[1]);
          res.push_back(std::make_pair(parentSmi, p[1]));
        }
      }
    }
  }
  return res;
}

ROMol *makeScaffoldGeneric(const ROMol &mol, bool doAtoms, bool doBonds) {
  RWMol *res = new RWMol(mol);
  if (doAtoms) {
    for (auto atom : res->atoms()) {
      atom->setAtomicNum(0);
      atom->setNumExplicitHs(0);
      atom->setNoImplicit(false);
    }
  }
  if (doBonds) {
    for (auto bond : res->bonds()) {
      bond->setBondType(Bond::SINGLE);
      bond->getBeginAtom()->setIsAromatic(false);
      bond->getEndAtom()->setIsAromatic(false);
    }
  }
  return static_cast<ROMol *>(res);
}
ROMol *removeAttachmentPoints(const ROMol &mol,
                              const ScaffoldNetworkParams &params) {
  RDUNUSED_PARAM(params);
  RWMol *res = new RWMol(mol);
  for (unsigned int i = 1; i <= mol.getNumAtoms(); ++i) {
    auto atom = res->getAtomWithIdx(mol.getNumAtoms() - i);
    if (!atom->getAtomicNum() && atom->getDegree() == 1) {
      res->removeAtom(atom);
    }
  }
  return static_cast<ROMol *>(res);
}
ROMol *pruneMol(const ROMol &mol, const ScaffoldNetworkParams &params) {
  RDUNUSED_PARAM(params);
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
namespace {
template <typename T, typename V>
size_t addEntryIfMissing(T &vect, const V &e) {
  auto viter = std::find(vect.begin(), vect.end(), e);
  size_t res;
  if (viter == vect.end()) {
    vect.push_back(e);
    res = vect.size() - 1;
  } else {
    res = viter - vect.begin();
  }
  return res;
}

}  // namespace
void addMolToNetwork(const ROMol &mol, ScaffoldNetwork &network,
                     const ScaffoldNetworkParams &params) {
  auto ismi = MolToSmiles(mol);
  boost::shared_ptr<ROMol> fmol(flattenMol(mol, params));
  if (params.pruneBeforeFragmenting) {
    boost::shared_ptr<ROMol> pmol(pruneMol(*fmol, params));
    fmol.swap(pmol);
  }
  auto fsmi = MolToSmiles(*fmol);
  if (ismi != fsmi) {
    auto iidx = addEntryIfMissing(network.nodes, ismi);
    auto fidx = addEntryIfMissing(network.nodes, fsmi);
    network.edges.push_back({iidx, fidx, EdgeType::Initialize});
  }

  auto frags = getMolFragments(*fmol, params);
  for (auto frag : frags) {
    auto smi = frag.first;
    auto fsmi = MolToSmiles(*frag.second);
    auto iidx = addEntryIfMissing(network.nodes, smi);
    auto fidx = addEntryIfMissing(network.nodes, fsmi);
    addEntryIfMissing(network.edges,
                      NetworkEdge({iidx, fidx, EdgeType::Fragment}));

    if (params.includeGenericScaffolds) {
      bool doAtoms = true;
      bool doBonds = false;
      std::unique_ptr<ROMol> gmol(
          makeScaffoldGeneric(*frag.second, doAtoms, doBonds));
      auto gsmi = MolToSmiles(*gmol);
      auto gidx = addEntryIfMissing(network.nodes, gsmi);
      addEntryIfMissing(network.edges,
                        NetworkEdge({fidx, gidx, EdgeType::Generic}));
      if (params.includeScaffoldsWithoutAttachments) {
        std::unique_ptr<ROMol> amol(removeAttachmentPoints(*gmol, params));
        auto asmi = MolToSmiles(*amol);
        auto aidx = addEntryIfMissing(network.nodes, asmi);
        addEntryIfMissing(
            network.edges,
            NetworkEdge({gidx, aidx, EdgeType::RemoveAttachment}));
      }
      if (params.includeGenericBondScaffolds) {
        bool doAtoms = true;
        bool doBonds = true;
        std::unique_ptr<ROMol> gbmol(
            makeScaffoldGeneric(*frag.second, doAtoms, doBonds));
        auto gbsmi = MolToSmiles(*gbmol);
        auto gbidx = addEntryIfMissing(network.nodes, gbsmi);
        addEntryIfMissing(network.edges,
                          NetworkEdge({fidx, gbidx, EdgeType::GenericBond}));
        if (params.includeScaffoldsWithoutAttachments) {
          std::unique_ptr<ROMol> amol(removeAttachmentPoints(*gbmol, params));
          auto asmi = MolToSmiles(*amol);
          auto aidx = addEntryIfMissing(network.nodes, asmi);
          addEntryIfMissing(
              network.edges,
              NetworkEdge({gbidx, aidx, EdgeType::RemoveAttachment}));
        }
      }
    }
    if (params.includeScaffoldsWithoutAttachments) {
      std::unique_ptr<ROMol> amol(removeAttachmentPoints(*frag.second, params));
      auto asmi = MolToSmiles(*amol);
      auto aidx = addEntryIfMissing(network.nodes, asmi);
      addEntryIfMissing(network.edges,
                        NetworkEdge({fidx, aidx, EdgeType::RemoveAttachment}));
    }
  }
}
}  // namespace detail

template <typename T>
void updateScaffoldNetwork(const T &mols, ScaffoldNetwork &network,
                           const ScaffoldNetworkParams &params) {
  for (const auto &mol : mols) {
    detail::addMolToNetwork(*mol, network, params);
  }
}

template RDKIT_SCAFFOLDNETWORK_EXPORT void updateScaffoldNetwork(
    const std::vector<ROMOL_SPTR> &ms, ScaffoldNetwork &network,
    const ScaffoldNetworkParams &params);

}  // namespace ScaffoldNetwork
}  // namespace RDKit