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
#include "detail.h"

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
    bondBreakersRxns.emplace_back(rxn);
  }
};

namespace detail {

void updateMolProps(RWMol &mol) { MolOps::sanitizeMol(mol); }

ROMol *makeScaffoldGeneric(const ROMol &mol, bool doAtoms, bool doBonds) {
  RWMol *res = new RWMol(mol);
  if (doAtoms) {
    for (auto atom : res->atoms()) {
      atom->setAtomicNum(0);
      atom->setNumExplicitHs(0);
      atom->setNoImplicit(false);
      atom->setIsotope(0);
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
ROMol *removeAttachmentPoints(const ROMol &mol, const ScaffoldNetworkParams &) {
  RWMol *res = new RWMol(mol);
  res->beginBatchEdit();
  for (const auto atom : res->atoms()) {
    if (!atom->getAtomicNum() && atom->getDegree() == 1) {
      // if we're removing a neighbor from an aromatic heteroatom,
      // don't forget to set the H count on that atom:
      auto nbri = res->getAtomNeighbors(atom);
      auto nbr = (*res)[*nbri.first];
      if (nbr->getIsAromatic() && nbr->getAtomicNum() > 0 &&
          nbr->getAtomicNum() != 6) {
        nbr->setNoImplicit(true);
        nbr->setNumExplicitHs(1);
      }
      res->removeAtom(atom);
    }
  }
  res->commitBatchEdit();
  return static_cast<ROMol *>(res);
}
ROMol *pruneMol(const ROMol &mol, const ScaffoldNetworkParams &) {
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
    if (params.flattenIsotopes) {
      atom->setIsotope(0);
    }
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
std::vector<std::pair<std::string, ROMOL_SPTR>> getMolFragments(
    const ROMol &mol, const ScaffoldNetworkParams &params) {
  std::vector<std::pair<std::string, ROMOL_SPTR>> res;
  std::deque<ROMOL_SPTR> stack;
  stack.push_back(ROMOL_SPTR(new ROMol(mol)));
  std::vector<std::string> seen;
  while (!stack.empty()) {
    auto wmol = stack.front();
    stack.pop_front();
    auto parentSmi = MolToSmiles(*wmol);
    for (auto rxn : params.bondBreakersRxns) {
      auto ps = rxn->runReactant(wmol, 0);
      for (auto p : ps) {
        if (!params.includeScaffoldsWithAttachments) {
          // we're only including scaffolds without attachments, so
          // go ahead and remove attachment points here
          p[0].reset(removeAttachmentPoints(*p[0], params));
        }

        updateMolProps(*static_cast<RWMol *>(p[0].get()));
        auto tsmi0 = MolToSmiles(*p[0]);
        if (std::find(seen.begin(), seen.end(), tsmi0) == seen.end()) {
          stack.push_back(p[0]);
          seen.push_back(tsmi0);
        }
        res.push_back(std::make_pair(parentSmi, p[0]));
        if (!params.keepOnlyFirstFragment) {
          if (!params.includeScaffoldsWithAttachments) {
            // we're only including scaffolds without attachments, so
            // go ahead and remove attachment points here
            p[1].reset(removeAttachmentPoints(*p[1], params));
          }
          updateMolProps(*static_cast<RWMol *>(p[1].get()));
          auto tsmi1 = MolToSmiles(*p[1]);
          if (std::find(seen.begin(), seen.end(), tsmi1) == seen.end()) {
            stack.push_back(p[1]);
            seen.push_back(tsmi1);
          }
          res.push_back(std::make_pair(parentSmi, p[1]));
        }
      }
    }
  }
  return res;
}
namespace {
template <typename T, typename V>
size_t addEntryIfMissing(T &vect, const V &e,
                         std::vector<unsigned> *counts = nullptr) {
  auto viter = std::find(vect.begin(), vect.end(), e);
  size_t res;
  if (viter == vect.end()) {
    vect.push_back(e);
    res = vect.size() - 1;
    if (counts) {
      counts->push_back(0);
    }
  } else {
    res = viter - vect.begin();
  }
  if (counts) {
    (*counts)[res] += 1;
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
  std::set<size_t> entriesFromMol;
  auto fsmi = MolToSmiles(*fmol);
  if (ismi != fsmi) {
    auto iidx = addEntryIfMissing(network.nodes, ismi, &network.counts);
    auto fidx = addEntryIfMissing(network.nodes, fsmi, &network.counts);
    if (params.collectMolCounts) {
      entriesFromMol.insert(iidx);
      entriesFromMol.insert(fidx);
    }

    addEntryIfMissing(network.edges,
                      NetworkEdge({iidx, fidx, EdgeType::Initialize}));

    if (params.includeGenericScaffolds) {
      bool doAtoms = true;
      bool doBonds = false;
      std::unique_ptr<ROMol> gmol(makeScaffoldGeneric(*fmol, doAtoms, doBonds));
      auto gsmi = MolToSmiles(*gmol);
      auto gidx = addEntryIfMissing(network.nodes, gsmi, &network.counts);
      if (params.collectMolCounts) {
        entriesFromMol.insert(gidx);
      }
      addEntryIfMissing(network.edges,
                        NetworkEdge({fidx, gidx, EdgeType::Generic}));
      if (params.includeGenericBondScaffolds) {
        bool doAtoms = true;
        bool doBonds = true;
        std::unique_ptr<ROMol> gbmol(
            makeScaffoldGeneric(*fmol, doAtoms, doBonds));
        auto gbsmi = MolToSmiles(*gbmol);
        auto gbidx = addEntryIfMissing(network.nodes, gbsmi, &network.counts);
        if (params.collectMolCounts) {
          entriesFromMol.insert(gbidx);
        }
        if (gidx != gbidx) {
          addEntryIfMissing(network.edges,
                            NetworkEdge({gidx, gbidx, EdgeType::GenericBond}));
        }
      }
    }
  } else {
    // add the base molecule to the network
    auto midx = addEntryIfMissing(network.nodes, fsmi, &network.counts);
    if (params.collectMolCounts) {
      entriesFromMol.insert(midx);
    }
  }

  auto frags = getMolFragments(*fmol, params);
  for (auto frag : frags) {
    auto smi = frag.first;
    ROMOL_SPTR fragMol = frag.second;

    // the first node is certain to already be in the network, so don't
    // update counts for it
    auto iidx = addEntryIfMissing(network.nodes, smi);
    auto lsmi = MolToSmiles(*fragMol);
    auto lidx = addEntryIfMissing(network.nodes, lsmi, &network.counts);
    if (params.collectMolCounts) {
      entriesFromMol.insert(lidx);
    }

    addEntryIfMissing(network.edges,
                      NetworkEdge({iidx, lidx, EdgeType::Fragment}));

    if (params.includeGenericScaffolds) {
      bool doAtoms = true;
      bool doBonds = false;
      std::unique_ptr<ROMol> gmol(
          makeScaffoldGeneric(*fragMol, doAtoms, doBonds));
      auto gsmi = MolToSmiles(*gmol);
      auto gidx = addEntryIfMissing(network.nodes, gsmi, &network.counts);
      if (params.collectMolCounts) {
        entriesFromMol.insert(gidx);
      }

      addEntryIfMissing(network.edges,
                        NetworkEdge({lidx, gidx, EdgeType::Generic}));
      if (params.includeGenericBondScaffolds) {
        bool doAtoms = true;
        bool doBonds = true;
        std::unique_ptr<ROMol> gbmol(
            makeScaffoldGeneric(*fragMol, doAtoms, doBonds));
        auto gbsmi = MolToSmiles(*gbmol);
        auto gbidx = addEntryIfMissing(network.nodes, gbsmi, &network.counts);
        if (params.collectMolCounts) {
          entriesFromMol.insert(gbidx);
        }

        if (gidx != gbidx) {
          addEntryIfMissing(network.edges,
                            NetworkEdge({gidx, gbidx, EdgeType::GenericBond}));
        }
      }
    }
    if (params.includeScaffoldsWithAttachments &&
        params.includeScaffoldsWithoutAttachments) {
      // we're including both scaffolds without attachments and those with.
      // do the one without now.
      std::unique_ptr<ROMol> amol(removeAttachmentPoints(*fragMol, params));
      auto asmi = MolToSmiles(*amol);
      auto aidx = addEntryIfMissing(network.nodes, asmi, &network.counts);
      if (params.collectMolCounts) {
        entriesFromMol.insert(aidx);
      }

      addEntryIfMissing(network.edges,
                        NetworkEdge({lidx, aidx, EdgeType::RemoveAttachment}));
      if (params.includeGenericScaffolds) {
        bool doAtoms = true;
        bool doBonds = false;
        std::unique_ptr<ROMol> gmol(
            makeScaffoldGeneric(*amol, doAtoms, doBonds));
        auto gsmi = MolToSmiles(*gmol);
        auto gidx = addEntryIfMissing(network.nodes, gsmi, &network.counts);
        if (params.collectMolCounts) {
          entriesFromMol.insert(gidx);
        }

        addEntryIfMissing(network.edges,
                          NetworkEdge({aidx, gidx, EdgeType::Generic}));
        if (params.includeGenericBondScaffolds) {
          bool doAtoms = true;
          bool doBonds = true;
          std::unique_ptr<ROMol> gbmol(
              makeScaffoldGeneric(*amol, doAtoms, doBonds));
          auto gbsmi = MolToSmiles(*gbmol);
          auto gbidx = addEntryIfMissing(network.nodes, gbsmi, &network.counts);
          if (params.collectMolCounts) {
            entriesFromMol.insert(gbidx);
          }

          if (gidx != gbidx) {
            addEntryIfMissing(
                network.edges,
                NetworkEdge({gidx, gbidx, EdgeType::GenericBond}));
          }
        }
      }
    }
  }
  if (params.collectMolCounts) {  // update the molCounts, entriesPerMol has
                                  // all the entries we just touched
    auto maxVal = *entriesFromMol.rbegin();
    if (network.molCounts.size() <= maxVal) {
      network.molCounts.resize(maxVal + 1, 0u);
    }
    for (auto idx : entriesFromMol) {
      network.molCounts[idx]++;
    }
  }
}
}  // namespace detail

template <typename T>
void updateScaffoldNetwork(const T &mols, ScaffoldNetwork &network,
                           const ScaffoldNetworkParams &params) {
  if (!params.includeScaffoldsWithAttachments &&
      !params.includeScaffoldsWithoutAttachments) {
    throw ValueErrorException(
        "must include at least one of scaffolds with attachments or scaffolds "
        "without attachments");
  }
  for (const auto &mol : mols) {
    if (!mol) {
      throw ValueErrorException(
          "updateScaffoldNetwork called with null molecule");
    }
    detail::addMolToNetwork(*mol, network, params);
  }
}

template RDKIT_SCAFFOLDNETWORK_EXPORT void updateScaffoldNetwork(
    const std::vector<ROMOL_SPTR> &ms, ScaffoldNetwork &network,
    const ScaffoldNetworkParams &params);
template RDKIT_SCAFFOLDNETWORK_EXPORT ScaffoldNetwork createScaffoldNetwork(
    const std::vector<ROMOL_SPTR> &ms, const ScaffoldNetworkParams &params);
template RDKIT_SCAFFOLDNETWORK_EXPORT void updateScaffoldNetwork(
    const std::vector<std::shared_ptr<ROMol>> &ms, ScaffoldNetwork &network,
    const ScaffoldNetworkParams &params);
template RDKIT_SCAFFOLDNETWORK_EXPORT ScaffoldNetwork
createScaffoldNetwork(const std::vector<std::shared_ptr<ROMol>> &ms,
                      const ScaffoldNetworkParams &params);

const std::vector<std::string> BRICSDefinitions = {
    "[$([C;D3]([#0,#6,#7,#8])(=O)):1]-;!@[$([O;D2]-;!@[#0,#6,#1]):2]>>[1*]-[*"
    ":1].[3*]-[*:2]",
    "[$([C;D3]([#0,#6,#7,#8])(=O)):1]-;!@[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#"
    "0;!#1]);!$([N;R]@[C;R]=O)]):2]>>[1*]-[*:1].[5*]-[*:2]",
    "[$([C;D3]([#0,#6,#7,#8])(=O)):1]-;!@[$([N;R;$(N(@C(=O))@[C,N,O,S])]):2]>"
    ">[1*]-[*:1].[10*]-[*:2]",
    "[$([O;D2]-;!@[#0,#6,#1]):1]-;!@[$([C;!D1;!$(C=*)]-;!@[#6]):2]>>[3*]-[*:"
    "1].[4*]-[*:2]",
    "[$([O;D2]-;!@[#0,#6,#1]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>"
    "[3*]-[*:1].[13*]-[*:2]",
    "[$([O;D2]-;!@[#0,#6,#1]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[3*]"
    "-[*:1].[14*]-[*:2]",
    "[$([O;D2]-;!@[#0,#6,#1]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[3*]-[*:1].["
    "15*]-[*:2]",
    "[$([O;D2]-;!@[#0,#6,#1]):1]-;!@[$([c;$(c(:c):c)]):2]>>[3*]-[*:1].[16*]-["
    "*:2]",
    "[$([C;!D1;!$(C=*)]-;!@[#6]):1]-;!@[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!"
    "#1]);!$([N;R]@[C;R]=O)]):2]>>[4*]-[*:1].[5*]-[*:2]",
    "[$([C;!D1;!$(C=*)]-;!@[#6]):1]-;!@[$([S;D2](-;!@[#0,#6])):2]>>[4*]-[*:1]"
    ".[11*]-[*:2]",
    "[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):1]-;!@[$("
    "[S;D4]([#6,#0])(=O)(=O)):2]>>[5*]-[*:1].[12*]-[*:2]",
    "[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):1]-;!@[$("
    "[c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[5*]-[*:1].[14*]-[*:2]",
    "[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):1]-;!@[$("
    "[c;$(c(:c):c)]):2]>>[5*]-[*:1].[16*]-[*:2]",
    "[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):1]-;!@[$("
    "[C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>[5*]-[*:1].[13*]-[*:2]",
    "[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):1]-;!@[$("
    "[C;$(C(-;@C)-;@C)]):2]>>[5*]-[*:1].[15*]-[*:2]",
    "[$([C;D3;!R](=O)-;!@[#0,#6,#7,#8]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,"
    "S])]):2]>>[6*]-[*:1].[13*]-[*:2]",
    "[$([C;D3;!R](=O)-;!@[#0,#6,#7,#8]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]"
    "):2]>>[6*]-[*:1].[14*]-[*:2]",
    "[$([C;D3;!R](=O)-;!@[#0,#6,#7,#8]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[6*]"
    "-[*:1].[15*]-[*:2]",
    "[$([C;D3;!R](=O)-;!@[#0,#6,#7,#8]):1]-;!@[$([c;$(c(:c):c)]):2]>>[6*]-[*:"
    "1].[16*]-[*:2]",
    "[$([C;D2,D3]-[#6]):1]=;!@[$([C;D2,D3]-[#6]):2]>>[7*]-[*:1].[7*]-[*:2]",
    "[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]):2]>>"
    "[8*]-[*:1].[9*]-[*:2]",
    "[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([N;R;$(N(@C(=O))@[C,N,O,S])]):2]>>[8*]-"
    "[*:1].[10*]-[*:2]",
    "[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>["
    "8*]-[*:1].[13*]-[*:2]",
    "[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[8*]-"
    "[*:1].[14*]-[*:2]",
    "[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[8*]-[*:1].[15*"
    "]-[*:2]",
    "[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([c;$(c(:c):c)]):2]>>[8*]-[*:1].[16*]-[*"
    ":2]",
    "[$([n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@["
    "N,O,S])]):2]>>[9*]-[*:1].[13*]-[*:2]",
    "[$([n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,"
    "s])]):2]>>[9*]-[*:1].[14*]-[*:2]",
    "[$([n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>["
    "9*]-[*:1].[15*]-[*:2]",
    "[$([n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]):1]-;!@[$([c;$(c(:c):c)]):2]>>[9*]-"
    "[*:1].[16*]-[*:2]",
    "[$([N;R;$(N(@C(=O))@[C,N,O,S])]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S]"
    ")]):2]>>[10*]-[*:1].[13*]-[*:2]",
    "[$([N;R;$(N(@C(=O))@[C,N,O,S])]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):"
    "2]>>[10*]-[*:1].[14*]-[*:2]",
    "[$([N;R;$(N(@C(=O))@[C,N,O,S])]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[10*]-"
    "[*:1].[15*]-[*:2]",
    "[$([N;R;$(N(@C(=O))@[C,N,O,S])]):1]-;!@[$([c;$(c(:c):c)]):2]>>[10*]-[*:"
    "1].[16*]-[*:2]",
    "[$([S;D2](-;!@[#0,#6])):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>["
    "11*]-[*:1].[13*]-[*:2]",
    "[$([S;D2](-;!@[#0,#6])):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[11*]"
    "-[*:1].[14*]-[*:2]",
    "[$([S;D2](-;!@[#0,#6])):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[11*]-[*:1].["
    "15*]-[*:2]",
    "[$([S;D2](-;!@[#0,#6])):1]-;!@[$([c;$(c(:c):c)]):2]>>[11*]-[*:1].[16*]-["
    "*:2]",
    "[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s]"
    ")]):2]>>[13*]-[*:1].[14*]-[*:2]",
    "[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>["
    "13*]-[*:1].[15*]-[*:2]",
    "[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):1]-;!@[$([c;$(c(:c):c)]):2]>>[13*]-"
    "[*:1].[16*]-[*:2]",
    "[$([c;$(c(:[c,n,o,s]):[n,o,s])]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):"
    "2]>>[14*]-[*:1].[14*]-[*:2]",
    "[$([c;$(c(:[c,n,o,s]):[n,o,s])]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[14*]-"
    "[*:1].[15*]-[*:2]",
    "[$([c;$(c(:[c,n,o,s]):[n,o,s])]):1]-;!@[$([c;$(c(:c):c)]):2]>>[14*]-[*:"
    "1].[16*]-[*:2]",
    "[$([C;$(C(-;@C)-;@C)]):1]-;!@[$([c;$(c(:c):c)]):2]>>[15*]-[*:1].[16*]-[*"
    ":2]",
    "[$([c;$(c(:c):c)]):1]-;!@[$([c;$(c(:c):c)]):2]>>[16*]-[*:1].[16*]-[*:"
    "2]"};

ScaffoldNetworkParams getBRICSNetworkParams() {
  ScaffoldNetworkParams res{BRICSDefinitions};
  res.keepOnlyFirstFragment = false;
  return res;
};

// const ScaffoldNetworkParams BRICSNetworkParams;

}  // namespace ScaffoldNetwork
}  // namespace RDKit
