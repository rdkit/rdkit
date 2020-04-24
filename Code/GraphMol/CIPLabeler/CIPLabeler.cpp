//
//
//  Copyright (C) 2020 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <algorithm>
#include <memory>

#include <GraphMol/RDKitBase.h>

#include "CIPLabeler.h"
#include "CIPMol.h"
#include "configs/Sp2Bond.h"
#include "configs/Tetrahedral.h"

#include "rules/Rules.hpp"
#include "rules/Rule1a.h"
#include "rules/Rule1b.h"
#include "rules/Rule2.h"
#include "rules/Rule3.h"
#include "rules/Rule4a.h"
#include "rules/Rule4b.h"
#include "rules/Rule4c.h"
#include "rules/Rule5New.h"
#include "rules/Rule6.h"

namespace RDKit {
namespace CIPLabeler {

namespace {

// constitutional rules
const Rules constitutional_rules({new Rule1a, new Rule1b, new Rule2});

// all rules (require aux calc)
const Rules all_rules({new Rule1a, new Rule1b, new Rule2, new Rule3, new Rule4a,
                       new Rule4b, new Rule4c, new Rule5New, new Rule6});

std::vector<std::unique_ptr<Configuration>> findConfigs(const CIPMol &mol) {
  std::vector<std::unique_ptr<Configuration>> configs;

  for (auto &atom : mol.atoms()) {
    auto chiraltag = atom->getChiralTag();
    if (chiraltag == Atom::CHI_TETRAHEDRAL_CW ||
        chiraltag == Atom::CHI_TETRAHEDRAL_CCW) {
      std::vector<Atom *> carriers;
      carriers.reserve(4);
      for (auto &bond : mol.getBonds(atom)) {
        auto nbr = bond->getOtherAtom(atom);
        carriers.push_back(nbr);
      }
      if (carriers.size() < 4) {
        // Implicit H -- use the central atom instead of a dummy H
        carriers.push_back(atom);
      }

      auto atom_cfg = -1;
      if (chiraltag == Atom::CHI_TETRAHEDRAL_CW) {
        atom_cfg = Tetrahedral::RIGHT;
      } else {
        atom_cfg = Tetrahedral::LEFT;
      }

      std::unique_ptr<Tetrahedral> cfg{
          new Tetrahedral(mol, atom, std::move(carriers), atom_cfg)};
      configs.push_back(std::move(cfg));
    }
  }

  for (auto i = 0u; i < mol.getNumBonds(); ++i) {
    auto bond = mol.getBond(i);

    int bond_cfg = -1;
    switch (bond->getStereo()) {
    case Bond::STEREOE:
    case Bond::STEREOTRANS:
      bond_cfg = Sp2Bond::OPPOSITE;
      break;
    case Bond::STEREOZ:
    case Bond::STEREOCIS:
      bond_cfg = Sp2Bond::TOGETHER;
      break;
    default:
      continue;
    }

    auto stereo_atoms = MolOps::findStereoAtoms(bond);
    std::vector<Atom *> anchors{
        {mol.getAtom(stereo_atoms[0]), mol.getAtom(stereo_atoms[1])}};
    std::vector<Atom *> atoms{{bond->getBeginAtom(), bond->getEndAtom()}};

    std::unique_ptr<Sp2Bond> cfg(
        new Sp2Bond(mol, bond, std::move(atoms), std::move(anchors), bond_cfg));
    configs.push_back(std::move(cfg));
  }

  return configs;
}

bool labelAux(const std::vector<std::unique_ptr<Configuration>> &configs,
              const Rules &rules,
              const std::unique_ptr<Configuration> &center) {
  using Node_Cfg_Pair = std::pair<Node *, Configuration *>;
  std::vector<Node_Cfg_Pair> aux;

  auto &digraph = center->getDigraph();
  for (const auto &config : configs) {
    if (config == center) {
      continue;
    }
    // FIXME: specific to each descriptor
    const auto &foci = config->getFoci();
    for (const auto &node : digraph.getNodes(foci[0])) {
      if (node->isDuplicate()) {
        continue;
      }
      auto low = node;
      if (foci.size() == 2) {
        for (const auto &edge : node->getEdges(foci[1])) {
          const auto &other_node = edge->getOther(node);
          if (other_node->getDistance() < node->getDistance())
            low = other_node;
        }
      }
      if (!low->isDuplicate()) {
        aux.emplace_back(low, config.get());
      }
    }
  }

  auto pair_cmp = [](const Node_Cfg_Pair &a, const Node_Cfg_Pair &b) {
    return a.first->getDistance() > b.first->getDistance();
  };
  std::sort(aux.begin(), aux.end(), pair_cmp);

  boost::unordered_map<Node *, Descriptor> queue;
  int prev = std::numeric_limits<int>::max();
  for (const auto &e : aux) {
    const auto &node = e.first;

    if (node->getDistance() < prev) {
      for (const auto &e2 : queue) {
        e2.first->setAux(e2.second);
      }
      queue.clear();
      prev = node->getDistance();
    }
    const auto &config = e.second;
    auto label = config->label(node, digraph, rules);
    queue.emplace(node, label);
  }

  for (const auto &e : queue) {
    e.first->setAux(e.second);
  }

  return true;
}

void label(CIPMol &mol,
           const std::vector<std::unique_ptr<Configuration>> &configs) {

  for (const auto &conf : configs) {
    auto desc = conf->label(constitutional_rules);
    if (desc != Descriptor::UNKNOWN) {
      conf->setPrimaryLabel(mol, desc);
    } else {
      if (labelAux(configs, all_rules, conf)) {
        desc = conf->label(all_rules);

        if (desc != Descriptor::UNKNOWN) {
          conf->setPrimaryLabel(mol, desc);
        }
      }
    }
  }
}

} // namespace

void assignCIPLabels(ROMol &mol) {
  auto cipmol = CIPMol(&mol);

  auto configs = findConfigs(cipmol);
  label(cipmol, configs);
}
} // namespace CIPLabeler
} // end of namespace RDKit
