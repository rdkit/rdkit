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
std::vector<Configuration *> findConfigs(CIPMol &mol) {
  std::vector<Configuration *> configs;

  for (auto &atom : mol.atoms()) {
    auto chiraltag = static_cast<int>(atom->getChiralTag());
    if (chiraltag > 0) {
      auto carriers = std::vector<Atom *>();
      carriers.reserve(4);
      for (auto &bond : mol.getBonds(atom)) {
        auto nbr = bond->getOtherAtom(atom);
        carriers.push_back(nbr);
      }
      if (carriers.size() < 4) {
        // Implicit H -- use the central atom instead of a dummy H
        carriers.push_back(atom);
      }

      // RDKit uses SMILES parity, and centres requires SDF parity, which
      // is opposite.
      auto parity = chiraltag ^ 0x3;

      auto cfg = new Tetrahedral(atom, std::move(carriers), parity);
      configs.push_back(cfg);
    }
  }

  for (auto i = 0u; i < mol.getNumBonds(); ++i) {
    auto bond = mol.getBond(i);

    if (bond->getBondType() == Bond::DOUBLE) {
      std::vector<Atom *> atoms{{bond->getBeginAtom(), bond->getEndAtom()}};
      std::vector<Atom *> anchors;
      auto previous_bond_dir = Bond::BondDir::NONE;
      auto bond_cfg = 1;
      for (const auto &atom : atoms) {
        for (auto &nbr_bnd : mol.getBonds(atom)) {
          auto nbr_bond_dir = nbr_bnd->getBondDir();
          if (nbr_bond_dir == Bond::BondDir::ENDDOWNRIGHT ||
              nbr_bond_dir == Bond::BondDir::ENDUPRIGHT) {
            anchors.push_back(nbr_bnd->getOtherAtom(atom));

            // Correct the BondDir to have the same reference
            // for the direction
            if (atom == nbr_bnd->getEndAtom()) {
              if (nbr_bond_dir == Bond::BondDir::ENDDOWNRIGHT) {
                nbr_bond_dir = Bond::BondDir::ENDUPRIGHT;
              } else {
                nbr_bond_dir = Bond::BondDir::ENDDOWNRIGHT;
              }
            }

            bond_cfg += nbr_bond_dir == previous_bond_dir;
            std::swap(previous_bond_dir, nbr_bond_dir);
            break;
          }
        }
      }

      if (anchors.size() == 2) {
        auto cfg =
            new Sp2Bond(bond, std::move(atoms), std::move(anchors), bond_cfg);
        configs.push_back(cfg);
      }
    }
  }

  return configs;
}

bool labelAux(const std::vector<Configuration *> &configs, const Rules *rules,
              const Configuration *center) {
  using Node_Cfg_Pair = std::pair<Node *, Configuration *>;
  auto aux = std::vector<Node_Cfg_Pair>{};

  auto digraph = center->getDigraph();
  for (const auto &config : configs) {
    if (config == center) {
      continue;
    }
    // FIXME: specific to each descriptor
    const auto &foci = config->getFoci();
    for (const auto &node : digraph->getNodes(foci[0])) {
      if (node->isDuplicate()) {
        continue;
      }
      auto low = node;
      if (foci.size() == 2) {
        for (const auto &edge : node->getEdges(foci[1])) {
          if (edge->getOther(node)->getDistance() < node->getDistance())
            low = edge->getOther(node);
        }
      }
      if (!low->isDuplicate()) {
        aux.emplace_back(low, config);
      }
    }
  }

  auto pair_cmp = [](const Node_Cfg_Pair &a, const Node_Cfg_Pair &b) {
    return a.first->getDistance() > b.first->getDistance();
  };
  std::sort(aux.begin(), aux.end(), pair_cmp);

  auto queue = boost::unordered_map<Node *, Descriptor>{};
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
    auto label = config->label(node, digraph.get(), rules);
    queue.emplace(node, label);
  }

  for (const auto &e : queue) {
    e.first->setAux(e.second);
  }

  return true;
}

void label(CIPMol *mol, const std::vector<Configuration *> &configs) {
  // constitutional rules
  const Rules begRules({new Rule1a(mol), new Rule1b(mol), new Rule2(mol)});

  // all rules (require aux calc)
  const Rules allRules({new Rule1a(mol), new Rule1b(mol), new Rule2(mol),
                        new Rule3(mol), new Rule4a(mol), new Rule4b(mol),
                        new Rule4c(mol), new Rule5New(mol), new Rule6(mol)});

  auto finalLabels = boost::unordered_map<Configuration *, Descriptor>{};
  for (const auto &conf : configs) {
    conf->setDigraph(new Digraph(mol));
    auto desc = conf->label(&begRules);
    if (desc != Descriptor::UNKNOWN) {
      conf->setPrimaryLabel(mol, desc);
    } else {
      if (labelAux(configs, &allRules, conf)) {
        desc = conf->label(&allRules);

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
  label(&cipmol, configs);

  // Cleanup
  for (auto &c : configs) {
    delete c;
  }
}
} // namespace CIPLabeler
} // end of namespace RDKit
