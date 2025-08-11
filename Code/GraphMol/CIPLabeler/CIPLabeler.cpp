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

#include <boost/algorithm/string.hpp>

#include "GraphMol/Chirality.h"
#include "GraphMol/RDKitBase.h"
#include <RDGeneral/ControlCHandler.h>

#include "CIPLabeler.h"
#include "CIPMol.h"
#include "configs/Sp2Bond.h"
#include "configs/Tetrahedral.h"
#include "configs/AtropisomerBond.h"

#include "rules/Rules.h"
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

std::vector<std::unique_ptr<Configuration>> findConfigs(
    CIPMol &mol, const boost::dynamic_bitset<> &atoms,
    const boost::dynamic_bitset<> &bonds) {
  std::vector<std::unique_ptr<Configuration>> configs;

  for (auto index = atoms.find_first(); index != boost::dynamic_bitset<>::npos;
       index = atoms.find_next(index)) {
    auto atom = mol.getAtom(index);
    auto chiraltag = atom->getChiralTag();
    if (chiraltag == Atom::CHI_TETRAHEDRAL_CW ||
        chiraltag == Atom::CHI_TETRAHEDRAL_CCW) {
      std::unique_ptr<Tetrahedral> cfg{new Tetrahedral(mol, atom)};
      configs.push_back(std::move(cfg));
    }
  }

  for (auto index = bonds.find_first(); index != boost::dynamic_bitset<>::npos;
       index = bonds.find_next(index)) {
    auto bond = mol.getBond(index);

    auto bond_cfg = bond->getStereo();
    switch (bond_cfg) {
      case Bond::STEREOE:
        bond_cfg = Bond::STEREOTRANS;
        break;
      case Bond::STEREOZ:
        bond_cfg = Bond::STEREOCIS;
        break;
      default:
        break;
    }
    switch (bond_cfg) {
      case Bond::STEREOTRANS:
      case Bond::STEREOCIS: {
        std::unique_ptr<Sp2Bond> cfg(new Sp2Bond(
            mol, bond, bond->getBeginAtom(), bond->getEndAtom(), bond_cfg));
        configs.push_back(std::move(cfg));
      } break;

      case Bond::STEREOATROPCCW:
      case Bond::STEREOATROPCW: {
        std::unique_ptr<AtropisomerBond> cfgAtrop(new AtropisomerBond(
            mol, bond, bond->getBeginAtom(), bond->getEndAtom(), bond_cfg));
        configs.push_back(std::move(cfgAtrop));
      } break;

      default:
        break;
    }
  }

  return configs;
}

bool labelAux(std::vector<std::unique_ptr<Configuration>> &configs,
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
          if (other_node->getDistance() < node->getDistance()) {
            low = other_node;
          }
        }
      }
      if (!low->isDuplicate()) {
        aux.emplace_back(low, config.get());
      }
    }
  }

  auto farthest = [](const Node_Cfg_Pair &a, const Node_Cfg_Pair &b) {
    return a.first->getDistance() > b.first->getDistance();
  };
  std::sort(aux.begin(), aux.end(), farthest);

  // Using a boost::unordered_map because it is more performant
  // than the STL version.
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

thread_local unsigned int remainingCallCount = 0;

// The chiral centers in current rdkit examples that can be resolved using only
// the constitutional rules average about 8 iterations (the highest count is
// 1039, in one of the examples in the CIP validation suite). We use 2000 as
// threshold to allow some margin.
constexpr unsigned int constitutionalRuleTimeout = 2000;

void label(std::vector<std::unique_ptr<Configuration>> &configs,
           unsigned int maxRecursiveIterations) {
  // First, if the specified number of iterations allows it, run all centers
  // through a fast pass with the constitutional rules allow easy stuff to be
  // resolved.
  for (auto &conf : configs) {
    // Make sure this stereo center has no label
    conf->getFocus()->clearProp(common_properties::_CIPCode);

    remainingCallCount = constitutionalRuleTimeout;
    try {
      auto desc = conf->label(constitutional_rules);
      if (desc != Descriptor::UNKNOWN) {
        conf->setPrimaryLabel(desc);
      }
    } catch (const MaxIterationsExceeded &) {
    }
  }

  // Now, retry everything that hasn't been solved with a more generous
  // threshold
  if (maxRecursiveIterations != 0) {
    remainingCallCount = maxRecursiveIterations;
  } else {
    remainingCallCount = UINT_MAX;  // really big - will never be hit
  }

  // try again on everything that hasn't been resolved yet
  for (const auto &conf : configs) {
    if (conf->getFocus()->hasProp(common_properties::_CIPCode)) {
      // already resolved!
      continue;
    }

    auto desc = conf->label(constitutional_rules);
    if (desc != Descriptor::UNKNOWN) {
      conf->setPrimaryLabel(desc);
    } else {
      if (labelAux(configs, all_rules, conf)) {
        desc = conf->label(all_rules);

        if (desc != Descriptor::UNKNOWN) {
          conf->setPrimaryLabel(desc);
        }
      }
    }
  }
}

}  // namespace

void assignCIPLabels(ROMol &mol, const boost::dynamic_bitset<> &atoms,
                     const boost::dynamic_bitset<> &bonds,
                     unsigned int maxRecursiveIterations) {
  ControlCHandler::reset();

  // reset the mark, for the case that this fails
  mol.clearProp(common_properties::_CIPComputed);
  CIPMol cipmol{mol};
  auto configs = findConfigs(cipmol, atoms, bonds);

  try {
    label(configs, maxRecursiveIterations);
  } catch (const ControlCCaught &) {
  }
  if (ControlCHandler::getGotSignal()) {
    BOOST_LOG(rdWarningLog)
        << "Interrupted, cancelling CIP label calculation" << std::endl;
    return;
  }

  const bool computed = true;
  mol.setProp(common_properties::_CIPComputed, true, computed);
}

void assignCIPLabels(ROMol &mol, unsigned int maxRecursiveIterations) {
  boost::dynamic_bitset<> atoms(mol.getNumAtoms());
  boost::dynamic_bitset<> bonds(mol.getNumBonds());
  atoms.set();
  bonds.set();
  assignCIPLabels(mol, atoms, bonds, maxRecursiveIterations);
}

}  // namespace CIPLabeler

namespace CIPLabeler_detail {

bool decrementRemainingCallCountAndCheck() {
  return (--CIPLabeler::remainingCallCount) > 0;
}

}  // namespace CIPLabeler_detail
}  // namespace RDKit
