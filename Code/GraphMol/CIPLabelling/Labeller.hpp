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

///
/// This is the core of the labelling code. The label() function
/// visits each configuration to be labelled, resolves it, and
/// assigns the proper descriptor.
///

#pragma once

#include <limits>
#include <sstream>
#include <vector>

#include <boost/unordered_map.hpp>

#include "Descriptor.hpp"
#include "Digraph.hpp"
#include "rules/Rule1a.hpp"
#include "rules/Rule1b.hpp"
#include "rules/Rule2.hpp"
#include "rules/Rule3.hpp"
#include "rules/Rule4a.hpp"
#include "rules/Rule4b.hpp"
#include "rules/Rule4c.hpp"
#include "rules/Rule5New.hpp"
#include "rules/Rule6.hpp"
#include "rules/Rules.hpp"

namespace RDKit {
namespace CIPLabelling {

template <typename A, typename B> class BaseMol;

template <typename A, typename B> class Configuration;

namespace {

template <typename A, typename B>
bool labelAux(const std::vector<Configuration<A, B> *> &configs,
              const Rules<A, B> *rules, const Configuration<A, B> *center) {
  using Node_Cfg_Pair = std::pair<Node<A, B> *, Configuration<A, B> *>;
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

  auto queue = boost::unordered_map<Node<A, B> *, Descriptor>{};
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

} // namespace

template <typename A, typename B>
void label(BaseMol<A, B> *mol,
           const std::vector<Configuration<A, B> *> configs) {
  // constitutional rules
  const Rules<A, B> begRules(
      {new Rule1a<A, B>(mol), new Rule1b<A, B>(mol), new Rule2<A, B>(mol)});

  // all rules (require aux calc)
  const Rules<A, B> allRules(
      {new Rule1a<A, B>(mol), new Rule1b<A, B>(mol), new Rule2<A, B>(mol),
       new Rule3<A, B>(mol), new Rule4a<A, B>(mol), new Rule4b<A, B>(mol),
       new Rule4c<A, B>(mol), new Rule5New<A, B>(mol), new Rule6<A, B>(mol)});

  auto finalLabels = boost::unordered_map<Configuration<A, B> *, Descriptor>{};
  for (const auto &conf : configs) {
    conf->setDigraph(new Digraph<A, B>(mol));
    try {
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
    } catch (const std::runtime_error &e) {
      throw;
    }
  }
}

} // namespace CIPLabelling
} // namespace RDKit