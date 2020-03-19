//
//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include <initializer_list>
#include <vector>

#include "Rule5.hpp"
#include "SequenceRule.hpp"

namespace RDKit {
namespace NewCIPLabelling {

namespace {
template <typename Base, typename T>
bool is_instance_of(const T* ptr) {
  // return std::is_base_of<Base, T>::value;
  return dynamic_cast<const Base*>(ptr) != nullptr;
}
}  // namespace

template <typename A, typename B>
class Rules : public SequenceRule<A, B> {
 private:
  static const bool SORT_BRANCHES_WITH_RULE5 = false;
  std::vector<const SequenceRule<A, B>*> rules;

 public:
  Rules() = delete;

  Rules(std::initializer_list<SequenceRule<A, B>*> rules)
      : SequenceRule<A, B>(nullptr) {
    for (auto& rule : rules) {
      add(rule);
    }
  }

  ~Rules() {
    for (auto& rule : rules) {
      delete rule;
    }
  }

  void add(SequenceRule<A, B>* rule) {
    if (rule == nullptr) {
      throw std::runtime_error("No sequence rule provided");
    }
    rules.push_back(rule);
    rule->setSorter(new Sort<A, B>(rules));
  }

  int getNumSubRules() const override { return rules.size(); }

  const Sort<A, B>* getSorter() const override {
    if (this->sorter == nullptr) {
      const_cast<Rules<A, B>*>(this)->setSorter(new Sort<A, B>(this));
    }
    return this->sorter.get();
  }

  const BaseMol<A, B>* getMol() const override {
    const BaseMol<A, B>* res = nullptr;
    for (const auto& rule : rules) {
      res = rule->getMol();
      if (res != nullptr) {
        break;
      }
    }
    return res;
  }

  int compare(const Edge<A, B>* o1, const Edge<A, B>* o2) const override {
    // Try using each rules. The rules will expand the search exhaustively
    // to all child ligands
    for (const auto& rule : rules) {
      // compare expands exhaustively across the whole graph
      int value = rule->recursiveCompare(o1, o2);
      if (value != 0) {
        return value;
      }
    }
    return 0;
  }

  int getComparision(const Edge<A, B>* a, const Edge<A, B>* b,
                     bool deep) const override {
    (void)deep;

    // Try using each rules. The rules will expand the search exhaustively
    // to all child ligands
    for (const auto& rule : rules) {
      // compare expands exhaustively across the whole graph
      int value = rule->recursiveCompare(a, b);

      if (value != 0) {
        return value;
      }
    }

    return 0;
  }
};

}  // namespace NewCIPLabelling
}  // namespace RDKit