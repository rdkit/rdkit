//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#pragma once
#include <map>
#include <vector>
#include <stdexcept>
#include <algorithm>

namespace RDKit {
namespace FMCS {
class DuplicatedSeedCache {
 public:
  typedef bool TValue;
  class TKey {
    std::vector<unsigned int> AtomIdx;  // sorted
    std::vector<unsigned int> BondIdx;  // sorted
   public:
    size_t getNumAtoms() const { return AtomIdx.size(); }
    size_t getNumBonds() const { return BondIdx.size(); }

    void addAtom(unsigned int i) {
      auto it = std::lower_bound(AtomIdx.begin(), AtomIdx.end(), i);
      AtomIdx.insert(it, i);
    }
    void addBond(unsigned int i) {
      auto it = std::lower_bound(BondIdx.begin(), BondIdx.end(), i);
      BondIdx.insert(it, i);
    }

    bool operator==(const TKey& right) const {  // opt.
      return AtomIdx.size() == right.AtomIdx.size() &&
             BondIdx.size() == right.BondIdx.size() &&
             0 == std::memcmp(&AtomIdx[0], &right.AtomIdx[0],
                              AtomIdx.size() * sizeof(unsigned int)) &&
             0 == std::memcmp(&BondIdx[0], &right.BondIdx[0],
                              BondIdx.size() * sizeof(unsigned int));
    }

    bool operator<(const TKey& right) const {
      if (AtomIdx.size() < right.AtomIdx.size()) {
        return true;
      }
      if (AtomIdx.size() > right.AtomIdx.size()) {
        return false;
      }

      if (BondIdx.size() < right.BondIdx.size()) {
        return true;
      }
      if (BondIdx.size() > right.BondIdx.size()) {
        return false;
      }

      // everything is equal -> perform straight comparison
      int diff;
      diff = std::memcmp(&AtomIdx[0], &right.AtomIdx[0],
                         AtomIdx.size() * sizeof(unsigned int));
      if (diff < 0) {
        return true;
      }
      if (diff > 0) {
        return false;
      }
      return std::memcmp(&BondIdx[0], &right.BondIdx[0],
                         BondIdx.size() * sizeof(unsigned int)) < 0;
    }
  };

 private:
  std::map<TKey, TValue> Index;
  size_t MaxAtoms{0};  // max key in the cache for fast failed find
 public:
  DuplicatedSeedCache() {}
  void clear() {
    Index.clear();
    MaxAtoms = 0;
  }

  bool find(const TKey& key, TValue& value) const {
    value = false;
    if (key.getNumAtoms() > MaxAtoms) {
      return false;  // fast check if key greater then max key in the cache
    }

    const auto entryit = Index.find(key);
    if (Index.end() != entryit) {
      value = entryit->second;
    }
    return Index.end() != entryit;
  }

  void add(const TKey& key, TValue found = true) {
    if (key.getNumAtoms() > MaxAtoms) {
      MaxAtoms = key.getNumAtoms();
    }

    Index.insert(std::pair<TKey, bool>(key, found));
  }

  size_t size() const {
    return Index.size();  // for statistics only
  }
};
}  // namespace FMCS
}  // namespace RDKit
