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
#include <list>
#include <vector>
#include <string>
#include <stdexcept>
#include "../RDKitBase.h"
#include "Graph.h"
#include "Seed.h"
#include "DebugTrace.h"  // algorithm filter definitions

namespace RDKit {
namespace FMCS {
class RDKIT_FMCS_EXPORT SubstructureCache {
 public:
#pragma pack(push, 1)
  struct KeyNumericMetrics {
    typedef unsigned long long TValue;
    TValue Value{0};

   public:
    KeyNumericMetrics() {}
  };
#pragma pack(pop)

  struct HashKey {
    KeyNumericMetrics NumericMetrics;

   public:
    void computeKey(const Seed &seed,
                    const std::vector<unsigned int> &queryAtomLabels,
                    const std::vector<unsigned int> &queryBondLabels) {
      computeMorganCodeHash(seed, queryAtomLabels, queryBondLabels);
    }

   private:
    void computeMorganCodeHash(
        const Seed &seed, const std::vector<unsigned int> &queryAtomLabels,
        const std::vector<unsigned int> &queryBondLabels) {
      size_t nv = seed.getNumAtoms();
      size_t ne = seed.getNumBonds();
      std::vector<unsigned long> currCodes(nv);
      std::vector<unsigned long> prevCodes(nv);
      size_t nIterations = seed.getNumBonds();
      if (nIterations > 5) {
        nIterations = 5;
      }

      for (unsigned int seedAtomIdx = 0; seedAtomIdx < seed.getNumAtoms();
           ++seedAtomIdx) {
        currCodes[seedAtomIdx] = queryAtomLabels.at(
            seed.MoleculeFragment.Atoms.at(seedAtomIdx)->getIdx());
      }

      for (size_t iter = 0; iter < nIterations; ++iter) {
        for (size_t i = 0; i < nv; ++i) {
          prevCodes[i] = currCodes[i];
        }

        for (size_t seedBondIdx = 0; seedBondIdx < ne; ++seedBondIdx) {
          const Bond *bond = seed.MoleculeFragment.Bonds[seedBondIdx];
          unsigned int order = queryBondLabels.at(
              seed.MoleculeFragment.Bonds.at(seedBondIdx)->getIdx());
          unsigned int atom1 = seed.MoleculeFragment.SeedAtomIdxMap
                                   .find(bond->getBeginAtomIdx())
                                   ->second;
          unsigned int atom2 =
              seed.MoleculeFragment.SeedAtomIdxMap.find(bond->getEndAtomIdx())
                  ->second;
          unsigned int v1 = prevCodes[atom1];
          unsigned int v2 = prevCodes[atom2];

          currCodes[atom1] += v2 * v2 + (v2 + 23) * (order + 1721);
          currCodes[atom2] += v1 * v1 + (v1 + 23) * (order + 1721);
        }
      }

      KeyNumericMetrics::TValue result = 0;
      for (unsigned int seedAtomIdx = 0; seedAtomIdx < nv; ++seedAtomIdx) {
        unsigned long code = currCodes[seedAtomIdx];
        result += code * (code + 6849) + 29;
      }

      NumericMetrics.Value = result;
    }

    // not implemented yet:
    /*
                void computeFingerprint(const Seed& seed)
                {
                unsigned int radius = seed.getNumBonds();
                if (radius > 5)
                    radius = 5;
                ExplicitBitVect *mf =
       RDKit::MorganFingerprints::getFingerprintAsBitVect(seed.GraphTopology,
       radius);   //SLOW !!!
                // ...
                delete mf;
                    NumericMetrics.Field.hasFingerprint = 1;
                }
    */
  };

  typedef HashKey TKey;
  typedef std::list<FMCS::Graph> TIndexEntry;  // hash-key is not unique key
 private:
  std::vector<TIndexEntry> ValueStorage;
  std::map<KeyNumericMetrics::TValue, size_t> NumericIndex;  // TIndexEntry
 public:
  // returns computed key, and pointer to index entry with a set of subgraphs
  // corresponding to the key if found.
  // then caller must find exactly matched subgraph in the result set with own
  // search algorithm,
  // including a resolving of collisions of hash key
  TIndexEntry *find(const Seed &seed,
                    const std::vector<unsigned int> &queryAtomLabels,
                    const std::vector<unsigned int> &queryBondLabels,
                    TKey &key) {  // compute key and find entry
    key.computeKey(seed, queryAtomLabels, queryBondLabels);
    const auto entryit = NumericIndex.find(key.NumericMetrics.Value);
    if (NumericIndex.end() != entryit) {
      return &ValueStorage[entryit->second];
    }
    return nullptr;  // not found
  }

  // if find() did not found any entry for this key of seed a new entry will be
  // created
  void add(const Seed &seed, TKey &key,
           TIndexEntry *entry) {  // "compute" value and store it in NEW entry
                                  // if not found
    if (!entry) {
      try {
        ValueStorage.emplace_back();
      } catch (...) {
        return;  // not enough memory room to add the item, but it's just a
                 // cache
      }
      entry = &ValueStorage.back();
    }
    entry->push_back(seed.Topology);

    if (!NumericIndex
             .insert(std::make_pair(key.NumericMetrics.Value,
                                    ValueStorage.size() - 1))
             .second) {
      return;  // not enough memory room to add the item, but it is just cache
    }
  }

  size_t keyssize() const {  // for statistics only
    return ValueStorage.size();
  }

  size_t fullsize() const {  // for statistics only
    return std::accumulate(
        ValueStorage.begin(), ValueStorage.end(), 0,
        [](const auto &acc, const auto &v) { return acc + v.size(); });
  }
};
}  // namespace FMCS
}  // namespace RDKit
