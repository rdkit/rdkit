//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once
#include <list>
#include <vector>
#include <string>
#include <stdexcept>
#include "../RDKitBase.h"
//#include "../Fingerprints/MorganFingerprints.h"
#include "Graph.h"
#include "Seed.h"
#include "DebugTrace.h" // algorithm filter definitions

namespace RDKit {
    namespace FMCS {
        class SubstructureCache {
        public:

#pragma pack(push,1)
            struct KeyNumericMetrics {
                typedef unsigned long long TValue;
                TValue Value;
            public:
                KeyNumericMetrics() : Value(0) {}
            };
#pragma pack(pop)

            struct HashKey {
                KeyNumericMetrics   NumericMetrics;
            public:
                void computeKey (const Seed& seed, const std::vector<unsigned>& queryAtomLabels, const std::vector<unsigned>& queryBondLabels) {
                    computeMorganCodeHash(seed, queryAtomLabels, queryBondLabels);
                }
            private:
                void computeMorganCodeHash(const Seed& seed
                                           , const std::vector<unsigned>& queryAtomLabels, const std::vector<unsigned>& queryBondLabels) {
                    size_t nv = seed.getNumAtoms();
                    size_t ne = seed.getNumBonds();
                    std::vector<unsigned long> currCodes(nv);
                    std::vector<unsigned long> prevCodes(nv);
                    size_t nIterations = seed.getNumBonds();
                    if (nIterations > 5)
                        nIterations = 5;

                    for(unsigned seedAtomIdx = 0; seedAtomIdx < seed.getNumAtoms(); seedAtomIdx++)
                        currCodes[seedAtomIdx] = queryAtomLabels[seed.MoleculeFragment.AtomsIdx[seedAtomIdx]];

                    for (size_t iter = 0; iter < nIterations; iter++) {
                        for (size_t i = 0; i < nv; i++)
                            prevCodes[i] = currCodes[i];

                        for (size_t seedBondIdx= 0; seedBondIdx< ne; seedBondIdx++) {
                            const Bond* bond = seed.MoleculeFragment.Bonds[seedBondIdx];
                            unsigned order =  queryBondLabels[seed.MoleculeFragment.BondsIdx[seedBondIdx]];
                            unsigned atom1= seed.MoleculeFragment.SeedAtomIdxMap.find(bond->getBeginAtomIdx())->second;
                            unsigned atom2= seed.MoleculeFragment.SeedAtomIdxMap.find(bond->getEndAtomIdx  ())->second;
                            unsigned v1 = prevCodes[atom1];
                            unsigned v2 = prevCodes[atom2];

                            currCodes[atom1] += v2*v2 + (v2 + 23) * (order + 1721);
                            currCodes[atom2] += v1*v1 + (v1 + 23) * (order + 1721);
                        }
                    }

                    KeyNumericMetrics::TValue result = 0;
                    for(unsigned seedAtomIdx = 0; seedAtomIdx < nv; seedAtomIdx++) {
                        unsigned long code = currCodes[seedAtomIdx];
                        result += code * (code + 6849) + 29;
                    }

                    NumericMetrics.Value = result;

                }

//not implemented yet:
                /*
                            void computeFingerprint(const Seed& seed)
                            {
                            unsigned int radius = seed.getNumBonds();
                            if (radius > 5)
                                radius = 5;
                            ExplicitBitVect *mf = RDKit::MorganFingerprints::getFingerprintAsBitVect(seed.GraphTopology, radius);   //SLOW !!!
                            // ...
                            delete mf;
                                NumericMetrics.Field.hasFingerprint = 1;
                            }
                */
            };

            typedef HashKey                TKey;
            typedef std::list<FMCS::Graph> TIndexEntry; // hash-key is not unique key
        private:
            std::vector<TIndexEntry>                    ValueStorage;
            std::map<KeyNumericMetrics::TValue, size_t> NumericIndex;  // TIndexEntry
        public:
            // returns computed key, and pointer to index entry with a set of subgraphs corresponding to the key if found.
            // then caller must find exactly matched subgraph in the result set with own search algorithm,
            // including a resolving of collisions of hash key
            TIndexEntry* find(const Seed& seed, const std::vector<unsigned>& queryAtomLabels
                              , const std::vector<unsigned>& queryBondLabels
                              , TKey& key) { // compute key and find entry
                key.computeKey(seed, queryAtomLabels, queryBondLabels);
                std::map<KeyNumericMetrics::TValue, size_t>::const_iterator
                entryit = NumericIndex.find(key.NumericMetrics.Value);
                if(NumericIndex.end() != entryit)
                    return & ValueStorage[entryit->second];
                return NULL;    // not found
            }

            //if find() did not found any entry for this key of seed a new entry will be created
            void add(const Seed& seed, TKey& key, TIndexEntry* entry) { // "compute" value and store it in NEW entry if not found
                if(!entry) {
                    try {
                        ValueStorage.push_back(TIndexEntry());
                    } catch(...) {
                        return;   // not enought memory room to add the item, but it's just a cache
                    }
                    entry = &ValueStorage.back();
                }
                entry->push_back(seed.Topology);

                if(!NumericIndex.insert( std::pair<KeyNumericMetrics::TValue, size_t>(key.NumericMetrics.Value, ValueStorage.size()-1) ).second)
                    return; // not enought memory room to add the item, but it is just cache
            }


            size_t keyssize() const { // for statistics only
                return ValueStorage.size();
            }

            size_t fullsize() const { // for statistics only
                size_t n = 0;
                for(std::vector<TIndexEntry>::const_iterator e = ValueStorage.begin(); e != ValueStorage.end(); e++)
                    n += e->size();
                return n;
            }
        };
    }
}
