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
#include <vector>
#include "FMCS.h"
#include "MatchTable.h"

namespace RDKit {
    namespace FMCS {
        struct TargetMatch {
            bool Empty;
            size_t  MatchedAtomSize;
            size_t  MatchedBondSize;
            std::vector<unsigned> TargetAtomIdx;
            std::vector<unsigned> TargetBondIdx;
            std::vector<bool> VisitedTargetBonds;
            std::vector<bool> VisitedTargetAtoms;   // for checking rings
        public:
            TargetMatch() : Empty(true), MatchedAtomSize(0), MatchedBondSize(0) {}
            TargetMatch(const TargetMatch& src) {
                *this = src;
            }
            TargetMatch& operator = (const TargetMatch& src) {
                Empty = src.Empty;
                if(!Empty) {
                    MatchedAtomSize = src.MatchedAtomSize;
                    MatchedBondSize = src.MatchedBondSize;
                    TargetAtomIdx.resize(src.TargetAtomIdx.size());
                    memcpy(&TargetAtomIdx[0], &src.TargetAtomIdx[0], sizeof(unsigned)*TargetAtomIdx.size());
                    TargetBondIdx.resize(src.TargetBondIdx.size());
                    memcpy(&TargetBondIdx[0], &src.TargetBondIdx[0], sizeof(unsigned)*TargetBondIdx.size());
                    VisitedTargetBonds = src.VisitedTargetBonds; // std::vector<bool> bitset
                    VisitedTargetAtoms = src.VisitedTargetAtoms; // std::vector<bool> bitset
                }
                return *this;
            }
            bool empty()const {
                return Empty;
            }
            void clear() {
                Empty = true;

                TargetAtomIdx.clear();
                TargetBondIdx.clear();
                VisitedTargetBonds.clear();
                VisitedTargetAtoms.clear();
            }
            void init(const Seed& seed, const match_V_t& match, const ROMol& query, const Target& target) {
                TargetAtomIdx.resize(query.getNumAtoms());
                TargetBondIdx.resize(query.getNumBonds());
                VisitedTargetBonds.resize(target.Molecule->getNumBonds());
                VisitedTargetAtoms.resize(target.Molecule->getNumAtoms());

                memset(&TargetAtomIdx[0], 0xFF, sizeof(unsigned)*TargetAtomIdx.size());
                memset(&TargetBondIdx[0], 0xFF, sizeof(unsigned)*TargetBondIdx.size());
                /*
                            for(size_t i = 0; i < TargetAtomIdx.size(); i++)
                                TargetAtomIdx[i] = -1;
                            for(size_t i = 0; i < TargetBondIdx.size(); i++)
                                TargetBondIdx[i] = -1;
                */
                for(size_t i = 0; i < VisitedTargetBonds.size(); i++)
                    VisitedTargetBonds[i] = false;
                for(size_t i = 0; i < VisitedTargetAtoms.size(); i++)
                    VisitedTargetAtoms[i] = false;

                MatchedAtomSize = match.size();
                for(match_V_t::const_iterator mit = match.begin(); mit != match.end(); mit++) {
                    TargetAtomIdx[seed.MoleculeFragment.AtomsIdx[mit->first]] = mit->second;
                    VisitedTargetAtoms[mit->second] = true;
                }

                MatchedBondSize = 0;
                for(std::vector<const Bond*>::const_iterator bond = seed.MoleculeFragment.Bonds.begin(); bond != seed.MoleculeFragment.Bonds.end(); bond++) {
                    unsigned i = (*bond)->getBeginAtomIdx();
                    unsigned j = (*bond)->getEndAtomIdx();
                    unsigned ti= TargetAtomIdx[i];
                    unsigned tj= TargetAtomIdx[j];
                    const Bond* tb = target.Molecule->getBondBetweenAtoms(ti, tj);
                    if(tb) {
                        MatchedBondSize++;
                        TargetBondIdx[(*bond)->getIdx()] = tb->getIdx(); // add
                        VisitedTargetBonds[tb->getIdx()] = true;
                    }
                }
                Empty = false;
            }
        };
    }
}
