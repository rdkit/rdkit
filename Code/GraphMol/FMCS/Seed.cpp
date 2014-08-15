//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MaximumCommonSubgraph.h"
#include "Composition2N.h"
#include "Seed.h"

#include "DebugTrace.h"
#include "../SmilesParse/SmilesWrite.h"

namespace RDKit {
    namespace FMCS {

        unsigned Seed::addAtom(const Atom* atom) {
            unsigned i = MoleculeFragment.AtomsIdx.size();
            unsigned aqi = atom->getIdx();
            MoleculeFragment.Atoms.push_back(atom);
            MoleculeFragment.AtomsIdx.push_back(aqi);
            MoleculeFragment.SeedAtomIdxMap[aqi] = i;
            Topology.addAtom(aqi);
#ifdef DUP_SUBSTRUCT_CACHE
            DupCacheKey.addAtom(aqi);
#endif
            return i;
        }

        unsigned Seed::addBond(const Bond* bond) {
            unsigned b = bond->getIdx();
            if(ExcludedBonds[b])
                throw -1;   //never, check the implementation
            ExcludedBonds[b] = true;
            MoleculeFragment.BondsIdx.push_back(b);
            MoleculeFragment.Bonds.push_back(bond);
            // remap idx to seed's indeces:
            unsigned i = MoleculeFragment.SeedAtomIdxMap[bond->getBeginAtomIdx()];
            unsigned j = MoleculeFragment.SeedAtomIdxMap[bond->getEndAtomIdx()];
            Topology.addBond(b, i, j);
#ifdef DUP_SUBSTRUCT_CACHE
            DupCacheKey.addBond(b);
#endif
            return getNumBonds();
        }

        void Seed::fillNewBonds(const ROMol& qmol) {
            std::vector<bool> excludedBonds = ExcludedBonds;
            for(unsigned srcAtomIdx = LastAddedAtomsBeginIdx; srcAtomIdx < getNumAtoms(); srcAtomIdx++) { // all atoms added on previous growing only
                const Atom* atom = MoleculeFragment.Atoms[srcAtomIdx];
                ROMol::OEDGE_ITER beg,end;
                for(boost::tie(beg,end) = qmol.getAtomBonds(atom); beg!=end; beg++) { // all bonds from MoleculeFragment.Atoms[srcAtomIdx]
                    const Bond* bond = &*(qmol[*beg]);
                    if( ! excludedBonds[bond->getIdx()]) { // already in the seed or in the NewBonds list from another atom in a RING
                        excludedBonds[bond->getIdx()] = true;
                        unsigned ai = (atom == bond->getBeginAtom()) ? bond->getEndAtomIdx() : bond->getBeginAtomIdx();
                        const Atom* end_atom = qmol.getAtomWithIdx(ai);
                        unsigned end_atom_idx = -1;
                        for(unsigned i=0; i < getNumAtoms(); i++)
                            if(end_atom == MoleculeFragment.Atoms[i]) {  // already exists in this seed
                                end_atom_idx = i;
                                break;
                            }
                        NewBonds.push_back(NewBond(srcAtomIdx, bond->getIdx(), ai, end_atom_idx, -1==end_atom_idx ? end_atom:0));
                    }
                }
            }
        }


        void Seed::grow(MaximumCommonSubgraph& mcs) const {
            const ROMol& qmol = mcs.getQueryMolecule();
            std::map<unsigned, unsigned> newAtomsMap; // map new added atoms to their seed's indeces

            if(!canGrowBiggerThan(mcs.getMaxNumberBonds(), mcs.getMaxNumberAtoms()) ) { // prune() parent
                GrowingStage = -1; //finished
#ifdef VERBOSE_STATISTICS_ON
                ++mcs.VerboseStatistics.RemainingSizeRejected;
#endif
                return;
            }

            if(0==GrowingStage) {
                // 0. Fill out list of all directly connected outgoing bonds
                ((Seed*)this)->fillNewBonds(qmol); // non const method, multistage growing optimisation

                if(NewBonds.empty()) {
                    GrowingStage = -1;  // finished
                    return;
                }

                // 1. Check and add the biggest child seed with all outgoing bonds added:
                // Add all bonds at first (build the biggest child seed). All new atoms are already in the seed
                Seed seed;
                seed.createFromParent(this);

                for(std::vector<NewBond>::const_iterator nbi = NewBonds.begin(); nbi != NewBonds.end(); nbi++) {
                    unsigned aIdx = nbi->EndAtomIdx;
                    if(-1 == aIdx) { // new atom
                        std::map<unsigned, unsigned>::const_iterator nai = newAtomsMap.find(nbi->NewAtomIdx);         // check RING
                        if(newAtomsMap.end() == nai) {
                            const Atom* end_atom = nbi->NewAtom;
                            aIdx = seed.addAtom(end_atom);
                            newAtomsMap[nbi->NewAtomIdx] = aIdx;    // store new possible ring end point
                        } else
                            aIdx = nai->second;
                    }
                    const Bond* src_bond = qmol.getBondWithIdx(nbi->BondIdx);
                    seed.addBond(src_bond);
                }
#ifdef VERBOSE_STATISTICS_ON
                {
                    ++mcs.VerboseStatistics.Seed;
                }
#endif
                seed.RemainingBonds = RemainingBonds - NewBonds.size();     // Added ALL !!!
                seed.RemainingAtoms = RemainingAtoms - newAtomsMap.size();  // new atoms added to seed

                // prune() Best Sizes
                if( ! seed.canGrowBiggerThan(mcs.getMaxNumberBonds(), mcs.getMaxNumberAtoms()) ) {
                    GrowingStage = -1;
#ifdef VERBOSE_STATISTICS_ON
                    ++mcs.VerboseStatistics.RemainingSizeRejected;
#endif
                    return; // the biggest possible subrgaph from this seed is too small for future growing. So, skip ALL children !
                }
                seed.MatchResult = MatchResult;
                bool allMatched = mcs.checkIfMatchAndAppend(seed);  // this seed + all extern bonds is a part of MCS


                GrowingStage = 1;
                if(allMatched && NewBonds.size() > 1)
                    return; // grow deep first. postpone next growing steps
            }
// 2. Check and add all 2^N-1-1 other possible seeds:
            if(1 == NewBonds.size()) {
                GrowingStage = -1;
                return; // everything has been done
            }
            // OPTIMISATION:
            // check each single bond first: if (this seed + single bond) does not exist in MCS, exclude this new bond from growing this seed.
            unsigned numErasedNewBonds = 0;
            for(std::vector<NewBond>::iterator nbi = NewBonds.begin(); nbi != NewBonds.end(); nbi++) {
#ifdef VERBOSE_STATISTICS_ON
                {
                    ++mcs.VerboseStatistics.Seed;
                }
#endif
                Seed seed;
                seed.createFromParent(this);
                newAtomsMap.clear();

                unsigned aIdx = nbi->EndAtomIdx;    // existed in this parent seed (ring) or -1
                if(-1 == aIdx) { // new atom
                    const Atom* end_atom = nbi->NewAtom;
                    aIdx = seed.addAtom(end_atom);
                }
                const Bond* src_bond = qmol.getBondWithIdx(nbi->BondIdx);
                seed.addBond(src_bond);
                seed.computeRemainingSize(qmol);

                if(seed.canGrowBiggerThan(mcs.getMaxNumberBonds(), mcs.getMaxNumberAtoms()) ) { // prune()
                    if(!MatchResult.empty())
                        seed.MatchResult = MatchResult;
                    if( ! mcs.checkIfMatchAndAppend(seed)) {
                        nbi->BondIdx = -1; // exclude this new bond from growing this seed - decrease 2^^N-1 to 2^^k-1, k<N.
                        ++numErasedNewBonds;
#ifdef VERBOSE_STATISTICS_ON
                        ++mcs.VerboseStatistics.SingleBondExcluded;
#endif
                    }
                } else { // seed too small
#ifdef VERBOSE_STATISTICS_ON
                    ++mcs.VerboseStatistics.RemainingSizeRejected;
#endif
                }
            }

            if(numErasedNewBonds > 0) {
                std::vector<NewBond> dirtyNewBonds;
                dirtyNewBonds.reserve(NewBonds.size());
                dirtyNewBonds.swap(NewBonds);
                for(std::vector<NewBond>::const_iterator nbi = dirtyNewBonds.begin(); nbi != dirtyNewBonds.end(); nbi++)
                    if(-1 != nbi->BondIdx)
                        NewBonds.push_back(*nbi);
            }

            // add all other from 2^k-1 possible seeds, where k=newBonds.size():
            if(NewBonds.size() > 1) { // if just one new bond, such seed has been already created
                if(sizeof(unsigned long long)*8 < NewBonds.size())
                    throw std::runtime_error("Max number of new external bonds of a seed more than 64");
                BitSet maxCompositionValue;
                Composition2N::compute2N(NewBonds.size(), maxCompositionValue);
                maxCompositionValue -= 1;   // 2^N-1
                Composition2N composition(maxCompositionValue, maxCompositionValue);

#ifdef EXCLUDE_WRONG_COMPOSITION
                std::vector<BitSet> failedCombinations;
                BitSet              failedCombinationsMask=0uLL;
#endif
                while(composition.generateNext()) {
                    if(composition.is2Power()) // exclude already processed single external bond combinations
                        continue;
                    if(0==numErasedNewBonds && composition.getBitSet() == maxCompositionValue)
                        continue;   // exclude already processed all external bonds combination 2N-1
#ifdef EXCLUDE_WRONG_COMPOSITION
// OPTIMISATION. reduce amount of generated seeds and match calls
// 2120 instead of 2208 match calls on small test. 43 wrongComp-s, 83 rejected
                    if(failedCombinationsMask & composition.getBitSet()) { // possibly exists in the list
                        bool compositionWrong = false;
                        for(std::vector<BitSet>::const_iterator failed = failedCombinations.begin();
                                failed != failedCombinations.end() ; failed++)
                            if(*failed == (*failed & composition.getBitSet())) { // combination includes failed combination
                                compositionWrong = true;
                                break;
                            }
                        if(compositionWrong) {
#ifdef VERBOSE_STATISTICS_ON
                            ++mcs.VerboseStatistics.WrongCompositionRejected;
#endif
                            continue;
                        }
                    }
#endif
#ifdef VERBOSE_STATISTICS_ON
                    {
                        ++mcs.VerboseStatistics.Seed;
                    }
#endif
                    Seed seed;
                    seed.createFromParent(this);
                    newAtomsMap.clear();

                    for(unsigned i=0; i<NewBonds.size(); i++)
                        if(composition.isSet(i)) {
                            const NewBond* nbi = & NewBonds[i];
                            unsigned aIdx = nbi->EndAtomIdx;    // existed in this parent seed (ring) or -1
                            if(-1 == aIdx) { // new atom
                                std::map<unsigned, unsigned>::const_iterator nai = newAtomsMap.find(nbi->NewAtomIdx);         // check RING
                                if(newAtomsMap.end() == nai) {
                                    const Atom* end_atom = nbi->NewAtom;//qmol.getAtomWithIdx(nbi->NewAtomIdx);
                                    aIdx = seed.addAtom(end_atom);
                                    newAtomsMap[nbi->NewAtomIdx] = aIdx;    // store new possible ring end point
                                } else
                                    aIdx = nai->second;
                            }

                            const Bond* src_bond = qmol.getBondWithIdx(nbi->BondIdx);
                            seed.addBond(src_bond);
                        }
                    seed.computeRemainingSize(qmol);
                    if( ! seed.canGrowBiggerThan(mcs.getMaxNumberBonds(), mcs.getMaxNumberAtoms()) ) { // prune(). // seed too small
#ifdef VERBOSE_STATISTICS_ON
                        ++mcs.VerboseStatistics.RemainingSizeRejected;
#endif
                    } else {
                        seed.MatchResult = MatchResult;
                        bool found = mcs.checkIfMatchAndAppend(seed);

                        if(!found) {
#ifdef EXCLUDE_WRONG_COMPOSITION  // if seed does not matched it is possible to exclude this mismatched combination for performance improvement
                            failedCombinations.push_back(composition.getBitSet());
                            failedCombinationsMask &= composition.getBitSet();
#ifdef VERBOSE_STATISTICS_ON
                            ++mcs.VerboseStatistics.WrongCompositionDetected;
#endif
#endif
                        }
                    }
                }
            }
            GrowingStage = -1; //finished
        }

        void Seed::computeRemainingSize(const ROMol& qmol) {
            RemainingBonds = RemainingAtoms = 0;

            std::vector<unsigned> end_atom_stack;
            std::vector<bool>     visitedBonds = ExcludedBonds;
            std::vector<bool>     visitedAtoms(qmol.getNumAtoms());

            for(size_t i = 0; i < visitedAtoms.size(); i++)
                visitedAtoms[i] = false;
            for(std::vector<unsigned>::const_iterator it = MoleculeFragment.AtomsIdx.begin(); it != MoleculeFragment.AtomsIdx.end(); it++)
                visitedAtoms[*it] = true;

            // SDF all paths
            // 1. direct neighbours
            for(unsigned seedAtomIdx = LastAddedAtomsBeginIdx; seedAtomIdx < getNumAtoms(); seedAtomIdx++) { // just now added new border vertices (candidates for future growing)
                const Atom* atom = MoleculeFragment.Atoms[seedAtomIdx];
                ROMol::OEDGE_ITER beg,end;
                for(boost::tie(beg,end) = qmol.getAtomBonds(atom); beg!=end; beg++) { // all bonds from MoleculeFragment.Atoms[srcAtomIdx]
                    const Bond& bond = *(qmol[*beg]);
                    if( ! visitedBonds[bond.getIdx()]) {
                        ++RemainingBonds;
                        visitedBonds[bond.getIdx()] = true;
                        unsigned end_atom_idx = (MoleculeFragment.AtomsIdx[seedAtomIdx] == bond.getBeginAtomIdx()) ? bond.getEndAtomIdx() : bond.getBeginAtomIdx();
                        if( ! visitedAtoms[end_atom_idx]) { // check RING/CYCLE
                            ++RemainingAtoms;
                            visitedAtoms[end_atom_idx] = true;
                            end_atom_stack.push_back(end_atom_idx);
                        }
                    }
                }
            }
            // 2. go deep
            while(!end_atom_stack.empty()) {
                unsigned ai = end_atom_stack.back();
                end_atom_stack.pop_back();
                const Atom* atom = qmol.getAtomWithIdx(ai);
                ROMol::OEDGE_ITER beg,end;
                for(boost::tie(beg,end) = qmol.getAtomBonds(atom); beg!=end; beg++) { // all bonds from end_atom
                    const Bond& bond = *(qmol[*beg]);
                    if( ! visitedBonds[bond.getIdx()]) {
                        ++RemainingBonds;
                        visitedBonds[bond.getIdx()] = true;
                        unsigned end_atom_idx = (ai == bond.getBeginAtomIdx()) ? bond.getEndAtomIdx() : bond.getBeginAtomIdx();
                        if( ! visitedAtoms[end_atom_idx]) { // check RING/CYCLE
                            ++RemainingAtoms;
                            visitedAtoms[end_atom_idx] = true;
                            end_atom_stack.push_back(end_atom_idx);
                        }
                    }
                }
            }
        }

    }
}
