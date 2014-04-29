// $Id$
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
#include <map>
#include "../RDKitBase.h"
#include "Graph.h"
#include "DuplicatedSeedCache.h"
#include "SubstructMatchCustom.h"
#include "DebugTrace.h" // algorithm optimisation definitions

namespace RDKit
{
 namespace FMCS
 {
    class  MaximumCommonSubgraph;
    struct TargetMatch;

    struct MolFragment  // Reference to a fragment of source molecule
    {
        std::vector<const Atom*>    Atoms;
        std::vector<const Bond*>    Bonds;
        std::vector<unsigned>    AtomsIdx;
        std::vector<unsigned>    BondsIdx;  // need for results and size() only !
        std::map<unsigned, unsigned> SeedAtomIdxMap;   // Full Query Molecule to Seed indeces backward conversion map
    };

    struct Seed
    {
        mutable unsigned    GrowingStage; // 0 new seed; -1 finished; n>0 in progress, exact stage of growing for SDF
    public:
        MolFragment         MoleculeFragment;   // Reference to a fragment of source molecule
        Graph               Topology;           // seed topology with references to source molecule

        std::vector<bool>   ExcludedBonds;
        unsigned            LastAddedAtomsBeginIdx; // in this subgraph for improving performance of future growing
        unsigned            LastAddedBondsBeginIdx; // in this subgraph for DEBUG ONLY
        unsigned            RemainingBonds;
        unsigned            RemainingAtoms;
#ifdef DUP_SUBSTRUCT_CACHE
        DuplicatedSeedCache::TKey DupCacheKey;
#endif
#ifdef FAST_INCREMENTAL_MATCH
        std::vector<TargetMatch> MatchResult;  // for each target
#endif // FAST_INCREMENTAL_MATCH
    public:
        Seed() : GrowingStage(0), LastAddedAtomsBeginIdx(0), LastAddedBondsBeginIdx(0)
               , RemainingBonds(-1), RemainingAtoms(-1) {}

        void createFromParent(const Seed* parent)
        {
//            *this = *parent; 
            MoleculeFragment = parent->MoleculeFragment;
            Topology         = parent->Topology;
            ExcludedBonds    = parent->ExcludedBonds;
            RemainingBonds   = parent->RemainingBonds;
            RemainingAtoms   = parent->RemainingAtoms;
#ifdef DUP_SUBSTRUCT_CACHE
            DupCacheKey      = parent->DupCacheKey;
#endif
            LastAddedAtomsBeginIdx = getNumAtoms();   // previous size
            LastAddedBondsBeginIdx = getNumBonds();   // previous size
            GrowingStage = 0;
        }

        unsigned getNumAtoms()const {return MoleculeFragment.AtomsIdx.size();}
        unsigned getNumBonds()const {return MoleculeFragment.BondsIdx.size();}

        void grow(MaximumCommonSubgraph& mcs, const ROMol& qmol)const;
        bool canGrowBiggerThan(unsigned maxBonds, unsigned maxAtoms)const  // prune()
        {
            return RemainingBonds + getNumBonds() >  maxBonds
                || RemainingAtoms + getNumAtoms() >= maxAtoms;
            
        }
        void computeRemainingSize(const ROMol& qmol);//, const std::vector<char>& excludedBonds);

        unsigned addAtom(const Atom* atom);
        unsigned addBond(const Bond* bond);
    };

}}
