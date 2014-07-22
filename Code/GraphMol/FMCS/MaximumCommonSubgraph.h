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
#include <string>
#include <stdexcept>
#include "../RDKitBase.h"
#include "FMCS.h"
#include "DebugTrace.h" // algorithm filter definitions
#include "SeedSet.h"
#include "Target.h"
#include "SubstructureCache.h"
#include "DuplicatedSeedCache.h"
#include "MatchTable.h"
#ifdef FAST_INCREMENTAL_MATCH
#include "TargetMatch.h"
#endif
#ifdef MULTI_THREAD
#include "ThreadPool.h"
#endif
#include "RingMatchTableSet.h"


namespace RDKit {
    namespace FMCS {
        class MaximumCommonSubgraph
#ifdef MULTI_THREAD
            : private Mutex
#endif
        {
            struct MCS { // current result. Reference to a fragment of source molecule
                std::vector<const Atom*> Atoms;
                std::vector<const Bond*> Bonds;
                std::vector<unsigned>    AtomsIdx;
                std::vector<unsigned>    BondsIdx;  // need for results and size() only !
                const ROMol*             QueryMolecule;
                std::vector<Target>      Targets;
            };

            unsigned long long To;
            MCSProgressData Stat;
            MCSParameters   Parameters;
            unsigned        ThresholdCount;   // min number of matching
            std::vector<const ROMol*> Molecules;
#ifdef FAST_SUBSTRUCT_CACHE
            std::vector<unsigned> QueryAtomLabels;  // for code Morgan. Value based on current functor and parameters
            std::vector<unsigned> QueryBondLabels;  // for code Morgan. Value based on current functor and parameters
            SubstructureCache     HashCache;
#ifdef PRECOMPUTED_TABLES_MATCH
            MatchTable          QueryAtomMatchTable;
            MatchTable          QueryBondMatchTable;
#endif
            RingMatchTableSet   RingMatchTables;
#endif
#ifdef DUP_SUBSTRUCT_CACHE
            DuplicatedSeedCache DuplicateCache;
#endif
            const ROMol*        QueryMolecule;
            unsigned            QueryMoleculeMatchedBonds;
            unsigned            QueryMoleculeMatchedAtoms;
            std::vector<Target> Targets;
            SeedSet             Seeds;
            MCS                 McsIdx;
        public:
#ifdef VERBOSE_STATISTICS_ON
            ExecStatistics VerboseStatistics;
#ifdef MULTI_THREAD
            Mutex   StatisticsMutex;
#endif
#endif

            MaximumCommonSubgraph(const MCSParameters* params);
            ~MaximumCommonSubgraph() {
                clear();
            }
            MCSResult find (const std::vector<ROMOL_SPTR>& mols);
            const ROMol& getQueryMolecule()const {
                return *QueryMolecule;
            }
            /* optional
               unsigned getMaxNumberBonds()
               #ifdef MULTI_THREAD
               {
               Guard lock(*this);
               #else
               const
               {
               #endif
            */        unsigned getMaxNumberBonds() const {
                return McsIdx.BondsIdx.size();
            }

            /* optional
               unsigned getMaxNumberAtoms()
               #ifdef MULTI_THREAD
               {
               Guard lock(*this);
               #else
               const
               {
               #endif
            */        unsigned getMaxNumberAtoms() const {
                return McsIdx.AtomsIdx.size();
            }
            //internal:
            bool checkIfMatchAndAppend(Seed& seed);//, const std::vector<char>& excludedBonds);   // REPLACE with Swap !!!
        private:
            void clear() {
                Targets.clear();
                Molecules.clear();
                To = nanoClock();
            }
            void init();
            void makeInitialSeeds();
            bool createSeedFromMCS(size_t newQueryTarget, Seed& seed);
            bool growSeeds();   //returns false if canceled
            std::string generateResultSMARTS(const MCS& McsIdx)const;

            bool match(Seed& seed);
            bool matchIncrementalFast(Seed& seed, unsigned itarget);
        };
    }
} // namespace RDKit
