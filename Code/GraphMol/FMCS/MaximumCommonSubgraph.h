#pragma once
#include <vector>
#include <string>
#include <stdexcept>
#include "../RDKitBase.h"
#include "FMCS.h"
#include "SeedSet.h"
#include "Target.h"
#include "SubstructureCache.h"
#include "DuplicatedSeedCache.h"
#include "MatchTable.h"
#include "RingMatchTableSet.h"
#include "DebugTrace.h" // algorithm filter definitions

namespace RDKit
{
 namespace FMCS
 {
    class MaximumCommonSubgraph
    {
        struct MCS  // current result. Reference to a fragment of source molecule
        {
            std::vector<const Atom*> Atoms;
            std::vector<const Bond*> Bonds;
            std::vector<unsigned>    AtomsIdx;
            std::vector<unsigned>    BondsIdx;  // need for results and size() only !
            const ROMol*             QueryMolecule;
            std::vector<Target>      Targets;
        };
        time_t          To;
        MCSProgressData Stat;
        MCSParameters   Parameters;
        unsigned        ThresholdCount;   // min number of matching
        std::vector<const ROMol*> Molecules;
#ifdef FAST_SUBSTRUCT_CACHE
        std::vector<unsigned> QueryAtomLabels;  // for code Morgan. Value based on current functor and parameters
        std::vector<unsigned> QueryBondLabels;  // for code Morgan. Value based on current functor and parameters
        SubstructureCache   HashCache;
#ifdef PRECOMPUTED_TABLES_MATCH
        MatchTable      QueryAtomMatchTable;
        MatchTable      QueryBondMatchTable;
#endif
        RingMatchTableSet   RingMatchTables;
#endif
#ifdef DUP_SUBSTRUCT_CACHE
        DuplicatedSeedCache DuplicateCache;
#endif
        const ROMol*        QueryMolecule;
        std::vector<Target> Targets;
        SeedSet             Seeds;
        MCS                 McsIdx;
    public:
#ifdef VERBOSE_STATISTICS_ON
        ExecStatistics VerboseStatistics;
#endif

        MaximumCommonSubgraph(const MCSParameters* params);
        ~MaximumCommonSubgraph() {clear();}
        MCSResult find (const std::vector<ROMOL_SPTR>& mols);
        unsigned getMaxNumberBonds()const {return McsIdx.BondsIdx.size();}  //Seeds.MaxBonds;}
        unsigned getMaxNumberAtoms()const {return McsIdx.AtomsIdx.size();}  //Seeds.MaxAtoms;}
    //internal:
        bool checkIfMatchAndAppend(Seed& seed);//, const std::vector<char>& excludedBonds);   // REPLACE with Swap !!!
    private:
        void clear()
        {
            Targets.clear();
            Molecules.clear();
        }
        void init();
        void makeInitialSeeds();
        bool createSeedFromMCS(size_t newQueryTarget, Seed& seed);
        bool growSeeds();   //returns false if canceled
        static std::string generateResultSMARTS(const MCS& McsIdx);

        bool match(Seed& seed);
        bool matchIncrementalFast(Seed& seed, unsigned itag); // use results of previous match stored in the seed
    };

}} // namespace RDKit
