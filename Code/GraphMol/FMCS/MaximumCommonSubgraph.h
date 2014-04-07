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

#ifdef SMILES_CACHE 
        CanonicalSMARTS_Set MatchedCanonicalSMARTS;
#ifdef SMILES_CACHE_NOT_MATCHED //opt. OPTIMIZATION ???:
        CanonicalSMARTS_Set DoNotMatchedCanonicalSMARTS;
#endif
#endif
        const ROMol*        QueryMolecule;
        std::vector<Target> Targets;
//    public:
        SeedSet Seeds;

    public:
        MaximumCommonSubgraph(const MCSParameters* params);
        ~MaximumCommonSubgraph() {clear();}
        MCSResult find (const std::vector<ROMOL_SPTR>& mols);
        unsigned getMaxNumberBonds()const {return Seeds.MaxBonds;}
        unsigned getMaxNumberAtoms()const {return Seeds.MaxAtoms;}
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
        void growSeeds(MolFragment& mcsIdx, MCSResult& res);
        std::string generateSMARTS(const MolFragment& mcsIdx);

        bool match(Seed& seed);
        bool matchIncrementalFast(Seed& seed, unsigned itag); // use results of previous match stored in the seed
    };

}} // namespace RDKit
