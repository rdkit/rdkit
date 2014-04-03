#pragma once
#include <map>
#include "../RDKitBase.h"
#include "Graph.h"
#include "DuplicatedSeedCache.h"
#include "DebugTrace.h" // algorithm optimisation definitions

namespace RDKit
{
 namespace FMCS
 {
    class MaximumCommonSubgraph;

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
        unsigned            LastAddedBondsBeginIdx; // in this subgraph for improving performance of future growing
        unsigned            RemainingBonds;
        unsigned            RemainingAtoms;
#ifdef DUP_SUBSTRUCT_CACHE
        DuplicatedSeedCache::TKey DupCacheKey;
#endif
#ifdef FAST_INCREMENTAL_MATCH
        std::vector<BondMatchSet> MatchResult;  // for each target
#endif // FAST_INCREMENTAL_MATCH
    public:
        Seed() : GrowingStage(0), LastAddedAtomsBeginIdx(0), RemainingBonds(-1), RemainingAtoms(-1) {}

        void createFromParent(const Seed* parent)
        {
            *this = *parent; 
            LastAddedAtomsBeginIdx = getNumAtoms();   // previous size
            LastAddedBondsBeginIdx = getNumBonds();   // previous size
            GrowingStage = 0;
        }

        unsigned getNumAtoms()const {return MoleculeFragment.AtomsIdx.size();}
        unsigned getNumBonds()const {return MoleculeFragment.BondsIdx.size();}

        void grow(MaximumCommonSubgraph& mcs, const ROMol& qmol)const;
        bool canGrowBiggerThan(unsigned maxBonds, unsigned maxAtoms)const  // prune()
        {
//TMP DEBUG
/*
if(0==MoleculeFragment.BondsIdx[0])
{
    std::cout<<"\n"
        <<(RemainingBonds + getNumBonds() > maxBonds || RemainingAtoms + getNumAtoms() >= maxAtoms ? "true":"FALSE")
        <<" ------------------------"<<(void*)this 
        <<": LastAddedAtomsBeginIdx = "<<LastAddedAtomsBeginIdx<<", LastAddedBondsBeginIdx = "<<LastAddedBondsBeginIdx<<"\n";
    for(size_t i = 0; i < MoleculeFragment.Bonds.size(); i++)
        std::cout << i << " "<<MoleculeFragment.Bonds[i]->getIdx()<<" : "
                        <<" "<<MoleculeFragment.Bonds[i]->getBeginAtom()->getIdx()
                        <<" "<<MoleculeFragment.Bonds[i]->getEndAtom()->getIdx()
                        <<"\n";
}
*/
//TMP DEBUG
/*
if(MoleculeFragment.BondsIdx.size() >=8
    && MoleculeFragment.BondsIdx[ 0]==0
    && MoleculeFragment.BondsIdx[ 1]==15
    && MoleculeFragment.BondsIdx[ 2]==1
    && MoleculeFragment.BondsIdx[ 3]==16
    && MoleculeFragment.BondsIdx[ 4]==2
    && MoleculeFragment.BondsIdx[ 5]==17
    && MoleculeFragment.BondsIdx[ 6]==13
    && MoleculeFragment.BondsIdx[ 7]==14
    )
;//return true;
*/
            return RemainingBonds + getNumBonds() >  maxBonds
                || RemainingAtoms + getNumAtoms() >= maxAtoms;
            
        }
        void computeRemainingSize(const ROMol& qmol);//, const std::vector<char>& excludedBonds);

        unsigned addAtom(const Atom* atom);
        unsigned addBond(const Bond* bond);
    };

}}
