#include "MaximumCommonSubgraph.h"
#include "Seed.h"

#include "DebugTrace.h"
#include "../SmilesParse/SmilesWrite.h"

namespace RDKit
{
namespace FMCS
{
    typedef unsigned long long BitSet;
//        template <class BitSet = unsigned long long>
    class Composition2N // generator of 2^^N-1 possible bit combinations
    {
        BitSet  Bits;
        BitSet  MaxValue;
    public:
        Composition2N(const BitSet& maxValue) : Bits(0), MaxValue(maxValue) {}

        static void compute2N(unsigned numBits, BitSet& maxValue)
        {
            maxValue = 0;
            while(0!=numBits)
                maxValue |= (1uLL << (--numBits));
        }
        BitSet getBitSet()const {return Bits;}
        bool generateNext()
        {
            return (++Bits) <= MaxValue;
        }
        bool is2Power() // one bit is set only
        {
            BitSet  bits = Bits;
            unsigned n = 0;
            while(0==(bits & 1uLL) && ++n < sizeof(bits)*8)    //find lowest bitwise 1
                bits >>= 1uLL;      //shift all zero lower bits
            if(0!=(bits & 1uLL))
                bits >>= 1uLL;      //shift first set bit
            return 0==bits;  //remained bits except lowest 1
        }
        bool nonZero() {return 0!=Bits;}
        bool isSet(unsigned bit) { return 0 != (Bits & (1uLL << bit));}
    };


unsigned Seed::addAtom(const Atom* atom)
{
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

unsigned Seed::addBond(const Bond* bond)
{
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


struct NewBond
{
    unsigned    SourceAtomIdx;  // index in the seed. Atom is already in the seed
    unsigned    BondIdx;        // index in qmol of new bond scheduled to be added into seed. This is outgoing bond from SourceAtomIdx
    unsigned    NewAtomIdx;     // index in qmol of new atom scheduled to be added into seed. Another end of new bond
    const Atom* NewAtom;        // pointer to qmol's new atom scheduled to be added into seed. Another end of new bond
    unsigned    EndAtomIdx;     // index in the seed. RING. "New" Atom on the another end of new bond is already exists in the seed.

    NewBond() : SourceAtomIdx(-1), BondIdx(-1), NewAtomIdx(-1), NewAtom(0), EndAtomIdx(-1) {}

    NewBond(unsigned from_atom, unsigned bond_idx, unsigned new_atom, unsigned to_atom, const Atom* a)
            : SourceAtomIdx(from_atom), BondIdx(bond_idx), NewAtomIdx(new_atom), NewAtom(a), EndAtomIdx(to_atom) {}
};


void Seed::grow(MaximumCommonSubgraph& mcs, const ROMol& qmol) const
{
    std::vector<bool>            excludedBonds = ExcludedBonds;
    std::vector<NewBond>         newBonds;    // all directly connected outgoing bonds
    std::map<unsigned, unsigned> newAtomsMap; // map new added atoms to their seed's indeces

//TMP DEBUG

//DEBUG BREAKPOINT test 504
// test 504
// MCS atoms=(1, 2, 0, 3, 17, 4, 18, 5, 16, 19, 6, 7, 20, 34, 8, 21, 9, 15, 22, 10, 14, 23, 11, 12, 24, 13, 25, 33, 26, 27, 32, 28, 31, 29, ) new atoms=[29, ]
//     bonds=(0, 15, 1, 16, 2, 17, 13, 14, 18, 3, 34, 19, 36, 4, 5, 20, 6, 35, 21, 7, 12, 32, 8, 9, 11, 22, 10, 30, 31, 23, 24, 37, 25, 29, 26, 28, ) 
// Size=34, 36 Remain=1, 1
// [CH2:34][N:19]([CH2:18][CH2:17][CH2:0][CH2:1][CH2:2][NH:3][C:4](=[O:16])[CH:5]1[CH2:6][CH:7]1[c:8]1[cH:15][cH:14][c:12]([Cl:13])[c:10]([Cl:11])[cH:9]1)[CH2:20][CH2:21][CH:22][NH:23][C:24](=[O:33])[NH:25][c:26]1[cH:32][cH:31][c:29][cH:28][cH:27]1
if(MoleculeFragment.BondsIdx.size() == 11+4 //+3
//common bonds
    && MoleculeFragment.BondsIdx[ 0]==0
    && MoleculeFragment.BondsIdx[ 1]==15
    && MoleculeFragment.BondsIdx[ 2]==1
    && MoleculeFragment.BondsIdx[ 3]==16
    && MoleculeFragment.BondsIdx[ 4]==2
    && MoleculeFragment.BondsIdx[ 5]==17
    && MoleculeFragment.BondsIdx[ 6]==13
    && MoleculeFragment.BondsIdx[ 7]==14
    && MoleculeFragment.BondsIdx[ 8]==18
  && MoleculeFragment.BondsIdx[ 9]==3
  && MoleculeFragment.BondsIdx[10]==34
    && MoleculeFragment.BondsIdx[11]==19
    && MoleculeFragment.BondsIdx[12]==36
    && MoleculeFragment.BondsIdx[13]==4
    && MoleculeFragment.BondsIdx[14]==5
// part of MCS
//  && MoleculeFragment.BondsIdx[15]==20
//  && MoleculeFragment.BondsIdx[16]==6
//  && MoleculeFragment.BondsIdx[17]==35
   )
{
    newAtomsMap.clear();
}

    for(unsigned srcAtomIdx = LastAddedAtomsBeginIdx; srcAtomIdx < getNumAtoms(); srcAtomIdx++)   // all atoms added on previous growing only
    {
        const Atom* atom = MoleculeFragment.Atoms[srcAtomIdx];
        ROMol::OEDGE_ITER beg,end;
        for(boost::tie(beg,end) = qmol.getAtomBonds(atom); beg!=end; beg++)  // all bonds from MoleculeFragment.Atoms[srcAtomIdx] 
        {
            const Bond* bond = &*(qmol[*beg]);
            if( ! excludedBonds[bond->getIdx()])   // already in the seed or in the newBonds list from another and of a RING
            {
                excludedBonds[bond->getIdx()] = true;
                unsigned ai = (atom == bond->getBeginAtom()) ? bond->getEndAtomIdx() : bond->getBeginAtomIdx();
                const Atom* end_atom = qmol.getAtomWithIdx(ai);
                unsigned end_atom_idx = -1;
                for(unsigned i=0; i < getNumAtoms(); i++)
                 if(end_atom == MoleculeFragment.Atoms[i])    // already exists in this seed
                {
                    end_atom_idx = i;
                    break;
                }
                newBonds.push_back(NewBond(srcAtomIdx, bond->getIdx(), ai, end_atom_idx, -1==end_atom_idx ? end_atom:0));
            }
        }
    }
    if(newBonds.empty())
    {
        GrowingStage = -1;  // finished
        return;
    }

if(0==GrowingStage)
{
    // 1. Check and add the biggest child seed with all outgoing bonds added:
    // Add all bonds at first (build the biggest child seed). All new atoms are already in the seed
    Seed seed;
    seed.createFromParent(this);

    for(std::vector<NewBond>::iterator nbi = newBonds.begin(); nbi != newBonds.end(); nbi++)
    {
        unsigned aIdx = nbi->EndAtomIdx;
        if(-1 == aIdx) // new atom
        {
            std::map<unsigned, unsigned>::const_iterator nai = newAtomsMap.find(nbi->NewAtomIdx);         // check RING
            if(newAtomsMap.end() == nai)
            {
                const Atom* end_atom = nbi->NewAtom;//qmol.getAtomWithIdx(nbi->NewAtomIdx);
                aIdx = seed.addAtom(end_atom);
                newAtomsMap[nbi->NewAtomIdx] = aIdx;    // store new possible ring end point
            }
            else
                aIdx = nai->second;
        }
        const Bond* src_bond = qmol.getBondWithIdx(nbi->BondIdx);
        seed.addBond(src_bond);
    }
#ifdef VERBOSE_STATISTICS_ON
    ++extstat.Seed;
#endif
    seed.RemainingBonds = RemainingBonds - newBonds.size();     // Added ALL !!!
    seed.RemainingAtoms = RemainingAtoms - newAtomsMap.size();  // new atoms added to seed
//TMP DEBUG
/*
 //TEST - PASSED
{
    unsigned rb = seed.RemainingBonds;
    unsigned ra = seed.RemainingAtoms;
    seed.computeRemainingSize(qmol);
    if(rb != seed.RemainingBonds || ra != seed.RemainingAtoms)
        printf("*** ERROR: rb != seed.RemainingBonds || ra != seed.RemainingAtoms\n");
}
*/
    // prune() Best Sizes
    if( ! seed.canGrowBiggerThan(mcs.getMaxNumberBonds(), mcs.getMaxNumberAtoms()) )
    {
#ifdef VERBOSE_STATISTICS_ON
        ++extstat.RemainingSizeRejected;
#endif
        GrowingStage = -1;
        return; // the biggest possible subrgaph from this seed is too small for future growing. So, skip ALL children !
    }

    bool allMatched = mcs.checkIfMatchAndAppend(seed);//, excludedBonds);   // this seed + all extern bonds is a part of MCS

    GrowingStage = 1;
    if(allMatched && newBonds.size() > 1)
        return; // grow deep first. postpone next growing steps
}
// 2. Check and add all 2^^N-1 - N other possible seeds:
    if(1 == newBonds.size())
    {
        GrowingStage = -1;
        return; // everything has been done
    }
    // OPTIMISATION:
    // check each single bond first: if (this seed + one bond) does not exist in MCS, exclude this new bond from growing this seed.
    unsigned numErasedNewBonds = 0;
    for(std::vector<NewBond>::iterator nbi = newBonds.begin(); nbi != newBonds.end(); nbi++)
    {
#ifdef VERBOSE_STATISTICS_ON
        ++extstat.Seed;
#endif
        Seed seed;
        seed.createFromParent(this);
        newAtomsMap.clear();

        unsigned aIdx = nbi->EndAtomIdx;    // existed in this parent seed (ring) or -1
        if(-1 == aIdx) // new atom
        {
            const Atom* end_atom = nbi->NewAtom;
            aIdx = seed.addAtom(end_atom);
        }
        const Bond* src_bond = qmol.getBondWithIdx(nbi->BondIdx);
        seed.addBond(src_bond);
        seed.computeRemainingSize(qmol);

        if(seed.canGrowBiggerThan(mcs.getMaxNumberBonds(), mcs.getMaxNumberAtoms()) )   // prune()
        {
            if( ! mcs.checkIfMatchAndAppend(seed))
            {
//DEBUG test 504
/*if(33!=nbi->BondIdx && 27!=nbi->BondIdx
   && seed.MoleculeFragment.BondsIdx.end()==find(seed.MoleculeFragment.BondsIdx.begin(),seed.MoleculeFragment.BondsIdx.end(),33)
   && seed.MoleculeFragment.BondsIdx.end()==find(seed.MoleculeFragment.BondsIdx.begin(),seed.MoleculeFragment.BondsIdx.end(),27))
{
    if( ! mcs.checkIfMatchAndAppend(seed))
        continue; //NEVER
}
*/
                nbi->BondIdx = -1; // exclude this new bond from growing this seed - decrease 2^^N-1 to 2^^k-1, k<N.
                ++numErasedNewBonds;
            }
        }
        else   // seed too small
        {
/*just this s eed is small
if(33!=nbi->BondIdx && 27!=nbi->BondIdx
   && seed.MoleculeFragment.BondsIdx.end()==find(seed.MoleculeFragment.BondsIdx.begin(),seed.MoleculeFragment.BondsIdx.end(),33)
   && seed.MoleculeFragment.BondsIdx.end()==find(seed.MoleculeFragment.BondsIdx.begin(),seed.MoleculeFragment.BondsIdx.end(),27))
{
    seed.computeRemainingSize(qmol);
    std::cout<<"REJECT atoms=(";
    for(unsigned seedAtomIdx = 0; seedAtomIdx < seed.getNumAtoms(); seedAtomIdx++)
    {
        const Atom* atom = seed.MoleculeFragment.Atoms[seedAtomIdx];
        std::cout<<atom->getIdx()<<", ";
    }
    std::cout<<") new atoms=[";
    for(unsigned seedAtomIdx = seed.LastAddedAtomsBeginIdx; seedAtomIdx < seed.getNumAtoms(); seedAtomIdx++)
    {
        const Atom* atom = seed.MoleculeFragment.Atoms[seedAtomIdx];
        std::cout<<atom->getIdx()<<", ";
    }
    std::cout<<"] Size="<< seed.getNumAtoms() <<", "<< seed.getNumBonds() <<" Remain=" << seed.RemainingAtoms <<", "<<seed.RemainingBonds<<" = ";
    std::cout<< MolFragmentToSmiles(qmol, *(const std::vector<int>*) &seed.MoleculeFragment.AtomsIdx, (const std::vector<int>*) &seed.MoleculeFragment.BondsIdx) <<"\n";  // unsigned
//    continue;
}
*/
#ifdef VERBOSE_STATISTICS_ON
            ++extstat.RemainingSizeRejected;
#endif
        }
    }

    if(numErasedNewBonds > 0)
    {
        std::vector<NewBond> dirtyNewBonds;
        dirtyNewBonds.reserve(newBonds.size());
        dirtyNewBonds.swap(newBonds);
        for(std::vector<NewBond>::const_iterator nbi = dirtyNewBonds.begin(); nbi != dirtyNewBonds.end(); nbi++)
            if(-1 != nbi->BondIdx)
                newBonds.push_back(*nbi);
/* // NEW with The SAME Result:
        std::vector<NewBond> clearNewBonds;
        for(std::vector<NewBond>::const_iterator nbi = newBonds.begin(); nbi != newBonds.end(); nbi++)
            if(-1 != nbi->BondIdx)
                clearNewBonds.push_back(*nbi);
        newBonds.clear();
        for(std::vector<NewBond>::const_iterator nbi = clearNewBonds.begin(); nbi != clearNewBonds.end(); nbi++)
            newBonds.push_back(*nbi);
*/
    }

    // add all other from 2^^k-1 possible seeds, where k=newBonds.size():
    if(newBonds.size() > 1) // if just one new bond, such seed has been already created
    {
        if(sizeof(unsigned long long)*8 < newBonds.size())
            throw std::runtime_error("Max number of new external bonds of a seed more than 64");
        BitSet maxCompositionValue;
        Composition2N::compute2N(newBonds.size(), maxCompositionValue);
        maxCompositionValue -= 1;   // 2^N-1
        if(0==numErasedNewBonds)
            maxCompositionValue -= 1;   // exclude already processed all external bonds combination 2N-2
        Composition2N composition(maxCompositionValue);

#ifdef EXCLUDE_WRONG_COMPOSITION
        std::vector<BitSet> failedCombinations;
        BitSet              failedCombinationsMask=0uLL;
#endif
        while(composition.generateNext())
        {
            if(composition.is2Power()) // exclude already processed single external bond combinations
                continue;
#ifdef EXCLUDE_WRONG_COMPOSITION
// OPTIMISATION. reduce amount of generated seeds and match calls
// 2120 instead of 2208 match calls on small test. 43 wrongComp-s, 83 rejected
            if(failedCombinationsMask & composition.getBitSet()) // possibly exists in the list
            {
                bool compositionWrong = false;
                for(std::vector<BitSet>::const_iterator failed = failedCombinations.begin(); 
                                                        failed != failedCombinations.end() ; failed++)
                 if(*failed == (*failed & composition.getBitSet())) // combination includes failed combination
                {
                    compositionWrong = true;
                    break;
                }
                if(compositionWrong)
                {
#ifdef VERBOSE_STATISTICS_ON
                    ++extstat.WrongCompositionRejected;
#endif
                    continue;
                }
            }
#endif
#ifdef VERBOSE_STATISTICS_ON
            ++extstat.Seed;
#endif
            Seed seed;
            seed.createFromParent(this);
            newAtomsMap.clear();

            for(unsigned i=0; i<newBonds.size(); i++)
             if(composition.isSet(i))
            {
                const NewBond* nbi = & newBonds[i];
                unsigned aIdx = nbi->EndAtomIdx;    // existed in this parent seed (ring) or -1
                if(-1 == aIdx) // new atom
                {
                    std::map<unsigned, unsigned>::const_iterator nai = newAtomsMap.find(nbi->NewAtomIdx);         // check RING
                    if(newAtomsMap.end() == nai)
                    {
                        const Atom* end_atom = nbi->NewAtom;//qmol.getAtomWithIdx(nbi->NewAtomIdx);
                        aIdx = seed.addAtom(end_atom);
                        newAtomsMap[nbi->NewAtomIdx] = aIdx;    // store new possible ring end point
                    }
                    else
                        aIdx = nai->second;
                }

                const Bond* src_bond = qmol.getBondWithIdx(nbi->BondIdx);
                seed.addBond(src_bond);
            }
            seed.computeRemainingSize(qmol);
            if( ! seed.canGrowBiggerThan(mcs.getMaxNumberBonds(), mcs.getMaxNumberAtoms()) )   // prune(). // seed too small
            {
#ifdef VERBOSE_STATISTICS_ON
                ++extstat.RemainingSizeRejected;
#endif
            }
            else
            {
//1                seedCorrupted = 
                bool found = mcs.checkIfMatchAndAppend(seed);

                if(!found)
                {
#ifdef EXCLUDE_WRONG_COMPOSITION  // if seed does not matched it is possible to exclude this FAILED combination for performance improvement
                    failedCombinations.push_back(composition.getBitSet());
                    failedCombinationsMask &= composition.getBitSet();
#ifdef VERBOSE_STATISTICS_ON
                    ++extstat.WrongCompositionDetected;
#endif
#endif
                }
            }
        }
    }
    GrowingStage = -1; //finished
}

void Seed::computeRemainingSize(const ROMol& qmol)//, const std::vector<char>& excludedBonds)
{
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
    for(unsigned seedAtomIdx = LastAddedAtomsBeginIdx; seedAtomIdx < getNumAtoms(); seedAtomIdx++) // just now added new border vertices (candidates for future growing)
    {
        const Atom* atom = MoleculeFragment.Atoms[seedAtomIdx];
        ROMol::OEDGE_ITER beg,end;
        for(boost::tie(beg,end) = qmol.getAtomBonds(atom); beg!=end; beg++)  // all bonds from MoleculeFragment.Atoms[srcAtomIdx] 
        {
            const Bond& bond = *(qmol[*beg]);
            if( ! visitedBonds[bond.getIdx()])
            {
                ++RemainingBonds;
                visitedBonds[bond.getIdx()] = true;
                unsigned end_atom_idx = (MoleculeFragment.AtomsIdx[seedAtomIdx] == bond.getBeginAtomIdx()) ? bond.getEndAtomIdx() : bond.getBeginAtomIdx();
                if( ! visitedAtoms[end_atom_idx]) // check RING/CYCLE
                {
                    ++RemainingAtoms;
                    visitedAtoms[end_atom_idx] = true;
                    end_atom_stack.push_back(end_atom_idx);
                }
            }
        }
    }
    // 2. go deep
    while(!end_atom_stack.empty())
    {
        unsigned ai = end_atom_stack.back();
                      end_atom_stack.pop_back();
        const Atom* atom = qmol.getAtomWithIdx(ai);
        ROMol::OEDGE_ITER beg,end;
        for(boost::tie(beg,end) = qmol.getAtomBonds(atom); beg!=end; beg++)  // all bonds from end_atom
        {
            const Bond& bond = *(qmol[*beg]);
            if( ! visitedBonds[bond.getIdx()])
            {
                ++RemainingBonds;
                visitedBonds[bond.getIdx()] = true;
                unsigned end_atom_idx = (ai == bond.getBeginAtomIdx()) ? bond.getEndAtomIdx() : bond.getBeginAtomIdx();
                if( ! visitedAtoms[end_atom_idx]) // check RING/CYCLE
                {
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
