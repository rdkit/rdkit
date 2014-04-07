#include <list>
#include <algorithm>
#include <math.h>
#include "../QueryAtom.h"
#include "../QueryBond.h"
#include "../SmilesParse/SmilesWrite.h"
#include "../SmilesParse/SmartsWrite.h"
//#include "../SmilesParse/SmilesParse.h"
//#include "../Substruct/SubstructMatch.h"
#include "SubstructMatchCustom.h"
#include "MaximumCommonSubgraph.h"

#ifdef VERBOSE_STATISTICS_ON
    ExecStatistics extstat;
#endif
namespace RDKit
{
 namespace FMCS
 {
    bool ConsoleOutputEnabled = false;
}}

namespace RDKit
{
 namespace FMCS
 {

    struct LabelDefinition
    {
        unsigned ItemIndex;   // item with this label value
        unsigned Value;
        LabelDefinition() : ItemIndex(-1), Value(-1) {}
        LabelDefinition(unsigned i, unsigned value) : ItemIndex(i), Value(value) {}
    };
        
MaximumCommonSubgraph::MaximumCommonSubgraph(const MCSParameters* params)
{
    Parameters = ( 0 != params ? *params : MCSParameters());
    if (Parameters.ProgressCallback == MCSProgressCallbackTimeout)
        Parameters.ProgressCallbackUserData = &To;
        To = time(0);
}

static
bool molPtr_NumBondLess (const ROMol* l, const ROMol* r)    // need for sorting the source molecules by size
{
    return l->getNumBonds() < r->getNumBonds();
}

void MaximumCommonSubgraph::init()
{
    QueryMolecule = Molecules.front();

    void* userData = Parameters.CompareFunctionsUserData;
    if(Parameters.BondCompareParameters.CompleteRingsOnly || Parameters.BondCompareParameters.RingMatchesRingOnly)
    {
#ifdef FAST_SUBSTRUCT_CACHE
        RingMatchTables.init(QueryMolecule);
        Parameters.CompareFunctionsUserData = &RingMatchTables;
#endif
    }
#ifdef FAST_SUBSTRUCT_CACHE
#ifdef PRECOMPUTED_TABLES_MATCH
        // fill out RingMatchTables for check cache Hash collision by checking match a part of Query to Query
        if(!userData // predefined functor - compute RingMatchTable for all targets
        && (Parameters.BondCompareParameters.CompleteRingsOnly || Parameters.BondCompareParameters.RingMatchesRingOnly))
            RingMatchTables.computeRingMatchTable(QueryMolecule, QueryMolecule, Parameters);

        // fill out match tables
        size_t nq = QueryMolecule->getNumAtoms();
        QueryAtomMatchTable.resize(nq, nq);
        for(size_t aj = 0; aj < nq; aj++)
            for(size_t ai = 0; ai < nq; ai++)
                QueryAtomMatchTable.set(ai, aj, Parameters.AtomTyper(Parameters.AtomCompareParameters, 
                    *QueryMolecule, ai, *QueryMolecule, aj, Parameters.CompareFunctionsUserData));

        nq = QueryMolecule->getNumBonds();
        QueryBondMatchTable.resize(nq, nq);
        for(size_t aj = 0; aj < nq; aj++)
            for(size_t ai = 0; ai < nq; ai++)
                QueryBondMatchTable.set(ai, aj, Parameters.BondTyper(Parameters.BondCompareParameters,
                    *QueryMolecule, ai, *QueryMolecule, aj, Parameters.CompareFunctionsUserData));
#endif // PRECOMPUTED_TABLES_MATCH
        // Compute label values based on current functor and parameters for code Morgan correct computation.
        unsigned currentLabelValue = 1;
        std::vector<LabelDefinition> labels;
        nq = QueryMolecule->getNumAtoms();
        QueryAtomLabels.resize(nq);
        for(size_t ai = 0; ai < nq; ai++)
        {
            if(MCSAtomCompareAny == Parameters.AtomTyper) // predefined functor without atom compare parameters
                QueryAtomLabels[ai] = 1;
            else
            {
                const Atom* atom = QueryMolecule->getAtomWithIdx(ai);
                if(MCSAtomCompareElements == Parameters.AtomTyper) // predefined functor without atom compare parameters
                    QueryAtomLabels[ai] = atom->getAtomicNum()|(Parameters.AtomCompareParameters.MatchValences ? (atom->getTotalValence()>>8):0);
                else if(MCSAtomCompareIsotopes == Parameters.AtomTyper) // predefined functor without atom compare parameters
                    QueryAtomLabels[ai] = atom->getAtomicNum()|(atom->getIsotope()>>8)|(Parameters.AtomCompareParameters.MatchValences ? (atom->getTotalValence()>>16):0);
                else // custom user defined functor
                {
                    QueryAtomLabels[ai] = -1;
                    for(size_t i = 0; i < labels.size(); i++)
                     if(Parameters.AtomTyper(Parameters.AtomCompareParameters, 
                         *QueryMolecule, labels[i].ItemIndex, *QueryMolecule, ai, userData)) // equal itoms
                    {
                        QueryAtomLabels[ai] = labels[i].Value;
                        break;
                    }
                    if(-1 == QueryAtomLabels[ai]) // not found -> create new label
                    {
                        QueryAtomLabels[ai] = ++currentLabelValue;
                        labels.push_back(LabelDefinition(ai, currentLabelValue));
                    }
                }
            }
        }
        labels.clear();
        currentLabelValue = 1;
        nq = QueryMolecule->getNumBonds();
        QueryBondLabels.resize(nq);
        for(size_t aj = 0; aj < nq; aj++)
         {
            const Bond* bond = QueryMolecule->getBondWithIdx(aj);
            unsigned ring = 0;
            if(Parameters.BondCompareParameters.CompleteRingsOnly || Parameters.BondCompareParameters.RingMatchesRingOnly)
            {
                ring = RingMatchTables.isQueryBondInRing(aj) ? 0 : 1;  // is bond in ring
            }
            if(MCSBondCompareAny == Parameters.BondTyper) // predefined functor without atom compare parameters
                QueryBondLabels[aj] = 1 | (ring>>8);
            else if(MCSBondCompareOrderExact == Parameters.BondTyper) // predefined functor without compare parameters
                QueryBondLabels[aj] = (bond->getBondType() + 1) | (ring>>8);
            else if(MCSBondCompareOrder == Parameters.BondTyper)     // predefined functor, ignore Aromatization
            {
                unsigned order = bond->getBondType();
                if(Bond::AROMATIC == order || Bond::ONEANDAHALF == order) // ignore Aromatization
                    order = Bond::SINGLE;
                else if(Bond::TWOANDAHALF == order)
                    order = Bond::DOUBLE;
                else if(Bond::THREEANDAHALF == order)
                    order = Bond::TRIPLE;
                else if(Bond::FOURANDAHALF == order)
                    order = Bond::QUADRUPLE;
                else if(Bond::FIVEANDAHALF == order)
                    order = Bond::QUINTUPLE;
                QueryBondLabels[aj] = (order + 1) | (ring>>8);
            }
            else // custom user defined functor
            {
                QueryBondLabels[aj] = -1;
                for(size_t i = 0; i < labels.size(); i++)
                 if(Parameters.BondTyper(Parameters.BondCompareParameters, 
                                        *QueryMolecule, labels[i].ItemIndex, *QueryMolecule, aj, userData)) // equal bonds + ring ...
                {
                    QueryBondLabels[aj] = labels[i].Value;
                    break;
                }
                if(-1 == QueryAtomLabels[aj]) // not found -> create new label
                {
                    QueryBondLabels[aj] = ++currentLabelValue;
                    labels.push_back( LabelDefinition(aj, currentLabelValue));
                }
            }
        }
#endif
    Targets.resize(Molecules.size()-1);
    size_t i=0;
    for(std::vector<const ROMol*>::iterator it = Molecules.begin()+1; it != Molecules.end(); it++, i++)
    {
        Targets[i].Molecule = *it;
        // build Target Topology ADD ATOMs
        size_t j=0;    // current item
        for(ROMol::ConstAtomIterator a = Targets[i].Molecule->beginAtoms(); a != Targets[i].Molecule->endAtoms(); a++, j++)
        {
            Targets[i].Topology.addAtom((*a)->getIdx());
#ifdef FAST_INCREMENTAL_MATCH
            Targets[i].AtomAdjacency.resize((*it)->getNumAtoms());
            ROMol::OEDGE_ITER beg,end;
            for(boost::tie(beg,end) = (*it)->getAtomBonds(*a); beg!=end; beg++)  // all bonds from the atom
            {
                AtomAdjacency aa;
                aa.Bond = (*(*it))[*beg].get();
                aa.BondIdx = aa.Bond->getIdx();
                aa.ConnectedAtomIdx = aa.Bond->getEndAtomIdx();
                if((*a)->getIdx() == aa.ConnectedAtomIdx)
                    aa.ConnectedAtomIdx = aa.Bond->getBeginAtomIdx();
                Targets[i].AtomAdjacency[j].push_back(aa);
            }
#endif // FAST_INCREMENTAL_MATCH
        }
        // build Target Topology ADD BONDs
        for(ROMol::ConstBondIterator b = Targets[i].Molecule->beginBonds(); b != Targets[i].Molecule->endBonds(); b++)
        {
            const Bond* bond = *b;
            unsigned ii = bond->getBeginAtomIdx();
            unsigned jj = bond->getEndAtomIdx();
            Targets[i].Topology.addBond((*b)->getIdx(), ii, jj);
        }

        // fill out RingMatchTables
        if(!userData // predefined functor - compute RingMatchTable for all targets
        && (Parameters.BondCompareParameters.CompleteRingsOnly || Parameters.BondCompareParameters.RingMatchesRingOnly))
        {
#ifdef FAST_SUBSTRUCT_CACHE
            RingMatchTables.addTargetBondRingsIndeces(Targets[i].Molecule);
            RingMatchTables.computeRingMatchTable(QueryMolecule, Targets[i].Molecule, Parameters);
#endif
        }

#ifdef PRECOMPUTED_TABLES_MATCH
        // fill out match tables
        size_t nq = QueryMolecule->getNumAtoms();
        size_t nt = (*it)->getNumAtoms();
        Targets[i].AtomMatchTable.resize(nq, nt);

        for(size_t aj = 0; aj < nt; aj++)
            for(size_t ai = 0; ai < nq; ai++)
                Targets[i].AtomMatchTable.set(ai, aj, Parameters.AtomTyper(Parameters.AtomCompareParameters, 
                    *QueryMolecule, ai, *Targets[i].Molecule, aj, Parameters.CompareFunctionsUserData));

        nq = QueryMolecule->getNumBonds();
        nt = (*it)->getNumBonds();
        Targets[i].BondMatchTable.resize(nq, nt);
        for(size_t aj = 0; aj < nt; aj++)
            for(size_t ai = 0; ai < nq; ai++)
                Targets[i].BondMatchTable.set(ai, aj, Parameters.BondTyper(Parameters.BondCompareParameters,
                    *QueryMolecule, ai, *Targets[i].Molecule, aj, Parameters.CompareFunctionsUserData));
#endif
    }

    Parameters.CompareFunctionsUserData = userData; // restore
}

void MaximumCommonSubgraph::makeInitialSeeds()
{
    // build a set of initial seeds as "all" single bonds from query molecule
    std::vector<bool> excludedBonds(QueryMolecule->getNumBonds());
    for(size_t i = 0; i < excludedBonds.size(); i++)
        excludedBonds[i] = false;

    Seeds.clear();

    //R1 additional performance OPTIMISATION
    //if(Parameters.BondCompareParameters.CompleteRingsOnly)
    // disable all mismatched rings, and do not generate initial seeds from such disabled bonds
    //  for(  rings .....) for(i......) 
    //   if(mismatched) excludedBonds[i.......] = true;

    for(RWMol::ConstBondIterator bi = QueryMolecule->beginBonds(); bi != QueryMolecule->endBonds(); bi++)
    {
        //R1 additional performance OPTIMISATION
        //if(excludedBonds[(*bi)->getIdx()])
        //    continue;
        Seed seed;
#ifdef FAST_INCREMENTAL_MATCH
        seed.MatchResult.resize(Targets.size());
#endif

#ifdef VERBOSE_STATISTICS_ON
        ++extstat.Seed;
        ++extstat.InitialSeed;
#endif
        seed.addAtom((*bi)->getBeginAtom());
        seed.addAtom((*bi)->getEndAtom());
        seed.ExcludedBonds = excludedBonds; // all bonds from first to current
        seed.addBond (*bi);
        excludedBonds[(*bi)->getIdx()] = true;

        seed.computeRemainingSize(*QueryMolecule);

        if( ! checkIfMatchAndAppend(seed))
        {
#ifdef VERBOSE_STATISTICS_ON
            ++extstat.MismatchedInitialSeed;
#endif
            // optionally remove all such bonds from all targets TOPOLOGY where it exists.
            //..........

            // disable (mark as already processed) mismatched bond in all seeds
            for(SeedSet::iterator si = Seeds.begin(); si != Seeds.end(); si++)
                si->ExcludedBonds[(*bi)->getIdx()] = true;

        }
    }
    Seeds.MaxBonds = 1;
    Seeds.MaxAtoms = 2;
}

void MaximumCommonSubgraph::growSeeds(MolFragment& mcsIdx, MCSResult& res)
{
    unsigned steps = 99999; // steps from last progress callback call. call it immediately in the begining

    // Find MCS -- SDF Seed growing OPTIMISATION (it works in 3 times faster)
    while(!Seeds.empty())
    {
        ++steps;
        extstat.TotalSteps++;
        SeedSet::iterator si = Seeds.begin();

        si->grow(*this, *QueryMolecule);

        const Seed& fs = Seeds.front();
        // bigger substructure found
        if((!Parameters.MaximizeBonds && (fs.getNumAtoms() > res.NumAtoms || (fs.getNumAtoms() == res.NumAtoms && fs.getNumBonds() > res.NumBonds)))
         ||( Parameters.MaximizeBonds && (fs.getNumBonds() > res.NumBonds || (fs.getNumBonds() == res.NumBonds && fs.getNumAtoms() > res.NumAtoms)))
         )
        {
            extstat.MCSFoundStep = extstat.TotalSteps;
            res.NumAtoms    = fs.getNumAtoms();
            res.NumBonds    = fs.getNumBonds();
            mcsIdx.Atoms    = fs.MoleculeFragment.Atoms;
            mcsIdx.Bonds    = fs.MoleculeFragment.Bonds;
            mcsIdx.AtomsIdx = fs.MoleculeFragment.AtomsIdx;
            mcsIdx.BondsIdx = fs.MoleculeFragment.BondsIdx;
//TMP DEBUG
/*    std::cout<<"MCS atoms=(";
    for(unsigned seedAtomIdx = 0; seedAtomIdx < mcs.getNumAtoms(); seedAtomIdx++)
    {
        const Atom* atom = mcs.MoleculeFragment.Atoms[seedAtomIdx];
        std::cout<<atom->getIdx()<<", ";
    }
    std::cout<<") new atoms=[";
    for(unsigned seedAtomIdx = mcs.LastAddedAtomsBeginIdx; seedAtomIdx < mcs.getNumAtoms(); seedAtomIdx++)
    {
        const Atom* atom = mcs.MoleculeFragment.Atoms[seedAtomIdx];
        std::cout<<atom->getIdx()<<", ";
    }
    std::cout<<"] bonds=(";
    for(int i=0; i<mcsIdx.BondsIdx.size(); i++)
        std::cout<<mcsIdx.BondsIdx[i]<<", ";
    std::cout<<") Size="<< mcs.getNumAtoms() <<", "<< mcs.getNumBonds() <<" Remain=" << mcs.RemainingAtoms <<", "<<mcs.RemainingBonds<<" = ";
    std::cout<< MolFragmentToSmiles(*QueryMolecule, *(const std::vector<int>*) &mcs.MoleculeFragment.AtomsIdx, (const std::vector<int>*) &mcs.MoleculeFragment.BondsIdx) <<"\n";  // unsigned
*/
/*
if(0==si->MoleculeFragment.BondsIdx[0])
{
    std::cout<<"\n"
        <<" -----------------"<<(void*)&*si 
        <<"Remaining "<<si->RemainingBonds<< ",  " << si->RemainingAtoms
        <<": LastAddedAtomsBeginIdx = "<<si->LastAddedAtomsBeginIdx<<", LastAddedBondsBeginIdx = "<<si->LastAddedBondsBeginIdx<<"\n";
    for(size_t i = 0; i < si->MoleculeFragment.Bonds.size(); i++)
        std::cout << i << " "<<si->MoleculeFragment.Bonds[i]->getIdx()<<" : "
                        <<" "<<si->MoleculeFragment.Bonds[i]->getBeginAtom()->getIdx()
                        <<" "<<si->MoleculeFragment.Bonds[i]->getEndAtom()->getIdx()
                        <<"\n";
}
*/            
/*
#ifdef xxxWIN32 // TEMP DEBUG !!!!!!!!!!!!!
  if(30 == res.NumAtoms)
  {
    std::cout<<"-------------\n";
      for(int i=0; i<mcs.getNumBonds(); i++)
        std::cout<<i<<" bIdx=" << mcs.MoleculeFragment.BondsIdx[i]//<<"\n";
            <<" =("<< QueryMolecule->getBondWithIdx(mcs.MoleculeFragment.BondsIdx[i])->getBeginAtomIdx()
            <<" , "<< QueryMolecule->getBondWithIdx(mcs.MoleculeFragment.BondsIdx[i])->getEndAtomIdx()<<")\n";
    std::cout<<"-------------\n";
    for(int i=0; i<30; i++)
    {
        const Atom* a = QueryMolecule->getAtomWithIdx(mcs.MoleculeFragment.AtomsIdx[i]);
        std::cout<<i<<" aIdx=" << mcs.MoleculeFragment.AtomsIdx[i]<<" : "<<a->getAtomicNum()<<"\n";
        if(a->getAtomicNum()==7)   // N
        {
            ROMol::OEDGE_ITER beg,end;
            for(boost::tie(beg,end) = QueryMolecule->getAtomBonds(a); beg!=end; beg++)
            {
                const Bond* bond = &*((*QueryMolecule)[*beg]);
                std::cout<<"  bond="<<bond->getIdx()<<" =("<< bond->getBeginAtomIdx()
                    <<" , "<< bond->getEndAtomIdx()<<")\n";
            }
        }
    }
    std::cout<<"-------------\n";
  }
#endif
*/
        }

        if(-1 == si->GrowingStage) //finished
            Seeds.erase(si);

        if(Parameters.ProgressCallback && (steps > 777))// || res.NumAtoms > Stat.NumAtoms))
        {
            steps = 0;
            Stat.NumAtoms = res.NumAtoms;
            Stat.NumBonds = res.NumBonds;
#ifdef VERBOSE_STATISTICS_ON
            Stat.SeedProcessed = extstat.Seed;
#endif
            if(!Parameters.ProgressCallback(Stat, Parameters, Parameters.ProgressCallbackUserData))
            {
                res.Canceled = true;
                break;
            }
        }
    }

}

    struct AtomMatch // for each seed atom (matched)
    {
        unsigned QueryAtomIdx;
        unsigned TargetAtomIdx;
        AtomMatch() : QueryAtomIdx(-1), TargetAtomIdx(-1) {}
    };
    typedef std::vector<AtomMatch> AtomMatchSet; 

std::string MaximumCommonSubgraph::generateSMARTS(const MolFragment& mcsIdx)
{
    // match the result MCS with all targets to check if it is exact match or template
    Seed seed; // result MCS
    seed.ExcludedBonds.resize(QueryMolecule->getNumBonds());
    for(size_t i = 0; i < seed.ExcludedBonds.size(); i++)
        seed.ExcludedBonds[i] = false;
    std::vector<AtomMatchSet> atomMatchResult(Targets.size());
    std::vector<unsigned> atomIdxMap(QueryMolecule->getNumAtoms());
    std::vector<std::map<unsigned, const Bond*> > bondMatchSet (mcsIdx.Bonds.size()); //key is unique BondType
    std::vector<std::map<unsigned, const Atom*> > atomMatchSet (mcsIdx.Atoms.size()); //key is unique atomic number
//isotope, mass, charge, Hs    std::vector<std::map<unsigned, const Atom*> > atomMatchSet (mcsIdx.Atoms.size()); //key is unique atomic number

    for(std::vector<const Atom*>::const_iterator atom = mcsIdx.Atoms.begin(); atom != mcsIdx.Atoms.end(); atom++)
    {
        atomIdxMap[(*atom)->getIdx()] = seed.getNumAtoms();
        seed.addAtom((*atom));
    }
    for(std::vector<const Bond*>::const_iterator bond = mcsIdx.Bonds.begin(); bond != mcsIdx.Bonds.end(); bond++)
        seed.addBond((*bond));

    unsigned itarget = 0;
    for(std::vector<Target>::const_iterator tag = Targets.begin(); tag != Targets.end(); tag++, itarget++)
    {
        match_V_t match;    // THERE IS NO Bonds match INFO !!!! 
        bool target_matched = 
#ifdef PRECOMPUTED_TABLES_MATCH
            SubstructMatchCustomTable(tag->Topology, seed.Topology, tag->AtomMatchTable, tag->BondMatchTable, &match);
#else //noticable slowly:
            SubstructMatchCustom(
            tag->Topology.GraphTopology
            , *tag->Molecule
            , seed.GraphTopology
            , *QueryMolecule
            , Parameters.AtomTyper, Parameters.BondTyper
            , Parameters.AtomCompareParameters, Parameters.BondCompareParameters, Parameters.CompareFunctionsUserData
            );
#endif
        atomMatchResult[itarget].resize(seed.getNumAtoms());
        for(match_V_t::const_iterator mit = match.begin(); target_matched && mit != match.end(); mit++)
        {
            unsigned ai = mit->first;  // SeedAtomIdx
            atomMatchResult[itarget][ai].QueryAtomIdx  = seed.Topology[mit->first];
            atomMatchResult[itarget][ai].TargetAtomIdx = tag->Topology[mit->second];
            const Atom* ta = tag->Molecule->getAtomWithIdx(tag->Topology[mit->second]);
            if(ta && ta->getAtomicNum() != seed.MoleculeFragment.Atoms[ai]->getAtomicNum())
                atomMatchSet[ai][ta->getAtomicNum()] = ta; // add
        }
        // AND BUILD BOND MATCH INFO
        unsigned bi=0;
        for(std::vector<const Bond*>::const_iterator bond = mcsIdx.Bonds.begin(); target_matched && bond != mcsIdx.Bonds.end(); bond++, bi++)
        {
            unsigned i = atomIdxMap[(*bond)->getBeginAtomIdx()];
            unsigned j = atomIdxMap[(*bond)->getEndAtomIdx()];
            unsigned ti= atomMatchResult[itarget][i].TargetAtomIdx;
            unsigned tj= atomMatchResult[itarget][j].TargetAtomIdx;
            const Bond* tBond = tag->Molecule->getBondBetweenAtoms(ti, tj);
            if(tBond && (*bond)->getBondType() != tBond->getBondType())
                bondMatchSet[bi] [tBond->getBondType()] = tBond; // add

        }
    }

    // Generate result's SMARTS

    RWMol mol;  // create molecule from MCS for MolToSmarts()
    unsigned ai = 0;  // SeedAtomIdx
    for(std::vector<const Atom*>::const_iterator atom = mcsIdx.Atoms.begin(); atom != mcsIdx.Atoms.end(); atom++, ai++)
    {
        QueryAtom* a = new QueryAtom(*(*atom));  // generate [#6] instead of C or c !
        //for all atomMatchSet[ai] items add atom query to template
        for(std::map<unsigned, const Atom*>::const_iterator am = atomMatchSet[ai].begin(); am != atomMatchSet[ai].end(); am++)
        {
            ATOM_OR_QUERY* a2 = new ATOM_OR_QUERY();
            a2->addChild(QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(makeAtomNumQuery((*atom)->getAtomicNum())));
            a2->addChild(QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(makeAtomNumQuery(am->second->getAtomicNum())));      
            a2->setDescription("AtomOr");
//ATOM_EQUALS_QUERY *makeAtomIsotopeQuery(int what)
//....
            a->setQuery(a2);
        }
        mol.addAtom(a, true);
    }
    unsigned bi = 0;  // Seed Idx
    for(std::vector<const Bond*>::const_iterator bond = mcsIdx.Bonds.begin(); bond != mcsIdx.Bonds.end(); bond++, bi++)
    {
        unsigned i = atomIdxMap[(*bond)->getBeginAtomIdx()];
        unsigned j = atomIdxMap[(*bond)->getEndAtomIdx()];
        QueryBond* b = new QueryBond(*(*bond));
        // add OR template if need
        for(std::map<unsigned, const Bond*>::const_iterator bm = bondMatchSet[bi].begin(); bm != bondMatchSet[bi].end(); bm++)
        {
            BOND_OR_QUERY* b2 = new BOND_OR_QUERY();
            b2->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery((*bond)->getBondType())));
            b2->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(bm->second->getBondType())));      
            b2->setDescription("BondOr");
            b->setQuery(b2);
        }
        b->setBeginAtomIdx(i);
        b->setEndAtomIdx(j);
        mol.addBond(b, true);
    }
    return MolToSmarts(mol, true);
}

MCSResult MaximumCommonSubgraph::find(const std::vector<ROMOL_SPTR>& src_mols)
{
    clear();
    MCSResult res;

    if(src_mols.size() < 2)
        throw std::runtime_error("FMCS. Invalid argument. mols.size() must be at least 2");
    if (Parameters.Threshold > 1.0)
        throw std::runtime_error("FMCS. Invalid argument. Parameter Threshold must be 1.0 or less.");

    ThresholdCount = (unsigned) ceil((src_mols.size()-1) * Parameters.Threshold);   // minimal required number of matched targets
    if (ThresholdCount < 1) // at least one target
        ThresholdCount = 1;
    if (ThresholdCount > src_mols.size()-1) // max all targets
        ThresholdCount = src_mols.size()-1;
    
    // Selecting CompleteRingsOnly option also enables --ring-matches-ring-only. ring--ring and chain bonds only match chain bonds.
    if(Parameters.BondCompareParameters.CompleteRingsOnly)
        Parameters.BondCompareParameters.RingMatchesRingOnly = true;

    for(std::vector<ROMOL_SPTR>::const_iterator it = src_mols.begin(); it != src_mols.end(); it++)
    {
        Molecules.push_back((*it).get());
        if(!Molecules.back()->getRingInfo()->isInitialized())
            Molecules.back()->getRingInfo()->initialize();  // but do not fill out !!!
    }

    // sort source set of molecules by their 'size' and assume the smallest molecule as a query
    std::stable_sort(Molecules.begin(), Molecules.end(), molPtr_NumBondLess);

    init();
    makeInitialSeeds();

/* // TEMP DEBUG !!!!!!!!!!!!!
  {
    std::cout<<"---Bonds----------\n";
      for(int i=0; i<QueryMolecule->getNumBonds(); i++)
        std::cout<<i<<" "
            <<" =("<< QueryMolecule->getBondWithIdx(i)->getBeginAtomIdx()
            <<" , "<< QueryMolecule->getBondWithIdx(i)->getEndAtomIdx()<<")\n";
    std::cout<<"---Atoms----------\n";
    for(int i=0; i<QueryMolecule->getNumAtoms(); i++)
    {
        const Atom* a = QueryMolecule->getAtomWithIdx(i);
        std::cout<<i<<" : elem="<<a->getAtomicNum()<<" bonds: ";
        ROMol::OEDGE_ITER beg,end;
        for(boost::tie(beg,end) = QueryMolecule->getAtomBonds(a); beg!=end; beg++)
        {
            const Bond* bond = &*((*QueryMolecule)[*beg]);
            std::cout<<bond->getIdx()<<"=("<< bond->getBeginAtomIdx()<<", "<<bond->getEndAtomIdx()<<") ";
        }
        std::cout<<"\n";
    }
    std::cout<<"-------------\n";
  }
// !!!!!!!!!!!!!!!!!!!!!!
*/
    MolFragment mcsIdx; // current MCS
    growSeeds(mcsIdx, res);

    res.SmartsString = generateSMARTS(mcsIdx);

/* // TMP DEBUG 
    std::cout<<"---MCS Bonds------\n";
    for(int i=0; i<mcsIdx.BondsIdx.size(); i++)
        std::cout<<i<<" bIdx=" << mcsIdx.BondsIdx[i]
            <<" =("<< QueryMolecule->getBondWithIdx(mcsIdx.BondsIdx[i])->getBeginAtomIdx()
            <<" , "<< QueryMolecule->getBondWithIdx(mcsIdx.BondsIdx[i])->getEndAtomIdx()<<")\n";
// TMP DEBUG 
//---------------------------
std::cout<<"Query      "<< MolToSmiles(*QueryMolecule) <<"\n";
std::cout<<"MCS Smiles "<< MolFragmentToSmiles(*QueryMolecule, *(const std::vector<int>*) &mcsIdx.AtomsIdx, (const std::vector<int>*) &mcsIdx.BondsIdx) <<"\n";  // unsigned
*/
/*
    unsigned itarget = 0;
    for(std::vector<Target>::const_iterator tag = Targets.begin(); tag != Targets.end(); tag++, itarget++)
    {
        MatchVectType match;
        bool target_matched = SubstructMatch(*tag->Molecule, *QueryMolecule, match);
        std::cout<<"Target "<< itarget+1 << (target_matched ? " matched" : " MISMATCHED") <<"\n";
    }
*/
//---------------------------

#ifdef VERBOSE_STATISTICS_ON
if(ConsoleOutputEnabled)
{
    std::cout << "STATISTICS:\n";
    std::cout << "Total Growing Steps  = " << extstat.TotalSteps<<",  MCS found on "<<extstat.MCSFoundStep<<" step\n";
    std::cout << "Initial   Seeds      = " << extstat.InitialSeed << ",  Mismatched " << extstat.MismatchedInitialSeed<<"\n";
    std::cout << "Inspected Seeds      = " << extstat.Seed<<"\n";
    std::cout << "Rejected by BestSize = " << extstat.RemainingSizeRejected << "\n";
#ifdef EXCLUDE_WRONG_COMPOSITION   
    std::cout << "Rejected by WrongComposition = " << extstat.WrongCompositionRejected
                        << " [ " << extstat.WrongCompositionDetected << " Detected ]\n";
#endif
    std::cout << "MatchCheck Seeds     = " << extstat.SeedCheck      <<"\n";
    std::cout //<< "\n"
              << "     MatchCalls = " << extstat.MatchCall      <<"\n"
              << "     MatchFound = " << extstat.MatchCallTrue  <<"\n";
#ifdef FAST_INCREMENTAL_MATCH
    std::cout << " fastMatchCalls = " << extstat.FastMatchCall <<"\n"
              << " fastMatchFound = " << extstat.FastMatchCallTrue <<"\n";
    std::cout << " slowMatchCalls = " << extstat.MatchCall     - extstat.FastMatchCallTrue <<"\n"
              << " slowMatchFound = " << extstat.MatchCallTrue - extstat.FastMatchCallTrue <<"\n";
#endif
#ifdef PRECOMPUTED_TABLES_MATCH
    std::cout << "--- USED PreComputed Match TABLES ! ---\n";
#endif
#ifdef VERBOSE_STATISTICS_FASTCALLS_ON
    std::cout << "AtomFunctorCalls = " << extstat.AtomFunctorCalls << "\n";
    std::cout << "BondCompareCalls = " << extstat.BondCompareCalls << "\n";
#endif
//    std::cout << "\n";
    std::cout << "  DupCacheFound = " << extstat.DupCacheFound 
        <<"   "<< extstat.DupCacheFoundMatch<<" matched, "
        <<extstat.DupCacheFound - extstat.DupCacheFoundMatch <<" mismatched\n";
#ifdef FAST_SUBSTRUCT_CACHE
    std::cout << "HashCache size  = " << HashCache.keyssize() << " keys\n";
    std::cout << "HashCache size  = " << HashCache.fullsize() << " entries\n";
    std::cout << "FindHashInCache = " << extstat.FindHashInCache << "\n";
    std::cout << "HashFoundInCache= " << extstat.HashKeyFoundInCache << "\n";
    std::cout << "ExactMatchCalls = " << extstat.ExactMatchCall  <<"\n"
              << "ExactMatchFound = " << extstat.ExactMatchCallTrue <<"\n";
#endif
}
#endif
//---------------------
    clear();
    return res;
}

bool MaximumCommonSubgraph::checkIfMatchAndAppend(Seed& seed)
{
#ifdef TRACE_ON
        TRACE() << "CHECK ";    // print out time
        for(std::vector<const Bond*>::const_iterator bi = seed.MoleculeFragment.Bonds.begin(); bi != seed.MoleculeFragment.Bonds.end(); bi++)
            TRACE(0) << (*bi)->getIdx() << " ";
        TRACE(0) << "\n";
#endif
#ifdef VERBOSE_STATISTICS_ON
    ++extstat.SeedCheck;
#endif

    bool foundInCache = false;
    bool foundInDupCache = false;

#ifdef DUP_SUBSTRUCT_CACHE
    if(DuplicateCache.find(seed.DupCacheKey, foundInCache))
    {
    // duplicate found. skip match() but store both seeds, because they will grow by different paths !!!
        #ifdef VERBOSE_STATISTICS_ON
            extstat.DupCacheFound++;
            extstat.DupCacheFoundMatch += foundInCache ? 1 : 0;
        #endif
        if(!foundInCache) // mismatched !!!
            return false;
    }
    foundInDupCache = foundInCache;
#endif
#ifdef FAST_SUBSTRUCT_CACHE
    SubstructureCache::HashKey      cacheKey;
    SubstructureCache::TIndexEntry* cacheEntry = 0;
    bool cacheEntryIsValid = false;
    if(!foundInCache)
    {
    #ifdef VERBOSE_STATISTICS_ON
        ++extstat.FindHashInCache;
    #endif
        cacheEntry = HashCache.find(seed, QueryAtomLabels, QueryBondLabels, cacheKey);
        cacheEntryIsValid = true;
        if(cacheEntry)   // possibly found. check for hash collision
        {
            #ifdef VERBOSE_STATISTICS_ON
                ++extstat.HashKeyFoundInCache;
            #endif
            // check hash collisions (time +3%):
            for(SubstructureCache::TIndexEntry::const_iterator g = cacheEntry->begin(); !foundInCache && g != cacheEntry->end(); g++)
            {
                if(g->m_vertices.size() != seed.getNumAtoms() || g->m_edges.size() != seed.getNumBonds())
                    continue;
            #ifdef VERBOSE_STATISTICS_ON
                ++extstat.ExactMatchCall;
            #endif
            // EXACT MATCH
            #ifdef PRECOMPUTED_TABLES_MATCH
                foundInCache = SubstructMatchCustomTable((*g), seed.Topology, QueryAtomMatchTable, QueryBondMatchTable);
            #else //..................
            #endif
            #ifdef VERBOSE_STATISTICS_ON
                if(foundInCache)
                    ++extstat.ExactMatchCallTrue;
            #endif
            }
        }
    }
    #endif
    bool found = foundInCache;

    if(!found) 
    {
        found = match(seed);
    }

    if(found)  // Store new generated seed, if found in cache or in all(- threshold) targets
    {
        Seed& new_seed = Seeds.add(seed);

    #ifdef DUP_SUBSTRUCT_CACHE
        if(!foundInDupCache && seed.getNumBonds() >= 3)  // only seed with a ring can be duplicated - do not store very small seed in cache
            DuplicateCache.add(seed.DupCacheKey, true);
    #endif
    #ifdef FAST_SUBSTRUCT_CACHE
        if(!foundInCache)
            HashCache.add(seed, cacheKey, cacheEntry);
    #endif
    }
    else
    {
    #ifdef DUP_SUBSTRUCT_CACHE
        if(seed.getNumBonds() > 3)
            DuplicateCache.add(seed.DupCacheKey, false);   //opt. cache mismatched duplicates too
    #endif
    }
    return found;  // new matched seed has been actualy added
}

bool MaximumCommonSubgraph::match(Seed& seed)
{
    unsigned max_miss = Targets.size() - ThresholdCount;
    unsigned missing  = 0;
    unsigned passed   = 0;
    unsigned itarget  = 0;

    for(std::vector<Target>::const_iterator tag = Targets.begin(); tag != Targets.end(); tag++, itarget++)
    {
#ifdef VERBOSE_STATISTICS_ON
    ++extstat.MatchCall;
#endif
        bool target_matched = false;
#ifdef FAST_INCREMENTAL_MATCH
        if(!seed.MatchResult[itarget].empty())
            target_matched = matchIncrementalFast(seed, itarget);
#endif
        if(!target_matched) // slow full match
        {
            match_V_t match;    // THERE IS NO Bonds match INFO !!!! 
            target_matched = 
#ifdef PRECOMPUTED_TABLES_MATCH
            SubstructMatchCustomTable(tag->Topology, seed.Topology, tag->AtomMatchTable, tag->BondMatchTable, &match);
#else //noticable slowly:
            SubstructMatchCustom(
            tag->Topology
            , *tag->Molecule
            , seed.GraphTopology
            , *QueryMolecule
            , Parameters.AtomTyper, Parameters.BondTyper
            , Parameters.AtomCompareParameters, Parameters.BondCompareParameters, Parameters.CompareFunctionsUserData
            );
#endif
#ifdef FAST_INCREMENTAL_MATCH   // save current match info
            seed.MatchResult[itarget].clear();//resize(target_matched ? match.size() : 0);
            for(match_V_t::const_iterator mit = match.begin(); target_matched && mit != match.end(); mit++)
            {
                BondMatch m;
                m.SeedAtomIdx   = mit->first;
                m.QueryAtomIdx  = seed.GraphTopology[mit->first];
                m.TargetAtomIdx = mit->second; //==target.Topology[it->second]);
                seed.MatchResult[itarget].push_back(m); // add for each matched bond
            }
            //AND RESTORE BOND MATCH INFO !!!!!!!!!!!!!!!!!!!!!
#endif
        }

        if(target_matched)
        {
            if(++passed >= ThresholdCount) // it's enought
                break;
        }
        else // mismatched
        {
            if(++missing > max_miss)
                break;
        }
    }
    if(missing <= max_miss)
    {
        #ifdef VERBOSE_STATISTICS_ON
            ++extstat.MatchCallTrue;
        #endif
        return true;
    }
    return false;
}

#ifdef FAST_INCREMENTAL_MATCH
// call it by target, if fail perform full match check
bool MaximumCommonSubgraph::matchIncrementalFast(Seed& seed, unsigned itarget)
{
    // use and update results of previous match stored in the seed
#ifdef VERBOSE_STATISTICS_ON
    ++extstat.FastMatchCall;
#endif
    Target& t = Targets[itarget];
    BondMatchSet& match = seed.MatchResult[itarget];
    size_t previousLen = match.size();
    for(unsigned newBondSeedIdx = seed.LastAddedBondsBeginIdx; newBondSeedIdx < seed.getNumBonds(); newBondSeedIdx++)
    {
        bool found=false;
        const Bond* bond = seed.MoleculeFragment.Bonds[newBondSeedIdx];
        unsigned newBondQueryIdx   = seed.MoleculeFragment.BondsIdx[newBondSeedIdx];
        unsigned i = seed.MoleculeFragment.SeedAtomIdxMap[bond->getBeginAtomIdx()];
        if(i >= seed.LastAddedAtomsBeginIdx) // old atom in the seed
            i = seed.MoleculeFragment.SeedAtomIdxMap[bond->getEndAtomIdx()];
        unsigned newBondSourceAtomSeedIdx   = i; // seed's index of atom from which new bond was added
        const BondMatch& ma = seed.MatchResult[itarget][i];
        unsigned newBondSourceAtomTargetIdx = ma.TargetAtomIdx; // corresponding atom in the target
unsigned newBondSourceAtomQueryIdx  = ma.QueryAtomIdx;
//const Atom* newBondSourceAtomTarget = t.Molecule->getAtomWithIdx(newBondSourceAtomTargetIdx);
        // check all unvisited bonds from newBondSourceAtomTargetIdx
        const AtomAdjacencyList& taAdj = t.AtomAdjacency[newBondSourceAtomTargetIdx];
        for(size_t itaAdji = 0; itaAdji < taAdj.size(); itaAdji++)
        {
            unsigned tbi = taAdj[itaAdji].BondIdx;
            if(t.BondMatchTable.at(newBondQueryIdx, tbi))
            {
/*
                // and check ending atom of the bond too 
                //..........
                bool visited = false;
                //match.find
                //MULTIPLE various bonds for one atom !!!
                for(size_t im = 0; !visited && im < match.size(); im++)
                    if(seed.MatchResult[itarget][im]. == )
                        visited = true;
                if(visited) // visited == targer bons already presents in the match of this seed
                    continue;
                // append to match result
                seed.MatchResult[itarget][]

                found = true;
                break;
*/
            }
        }
    }
/*/-----------------------
//Seed:
unsigned LastAddedBondBeginIdx;
//history for each tag
seedBond, tagBond,
both (seedAtom, tagAtom) of

        Graph q = seed.GraphTopology;
        for(unsigned ib = seed.LastAddedBondBeginIdx; ib < seed.getNumBonds(); ib++)   // all added bonds

            ib
            endAtom // of ib
            unsigned tagBeginAtomIdx = ...; // already verified. matched with seed.
            for() all outgoing bonds except already matched to seed
                check match of bond to ib and endAtom to ia if ia is not already verified (ring back to seed)
        //        for(unsigned srcAtomIdx = seed.LastAddedAtomsBeginIdx; srcAtomIdx < seed.getNumAtoms(); srcAtomIdx++)   // all added atoms

//-----------------------*/
    if(false // match_ok
        )
    {
    #ifdef VERBOSE_STATISTICS_ON
        ++extstat.FastMatchCallTrue;
    #endif
        return true;
    }
    match.resize(previousLen);  // clean up - remove new items only
    return false;
}
#endif // FAST_INCREMENTAL_MATCH

}}   // namespace RDKit
