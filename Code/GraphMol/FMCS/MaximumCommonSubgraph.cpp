//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <list>
#include <algorithm>
#include <cmath>
#include "../QueryAtom.h"
#include "../QueryBond.h"
#include "../SmilesParse/SmilesWrite.h"
#include "../SmilesParse/SmartsWrite.h"
#include "../SmilesParse/SmilesParse.h"
#include "../Substruct/SubstructMatch.h"
#include "SubstructMatchCustom.h"
#include "MaximumCommonSubgraph.h"
#include <boost/graph/adjacency_list.hpp>

namespace RDKit {
namespace FMCS {

struct LabelDefinition {
  unsigned ItemIndex;  // item with this label value
  unsigned Value;
  LabelDefinition() : ItemIndex(NotSet), Value(NotSet) {}
  LabelDefinition(unsigned i, unsigned value) : ItemIndex(i), Value(value) {}
};

MaximumCommonSubgraph::MaximumCommonSubgraph(const MCSParameters* params) {
  Parameters = (nullptr != params ? *params : MCSParameters());
  if (!Parameters.ProgressCallback) {
    Parameters.ProgressCallback = MCSProgressCallbackTimeout;
    Parameters.ProgressCallbackUserData = &To;
  }
  if ((Parameters.AtomCompareParameters.MatchChiralTag ||
       Parameters.BondCompareParameters.MatchFusedRings ||
       Parameters.BondCompareParameters.MatchFusedRingsStrict) &&
      nullptr == Parameters.FinalMatchChecker) {
    Parameters.FinalMatchChecker = FinalMatchCheckFunction;
    if (Parameters.AtomCompareParameters.MatchChiralTag) {
      Parameters.BondCompareParameters.MatchStereo = true;
    }
  }
  To = nanoClock();
}

static bool molPtr_NumBondLess(
    const ROMol* l,
    const ROMol* r) {  // need for sorting the source molecules by size
  return l->getNumBonds() < r->getNumBonds();
}

void MaximumCommonSubgraph::init() {
  QueryMolecule = Molecules.front();

  Targets.clear();
#ifdef FAST_SUBSTRUCT_CACHE
  QueryAtomLabels.clear();
  QueryBondLabels.clear();
  QueryAtomMatchTable.clear();
  QueryBondMatchTable.clear();
  RingMatchTables.clear();
#endif
#ifdef DUP_SUBSTRUCT_CACHE
  DuplicateCache.clear();
#endif
  void* userData = Parameters.CompareFunctionsUserData;

  size_t nq = 0;
#ifdef FAST_SUBSTRUCT_CACHE
  // fill out RingMatchTables to check cache Hash collision by checking match a
  // part of Query to Query
  if (!userData  // predefined functor - compute RingMatchTable for all targets
      && (Parameters.BondCompareParameters.CompleteRingsOnly ||
          Parameters.BondCompareParameters.RingMatchesRingOnly ||
          Parameters.AtomCompareParameters.RingMatchesRingOnly)) {
    RingMatchTables.init(QueryMolecule);
    Parameters.CompareFunctionsUserData = &RingMatchTables;
    RingMatchTables.computeRingMatchTable(QueryMolecule, QueryMolecule,
                                          Parameters);
  }

  // fill out match tables
  nq = QueryMolecule->getNumAtoms();
  QueryAtomMatchTable.resize(nq, nq);
  for (size_t aj = 0; aj < nq; aj++) {
    for (size_t ai = 0; ai < nq; ai++) {
      QueryAtomMatchTable.set(
          ai, aj,
          Parameters.AtomTyper(Parameters.AtomCompareParameters, *QueryMolecule,
                               ai, *QueryMolecule, aj,
                               Parameters.CompareFunctionsUserData));
    }
  }
  nq = QueryMolecule->getNumBonds();
  QueryBondMatchTable.resize(nq, nq);
  for (size_t aj = 0; aj < nq; aj++) {
    for (size_t ai = 0; ai < nq; ai++) {
      QueryBondMatchTable.set(
          ai, aj,
          Parameters.BondTyper(Parameters.BondCompareParameters, *QueryMolecule,
                               ai, *QueryMolecule, aj,
                               Parameters.CompareFunctionsUserData));
    }
  }
  // Compute label values based on current functor and parameters for code
  // Morgan correct computation.
  unsigned currentLabelValue = 1;
  std::vector<LabelDefinition> labels;
  nq = QueryMolecule->getNumAtoms();
  QueryAtomLabels.resize(nq);
  for (size_t ai = 0; ai < nq; ai++) {
    if (MCSAtomCompareAny ==
        Parameters.AtomTyper) {  // predefined functor without atom compare
                                 // parameters
      QueryAtomLabels[ai] = 1;
    } else {
      const Atom* atom = QueryMolecule->getAtomWithIdx(ai);
      if (MCSAtomCompareElements ==
          Parameters.AtomTyper) {  // predefined functor without atom compare
                                   // parameters
        QueryAtomLabels[ai] = atom->getAtomicNum() |
                              (Parameters.AtomCompareParameters.MatchValences
                                   ? (atom->getTotalValence() >> 8)
                                   : 0);
      } else if (MCSAtomCompareIsotopes ==
                 Parameters.AtomTyper) {  // predefined functor without atom
                                          // compare parameters
        QueryAtomLabels[ai] = atom->getAtomicNum() | (atom->getIsotope() >> 8) |
                              (Parameters.AtomCompareParameters.MatchValences
                                   ? (atom->getTotalValence() >> 16)
                                   : 0);
      } else {  // custom user defined functor
        QueryAtomLabels[ai] = NotSet;
        for (auto& label : labels) {
          if (Parameters.AtomTyper(Parameters.AtomCompareParameters,
                                   *QueryMolecule, label.ItemIndex,
                                   *QueryMolecule, ai,
                                   userData)) {  // equal items
            QueryAtomLabels[ai] = label.Value;
            break;
          }
        }
        if (NotSet == QueryAtomLabels[ai]) {  // not found -> create new label
          QueryAtomLabels[ai] = ++currentLabelValue;
          labels.emplace_back(ai, currentLabelValue);
        }
      }
    }
  }
  labels.clear();
  currentLabelValue = 1;
  nq = QueryMolecule->getNumBonds();
  QueryBondLabels.resize(nq);
  for (size_t aj = 0; aj < nq; aj++) {
    const Bond* bond = QueryMolecule->getBondWithIdx(aj);
    unsigned ring = 0;
    if (!userData && (Parameters.BondCompareParameters.CompleteRingsOnly ||
                      Parameters.BondCompareParameters.RingMatchesRingOnly)) {
      ring = RingMatchTables.isQueryBondInRing(aj) ? 0 : 1;  // is bond in ring
    }
    if (MCSBondCompareAny ==
        Parameters.BondTyper) {  // predefined functor without atom compare
                                 // parameters
      QueryBondLabels[aj] = 1 | (ring >> 8);
    } else if (MCSBondCompareOrderExact ==
               Parameters.BondTyper) {  // predefined functor without compare
                                        // parameters
      QueryBondLabels[aj] = (bond->getBondType() + 1) | (ring >> 8);
    } else if (MCSBondCompareOrder ==
               Parameters
                   .BondTyper) {  // predefined functor, ignore Aromatization
      unsigned order = bond->getBondType();
      if (Bond::AROMATIC == order ||
          Bond::ONEANDAHALF == order) {  // ignore Aromatization
        order = Bond::SINGLE;
      } else if (Bond::TWOANDAHALF == order) {
        order = Bond::DOUBLE;
      } else if (Bond::THREEANDAHALF == order) {
        order = Bond::TRIPLE;
      } else if (Bond::FOURANDAHALF == order) {
        order = Bond::QUADRUPLE;
      } else if (Bond::FIVEANDAHALF == order) {
        order = Bond::QUINTUPLE;
      }
      QueryBondLabels[aj] = (order + 1) | (ring >> 8);
    } else {  // custom user defined functor
      QueryBondLabels[aj] = NotSet;
      for (auto& label : labels) {
        if (Parameters.BondTyper(Parameters.BondCompareParameters,
                                 *QueryMolecule, label.ItemIndex,
                                 *QueryMolecule, aj,
                                 userData)) {  // equal bonds + ring ...
          QueryBondLabels[aj] = label.Value;
          break;
        }
      }
      if (NotSet == QueryAtomLabels[aj]) {  // not found -> create new label
        QueryBondLabels[aj] = ++currentLabelValue;
        labels.emplace_back(aj, currentLabelValue);
      }
    }
  }
#endif
  Targets.resize(Molecules.size() - 1);
  size_t i = 0;
  for (auto it = Molecules.begin() + 1; it != Molecules.end(); it++, i++) {
    Targets[i].Molecule = *it;
    // build Target Topology ADD ATOMs
    size_t j = 0;  // current item
    for (ROMol::ConstAtomIterator a = Targets[i].Molecule->beginAtoms();
         a != Targets[i].Molecule->endAtoms(); a++, j++) {
      Targets[i].Topology.addAtom((*a)->getIdx());
    }
    // build Target Topology ADD BONDs
    for (ROMol::ConstBondIterator b = Targets[i].Molecule->beginBonds();
         b != Targets[i].Molecule->endBonds(); b++) {
      const Bond* bond = *b;
      unsigned ii = bond->getBeginAtomIdx();
      unsigned jj = bond->getEndAtomIdx();
      Targets[i].Topology.addBond((*b)->getIdx(), ii, jj);
    }

#ifdef FAST_SUBSTRUCT_CACHE
    // fill out RingMatchTables
    if (!userData  // predefined functor - compute RingMatchTable for all
                   // targets
        && (Parameters.BondCompareParameters.CompleteRingsOnly ||
            Parameters.BondCompareParameters.RingMatchesRingOnly)) {
      RingMatchTables.addTargetBondRingsIndeces(Targets[i].Molecule);
      RingMatchTables.computeRingMatchTable(QueryMolecule, Targets[i].Molecule,
                                            Parameters);
    }
#endif

    // fill out match tables
    size_t nq = QueryMolecule->getNumAtoms();
    size_t nt = (*it)->getNumAtoms();
    Targets[i].AtomMatchTable.resize(nq, nt);

    for (size_t aj = 0; aj < nt; aj++) {
      for (size_t ai = 0; ai < nq; ai++) {
        Targets[i].AtomMatchTable.set(
            ai, aj,
            Parameters.AtomTyper(Parameters.AtomCompareParameters,
                                 *QueryMolecule, ai, *Targets[i].Molecule, aj,
                                 Parameters.CompareFunctionsUserData));
      }
    }
    nq = QueryMolecule->getNumBonds();
    nt = (*it)->getNumBonds();
    Targets[i].BondMatchTable.resize(nq, nt);
    for (size_t aj = 0; aj < nt; aj++) {
      for (size_t ai = 0; ai < nq; ai++) {
        Targets[i].BondMatchTable.set(
            ai, aj,
            Parameters.BondTyper(Parameters.BondCompareParameters,
                                 *QueryMolecule, ai, *Targets[i].Molecule, aj,
                                 Parameters.CompareFunctionsUserData));
      }
    }
  }
  Parameters.CompareFunctionsUserData = userData;  // restore
}

struct QueryRings {
  std::vector<unsigned> BondRings;  // amount of rings
  std::vector<unsigned> AtomRings;  // amount of rings

  QueryRings(const ROMol* query)
      : BondRings(query->getNumBonds()), AtomRings(query->getNumAtoms()) {
    for (unsigned int& BondRing : BondRings) {
      BondRing = 0;
    }
    const RingInfo::VECT_INT_VECT& brings = query->getRingInfo()->bondRings();
    for (const auto& ring : brings) {
      for (int ri : ring) {
        ++BondRings[ri];
      }
    }
    for (unsigned int& AtomRing : AtomRings) {
      AtomRing = 0;
    }
    const RingInfo::VECT_INT_VECT& arings = query->getRingInfo()->atomRings();
    for (const auto& ring : arings) {
      for (int ri : ring) {
        ++AtomRings[ri];
      }
    }
  }

  inline unsigned getNumberRings(const Bond* bond) const {
    return BondRings[bond->getIdx()];
  }

  inline unsigned getNumberRings(const Atom* atom) const {
    return AtomRings[atom->getIdx()];
  }
};  // namespace RDKit

struct WeightedBond {
  const Bond* BondPtr{nullptr};
  unsigned Weight{0};
  WeightedBond()  {}
  WeightedBond(const Bond* bond, const QueryRings& r)
      : BondPtr(bond), Weight(0) {
    // score ((bond.is_in_ring + atom1.is_in_ring + atom2.is_in_ring)
    if (r.getNumberRings(bond)) {
      Weight += 1;
    }
    if (r.getNumberRings(bond->getBeginAtom())) {
      Weight += 1;
    }
    if (r.getNumberRings(bond->getEndAtom())) {
      Weight += 1;
    }
  }
  bool operator<(const WeightedBond& r) {
    return Weight >= r.Weight;  // sort in Z-A order (Rings first)
  }
};

void MaximumCommonSubgraph::makeInitialSeeds() {
  // build a set of initial seeds as "all" single bonds from query
  // molecule
  std::vector<bool> excludedBonds(QueryMolecule->getNumBonds(), false);

  Seeds.clear();
  QueryMoleculeMatchedBonds = 0;
  QueryMoleculeMatchedAtoms = 0;
  QueryMoleculeSingleMatchedAtom = nullptr;
  if (!Parameters.InitialSeed.empty()) {  // make user defined seed
    std::unique_ptr<const ROMol> initialSeedMolecule(
        (const ROMol*)SmartsToMol(Parameters.InitialSeed));
    // make a set of of seed as indices and pointers to current query
    // molecule items based on matching results
    std::vector<MatchVectType> matching_substructs;
    SubstructMatch(*QueryMolecule, *initialSeedMolecule, matching_substructs);
    // loop throw all fragments of Query matched to initial seed
    for (std::vector<MatchVectType>::const_iterator ms =
             matching_substructs.begin();
         ms != matching_substructs.end(); ms++) {
      Seed seed;
      seed.ExcludedBonds = excludedBonds;
      seed.MatchResult.resize(Targets.size());
#ifdef VERBOSE_STATISTICS_ON
      {
        ++VerboseStatistics.Seed;
        ++VerboseStatistics.InitialSeed;
      }
#endif
      // add all matched atoms of the matched query fragment
      std::map<unsigned, unsigned> initialSeedToQueryAtom;
      for (const auto& msb : *ms) {
        unsigned qai = msb.second;
        unsigned sai = msb.first;
        seed.addAtom(QueryMolecule->getAtomWithIdx(qai));
        initialSeedToQueryAtom[sai] = qai;
      }
      // add all bonds (existed in initial seed !!!) between all matched
      // atoms in query
      for (const auto& msb : *ms) {
        const Atom* atom = initialSeedMolecule->getAtomWithIdx(msb.first);
	for(auto *initialBond: atom->bonds()) {
          unsigned qai1 =
              initialSeedToQueryAtom.find(initialBond->getBeginAtomIdx())
                  ->second;
          unsigned qai2 =
              initialSeedToQueryAtom.find(initialBond->getEndAtomIdx())->second;

          const Bond* b = QueryMolecule->getBondBetweenAtoms(qai1, qai2);
          if (!seed.ExcludedBonds[b->getIdx()]) {
            seed.addBond(b);
            seed.ExcludedBonds[b->getIdx()] = true;
          }
        }
      }
      seed.computeRemainingSize(*QueryMolecule);

      if (checkIfMatchAndAppend(seed)) {
        QueryMoleculeMatchedBonds = seed.getNumBonds();
      }
    }
  } else {  // create a set of seeds from each query bond
    // R1 additional performance OPTIMISATION
    // if(Parameters.BondCompareParameters.CompleteRingsOnly)
    // disable all mismatched rings, and do not generate initial seeds
    // from such disabled bonds
    //  for(  rings .....) for(i......)
    //   if(mismatched) excludedBonds[i.......] = true;
    QueryRings r(QueryMolecule);
    std::vector<WeightedBond> wb;
    wb.reserve(QueryMolecule->getNumBonds());
    for (RWMol::ConstBondIterator bi = QueryMolecule->beginBonds();
         bi != QueryMolecule->endBonds(); bi++) {
      wb.emplace_back(*bi, r);
    }

    for (std::vector<WeightedBond>::const_iterator bi = wb.begin();
         bi != wb.end(); bi++) {
      // R1 additional performance OPTIMISATION
      // if(excludedBonds[(*bi)->getIdx()])
      //    continue;
      Seed seed;
      seed.MatchResult.resize(Targets.size());

#ifdef VERBOSE_STATISTICS_ON
      {
        ++VerboseStatistics.Seed;
        ++VerboseStatistics.InitialSeed;
      }
#endif
      seed.addAtom(bi->BondPtr->getBeginAtom());
      seed.addAtom(bi->BondPtr->getEndAtom());
      seed.ExcludedBonds = excludedBonds;  // all bonds from first to current
      seed.addBond(bi->BondPtr);
      excludedBonds[bi->BondPtr->getIdx()] = true;

      seed.computeRemainingSize(*QueryMolecule);

      if (checkIfMatchAndAppend(seed)) {
        ++QueryMoleculeMatchedBonds;
      } else {
        // optionally remove all such bonds from all targets TOPOLOGY
        // where it exists.
        //..........

        // disable (mark as already processed) mismatched bond in all
        // seeds
        for (auto& Seed : Seeds) {
          Seed.ExcludedBonds[bi->BondPtr->getIdx()] = true;
        }

#ifdef VERBOSE_STATISTICS_ON
        ++VerboseStatistics.MismatchedInitialSeed;
#endif
      }
    }
  }
  size_t nq = QueryMolecule->getNumAtoms();
  for (size_t i = 0; i < nq; i++) {  // all query's atoms
    const Atom* queryMolAtom = QueryMolecule->getAtomWithIdx(i);
    bool isQueryMolAtomInRing = queryIsAtomInRing(queryMolAtom);
    unsigned matched = 0;
    for (std::vector<Target>::const_iterator tag = Targets.begin();
         tag != Targets.end(); tag++) {
      size_t nt = tag->Molecule->getNumAtoms();
      for (size_t aj = 0; aj < nt; aj++) {
        if (tag->AtomMatchTable.at(i, aj)) {
          const Atom* targetMolAtom = tag->Molecule->getAtomWithIdx(aj);
          bool isTargetMolAtomInRing = queryIsAtomInRing(targetMolAtom);
          ++matched;
          if (!(Parameters.BondCompareParameters.CompleteRingsOnly &&
                (isQueryMolAtomInRing || isTargetMolAtomInRing))) {
            if (!QueryMoleculeSingleMatchedAtom) {
              QueryMoleculeSingleMatchedAtom = queryMolAtom;
            } else {
              QueryMoleculeSingleMatchedAtom = (std::max)(
                  queryMolAtom, QueryMoleculeSingleMatchedAtom,
                  [](const Atom* a, const Atom* b) {
                    if (a->getDegree() != b->getDegree()) {
                      return (a->getDegree() < b->getDegree());
                    } else if (a->getFormalCharge() != b->getFormalCharge()) {
                      return (a->getFormalCharge() < b->getFormalCharge());
                    } else if (a->getAtomicNum() != b->getAtomicNum()) {
                      return (a->getAtomicNum() < b->getAtomicNum());
                    }
                    return (a->getIdx() < b->getIdx());
                  });
            }
          }
          break;
        }
      }
    }
    if (matched >= ThresholdCount) {
      ++QueryMoleculeMatchedAtoms;
    }
  }
}

namespace {

void _DFS(const Graph& g, const boost::dynamic_bitset<>& ringBonds,
          boost::dynamic_bitset<>& openBonds, unsigned vertex,
          std::vector<unsigned> apath, std::vector<unsigned>& colors,
          unsigned fromVertex) {
  std::pair<Graph::out_edge_iterator, Graph::out_edge_iterator> bonds =
      boost::out_edges(vertex, g);
  colors[vertex] = 1;
  while (bonds.first != bonds.second) {
    unsigned bIdx = g[*(bonds.first)];
    if (!ringBonds[bIdx]) {
      ++bonds.first;
      continue;
    }
    // find the other vertex:
    unsigned oVertex = boost::source(*(bonds.first), g);
    if (oVertex == vertex) {
      oVertex = boost::target((*bonds.first), g);
    }
    ++bonds.first;

    if (oVertex == fromVertex) {
      continue;
    }
    if (colors[oVertex] == 1) {
      // closing a cycle
      for (auto ai = apath.rbegin(); ai != apath.rend() && *ai != oVertex;
           ++ai) {
        auto epair = boost::edge(*ai, *(ai + 1), g);
        CHECK_INVARIANT(epair.second, "edge not found");
        openBonds.reset(g[epair.first]);
      }
      auto epair = boost::edge(oVertex, *apath.rbegin(), g);
      CHECK_INVARIANT(epair.second, "edge not found");
      openBonds.reset(g[epair.first]);
    } else if (colors[oVertex] == 0) {
      std::vector<unsigned> napath = apath;
      napath.push_back(oVertex);
      _DFS(g, ringBonds, openBonds, oVertex, napath, colors, vertex);
    }
  }
  colors[vertex] = 2;
  return;
}

bool checkIfRingsAreClosed(const Seed& fs) {
  if (!fs.MoleculeFragment.Bonds.size()) {
    return true;
  }

  bool res = true;
  const auto& om = fs.MoleculeFragment.Bonds[0]->getOwningMol();
  const auto ri = om.getRingInfo();
  boost::dynamic_bitset<> ringBonds(om.getNumBonds());
  for (const auto bond : fs.MoleculeFragment.Bonds) {
    if (ri->numBondRings(bond->getIdx())) {
      ringBonds.set(bond->getIdx());
    }
  }
  boost::dynamic_bitset<> openBonds = ringBonds;
  if (openBonds.any()) {
    std::vector<unsigned> colors(om.getNumAtoms());
    for (unsigned bi = 0; bi < openBonds.size(); ++bi) {
      if (!openBonds[bi]) {
        continue;
      }
      std::pair<Graph::edge_iterator, Graph::edge_iterator> bonds =
          boost::edges(fs.Topology);
      while (bonds.first != bonds.second && fs.Topology[*(bonds.first)] != bi) {
        ++bonds.first;
      }
      CHECK_INVARIANT(bonds.first != bonds.second, "bond not found");
      unsigned startVertex = boost::source(*(bonds.first), fs.Topology);
      std::vector<unsigned> apath = {startVertex};
      _DFS(fs.Topology, ringBonds, openBonds, startVertex, apath, colors,
           om.getNumAtoms() + 1);
    }
    res = openBonds.none();
  }
  return res;
}

}  // namespace
bool MaximumCommonSubgraph::growSeeds() {
  bool mcsFound = false;
  bool canceled = false;
  unsigned steps = 99999;  // steps from last progress callback call. call
                           // it immediately in the beginning

  // Find MCS -- SDF Seed growing OPTIMISATION (it works in 3 times
  // faster)
  while (!Seeds.empty()) {
    if (getMaxNumberBonds() == QueryMoleculeMatchedBonds) {  // MCS == Query
      break;
    }
    ++steps;
#ifdef VERBOSE_STATISTICS_ON
    VerboseStatistics.TotalSteps++;
#endif
    auto si = Seeds.begin();

    si->grow(*this);
    {
      const Seed& fs = Seeds.front();
      // bigger substructure found
      if (fs.CopyComplete) {
        bool possibleMCS = false;
        if ((!Parameters.MaximizeBonds &&
             (fs.getNumAtoms() > getMaxNumberAtoms() ||
              (fs.getNumAtoms() == getMaxNumberAtoms() &&
               fs.getNumBonds() > getMaxNumberBonds())))) {
          possibleMCS = true;
        } else if (Parameters.MaximizeBonds &&
                   (fs.getNumBonds() > getMaxNumberBonds() ||
                    (fs.getNumBonds() == getMaxNumberBonds() &&
                     fs.getNumAtoms() > getMaxNumberAtoms()))) {
          possibleMCS = true;
        }
        // #945: test here to see if the MCS actually has all rings closed
        if (possibleMCS && Parameters.BondCompareParameters.CompleteRingsOnly) {
          possibleMCS = checkIfRingsAreClosed(fs);
        }
        if (possibleMCS) {
          mcsFound = true;
#ifdef VERBOSE_STATISTICS_ON
          VerboseStatistics.MCSFoundStep = VerboseStatistics.TotalSteps;
          VerboseStatistics.MCSFoundTime = nanoClock();
#endif
          McsIdx.Atoms = fs.MoleculeFragment.Atoms;
          McsIdx.Bonds = fs.MoleculeFragment.Bonds;
          McsIdx.AtomsIdx = fs.MoleculeFragment.AtomsIdx;
          McsIdx.BondsIdx = fs.MoleculeFragment.BondsIdx;
          if (Parameters.Verbose) {
            std::cout << VerboseStatistics.TotalSteps
                      << " Seeds:" << Seeds.size() << " MCS "
                      << McsIdx.Atoms.size() << " atoms, "
                      << McsIdx.Bonds.size() << " bonds";
            printf(" for %.4lf seconds. bond[0]=%u\n",
                   double(VerboseStatistics.MCSFoundTime - To) / 1000000.,
                   McsIdx.BondsIdx[0]);
          }
        }
      }
    }
    if (NotSet == si->GrowingStage) {  // finished
      Seeds.erase(si);
    }
    if (Parameters.ProgressCallback) {
      steps = 0;
      Stat.NumAtoms = getMaxNumberAtoms();
      Stat.NumBonds = getMaxNumberBonds();
      if (!Parameters.ProgressCallback(Stat, Parameters,
                                       Parameters.ProgressCallbackUserData)) {
        canceled = true;
        break;
      }
    }
  }

  if (mcsFound) {  // postponed copy of current set of molecules for
                   // threshold < 1.
    McsIdx.QueryMolecule = QueryMolecule;
    McsIdx.Targets = Targets;
  }
  return !canceled;
}  // namespace FMCS

struct AtomMatch {  // for each seed atom (matched)
  unsigned QueryAtomIdx;
  unsigned TargetAtomIdx;
  AtomMatch() : QueryAtomIdx(NotSet), TargetAtomIdx(NotSet) {}
};
typedef std::vector<AtomMatch> AtomMatchSet;

std::pair<std::string, RWMol*>
MaximumCommonSubgraph::generateResultSMARTSAndQueryMol(
    const MCS& mcsIdx) const {
  // match the result MCS with all targets to check if it is exact match
  // or template
  Seed seed;  // result MCS
  seed.ExcludedBonds.resize(mcsIdx.QueryMolecule->getNumBonds(), false);
  std::vector<AtomMatchSet> atomMatchResult(mcsIdx.Targets.size());
  std::vector<unsigned> atomIdxMap(mcsIdx.QueryMolecule->getNumAtoms());
  std::vector<std::map<unsigned, const Bond*>> bondMatchSet(
      mcsIdx.Bonds.size());  // key is unique BondType
  std::vector<std::map<unsigned, const Atom*>> atomMatchSet(
      mcsIdx.Atoms.size());  // key is unique atomic number

  for (auto atom : mcsIdx.Atoms) {
    atomIdxMap[atom->getIdx()] = seed.getNumAtoms();
    seed.addAtom(atom);
  }
  for (auto bond : mcsIdx.Bonds) {
    seed.addBond(bond);
  }

  if (!mcsIdx.Bonds.empty()) {
    unsigned itarget = 0;
    for (auto tag = mcsIdx.Targets.begin(); tag != mcsIdx.Targets.end();
         tag++, itarget++) {
      match_V_t match;  // THERE IS NO Bonds match INFO !!!!
      bool target_matched = SubstructMatchCustomTable(
          tag->Topology, *tag->Molecule, seed.Topology, *QueryMolecule,
          tag->AtomMatchTable, tag->BondMatchTable, &Parameters, &match);
      if (!target_matched) {
        continue;
      }
      atomMatchResult[itarget].resize(seed.getNumAtoms());
      for (match_V_t::const_iterator mit = match.begin(); mit != match.end();
           mit++) {
        unsigned ai = mit->first;  // SeedAtomIdx
        atomMatchResult[itarget][ai].QueryAtomIdx = seed.Topology[mit->first];
        atomMatchResult[itarget][ai].TargetAtomIdx = tag->Topology[mit->second];
        const Atom* ta =
            tag->Molecule->getAtomWithIdx(tag->Topology[mit->second]);
        if (ta && ta->getAtomicNum() !=
                      seed.MoleculeFragment.Atoms[ai]->getAtomicNum()) {
          atomMatchSet[ai][ta->getAtomicNum()] = ta;  // add
        }
      }
      // AND BUILD BOND MATCH INFO
      unsigned bi = 0;
      for (auto bond = mcsIdx.Bonds.begin(); bond != mcsIdx.Bonds.end();
           bond++, bi++) {
        unsigned i = atomIdxMap[(*bond)->getBeginAtomIdx()];
        unsigned j = atomIdxMap[(*bond)->getEndAtomIdx()];
        unsigned ti = atomMatchResult[itarget][i].TargetAtomIdx;
        unsigned tj = atomMatchResult[itarget][j].TargetAtomIdx;
        const Bond* tb = tag->Molecule->getBondBetweenAtoms(ti, tj);
        if (tb && (*bond)->getBondType() != tb->getBondType()) {
          bondMatchSet[bi][tb->getBondType()] = tb;  // add
        }
      }
    }
  }

  // Generate result's SMARTS

  // create molecule from MCS for MolToSmarts()
  auto* mol = new RWMol();
  const RingInfo* ri = mcsIdx.QueryMolecule->getRingInfo();
  unsigned ai = 0;  // SeedAtomIdx
  for (auto atom = mcsIdx.Atoms.begin(); atom != mcsIdx.Atoms.end();
       atom++, ai++) {
    QueryAtom a;
    if (Parameters.AtomTyper == MCSAtomCompareIsotopes ||
        Parameters.AtomCompareParameters
            .MatchIsotope) {  // do '[0*]-[0*]-[13*]' for CC[13NH2]
      a.setQuery(makeAtomIsotopeQuery((int)(*atom)->getIsotope()));
    } else {
      // generate [#6] instead of C or c !
      a.setQuery(makeAtomNumQuery((*atom)->getAtomicNum()));
      // for all atomMatchSet[ai] items add atom query to template like
      // [#6,#17,#9, ... ]
      for (std::map<unsigned, const Atom*>::const_iterator am =
               atomMatchSet[ai].begin();
           am != atomMatchSet[ai].end(); am++) {
        a.expandQuery(makeAtomNumQuery(am->second->getAtomicNum()),
                      Queries::COMPOSITE_OR);
        if (Parameters.AtomCompareParameters.MatchChiralTag &&
            (am->second->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
             am->second->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW)) {
          a.setChiralTag(am->second->getChiralTag());
        }
      }
    }
    if (Parameters.AtomCompareParameters.RingMatchesRingOnly) {
      ATOM_EQUALS_QUERY* q = makeAtomInRingQuery();
      q->setNegation(!ri->numAtomRings((*atom)->getIdx()));
      a.expandQuery(q, Queries::COMPOSITE_AND, true);
    }
    mol->addAtom(&a, true, false);
  }
  unsigned bi = 0;  // Seed Idx
  for (auto bond = mcsIdx.Bonds.begin(); bond != mcsIdx.Bonds.end();
       bond++, bi++) {
    QueryBond b;
    unsigned beginAtomIdx = atomIdxMap[(*bond)->getBeginAtomIdx()];
    unsigned endAtomIdx = atomIdxMap[(*bond)->getEndAtomIdx()];
    b.setBeginAtomIdx(beginAtomIdx);
    b.setEndAtomIdx(endAtomIdx);
    b.setQuery(makeBondOrderEqualsQuery((*bond)->getBondType()));
    // add OR template if need
    for (std::map<unsigned, const Bond*>::const_iterator bm =
             bondMatchSet[bi].begin();
         bm != bondMatchSet[bi].end(); bm++) {
      b.expandQuery(makeBondOrderEqualsQuery(bm->second->getBondType()),
                    Queries::COMPOSITE_OR);
      if (Parameters.BondCompareParameters.MatchStereo &&
          bm->second->getStereo() > Bond::STEREOANY) {
        b.setStereo(bm->second->getStereo());
      }
    }
    if (Parameters.BondCompareParameters.RingMatchesRingOnly) {
      BOND_EQUALS_QUERY* q = makeBondIsInRingQuery();
      q->setNegation(!ri->numBondRings((*bond)->getIdx()));
      b.expandQuery(q, Queries::COMPOSITE_AND, true);
    }
    mol->addBond(&b, false);
  }

  return std::make_pair(MolToSmarts(*mol, true), mol);
}

bool MaximumCommonSubgraph::addFusedBondQueries(const MCS& mcsIdx,
                                                RWMol* rwMol) const {
  const RingInfo* ri = mcsIdx.QueryMolecule->getRingInfo();
  unsigned bi = 0;  // Seed Idx
  bool haveFusedBondQuery = false;
  std::map<int, std::set<size_t>> bondRingMembership;
  const VECT_INT_VECT& br = ri->bondRings();
  std::vector<size_t> mcsRingSize(br.size(), 0);
  for (size_t ringIdx = 0; ringIdx < br.size(); ++ringIdx) {
    for (int bondIdx : br[ringIdx]) {
      bondRingMembership[bondIdx].insert(ringIdx);
    }
  }
  for (auto Bond : mcsIdx.Bonds) {
    unsigned int bondIdx = Bond->getIdx();
    if (!ri->numBondRings(bondIdx)) {
      continue;
    }
    for (size_t ringIdx : bondRingMembership[bondIdx]) {
      ++mcsRingSize[ringIdx];
    }
  }
  for (auto bond = mcsIdx.Bonds.begin(); bond != mcsIdx.Bonds.end();
       ++bond, ++bi) {
    unsigned int bondIdx = (*bond)->getIdx();
    if (!ri->numBondRings(bondIdx)) {
      continue;
    }
    haveFusedBondQuery = true;
    for (int ringIdx : bondRingMembership[bondIdx]) {
      if (mcsRingSize[ringIdx] < br[ringIdx].size()) {
        continue;
      }
      rwMol->getBondWithIdx(bi)->expandQuery(
          makeBondInRingOfSizeQuery(br[ringIdx].size()), Queries::COMPOSITE_AND,
          true);
    }
  }
  return haveFusedBondQuery;
}

bool MaximumCommonSubgraph::createSeedFromMCS(size_t newQueryTarget,
                                              Seed& newSeed) {
  Seed mcs;
  mcs.ExcludedBonds.resize(McsIdx.QueryMolecule->getNumBonds(), false);
  std::vector<unsigned> mcsAtomIdxMap(McsIdx.QueryMolecule->getNumAtoms());

  for (std::vector<const Atom*>::const_iterator atom = McsIdx.Atoms.begin();
       atom != McsIdx.Atoms.end(); atom++) {
    mcsAtomIdxMap[(*atom)->getIdx()] = mcs.addAtom((*atom));
  }
  for (std::vector<const Bond*>::const_iterator bond = McsIdx.Bonds.begin();
       bond != McsIdx.Bonds.end(); bond++) {
    mcs.addBond((*bond));
  }

  const Target& newQuery = McsIdx.Targets[newQueryTarget];

  match_V_t match;
  bool target_matched = SubstructMatchCustomTable(
      newQuery.Topology, *newQuery.Molecule, mcs.Topology,
      *McsIdx.QueryMolecule, newQuery.AtomMatchTable, newQuery.BondMatchTable,
      &Parameters, &match);
  if (!target_matched) {
    return false;
  }

  AtomMatchSet atomMatchResult(mcs.getNumAtoms());

  newSeed.ExcludedBonds.resize(newQuery.Molecule->getNumBonds(), false);

  for (match_V_t::const_iterator mit = match.begin(); mit != match.end();
       mit++) {
    unsigned ai = mit->first;  // SeedAtomIdx in mcs seed
    atomMatchResult[ai].QueryAtomIdx = mcs.Topology[mit->first];
    atomMatchResult[ai].TargetAtomIdx = newQuery.Topology[mit->second];
    const Atom* ta =
        newQuery.Molecule->getAtomWithIdx(newQuery.Topology[mit->second]);
    newSeed.addAtom(ta);
  }

  for (std::vector<const Bond*>::const_iterator bond = McsIdx.Bonds.begin();
       bond != McsIdx.Bonds.end(); bond++) {
    unsigned i = mcsAtomIdxMap[(*bond)->getBeginAtomIdx()];
    unsigned j = mcsAtomIdxMap[(*bond)->getEndAtomIdx()];
    unsigned ti = atomMatchResult[i].TargetAtomIdx;
    unsigned tj = atomMatchResult[j].TargetAtomIdx;
    const Bond* tb = newQuery.Molecule->getBondBetweenAtoms(ti, tj);
    newSeed.addBond(tb);
  }
  newSeed.computeRemainingSize(*newQuery.Molecule);
  return true;
}

MCSResult MaximumCommonSubgraph::find(const std::vector<ROMOL_SPTR>& src_mols) {
  clear();
  MCSResult res;

  if (src_mols.size() < 2) {
    throw std::runtime_error(
        "FMCS. Invalid argument. mols.size() must be at least 2");
  }
  if (Parameters.Threshold > 1.0) {
    throw std::runtime_error(
        "FMCS. Invalid argument. Parameter Threshold must be 1.0 or "
        "less.");
  }

  ThresholdCount = (unsigned)ceil((src_mols.size()) * Parameters.Threshold) -
                   1;        // minimal required number of matched targets
  if (ThresholdCount < 1) {  // at least one target
    ThresholdCount = 1;
  }
  if (ThresholdCount > src_mols.size() - 1) {  // max all targets
    ThresholdCount = src_mols.size() - 1;
  }

  // Selecting CompleteRingsOnly option also enables
  // --ring-matches-ring-only. ring--ring and chain bonds only match chain
  // bonds.
  if (Parameters.BondCompareParameters.CompleteRingsOnly) {
    Parameters.BondCompareParameters.RingMatchesRingOnly = true;
  }

  for (const auto& src_mol : src_mols) {
    Molecules.push_back(src_mol.get());
    if (!Molecules.back()->getRingInfo()->isInitialized()) {
      Molecules.back()->getRingInfo()->initialize();  // but do not fill out !!!
    }
  }

  // sort source set of molecules by their 'size' and assume the smallest
  // molecule as a query
  std::stable_sort(Molecules.begin(), Molecules.end(), molPtr_NumBondLess);

  bool areSeedsEmpty = false;
  for (size_t i = 0;
       i < Molecules.size() - ThresholdCount && !areSeedsEmpty && !res.Canceled;
       ++i) {
    init();
    if (Targets.empty()) {
      break;
    }
    MCSFinalMatchCheckFunction tff = Parameters.FinalMatchChecker;
    if (FinalMatchCheckFunction == Parameters.FinalMatchChecker) {
      Parameters.FinalMatchChecker = nullptr;  // skip final match check for
    }
    // initial seed to allow future growing
    // of it
    makeInitialSeeds();
    Parameters.FinalMatchChecker = tff;  // restore final functor

    if (Parameters.Verbose) {
      std::cout << "Query " << MolToSmiles(*QueryMolecule) << " "
                << QueryMolecule->getNumAtoms() << "("
                << QueryMoleculeMatchedAtoms << ") atoms, "
                << QueryMolecule->getNumBonds() << "("
                << QueryMoleculeMatchedBonds << ") bonds\n";
    }

    areSeedsEmpty = Seeds.empty();
    res.Canceled = areSeedsEmpty || growSeeds() ? false : true;
    // verify what MCS is equal to one of initial seed for chirality match
    if ((FinalMatchCheckFunction == Parameters.FinalMatchChecker &&
         1 == getMaxNumberBonds()) ||
        0 == getMaxNumberBonds()) {
      McsIdx = MCS();      // clear
      makeInitialSeeds();  // check all possible initial seeds
      if (!areSeedsEmpty) {
        const Seed& fs = Seeds.front();
        if (1 == getMaxNumberBonds() ||
            !(Parameters.BondCompareParameters.CompleteRingsOnly &&
              fs.MoleculeFragment.Bonds.size() == 1 &&
              queryIsBondInRing(fs.MoleculeFragment.Bonds.front()))) {
          McsIdx.QueryMolecule = QueryMolecule;
          McsIdx.Atoms = fs.MoleculeFragment.Atoms;
          McsIdx.Bonds = fs.MoleculeFragment.Bonds;
          McsIdx.AtomsIdx = fs.MoleculeFragment.AtomsIdx;
          McsIdx.BondsIdx = fs.MoleculeFragment.BondsIdx;
        }
      }
      if (!McsIdx.QueryMolecule && QueryMoleculeSingleMatchedAtom) {
        McsIdx.QueryMolecule = QueryMolecule;
        McsIdx.Atoms = std::vector<const Atom*>{QueryMoleculeSingleMatchedAtom};
        McsIdx.Bonds = std::vector<const Bond*>();
        McsIdx.AtomsIdx = std::vector<unsigned>{0};
        McsIdx.BondsIdx = std::vector<unsigned>();
      }

    } else if (i + 1 < Molecules.size() - ThresholdCount) {
      Seed seed;
      if (createSeedFromMCS(i, seed)) {  // MCS is matched with new query
        Seeds.push_back(seed);
      }
      std::swap(Molecules[0],
                Molecules[i + 1]);  // change query molecule for threshold < 1.
    }
  }

  res.NumAtoms = getMaxNumberAtoms();
  if (!res.NumAtoms && QueryMoleculeSingleMatchedAtom) {
    res.NumAtoms = 1;
  }
  res.NumBonds = getMaxNumberBonds();

  if (res.NumBonds > 0 || QueryMoleculeSingleMatchedAtom) {
    std::pair<std::string, RWMol*> smartsQueryMolPair =
        generateResultSMARTSAndQueryMol(McsIdx);
    res.SmartsString = smartsQueryMolPair.first;
    res.QueryMol = ROMOL_SPTR(smartsQueryMolPair.second);
    if (Parameters.BondCompareParameters.MatchFusedRingsStrict) {
      addFusedBondQueries(McsIdx, smartsQueryMolPair.second);
    }
  }

#ifdef VERBOSE_STATISTICS_ON
  if (Parameters.Verbose) {
    unsigned itarget = 0;
    for (std::vector<Target>::const_iterator tag = Targets.begin();
         res.NumAtoms > 0 && tag != Targets.end(); tag++, itarget++) {
      MatchVectType match;

      bool target_matched =
          SubstructMatch(*tag->Molecule, *res.QueryMol, match);
      if (!target_matched) {
        std::cout << "Target " << itarget + 1
                  << (target_matched ? " matched " : " MISMATCHED ")
                  << MolToSmiles(*tag->Molecule) << "\n";
      }
    }

    std::cout << "STATISTICS:\n";
    std::cout << "Total Growing Steps  = " << VerboseStatistics.TotalSteps
              << ", MCS found on " << VerboseStatistics.MCSFoundStep << " step";
    if (VerboseStatistics.MCSFoundTime - To > 0) {
      printf(", for %.4lf seconds\n",
             double(VerboseStatistics.MCSFoundTime - To) / 1000000.);
    } else {
      std::cout << ", for less than 1 second\n";
    }
    std::cout << "Initial   Seeds      = " << VerboseStatistics.InitialSeed
              << ",  Mismatched " << VerboseStatistics.MismatchedInitialSeed
              << "\n";
    std::cout << "Inspected Seeds      = " << VerboseStatistics.Seed << "\n";
    std::cout << "Rejected by BestSize = "
              << VerboseStatistics.RemainingSizeRejected << "\n";
    std::cout << "SingleBondExcluded   = "
              << VerboseStatistics.SingleBondExcluded << "\n";
#ifdef EXCLUDE_WRONG_COMPOSITION
    std::cout << "Rejected by WrongComposition = "
              << VerboseStatistics.WrongCompositionRejected << " [ "
              << VerboseStatistics.WrongCompositionDetected << " Detected ]\n";
#endif
    std::cout << "MatchCheck Seeds     = " << VerboseStatistics.SeedCheck
              << "\n";
    std::cout  //<< "\n"
        << "     MatchCalls = " << VerboseStatistics.MatchCall << "\n"
        << "     MatchFound = " << VerboseStatistics.MatchCallTrue << "\n";
    std::cout << " fastMatchCalls = " << VerboseStatistics.FastMatchCall << "\n"
              << " fastMatchFound = " << VerboseStatistics.FastMatchCallTrue
              << "\n";
    std::cout << " slowMatchCalls = "
              << VerboseStatistics.MatchCall -
                     VerboseStatistics.FastMatchCallTrue
              << "\n"
              << " slowMatchFound = " << VerboseStatistics.SlowMatchCallTrue
              << "\n";

#ifdef VERBOSE_STATISTICS_FASTCALLS_ON
    std::cout << "AtomFunctorCalls = " << VerboseStatistics.AtomFunctorCalls
              << "\n";
    std::cout << "BondCompareCalls = " << VerboseStatistics.BondCompareCalls
              << "\n";
#endif
    std::cout << "  DupCacheFound = " << VerboseStatistics.DupCacheFound
              << "   " << VerboseStatistics.DupCacheFoundMatch << " matched, "
              << VerboseStatistics.DupCacheFound -
                     VerboseStatistics.DupCacheFoundMatch
              << " mismatched\n";
#ifdef FAST_SUBSTRUCT_CACHE
    std::cout << "HashCache size  = " << HashCache.keyssize() << " keys\n";
    std::cout << "HashCache size  = " << HashCache.fullsize() << " entries\n";
    std::cout << "FindHashInCache = " << VerboseStatistics.FindHashInCache
              << "\n";
    std::cout << "HashFoundInCache= " << VerboseStatistics.HashKeyFoundInCache
              << "\n";
    std::cout << "ExactMatchCalls = " << VerboseStatistics.ExactMatchCall
              << "\n"
              << "ExactMatchFound = " << VerboseStatistics.ExactMatchCallTrue
              << "\n";
#endif
  }
#endif

  clear();
  return res;
}

bool MaximumCommonSubgraph::checkIfMatchAndAppend(Seed& seed) {
#ifdef VERBOSE_STATISTICS_ON
  ++VerboseStatistics.SeedCheck;
#endif
#ifdef FAST_SUBSTRUCT_CACHE
  SubstructureCache::HashKey cacheKey;
  SubstructureCache::TIndexEntry* cacheEntry = nullptr;
  bool cacheEntryIsValid = false;
  RDUNUSED_PARAM(cacheEntryIsValid);  // unused var
#endif

  bool foundInCache = false;
  bool foundInDupCache = false;

  {
#ifdef DUP_SUBSTRUCT_CACHE
    if (DuplicateCache.find(seed.DupCacheKey, foundInCache)) {
// duplicate found. skip match() but store both seeds, because they will grow by
// different paths !!!
#ifdef VERBOSE_STATISTICS_ON
      VerboseStatistics.DupCacheFound++;
      VerboseStatistics.DupCacheFoundMatch += foundInCache ? 1 : 0;
#endif
      if (!foundInCache) {  // mismatched !!!
        return false;
      }
    }
    foundInDupCache = foundInCache;
#endif
#ifdef FAST_SUBSTRUCT_CACHE
    if (!foundInCache) {
#ifdef VERBOSE_STATISTICS_ON
      ++VerboseStatistics.FindHashInCache;
#endif
      cacheEntry =
          HashCache.find(seed, QueryAtomLabels, QueryBondLabels, cacheKey);
      cacheEntryIsValid = true;
      if (cacheEntry) {  // possibly found. check for hash collision
#ifdef VERBOSE_STATISTICS_ON
        ++VerboseStatistics.HashKeyFoundInCache;
#endif
        // check hash collisions (time +3%):
        for (SubstructureCache::TIndexEntry::const_iterator g =
                 cacheEntry->begin();
             !foundInCache && g != cacheEntry->end(); g++) {
          if (g->m_vertices.size() != seed.getNumAtoms() ||
              g->m_edges.size() != seed.getNumBonds()) {
            continue;
          }
#ifdef VERBOSE_STATISTICS_ON
          ++VerboseStatistics.ExactMatchCall;
#endif
          // EXACT MATCH
          foundInCache = SubstructMatchCustomTable(
              (*g), *QueryMolecule, seed.Topology, *QueryMolecule,
              QueryAtomMatchTable, QueryBondMatchTable, &Parameters);
#ifdef VERBOSE_STATISTICS_ON
          if (foundInCache) {
            ++VerboseStatistics.ExactMatchCallTrue;
          }
#endif
        }
      }
    }
#endif
  }
  bool found = foundInCache;

  if (!found) {
    found = match(seed);
  }

  Seed* newSeed = nullptr;

  {
    if (found) {  // Store new generated seed, if found in cache or in
                  // all(- threshold) targets
      {
        newSeed = &Seeds.add(seed);
        newSeed->CopyComplete = false;
      }

#ifdef DUP_SUBSTRUCT_CACHE
      if (!foundInDupCache &&
          seed.getNumBonds() >= 3) {  // only seed with a ring
                                      // can be duplicated -
                                      // do not store very
                                      // small seed in cache
        DuplicateCache.add(seed.DupCacheKey, true);
      }
#endif
#ifdef FAST_SUBSTRUCT_CACHE
      if (!foundInCache) {
        HashCache.add(seed, cacheKey, cacheEntry);
      }
#endif
    } else {
#ifdef DUP_SUBSTRUCT_CACHE
      if (seed.getNumBonds() > 3) {
        DuplicateCache.add(seed.DupCacheKey,
                           false);  // opt. cache mismatched duplicates too
      }
#endif
    }
  }
  if (newSeed) {
    *newSeed = seed;  // non-blocking copy for MULTI_THREAD and best CPU
                      // utilization
  }

  return found;  // new matched seed has been actually added
}

bool MaximumCommonSubgraph::match(Seed& seed) {
  unsigned max_miss = Targets.size() - ThresholdCount;
  unsigned missing = 0;
  unsigned passed = 0;
  unsigned itarget = 0;

  for (std::vector<Target>::const_iterator tag = Targets.begin();
       tag != Targets.end(); tag++, itarget++) {
#ifdef VERBOSE_STATISTICS_ON
    { ++VerboseStatistics.MatchCall; }
#endif
    bool target_matched = false;
    if (!seed.MatchResult.empty() && !seed.MatchResult[itarget].empty()) {
      target_matched = matchIncrementalFast(seed, itarget);
    }
    if (!target_matched) {  // slow full match
      match_V_t match;      // THERE IS NO Bonds match INFO !!!!
      target_matched = SubstructMatchCustomTable(
          tag->Topology, *tag->Molecule, seed.Topology, *QueryMolecule,
          tag->AtomMatchTable, tag->BondMatchTable, &Parameters, &match);
      // save current match info
      if (target_matched) {
        if (seed.MatchResult.empty()) {
          seed.MatchResult.resize(Targets.size());
        }
        seed.MatchResult[itarget].init(seed, match, *QueryMolecule, *tag);
      } else if (!seed.MatchResult.empty()) {
        seed.MatchResult[itarget].clear();  //.Empty = true; // == fast clear();
      }
#ifdef VERBOSE_STATISTICS_ON
      if (target_matched) {
        ++VerboseStatistics.SlowMatchCallTrue;
      }
#endif
    }

    if (target_matched) {
      if (++passed >= ThresholdCount) {  // it's enough
        break;
      }
    } else {  // mismatched
      if (++missing > max_miss) {
        break;
      }
    }
  }
  if (missing <= max_miss) {
#ifdef VERBOSE_STATISTICS_ON
    ++VerboseStatistics.MatchCallTrue;
#endif
    return true;
  }
  return false;
}

// call it for each target, if failed perform full match check
bool MaximumCommonSubgraph::matchIncrementalFast(Seed& seed, unsigned itarget) {
// use and update results of previous match stored in the seed
#ifdef VERBOSE_STATISTICS_ON
  { ++VerboseStatistics.FastMatchCall; }
#endif
  const Target& target = Targets[itarget];
  TargetMatch& match = seed.MatchResult[itarget];
  if (match.empty()) {
    return false;
  }
  /*
  // CHIRALITY: FinalMatchCheck:
  if(Parameters.AtomCompareParameters.MatchChiralTag ||
  Parameters.FinalMatchChecker) {   // TEMP
          match.clear();
          return false;
  }
  */
  bool matched = false;
  for (unsigned newBondSeedIdx = match.MatchedBondSize;
       newBondSeedIdx < seed.getNumBonds(); newBondSeedIdx++) {
    matched = false;
    bool atomAdded = false;
    const Bond* newBond = seed.MoleculeFragment.Bonds[newBondSeedIdx];
    unsigned newBondQueryIdx = seed.MoleculeFragment.BondsIdx[newBondSeedIdx];

    unsigned newBondSourceAtomSeedIdx;   // seed's index of atom from which
                                         // new bond was added
    unsigned newBondAnotherAtomSeedIdx;  // seed's index of atom on
                                         // another end of the bond
    unsigned i =
        seed.MoleculeFragment.SeedAtomIdxMap[newBond->getBeginAtomIdx()];
    unsigned j = seed.MoleculeFragment.SeedAtomIdxMap[newBond->getEndAtomIdx()];
    if (i >= match.MatchedAtomSize) {  // this is new atom in the seed
      newBondSourceAtomSeedIdx = j;
      newBondAnotherAtomSeedIdx = i;
    } else {
      newBondSourceAtomSeedIdx = i;
      newBondAnotherAtomSeedIdx = j;
    }
    unsigned newBondAnotherAtomQueryIdx =
        seed.MoleculeFragment.AtomsIdx[newBondAnotherAtomSeedIdx];
    unsigned newBondSourceAtomQueryIdx =
        seed.MoleculeFragment.AtomsIdx[newBondSourceAtomSeedIdx];
    unsigned newBondSourceAtomTargetIdx =
        match.TargetAtomIdx
            [newBondSourceAtomQueryIdx];  // matched to
                                          // newBondSourceAtomSeedIdx

    const Bond* tb = nullptr;
    unsigned newBondAnotherAtomTargetIdx = NotSet;

    if (newBondAnotherAtomSeedIdx <
        match.MatchedAtomSize) {  // new bond between old atoms - both are
                                  // matched to exact atoms in the target
      newBondAnotherAtomTargetIdx =
          match.TargetAtomIdx[newBondAnotherAtomQueryIdx];
      tb = target.Molecule->getBondBetweenAtoms(
          newBondSourceAtomTargetIdx,
          newBondAnotherAtomTargetIdx);  // target bond between Source and
                                         // Another atoms
      if (tb) {  // bond exists, check match with query molecule
        unsigned tbi = tb->getIdx();
        unsigned qbi = seed.MoleculeFragment.BondsIdx[newBondSeedIdx];
        if (!match.VisitedTargetBonds[tbi]) {  // false if target bond is
                                               // already matched
          matched = target.BondMatchTable.at(qbi, tbi);
        }
      }
    } else {  // enumerate all bonds from source atom in the target
      const Atom* atom =
          target.Molecule->getAtomWithIdx(newBondSourceAtomTargetIdx);
      ROMol::OEDGE_ITER beg, end;
      for (boost::tie(beg, end) = target.Molecule->getAtomBonds(atom);
           beg != end; beg++) {
        tb = *beg;
        if (!match.VisitedTargetBonds[tb->getIdx()]) {
          newBondAnotherAtomTargetIdx = tb->getBeginAtomIdx();
          if (newBondSourceAtomTargetIdx == newBondAnotherAtomTargetIdx) {
            newBondAnotherAtomTargetIdx = tb->getEndAtomIdx();
          }

          if (newBondAnotherAtomSeedIdx <
                  seed.LastAddedAtomsBeginIdx  // RING: old atom, new atom
                                               // in matched substructure
                                               // must be new in seed
              || std::find(seed.MoleculeFragment.AtomsIdx.begin() +
                               seed.LastAddedAtomsBeginIdx,
                           seed.MoleculeFragment.AtomsIdx.begin() +
                               newBondAnotherAtomSeedIdx,
                           newBondAnotherAtomQueryIdx) ==
                     seed.MoleculeFragment.AtomsIdx.end()) {
            if (!match.VisitedTargetAtoms[newBondAnotherAtomTargetIdx]) {
              continue;
            }
          } else {
            if (match.VisitedTargetAtoms[newBondAnotherAtomTargetIdx]) {
              continue;
            }
          }

          // check AnotherAtom and bond
          matched =
              target.AtomMatchTable.at(newBondAnotherAtomQueryIdx,
                                       newBondAnotherAtomTargetIdx) &&
              target.BondMatchTable.at(
                  seed.MoleculeFragment.BondsIdx[newBondSeedIdx], tb->getIdx());

          if (matched) {
            atomAdded = true;
            break;
          }
        }
      }
    }

    if (matched) {      // update match history
      if (atomAdded) {  // new atom has been added
        match.MatchedAtomSize++;
        match.TargetAtomIdx[newBondAnotherAtomQueryIdx] =
            newBondAnotherAtomTargetIdx;
        match.VisitedTargetAtoms[newBondAnotherAtomTargetIdx] = true;
      }
      match.MatchedBondSize++;
      match.TargetBondIdx[newBondQueryIdx] = tb->getIdx();
      match.VisitedTargetBonds[tb->getIdx()] = true;
    } else {
      match.clear();
      return false;
    }
  }

  if (match.MatchedAtomSize != seed.getNumAtoms() ||
      match.MatchedBondSize !=
          seed.getNumBonds()) {  // number of unique items !!!
    match.clear();
    return false;
  }
  // CHIRALITY: FinalMatchCheck
  if (matched && Parameters.FinalMatchChecker) {
    std::vector<std::uint32_t> c1;
    c1.reserve(4096);
    std::vector<std::uint32_t> c2;
    c2.reserve(4096);
    for (unsigned si = 0; si < seed.getNumAtoms();
         si++) {  // index in the seed topology
      c1.push_back(si);
      c2.push_back(match.TargetAtomIdx[seed.Topology[si]]);
    }
    matched = Parameters.FinalMatchChecker(c1.data(), c2.data(), *QueryMolecule,
                                           seed.Topology, *target.Molecule,
                                           target.Topology,
                                           &Parameters);  // check CHIRALITY
    if (!matched) {
      match.clear();
    }
  }
#ifdef VERBOSE_STATISTICS_ON
  if (matched) {
#ifdef MULTI_THREAD
    Guard statlock(StatisticsMutex);
#endif
    ++VerboseStatistics.FastMatchCallTrue;
  }
#endif
  return matched;
}
}  // namespace FMCS
}  // namespace RDKit
