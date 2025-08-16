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
#include <GraphMol/FMCS/TwoMolMCSS.h>
#include "SubstructMatchCustom.h"
#include "MaximumCommonSubgraph.h"
#include <RDGeneral/BoostStartInclude.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {
namespace FMCS {

struct LabelDefinition {
  unsigned int ItemIndex;  // item with this label value
  unsigned int Value;
  LabelDefinition() : ItemIndex(NotSet), Value(NotSet) {}
  LabelDefinition(unsigned int i, unsigned int value)
      : ItemIndex(i), Value(value) {}
};

MaximumCommonSubgraph::MaximumCommonSubgraph(const MCSParameters *params) {
  Parameters = detail::MCSParametersInternal(params);
  if (!Parameters.ProgressCallback) {
    Parameters.ProgressCallback = MCSProgressCallbackTimeout;
    Parameters.ProgressCallbackUserData = &To;
  }
  if (Parameters.AtomCompareParameters.MatchChiralTag) {
    Parameters.BondCompareParameters.MatchStereo = true;
  }
  QueryMoleculeMatchedBonds = 0;
  QueryMoleculeMatchedAtoms = 0;
  QueryMoleculeSingleMatchedAtom = nullptr;
  To = nanoClock();
}

static bool molPtr_NumBondLess(
    const ROMol *l,
    const ROMol *r) {  // need for sorting the source molecules by size
  return l->getNumBonds() < r->getNumBonds();
}

void MaximumCommonSubgraph::init(size_t startIdx) {
  QueryMolecule = Molecules.at(startIdx);

  Targets.clear();
#ifdef FAST_SUBSTRUCT_CACHE
  QueryAtomLabels.clear();
  QueryBondLabels.clear();
  QueryAtomMatchTable.clear();
  QueryBondMatchTable.clear();
#endif
#ifdef DUP_SUBSTRUCT_CACHE
  DuplicateCache.clear();
#endif

  size_t nq = 0;
#ifdef FAST_SUBSTRUCT_CACHE
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
  unsigned int currentLabelValue = 1;
  std::vector<LabelDefinition> labels;
  nq = QueryMolecule->getNumAtoms();
  QueryAtomLabels.resize(nq, NotSet);
  for (size_t ai = 0; ai < nq; ++ai) {
    if (MCSAtomCompareAny ==
        Parameters.AtomTyper) {  // predefined functor without atom compare
                                 // parameters
      QueryAtomLabels[ai] = 1;
    } else {
      const auto atom = QueryMolecule->getAtomWithIdx(ai);
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
        for (auto &label : labels) {
          if (Parameters.AtomTyper(
                  Parameters.AtomCompareParameters, *QueryMolecule,
                  label.ItemIndex, *QueryMolecule, ai,
                  Parameters.CompareFunctionsUserData)) {  // equal items
            QueryAtomLabels[ai] = label.Value;
            break;
          }
        }
        if (NotSet ==
            QueryAtomLabels.at(ai)) {  // not found -> create new label
          QueryAtomLabels[ai] = ++currentLabelValue;
          labels.emplace_back(ai, currentLabelValue);
        }
      }
    }
  }
  labels.clear();
  currentLabelValue = 1;
  nq = QueryMolecule->getNumBonds();
  QueryBondLabels.resize(nq, NotSet);
  for (size_t aj = 0; aj < nq; ++aj) {
    const Bond *bond = QueryMolecule->getBondWithIdx(aj);
    unsigned ring = 0;
    if (!Parameters.CompareFunctionsUserData &&
        (Parameters.BondCompareParameters.CompleteRingsOnly ||
         Parameters.BondCompareParameters.RingMatchesRingOnly)) {
      // is bond in ring
      ring = QueryMolecule->getRingInfo()->numBondRings(aj) ? 0 : 1;
    }
    if (MCSBondCompareAny == Parameters.BondTyper) {
      QueryBondLabels[aj] = 1 | (ring >> 8);
    } else if (MCSBondCompareOrderExact == Parameters.BondTyper) {
      QueryBondLabels[aj] = (bond->getBondType() + 1) | (ring >> 8);
    } else if (MCSBondCompareOrder == Parameters.BondTyper) {
      auto order = bond->getBondType();
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
      for (const auto &label : labels) {
        if (Parameters.BondTyper(
                Parameters.BondCompareParameters, *QueryMolecule,
                label.ItemIndex, *QueryMolecule, aj,
                Parameters
                    .CompareFunctionsUserData)) {  // equal bonds + ring ...
          QueryBondLabels[aj] = label.Value;
          break;
        }
      }
      if (NotSet == QueryBondLabels.at(aj)) {  // not found -> create new label
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
    for (const auto &a : Targets.at(i).Molecule->atoms()) {
      Targets[i].Topology.addAtom(a->getIdx());
    }
    // build Target Topology ADD BONDs
    for (const auto &b : Targets.at(i).Molecule->bonds()) {
      auto ii = b->getBeginAtomIdx();
      auto jj = b->getEndAtomIdx();
      Targets[i].Topology.addBond(b->getIdx(), ii, jj);
    }

    // fill out match tables
    size_t nq = QueryMolecule->getNumAtoms();
    size_t nt = Targets.at(i).Molecule->getNumAtoms();
    Targets[i].AtomMatchTable.resize(nq, nt);

    for (size_t aj = 0; aj < nt; aj++) {
      for (size_t ai = 0; ai < nq; ai++) {
        Targets[i].AtomMatchTable.set(
            ai, aj,
            Parameters.AtomTyper(Parameters.AtomCompareParameters,
                                 *QueryMolecule, ai, *Targets.at(i).Molecule,
                                 aj, Parameters.CompareFunctionsUserData));
      }
    }
    nq = QueryMolecule->getNumBonds();
    nt = Targets.at(i).Molecule->getNumBonds();
    Targets[i].BondMatchTable.resize(nq, nt);
    for (size_t aj = 0; aj < nt; aj++) {
      for (size_t ai = 0; ai < nq; ai++) {
        Targets[i].BondMatchTable.set(
            ai, aj,
            Parameters.BondTyper(Parameters.BondCompareParameters,
                                 *QueryMolecule, ai, *Targets.at(i).Molecule,
                                 aj, Parameters.CompareFunctionsUserData));
      }
    }
  }
}

struct WeightedBond {
  const Bond *BondPtr{nullptr};
  unsigned int Weight{0};
  WeightedBond() {}
  WeightedBond(const Bond *bond) : BondPtr(bond), Weight(0) {
    const auto ringInfo = bond->getOwningMol().getRingInfo();
    // score ((bond.is_in_ring + atom1.is_in_ring + atom2.is_in_ring)
    if (ringInfo->numBondRings(bond->getIdx())) {
      ++Weight;
    }
    if (ringInfo->numAtomRings(bond->getBeginAtomIdx())) {
      ++Weight;
    }
    if (ringInfo->numAtomRings(bond->getEndAtomIdx())) {
      ++Weight;
    }
  }
  bool operator<(const WeightedBond &r) {
    return Weight >= r.Weight;  // sort in Z-A order (Rings first)
  }
};

void MaximumCommonSubgraph::makeInitialSeeds() {
  // build a set of initial seeds as "all" single bonds from query
  // molecule
  boost::dynamic_bitset<> excludedBonds(QueryMolecule->getNumBonds());

  Seeds.clear();
  QueryMoleculeMatchedBonds = 0;
  QueryMoleculeMatchedAtoms = 0;
  QueryMoleculeSingleMatchedAtom = nullptr;
  std::vector<MatchVectType> matching_substructs;
  std::unique_ptr<const ROMol> initialSeedMolecule;

  if (Parameters.InitialSeed.empty() && Parameters.FastInitialSeed &&
      Targets.size() == 1) {
    // We can do the clique detection MCS to create seeds which should
    // be the final answer.  This just plumbs it in as a check and so that
    // all the downstream results preparation is the same.
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> maxCliques;
    TwoMolMCSS(*QueryMolecule, *Targets[0].Molecule, Parameters.MinMCSSSize,
               Targets[0].AtomMatchTable, Targets[0].BondMatchTable,
               maxCliques);
    if (maxCliques.empty()) {
      return;
    }
    // Use a trimmed down copy of QueryMolecule as the initialSeedMolecule, and
    // the clique as the match of it onto the QueryMolecule.  Only use the first
    // clique, as they should all be the same barring symmetry.
    std::unique_ptr<RWMol> qmol(new RWMol(*QueryMolecule));

    const auto &clique = maxCliques.front();
    std::unordered_set<unsigned int> atomsToKeep;
    std::transform(clique.begin(), clique.end(),
                   std::inserter(atomsToKeep, atomsToKeep.begin()),
                   [](const auto &a) -> unsigned int { return a.first; });
    for (auto atom : qmol->atoms()) {
      atom->setProp<int>("OrigIdx", atom->getIdx());
    }
    qmol->beginBatchEdit();
    for (auto atom : qmol->atoms()) {
      if (atomsToKeep.find(atom->getIdx()) == atomsToKeep.end()) {
        qmol->removeAtom(atom);
      }
    }
    qmol->commitBatchEdit();
    MatchVectType nextMatch;
    nextMatch.reserve(qmol->getNumAtoms());
    for (auto atom : qmol->atoms()) {
      nextMatch.emplace_back(
          std::pair<int, int>(atom->getIdx(), atom->getProp<int>("OrigIdx")));
    }
    matching_substructs.emplace_back(std::move(nextMatch));
    initialSeedMolecule.reset(static_cast<const ROMol *>(qmol.release()));
  } else if (!Parameters.InitialSeed.empty()) {
    // make user defined seed
    initialSeedMolecule.reset(
        static_cast<const ROMol *>(SmartsToMol(Parameters.InitialSeed)));
    // make a set of seeds as indices and pointers to current query
    // molecule items based on matching results
    std::vector<MatchVectType> matching_substructs;
    SubstructMatch(*QueryMolecule, *initialSeedMolecule, matching_substructs);
  }

  if (initialSeedMolecule) {
    std::cout << "Initial seed : " << MolToSmarts(*initialSeedMolecule) << " : "
              << initialSeedMolecule->getNumAtoms() << " and "
              << initialSeedMolecule->getNumBonds() << std::endl;
    // loop throw all fragments of Query matched to initial seed
    for (const auto &ms : matching_substructs) {
      Seed seed;
      seed.setStoreAllDegenerateMCS(Parameters.StoreAll);
      seed.ExcludedBonds = excludedBonds;
      seed.MatchResult.resize(Targets.size());
#ifdef VERBOSE_STATISTICS_ON
      {
        ++VerboseStatistics.Seed;
        ++VerboseStatistics.InitialSeed;
      }
#endif
      // add all matched atoms of the matched query fragment
      std::map<unsigned int, unsigned int> initialSeedToQueryAtom;
      for (const auto &msb : ms) {
        unsigned int qai = msb.second;
        unsigned int sai = msb.first;
        seed.addAtom(QueryMolecule->getAtomWithIdx(qai));
        initialSeedToQueryAtom[sai] = qai;
      }
      // add all bonds (existed in initial seed !!!) between all matched
      // atoms in query
      for (const auto &msb : ms) {
        const auto atom = initialSeedMolecule->getAtomWithIdx(msb.first);
        for (const auto &nbri : boost::make_iterator_range(
                 initialSeedMolecule->getAtomBonds(atom))) {
          const auto initialBond = (*initialSeedMolecule)[nbri];
          unsigned int qai1 =
              initialSeedToQueryAtom.at(initialBond->getBeginAtomIdx());
          unsigned int qai2 =
              initialSeedToQueryAtom.at(initialBond->getEndAtomIdx());

          const auto b = QueryMolecule->getBondBetweenAtoms(qai1, qai2);
          CHECK_INVARIANT(b, "bond must not be NULL");
          if (!seed.ExcludedBonds.test(b->getIdx())) {
            seed.addBond(b);
            seed.ExcludedBonds.set(b->getIdx());
          }
        }
      }
      seed.computeRemainingSize(*QueryMolecule);

      if (checkIfMatchAndAppend(seed)) {
        QueryMoleculeMatchedBonds = seed.getNumBonds();
      }
    }
    if (Seeds.empty()) {
      BOOST_LOG(rdWarningLog)
          << "The provided InitialSeed is not an MCS and will be ignored : "
          << Parameters.InitialSeed << " : " << MolToSmiles(*QueryMolecule)
          << " vs " << MolToSmiles(*Targets[0].Molecule) << " with "
          << MolToSmarts(*initialSeedMolecule) << std::endl;
    }
  }
  if (Seeds.empty()) {  // create a set of seeds from each query bond
    // R1 additional performance OPTIMISATION
    // if(Parameters.BondCompareParameters.CompleteRingsOnly)
    // disable all mismatched rings, and do not generate initial seeds
    // from such disabled bonds
    //  for(  rings .....) for(i......)
    //   if(mismatched) excludedBonds[i.......] = true;
    std::vector<WeightedBond> wbVec;
    wbVec.reserve(QueryMolecule->getNumBonds());
    for (const auto &bond : QueryMolecule->bonds()) {
      wbVec.emplace_back(bond);
    }

    for (const auto &wb : wbVec) {
      // R1 additional performance OPTIMISATION
      // if(excludedBonds[(*bi)->getIdx()])
      //    continue;
      Seed seed;
      seed.setStoreAllDegenerateMCS(Parameters.StoreAll);
      seed.MatchResult.resize(Targets.size());

#ifdef VERBOSE_STATISTICS_ON
      {
        ++VerboseStatistics.Seed;
        ++VerboseStatistics.InitialSeed;
      }
#endif
      seed.addAtom(wb.BondPtr->getBeginAtom());
      seed.addAtom(wb.BondPtr->getEndAtom());
      seed.ExcludedBonds = excludedBonds;  // all bonds from first to current
      seed.addBond(wb.BondPtr);
      excludedBonds.set(wb.BondPtr->getIdx());

      seed.computeRemainingSize(*QueryMolecule);

      if (checkIfMatchAndAppend(seed)) {
        ++QueryMoleculeMatchedBonds;
      } else {
        // optionally remove all such bonds from all targets TOPOLOGY
        // where it exists.
        //..........

        // disable (mark as already processed) mismatched bond in all
        // seeds
        for (auto &Seed : Seeds) {
          Seed.ExcludedBonds.set(wb.BondPtr->getIdx());
        }

#ifdef VERBOSE_STATISTICS_ON
        ++VerboseStatistics.MismatchedInitialSeed;
#endif
      }
    }
  }
  auto nq = QueryMolecule->getNumAtoms();
  MatchVectType singleAtomPairMatch(1);
  MatchVectType emptyBondPairMatch;
  for (size_t i = 0; i < nq; i++) {  // all query's atoms
    const auto queryMolAtom = QueryMolecule->getAtomWithIdx(i);
    bool isQueryMolAtomInRing = queryIsAtomInRing(queryMolAtom);
    unsigned int matched = 0;
    const Atom *candQueryMoleculeSingleMatchedAtom = nullptr;
    for (const auto &tag : Targets) {
      auto nt = tag.Molecule->getNumAtoms();
      for (size_t aj = 0; aj < nt; aj++) {
        if (tag.AtomMatchTable.at(i, aj)) {
          const auto targetMolAtom = tag.Molecule->getAtomWithIdx(aj);
          bool isTargetMolAtomInRing = queryIsAtomInRing(targetMolAtom);
          ++matched;
          if (!(Parameters.BondCompareParameters.CompleteRingsOnly &&
                (isQueryMolAtomInRing || isTargetMolAtomInRing))) {
            bool shouldAccept = !Parameters.ShouldAcceptMCS;
            if (!shouldAccept) {
              singleAtomPairMatch[0] = std::make_pair(i, aj);
              shouldAccept = Parameters.ShouldAcceptMCS(
                  *QueryMolecule, *tag.Molecule, singleAtomPairMatch,
                  emptyBondPairMatch, &Parameters);
            }
            if (shouldAccept) {
              candQueryMoleculeSingleMatchedAtom = queryMolAtom;
            }
          }
          break;
        }
      }
    }
    if (matched && matched >= ThresholdCount) {
      ++QueryMoleculeMatchedAtoms;
      if (candQueryMoleculeSingleMatchedAtom) {
        if (!QueryMoleculeSingleMatchedAtom) {
          QueryMoleculeSingleMatchedAtom = candQueryMoleculeSingleMatchedAtom;
        } else {
          QueryMoleculeSingleMatchedAtom =
              (std::max)(candQueryMoleculeSingleMatchedAtom,
                         QueryMoleculeSingleMatchedAtom,
                         [](const Atom *a, const Atom *b) {
                           if (a->getDegree() != b->getDegree()) {
                             return (a->getDegree() < b->getDegree());
                           } else if (a->getFormalCharge() !=
                                      b->getFormalCharge()) {
                             return (a->getFormalCharge() <
                                     b->getFormalCharge());
                           } else if (a->getAtomicNum() != b->getAtomicNum()) {
                             return (a->getAtomicNum() < b->getAtomicNum());
                           }
                           return (a->getIdx() < b->getIdx());
                         });
        }
      }
    }
  }
}

namespace {

bool checkIfRingsAreClosed(const Seed &fs, bool noLoneRingAtoms) {
  if (fs.MoleculeFragment.Bonds.empty() && fs.MoleculeFragment.Atoms.empty()) {
    return true;
  }
  const auto &om = fs.MoleculeFragment.Atoms.front()->getOwningMol();
  const auto ri = om.getRingInfo();
  if (!ri->numRings()) {
    return true;
  }
  boost::dynamic_bitset<> mcsBonds(om.getNumBonds());
  boost::dynamic_bitset<> mcsNonFusedRings(ri->numRings());
  boost::dynamic_bitset<> mcsFusedRings(ri->numRings());
  for (const auto &bond : fs.MoleculeFragment.Bonds) {
    auto bi = bond->getIdx();
    mcsBonds.set(bi);
    if (ri->numBondRings(bi) == 1) {
      mcsNonFusedRings.set(ri->bondMembers(bi).front());
    }
  }
  for (unsigned int ringIdx = 0; ringIdx < mcsNonFusedRings.size(); ++ringIdx) {
    if (!mcsNonFusedRings.test(ringIdx)) {
      continue;
    }
    for (const auto &bi : ri->bondRings().at(ringIdx)) {
      bool keepBond = false;
      for (unsigned int memberOf : ri->bondMembers(bi)) {
        if (memberOf == ringIdx) {
          keepBond = true;
        } else if (mcsNonFusedRings.test(memberOf)) {
          keepBond = false;
          break;
        }
      }
      if (keepBond && !mcsBonds.test(bi)) {
        return false;
      }
    }
  }
  if (noLoneRingAtoms) {
    for (const auto &atom : fs.MoleculeFragment.Atoms) {
      auto ai = atom->getIdx();
      const auto &ringIndices = ri->atomMembers(ai);
      if (!ringIndices.empty() &&
          !std::any_of(ringIndices.begin(), ringIndices.end(),
                       [&mcsNonFusedRings](const auto &ringIdx) {
                         return mcsNonFusedRings.test(ringIdx);
                       })) {
        return false;
      }
    }
  }
  if (mcsNonFusedRings.none()) {
    for (const auto &bond : fs.MoleculeFragment.Bonds) {
      auto bi = bond->getIdx();
      if (ri->numBondRings(bi) > 1) {
        for (auto ringIdx : ri->bondMembers(bi)) {
          mcsFusedRings.set(ringIdx);
        }
      }
    }
  }
  if (mcsFusedRings.any()) {
    for (unsigned int ringIdx = 0; ringIdx < mcsFusedRings.size(); ++ringIdx) {
      if (!mcsFusedRings.test(ringIdx)) {
        continue;
      }
      const auto &ringBondIndices = ri->bondRings().at(ringIdx);
      if (std::all_of(
              ringBondIndices.begin(), ringBondIndices.end(),
              [&mcsBonds](const auto &bi) { return mcsBonds.test(bi); })) {
        return true;
      }
    }
    return false;
  }
  return true;
}

bool checkIfShouldAcceptMCS(const FMCS::MolFragment &f, const ROMol &query,
                            const std::vector<Target> &targets,
                            const MCSParameters &p) {
  if (!p.ShouldAcceptMCS || f.Bonds.empty()) {
    return true;
  }
  Seed seed;  // result MCS
  seed.ExcludedBonds.resize(query.getNumBonds(), false);

  for (const auto &atom : f.Atoms) {
    seed.addAtom(atom);
  }
  for (const auto &bond : f.Bonds) {
    seed.addBond(bond);
  }
  for (const auto &target : targets) {
    match_V_t match;
    bool targetMatched = SubstructMatchCustomTable(
        target.Topology, *target.Molecule, seed.Topology, query,
        target.AtomMatchTable, target.BondMatchTable, &p, &match);
    // it is OK to skip non-matches here as threshold could be < 1.0;
    // we only want to triage matches
    if (!targetMatched) {
      continue;
    }
    MatchVectType atomPairMatch;
    atomPairMatch.reserve(match.size());
    std::vector<unsigned int> queryToTargetAtomIdx(query.getNumAtoms(), NotSet);
    for (const auto &m : match) {
      auto queryAtomIdx = seed.Topology[m.first];
      auto targetAtomIdx = target.Topology[m.second];
      atomPairMatch.emplace_back(queryAtomIdx, targetAtomIdx);
      queryToTargetAtomIdx[queryAtomIdx] = targetAtomIdx;
    }
    MatchVectType bondPairMatch;
    auto numMcsBonds = boost::num_edges(seed.Topology);
    bondPairMatch.reserve(numMcsBonds);
    for (const auto &it :
         boost::make_iterator_range(boost::edges(seed.Topology))) {
      int queryBeginAtomIdx = seed.Topology[boost::source(it, seed.Topology)];
      int queryEndAtomIdx = seed.Topology[boost::target(it, seed.Topology)];
      const auto queryBond =
          query.getBondBetweenAtoms(queryBeginAtomIdx, queryEndAtomIdx);
      CHECK_INVARIANT(queryBond, "");
      const auto targetBeginAtomIdx =
          queryToTargetAtomIdx.at(queryBeginAtomIdx);
      CHECK_INVARIANT(targetBeginAtomIdx != NotSet, "");
      const auto targetEndAtomIdx = queryToTargetAtomIdx.at(queryEndAtomIdx);
      CHECK_INVARIANT(targetEndAtomIdx != NotSet, "");
      const auto targetBond = target.Molecule->getBondBetweenAtoms(
          targetBeginAtomIdx, targetEndAtomIdx);
      CHECK_INVARIANT(targetBond, "");
      bondPairMatch.emplace_back(queryBond->getIdx(), targetBond->getIdx());
    }
    if (!p.ShouldAcceptMCS(query, *target.Molecule, atomPairMatch,
                           bondPairMatch, &p)) {
      return false;
    }
  }
  return true;
}

}  // namespace
bool MaximumCommonSubgraph::growSeeds() {
  bool mcsFound = false;
  bool canceled = false;
  // Find MCS -- SDF Seed growing OPTIMISATION (it works in 3 times
  // faster)
  while (!Seeds.empty()) {
    if (getMaxNumberBonds() == QueryMoleculeMatchedBonds) {  // MCS == Query
      break;
    }
#ifdef VERBOSE_STATISTICS_ON
    VerboseStatistics.TotalSteps++;
#endif
    auto si = Seeds.begin();

    si->grow(*this);
    {
      const Seed &fs = Seeds.front();
      // bigger substructure found
      if (fs.CopyComplete) {
        bool possibleMCS = false;
        if (!Parameters.MaximizeBonds) {
          possibleMCS = (fs.getNumAtoms() > getMaxNumberAtoms() ||
                         (fs.getNumAtoms() == getMaxNumberAtoms() &&
                          fs.getNumBonds() > getMaxNumberBonds()));
        } else {
          possibleMCS = (fs.getNumBonds() > getMaxNumberBonds() ||
                         (fs.getNumBonds() == getMaxNumberBonds() &&
                          fs.getNumAtoms() > getMaxNumberAtoms()));
        }
        bool isDegenerateMCS = (fs.getNumBonds() == getMaxNumberBonds() &&
                                fs.getNumAtoms() == getMaxNumberAtoms());
        if (!possibleMCS && Parameters.StoreAll) {
          possibleMCS = isDegenerateMCS;
        }
        // #945: test here to see if the MCS actually has all rings closed
        if (possibleMCS && Parameters.BondCompareParameters.CompleteRingsOnly) {
          possibleMCS = checkIfRingsAreClosed(
              fs, Parameters.AtomCompareParameters.CompleteRingsOnly);
        }
        if (possibleMCS) {
          possibleMCS = checkIfShouldAcceptMCS(
              fs.MoleculeFragment, *QueryMolecule, Targets, Parameters);
        }
        if (possibleMCS) {
          mcsFound = true;
#ifdef VERBOSE_STATISTICS_ON
          VerboseStatistics.MCSFoundStep = VerboseStatistics.TotalSteps;
          VerboseStatistics.MCSFoundTime = nanoClock();
#endif
          McsIdx.Atoms = fs.MoleculeFragment.Atoms;
          McsIdx.Bonds = fs.MoleculeFragment.Bonds;
          if (Parameters.Verbose) {
            std::cout << VerboseStatistics.TotalSteps
                      << " Seeds:" << Seeds.size() << " MCS "
                      << McsIdx.Atoms.size() << " atoms, "
                      << McsIdx.Bonds.size() << " bonds";
            printf(" for %.4lf seconds. bond[0]=%u\n",
                   double(VerboseStatistics.MCSFoundTime - To) / 1000000.,
                   McsIdx.Bonds.front()->getIdx());
          }
          if (Parameters.StoreAll) {
            if (!isDegenerateMCS) {
              DegenerateMcsMap.clear();
            }
            std::vector<unsigned int> key(McsIdx.Bonds.size());
            std::transform(McsIdx.Bonds.begin(), McsIdx.Bonds.end(),
                           key.begin(),
                           [](const auto bond) { return bond->getIdx(); });
            std::sort(key.begin(), key.end());
            MCS value(McsIdx);
            value.QueryMolecule = QueryMolecule;
            value.Targets = Targets;
            DegenerateMcsMap[key] = value;
          }
        }
      }
    }
    if (NotSet == si->GrowingStage) {  // finished
      Seeds.erase(si);
    }
    if (Parameters.ProgressCallback) {
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
  unsigned int QueryAtomIdx;
  unsigned int TargetAtomIdx;
  AtomMatch() : QueryAtomIdx(NotSet), TargetAtomIdx(NotSet) {}
};
typedef std::vector<AtomMatch> AtomMatchSet;

std::pair<std::string, ROMOL_SPTR>
MaximumCommonSubgraph::generateResultSMARTSAndQueryMol(
    const MCS &mcsIdx) const {
  // match the result MCS with all targets to check if it is exact match
  // or template
  Seed seed;  // result MCS
  seed.setStoreAllDegenerateMCS(Parameters.StoreAll);
  seed.ExcludedBonds.resize(mcsIdx.QueryMolecule->getNumBonds(), false);
  std::vector<AtomMatchSet> atomMatchResult(mcsIdx.Targets.size());
  std::vector<unsigned int> atomIdxMap(mcsIdx.QueryMolecule->getNumAtoms());
  std::vector<std::map<unsigned int, const Bond *>> bondMatchSet(
      mcsIdx.Bonds.size());  // key is unique BondType
  std::vector<std::map<unsigned int, const Atom *>> atomMatchSet(
      mcsIdx.Atoms.size());  // key is unique atomic number

  for (const auto &atom : mcsIdx.Atoms) {
    atomIdxMap[atom->getIdx()] = seed.getNumAtoms();
    seed.addAtom(atom);
  }
  for (const auto &bond : mcsIdx.Bonds) {
    seed.addBond(bond);
  }

  std::vector<unsigned int> matchedTargetIndices;
  if (!mcsIdx.Bonds.empty()) {
    for (const auto &tag : mcsIdx.Targets) {
      match_V_t match;  // THERE IS NO Bonds match INFO !!!!
      bool target_matched = SubstructMatchCustomTable(
          tag.Topology, *tag.Molecule, seed.Topology, *QueryMolecule,
          tag.AtomMatchTable, tag.BondMatchTable, &Parameters, &match);
      if (!target_matched) {
        continue;
      }
      unsigned int itarget = &tag - &mcsIdx.Targets.front();
      matchedTargetIndices.push_back(itarget);
      atomMatchResult.at(itarget).resize(seed.getNumAtoms());
      for (const auto &m : match) {
        const auto ai = m.first;  // SeedAtomIdx
        atomMatchResult.at(itarget).at(ai).QueryAtomIdx =
            seed.Topology[m.first];
        atomMatchResult.at(itarget).at(ai).TargetAtomIdx =
            tag.Topology[m.second];
        const auto ta = tag.Molecule->getAtomWithIdx(tag.Topology[m.second]);
        if (ta->getAtomicNum() !=
            seed.MoleculeFragment.Atoms.at(ai)->getAtomicNum()) {
          atomMatchSet[ai][ta->getAtomicNum()] = ta;  // add
        }
      }
      // AND BUILD BOND MATCH INFO
      for (const auto &bond : mcsIdx.Bonds) {
        const auto bi = &bond - &mcsIdx.Bonds.front();
        const auto i = atomIdxMap.at(bond->getBeginAtomIdx());
        const auto j = atomIdxMap.at(bond->getEndAtomIdx());
        const auto ti = atomMatchResult.at(itarget).at(i).TargetAtomIdx;
        const auto tj = atomMatchResult.at(itarget).at(j).TargetAtomIdx;
        const auto tb = tag.Molecule->getBondBetweenAtoms(ti, tj);
        if (tb && bond->getBondType() != tb->getBondType()) {
          bondMatchSet[bi][tb->getBondType()] = tb;  // add
        }
      }
    }
  }

  // Generate result's SMARTS

  // create molecule from MCS for MolToSmarts()
  auto mol = new RWMol();
  ROMOL_SPTR molSptr(mol);
  const auto ri = mcsIdx.QueryMolecule->getRingInfo();
  boost::dynamic_bitset<> mcsRingIsComplete;
  bool needAtomRingQueries =
      (Parameters.AtomCompareParameters.RingMatchesRingOnly ||
       Parameters.BondCompareParameters.MatchFusedRingsStrict);
  if (needAtomRingQueries) {
    mcsRingIsComplete.resize(ri->numRings(), true);
    boost::dynamic_bitset<> queryBondInMcs(mcsIdx.QueryMolecule->getNumBonds());
    for (const auto &bond : mcsIdx.Bonds) {
      queryBondInMcs.set(bond->getIdx());
    }
    const auto &bondRings = ri->bondRings();
    for (const auto &bondRing : bondRings) {
      auto ringIdx = &bondRing - &bondRings.front();
      for (const auto &bondIdx : bondRing) {
        if (!queryBondInMcs.test(bondIdx)) {
          mcsRingIsComplete.reset(ringIdx);
          break;
        }
      }
    }
  }
  for (const auto &atom : mcsIdx.Atoms) {
    auto queryAtomIdx = atom->getIdx();
    auto numAtomRings = ri->numAtomRings(queryAtomIdx);
    QueryAtom a;
    const auto ai = &atom - &mcsIdx.Atoms.front();
    if (Parameters.AtomTyper == MCSAtomCompareIsotopes ||
        Parameters.AtomCompareParameters
            .MatchIsotope) {  // do '[0*]-[0*]-[13*]' for CC[13NH2]
      a.setQuery(makeAtomIsotopeQuery(static_cast<int>(atom->getIsotope())));
    } else {
      // generate [#6] instead of C or c !
      a.setQuery(makeAtomNumQuery(atom->getAtomicNum()));
      // for all atomMatchSet[ai] items add atom query to template like
      // [#6,#17,#9, ... ]
      for (const auto &am : atomMatchSet.at(ai)) {
        a.expandQuery(makeAtomNumQuery(am.second->getAtomicNum()),
                      Queries::COMPOSITE_OR);
        if (Parameters.AtomCompareParameters.MatchChiralTag &&
            (am.second->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
             am.second->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW)) {
          a.setChiralTag(am.second->getChiralTag());
        }
      }
    }
    if (needAtomRingQueries) {
      const auto &ringIndicesAtomIsMemberOf = ri->atomMembers(queryAtomIdx);
      auto numCompleteRings = std::count_if(
          ringIndicesAtomIsMemberOf.begin(), ringIndicesAtomIsMemberOf.end(),
          [&mcsRingIsComplete](const auto &ringIdx) {
            return mcsRingIsComplete.test(ringIdx);
          });
      if (Parameters.AtomCompareParameters.RingMatchesRingOnly &&
          !numCompleteRings) {
        auto q = makeAtomInRingQuery();
        q->setNegation(!numAtomRings);
        a.expandQuery(q, Queries::COMPOSITE_AND, true);
      } else if (Parameters.BondCompareParameters.MatchFusedRingsStrict &&
                 numAtomRings == 1 && numCompleteRings == 1) {
        auto ringSize =
            ri->atomRings().at(ringIndicesAtomIsMemberOf.front()).size();
        auto q = new ATOM_OR_QUERY;
        q->setDescription("AtomOr");
        q->addChild(QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(
            makeAtomMinRingSizeQuery(ringSize)));
        auto q2 = makeAtomInNRingsQuery(1);
        q2->setNegation(true);
        q->addChild(QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(q2));
        a.expandQuery(q, Queries::COMPOSITE_AND, true);
      }
    }
    mol->addAtom(&a, true, false);
  }
  for (const auto &bond : mcsIdx.Bonds) {
    QueryBond b;
    const auto bi = &bond - &mcsIdx.Bonds.front();
    const auto beginAtomIdx = atomIdxMap.at(bond->getBeginAtomIdx());
    const auto endAtomIdx = atomIdxMap.at(bond->getEndAtomIdx());
    b.setBeginAtomIdx(beginAtomIdx);
    b.setEndAtomIdx(endAtomIdx);
    b.setQuery(makeBondOrderEqualsQuery(bond->getBondType()));
    // add OR template if need
    for (const auto &bm : bondMatchSet.at(bi)) {
      b.expandQuery(makeBondOrderEqualsQuery(bm.second->getBondType()),
                    Queries::COMPOSITE_OR);
      if (Parameters.BondCompareParameters.MatchStereo &&
          bm.second->getStereo() > Bond::STEREOANY) {
        b.setStereo(bm.second->getStereo());
      }
    }
    if (Parameters.BondCompareParameters.RingMatchesRingOnly ||
        Parameters.BondCompareParameters.MatchFusedRingsStrict) {
      const auto numBondRings = ri->numBondRings(bond->getIdx());
      auto q = makeBondIsInRingQuery();
      q->setNegation(!numBondRings);
      b.expandQuery(q, Queries::COMPOSITE_AND, true);
    }
    mol->addBond(&b, false);
  }
  return std::make_pair(MolToSmarts(*mol, true), molSptr);
}

bool MaximumCommonSubgraph::createSeedFromMCS(size_t newQueryTarget,
                                              Seed &newSeed) {
  Seed mcs;
  mcs.setStoreAllDegenerateMCS(Parameters.StoreAll);
  mcs.ExcludedBonds.resize(McsIdx.QueryMolecule->getNumBonds(), false);
  std::vector<unsigned int> mcsAtomIdxMap(McsIdx.QueryMolecule->getNumAtoms());

  for (const auto &atom : McsIdx.Atoms) {
    mcsAtomIdxMap[atom->getIdx()] = mcs.addAtom(atom);
  }
  for (const auto &bond : McsIdx.Bonds) {
    mcs.addBond(bond);
  }

  const Target &newQuery = McsIdx.Targets.at(newQueryTarget);

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

  for (const auto &m : match) {
    unsigned int ai = m.first;  // SeedAtomIdx in mcs seed
    atomMatchResult[ai].QueryAtomIdx = mcs.Topology[m.first];
    atomMatchResult[ai].TargetAtomIdx = newQuery.Topology[m.second];
    const auto ta =
        newQuery.Molecule->getAtomWithIdx(newQuery.Topology[m.second]);
    newSeed.addAtom(ta);
  }

  for (const auto &bond : McsIdx.Bonds) {
    unsigned int i = mcsAtomIdxMap.at(bond->getBeginAtomIdx());
    unsigned int j = mcsAtomIdxMap.at(bond->getEndAtomIdx());
    unsigned int ti = atomMatchResult.at(i).TargetAtomIdx;
    unsigned int tj = atomMatchResult.at(j).TargetAtomIdx;
    const auto tb = newQuery.Molecule->getBondBetweenAtoms(ti, tj);
    CHECK_INVARIANT(tb, "tb most not be NULL");
    newSeed.addBond(tb);
  }
  newSeed.computeRemainingSize(*newQuery.Molecule);
  return true;
}

MCSResult MaximumCommonSubgraph::find(const std::vector<ROMOL_SPTR> &src_mols) {
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

  // minimal required number of matched targets:
  // at least one target, max all targets
  ThresholdCount = static_cast<unsigned int>(std::min(
      static_cast<int>(src_mols.size()) - 1,
      std::max(1, static_cast<int>(ceil(static_cast<double>(src_mols.size()) *
                                        Parameters.Threshold)) -
                      1)));

  // AtomCompareParameters.CompleteRingsOnly implies
  // BondCompareParameters.CompleteRingsOnly
  if (Parameters.AtomCompareParameters.CompleteRingsOnly) {
    Parameters.BondCompareParameters.CompleteRingsOnly = true;
  }

  // Selecting CompleteRingsOnly option also enables
  // --ring-matches-ring-only. ring--ring and chain bonds only match chain
  // bonds.
  if (Parameters.BondCompareParameters.CompleteRingsOnly) {
    Parameters.BondCompareParameters.RingMatchesRingOnly = true;
  }
  if (Parameters.AtomCompareParameters.CompleteRingsOnly) {
    Parameters.AtomCompareParameters.RingMatchesRingOnly = true;
  }

  unsigned int i = 0;
  boost::dynamic_bitset<> faked_ring_info(src_mols.size());
  for (const auto &src_mol : src_mols) {
    Molecules.push_back(src_mol.get());
    if (!Molecules.back()->getRingInfo()->isInitialized()) {
      Molecules.back()->getRingInfo()->initialize();  // but do not fill out !!!
      faked_ring_info.set(i);
    }
    ++i;
  }
  // sort source set of molecules by their 'size' and assume the smallest
  // molecule as a query
  std::stable_sort(Molecules.begin(), Molecules.end(), molPtr_NumBondLess);
  size_t startIdx = 0;
  size_t endIdx = Molecules.size() - ThresholdCount;
  while (startIdx < endIdx && !Molecules.at(startIdx)->getNumAtoms()) {
    ++startIdx;
  }
  bool areSeedsEmpty = false;
  for (size_t i = startIdx; i < endIdx && !areSeedsEmpty && !res.Canceled;
       ++i) {
    init(startIdx);
    if (Targets.empty()) {
      break;
    }
    MCSFinalMatchCheckFunction tff = Parameters.FinalMatchChecker;
    // skip final match check for initial seed to allow future growing
    Parameters.FinalMatchChecker = nullptr;
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
    res.Canceled = !(areSeedsEmpty || growSeeds());
    // verify what MCS is equal to one of initial seed for chirality match
    if (getMaxNumberBonds() == 0) {
      McsIdx = MCS();      // clear
      makeInitialSeeds();  // check all possible initial seeds
      if (!areSeedsEmpty) {
        const Seed &fs = Seeds.front();
        if ((1 == getMaxNumberBonds() ||
             !(Parameters.BondCompareParameters.CompleteRingsOnly &&
               fs.MoleculeFragment.Bonds.size() == 1 &&
               queryIsBondInRing(fs.MoleculeFragment.Bonds.front()))) &&
            checkIfShouldAcceptMCS(fs.MoleculeFragment, *QueryMolecule, Targets,
                                   Parameters)) {
          McsIdx.QueryMolecule = QueryMolecule;
          McsIdx.Targets = Targets;
          McsIdx.Atoms = fs.MoleculeFragment.Atoms;
          McsIdx.Bonds = fs.MoleculeFragment.Bonds;
        }
      }
      if (!McsIdx.QueryMolecule && QueryMoleculeSingleMatchedAtom) {
        McsIdx.QueryMolecule = QueryMolecule;
        McsIdx.Targets = Targets;
        McsIdx.Atoms =
            std::vector<const Atom *>{QueryMoleculeSingleMatchedAtom};
        McsIdx.Bonds = std::vector<const Bond *>();
      }
    } else if (i + 1 < endIdx) {
      Seed seed;
      if (createSeedFromMCS(i, seed)) {  // MCS is matched with new query
        Seeds.push_back(seed);
      }
      std::swap(
          Molecules.at(startIdx),
          Molecules.at(i + 1));  // change query molecule for threshold < 1.
    }
  }

  res.NumAtoms = getMaxNumberAtoms();
  if (!res.NumAtoms && QueryMoleculeSingleMatchedAtom) {
    res.NumAtoms = 1;
  }
  res.NumBonds = getMaxNumberBonds();

  if (res.NumBonds > 0 || QueryMoleculeSingleMatchedAtom) {
    if (!Parameters.StoreAll) {
      auto smartsQueryMolPair = generateResultSMARTSAndQueryMol(McsIdx);
      res.SmartsString = std::move(smartsQueryMolPair.first);
      res.QueryMol = std::move(smartsQueryMolPair.second);
    } else {
      std::transform(DegenerateMcsMap.begin(), DegenerateMcsMap.end(),
                     std::inserter(res.DegenerateSmartsQueryMolDict,
                                   res.DegenerateSmartsQueryMolDict.end()),
                     [this](const auto &pair) {
                       return generateResultSMARTSAndQueryMol(pair.second);
                     });
    }
  }

#ifdef VERBOSE_STATISTICS_ON
  if (Parameters.Verbose && res.NumAtoms > 0) {
    for (const auto &tag : Targets) {
      unsigned int itarget = &tag - &Targets.front();
      MatchVectType match;

      bool target_matched = SubstructMatch(*tag.Molecule, *res.QueryMol, match);
      if (!target_matched) {
        std::cout << "Target " << itarget + 1
                  << (target_matched ? " matched " : " MISMATCHED ")
                  << MolToSmiles(*tag.Molecule) << "\n";
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
    std::cout << "IndividualBondExcluded   = "
              << VerboseStatistics.IndividualBondExcluded << "\n";
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

  auto pos = faked_ring_info.find_first();
  while (pos != boost::dynamic_bitset<>::npos) {
    src_mols[pos]->getRingInfo()->reset();
    pos = faked_ring_info.find_next(pos);
  }

  clear();
  return res;
}

bool MaximumCommonSubgraph::checkIfMatchAndAppend(Seed &seed) {
#ifdef VERBOSE_STATISTICS_ON
  ++VerboseStatistics.SeedCheck;
#endif
#ifdef FAST_SUBSTRUCT_CACHE
  SubstructureCache::HashKey cacheKey;
  SubstructureCache::TIndexEntry *cacheEntry = nullptr;
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
      if (cacheEntry) {  // possibly found. check for hash collision
#ifdef VERBOSE_STATISTICS_ON
        ++VerboseStatistics.HashKeyFoundInCache;
#endif
        // check hash collisions (time +3%):
        for (const auto &g : *cacheEntry) {
          if (g.m_vertices.size() != seed.getNumAtoms() ||
              g.m_edges.size() != seed.getNumBonds()) {
            continue;
          }
#ifdef VERBOSE_STATISTICS_ON
          ++VerboseStatistics.ExactMatchCall;
#endif
          // EXACT MATCH
          foundInCache = SubstructMatchCustomTable(
              g, *QueryMolecule, seed.Topology, *QueryMolecule,
              QueryAtomMatchTable, QueryBondMatchTable, &Parameters);
#ifdef VERBOSE_STATISTICS_ON
          if (foundInCache) {
            ++VerboseStatistics.ExactMatchCallTrue;
          } else {
            break;
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

  Seed *newSeed = nullptr;

  {
    if (found) {  // Store new generated seed, if found in cache or in
                  // all(- threshold) targets
      newSeed = &Seeds.add(seed);
      newSeed->CopyComplete = false;

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

bool MaximumCommonSubgraph::match(Seed &seed) {
  unsigned int max_miss = Targets.size() - ThresholdCount;
  unsigned int missing = 0;
  unsigned int passed = 0;

  for (const auto &tag : Targets) {
    unsigned int itarget = &tag - &Targets.front();
#ifdef VERBOSE_STATISTICS_ON
    {
      ++VerboseStatistics.MatchCall;
    }
#endif
    bool target_matched = false;
    if (!seed.MatchResult.empty() && !seed.MatchResult.at(itarget).empty()) {
      target_matched = matchIncrementalFast(seed, itarget);
    }
    if (!target_matched) {  // slow full match
      match_V_t match;      // THERE IS NO Bonds match INFO !!!!
      target_matched = SubstructMatchCustomTable(
          tag.Topology, *tag.Molecule, seed.Topology, *QueryMolecule,
          tag.AtomMatchTable, tag.BondMatchTable, &Parameters, &match);
      // save current match info
      if (target_matched) {
        if (seed.MatchResult.empty()) {
          seed.MatchResult.resize(Targets.size());
        }
        seed.MatchResult[itarget].init(seed, match, *QueryMolecule, tag);
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
bool MaximumCommonSubgraph::matchIncrementalFast(Seed &seed,
                                                 unsigned int itarget) {
// use and update results of previous match stored in the seed
#ifdef VERBOSE_STATISTICS_ON
  {
    ++VerboseStatistics.FastMatchCall;
  }
#endif
  const auto &target = Targets.at(itarget);
  auto &match = seed.MatchResult.at(itarget);
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
  for (unsigned int newBondSeedIdx = match.MatchedBondSize;
       newBondSeedIdx < seed.getNumBonds(); newBondSeedIdx++) {
    matched = false;
    bool atomAdded = false;
    const auto newBond = seed.MoleculeFragment.Bonds.at(newBondSeedIdx);
    unsigned int newBondQueryIdx = newBond->getIdx();

    // seed's index of atom from which new bond was added
    unsigned int newBondSourceAtomSeedIdx;
    // seed's index of atom on other end of the bond
    unsigned int newBondOtherAtomSeedIdx;
    unsigned int i =
        seed.MoleculeFragment.SeedAtomIdxMap.at(newBond->getBeginAtomIdx());
    unsigned int j =
        seed.MoleculeFragment.SeedAtomIdxMap.at(newBond->getEndAtomIdx());
    if (i >= match.MatchedAtomSize) {
      // this is new atom in the seed
      newBondSourceAtomSeedIdx = j;
      newBondOtherAtomSeedIdx = i;
    } else {
      newBondSourceAtomSeedIdx = i;
      newBondOtherAtomSeedIdx = j;
    }
    unsigned int newBondOtherAtomQueryIdx =
        seed.MoleculeFragment.Atoms.at(newBondOtherAtomSeedIdx)->getIdx();
    unsigned int newBondSourceAtomQueryIdx =
        seed.MoleculeFragment.Atoms.at(newBondSourceAtomSeedIdx)->getIdx();
    // matched to newBondSourceAtomSeedIdx
    unsigned int newBondSourceAtomTargetIdx =
        match.TargetAtomIdx.at(newBondSourceAtomQueryIdx);
    const Bond *tb = nullptr;
    unsigned int newBondOtherAtomTargetIdx = NotSet;

    if (newBondOtherAtomSeedIdx < match.MatchedAtomSize) {
      // new bond between old atoms - both are
      // matched to exact atoms in the target
      newBondOtherAtomTargetIdx =
          match.TargetAtomIdx.at(newBondOtherAtomQueryIdx);
      // target bond between Source and Other atom
      tb = target.Molecule->getBondBetweenAtoms(newBondSourceAtomTargetIdx,
                                                newBondOtherAtomTargetIdx);
      if (tb) {
        // bond exists, check match with query molecule
        unsigned int tbi = tb->getIdx();
        unsigned int qbi =
            seed.MoleculeFragment.Bonds.at(newBondSeedIdx)->getIdx();
        if (!match.VisitedTargetBonds.test(tbi)) {
          // false if target bond is already matched
          matched = target.BondMatchTable.at(qbi, tbi);
        }
      }
    } else {
      // enumerate all bonds from source atom in the target
      const auto atom =
          target.Molecule->getAtomWithIdx(newBondSourceAtomTargetIdx);
      for (const auto &nbri :
           boost::make_iterator_range(target.Molecule->getAtomBonds(atom))) {
        tb = (*target.Molecule)[nbri];
        if (match.VisitedTargetBonds.test(tb->getIdx())) {
          continue;
        }
        newBondOtherAtomTargetIdx = tb->getBeginAtomIdx();
        if (newBondSourceAtomTargetIdx == newBondOtherAtomTargetIdx) {
          newBondOtherAtomTargetIdx = tb->getEndAtomIdx();
        }
        if (match.VisitedTargetAtoms.test(newBondOtherAtomTargetIdx)) {
          continue;
        }
        // check OtherAtom and bond
        matched = target.AtomMatchTable.at(newBondOtherAtomQueryIdx,
                                           newBondOtherAtomTargetIdx) &&
                  target.BondMatchTable.at(
                      seed.MoleculeFragment.Bonds.at(newBondSeedIdx)->getIdx(),
                      tb->getIdx());
        if (matched) {
          atomAdded = true;
          break;
        }
      }
    }

    if (matched) {      // update match history
      if (atomAdded) {  // new atom has been added
        match.MatchedAtomSize++;
        match.TargetAtomIdx[newBondOtherAtomQueryIdx] =
            newBondOtherAtomTargetIdx;
        match.VisitedTargetAtoms.set(newBondOtherAtomTargetIdx);
      }
      match.MatchedBondSize++;
      match.TargetBondIdx[newBondQueryIdx] = tb->getIdx();
      match.VisitedTargetBonds.set(tb->getIdx());
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
    c1.reserve(seed.getNumAtoms());
    std::vector<std::uint32_t> c2;
    c2.reserve(target.Molecule->getNumAtoms());
    for (unsigned int si = 0; si < seed.getNumAtoms(); ++si) {
      // index in the seed topology
      c1.push_back(si);
      c2.push_back(match.TargetAtomIdx.at(seed.Topology[si]));
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
