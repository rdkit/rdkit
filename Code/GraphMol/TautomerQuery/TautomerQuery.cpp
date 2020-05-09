//
// Created by Gareth Jones on 5/7/2020.
//
// Copyright 2020 Schrodinger, Inc
//

#include "TautomerQuery.h"
#include <boost/smart_ptr.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/Tautomer.h>
#include <GraphMol/Bond.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/Substruct/SubstructUtils.h>
#include <GraphMol/Fingerprints/Fingerprints.h>

namespace {

int getTargetIdx(int queryIdx, const RDKit::MatchVectType &match) {
  const auto pair = match[queryIdx];
  // is this always true ?
  assert(pair.first == static_cast<int>(queryIdx));
  return pair.second;
}

}  // namespace

namespace RDKit {

TautomerQuery::TautomerQuery(const ROMol &query,
                             const std::vector<ROMOL_SPTR> tautomers,
                             const ROMOL_SPTR templateMolecule,
                             const std::vector<size_t> modifiedAtoms,
                             const std::vector<size_t> modifiedBonds)
    : query(query),
      tautomers(tautomers),
      templateMolecule(templateMolecule),
      modifiedAtoms(modifiedAtoms),
      modifiedBonds(modifiedBonds) {}

TautomerQuery::~TautomerQuery() {}

boost::shared_ptr<TautomerQuery> TautomerQuery::fromMol(
    const ROMol &query, const std::string &tautomerTransformFile) {
  auto tautomerFile = !tautomerTransformFile.empty()
                          ? tautomerTransformFile
                          : std::string(getenv("RDBASE")) +
                                "/Data/MolStandardize/tautomerTransforms.in";
  auto tautomerParams = std::unique_ptr<MolStandardize::TautomerCatalogParams>(
      new MolStandardize::TautomerCatalogParams(tautomerFile));
  MolStandardize::TautomerEnumerator tautomerEnumerator(
      new MolStandardize::TautomerCatalog(tautomerParams.get()));
  boost::dynamic_bitset<> modifiedAtomsBitSet(query.getNumAtoms());
  boost::dynamic_bitset<> modifiedBondsBitSet(query.getNumBonds());
  const auto t2 = tautomerEnumerator.enumerate(query);
  const auto tautomers = tautomerEnumerator.enumerate(
      query, &modifiedAtomsBitSet, &modifiedBondsBitSet);

  std::vector<size_t> modifiedAtoms;
  modifiedAtoms.reserve(modifiedAtomsBitSet.count());
  for (size_t i = 0; i < query.getNumAtoms(); i++) {
    if (modifiedAtomsBitSet[i]) {
      modifiedAtoms.push_back(i);
    }
  }
  std::vector<size_t> modifiedBonds;
  modifiedBonds.reserve(modifiedBondsBitSet.count());
  for (size_t i = 0; i < query.getNumBonds(); i++) {
    if (modifiedBondsBitSet[i]) {
      modifiedBonds.push_back(i);
    }
  }

  auto templateMolecule = boost::make_shared<RWMol>(query);
  for (auto idx : modifiedAtoms) {
    const auto atom = templateMolecule->getAtomWithIdx(idx);
    // I think this should match on atom number only (no formal charge
    // matching), and also fingerprint properly
    const auto queryAtom = new QueryAtom(*atom);
    // auto atomQuery = makeAtomNumQuery(atom->getAtomicNum());
    // queryAtom->setQuery(atomQuery);
    const bool updateLabel = false;
    const bool preserveProps = true;
    templateMolecule->replaceAtom(idx, queryAtom, updateLabel, preserveProps);
    delete queryAtom;
  }
  for (auto idx : modifiedBonds) {
    auto bondQuery1 = makeSingleOrAromaticBondQuery();
    auto bondQuery2 = makeBondOrderEqualsQuery(Bond::DOUBLE);
    auto queryBond = new QueryBond();
    queryBond->setQuery(bondQuery1);
    queryBond->expandQuery(bondQuery2, Queries::COMPOSITE_OR);
    templateMolecule->replaceBond(idx, queryBond, true);
    delete queryBond;
  }

  boost::shared_ptr<TautomerQuery> tautomerQuery(new TautomerQuery(
      query, tautomers, boost::make_shared<ROMol>(*templateMolecule),
      modifiedAtoms, modifiedBonds));
  return tautomerQuery;
}

bool TautomerQuery::matchTautomer(
    const ROMol &mol, const ROMol &tautomer, const MatchVectType &match,
    const SubstructMatchParameters &params) const {
  for (auto idx : modifiedAtoms) {
    const auto queryAtom = tautomer.getAtomWithIdx(idx);
    const auto targetAtom = mol.getAtomWithIdx(getTargetIdx(idx, match));
    if (!atomCompat(queryAtom, targetAtom, params)) return false;
  }

  for (auto idx : modifiedBonds) {
    const auto queryBond = tautomer.getBondWithIdx(idx);
    const auto beginIdx = queryBond->getBeginAtomIdx();
    const auto endIdx = queryBond->getEndAtomIdx();
    const auto targetBeginIdx = getTargetIdx(beginIdx, match);
    const auto targetEndIdx = getTargetIdx(endIdx, match);
    const auto targetBond =
        mol.getBondBetweenAtoms(targetBeginIdx, targetEndIdx);
    if (!bondCompat(queryBond, targetBond, params)) {
      return false;
    }
  }

  return true;
}

std::vector<MatchVectType> TautomerQuery::SubstructMatch(
    const ROMol &mol, const SubstructMatchParameters &params,
    std::vector<ROMOL_SPTR> *matchingTautomers) const {
  if (matchingTautomers) {
    matchingTautomers->clear();
  }
  std::vector<MatchVectType> matches;

  const auto templateMatches =
      RDKit::SubstructMatch(mol, *templateMolecule, params);

  for (auto templateMatch : templateMatches) {
    for (auto tautomer : tautomers) {
      if (matchTautomer(mol, *tautomer, templateMatch, params)) {
        matches.push_back(templateMatch);
        if (matchingTautomers) {
          matchingTautomers->push_back(tautomer);
        }
        break;
      }
    }
  }

  return matches;
}

ExplicitBitVect *TautomerQuery::patternFingerprintTemplate(uint fpSize) {
  return PatternFingerprintMol(*templateMolecule, fpSize, nullptr, nullptr,
                               true);
}

std::vector<MatchVectType> SubstructMatch(
    const ROMol &mol, const TautomerQuery &query,
    const SubstructMatchParameters &params) {
  return query.SubstructMatch(mol, params);
}

}  // namespace RDKit
