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
  return pair.second;
}

}  // namespace

namespace RDKit {

TautomerQuery::TautomerQuery(const std::vector<ROMOL_SPTR> &tautomers,
                             const ROMol *const templateMolecule,
                             const std::vector<size_t> &modifiedAtoms,
                             const std::vector<size_t> &modifiedBonds)
    : d_tautomers(tautomers),
      d_templateMolecule(templateMolecule),
      d_modifiedAtoms(modifiedAtoms),
      d_modifiedBonds(modifiedBonds) {}

TautomerQuery::~TautomerQuery() {
   delete d_templateMolecule;
}

TautomerQuery *TautomerQuery::fromMol(
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

  auto templateMolecule = new RWMol(query);
  for (auto idx : modifiedAtoms) {
    const auto atom = templateMolecule->getAtomWithIdx(idx);
    const auto queryAtom = new QueryAtom(atom->getAtomicNum());
    const bool updateLabel = false;
    const bool preserveProps = true;
    templateMolecule->replaceAtom(idx, queryAtom, updateLabel, preserveProps);
    delete queryAtom;
  }
  for (auto idx : modifiedBonds) {
    auto bondQuery = makeSingleOrDoubleOrAromaticBondQuery();
    auto queryBond = new QueryBond();
    queryBond->setQuery(bondQuery);
    templateMolecule->replaceBond(idx, queryBond, true);
    delete queryBond;
  }

  return new TautomerQuery(tautomers, (ROMol *)templateMolecule, modifiedAtoms,
                           modifiedBonds);
}

bool TautomerQuery::matchTautomer(
    const ROMol &mol, const ROMol &tautomer, const MatchVectType &match,
    const SubstructMatchParameters &params) const {
  for (auto idx : d_modifiedAtoms) {
    const auto queryAtom = tautomer.getAtomWithIdx(idx);
    const auto targetAtom = mol.getAtomWithIdx(getTargetIdx(idx, match));
    if (!atomCompat(queryAtom, targetAtom, params)) return false;
  }

  for (auto idx : d_modifiedBonds) {
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

std::vector<MatchVectType> TautomerQuery::substructOf(
    const ROMol &mol, const SubstructMatchParameters &params,
    std::vector<ROMOL_SPTR> *matchingTautomers) const {
  if (matchingTautomers) {
    matchingTautomers->clear();
  }
  std::vector<MatchVectType> matches;

  const auto templateMatches =
      RDKit::SubstructMatch(mol, *d_templateMolecule, params);

  for (auto templateMatch : templateMatches) {
    for (auto tautomer : d_tautomers) {
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

bool TautomerQuery::isSubstructOf(const ROMol &mol,
                                  const SubstructMatchParameters &params) {
  auto matches = substructOf(mol, params);
  return matches.size() > 0;
}

ExplicitBitVect *TautomerQuery::patternFingerprintTemplate(
    unsigned int fpSize) {
  return PatternFingerprintMol(*d_templateMolecule, fpSize, nullptr, nullptr,
                               true);
}

ExplicitBitVect *TautomerQuery::patternFingerprintTarget(const ROMol &target,
                                                         unsigned int fpSize) {
  return PatternFingerprintMol(target, fpSize, nullptr, nullptr, true);
}

std::vector<MatchVectType> SubstructMatch(
    const ROMol &mol, const TautomerQuery &query,
    const SubstructMatchParameters &params) {
  return query.substructOf(mol, params);
}

}  // namespace RDKit