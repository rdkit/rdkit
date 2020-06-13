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

// #define VERBOSE

#ifdef VERBOSE
#include <GraphMol/SmilesParse/SmilesWrite.h>
#endif

namespace {

using namespace RDKit;

int getTargetIdx(int queryIdx, const MatchVectType &match) {
  const auto pair = match[queryIdx];
  return pair.second;
}

// Adapted from Code/GraphMol/Substruct/SubstructUtils.cpp#removeDuplicates
void removeTautomerDuplicates(std::vector<MatchVectType> &matches,
                              std::vector<ROMOL_SPTR> *matchingTautomers,
                              int nAtoms) {
  //
  //  This works by tracking the indices of the atoms in each match vector.
  //  This can lead to unexpected behavior when looking at rings and queries
  //  that don't specify bond orders.  For example querying this molecule:
  //    C1CCC=1
  //  with the pattern constructed from SMARTS C~C~C~C will return a
  //  single match, despite the fact that there are 4 different paths
  //  when valence is considered.  The defense of this behavior is
  //  that the 4 paths are equivalent in the semantics of the query.
  //  Also, OELib returns the same results
  //

  std::vector<boost::dynamic_bitset<>> seen;
  std::vector<MatchVectType> res;
  for (size_t i = 0; i < matches.size(); i++) {
    auto match = matches[i];
    boost::dynamic_bitset<> val(nAtoms);
    for (const auto &ci : match) {
      val.set(ci.second);
    }
    if (std::find(seen.begin(), seen.end(), val) == seen.end()) {
      // it's something new
      res.push_back(match);
      seen.push_back(val);
    } else if (matchingTautomers) {
      int position = res.size();
      matchingTautomers->erase(matchingTautomers->begin() + position);
    }
  }

  matches = res;
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

TautomerQuery::~TautomerQuery() { delete d_templateMolecule; }

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
#ifdef VERBOSE
    std::cout << "Query atom " << queryAtom->getSymbol() << " target atom "
              << targetAtom->getSymbol() << std::endl;
#endif
    if (!atomCompat(queryAtom, targetAtom, params)) {
#ifdef VERBOSE
      std::cout << "Atom incompatibility" << std::endl;
#endif
      return false;
    }
  }

  for (auto idx : d_modifiedBonds) {
    const auto queryBond = tautomer.getBondWithIdx(idx);
    const auto beginIdx = queryBond->getBeginAtomIdx();
    const auto endIdx = queryBond->getEndAtomIdx();
    const auto targetBeginIdx = getTargetIdx(beginIdx, match);
    const auto targetEndIdx = getTargetIdx(endIdx, match);
    const auto targetBond =
        mol.getBondBetweenAtoms(targetBeginIdx, targetEndIdx);
#ifdef VERBOSE
    std::cout << "Query bond " << queryBond->getBondType() << " target bond "
              << targetBond->getBondType() << std::endl;
#endif
    if (!bondCompat(queryBond, targetBond, params)) {
#ifdef VERBOSE
      std::cout << "Bond incompatibility" << std::endl;
#endif
      return false;
    }
  }

#ifdef VERBOSE
  std::cout << "Tautomer match" << std::endl;
#endif
  return true;
}

std::vector<MatchVectType> TautomerQuery::substructOf(
    const ROMol &mol, const SubstructMatchParameters &params,
    std::vector<ROMOL_SPTR> *matchingTautomers) const {
  if (matchingTautomers) {
    matchingTautomers->clear();
  }
  std::vector<MatchVectType> matches;

#ifdef VERBOSE
  std::cout << "Tautomer search with query " << MolToSmiles(*d_templateMolecule)
            << " on " << MolToSmiles(mol) << " max matches "
            << params.maxMatches << std::endl;
#endif
  SubstructMatchParameters templateParams(params);
  templateParams.maxMatches = 1000;
  templateParams.uniquify = false;
  const auto templateMatches =
      RDKit::SubstructMatch(mol, *d_templateMolecule, templateParams);
#ifdef VERBOSE
  std::cout << "Number of template matches " << templateMatches.size()
            << std::endl;
#endif

  // TODO create a functor so that I dont have to get all template matches
  // before evaluating tautomer matches.
  for (auto templateMatch : templateMatches) {
#ifdef VERBOSE
    std::cout << "Checking template match" << std::endl;
#endif
    for (auto tautomer : d_tautomers) {
#ifdef VERBOSE
      std::cout << "Checking Tautomer " << MolToSmiles(*tautomer) << std::endl;
#endif
      if (matchTautomer(mol, *tautomer, templateMatch, params)) {
#ifdef VERBOSE
        std::cout << "Got Match " << std::endl;
#endif
        matches.push_back(templateMatch);
        if (matchingTautomers) {
          matchingTautomers->push_back(tautomer);
        }
        if (matches.size() == params.maxMatches) goto searchFinished;
        break;
      }
    }
  }

searchFinished:
#ifdef VERBOSE
  std::cout << "Found " << matches.size() << " matches " << std::endl;
#endif

  if (params.uniquify && matches.size() > 1) {
    removeTautomerDuplicates(matches, matchingTautomers, mol.getNumAtoms());
#ifdef VERBOSE
    std::cout << "After removing duplicates " << matches.size() << " matches " << std::endl;
#endif
  }
  return matches;
}

bool TautomerQuery::isSubstructOf(const ROMol &mol,
                                  const SubstructMatchParameters &params) {
  SubstructMatchParameters params2(params);
  params2.maxMatches = 1;
  auto matches = substructOf(mol, params2);
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