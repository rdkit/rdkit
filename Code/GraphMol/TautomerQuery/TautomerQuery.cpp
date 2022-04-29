//
// Created by Gareth Jones on 5/7/2020.
//
// Copyright 2020 Schrodinger, Inc
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include "TautomerQuery.h"
#include <functional>
#include <set>
#include <utility>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/Tautomer.h>
#include <GraphMol/Bond.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/Substruct/SubstructUtils.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
// #define VERBOSE

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/archive_exception.hpp>
#include <RDGeneral/BoostEndInclude.h>
#endif

#ifdef VERBOSE
#include <GraphMol/SmilesParse/SmilesWrite.h>
#endif

namespace {

using namespace RDKit;

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

  std::set<boost::dynamic_bitset<>> seen;
  std::vector<MatchVectType> res;
  res.reserve(matches.size());
  for (auto &&match : matches) {
    boost::dynamic_bitset<> val(nAtoms);
    for (const auto &ci : match) {
      val.set(ci.second);
    }
    auto pos = seen.lower_bound(val);
    if (pos == seen.end() || *pos != val) {
      res.push_back(std::move(match));
      seen.insert(pos, std::move(val));
    } else if (matchingTautomers) {
      auto position = res.size();
      matchingTautomers->erase(matchingTautomers->begin() + position);
    }
  }
  res.shrink_to_fit();
  matches = std::move(res);
}

}  // namespace

namespace RDKit {

bool TautomerQueryCanSerialize() {
#ifdef RDK_USE_BOOST_SERIALIZATION
  return true;
#else
  return false;
#endif
}

class TautomerQueryMatcher {
 private:
  const TautomerQuery &d_tautomerQuery;
  const SubstructMatchParameters &d_params;
  std::vector<ROMOL_SPTR> *d_matchingTautomers;

 public:
  TautomerQueryMatcher(const TautomerQuery &tautomerQuery,
                       const SubstructMatchParameters &params,
                       std::vector<ROMOL_SPTR> *matchingTautomers)
      : d_tautomerQuery(tautomerQuery),
        d_params(params),
        d_matchingTautomers(matchingTautomers) {}

  bool match(const ROMol &mol, const std::vector<unsigned int> &match) {
#ifdef VERBOSE
    std::cout << "Checking template match" << std::endl;
#endif

    for (auto tautomer : d_tautomerQuery.getTautomers()) {
#ifdef VERBOSE
      std::cout << "Checking Tautomer " << MolToSmiles(*tautomer) << std::endl;
#endif
      if (d_tautomerQuery.matchTautomer(mol, *tautomer, match, d_params)) {
        auto matchingTautomer = d_params.extraFinalCheck
                                    ? d_params.extraFinalCheck(mol, match)
                                    : true;
        if (matchingTautomer) {
#ifdef VERBOSE
          std::cout << "Got Match " << std::endl;
#endif
          if (d_matchingTautomers) {
            d_matchingTautomers->push_back(tautomer);
          }
        }
        return matchingTautomer;
      }
    }
    return false;
  }
};

TautomerQuery::TautomerQuery(std::vector<ROMOL_SPTR> tautomers,
                             const ROMol *const templateMolecule,
                             std::vector<size_t> modifiedAtoms,
                             std::vector<size_t> modifiedBonds)
    : d_tautomers(std::move(tautomers)),
      d_templateMolecule(templateMolecule),
      d_modifiedAtoms(std::move(modifiedAtoms)),
      d_modifiedBonds(std::move(modifiedBonds)) {}

TautomerQuery *TautomerQuery::fromMol(
    const ROMol &query, const std::string &tautomerTransformFile) {
  auto tautomerParams = std::unique_ptr<MolStandardize::TautomerCatalogParams>(
      new MolStandardize::TautomerCatalogParams(tautomerTransformFile));
  MolStandardize::TautomerEnumerator tautomerEnumerator(
      new MolStandardize::TautomerCatalog(tautomerParams.get()));
  const auto res = tautomerEnumerator.enumerate(query);

  std::vector<size_t> modifiedAtoms;
  modifiedAtoms.reserve(res.modifiedAtoms().count());
  for (size_t i = 0; i < query.getNumAtoms(); i++) {
    if (res.modifiedAtoms().test(i)) {
      modifiedAtoms.push_back(i);
    }
  }
  std::vector<size_t> modifiedBonds;
  modifiedBonds.reserve(res.modifiedBonds().count());
  for (size_t i = 0; i < query.getNumBonds(); i++) {
    if (res.modifiedBonds().test(i)) {
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

  return new TautomerQuery(res.tautomers(),
                           static_cast<ROMol *>(templateMolecule),
                           modifiedAtoms, modifiedBonds);
}

bool TautomerQuery::matchTautomer(
    const ROMol &mol, const ROMol &tautomer,
    const std::vector<unsigned int> &match,
    const SubstructMatchParameters &params) const {
  for (auto idx : d_modifiedAtoms) {
    const auto queryAtom = tautomer.getAtomWithIdx(idx);
    const auto targetAtom = mol.getAtomWithIdx(match[idx]);
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
    const auto targetBeginIdx = match[beginIdx];
    const auto targetEndIdx = match[endIdx];
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

#ifdef VERBOSE
  std::cout << "Tautomer search with query " << MolToSmiles(*d_templateMolecule)
            << " on " << MolToSmiles(mol) << " max matches "
            << params.maxMatches << std::endl;
#endif
  SubstructMatchParameters templateParams(params);
  // need to check all mappings of template to target
  templateParams.uniquify = false;

  TautomerQueryMatcher tautomerQueryMatcher(*this, params, matchingTautomers);
  // use this functor as a final check to see if any tautomer matches the target
  auto checker = [&tautomerQueryMatcher](
                     const ROMol &mol,
                     const std::vector<unsigned int> &match) mutable {
    return tautomerQueryMatcher.match(mol, match);
  };
  templateParams.extraFinalCheck = checker;

  auto matches =
      RDKit::SubstructMatch(mol, *d_templateMolecule, templateParams);

#ifdef VERBOSE
  std::cout << "Found " << matches.size() << " matches " << std::endl;
#endif

  if (params.uniquify && matches.size() > 1) {
    removeTautomerDuplicates(matches, matchingTautomers, mol.getNumAtoms());
#ifdef VERBOSE
    std::cout << "After removing duplicates " << matches.size() << " matches "
              << std::endl;
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
    unsigned int fpSize) const {
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

void TautomerQuery::toStream(std::ostream &ss) const {
#ifndef RDK_USE_BOOST_SERIALIZATION
  PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
  boost::archive::text_oarchive ar(ss);
  ar << *this;
#endif
}

std::string TautomerQuery::serialize() const {
  std::stringstream ss;
  toStream(ss);
  return ss.str();
}

void TautomerQuery::initFromStream(std::istream &ss) {
#ifndef RDK_USE_BOOST_SERIALIZATION
  PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
  boost::archive::text_iarchive ar(ss);
  ar >> *this;
#endif
}

void TautomerQuery::initFromString(const std::string &text) {
  std::stringstream ss(text);
  initFromStream(ss);
}

}  // namespace RDKit
