//
//  Copyright (C) 2018-2021 Susan H. Leung and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Tautomer.h"
#include "Fragment.h"
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentCatalogUtils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>
#include <limits>

#include <utility>

// #define VERBOSE_ENUMERATION 1

#ifdef VERBOSE_ENUMERATION
#include <GraphMol/SmilesParse/SmartsWrite.h>
#endif

namespace RDKit {

namespace MolStandardize {

namespace TautomerScoringFunctions {
int scoreRings(const ROMol &mol) {
  int score = 0;
  auto ringInfo = mol.getRingInfo();
  // For aromatic ring scoring, we only need SSSR (not symmetrized SSSR).
  // Aromatic rings are always essential to the ring system and will be in SSSR.
  // Avoiding symmetrizeSSSR and the molecule copy saves significant time
  // since scoreRings is called for each tautomer during canonicalization.
  if (!ringInfo->isSssrOrBetter()) {
    MolOps::findSSSR(mol);
    // ringInfo pointer remains valid - just refresh from the mol
    ringInfo = mol.getRingInfo();
  }
  boost::dynamic_bitset<> isArom(mol.getNumBonds());
  boost::dynamic_bitset<> bothCarbon(mol.getNumBonds());
  for (const auto &bnd : mol.bonds()) {
    if (bnd->getIsAromatic()) {
      isArom.set(bnd->getIdx());
      if (bnd->getBeginAtom()->getAtomicNum() == 6 &&
          bnd->getEndAtom()->getAtomicNum() == 6) {
        bothCarbon.set(bnd->getIdx());
      }
    }
  }
  for (const auto &bring : ringInfo->bondRings()) {
    bool allC = true;
    bool allAromatic = true;
    for (const auto bidx : bring) {
      if (!isArom[bidx]) {
        allAromatic = false;
        break;
      }
      if (!bothCarbon[bidx]) {
        allC = false;
      }
    }
    if (allAromatic) {
      score += 100;
      if (allC) {
        score += 150;
      }
    }
  }
  return score;
};

SubstructTerm::SubstructTerm(std::string aname, std::string asmarts, int ascore,
                             std::vector<int> reqElements,
                             std::string connSmarts)
    : name(std::move(aname)),
      smarts(std::move(asmarts)),
      score(ascore),
      requiredElements(std::move(reqElements)),
      connectivitySmarts(std::move(connSmarts)) {
  std::unique_ptr<ROMol> pattern(SmartsToMol(smarts));
  if (pattern) {
    matcher = std::move(*pattern);
  }
  // Initialize connectivity matcher if provided
  if (!connectivitySmarts.empty()) {
    std::unique_ptr<ROMol> connPattern(SmartsToMol(connectivitySmarts));
    if (connPattern) {
      connectivityMatcher = std::move(*connPattern);
    }
  }
}

const std::vector<SubstructTerm> &getDefaultTautomerScoreSubstructs() {
  // Each term specifies:
  //   - name, SMARTS, score
  //   - requiredElements: atomic numbers that must be present (empty = no filter)
  //   - connectivitySmarts: bond-order-agnostic pattern for pre-screening
  // Since tautomerization only moves H and changes bond orders (never creates/
  // destroys heavy-atom bonds), we can skip patterns whose connectivity
  // prerequisites aren't met by the input molecule.
  static std::vector<SubstructTerm> substructureTerms{
      {"benzoquinone", "[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]", 25, {6}, "[#6]1(~[#6]~[#6]~[#6](~[#6]~[#6]~1)~[N,S,O])~[N,S,O]"},
      {"oxim", "[#6]=[N][OH]", 4, {6, 7, 8}, "[#6]~[#7]~[#8]"},
      {"C=O", "[#6]=,:[#8]", 2, {6, 8}, "[#6]~[#8]"},
      {"N=O", "[#7]=,:[#8]", 2, {7, 8}, "[#7]~[#8]"},
      {"P=O", "[#15]=,:[#8]", 2, {15, 8}, "[#15]~[#8]"},
      {"C=hetero", "[C]=[!#1;!#6]", 1, {6}, "[C]~[!#1;!#6]"},
      {"C(=hetero)-hetero", "[C](=[!#1;!#6])[!#1;!#6]", 2, {6}, "[C](~[!#1;!#6])~[!#1;!#6]"},
      {"aromatic C = exocyclic N", "[c]=!@[N]", -1, {6, 7}, "[c]~[N]"},
      {"methyl", "[CX4H3]", 1, {6}, ""},
      {"guanidine terminal=N", "[#7]C(=[NR0])[#7H0]", 1, {6, 7}, "[#7]~[#6]~[#7]"},
      {"guanidine endocyclic=N", "[#7;R][#6;R]([N])=[#7;R]", 2, {6, 7}, "[#7]~[#6]~[#7]"},
      {"aci-nitro", "[#6]=[N+]([O-])[OH]", -4, {6, 7, 8}, "[#6]~[#7]~[#8]"}};
  return substructureTerms;
}

std::vector<size_t> getRelevantSubstructTermIndices(
    const ROMol &mol, const std::vector<SubstructTerm> &terms) {
  // Prepare relevant SubstructTerms for this molecule by filtering in two
  // stages:
  //   1. Element check: skip terms requiring elements not in the molecule
  //   2. Connectivity check: skip terms whose bond-order-agnostic pattern
  //      doesn't match (since tautomerization doesn't create/destroy bonds)

  std::unordered_set<int> presentElements;
  for (const auto atom : mol.atoms()) {
    presentElements.insert(atom->getAtomicNum());
  }

  std::vector<size_t> relevantIndices;
  relevantIndices.reserve(terms.size());

  SubstructMatchParameters params;
  params.maxMatches = 1;

  for (size_t i = 0; i < terms.size(); ++i) {
    const auto &term = terms[i];

    bool hasAllElements = true;
    for (int elem : term.requiredElements) {
      if (presentElements.find(elem) == presentElements.end()) {
        hasAllElements = false;
        break;
      }
    }
    if (!hasAllElements) {
      continue;
    }
    if (term.connectivityMatcher.getNumAtoms() > 0) {
      auto matches = SubstructMatch(mol, term.connectivityMatcher, params);
      if (matches.empty()) {
        continue;
      }
    }
    relevantIndices.push_back(i);
  }

  return relevantIndices;
}

int scoreSubstructsFiltered(const ROMol &mol,
                            const std::vector<SubstructTerm> &terms,
                            const std::vector<size_t> &relevantIndices) {
  int score = 0;
  SubstructMatchParameters params;
  for (size_t idx : relevantIndices) {
    const auto &term = terms[idx];
    if (!term.matcher.getNumAtoms()) {
      continue;
    }
    const auto nMatches = SubstructMatchCount(mol, term.matcher, params);
    score += static_cast<int>(nMatches) * term.score;
  }
  return score;
}

int scoreSubstructs(const ROMol &mol,
                    const std::vector<SubstructTerm> &substructureTerms) {
  int score = 0;
  SubstructMatchParameters params;
  for (const auto &term : substructureTerms) {
    if (!term.matcher.getNumAtoms()) {
      BOOST_LOG(rdErrorLog) << " matcher for term " << term.name
                            << " is invalid, ignoring it." << std::endl;
      continue;
    }
    const auto nMatches = SubstructMatchCount(mol, term.matcher, params);
    score += static_cast<int>(nMatches) * term.score;
  }
  return score;
}

int scoreHeteroHs(const ROMol &mol) {
  int score = 0;
  for (const auto &at : mol.atoms()) {
    int anum = at->getAtomicNum();
    if (anum == 15 || anum == 16 || anum == 34 || anum == 52) {
      score -= at->getTotalNumHs();
    }
  }
  return score;
}
}  // namespace TautomerScoringFunctions

TautomerEnumerator::TautomerEnumerator(const CleanupParameters &params)
    : d_maxTautomers(params.maxTautomers),
      d_maxTransforms(params.maxTransforms),
      d_removeSp3Stereo(params.tautomerRemoveSp3Stereo),
      d_removeBondStereo(params.tautomerRemoveBondStereo),
      d_removeIsotopicHs(params.tautomerRemoveIsotopicHs),
      d_reassignStereo(params.tautomerReassignStereo) {
  std::unique_ptr<TautomerCatalogParams> tautParams;
  if (params.tautomerTransformData.empty()) {
    tautParams.reset(new TautomerCatalogParams(params.tautomerTransforms));
  } else {
    tautParams.reset(new TautomerCatalogParams(params.tautomerTransformData));
  }
  dp_catalog.reset(new TautomerCatalog(tautParams.get()));
}

bool TautomerEnumerator::setTautomerStereoAndIsoHs(
    const ROMol &mol, ROMol &taut, const TautomerEnumeratorResult &res) const {
  bool modified = false;
  // Iterate only the atoms/bonds actually modified by transforms.
  for (auto atomIdx = res.d_modifiedAtoms.find_first();
       atomIdx != boost::dynamic_bitset<>::npos;
       atomIdx = res.d_modifiedAtoms.find_next(atomIdx)) {
    const auto atom = mol.getAtomWithIdx(static_cast<unsigned int>(atomIdx));
    auto tautAtom = taut.getAtomWithIdx(atomIdx);
    // clear chiral tag on sp2 atoms (also sp3 if d_removeSp3Stereo is true)
    if (tautAtom->getHybridization() == Atom::SP2 || d_removeSp3Stereo) {
      modified |= (tautAtom->getChiralTag() != Atom::CHI_UNSPECIFIED);
      tautAtom->setChiralTag(Atom::CHI_UNSPECIFIED);
      if (tautAtom->hasProp(common_properties::_CIPCode)) {
        tautAtom->clearProp(common_properties::_CIPCode);
      }
    } else {
      modified |= (tautAtom->getChiralTag() != atom->getChiralTag());
      tautAtom->setChiralTag(atom->getChiralTag());
      if (atom->hasProp(common_properties::_CIPCode)) {
        tautAtom->setProp(
            common_properties::_CIPCode,
            atom->getProp<std::string>(common_properties::_CIPCode));
      }
    }
    // remove isotopic Hs if present (and if d_removeIsotopicHs is true)
    if (tautAtom->hasProp(common_properties::_isotopicHs) &&
        (d_removeIsotopicHs || !tautAtom->getTotalNumHs())) {
      tautAtom->clearProp(common_properties::_isotopicHs);
    }
  }
  // remove stereochemistry on bonds that are part of a tautomeric path
  for (auto bondIdx = res.d_modifiedBonds.find_first();
       bondIdx != boost::dynamic_bitset<>::npos;
       bondIdx = res.d_modifiedBonds.find_next(bondIdx)) {
    const auto bond = mol.getBondWithIdx(static_cast<unsigned int>(bondIdx));
    std::vector<unsigned int> bondsToClearDirs;
    if (bond->getBondType() == Bond::DOUBLE &&
        bond->getStereo() > Bond::STEREOANY) {
      // look around the beginning and end atoms and check for bonds with
      // direction set
      for (auto atom : {bond->getBeginAtom(), bond->getEndAtom()}) {
        for (const auto &nbri :
             boost::make_iterator_range(mol.getAtomBonds(atom))) {
          const auto &obnd = mol[nbri];
          if (obnd->getBondDir() == Bond::ENDDOWNRIGHT ||
              obnd->getBondDir() == Bond::ENDUPRIGHT) {
            bondsToClearDirs.push_back(obnd->getIdx());
          }
        }
      }
    }
    auto tautBond = taut.getBondWithIdx(static_cast<unsigned int>(bondIdx));
    if (tautBond->getBondType() != Bond::DOUBLE || d_removeBondStereo) {
      // When bond stereo is being removed for bonds involved in tautomerism,
      // use STEREOANY (for *non-ring* double bonds) instead of STEREONONE.
      // This prevents downstream tools (notably InChI) from inferring a specific E/Z
      // assignment from 2D coordinates after bond orders have been changed.
      bool isRingBond = false;
      if (tautBond->getBondType() == Bond::DOUBLE) {
        auto ringInfo = taut.getRingInfo();
        // We only need to know whether the bond is in a ring; avoid forcing
        // SymmSSSR here since that has a performance cost.
        // this might be not precise enough for large rings like 9-or-more macrocycles
        if (!ringInfo || !ringInfo->isFindFastOrBetter()) {
          MolOps::fastFindRings(taut);
          ringInfo = taut.getRingInfo();
        }
        if (ringInfo) {
          const auto tbidx = tautBond->getIdx();
          const bool bondInRing = ringInfo->numBondRings(tbidx) != 0;
          const bool beginInRing =
              ringInfo->numAtomRings(tautBond->getBeginAtomIdx()) != 0;
          const bool endInRing =
              ringInfo->numAtomRings(tautBond->getEndAtomIdx()) != 0;
          // Only set STEREOANY on clearly non-ring bonds. Using atom-in-ring as
          // a fallback avoids any edge cases where bond indices or ring data
          // aren't aligned during enumeration.
          isRingBond = bondInRing || (beginInRing && endInRing);
        }
      }
      const auto targetStereo = (tautBond->getBondType() == Bond::DOUBLE &&
                                 !isRingBond)
                                    ? Bond::STEREOANY
                                    : Bond::STEREONONE;
      modified |= (tautBond->getStereo() != targetStereo);
      tautBond->setStereo(targetStereo);
      tautBond->getStereoAtoms().clear();
      for (auto bi : bondsToClearDirs) {
        taut.getBondWithIdx(bi)->setBondDir(Bond::NONE);
      }
    } else {
      const INT_VECT &sa = bond->getStereoAtoms();
      modified |= (tautBond->getStereo() != bond->getStereo() ||
                   sa.size() != tautBond->getStereoAtoms().size());
      if (sa.size() == 2) {
        tautBond->setStereoAtoms(sa.front(), sa.back());
      }
      tautBond->setStereo(bond->getStereo());
      for (auto bi : bondsToClearDirs) {
        taut.getBondWithIdx(bi)->setBondDir(
            mol.getBondWithIdx(bi)->getBondDir());
      }
    }
  }
  if (d_reassignStereo) {
    static const bool cleanIt = true;
    static const bool force = true;
    MolOps::assignStereochemistry(taut, cleanIt, force);

    // assignStereochemistry() can overwrite the explicit "undefined" bond
    // stereo (STEREOANY) that we set above in order to prevent downstream
    // coordinate-based E/Z inference. If bond stereo removal is enabled,
    // re-apply our contract to the bonds involved in tautomerism.
    if (d_removeBondStereo) {
      auto ringInfo = taut.getRingInfo();
      if (!ringInfo || !ringInfo->isFindFastOrBetter()) {
        // might prefer more expensive calc here to better support 9-or-more macrocycles
        MolOps::fastFindRings(taut);
        ringInfo = taut.getRingInfo();
      }
      for (auto bond : taut.bonds()) {
        const auto bondIdx = bond->getIdx();
        if (!res.d_modifiedBonds.test(bondIdx)) {
          continue;
        }
        if (bond->getBondType() != Bond::DOUBLE) {
          bond->setStereo(Bond::STEREONONE);
          bond->getStereoAtoms().clear();
          continue;
        }

        const bool bondInRing = ringInfo && ringInfo->numBondRings(bondIdx);
        const bool beginInRing =
            ringInfo && ringInfo->numAtomRings(bond->getBeginAtomIdx());
        const bool endInRing =
            ringInfo && ringInfo->numAtomRings(bond->getEndAtomIdx());
        const bool isRingBond = bondInRing || (beginInRing && endInRing);

        bond->setStereo(isRingBond ? Bond::STEREONONE : Bond::STEREOANY);
        bond->getStereoAtoms().clear();
      }
    }
  } else {
    taut.setProp(common_properties::_StereochemDone, 1);
  }
  return modified;
}

std::vector<ROMOL_SPTR> TautomerEnumerator::enumerate(
    const ROMol &mol, boost::dynamic_bitset<> *modifiedAtoms,
    boost::dynamic_bitset<> *modifiedBonds) const {
  TautomerEnumeratorResult tresult = enumerate(mol);
  if (modifiedAtoms) {
    *modifiedAtoms = tresult.modifiedAtoms();
  }
  if (modifiedBonds) {
    *modifiedBonds = tresult.modifiedBonds();
  }
  return tresult.tautomers();
}

TautomerEnumeratorResult TautomerEnumerator::enumerate(const ROMol &mol) const {
#ifdef VERBOSE_ENUMERATION
  std::cout << "**********************************" << std::endl;
#endif
  PRECONDITION(dp_catalog, "no catalog!");
  const TautomerCatalogParams *tautparams = dp_catalog->getCatalogParams();
  PRECONDITION(tautparams, "");

  TautomerEnumeratorResult res;

  const std::vector<TautomerTransform> &transforms =
      tautparams->getTransforms();

  // Enumerate all possible tautomers and return them as a vector.
  // smi is the input molecule SMILES
  std::string smi = MolToSmiles(mol, true);
  // taut is a copy of the input molecule
  ROMOL_SPTR taut(new ROMol(mol));
  // do whatever sanitization bits are required
  if (taut->needsUpdatePropertyCache()) {
    taut->updatePropertyCache(false);
  }
  if (!taut->getRingInfo()->isSymmSssr()) {
    MolOps::symmetrizeSSSR(*taut);
  }

  // Kekulized form will be created lazily when needed for transform matching
  res.d_tautomers = {{smi, Tautomer(taut, 0, 0)}};
  res.d_modifiedAtoms.resize(mol.getNumAtoms());
  res.d_modifiedBonds.resize(mol.getNumBonds());

  // Keep running counts of modified atoms/bonds.
  // `boost::dynamic_bitset<>::count()` is O(n) in the number of blocks, and we
  // were previously calling it once per new tautomer, which is avoidable.
  size_t numModifiedAtoms = 0;
  size_t numModifiedBonds = 0;
  const auto markAtomModified = [&res, &numModifiedAtoms](unsigned int idx) {
    if (!res.d_modifiedAtoms.test(idx)) {
      res.d_modifiedAtoms.set(idx);
      ++numModifiedAtoms;
    }
  };
  const auto markBondModified = [&res, &numModifiedBonds](unsigned int idx) {
    if (!res.d_modifiedBonds.test(idx)) {
      res.d_modifiedBonds.set(idx);
      ++numModifiedBonds;
    }
  };
  bool completed = false;
  bool bailOut = false;
  unsigned int nTransforms = 0;
  static const std::array<const char *, 4> statusMsg{
      "completed", "max tautomers reached", "max transforms reached",
      "canceled"};

  while (!completed && !bailOut) {
    // std::map automatically sorts res.d_tautomers into alphabetical order
    // (SMILES)
    for (auto &smilesTautomerPair : res.d_tautomers) {
#ifdef VERBOSE_ENUMERATION
      std::cout << "Current tautomers: " << std::endl;
      for (const auto &smilesTautomerPair : res.d_tautomers) {
        std::cout << smilesTautomerPair.first << " done "
                  << smilesTautomerPair.second.d_done << std::endl;
      }
#endif
      std::string tsmiles;
      if (smilesTautomerPair.second.d_done) {
#ifdef VERBOSE_ENUMERATION
        std::cout << "Skipping " << smilesTautomerPair.first
                  << " as already done" << std::endl;
#endif
        continue;
      }
#ifdef VERBOSE_ENUMERATION
      std::cout << "Looking at tautomer: " << smilesTautomerPair.first
                << std::endl;
#endif
      // tautomer not yet done
      for (const auto &transform : transforms) {
        if (bailOut) {
          break;
        }
        // kmol is the kekulized version of the tautomer (created lazily)
        const auto &kmol = smilesTautomerPair.second.getKekulized();
        std::vector<MatchVectType> matches;
        unsigned int matched = SubstructMatch(*kmol, *(transform.Mol), matches);

        if (!matched) {
          continue;
        }
        ++nTransforms;
#ifdef VERBOSE_ENUMERATION
        std::string name;
        (transform.Mol)->getProp(common_properties::_Name, name);
        SmilesWriteParams smilesWriteParams;
        smilesWriteParams.allBondsExplicit = true;
        std::cout << "kmol for " << smilesTautomerPair.first << " : "
                  << MolToSmiles(*kmol, smilesWriteParams) << std::endl;
        std::cout << "transform mol: " << MolToSmarts(*(transform.Mol))
                  << std::endl;

        std::cout << "Matched: " << name << std::endl;
#endif
        // loop over transform matches
        for (const auto &match : matches) {
          if (nTransforms >= d_maxTransforms) {
            res.d_status = TautomerEnumeratorStatus::MaxTransformsReached;
            bailOut = true;
          } else if (res.d_tautomers.size() >= d_maxTautomers) {
            res.d_status = TautomerEnumeratorStatus::MaxTautomersReached;
            bailOut = true;
          } else if (d_callback.get() && !(*d_callback)(mol, res)) {
            res.d_status = TautomerEnumeratorStatus::Canceled;
            bailOut = true;
          }
          if (bailOut) {
            break;
          }
          // Create a copy of in the input molecule so we can modify it
          // Use kekule form so bonds are explicitly single/double instead of
          // aromatic
          RWMOL_SPTR product(new RWMol(*kmol, true));
          // Remove a hydrogen from the first matched atom and add one to the
          // last
          int firstIdx = match.front().second;
          int lastIdx = match.back().second;
          Atom *first = product->getAtomWithIdx(firstIdx);
          Atom *last = product->getAtomWithIdx(lastIdx);
          markAtomModified(static_cast<unsigned int>(firstIdx));
          markAtomModified(static_cast<unsigned int>(lastIdx));
          first->setNumExplicitHs(
              std::max(0, static_cast<int>(first->getTotalNumHs()) - 1));
          last->setNumExplicitHs(last->getTotalNumHs() + 1);
          // Remove any implicit hydrogens from the first and last atoms
          // now we have set the count explicitly
          first->setNoImplicit(true);
          last->setNoImplicit(true);
          // Adjust bond orders
          unsigned int bi = 0;
          for (size_t i = 0; i < transform.Mol->getNumBonds(); ++i) {
            const auto tbond = transform.Mol->getBondWithIdx(i);
            Bond *bond = product->getBondBetweenAtoms(
                match[tbond->getBeginAtomIdx()].second,
                match[tbond->getEndAtomIdx()].second);
            ASSERT_INVARIANT(bond, "required bond not found");
            // check if bonds is specified in tautomer.in file
            if (!transform.BondTypes.empty()) {
              bond->setBondType(transform.BondTypes[bi]);
              ++bi;
            } else {
              Bond::BondType bondtype = bond->getBondType();
#ifdef VERBOSE_ENUMERATION
              std::cout << "Bond as double: " << bond->getBondTypeAsDouble()
                        << std::endl;
              std::cout << bondtype << std::endl;
#endif
              if (bondtype == Bond::SINGLE) {
                bond->setBondType(Bond::DOUBLE);
#ifdef VERBOSE_ENUMERATION
                std::cout << "Set bond to double" << std::endl;
#endif
              }
              if (bondtype == Bond::DOUBLE) {
                bond->setBondType(Bond::SINGLE);
#ifdef VERBOSE_ENUMERATION
                std::cout << "Set bond to single" << std::endl;
#endif
              }
            }
            markBondModified(bond->getIdx());
          }
          // TODO adjust charges
          if (!transform.Charges.empty()) {
            unsigned int ci = 0;
            for (const auto &pair : match) {
              Atom *atom = product->getAtomWithIdx(pair.second);
              atom->setFormalCharge(atom->getFormalCharge() +
                                    transform.Charges[ci++]);
            }
          }
#ifdef VERBOSE_ENUMERATION
          {
            SmilesWriteParams smilesWriteParams;
            smilesWriteParams.allBondsExplicit = true;
            std::cout << "pre-sanitize: "
                      << MolToSmiles(*product, smilesWriteParams) << std::endl;
          }
#endif

          try {
            // We only change bond orders/H counts/charges; the molecular graph
            // (and therefore ring topology) is unchanged.
            // `sanitizeMol()` always calls `clearComputedProps()` which resets
            // ring info and forces ring-finding for each generated tautomer.
            // Avoid that by clearing computed props without touching rings,
            // then running the specific sanitize steps we need.
            product->clearComputedProps(false);
            product->updatePropertyCache(false);
            MolOps::Kekulize(*product);
            MolOps::setAromaticity(*product);
            MolOps::setConjugation(*product);
            MolOps::setHybridization(*product);
            MolOps::adjustHs(*product);
          } catch (const KekulizeException &) {
            continue;
          }
#ifdef VERBOSE_ENUMERATION
          SmilesWriteParams smilesWriteParams;
          smilesWriteParams.allBondsExplicit = true;
          std::cout << "pre-setTautomerStereo: "
                    << MolToSmiles(*product, smilesWriteParams) << std::endl;
#endif
          setTautomerStereoAndIsoHs(mol, *product, res);
          tsmiles = MolToSmiles(*product, true);
#ifdef VERBOSE_ENUMERATION
          (transform.Mol)->getProp(common_properties::_Name, name);
          std::cout << "Applied rule: " << name << " to "
                    << smilesTautomerPair.first << std::endl;
#endif
          if (res.d_tautomers.find(tsmiles) != res.d_tautomers.end()) {
#ifdef VERBOSE_ENUMERATION
            std::cout << "Previous tautomer produced again: " << tsmiles
                      << std::endl;
#endif
            continue;
          }
          // in addition to the above transformations, sanitization may modify
          // bonds, e.g. Cc1nc2ccccc2[nH]1
          for (size_t i = 0; i < mol.getNumBonds(); i++) {
            auto molBondType = mol.getBondWithIdx(i)->getBondType();
            auto tautBondType = product->getBondWithIdx(i)->getBondType();
            if (molBondType != tautBondType && !res.d_modifiedBonds.test(i)) {
#ifdef VERBOSE_ENUMERATION
              std::cout << "Sanitization has modified bond " << i << std::endl;
#endif
              markBondModified(static_cast<unsigned int>(i));
            }
          }
          // Kekulized form will be created lazily when needed
#ifdef VERBOSE_ENUMERATION
          auto it = res.d_tautomers.find(tsmiles);
          if (it == res.d_tautomers.end()) {
            std::cout << "New tautomer added as ";
          } else {
            std::cout << "New tautomer replaced for ";
          }
          std::cout << tsmiles << ", taut: " << MolToSmiles(*product)
                    << std::endl;
#endif
          // BOOST_LOG(rdInfoLog)
          //     << "Tautomer transform "
          //     <<
          //     transform.Mol->getProp<std::string>(common_properties::_Name)
          //     << " produced tautomer " << tsmiles << std::endl;
          res.d_tautomers[tsmiles] = Tautomer(
              std::move(product),
              numModifiedAtoms, numModifiedBonds);
        }
      }
      smilesTautomerPair.second.d_done = true;
    }
    completed = true;
    size_t maxNumModifiedAtoms = numModifiedAtoms;
    size_t maxNumModifiedBonds = numModifiedBonds;
    for (auto it = res.d_tautomers.begin(); it != res.d_tautomers.end();) {
      auto &taut = it->second;
      if (!taut.d_done) {
        completed = false;
      }
      if ((taut.d_numModifiedAtoms < maxNumModifiedAtoms ||
           taut.d_numModifiedBonds < maxNumModifiedBonds) &&
          setTautomerStereoAndIsoHs(mol, *taut.tautomer, res)) {
        Tautomer tautStored = std::move(taut);
        it = res.d_tautomers.erase(it);
        tautStored.d_numModifiedAtoms = maxNumModifiedAtoms;
        tautStored.d_numModifiedBonds = maxNumModifiedBonds;
        auto insertRes = res.d_tautomers.insert(std::make_pair(
            MolToSmiles(*tautStored.tautomer), std::move(tautStored)));
        if (insertRes.second) {
          it = insertRes.first;
        }
      } else {
        ++it;
      }
    }
    if (bailOut && res.d_tautomers.size() < d_maxTautomers &&
        res.d_status == TautomerEnumeratorStatus::MaxTautomersReached) {
      res.d_status = TautomerEnumeratorStatus::Completed;
      bailOut = false;
    }
  }  // while
  res.fillTautomersItVec();
  if (!completed) {
    BOOST_LOG(rdWarningLog)
        << "Tautomer enumeration stopped at " << res.d_tautomers.size()
        << " tautomers: " << statusMsg.at(static_cast<size_t>(res.d_status))
        << std::endl;
  }

  return res;
}

// pickCanonical non-templated overload that avoids recomputing SMILES
ROMol *TautomerEnumerator::pickCanonical(
    const TautomerEnumeratorResult &tautRes,
    boost::function<int(const ROMol &mol)> scoreFunc) const {
  ROMOL_SPTR bestMol;
  if (tautRes.d_tautomers.size() == 1) {
    bestMol = tautRes.d_tautomers.begin()->second.tautomer;
  } else {
    // Calculate score for each tautomer
    int bestScore = std::numeric_limits<int>::min();
    std::string bestSmiles = "";
    for (const auto &t : tautRes.d_tautomers) {
      auto score = scoreFunc(*t.second.tautomer);
#ifdef VERBOSE_ENUMERATION
      std::cerr << "  " << t.first << " " << score << std::endl;
#endif
      if (score > bestScore) {
        bestScore = score;
        bestSmiles = t.first;
        bestMol = t.second.tautomer;
      } else if (score == bestScore) {
        if (t.first < bestSmiles) {
          bestSmiles = t.first;
          bestMol = t.second.tautomer;
        }
      }
    }
  }
  ROMol *res = new ROMol(*bestMol);
  static const bool cleanIt = true;
  static const bool force = true;
  MolOps::assignStereochemistry(*res, cleanIt, force);

  return res;
}

ROMol *TautomerEnumerator::canonicalize(
    const ROMol &mol, boost::function<int(const ROMol &mol)> scoreFunc) const {
  auto thisCopy = TautomerEnumerator(*this);
  thisCopy.setReassignStereo(false);
  auto res = thisCopy.enumerate(mol);
  if (res.empty()) {
    BOOST_LOG(rdWarningLog)
        << "no tautomers found, returning input molecule" << std::endl;
    return new ROMol(mol);
  }
  // When no custom scorer provided, use optimized scoring that pre-filters
  // SubstructTerm patterns once for the input molecule rather than evaluating
  // all 12 substructure matches per tautomer. This is safe because
  // tautomerization only moves H and changes bond orders, never creates or
  // destroys heavy-atom bonds.
  if (!scoreFunc) {
    scoreFunc = TautomerScoringFunctions::makeOptimizedScorer(mol);
  }
  return pickCanonical(res, scoreFunc);
}

void TautomerEnumerator::canonicalizeInPlace(
    RWMol &mol, boost::function<int(const ROMol &mol)> scoreFunc) const {
  auto thisCopy = TautomerEnumerator(*this);
  thisCopy.setReassignStereo(false);
  auto res = thisCopy.enumerate(mol);
  if (res.empty()) {
    BOOST_LOG(rdWarningLog)
        << "no tautomers found, molecule unchanged" << std::endl;
    return;
  }
  // When no custom scorer provided, use optimized scoring that pre-filters
  // SubstructTerm patterns once for the input molecule.
  if (!scoreFunc) {
    scoreFunc = TautomerScoringFunctions::makeOptimizedScorer(mol);
  }
  std::unique_ptr<ROMol> tmp{pickCanonical(res, scoreFunc)};

  TEST_ASSERT(tmp->getNumAtoms() == mol.getNumAtoms());
  TEST_ASSERT(tmp->getNumBonds() == mol.getNumBonds());
  // now copy the info from the canonical tautomer over to the input molecule
  for (const auto tmpAtom : tmp->atoms()) {
    auto atom = mol.getAtomWithIdx(tmpAtom->getIdx());
    TEST_ASSERT(tmpAtom->getAtomicNum() == atom->getAtomicNum());
    atom->setFormalCharge(tmpAtom->getFormalCharge());
    atom->setNoImplicit(tmpAtom->getNoImplicit());
    atom->setIsAromatic(tmpAtom->getIsAromatic());
    atom->setNumExplicitHs(tmpAtom->getNumExplicitHs());
    atom->setNumRadicalElectrons(tmpAtom->getNumRadicalElectrons());
    atom->setChiralTag(tmpAtom->getChiralTag());
  }
  for (const auto tmpBond : tmp->bonds()) {
    auto bond = mol.getBondWithIdx(tmpBond->getIdx());
    TEST_ASSERT(tmpBond->getBeginAtomIdx() == bond->getBeginAtomIdx());
    TEST_ASSERT(tmpBond->getEndAtomIdx() == bond->getEndAtomIdx());
    bond->setBondType(tmpBond->getBondType());
    bond->setBondDir(tmpBond->getBondDir());
    bond->setIsAromatic(tmpBond->getIsAromatic());
    bond->setIsConjugated(tmpBond->getIsConjugated());
    if (tmpBond->getStereoAtoms().size() == 2) {
      bond->setStereoAtoms(tmpBond->getStereoAtoms()[0],
                           tmpBond->getStereoAtoms()[1]);
    }
    bond->setStereo(tmpBond->getStereo());
  }
  mol.updatePropertyCache(false);
}

}  // namespace MolStandardize
}  // namespace RDKit
