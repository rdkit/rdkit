//
//  Copyright (C) 2018-2020 Susan H. Leung and Greg Landrum
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

#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <utility>

// #define VERBOSE_ENUMERATION 1

#ifdef VERBOSE_ENUMERATION
#include <GraphMol/SmilesParse/SmartsWrite.h>
#endif

using namespace RDKit;

namespace RDKit {

namespace MolStandardize {

namespace TautomerScoringFunctions {
int scoreRings(const ROMol &mol) {
  int score = 0;
  auto ringInfo = mol.getRingInfo();
  std::unique_ptr<ROMol> cp;
  if (!ringInfo->isInitialized()) {
    cp.reset(new ROMol(mol));
    MolOps::symmetrizeSSSR(*cp);
    ringInfo = cp->getRingInfo();
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

struct smarts_mol_holder {
  std::string d_smarts;
  ROMOL_SPTR dp_mol;
  smarts_mol_holder(const std::string &smarts) : d_smarts(smarts) {
    dp_mol.reset(SmartsToMol(smarts));
  }
};

typedef boost::flyweight<
    boost::flyweights::key_value<std::string, smarts_mol_holder>,
    boost::flyweights::no_tracking>
    smarts_mol_flyweight;

struct SubstructTerm {
  std::string name;
  std::string smarts;
  int score;
  ROMOL_SPTR matcher;
  SubstructTerm(std::string aname, std::string asmarts, int ascore)
      : name(std::move(aname)), smarts(std::move(asmarts)), score(ascore) {
    matcher = smarts_mol_flyweight(smarts).get().dp_mol;
  };
};

int scoreSubstructs(const ROMol &mol) {
  // a note on efficiency here: we'll construct the SubstructTerm objects here
  // repeatedly, but the SMARTS parsing for each entry will only be done once
  // since we're using the boost::flyweights above to cache them
  const std::vector<SubstructTerm> substructureTerms{
      {"benzoquinone", "[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]",
       25},
      {"oxim", "[#6]=[N][OH]", 4},
      {"C=O", "[#6]=,:[#8]", 2},
      {"N=O", "[#7]=,:[#8]", 2},
      {"P=O", "[#15]=,:[#8]", 2},
      {"C=hetero", "[C]=[!#1;!#6]", 1},
      {"aromatic C = exocyclic N", "[c]=!@[N]", -1},
      {"methyl", "[CX4H3]", 1},
      {"guanidine terminal=N", "[#7]C(=[NR0])[#7H0]", 1},
      {"guanidine endocyclic=N", "[#7;R][#6;R]([N])=[#7;R]", 2},
      {"aci-nitro", "[#6]=[N+]([O-])[OH]", -4}};
  int score = 0;
  for (const auto &term : substructureTerms) {
    if (!term.matcher) {
      BOOST_LOG(rdErrorLog) << " matcher for term " << term.name
                            << " is invalid, ignoring it." << std::endl;
      continue;
    }
    SubstructMatchParameters params;
    const auto matches = SubstructMatch(mol, *term.matcher, params);
    // if (!matches.empty()) {
    //   std::cerr << " " << matches.size() << " matches to " << term.name
    //             << std::endl;
    // }
    score += static_cast<int>(matches.size()) * term.score;
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
      d_reassignStereo(params.tautomerReassignStereo) {
  TautomerCatalogParams tautParams(params.tautomerTransforms);
  dp_catalog.reset(new TautomerCatalog(&tautParams));
}

bool TautomerEnumerator::setTautomerStereo(
    const ROMol &mol, ROMol &taut, const TautomerEnumeratorResult &res) const {
  bool modified = false;
  for (auto atom : mol.atoms()) {
    auto atomIdx = atom->getIdx();
    if (!res.d_modifiedAtoms.test(atomIdx)) {
      continue;
    }
    auto tautAtom = taut.getAtomWithIdx(atomIdx);
    // clear chiral tag on sp2 atoms
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
  }
  // remove stereochemistry on bonds that are part of a tautomeric path
  for (auto bond : mol.bonds()) {
    auto bondIdx = bond->getIdx();
    if (!res.d_modifiedBonds.test(bondIdx)) {
      continue;
    }
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
    auto tautBond = taut.getBondWithIdx(bondIdx);
    if (tautBond->getBondType() != Bond::DOUBLE || d_removeBondStereo) {
      modified |= (tautBond->getStereo() != Bond::STEREONONE);
      tautBond->setStereo(Bond::STEREONONE);
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
  // Create a kekulized form of the molecule to match the SMARTS against
  RWMOL_SPTR kekulized(new RWMol(mol));
  MolOps::Kekulize(*kekulized, false);
  res.d_tautomers = {{smi, Tautomer(taut, kekulized, 0, 0)}};
  std::set<std::string> done;
  res.d_modifiedAtoms.resize(mol.getNumAtoms());
  res.d_modifiedBonds.resize(mol.getNumBonds());
  bool completed = false;
  bool bailOut = false;
  unsigned int nTransforms = 0;
  static const std::array<const char *, 4> statusMsg{
      "completed", "max tautomers reached", "max transforms reached",
      "canceled"};

  while (!completed && !bailOut) {
    // std::map automatically sorts res.d_tautomers into alphabetical order
    // (SMILES)
    for (const auto &tautomer : res.d_tautomers) {
#ifdef VERBOSE_ENUMERATION
      std::cout << "Done : " << std::endl;
      for (const auto d : done) {
        std::cout << d << std::endl;
      }
      std::cout << "Looking at tautomer: " << tautomer.first << std::endl;
#endif
      std::string tsmiles;
      if (done.count(tautomer.first)) {
        continue;
      } else {
        // done does not contain tautomer
        for (const auto &transform : transforms) {
          if (bailOut) {
            break;
          }
          // kmol is the kekulized version of the tautomer
          const auto &kmol = tautomer.second.kekulized;
          std::vector<MatchVectType> matches;
          unsigned int matched =
              SubstructMatch(*kmol, *(transform.Mol), matches);

          if (!matched) {
            continue;
          }
          ++nTransforms;
#ifdef VERBOSE_ENUMERATION
          std::string name;
          (transform.Mol)->getProp(common_properties::_Name, name);
          std::cout << "kmol: " << kmol->first << std::endl;
          std::cout << MolToSmiles(*(kmol->second)) << std::endl;
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
            ROMOL_SPTR product(new ROMol(*kmol));
            // Remove a hydrogen from the first matched atom and add one to the
            // last
            int firstIdx = match.front().second;
            int lastIdx = match.back().second;
            Atom *first = product->getAtomWithIdx(firstIdx);
            Atom *last = product->getAtomWithIdx(lastIdx);
            res.d_modifiedAtoms.set(firstIdx);
            res.d_modifiedAtoms.set(lastIdx);
            first->setNumExplicitHs(std::max(0U, first->getTotalNumHs() - 1));
            last->setNumExplicitHs(last->getTotalNumHs() + 1);
            // Remove any implicit hydrogens from the first and last atoms
            // now we have set the count explicitly
            first->setNoImplicit(true);
            last->setNoImplicit(true);
            // Adjust bond orders
            unsigned int bi = 0;
            for (size_t i = 0; i < match.size() - 1; ++i) {
              Bond *bond = product->getBondBetweenAtoms(match[i].second,
                                                        match[i + 1].second);
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
              res.d_modifiedBonds.set(bond->getIdx());
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

            RWMOL_SPTR wproduct(new RWMol(*product));
#ifdef VERBOSE_ENUMERATION
            std::cout << "pre-sanitization: "
                      << MolToSmiles(*wproduct, true, true) << std::endl;
#endif
            unsigned int failedOp;
            MolOps::sanitizeMol(*wproduct, failedOp,
                                MolOps::SANITIZE_KEKULIZE |
                                    MolOps::SANITIZE_SETAROMATICITY |
                                    MolOps::SANITIZE_SETCONJUGATION |
                                    MolOps::SANITIZE_SETHYBRIDIZATION |
                                    MolOps::SANITIZE_ADJUSTHS);
            setTautomerStereo(mol, *wproduct, res);
            tsmiles = MolToSmiles(*wproduct, true);
#ifdef VERBOSE_ENUMERATION
            std::string name;
            (transform.Mol)->getProp(common_properties::_Name, name);
            std::cout << "Applied rule: " << name << " to " << tautomer.first
                      << std::endl;
#endif
            if (res.d_tautomers.find(tsmiles) != res.d_tautomers.end()) {
#ifdef VERBOSE_ENUMERATION
              std::cout << "Previous tautomer produced again: " << tsmiles
                        << std::endl;
#endif
              continue;
            }
#ifdef VERBOSE_ENUMERATION
            std::cout << "New tautomer produced: " << tsmiles << std::endl;
#endif
            // in addition to the above transformations, sanitzation may modify
            // bonds, e.g. Cc1nc2ccccc2[nH]1
            for (size_t i = 0; i < mol.getNumBonds(); i++) {
              auto molBondType = mol.getBondWithIdx(i)->getBondType();
              auto tautBondType = wproduct->getBondWithIdx(i)->getBondType();
              if (molBondType != tautBondType && !res.d_modifiedBonds.test(i)) {
#ifdef VERBOSE_ENUMERATION
                std::cout << "Sanitization has modified bond " << i
                          << std::endl;
#endif
                res.d_modifiedBonds.set(i);
              }
            }
            RWMOL_SPTR kekulized_product(new RWMol(*wproduct));
            MolOps::Kekulize(*kekulized_product, false);
            res.d_tautomers[tsmiles] = Tautomer(
                std::move(wproduct), std::move(kekulized_product),
                res.d_modifiedAtoms.count(), res.d_modifiedBonds.count());

#ifdef VERBOSE_ENUMERATION
            std::cout << "Now completed: " << std::endl;
            for (const auto &tautomer : res.d_tautomers) {
              std::cout << tautomer.first << std::endl;
            }
#endif
          }
        }
      }
      done.insert(tautomer.first);
    }
    completed = (res.d_tautomers.size() <= done.size());
    size_t maxNumModifiedAtoms = res.d_modifiedAtoms.count();
    size_t maxNumModifiedBonds = res.d_modifiedBonds.count();
    for (auto it = res.d_tautomers.begin(); it != res.d_tautomers.end();) {
      auto &taut = it->second;
      if ((taut.d_numModifiedAtoms < maxNumModifiedAtoms ||
           taut.d_numModifiedBonds < maxNumModifiedBonds) &&
          setTautomerStereo(mol, *taut.tautomer, res)) {
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
  return pickCanonical(res, scoreFunc);
}

}  // namespace MolStandardize
}  // namespace RDKit
