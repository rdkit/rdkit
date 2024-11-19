//
//  Copyright (C) 2013-2020 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/QueryOps.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include "Fingerprints.h"
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Invariant.h>
#include <boost/random.hpp>
#include <climits>
#include <cstdint>
#include <RDGeneral/hash/hash.hpp>
#include <RDGeneral/types.h>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <RDGeneral/BoostEndInclude.h>

// #define VERBOSE_FINGERPRINTING 1

namespace {
class ss_matcher {
 public:
  ss_matcher() {};
  ss_matcher(const std::string &pattern) {
    RDKit::RWMol *p = RDKit::SmartsToMol(pattern);
    TEST_ASSERT(p);
    m_matcher.reset(p);
  };

  // const RDKit::ROMOL_SPTR &getMatcher() const { return m_matcher; };
  const RDKit::ROMol *getMatcher() const { return m_matcher.get(); };

 private:
  RDKit::ROMOL_SPTR m_matcher;
};
}  // namespace

namespace RDKit {
const char *pqs[] = {
    "[*]~[*]", "[*]~[*]~[*]", "[R]~1~[R]~[R]~1",
    //"[*]~[*]~[*]~[*]",
    "[*]~[*](~[*])~[*]",
    //"[*]~[R]~1[R]~[R]~1",
    "[R]~1[R]~[R]~[R]~1",
    //"[*]~[*]~[*]~[*]~[*]",
    "[*]~[*]~[*](~[*])~[*]",
    //"[*]~[R]~1[R]~[R]~1~[*]",
    "[R]~1~[R]~[R]~[R]~[R]~1", "[R]~1~[R]~[R]~[R]~[R]~[R]~1",
    //"[R2]~[R1]~[R2]", Github #151: can't have ring counts in an SSS pattern
    //"[R2]~[R1]~[R1]~[R2]",  Github #151: can't have ring counts in an SSS
    // pattern
    "[R](@[R])(@[R])~[R]~[R](@[R])(@[R])",
    "[R](@[R])(@[R])~[R]@[R]~[R](@[R])(@[R])",

    //"[*]!@[R]~[R]!@[*]",  Github #151: can't have !@ in an SSS pattern
    //"[*]!@[R]~[R]~[R]!@[*]", Github #151: can't have !@ in an SSS pattern
    "[*]~[R](@[R])@[R](@[R])~[*]", "[*]~[R](@[R])@[R]@[R](@[R])~[*]",
    "[*]",  // single atom fragment
    ""};
typedef boost::flyweight<boost::flyweights::key_value<std::string, ss_matcher>,
                         boost::flyweights::no_tracking>
    pattern_flyweight;

namespace detail {
void getAtomNumbers(const Atom *a, std::vector<int> &atomNums) {
  atomNums.clear();
  if (!a->hasQuery()) {
    atomNums.push_back(a->getAtomicNum());
    return;
  }
  // negated things are always complex:
  if (a->getQuery()->getNegation()) {
    return;
  }
  std::string descr = a->getQuery()->getDescription();
  if (descr == "AtomAtomicNum") {
    atomNums.push_back(
        static_cast<ATOM_EQUALS_QUERY *>(a->getQuery())->getVal());
  } else if (descr == "AtomXor") {
    return;
  } else if (descr == "AtomAnd") {
    auto childIt = a->getQuery()->beginChildren();
    if ((*childIt)->getDescription() == "AtomAtomicNum" &&
        ((*(childIt + 1))->getDescription() == "AtomIsAliphatic" ||
         (*(childIt + 1))->getDescription() == "AtomIsAromatic") &&
        (childIt + 2) == a->getQuery()->endChildren()) {
      atomNums.push_back(
          static_cast<ATOM_EQUALS_QUERY *>((*childIt).get())->getVal());
      return;
    }
  } else if (descr == "AtomOr") {
    auto childIt = a->getQuery()->beginChildren();
    while (childIt != a->getQuery()->endChildren()) {
      if ((*childIt)->getDescription() == "AtomAtomicNum") {
        atomNums.push_back(
            static_cast<ATOM_EQUALS_QUERY *>((*childIt).get())->getVal());
      } else if ((*childIt)->getDescription() == "AtomAnd") {
        auto childIt2 = (*childIt)->beginChildren();
        if ((*childIt2)->getDescription() == "AtomAtomicNum" &&
            ((*(childIt2 + 1))->getDescription() == "AtomIsAliphatic" ||
             (*(childIt2 + 1))->getDescription() == "AtomIsAromatic") &&
            (childIt2 + 2) == (*childIt)->endChildren()) {
          atomNums.push_back(
              static_cast<ATOM_EQUALS_QUERY *>((*childIt2).get())->getVal());
        } else {
          atomNums.clear();
          return;
        }
      } else {
        atomNums.clear();
        return;
      }
      ++childIt;
    }
  }
  return;
}
}  // namespace detail

namespace {

bool isPatternComplexQuery(const Bond *b) {
  if (!b->hasQuery()) {
    return false;
  }
  // negated things are always complex:
  if (b->getQuery()->getNegation()) {
    return true;
  }
  std::string descr = b->getQuery()->getDescription();
  // std::cerr<<"   !!!!!! "<<b->getIdx()<<"
  // "<<b->getBeginAtomIdx()<<"-"<<b->getEndAtomIdx()<<" "<<descr<<std::endl;
  return descr != "BondOrder";
}

bool isTautomerBondQuery(const Bond *b) {
  // assumes we have already tested true for isPatternComplexQuery
  auto description = b->getQuery()->getDescription();
  return description == "SingleOrDoubleOrAromaticBond" ||
         description == "SingleOrAromaticBond";
}

void updatePatternFingerprint(const ROMol &mol, ExplicitBitVect &fp,
                              unsigned int fpSize,
                              std::vector<unsigned int> *atomCounts,
                              ExplicitBitVect *setOnlyBits,
                              bool tautomericFingerprint) {
  PRECONDITION(fpSize != 0, "fpSize==0");
  PRECONDITION(!atomCounts || atomCounts->size() >= mol.getNumAtoms(),
               "bad atomCounts size");
  PRECONDITION(!setOnlyBits || setOnlyBits->getNumBits() == fpSize,
               "bad setOnlyBits size");

  std::vector<const ROMol *> patts;
  patts.reserve(10);
  unsigned int idx = 0;
  while (1) {
    std::string pq = pqs[idx];
    if (pq == "") {
      break;
    }
    ++idx;
    const ROMol *matcher = pattern_flyweight(pq).get().getMatcher();
    CHECK_INVARIANT(matcher, "bad smarts");
    patts.push_back(matcher);
  }

  if (!mol.getRingInfo()->isFindFastOrBetter()) {
    MolOps::fastFindRings(mol);
  }

  boost::dynamic_bitset<> isQueryAtom(mol.getNumAtoms()),
      isQueryBond(mol.getNumBonds()), isTautomerBond(mol.getNumBonds());
  for (const auto at : mol.atoms()) {
    // isComplexQuery() no longer considers "AtomNull" to be complex, but for
    // the purposes of the pattern FP, it definitely needs to be treated as a
    // query feature.
    if (at->hasQuery() && (at->getQuery()->getDescription() == "AtomNull" ||
                           isComplexQuery(at))) {
      isQueryAtom.set(at->getIdx());
    }
  }

  for (const auto bond : mol.bonds()) {
    if (isPatternComplexQuery(bond)) {
      isQueryBond.set(bond->getIdx());
      if (tautomericFingerprint && isTautomerBondQuery(bond)) {
        isTautomerBond.set(bond->getIdx());
      }
    }
  }

  unsigned int pIdx = 0;
  for (const auto patt : patts) {
    ++pIdx;
    std::vector<MatchVectType> matches;
    // uniquify matches?
    //   time for 10K molecules w/ uniquify: 5.24s
    //   time for 10K molecules w/o uniquify: 4.87s

    SubstructMatchParameters params;
    params.uniquify = false;
    // raise maxMatches really high. This was the cause for github #2614.
    // if we end up with more matches than this, we're completely hosed: :-)
    params.maxMatches = 100000000;
    matches = SubstructMatch(mol, *patt, params);

    std::uint32_t mIdx = pIdx + patt->getNumAtoms() + patt->getNumBonds();
    for (const auto &mv : matches) {
#ifdef VERBOSE_FINGERPRINTING
      std::cerr << "\nPatt: " << pIdx << " | ";
#endif
      // collect bits counting the number of occurrences of the pattern:
      gboost::hash_combine(mIdx, 0xBEEF);
      fp.setBit(mIdx % fpSize);
#ifdef VERBOSE_FINGERPRINTING
      std::cerr << "count: " << mIdx % fpSize << " | ";
#endif

      bool isQuery = false;
      std::uint32_t bitId = pIdx;
      std::vector<unsigned int> amap(mv.size(), 0);
      for (const auto &p : mv) {
#ifdef VERBOSE_FINGERPRINTING
        std::cerr << p.second << " ";
#endif
        if (isQueryAtom[p.second]) {
          isQuery = true;
#ifdef VERBOSE_FINGERPRINTING
          std::cerr << "atom query.";
#endif
          break;
        }
        gboost::hash_combine(bitId,
                             mol.getAtomWithIdx(p.second)->getAtomicNum());
        amap[p.first] = p.second;
      }
      if (isQuery) {
        continue;
      }
      auto tautomerBitId = bitId;
      auto tautomerQuery = false;
      ROMol::EDGE_ITER firstB, lastB;
      boost::tie(firstB, lastB) = patt->getEdges();
#ifdef VERBOSE_FINGERPRINTING
      std::cerr << " bs:|| ";
#endif
      while (!isQuery && firstB != lastB) {
        const Bond *pbond = (*patt)[*firstB];
        ++firstB;
        const Bond *mbond = mol.getBondBetweenAtoms(
            amap[pbond->getBeginAtomIdx()], amap[pbond->getEndAtomIdx()]);
        const auto bondIdx = mbond->getIdx();

        if (isQueryBond[bondIdx]) {
          isQuery = true;
          if (isTautomerBond[bondIdx]) {
            isQuery = false;
            tautomerQuery = true;
#ifdef VERBOSE_FINGERPRINTING
            std::cerr << "tautomer query: " << mbond->getIdx();
#endif
          }
          if (isQuery) {
#ifdef VERBOSE_FINGERPRINTING
            std::cerr << "bond query: " << mbond->getIdx();
#endif
            break;
          }
        }

        if (tautomericFingerprint) {
          if (isTautomerBond[bondIdx] || mbond->getIsAromatic() ||
              mbond->getBondType() == Bond::SINGLE ||
              mbond->getBondType() == Bond::DOUBLE ||
              mbond->getBondType() == Bond::AROMATIC) {
            gboost::hash_combine(tautomerBitId, -1);
#ifdef VERBOSE_FINGERPRINTING
            std::cerr << "T ";
#endif
          }
        }

        if (!tautomerQuery) {
          if (!mbond->getIsAromatic()) {
            gboost::hash_combine(bitId, (std::uint32_t)mbond->getBondType());
#ifdef VERBOSE_FINGERPRINTING
            std::cerr << mbond->getBondType() << " ";
#endif
          } else {
            gboost::hash_combine(bitId, (std::uint32_t)Bond::AROMATIC);
#ifdef VERBOSE_FINGERPRINTING
            std::cerr << Bond::AROMATIC << " ";
#endif
          }
        }
      }

      if (!isQuery) {
        if (!tautomerQuery) {
#ifdef VERBOSE_FINGERPRINTING
          std::cerr << " set: " << bitId << " " << bitId % fpSize;
#endif
          fp.setBit(bitId % fpSize);
        }
        if (tautomericFingerprint) {
#ifdef VERBOSE_FINGERPRINTING
          std::cerr << " tset: " << tautomerBitId << " "
                    << tautomerBitId % fpSize;
#endif
          fp.setBit(tautomerBitId % fpSize);
        }
      }
    }
  }
}

}  // namespace

// caller owns the result, it must be deleted
ExplicitBitVect *PatternFingerprintMol(const ROMol &mol, unsigned int fpSize,
                                       std::vector<unsigned int> *atomCounts,
                                       ExplicitBitVect *setOnlyBits,
                                       bool tautomericFingerprint) {
  PRECONDITION(fpSize != 0, "fpSize==0");
  PRECONDITION(!atomCounts || atomCounts->size() >= mol.getNumAtoms(),
               "bad atomCounts size");
  PRECONDITION(!setOnlyBits || setOnlyBits->getNumBits() == fpSize,
               "bad setOnlyBits size");
  auto *res = new ExplicitBitVect(fpSize);
  updatePatternFingerprint(mol, *res, fpSize, atomCounts, setOnlyBits,
                           tautomericFingerprint);
  return res;
}

// caller owns the result, it must be deleted
ExplicitBitVect *PatternFingerprintMol(const MolBundle &bundle,
                                       unsigned int fpSize,
                                       ExplicitBitVect *setOnlyBits,
                                       bool tautomericFingerprint) {
  PRECONDITION(fpSize != 0, "fpSize==0");
  PRECONDITION(!setOnlyBits || setOnlyBits->getNumBits() == fpSize,
               "bad setOnlyBits size");
  ExplicitBitVect *res = nullptr;
  for (const auto &molp : bundle.getMols()) {
    ExplicitBitVect molfp(fpSize);
    updatePatternFingerprint(*molp, molfp, fpSize, nullptr, setOnlyBits,
                             tautomericFingerprint);
    if (!res) {
      res = new ExplicitBitVect(molfp);
    } else {
      (*res) &= molfp;
    }
  }
  return res;
}

}  // namespace RDKit
