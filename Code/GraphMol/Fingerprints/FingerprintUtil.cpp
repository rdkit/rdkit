//
//  Copyright (C) 2018 Boran Adas, Google Summer of Code
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/FingerprintUtil.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <RDGeneral/hash/hash.hpp>
#include <boost/dynamic_bitset.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <GraphMol/QueryOps.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/random.hpp>
#include <cstdint>
#include <RDGeneral/BoostEndInclude.h>
#include <climits>
#include <RDGeneral/types.h>

namespace RDKit {
namespace AtomPairs {

std::uint32_t getAtomCode(const Atom *atom, unsigned int branchSubtract,
                          bool includeChirality) {
  PRECONDITION(atom, "no atom");
  std::uint32_t code;

  unsigned int numBranches = 0;
  if (atom->getDegree() > branchSubtract) {
    numBranches = atom->getDegree() - branchSubtract;
  }

  code = numBranches % maxNumBranches;
  unsigned int nPi = numPiElectrons(*atom) % maxNumPi;
  code |= nPi << numBranchBits;

  unsigned int typeIdx = 0;
  unsigned int nTypes = 1 << numTypeBits;
  while (typeIdx < nTypes) {
    if (atomNumberTypes[typeIdx] ==
        static_cast<unsigned int>(atom->getAtomicNum())) {
      break;
    } else if (atomNumberTypes[typeIdx] >
               static_cast<unsigned int>(atom->getAtomicNum())) {
      typeIdx = nTypes;
      break;
    }
    ++typeIdx;
  }
  if (typeIdx == nTypes) {
    --typeIdx;
  }
  code |= typeIdx << (numBranchBits + numPiBits);
  if (includeChirality) {
    std::string cipCode;
    if (atom->getPropIfPresent(common_properties::_CIPCode, cipCode)) {
      std::uint32_t offset = numBranchBits + numPiBits + numTypeBits;
      if (cipCode == "R") {
        code |= 1 << offset;
      } else if (cipCode == "S") {
        code |= 2 << offset;
      }
    }
  }
  POSTCONDITION(code < static_cast<std::uint32_t>(
                           1 << (codeSize + (includeChirality ? 2 : 0))),
                "code exceeds number of bits");
  return code;
};

std::uint32_t getAtomPairCode(std::uint32_t codeI, std::uint32_t codeJ,
                              unsigned int dist, bool includeChirality) {
  PRECONDITION(dist < maxPathLen, "dist too long");
  std::uint32_t res = dist;
  res |= std::min(codeI, codeJ) << numPathBits;
  res |= std::max(codeI, codeJ)
         << (numPathBits + codeSize + (includeChirality ? numChiralBits : 0));
  return res;
}

std::uint64_t getTopologicalTorsionCode(
    const std::vector<std::uint32_t> &pathCodes, bool includeChirality) {
  bool reverseIt = false;
  unsigned int i = 0;
  unsigned int j = pathCodes.size() - 1;
  while (i < j) {
    if (pathCodes[i] > pathCodes[j]) {
      reverseIt = true;
      break;
    } else if (pathCodes[i] < pathCodes[j]) {
      break;
    }
    ++i;
    --j;
  }

  int shiftSize = codeSize + (includeChirality ? numChiralBits : 0);
  std::uint64_t res = 0;
  if (reverseIt) {
    for (unsigned int i = 0; i < pathCodes.size(); ++i) {
      res |= static_cast<std::uint64_t>(pathCodes[pathCodes.size() - i - 1])
             << (shiftSize * i);
    }
  } else {
    for (unsigned int i = 0; i < pathCodes.size(); ++i) {
      res |= static_cast<std::uint64_t>(pathCodes[i]) << (shiftSize * i);
    }
  }
  return res;
}

std::uint32_t getTopologicalTorsionHash(
    const std::vector<std::uint32_t> &pathCodes) {
  bool reverseIt = false;
  unsigned int i = 0;
  unsigned int j = pathCodes.size() - 1;
  while (i < j) {
    if (pathCodes[i] > pathCodes[j]) {
      reverseIt = true;
      break;
    } else if (pathCodes[i] < pathCodes[j]) {
      break;
    }
    ++i;
    --j;
  }

  std::uint32_t res = 0;
  if (reverseIt) {
    for (unsigned int i = 0; i < pathCodes.size(); ++i) {
      gboost::hash_combine(res, pathCodes[pathCodes.size() - i - 1]);
    }
  } else {
    for (unsigned int pathCode : pathCodes) {
      gboost::hash_combine(res, pathCode);
    }
  }
  return res;
}
}  // namespace AtomPairs

namespace MorganFingerprints {

// Definitions for feature points adapted from:
// Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998)
const char *smartsPatterns[6] = {
    "[$([N;!H0;v3,v4&+1]),\
$([O,S;H1;+0]),\
n&H1&+0]",                                                  // Donor
    "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),\
$([O,S;H0;v2]),\
$([O,S;-]),\
$([N;v3;!$(N-*=[O,N,P,S])]),\
n&H0&+0,\
$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]",                    // Acceptor
    "[a]",                                                  // Aromatic
    "[F,Cl,Br,I]",                                          // Halogen
    "[#7;+,\
$([N;H2&+0][$([C,a]);!$([C,a](=O))]),\
$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),\
$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]",  // Basic
    "[$([C,S](=[O,S,P])-[O;H1,-1])]"                        // Acidic
};

const RDKit::ROMol *ss_matcher::getMatcher() const { return m_matcher.get(); }

ss_matcher::ss_matcher(){};
ss_matcher::ss_matcher(const std::string &pattern) {
  RDKit::RWMol *p = RDKit::SmartsToMol(pattern);
  TEST_ASSERT(p);
  m_matcher.reset(p);
};

typedef boost::flyweight<boost::flyweights::key_value<std::string, ss_matcher>,
                         boost::flyweights::no_tracking>
    pattern_flyweight;

std::vector<std::string> defaultFeatureSmarts(smartsPatterns,
                                              smartsPatterns + 6);
typedef boost::flyweight<boost::flyweights::key_value<std::string, ss_matcher>,
                         boost::flyweights::no_tracking>
    pattern_flyweight;
void getFeatureInvariants(const ROMol &mol, std::vector<uint32_t> &invars,
                          std::vector<const ROMol *> *patterns) {
  unsigned int nAtoms = mol.getNumAtoms();
  PRECONDITION(invars.size() >= nAtoms, "vector too small");

  std::vector<const ROMol *> featureMatchers;
  if (!patterns) {
    featureMatchers.reserve(defaultFeatureSmarts.size());
    for (std::vector<std::string>::const_iterator smaIt =
             defaultFeatureSmarts.begin();
         smaIt != defaultFeatureSmarts.end(); ++smaIt) {
      const ROMol *matcher = pattern_flyweight(*smaIt).get().getMatcher();
      CHECK_INVARIANT(matcher, "bad smarts");
      featureMatchers.push_back(matcher);
    }
    patterns = &featureMatchers;
  }
  std::fill(invars.begin(), invars.end(), 0);
  for (unsigned int i = 0; i < patterns->size(); ++i) {
    unsigned int mask = 1 << i;
    std::vector<MatchVectType> matchVect;
    // to maintain thread safety, we have to copy the pattern
    // molecules:
    SubstructMatch(mol, ROMol(*(*patterns)[i], true), matchVect);
    for (std::vector<MatchVectType>::const_iterator mvIt = matchVect.begin();
         mvIt != matchVect.end(); ++mvIt) {
      for (const auto &mIt : *mvIt) {
        invars[mIt.second] |= mask;
      }
    }
  }
}  // end of getFeatureInvariants()

void getConnectivityInvariants(const ROMol &mol, std::vector<uint32_t> &invars,
                               bool includeRingMembership) {
  unsigned int nAtoms = mol.getNumAtoms();
  PRECONDITION(invars.size() >= nAtoms, "vector too small");
  gboost::hash<std::vector<uint32_t>> vectHasher;
  for (unsigned int i = 0; i < nAtoms; ++i) {
    Atom const *atom = mol.getAtomWithIdx(i);
    std::vector<uint32_t> components;
    components.push_back(atom->getAtomicNum());
    components.push_back(atom->getTotalDegree());
    components.push_back(atom->getTotalNumHs(true));
    components.push_back(atom->getFormalCharge());
    int deltaMass = static_cast<int>(
        atom->getMass() -
        PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
    components.push_back(deltaMass);

    if (includeRingMembership &&
        atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx())) {
      components.push_back(1);
    }
    invars[i] = vectHasher(components);
  }
}  // end of getConnectivityInvariants()

}  // namespace MorganFingerprints

namespace RDKitFPUtils {

void buildDefaultRDKitFingerprintAtomInvariants(
    const ROMol &mol, std::vector<std::uint32_t> &lAtomInvariants) {
  lAtomInvariants.clear();
  lAtomInvariants.reserve(mol.getNumAtoms());
  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    unsigned int aHash = ((*atomIt)->getAtomicNum() % 128) << 1 |
                         static_cast<unsigned int>((*atomIt)->getIsAromatic());
    lAtomInvariants.push_back(aHash);
  }
}

void enumerateAllPaths(const ROMol &mol, INT_PATH_LIST_MAP &allPaths,
                       const std::vector<std::uint32_t> *fromAtoms,
                       bool branchedPaths, bool useHs, unsigned int minPath,
                       unsigned int maxPath) {
  if (!fromAtoms) {
    if (branchedPaths) {
      allPaths = findAllSubgraphsOfLengthsMtoN(mol, minPath, maxPath, useHs);
    } else {
      allPaths = findAllPathsOfLengthsMtoN(mol, minPath, maxPath, true, useHs);
    }
  } else {
    for (auto aidx : *fromAtoms) {
      INT_PATH_LIST_MAP tPaths;
      if (branchedPaths) {
        tPaths =
            findAllSubgraphsOfLengthsMtoN(mol, minPath, maxPath, useHs, aidx);
      } else {
        tPaths =
            findAllPathsOfLengthsMtoN(mol, minPath, maxPath, true, useHs, aidx);
      }
      for (INT_PATH_LIST_MAP::const_iterator tpit = tPaths.begin();
           tpit != tPaths.end(); ++tpit) {
#ifdef VERBOSE_FINGERPRINTING
        std::cerr << "paths from " << aidx << " size: " << tpit->first
                  << std::endl;
        for (auto path : tpit->second) {
          std::cerr << " path: ";
          std::copy(path.begin(), path.end(),
                    std::ostream_iterator<int>(std::cerr, ", "));
          std::cerr << std::endl;
        }
#endif

        allPaths[tpit->first].insert(allPaths[tpit->first].begin(),
                                     tpit->second.begin(), tpit->second.end());
      }
    }
  }
}

void identifyQueryBonds(const ROMol &mol, std::vector<const Bond *> &bondCache,
                        std::vector<short> &isQueryBond) {
  bondCache.resize(mol.getNumBonds());
  ROMol::EDGE_ITER firstB, lastB;
  boost::tie(firstB, lastB) = mol.getEdges();
  while (firstB != lastB) {
    const Bond *bond = mol[*firstB];
    isQueryBond[bond->getIdx()] = 0x0;
    bondCache[bond->getIdx()] = bond;
    if (isComplexQuery(bond)) {
      isQueryBond[bond->getIdx()] = 0x1;
    }
    if (isComplexQuery(bond->getBeginAtom())) {
      isQueryBond[bond->getIdx()] |= 0x2;
    }
    if (isComplexQuery(bond->getEndAtom())) {
      isQueryBond[bond->getIdx()] |= 0x4;
    }
    ++firstB;
  }
}

std::vector<unsigned int> generateBondHashes(
    const ROMol &mol, boost::dynamic_bitset<> &atomsInPath,
    const std::vector<const Bond *> &bondCache,
    const std::vector<short> &isQueryBond, const PATH_TYPE &path,
    bool useBondOrder, const std::vector<std::uint32_t> *atomInvariants) {
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");

  std::vector<unsigned int> bondHashes;
  atomsInPath.reset();
  bool queryInPath = false;
  std::vector<unsigned int> atomDegrees(mol.getNumAtoms(), 0);
  for (unsigned int i = 0; i < path.size() && !queryInPath; ++i) {
    const Bond *bi = bondCache[path[i]];
    CHECK_INVARIANT(bi, "bond not in cache");
    atomDegrees[bi->getBeginAtomIdx()]++;
    atomDegrees[bi->getEndAtomIdx()]++;
    atomsInPath.set(bi->getBeginAtomIdx());
    atomsInPath.set(bi->getEndAtomIdx());
    if (isQueryBond[path[i]]) {
      queryInPath = true;
    }
  }
  if (queryInPath) {
    return bondHashes;
  }

  // -----------------
  // calculate the bond hashes:
  std::vector<unsigned int> bondNbrs(path.size(), 0);
  bondHashes.reserve(path.size() + 1);

  for (unsigned int i = 0; i < path.size(); ++i) {
    const Bond *bi = bondCache[path[i]];
#ifdef REPORT_FP_STATS
    if (std::find(atomsToUse.begin(), atomsToUse.end(),
                  bi->getBeginAtomIdx()) == atomsToUse.end()) {
      atomsToUse.push_back(bi->getBeginAtomIdx());
    }
    if (std::find(atomsToUse.begin(), atomsToUse.end(), bi->getEndAtomIdx()) ==
        atomsToUse.end()) {
      atomsToUse.push_back(bi->getEndAtomIdx());
    }
#endif
    for (unsigned int j = i + 1; j < path.size(); ++j) {
      const Bond *bj = bondCache[path[j]];
      if (bi->getBeginAtomIdx() == bj->getBeginAtomIdx() ||
          bi->getBeginAtomIdx() == bj->getEndAtomIdx() ||
          bi->getEndAtomIdx() == bj->getBeginAtomIdx() ||
          bi->getEndAtomIdx() == bj->getEndAtomIdx()) {
        ++bondNbrs[i];
        ++bondNbrs[j];
      }
    }
#ifdef VERBOSE_FINGERPRINTING
    std::cerr << "   bond(" << i << "):" << bondNbrs[i] << std::endl;
#endif
    // we have the count of neighbors for bond bi, compute its hash:
    unsigned int a1Hash = (*atomInvariants)[bi->getBeginAtomIdx()];
    unsigned int a2Hash = (*atomInvariants)[bi->getEndAtomIdx()];
    unsigned int deg1 = atomDegrees[bi->getBeginAtomIdx()];
    unsigned int deg2 = atomDegrees[bi->getEndAtomIdx()];
    if (a1Hash < a2Hash) {
      std::swap(a1Hash, a2Hash);
      std::swap(deg1, deg2);
    } else if (a1Hash == a2Hash && deg1 < deg2) {
      std::swap(deg1, deg2);
    }
    unsigned int bondHash = 1;
    if (useBondOrder) {
      if (bi->getIsAromatic() || bi->getBondType() == Bond::AROMATIC) {
        // makes sure aromatic bonds always hash as aromatic
        bondHash = Bond::AROMATIC;
      } else {
        bondHash = bi->getBondType();
      }
    }
    std::uint32_t ourHash = bondNbrs[i];
    gboost::hash_combine(ourHash, bondHash);
    gboost::hash_combine(ourHash, a1Hash);
    gboost::hash_combine(ourHash, deg1);
    gboost::hash_combine(ourHash, a2Hash);
    gboost::hash_combine(ourHash, deg2);
    bondHashes.push_back(ourHash);
    // std::cerr<<"    "<<bi->getIdx()<<"
    // "<<a1Hash<<"("<<deg1<<")"<<"-"<<a2Hash<<"("<<deg2<<")"<<" "<<bondHash<<"
    // -> "<<ourHash<<std::endl;
  }
  return bondHashes;
}

}  // namespace RDKitFPUtils

}  // namespace RDKit
