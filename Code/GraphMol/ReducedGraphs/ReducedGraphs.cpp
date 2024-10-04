//
//  Copyright (C) 2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <Numerics/Vector.h>

#include "ReducedGraphs.h"

#include <RDGeneral/BoostStartInclude.h>
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <boost/ref.hpp>
#include <RDGeneral/BoostEndInclude.h>

// #define VERBOSE_FINGERPRINTING 1

namespace RDKit {

namespace {
// FIX: this is duplicated here and in the MorganFingerprints code
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

namespace ReducedGraphs {
// Definitions for feature points adapted from:
// Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998)
const int nFeatures = 5;
const char *smartsPatterns[nFeatures] = {
    "[$([N;!H0;v3,v4&+1]),\
$([O,S;H1;+0]),\
n&H1&+0]",                                                  // Donor
    "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),\
$([O;H0;v2]),\
$([O,S;v1;-]),\
$([N;v3;!$(N-*=[O,N,P,S])]),\
n&H0&+0,\
$([o;+0;!$([o]:n);!$([o]:c:n)])]",                          // Acceptor
    "[#7;+,\
$([N;H2&+0][$([C,a]);!$([C,a](=O))]),\
$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),\
$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]",  // Positive
    "[$([C,S](=[O,S,P])-[O;H1,-1])]",                       // Negative
    "[$([C;D3,D4](-[CH3])-[CH3]),$([S;D2](-C)-C)]"  // Hydrophobic, note that
                                                    // this needs to be last.
                                                    // // FIX: not at all sure
                                                    // that this is right
};
std::vector<std::string> defaultFeatureSmarts(smartsPatterns,
                                              smartsPatterns + nFeatures);
typedef boost::flyweight<boost::flyweights::key_value<std::string, ss_matcher>,
                         boost::flyweights::no_tracking>
    pattern_flyweight;

void getErGAtomTypes(const ROMol &mol,
                     std::vector<boost::dynamic_bitset<>> &types,
                     std::vector<const ROMol *> *patterns = nullptr) {
  unsigned int nAtoms = mol.getNumAtoms();

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

  types.resize(patterns->size());

  for (unsigned int i = 0; i < patterns->size(); ++i) {
    types[i].resize(nAtoms);
    types[i].reset();
    // unsigned int mask=1<<i;
    std::vector<MatchVectType> matchVect;
    // to maintain thread safety, we have to copy the pattern
    // molecules:
    SubstructMatch(mol, ROMol(*(*patterns)[i], true), matchVect);
    for (std::vector<MatchVectType>::const_iterator mvIt = matchVect.begin();
         mvIt != matchVect.end(); ++mvIt) {
      types[i].set((*mvIt)[0].second);
    }
  }
}  // end of getAtomTypes;

RDNumeric::DoubleVector *getErGFingerprint(
    const ROMol &mol, std::vector<boost::dynamic_bitset<>> *atomTypes,
    double fuzzIncrement, unsigned int minPath, unsigned int maxPath) {
  ROMol *rg = generateMolExtendedReducedGraph(mol, atomTypes);
#ifdef VERBOSE_FINGERPRINTING
  rg->updatePropertyCache(false);
  std::cerr << " reduced graph smiles: "
            << MolToSmiles(*rg, false, false, -1, false) << std::endl;
#endif
  RDNumeric::DoubleVector *res = generateErGFingerprintForReducedGraph(
      *rg, atomTypes, fuzzIncrement, minPath, maxPath);
  delete rg;
  return res;
}

RDNumeric::DoubleVector *generateErGFingerprintForReducedGraph(
    const ROMol &mol, std::vector<boost::dynamic_bitset<>> *atomTypes,
    double fuzzIncrement, unsigned int minPath, unsigned int maxPath) {
  PRECONDITION(maxPath > minPath, "maxPath<=minPath");
  // FIX: this isn't doing the special handling for flip/flop bits
  unsigned int nTypes = nFeatures;
  if (atomTypes) {
    nTypes = atomTypes->size();
  }
  nTypes += 1;  // need the extra one for aromatic features
  unsigned int nBins = maxPath - minPath;

  unsigned int vSize = (nTypes * (nTypes + 1)) / 2 * (maxPath - minPath + 1);
  auto *res = new RDNumeric::DoubleVector(vSize, 0.0);

  // we need the topological distance matrix:
  double *dm = MolOps::getDistanceMat(mol);

  // cache the atom type vectors:
  std::vector<std::vector<int>> tvs;
  tvs.reserve(mol.getNumAtoms());
  for (ROMol::ConstAtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
       ++atIt) {
    const std::vector<int> &tv =
        (*atIt)->getProp<std::vector<int>>("_ErGAtomTypes");
    tvs.push_back(tv);
  }

  // now loop and set the bits between features
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    if (tvs[i].empty()) {
      continue;
    }
#ifdef VERBOSE_FINGERPRINTING
    std::cerr << " atom: " << i << " ";
    std::copy(tvs[i].begin(), tvs[i].end(),
              std::ostream_iterator<int>(std::cerr, ", "));
    std::cerr << std::endl;
#endif
    for (unsigned int j = i + 1; j < mol.getNumAtoms(); ++j) {
      if (tvs[j].empty()) {
        continue;
      }
      int dist = int(dm[i * mol.getNumAtoms() + j]);
      if (dist < rdcast<int>(minPath) || dist > rdcast<int>(maxPath)) {
        continue;
      }
      for (auto ti : tvs[i]) {
        for (auto tj : tvs[j]) {
          int ijMin = std::min(ti, tj);
          int ijMax = std::max(ti, tj);
          int block;
          if (ijMin > 0) {
            block = ijMin * nTypes - (ijMin * (ijMin + 1)) / 2;
          } else {
            block = 0;
          }
          block += ijMax;
          unsigned int bin = block * nBins + (dist - minPath);
          (*res)[bin] += 1;
          if (dist > rdcast<int>(minPath) && fuzzIncrement > 0) {
            (*res)[bin - 1] += fuzzIncrement;
          }
          if (dist < rdcast<int>(maxPath) && fuzzIncrement > 0) {
            (*res)[bin + 1] += fuzzIncrement;
          }
        }
      }
    }
  }
  return res;
}

ROMol *generateMolExtendedReducedGraph(
    const ROMol &mol, std::vector<boost::dynamic_bitset<>> *atomTypes) {
  std::vector<boost::dynamic_bitset<>> *latomTypes = nullptr;
  if (!atomTypes) {
    latomTypes = new std::vector<boost::dynamic_bitset<>>();
    atomTypes = latomTypes;
    getErGAtomTypes(mol, *atomTypes);
  }
  auto *res = new RWMol(mol);

  const int aliphaticFlag = atomTypes->size() - 1;  // the last type
  const int aromaticFlag = atomTypes->size();

  for (ROMol::AtomIterator atIt = res->beginAtoms(); atIt != res->endAtoms();
       ++atIt) {
    std::vector<int> tv;
    tv.clear();
    for (unsigned int i = 0; i < atomTypes->size(); ++i) {
      if ((*atomTypes)[i][(*atIt)->getIdx()]) {
        tv.push_back(i);
      }
    }
    (*atIt)->setProp("_ErGAtomTypes", tv);
  }

  // start by adding dummies at the ring centroids
  for (const auto &ring : mol.getRingInfo()->atomRings()) {
    if (ring.size() < 8) {
      int nIdx = res->addAtom(new Atom(0), false, true);
      int nAromatic = 0, nSP2 = 0;
      for (auto idx : ring) {
        res->addBond(idx, nIdx, Bond::SINGLE);
        if (mol.getAtomWithIdx(idx)->getIsAromatic()) {
          ++nAromatic;
        } else if (mol.getAtomWithIdx(idx)->getHybridization() == Atom::SP2) {
          ++nSP2;
        }
      }
      std::vector<int> tv;
      if (nAromatic >= 2 || nSP2 >= rdcast<int>(ring.size() / 2)) {
        tv.push_back(aromaticFlag);
      } else {
        tv.push_back(aliphaticFlag);
      }
      res->getAtomWithIdx(nIdx)->setProp("_ErGAtomTypes", tv);
    }
  }

  // now remove any degree-two ring atoms that have no features:
  res->beginBatchEdit();
  for (int i = mol.getNumAtoms() - 1; i >= 0; --i) {
    if (mol.getRingInfo()->numAtomRings(i) &&
        mol.getAtomWithIdx(i)->getDegree() == 2 &&
        res->getAtomWithIdx(i)
            ->getProp<std::vector<int>>("_ErGAtomTypes")
            .empty()) {
      res->removeAtom(i);
    }
  }
  res->commitBatchEdit();

  // FIX: still need to do the "highly fused rings" simplification for things
  // like adamantane

  if (latomTypes) {
    delete latomTypes;
  }
  return res;
}
}  // end of namespace ReducedGraphs
}  // end of namespace RDKit
