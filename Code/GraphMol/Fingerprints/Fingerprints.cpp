//
//  Copyright (C) 2003-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include "Fingerprints.h"
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/random.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <limits>
#include <RDGeneral/hash/hash.hpp>
#include <RDGeneral/types.h>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <GraphMol/Fingerprints/FingerprintUtil.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>

//#define VERBOSE_FINGERPRINTING 1
//#define REPORT_FP_STATS 1
#ifdef REPORT_FP_STATS
#include <GraphMol/SmilesParse/SmilesWrite.h>
#endif

namespace RDKit {
namespace Fingerprints {
namespace detail {}  // namespace detail
}  // namespace Fingerprints
namespace {
/*
std::uint32_t hashBond(const Bond *bnd,const std::vector<std::uint32_t>
&atomInvariants,
              const std::vector<std::uint32_t> &atomDegrees,std::uint32_t
bondDegree,
              bool useBondOrder){
PRECONDITION(bnd,"bad bond");
std::uint32_t res;
if(useBondOrder) {
if(bnd->getIsAromatic()){
  res = Bond::AROMATIC;
} else {
  res=bnd->getBondType();
}
} else {
res = 1;
}
std::uint32_t iv1=atomInvariants[bnd->getBeginAtomIdx()];
std::uint32_t iv2=atomInvariants[bnd->getEndAtomIdx()];
std::uint32_t deg1=atomDegrees[bnd->getBeginAtomIdx()];
std::uint32_t deg2=atomDegrees[bnd->getEndAtomIdx()];

if(iv1>iv2){
std::swap(iv1,iv2);
std::swap(deg1,deg2);
} else if(iv1==iv2){
if(deg1>deg2){
  std::swap(deg1,deg2);
}
}

res = (res%8) | (iv1%128)<<3 | (iv2%128)<<10 | (deg1%8)<<17 | (deg2%8)<<20 |
(bondDegree%8)<<23 ;
//std::cerr<<"---->("<<bnd->getIdx()<<")
"<<bnd->getBeginAtomIdx()<<"-"<<bnd->getEndAtomIdx()<<" "<<res<<"
"<<iv1<<"-"<<iv2<<":"<<deg1<<"-"<<deg2<<std::endl;
return res;
}

std::uint32_t canonicalPathHash(const PATH_TYPE &path,
                       const ROMol &mol,
                       const std::vector<const Bond *> &bondCache,
                       const std::vector<std::uint32_t> &bondHashes){
std::deque< std::pair<unsigned int,boost::dynamic_bitset<> > > stack;
std::uint32_t best;
//std::cerr<<" hash: ";
//std::copy(path.begin(),path.end(),std::ostream_iterator<int>(std::cerr,", "));

for(unsigned int i=0;i<path.size();++i){
//std::cerr<<"
"<<bondCache[path[i]]->getBeginAtomIdx()<<"-"<<bondCache[path[i]]->getEndAtomIdx()<<"
"<<bondHashes[i];
if(i==0){
  boost::dynamic_bitset<> bs(mol.getNumBonds());
  bs.set(path[i]);
  stack.push_back(std::make_pair(i,bs));
  best=bondHashes[i];
} else {
  if(bondHashes[i]<=best){
    if(bondHashes[i]<best){
      stack.clear();
      best = bondHashes[i];
    }
    boost::dynamic_bitset<> bs(mol.getNumBonds());
    bs.set(path[i]);
    stack.push_back(std::make_pair(i,bs));
  }
}
}
//std::cerr<<std::endl;

std::uint32_t res=best;
//std::cerr<<"  best: "<<best<<std::endl;
if(path.size()==1) return res;
best = std::numeric_limits<std::uint32_t>::max();
std::deque< std::pair<unsigned int,boost::dynamic_bitset<> > > newStack;
while(!stack.empty()){
// assumption: each element of the stack corresponds to
// the last point of a traversal of the path
// res has been updated with all elements already traversed

unsigned int i;
boost::dynamic_bitset<> bondsThere;
boost::tie(i,bondsThere)=stack.front();

//std::cerr<<" "<<path[i]<<"("<<bondsThere<<")";

const Bond *bnd=bondCache[path[i]];
for(unsigned int j=0;j<path.size();++j){
  //std::cerr<<" c:"<<path[j];
  if(bondsThere[path[j]]) {
    //std::cerr<<"x";
    continue;
  }
  const Bond *obnd=bondCache[path[j]];
  if(bondHashes[j]>best) continue;
  if(obnd->getBeginAtomIdx()==bnd->getBeginAtomIdx() ||
     obnd->getBeginAtomIdx()==bnd->getEndAtomIdx() ||
     obnd->getEndAtomIdx()==bnd->getBeginAtomIdx() ||
     obnd->getEndAtomIdx()==bnd->getEndAtomIdx() ){
    // it's a neighbor and the hash is at least as good as what we've seen so
far
    if(bondHashes[j]<best){
      newStack.clear();
      best=bondHashes[j];
    }
    boost::dynamic_bitset<> bs(bondsThere);
    bs.set(path[j]);
    newStack.push_back(std::make_pair(j,bs));
    //std::cerr<<"  "<<path[j];
  }
}

stack.pop_front();
if(stack.empty()){
  //std::cerr<<"\n     new round "<<" best: "<<best<<" res: "<<res<<" sz:
"<<newStack.size();
  // at the end of this round, start the next one
  gboost::hash_combine(res,best);
  //std::cerr<<" nres: "<<res<<std::endl;
  //stack=newStack;
  std::swap(stack,newStack);
  best = std::numeric_limits<std::uint32_t>::max();
  newStack.clear();
}
}
gboost::hash_combine(res,path.size());
return res;
}
*/

}  // end of anonymous namespace

// caller owns the result, it must be deleted
ExplicitBitVect *RDKFingerprintMol(
    const ROMol &mol, unsigned int minPath, unsigned int maxPath,
    unsigned int fpSize, unsigned int nBitsPerHash, bool useHs,
    double tgtDensity, unsigned int minSize, bool branchedPaths,
    bool useBondOrder, std::vector<std::uint32_t> *atomInvariants,
    const std::vector<std::uint32_t> *fromAtoms,
    std::vector<std::vector<std::uint32_t>> *atomBits,
    std::map<std::uint32_t, std::vector<std::vector<int>>> *bitInfo) {
  PRECONDITION(minPath != 0, "minPath==0");
  PRECONDITION(maxPath >= minPath, "maxPath<minPath");
  PRECONDITION(fpSize != 0, "fpSize==0");
  PRECONDITION(nBitsPerHash != 0, "nBitsPerHash==0");
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");
  PRECONDITION(!atomBits || atomBits->size() >= mol.getNumAtoms(),
               "bad atomBits size");

  std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen(
      RDKit::RDKitFP::getRDKitFPGenerator<std::uint32_t>(
          minPath, maxPath, useHs, branchedPaths, useBondOrder));
  fpgen->getOptions()->d_fpSize = fpSize;
  fpgen->getOptions()->d_numBitsPerFeature = nBitsPerHash;

  FingerprintFuncArguments args;
  args.customAtomInvariants = atomInvariants;
  args.fromAtoms = fromAtoms;

  AdditionalOutput ao;
  if (atomBits) {
    args.additionalOutput = &ao;
    ao.allocateAtomToBits();
  }
  if (bitInfo) {
    args.additionalOutput = &ao;
    ao.allocateBitPaths();
  }

  auto res = fpgen->getFingerprint(mol, args).release();

  if (atomBits) {
    atomBits->clear();
    for (const auto &abl : *ao.atomToBits) {
      std::vector<std::uint32_t> uv;
      uv.reserve(abl.size());
      for (auto l : abl) {
        uv.push_back(static_cast<std::uint32_t>(l));
      }
      atomBits->emplace_back(std::move(uv));
    }
  }

  if (bitInfo) {
    bitInfo->clear();
    for (const auto &abl : *ao.bitPaths) {
      (*bitInfo)[static_cast<std::uint32_t>(abl.first)] = abl.second;
    }
  }

  // EFF: this could be faster by folding by more than a factor
  // of 2 each time, but we're not going to be spending much
  // time here anyway
  if (tgtDensity > 0.0) {
    while (static_cast<double>(res->getNumOnBits()) / res->getNumBits() <
               tgtDensity &&
           res->getNumBits() >= 2 * minSize) {
      ExplicitBitVect *tmpV = FoldFingerprint(*res, 2);
      delete res;
      res = tmpV;
    }
  }
  return res;
}

// caller owns the result, it must be deleted
ExplicitBitVect *LayeredFingerprintMol(
    const ROMol &mol, unsigned int layerFlags, unsigned int minPath,
    unsigned int maxPath, unsigned int fpSize,
    std::vector<unsigned int> *atomCounts, ExplicitBitVect *setOnlyBits,
    bool branchedPaths, const std::vector<std::uint32_t> *fromAtoms) {
  PRECONDITION(minPath != 0, "minPath==0");
  PRECONDITION(maxPath >= minPath, "maxPath<minPath");
  PRECONDITION(fpSize != 0, "fpSize==0");
  PRECONDITION(!atomCounts || atomCounts->size() >= mol.getNumAtoms(),
               "bad atomCounts size");
  PRECONDITION(!setOnlyBits || setOnlyBits->getNumBits() == fpSize,
               "bad setOnlyBits size");

  if (!mol.getRingInfo()->isInitialized()) {
    MolOps::findSSSR(mol);
  }

  std::vector<const Bond *> bondCache;
  bondCache.resize(mol.getNumBonds());
  std::vector<short> isQueryBond(mol.getNumBonds(), 0);
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

  std::vector<bool> aromaticAtoms(mol.getNumAtoms(), false);
  std::vector<int> anums(mol.getNumAtoms(), 0);
  ROMol::VERTEX_ITER firstA, lastA;
  boost::tie(firstA, lastA) = mol.getVertices();
  while (firstA != lastA) {
    const Atom *atom = mol[*firstA];
    if (isAtomAromatic(atom)) {
      aromaticAtoms[atom->getIdx()] = true;
    }
    anums[atom->getIdx()] = atom->getAtomicNum();
    ++firstA;
  }

  auto *res = new ExplicitBitVect(fpSize);

  INT_PATH_LIST_MAP allPaths;
  if (!fromAtoms) {
    if (branchedPaths) {
      allPaths = findAllSubgraphsOfLengthsMtoN(mol, minPath, maxPath, false);
    } else {
      allPaths = findAllPathsOfLengthsMtoN(mol, minPath, maxPath, false);
    }
  } else {
    for (auto aidx : *fromAtoms) {
      INT_PATH_LIST_MAP tPaths;
      if (branchedPaths) {
        tPaths =
            findAllSubgraphsOfLengthsMtoN(mol, minPath, maxPath, false, aidx);
      } else {
        tPaths =
            findAllPathsOfLengthsMtoN(mol, minPath, maxPath, true, false, aidx);
      }
      for (INT_PATH_LIST_MAP::const_iterator tpit = tPaths.begin();
           tpit != tPaths.end(); ++tpit) {
        allPaths[tpit->first].insert(allPaths[tpit->first].begin(),
                                     tpit->second.begin(), tpit->second.end());
      }
    }
  }

  boost::dynamic_bitset<> atomsInPath(mol.getNumAtoms());
  boost::dynamic_bitset<> bondsInPath(mol.getNumBonds());
  for (INT_PATH_LIST_MAP_CI paths = allPaths.begin(); paths != allPaths.end();
       ++paths) {
    for (const auto &path : paths->second) {
#ifdef VERBOSE_FINGERPRINTING
      std::cerr << "Path: ";
      std::copy(path.begin(), path.end(),
                std::ostream_iterator<int>(std::cerr, ", "));
      std::cerr << std::endl;
#endif

      std::vector<std::vector<unsigned int>> hashLayers(maxFingerprintLayers);
      for (unsigned int i = 0; i < maxFingerprintLayers; ++i) {
        if (layerFlags & (0x1 << i)) {
          hashLayers[i].reserve(maxPath);
        }
      }

      // details about what kinds of query features appear on the path:
      unsigned int pathQueries = 0;
      // std::cerr<<" path: ";
      for (int pIt : path) {
        pathQueries |= isQueryBond[pIt];
        // std::cerr<< *pIt <<"("<<isQueryBond[*pIt]<<") ";
      }
      // std::cerr<<" : "<<pathQueries<<std::endl;

      // calculate the number of neighbors each bond has in the path:
      std::vector<unsigned int> bondNbrs(path.size(), 0);
      atomsInPath.reset();

      std::vector<unsigned int> atomDegrees(mol.getNumAtoms(), 0);
      for (int i : path) {
        const Bond *bi = bondCache[i];
        atomDegrees[bi->getBeginAtomIdx()]++;
        atomDegrees[bi->getEndAtomIdx()]++;
        atomsInPath.set(bi->getBeginAtomIdx());
        atomsInPath.set(bi->getEndAtomIdx());
      }

      for (unsigned int i = 0; i < path.size(); ++i) {
        const Bond *bi = bondCache[path[i]];
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
        // we have the count of neighbors for bond bi, compute its hash layers:
        unsigned int ourHash = 0;

        if (layerFlags & 0x1) {
          // layer 1: straight topology
          unsigned int a1Deg, a2Deg;
          a1Deg = atomDegrees[bi->getBeginAtomIdx()];
          a2Deg = atomDegrees[bi->getEndAtomIdx()];
          if (a1Deg < a2Deg) {
            std::swap(a1Deg, a2Deg);
          }
          ourHash = bondNbrs[i] % 8;  // 3 bits here
          ourHash |= (a1Deg % 8) << 3;
          ourHash |= (a2Deg % 8) << 6;
          hashLayers[0].push_back(ourHash);
        }
        if (layerFlags & 0x2 && !(pathQueries & 0x1)) {
          // layer 2: include bond orders:
          unsigned int bondHash;
          // makes sure aromatic bonds and single bonds  always hash the same:
          if (!bi->getIsAromatic() && bi->getBondType() != Bond::SINGLE &&
              bi->getBondType() != Bond::AROMATIC) {
            bondHash = bi->getBondType();
          } else {
            bondHash = Bond::SINGLE;
          }
          unsigned int a1Deg, a2Deg;
          a1Deg = atomDegrees[bi->getBeginAtomIdx()];
          a2Deg = atomDegrees[bi->getEndAtomIdx()];
          if (a1Deg < a2Deg) {
            std::swap(a1Deg, a2Deg);
          }
          ourHash = bondHash % 8;
          ourHash |= (bondNbrs[i] % 8) << 3;
          ourHash |= (a1Deg % 8) << 6;
          ourHash |= (a2Deg % 8) << 9;

          hashLayers[1].push_back(ourHash);
        }
        if (layerFlags & 0x4 && !(pathQueries & 0x6)) {
          // std::cerr<<" consider: "<<bi->getBeginAtomIdx()<<" - "
          // <<bi->getEndAtomIdx()<<std::endl;
          // layer 3: include atom types:
          unsigned int a1Hash, a2Hash;
          a1Hash = (anums[bi->getBeginAtomIdx()] % 128);
          a2Hash = (anums[bi->getEndAtomIdx()] % 128);
          unsigned int a1Deg, a2Deg;
          a1Deg = atomDegrees[bi->getBeginAtomIdx()];
          a2Deg = atomDegrees[bi->getEndAtomIdx()];
          if (a1Hash < a2Hash) {
            std::swap(a1Hash, a2Hash);
            std::swap(a1Deg, a2Deg);
          } else if (a1Hash == a2Hash && a1Deg < a2Deg) {
            std::swap(a1Deg, a2Deg);
          }
          ourHash = a1Hash;
          ourHash |= a2Hash << 7;
          ourHash |= (a1Deg % 8) << 14;
          ourHash |= (a2Deg % 8) << 17;
          ourHash |= (bondNbrs[i] % 8) << 20;
          hashLayers[2].push_back(ourHash);
        }
        if (layerFlags & 0x8 && !(pathQueries & 0x6)) {
          // layer 4: include ring information
          if (queryIsBondInRing(bi)) {
            hashLayers[3].push_back(1);
          }
        }
        if (layerFlags & 0x10 && !(pathQueries & 0x6)) {
          // layer 5: include ring size information
          ourHash = (queryBondMinRingSize(bi) % 8);
          hashLayers[4].push_back(ourHash);
        }
        if (layerFlags & 0x20 && !(pathQueries & 0x6)) {
          // std::cerr<<" consider: "<<bi->getBeginAtomIdx()<<" - "
          // <<bi->getEndAtomIdx()<<std::endl;
          // layer 6: aromaticity:
          bool a1Hash = aromaticAtoms[bi->getBeginAtomIdx()];
          bool a2Hash = aromaticAtoms[bi->getEndAtomIdx()];

          if ((!a1Hash) && a2Hash) {
            std::swap(a1Hash, a2Hash);
          }
          ourHash = a1Hash;
          ourHash |= a2Hash << 1;
          ourHash |= (bondNbrs[i] % 8) << 5;
          hashLayers[5].push_back(ourHash);
        }
      }
      unsigned int l = 0;
      bool flaggedPath = false;
      for (auto layerIt = hashLayers.begin(); layerIt != hashLayers.end();
           ++layerIt, ++l) {
        if (!layerIt->size()) {
          continue;
        }
        // ----
        std::sort(layerIt->begin(), layerIt->end());

        // finally, we will add the number of distinct atoms in the path at the
        // end
        // of the vect. This allows us to distinguish C1CC1 from CC(C)C
        layerIt->push_back(static_cast<unsigned int>(atomsInPath.count()));

        layerIt->push_back(l + 1);

        // hash the path to generate a seed:
        unsigned long seed =
            gboost::hash_range(layerIt->begin(), layerIt->end());

#ifdef VERBOSE_FINGERPRINTING
        std::cerr << " hash: " << seed << std::endl;
#endif
        unsigned int bitId = seed % fpSize;
#ifdef VERBOSE_FINGERPRINTING
        std::cerr << "   bit: " << bitId << std::endl;
#endif
        if (!setOnlyBits || (*setOnlyBits)[bitId]) {
          res->setBit(bitId);
          if (atomCounts && !flaggedPath) {
            for (unsigned int aIdx = 0; aIdx < atomsInPath.size(); ++aIdx) {
              if (atomsInPath[aIdx]) {
                (*atomCounts)[aIdx] += 1;
              }
            }
            flaggedPath = true;
          }
        }
      }
    }
  }
  return res;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// caller owns the result, it must be deleted
SparseIntVect<boost::uint64_t> *getUnfoldedRDKFingerprintMol(
    const ROMol &mol, unsigned int minPath, unsigned int maxPath, bool useHs,
    bool branchedPaths, bool useBondOrder,
    std::vector<std::uint32_t> *atomInvariants,
    const std::vector<std::uint32_t> *fromAtoms,
    std::vector<std::vector<boost::uint64_t>> *atomBits,
    std::map<boost::uint64_t, std::vector<std::vector<int>>> *bitInfo) {
  PRECONDITION(minPath != 0, "minPath==0");
  PRECONDITION(maxPath >= minPath, "maxPath<minPath");
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");
  PRECONDITION(!atomBits || atomBits->size() >= mol.getNumAtoms(),
               "bad atomBits size");

  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpgen(
      RDKit::RDKitFP::getRDKitFPGenerator<std::uint64_t>(
          minPath, maxPath, useHs, branchedPaths, useBondOrder));
  fpgen->getOptions()->d_numBitsPerFeature = 1;
  FingerprintFuncArguments args;
  args.customAtomInvariants = atomInvariants;
  args.fromAtoms = fromAtoms;

  AdditionalOutput ao;
  if (atomBits) {
    args.additionalOutput = &ao;
    ao.allocateAtomToBits();
  }
  if (bitInfo) {
    args.additionalOutput = &ao;
    ao.allocateBitPaths();
  }

  auto fp = fpgen->getSparseCountFingerprint(mol, args);

  if (atomBits) {
    *atomBits = *ao.atomToBits;
  }

  if (bitInfo) {
    *bitInfo = *ao.bitPaths;
  }

  return fp.release();
}
}  // namespace RDKit
