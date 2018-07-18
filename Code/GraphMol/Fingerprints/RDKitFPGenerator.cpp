
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <RDGeneral/hash/hash.hpp>
#include <cstdint>

#include <GraphMol/QueryOps.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/random.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <limits.h>
#include <RDGeneral/hash/hash.hpp>
#include <RDGeneral/types.h>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>

namespace RDKit {
namespace RDKitFP {

namespace utils {

void enumerateAllPaths(const ROMol &mol, INT_PATH_LIST_MAP &allPaths,
                       const std::vector<boost::uint32_t> *fromAtoms,
                       bool branchedPaths, bool useHs, unsigned int minPath,
                       unsigned int maxPath) {
  if (!fromAtoms) {
    if (branchedPaths) {
      allPaths = findAllSubgraphsOfLengthsMtoN(mol, minPath, maxPath, useHs);
    } else {
      allPaths = findAllPathsOfLengthsMtoN(mol, minPath, maxPath, true, useHs);
    }
  } else {
    BOOST_FOREACH (boost::uint32_t aidx, *fromAtoms) {
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
        BOOST_FOREACH (PATH_TYPE path, tpit->second) {
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
    bool useBondOrder, const std::vector<boost::uint32_t> *atomInvariants) {
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
    if (isQueryBond[path[i]]) queryInPath = true;
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
    boost::uint32_t ourHash = bondNbrs[i];
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

}  // namespace utils

std::vector<std::uint32_t> *RDKitFPAtomInvGenerator::getAtomInvariants(
    const ROMol &mol) const {
  std::vector<std::uint32_t> *result = new std::vector<std::uint32_t>();
  result->reserve(mol.getNumAtoms());
  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    unsigned int aHash = ((*atomIt)->getAtomicNum() % 128) << 1 |
                         static_cast<unsigned int>((*atomIt)->getIsAromatic());
    result->push_back(aHash);
  }
  return result;
}

std::string RDKitFPAtomInvGenerator::infoString() const {
  return "RDKitFPAtomInvGenerator";
}

std::uint64_t RDKitFPArguments::getResultSize() const {
  return std::numeric_limits<uint32_t>::max();
}

std::string RDKitFPArguments::infoString() const {
  return "RDKitFPArguments minPath=" + std::to_string(d_minPath) +
         " maxPath=" + std::to_string(d_maxPath) +
         " useHs=" + std::to_string(df_useHs) +
         " branchedPaths=" + std::to_string(df_branchedPaths) +
         " useBondOrder=" + std::to_string(df_useBondOrder);
}

RDKitFPArguments::RDKitFPArguments(unsigned int minPath, unsigned int maxPath,
                                   bool useHs, bool branchedPaths,
                                   bool useBondOrder,
                                   const bool countSimulation,
                                   const std::vector<std::uint32_t> countBounds,
                                   const std::uint32_t foldedSize)
    : FingerprintArguments(countSimulation, countBounds, foldedSize),
      d_minPath(minPath),
      d_maxPath(maxPath),
      df_useHs(useHs),
      df_branchedPaths(branchedPaths),
      df_useBondOrder(useBondOrder) {
  PRECONDITION(minPath != 0, "minPath==0");
  PRECONDITION(maxPath >= minPath, "maxPath<minPath");
}

std::uint32_t RDKitFPAtomEnv::getBitId(
    FingerprintArguments *arguments,
    const std::vector<std::uint32_t> *atomInvariants,
    const std::vector<std::uint32_t> *bondInvariants,
    const AdditionalOutput *additionalOutput) const {
  // todo set additional outputs
  return d_bitId;
}

RDKitFPAtomEnv::RDKitFPAtomEnv(const std::uint32_t bitId,
                               const boost::dynamic_bitset<> atomsInPath)
    : d_bitId(bitId), d_atomsInPath(atomsInPath) {}

std::string RDKitFPEnvGenerator::infoString() const {
  return "RDKitFPEnvGenerator";
}

std::vector<AtomEnvironment *> RDKitFPEnvGenerator::getEnvironments(
    const ROMol &mol, FingerprintArguments *arguments,
    const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *atomInvariants,
    const std::vector<std::uint32_t> *bondInvariants) const {
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");

  RDKitFPArguments *rDKitFPArguments =
      dynamic_cast<RDKitFPArguments *>(arguments);

  std::vector<AtomEnvironment *> result;

  // get all paths
  INT_PATH_LIST_MAP allPaths;
  utils::enumerateAllPaths(
      mol, allPaths, fromAtoms, rDKitFPArguments->df_branchedPaths,
      rDKitFPArguments->df_useHs, rDKitFPArguments->d_minPath,
      rDKitFPArguments->d_maxPath);

  // identify query bonds
  std::vector<short> isQueryBond(mol.getNumBonds(), 0);
  std::vector<const Bond *> bondCache;
  utils::identifyQueryBonds(mol, bondCache, isQueryBond);

  boost::dynamic_bitset<> atomsInPath(mol.getNumAtoms());
  for (INT_PATH_LIST_MAP_CI paths = allPaths.begin(); paths != allPaths.end();
       paths++) {
    BOOST_FOREACH (const PATH_TYPE &path, paths->second) {
      // the bond hashes of the path
      std::vector<unsigned int> bondHashes = utils::generateBondHashes(
          mol, atomsInPath, bondCache, isQueryBond, path,
          rDKitFPArguments->df_useBondOrder, atomInvariants);
      if (!bondHashes.size()) {
        continue;
      }

      // hash the path to generate a seed:
      unsigned long seed;
      if (path.size() > 1) {
        std::sort(bondHashes.begin(), bondHashes.end());

        // finally, we will add the number of distinct atoms in the path at the
        // end
        // of the vect. This allows us to distinguish C1CC1 from CC(C)C
        bondHashes.push_back(static_cast<unsigned int>(atomsInPath.count()));
        seed = gboost::hash_range(bondHashes.begin(), bondHashes.end());
      } else {
        seed = bondHashes[0];
      }

      result.push_back(new RDKitFPAtomEnv(seed, atomsInPath));
    }
  }

  return result;
}

FingerprintGenerator *getRDKitFPGenerator(
    const unsigned int minPath, const unsigned int maxPath, const bool useHs,
    const bool branchedPaths, const bool useBondOrder,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator,
    const bool countSimulation, const std::vector<std::uint32_t> countBounds,
    const std::uint32_t foldedSize) {
  AtomEnvironmentGenerator *envGenerator = new RDKitFPEnvGenerator();
  FingerprintArguments *arguments =
      new RDKitFPArguments(minPath, maxPath, useHs, branchedPaths, useBondOrder,
                           countSimulation, countBounds, foldedSize);

  bool ownsAtomInvGenerator = false;
  if (!atomInvariantsGenerator) {
    atomInvariantsGenerator = new RDKitFPAtomInvGenerator();
    ownsAtomInvGenerator = true;
  }

  bool ownsBondInvGenerator = false;

  return new FingerprintGenerator(
      envGenerator, arguments, atomInvariantsGenerator, bondInvariantsGenerator,
      ownsAtomInvGenerator, ownsBondInvGenerator);
}

}  // namespace RDKitFP
}  // namespace RDKit
