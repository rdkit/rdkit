#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>

namespace RDKit {
namespace AtomPair {

unsigned int numPiElectrons(const Atom *atom) {
  PRECONDITION(atom, "no atom");
  unsigned int res = 0;
  if (atom->getIsAromatic()) {
    res = 1;
  } else if (atom->getHybridization() != Atom::SP3) {
    unsigned int val = static_cast<unsigned int>(atom->getExplicitValence());
    val -= atom->getNumExplicitHs();
    CHECK_INVARIANT(val >= atom->getDegree(),
                    "explicit valence exceeds atom degree");
    res = val - atom->getDegree();
  }
  return res;
}

boost::uint32_t getAtomCode(const Atom *atom, unsigned int branchSubtract,
                            bool includeChirality) {
  PRECONDITION(atom, "no atom");
  boost::uint32_t code;

  unsigned int numBranches = 0;
  if (atom->getDegree() > branchSubtract) {
    numBranches = atom->getDegree() - branchSubtract;
  }

  code = numBranches % maxNumBranches;
  unsigned int nPi = numPiElectrons(atom) % maxNumPi;
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
  if (typeIdx == nTypes) --typeIdx;
  code |= typeIdx << (numBranchBits + numPiBits);
  if (includeChirality) {
    std::string cipCode;
    if (atom->getPropIfPresent(common_properties::_CIPCode, cipCode)) {
      boost::uint32_t offset = numBranchBits + numPiBits + numTypeBits;
      if (cipCode == "R") {
        code |= 1 << offset;
      } else if (cipCode == "S") {
        code |= 2 << offset;
      }
    }
  }
  POSTCONDITION(code < static_cast<boost::uint32_t>(
                           1 << (codeSize + (includeChirality ? 2 : 0))),
                "code exceeds number of bits");
  return code;
};

boost::uint32_t getAtomPairCode(boost::uint32_t codeI, boost::uint32_t codeJ,
                                unsigned int dist, bool includeChirality) {
  PRECONDITION(dist < maxPathLen, "dist too long");
  boost::uint32_t res = dist;
  res |= std::min(codeI, codeJ) << numPathBits;
  res |= std::max(codeI, codeJ)
         << (numPathBits + codeSize + (includeChirality ? numChiralBits : 0));
  return res;
}

// AtomPairArguments

unsigned int AtomPairArguments::getResultSize() const {
  return (1 << (numAtomPairFingerprintBits + 2 * (includeChirality ? 2 : 0)));
}

AtomPairArguments::AtomPairArguments(const bool countSimulation,
                                     const bool includeChirality,
                                     const bool use2D,
                                     const unsigned int minDistance,
                                     const unsigned int maxDistance)
    : FingerprintArguments(countSimulation),
      includeChirality(includeChirality),
      use2D(use2D),
      minDistance(minDistance),
      maxDistance(maxDistance) {
  PRECONDITION(minDistance <= maxDistance, "bad distances provided");
};

// AtomPairAtomEnv

boost::uint32_t AtomPairAtomEnv::getBitId(
    AtomPairArguments arguments,
    const std::vector<boost::uint32_t> *atomInvariants,
    const std::vector<boost::uint32_t> *bondInvariants) const {
  return getAtomPairCode(atomCodeCache->at(atomIdFirst),
                         atomCodeCache->at(atomIdSecond), distance,
                         arguments.includeChirality);
}

const std::vector<boost::uint32_t> *AtomPairAtomEnv::getAtomCodeCache() const {
  return atomCodeCache;
}

AtomPairAtomEnv::AtomPairAtomEnv(
    const std::vector<boost::uint32_t> *atomCodeCache,
    const unsigned int atomIdFirst, const unsigned int atomIdSecond,
    const unsigned int distance)
    : atomCodeCache(atomCodeCache),
      atomIdFirst(atomIdFirst),
      atomIdSecond(atomIdSecond),
      distance(distance) {}

// AtomPairEnvGenerator

std::vector<AtomPairAtomEnv> AtomPairEnvGenerator::getEnvironments(
    const ROMol &mol, const AtomPairArguments arguments,
    const std::vector<boost::uint32_t> *fromAtoms,
    const std::vector<boost::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<boost::uint32_t> *atomInvariants,
    const std::vector<boost::uint32_t> *bondInvariants) const {
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");

  std::vector<AtomPairAtomEnv> result = std::vector<AtomPairAtomEnv>();
  const double *distanceMatrix;
  if (arguments.use2D) {
    distanceMatrix = MolOps::getDistanceMat(mol);
  } else {
    distanceMatrix = MolOps::get3DDistanceMat(mol, confId);
  }

  const unsigned int atomCount = mol.getNumAtoms();

  std::vector<boost::uint32_t> *atomCodeCache =
      new std::vector<boost::uint32_t>();
  for (ROMol::ConstAtomIterator atomItI = mol.beginAtoms();
       atomItI != mol.endAtoms(); ++atomItI) {
    if (!atomInvariants) {
      atomCodeCache->push_back(
          getAtomCode(*atomItI, 0, arguments.includeChirality));
    } else {
      atomCodeCache->push_back((*atomInvariants)[(*atomItI)->getIdx()] %
                               ((1 << codeSize) - 1));
    }

    for (ROMol::ConstAtomIterator atomItI = mol.beginAtoms();
         atomItI != mol.endAtoms(); ++atomItI) {
      unsigned int i = (*atomItI)->getIdx();
      if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(),
                                   i) != ignoreAtoms->end()) {
        continue;
      }

      for (ROMol::ConstAtomIterator atomItJ = atomItI + 1;
           atomItJ != mol.endAtoms(); ++atomItJ) {
        unsigned int j = (*atomItJ)->getIdx();
        if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(),
                                     j) != ignoreAtoms->end()) {
          continue;
        }

        if (fromAtoms &&
            (std::find(fromAtoms->begin(), fromAtoms->end(), i) ==
             fromAtoms->end()) &&
            (std::find(fromAtoms->begin(), fromAtoms->end(), j) ==
             fromAtoms->end())) {
          continue;
        }
        unsigned int distance =
            static_cast<unsigned int>(floor(distanceMatrix[i * atomCount + j]));

        if (distance >= arguments.minDistance &&
            distance <= arguments.maxDistance) {
          result.push_back(AtomPairAtomEnv(atomCodeCache, i, j, distance));
        }
      }
    }
  }

  if (result.empty()) {
    // cleanUpEnvironments will not be able to clear shared resources since the
    // result is empty so we do it here
    delete atomCodeCache;
  }

  return std::vector<AtomPairAtomEnv>();
};

void AtomPairEnvGenerator::cleanUpEnvironments(
    const std::vector<AtomPairAtomEnv> &atomEnvironments) const {
  if (!atomEnvironments.empty()) {
    delete atomEnvironments[0].getAtomCodeCache();
  }
}

FingerprintGenerator<AtomPairEnvGenerator, AtomPairArguments, AtomPairAtomEnv>
getAtomPairGenerator(const ROMol &mol, const unsigned int minDistance,
                     const unsigned int maxDistance,
                     const bool includeChirality, const bool use2D,
                     const bool useCountSimulation,
                     AtomInvariantsGenerator *atomInvariantsGenerator = 0,
                     BondInvariantsGenerator *bondInvariantsGenerator = 0) {
  AtomPairArguments atomPairArguments = AtomPairArguments(
      useCountSimulation, includeChirality, use2D, minDistance, maxDistance);
  AtomPairEnvGenerator atomPairEnvGenerator = AtomPairEnvGenerator();

  return FingerprintGenerator<AtomPairEnvGenerator, AtomPairArguments,
                              AtomPairAtomEnv>(
      atomPairEnvGenerator, atomPairArguments, atomInvariantsGenerator,
      bondInvariantsGenerator);
}
}  // namespace AtomPair

template<>
FingerprintGenerator<AtomPair::AtomPairEnvGenerator,
                     AtomPair::AtomPairArguments, AtomPair::AtomPairAtomEnv>::
    FingerprintGenerator(
        AtomPair::AtomPairEnvGenerator atomEnvironmentGenerator,
        AtomPair::AtomPairArguments fingerprintArguments,
        AtomInvariantsGenerator *atomInvariantsGenerator,
        BondInvariantsGenerator *bondInvariantsGenerator) {
  this->atomEnvironmentGenerator = atomEnvironmentGenerator;
  this->fingerprintArguments = fingerprintArguments;
  this->asCommonArguments = &this->fingerprintArguments;
  this->atomInvariantsGenerator = atomInvariantsGenerator;
  this->bondInvariantsGenerator = bondInvariantsGenerator;
}
}  // namespace RDKit