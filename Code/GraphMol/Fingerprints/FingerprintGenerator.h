#ifndef RD_FINGERPRINTGEN_H_2018_05
#define RD_FINGERPRINTGEN_H_2018_05

#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseBitVect.h>
#include <cstdint>

namespace RDKit {
class ROMol;

struct AdditionalOutput {
  // will review this structure once more fignerprint types are implemented

  std::vector<std::vector<std::uint64_t>> *atomToBits;

  std::map<std::uint32_t, std::vector<std::pair<std::uint32_t, std::uint32_t>>>
      *bitInfoMap;
  // morgan fp
  // maps bitId -> vector of (atomId, radius)

  std::pair<std::vector<std::vector<std::uint32_t>>,
            std::map<std::uint32_t, std::vector<std::vector<int>>>> *bitInfo;
  // rdkit fp
  // first part, vector of bits set for each atom, must have the same size as
  // atom count for molecule
  // second part, maps bitId -> vector of paths

  std::vector<unsigned int> *atomCounts;
  // number of paths that set bits for each atom, must have the same size as
  // atom count for molecule
};

/*!
  /brief Abstract base class that holds molecule independent arguments that are
  common amongst all fingerprint types and classes inherited from this would
  hold fingerprint type specific arguments

 */
template <typename OutputType>
class FingerprintArguments : private boost::noncopyable {
 public:
  FingerprintArguments(const bool countSimulation,
                       const std::vector<std::uint32_t> countBounds,
                       const std::uint32_t foldedSize);
  const bool d_countSimulation;
  const std::vector<std::uint32_t> d_countBounds;
  const std::uint32_t d_foldedSize;

  /*!
    /brief Returns the size of the fingerprint based on arguments

    /return OutputType size of the fingerprint
   */
  virtual OutputType getResultSize() const = 0;

  /**
   /brief method that returns information string about the fingerprint specific
   argument set and the arguments themselves

   /return std::string information string
   */
  virtual std::string infoString() const = 0;

  /**
   /brief method that returns information string about common fingerprinting
   arguments' values

   /return std::string information string
   */
  std::string commonArgumentsString() const;

  virtual ~FingerprintArguments() = 0;
};

/*!
  /brief abstract base class that holds atom-environments that will be hashed to
  generate the fingerprint

 */
template <typename OutputType>
class AtomEnvironment : private boost::noncopyable {
 public:
  /*!
    /brief calculates and returns the bit id to be set for this atom-environment

    /param arguments         Fingerprinting type specific molecule independent
    arguments
    /param atomInvariants    Atom-invariants to be used during hashing
    /param bondInvariants    Bond-invariants to be used during hashing
    /param hashResults   if set results will be ready to be folded by modding

    /return OutputType  calculated bit id for this environment
   */
  virtual OutputType getBitId(FingerprintArguments<OutputType> *arguments,
                              const std::vector<std::uint32_t> *atomInvariants,
                              const std::vector<std::uint32_t> *bondInvariants,
                              const AdditionalOutput *AdditionalOutput,
                              const bool hashResults = false) const = 0;

  virtual ~AtomEnvironment() = 0;
};

/*!
  /brief abstract base class that generates atom-environments from a molecule

 */
template <typename OutputType>
class AtomEnvironmentGenerator : private boost::noncopyable {
 public:
  /*!
    /brief generate and return all atom-envorinments from a molecule

    /param mol               molecule to generate the atom-environments from
    /param arguments         fingerprint type specific molecule independent
    arguments
    /param fromAtoms         atoms to be used during environment generation,
    usage of this parameter depends on the implementation of different
    fingerprint types
    /param ignoreAtoms      atoms to be ignored during environment generation,
    usage of this parameter depends on the implementation of different
    fingerprint types
    /param confId           which conformation to use during environment
    generation, needed for some fingerprint types
    /param additionalOutput contains pointers for additional outputs of
    fingerprinting operation, usage depends on implementation of the fingerprint
    type
    /param atomInvariants   atom invariants to be used during environment
    generation, in some cases some of the hashing can be done during environment
    generation so it is also passed here
    /param bondInvariants   bond invariants to be used during environment
    generation, same as atomInvariants it might be needed
    /param hashResults   if set results will be ready to be folded by modding

    /return std::vector<AtomEnvironment *>  atom-environments generated from
    this molecule
   */
  virtual std::vector<AtomEnvironment<OutputType> *> getEnvironments(
      const ROMol &mol, FingerprintArguments<OutputType> *arguments,
      const std::vector<std::uint32_t> *fromAtoms = nullptr,
      const std::vector<std::uint32_t> *ignoreAtoms = nullptr,
      const int confId = -1, const AdditionalOutput *additionalOutput = nullptr,
      const std::vector<std::uint32_t> *atomInvariants = nullptr,
      const std::vector<std::uint32_t> *bondInvariants = nullptr,
      const bool hashResults = false) const = 0;

  /**
   /brief method that returns information about this /c AtomEnvironmentGenerator
   and its arguments if any

   /return std::string information string
   */
  virtual std::string infoString() const = 0;

  virtual ~AtomEnvironmentGenerator() = 0;
};

/*!
  /brief abstract base class for atom invariants generators

 */
class AtomInvariantsGenerator : private boost::noncopyable {
 public:
  /*!
    /brief get atom invariants from a molecule

    /param mol                   molecule to generate the atom invariants for

    /return std::vector<std::uint32_t> atom invariants generated for the given
    molecule
   */
  virtual std::vector<std::uint32_t> *getAtomInvariants(
      const ROMol &mol) const = 0;

  /**
   /brief method that returns information about this /c AtomInvariantsGenerator
   and its arguments

   /return std::string information string
   */
  virtual std::string infoString() const = 0;

  virtual ~AtomInvariantsGenerator() = 0;
  virtual AtomInvariantsGenerator *clone() const = 0;
};

/*!
  /brief abstract base class for bond invariants generators

 */
class BondInvariantsGenerator : private boost::noncopyable {
 public:
  /*!
    /brief get bond invariants from a molecule

    /param mol                   molecule to generate the bond invariants for

    /return std::vector<std::uint32_t> bond invariants generated for the given
    molecule
   */
  virtual std::vector<std::uint32_t> *getBondInvariants(
      const ROMol &mol) const = 0;

  /**
 /brief method that returns information about this /c BondInvariantsGenerator
 and its arguments

 /return std::string information string
 */
  virtual std::string infoString() const = 0;

  virtual ~BondInvariantsGenerator() = 0;
  virtual BondInvariantsGenerator *clone() const = 0;
};

/*!
  /brief class that generates same fingerprint style for different output
  formats

 */
template <typename OutputType>
class FingerprintGenerator : private boost::noncopyable {
  FingerprintArguments<OutputType> *dp_fingerprintArguments;
  AtomEnvironmentGenerator<OutputType> *dp_atomEnvironmentGenerator;
  AtomInvariantsGenerator *dp_atomInvariantsGenerator;
  BondInvariantsGenerator *dp_bondInvariantsGenerator;
  const bool df_ownsAtomInvGenerator;
  const bool df_ownsBondInvGenerator;

 public:
  FingerprintGenerator(
      AtomEnvironmentGenerator<OutputType> *atomEnvironmentGenerator,
      FingerprintArguments<OutputType> *fingerprintArguments,
      AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
      BondInvariantsGenerator *bondInvariantsGenerator = nullptr,
      bool ownsAtomInvGenerator = false, bool ownsBondInvGenerator = false);

  ~FingerprintGenerator();

  SparseIntVect<OutputType> *getFingerprint(
      const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms = nullptr,
      const std::vector<std::uint32_t> *ignoreAtoms = nullptr,
      const int confId = -1, const AdditionalOutput *additionalOutput = nullptr,
      const std::vector<std::uint32_t> *customAtomInvariants = nullptr,
      const std::vector<std::uint32_t> *customBondInvariants = nullptr) const;

  SparseBitVect *getFingerprintAsBitVect(
      const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms = nullptr,
      const std::vector<std::uint32_t> *ignoreAtoms = nullptr,
      const int confId = -1, const AdditionalOutput *additionalOutput = nullptr,
      const std::vector<std::uint32_t> *customAtomInvariants = nullptr,
      const std::vector<std::uint32_t> *customBondInvariants = nullptr) const;

  SparseIntVect<OutputType> *getFoldedFingerprint(
      const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms = nullptr,
      const std::vector<std::uint32_t> *ignoreAtoms = nullptr,
      const int confId = -1, const AdditionalOutput *additionalOutput = nullptr,
      const std::vector<std::uint32_t> *customAtomInvariants = nullptr,
      const std::vector<std::uint32_t> *customBondInvariants = nullptr) const;

  ExplicitBitVect *getFoldedFingerprintAsBitVect(
      const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms = nullptr,
      const std::vector<std::uint32_t> *ignoreAtoms = nullptr,
      const int confId = -1, const AdditionalOutput *additionalOutput = nullptr,
      const std::vector<std::uint32_t> *customAtomInvariants = nullptr,
      const std::vector<std::uint32_t> *customBondInvariants = nullptr) const;

  std::string infoString() const;
};
}  // namespace RDKit

#include <RDGeneral/hash/hash.hpp>
namespace RDKit {

template <typename OutputType>
FingerprintArguments<OutputType>::FingerprintArguments(
    const bool countSimulation, const std::vector<std::uint32_t> countBounds,
    std::uint32_t foldedSize)
    : d_countSimulation(countSimulation),
      d_countBounds(countBounds),
      d_foldedSize(foldedSize) {
  PRECONDITION(!countSimulation || !countBounds.empty(),
               "bad count bounds provided");
}

template <typename OutputType>
std::string FingerprintArguments<OutputType>::commonArgumentsString() const {
  return "Common arguments\tcountSimulation=" +
         std::to_string(d_countSimulation) +
         " foldedSize=" + std::to_string(d_foldedSize);
}

template <typename OutputType>
FingerprintArguments<OutputType>::~FingerprintArguments() {}

template <typename OutputType>
AtomEnvironmentGenerator<OutputType>::~AtomEnvironmentGenerator() {}

template <typename OutputType>
AtomEnvironment<OutputType>::~AtomEnvironment() {}

inline AtomInvariantsGenerator::~AtomInvariantsGenerator() {}
inline AtomInvariantsGenerator *AtomInvariantsGenerator::clone() const {}

inline BondInvariantsGenerator::~BondInvariantsGenerator() {}
inline BondInvariantsGenerator *BondInvariantsGenerator::clone() const {}

template <typename OutputType>
FingerprintGenerator<OutputType>::FingerprintGenerator(
    AtomEnvironmentGenerator<OutputType> *atomEnvironmentGenerator,
    FingerprintArguments<OutputType> *fingerprintArguments,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator, bool ownsAtomInvGenerator,
    bool ownsBondInvGenerator)
    : df_ownsAtomInvGenerator(ownsAtomInvGenerator),
      df_ownsBondInvGenerator(ownsBondInvGenerator) {
  this->dp_atomEnvironmentGenerator = atomEnvironmentGenerator;
  this->dp_fingerprintArguments = fingerprintArguments;
  this->dp_atomInvariantsGenerator = atomInvariantsGenerator;
  this->dp_bondInvariantsGenerator = bondInvariantsGenerator;
}

template <typename OutputType>
FingerprintGenerator<OutputType>::~FingerprintGenerator() {
  delete dp_atomEnvironmentGenerator;
  delete dp_fingerprintArguments;
  if (df_ownsAtomInvGenerator) {
    delete dp_atomInvariantsGenerator;
  }
  if (df_ownsBondInvGenerator) {
    delete dp_bondInvariantsGenerator;
  }
}

template <typename OutputType>
std::string FingerprintGenerator<OutputType>::infoString() const {
  std::string seperator = " : ";
  return dp_fingerprintArguments->commonArgumentsString() + seperator +
         dp_fingerprintArguments->infoString() + seperator +
         dp_atomEnvironmentGenerator->infoString() + seperator +
         (dp_atomInvariantsGenerator
              ? (dp_atomInvariantsGenerator->infoString() + seperator)
              : "No atom invariants generator") +
         (dp_bondInvariantsGenerator
              ? (dp_bondInvariantsGenerator->infoString() + seperator)
              : "No bond invariants generator");
}

template <typename OutputType>
SparseIntVect<OutputType> *FingerprintGenerator<OutputType>::getFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *customAtomInvariants,
    const std::vector<std::uint32_t> *customBondInvariants) const {
  // create the atom and bond invariants using the generators if there are any,
  // created invariants will be passed to each atom environment's getBitId call
  std::vector<std::uint32_t> *atomInvariants = nullptr;
  if (customAtomInvariants) {
    atomInvariants = new std::vector<std::uint32_t>(*customAtomInvariants);
  } else if (dp_atomInvariantsGenerator) {
    atomInvariants = dp_atomInvariantsGenerator->getAtomInvariants(mol);
  }

  std::vector<std::uint32_t> *bondInvariants = nullptr;
  if (customBondInvariants) {
    bondInvariants = new std::vector<std::uint32_t>(*customBondInvariants);
  } else if (dp_bondInvariantsGenerator) {
    bondInvariants = dp_bondInvariantsGenerator->getBondInvariants(mol);
  }

  // create all atom environments that will generate the bit-ids that will make
  // up the fingerprint
  std::vector<AtomEnvironment<OutputType> *> atomEnvironments =
      dp_atomEnvironmentGenerator->getEnvironments(
          mol, dp_fingerprintArguments, fromAtoms, ignoreAtoms, confId,
          additionalOutput, atomInvariants, bondInvariants);

  // allocate the result
  SparseIntVect<OutputType> *res =
      new SparseIntVect<OutputType>(dp_fingerprintArguments->getResultSize());

  // iterate over every atom environment and generate bit-ids that will make up
  // the fingerprint
  for (auto it = atomEnvironments.begin(); it != atomEnvironments.end(); it++) {
    OutputType bitId = (*it)->getBitId(dp_fingerprintArguments, atomInvariants,
                                       bondInvariants, additionalOutput);
    if (dp_fingerprintArguments->d_countSimulation) {
      // keep the occurrence count for every bit generated
      res->setVal(bitId, res->getVal(bitId) + 1);
    } else {
      // do not keep the count, just set to 1
      res->setVal(bitId, 1);
    }

    delete (*it);
  }

  delete atomInvariants;
  delete bondInvariants;

  return res;
}

template <typename OutputType>
SparseBitVect *FingerprintGenerator<OutputType>::getFingerprintAsBitVect(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *customAtomInvariants,
    const std::vector<std::uint32_t> *customBondInvariants) const {
  // generate fingerprint using getFingerprint and convert it to SparseBitVect
  SparseIntVect<OutputType> *tempResult =
      getFingerprint(mol, fromAtoms, ignoreAtoms, confId, additionalOutput,
                     customAtomInvariants, customBondInvariants);
  OutputType countBitsPerBit = dp_fingerprintArguments->d_countBounds.size();

  SparseBitVect *result;

  if (dp_fingerprintArguments->d_countSimulation) {
    OutputType sizeWithCount =
        dp_fingerprintArguments->getResultSize() * countBitsPerBit;
    result = new SparseBitVect(sizeWithCount);
  } else {
    result = new SparseBitVect(dp_fingerprintArguments->getResultSize());
  }

  BOOST_FOREACH (auto val, tempResult->getNonzeroElements()) {
    if (dp_fingerprintArguments->d_countSimulation) {
      for (unsigned int i = 0; i < countBitsPerBit; ++i) {
        // for every bound in the d_countBounds in dp_fingerprintArguments, set
        // a bit if the occurrence count is equal or higher than the bound for
        // that bit
        if (val.second >= dp_fingerprintArguments->d_countBounds[i]) {
          result->setBit(val.first * countBitsPerBit + i);
        }
      }
    } else {
      result->setBit(val.first);
    }
  }

  delete tempResult;
  return result;
}

template <typename OutputType>
SparseIntVect<OutputType>
    *FingerprintGenerator<OutputType>::getFoldedFingerprint(
        const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
        const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
        const AdditionalOutput *additionalOutput,
        const std::vector<std::uint32_t> *customAtomInvariants,
        const std::vector<std::uint32_t> *customBondInvariants) const {
  // create the atom and bond invariants using the generators if there are any,
  // created invariants will be passed to each atom environment's getBitId call
  std::vector<std::uint32_t> *atomInvariants = nullptr;
  if (customAtomInvariants) {
    atomInvariants = new std::vector<std::uint32_t>(*customAtomInvariants);
  } else if (dp_atomInvariantsGenerator) {
    atomInvariants = dp_atomInvariantsGenerator->getAtomInvariants(mol);
  }

  std::vector<std::uint32_t> *bondInvariants = nullptr;
  if (customBondInvariants) {
    bondInvariants = new std::vector<std::uint32_t>(*customBondInvariants);
  } else if (dp_bondInvariantsGenerator) {
    bondInvariants = dp_bondInvariantsGenerator->getBondInvariants(mol);
  }

  // create all atom environments that will generate the bit-ids that will make
  // up the fingerprint
  std::vector<AtomEnvironment<OutputType> *> atomEnvironments =
      dp_atomEnvironmentGenerator->getEnvironments(
          mol, dp_fingerprintArguments, fromAtoms, ignoreAtoms, confId,
          additionalOutput, atomInvariants, bondInvariants, true);

  // allocate the result
  SparseIntVect<OutputType> *res =
      new SparseIntVect<OutputType>(dp_fingerprintArguments->d_foldedSize);

  // iterate over every atom environment and generate bit-ids that will make up
  // the fingerprint
  for (auto it = atomEnvironments.begin(); it != atomEnvironments.end(); it++) {
    OutputType bitId = (*it)->getBitId(dp_fingerprintArguments, atomInvariants,
                                       bondInvariants, additionalOutput, true);

    bitId = bitId % dp_fingerprintArguments->d_foldedSize;

    if (dp_fingerprintArguments->d_countSimulation) {
      // keep the occurrence count for every bit generated
      res->setVal(bitId, res->getVal(bitId) + 1);
    } else {
      // do not keep the count, just set to 1
      res->setVal(bitId, 1);
    }

    delete (*it);
  }

  delete atomInvariants;
  delete bondInvariants;

  return res;
}

template <typename OutputType>
ExplicitBitVect *
FingerprintGenerator<OutputType>::getFoldedFingerprintAsBitVect(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *customAtomInvariants,
    const std::vector<std::uint32_t> *customBondInvariants) const {
  // same idea from getFingerprintAsBitVect, generate the fingerprint using
  // getFoldedFingerprint and convert it to a ExplicitBitVect
  SparseIntVect<OutputType> *tempResult = getFoldedFingerprint(
      mol, fromAtoms, ignoreAtoms, confId, additionalOutput,
      customAtomInvariants, customBondInvariants);
  OutputType countBitsPerBit = dp_fingerprintArguments->d_countBounds.size();

  ExplicitBitVect *result;

  if (dp_fingerprintArguments->d_countSimulation) {
    OutputType sizeWithCount =
        dp_fingerprintArguments->d_foldedSize * countBitsPerBit;
    result = new ExplicitBitVect(sizeWithCount);
  } else {
    result = new ExplicitBitVect(dp_fingerprintArguments->d_foldedSize);
  }

  BOOST_FOREACH (auto val, tempResult->getNonzeroElements()) {
    // same count simulation logic used in getFingerprintAsBitVect
    if (dp_fingerprintArguments->d_countSimulation) {
      for (unsigned int i = 0; i < countBitsPerBit; ++i) {
        if (val.second >= dp_fingerprintArguments->d_countBounds[i]) {
          result->setBit(val.first * countBitsPerBit + i);
        }
      }
    } else {
      result->setBit(val.first);
    }
  }

  delete tempResult;
  return result;
}

}  // namespace RDKit

#endif
