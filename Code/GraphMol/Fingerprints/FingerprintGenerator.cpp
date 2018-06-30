
#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseBitVect.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <RDGeneral/hash/hash.hpp>
#include <cstdint>

namespace RDKit {

FingerprintArguments::FingerprintArguments(
    const bool countSimulation, const std::vector<std::uint32_t> countBounds,
    std::uint32_t foldedSize)
    : d_countSimulation(countSimulation),
      d_countBounds(countBounds),
      d_foldedSize(foldedSize) {
  PRECONDITION(!countSimulation || !countBounds.empty(),
               "bad count bounds provided");
}

std::string FingerprintArguments::commonArgumentsString() const {
  return "Common arguments\tcountSimulation=" +
         std::to_string(d_countSimulation) +
         " foldedSize=" + std::to_string(d_foldedSize);
}

FingerprintArguments::~FingerprintArguments() {}

AtomEnvironmentGenerator::~AtomEnvironmentGenerator() {}

AtomEnvironment::~AtomEnvironment() {}

AtomInvariantsGenerator::~AtomInvariantsGenerator() {}

BondInvariantsGenerator::~BondInvariantsGenerator() {}

FingerprintGenerator::FingerprintGenerator(
    AtomEnvironmentGenerator *atomEnvironmentGenerator,
    FingerprintArguments *fingerprintArguments,
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

FingerprintGenerator::~FingerprintGenerator() {
  delete dp_atomEnvironmentGenerator;
  delete dp_fingerprintArguments;
  if (df_ownsAtomInvGenerator) {
    delete dp_atomInvariantsGenerator;
  }
  if (df_ownsBondInvGenerator) {
    delete dp_bondInvariantsGenerator;
  }
}

std::string FingerprintGenerator::infoString() const {
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

SparseIntVect<std::uint32_t> *FingerprintGenerator::getFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput) const {
  // create the atom and bond invariants using the generators if there are any,
  // created invariants will be passed to each atom environment's getBitId call
  std::vector<std::uint32_t> *atomInvariants =
      dp_atomInvariantsGenerator
          ? dp_atomInvariantsGenerator->getAtomInvariants(mol)
          : nullptr;
  std::vector<std::uint32_t> *bondInvariants =
      dp_bondInvariantsGenerator
          ? dp_bondInvariantsGenerator->getBondInvariants(mol)
          : nullptr;

  // create all atom environments that will generate the bit-ids that will make
  // up the fingerprint
  std::vector<AtomEnvironment *> atomEnvironments =
      dp_atomEnvironmentGenerator->getEnvironments(
          mol, dp_fingerprintArguments, fromAtoms, ignoreAtoms, confId,
          additionalOutput, atomInvariants, bondInvariants);

  // allocate the result
  SparseIntVect<std::uint32_t> *res = new SparseIntVect<std::uint32_t>(
      dp_fingerprintArguments->getResultSize());

  // iterate over every atom environment and generate bit-ids that will make up
  // the fingerprint
  for (std::vector<AtomEnvironment *>::iterator it = atomEnvironments.begin();
       it != atomEnvironments.end(); it++) {
    std::uint32_t bitId =
        (*it)->getBitId(dp_fingerprintArguments, atomInvariants, bondInvariants,
                        additionalOutput);

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

SparseBitVect *FingerprintGenerator::getFingerprintAsBitVect(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput) const {
  // generate fingerprint using getFingerprint and convert it to SparseBitVect
  SparseIntVect<std::uint32_t> *tempResult =
      getFingerprint(mol, fromAtoms, ignoreAtoms, confId, additionalOutput);
  std::uint32_t countBitsPerBit = dp_fingerprintArguments->d_countBounds.size();

  SparseBitVect *result;

  if (dp_fingerprintArguments->d_countSimulation) {
    std::uint32_t sizeWithCount =
        dp_fingerprintArguments->getResultSize() * countBitsPerBit;
    result = new SparseBitVect(sizeWithCount);
  } else {
    result = new SparseBitVect(dp_fingerprintArguments->getResultSize());
  }

  BOOST_FOREACH (SparseIntVect<std::uint32_t>::StorageType::value_type val,
                 tempResult->getNonzeroElements()) {
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

SparseIntVect<std::uint32_t> *FingerprintGenerator::getFoldedFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput) const {
  // generate fingerprint using getFingerprint and fold it to reduce the size
  SparseIntVect<std::uint32_t> *tempResult =
      getFingerprint(mol, fromAtoms, ignoreAtoms, confId, additionalOutput);
  SparseIntVect<std::uint32_t> *result =
      new SparseIntVect<std::uint32_t>(dp_fingerprintArguments->d_foldedSize);

  BOOST_FOREACH (SparseIntVect<std::uint32_t>::StorageType::value_type val,
                 tempResult->getNonzeroElements()) {
    boost::uint32_t foldedBitId = 0;
    // hash the bit-id before modding to distrubute bits set in the bit-id more
    // evenly
    gboost::hash_combine(foldedBitId, val.first);

    // mod the bit-id to limit the size to the desired amount (d_foldedSize in
    // fingerprint arguments)
    result->setVal(foldedBitId % dp_fingerprintArguments->d_foldedSize,
                   val.second);
  }

  delete tempResult;
  return result;
}

ExplicitBitVect *FingerprintGenerator::getFoldedFingerprintAsBitVect(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput) const {
  // same idea from getFingerprintAsBitVect, generate the fingerprint using
  // getFoldedFingerprint and convert it to a ExplicitBitVect
  SparseIntVect<std::uint32_t> *tempResult = getFoldedFingerprint(
      mol, fromAtoms, ignoreAtoms, confId, additionalOutput);
  std::uint32_t countBitsPerBit = dp_fingerprintArguments->d_countBounds.size();

  ExplicitBitVect *result;

  if (dp_fingerprintArguments->d_countSimulation) {
    std::uint32_t sizeWithCount =
        dp_fingerprintArguments->d_foldedSize * countBitsPerBit;
    result = new ExplicitBitVect(sizeWithCount);
  } else {
    result = new ExplicitBitVect(dp_fingerprintArguments->d_foldedSize);
  }

  BOOST_FOREACH (SparseIntVect<std::uint32_t>::StorageType::value_type val,
                 tempResult->getNonzeroElements()) {
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