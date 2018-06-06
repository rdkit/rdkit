
#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseBitVect.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <cstdint>

namespace RDKit {

FingerprintArguments::FingerprintArguments(const bool countSimulation)
    : d_countSimulation(countSimulation) {}

FingerprintArguments::~FingerprintArguments() {}

AtomEnvironmentGenerator::~AtomEnvironmentGenerator() {}

AtomEnvironment::~AtomEnvironment() {}

FingerprintGenerator::FingerprintGenerator(
    AtomEnvironmentGenerator *atomEnvironmentGenerator,
    FingerprintArguments *fingerprintArguments,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator) {
  this->dp_atomEnvironmentGenerator = atomEnvironmentGenerator;
  this->dp_fingerprintArguments = fingerprintArguments;
  this->dp_atomInvariantsGenerator = atomInvariantsGenerator;
  this->dp_bondInvariantsGenerator = bondInvariantsGenerator;
}

FingerprintGenerator::~FingerprintGenerator() {
  delete dp_atomEnvironmentGenerator;
  delete dp_fingerprintArguments;
  dp_atomEnvironmentGenerator = nullptr;
  dp_fingerprintArguments = nullptr;
}

SparseIntVect<std::uint32_t> *FingerprintGenerator::getFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput) const {
  SparseIntVect<std::uint32_t> *res = new SparseIntVect<std::uint32_t>(
      dp_fingerprintArguments->getResultSize());

  std::vector<AtomEnvironment *> atomEnvironments =
      dp_atomEnvironmentGenerator->getEnvironments(mol, dp_fingerprintArguments,
                                                   fromAtoms, ignoreAtoms,
                                                   confId, additionalOutput);

  std::vector<std::uint32_t> *atomInvariantsAsPointer = nullptr;
  std::vector<std::uint32_t> atomInvariants;
  if (dp_atomInvariantsGenerator) {
    atomInvariants = dp_atomInvariantsGenerator->getAtomInvariants(
        mol, dp_fingerprintArguments);
    atomInvariantsAsPointer = &atomInvariants;
  }

  std::vector<std::uint32_t> *bondInvariantsAsPointer = nullptr;
  std::vector<std::uint32_t> bondInvariants;
  if (dp_bondInvariantsGenerator) {
    bondInvariants = dp_bondInvariantsGenerator->getBondInvariants(
        mol, dp_fingerprintArguments);
    bondInvariantsAsPointer = &bondInvariants;
  }

  for (std::vector<AtomEnvironment *>::iterator it = atomEnvironments.begin();
       it != atomEnvironments.end(); it++) {
    std::uint32_t bitId =
        (*it)->getBitId(dp_fingerprintArguments, atomInvariantsAsPointer,
                        bondInvariantsAsPointer, additionalOutput);

    if (dp_fingerprintArguments->d_countSimulation) {
      res->setVal(bitId, res->getVal(bitId) + 1);
    } else {
      res->setVal(bitId, 1);
    }
  }

  return res;
}

SparseBitVect *FingerprintGenerator::getFingerprintAsBitVect(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput) const {
  return 0;
}

SparseIntVect<std::uint32_t> *FingerprintGenerator::getFoldedFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms,
    const AdditionalOutput *additionalOutput) const {
  return 0;
}

ExplicitBitVect *FingerprintGenerator::getFoldedFingerprintAsBitVect(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms,
    const AdditionalOutput *additionalOutput) const {
  return 0;
}

}  // namespace RDKit