
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

SparseIntVect<std::uint32_t> *FingerprintGenerator::getFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput) const {
  std::vector<AtomEnvironment *> atomEnvironments =
      dp_atomEnvironmentGenerator->getEnvironments(mol, dp_fingerprintArguments,
                                                   fromAtoms, ignoreAtoms,
                                                   confId, additionalOutput);

  std::vector<std::uint32_t> *atomInvariants =
      dp_atomInvariantsGenerator
          ? dp_atomInvariantsGenerator->getAtomInvariants(mol)
          : nullptr;

  std::vector<std::uint32_t> *bondInvariants =
      dp_bondInvariantsGenerator
          ? dp_bondInvariantsGenerator->getBondInvariants(mol)
          : nullptr;

  SparseIntVect<std::uint32_t> *res = new SparseIntVect<std::uint32_t>(
      dp_fingerprintArguments->getResultSize());

  for (std::vector<AtomEnvironment *>::iterator it = atomEnvironments.begin();
       it != atomEnvironments.end(); it++) {
    std::uint32_t bitId =
        (*it)->getBitId(dp_fingerprintArguments, atomInvariants, bondInvariants,
                        additionalOutput);

    if (dp_fingerprintArguments->d_countSimulation) {
      res->setVal(bitId, res->getVal(bitId) + 1);
    } else {
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
  return nullptr;
}

SparseIntVect<std::uint32_t> *FingerprintGenerator::getFoldedFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms,
    const AdditionalOutput *additionalOutput) const {
  return nullptr;
}

ExplicitBitVect *FingerprintGenerator::getFoldedFingerprintAsBitVect(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms,
    const AdditionalOutput *additionalOutput) const {
  return nullptr;
}

}  // namespace RDKit