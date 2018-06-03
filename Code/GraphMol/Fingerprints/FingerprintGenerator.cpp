
#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseBitVect.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>

namespace RDKit {

FingerprintArguments::FingerprintArguments(const bool countSimulation)
    : countSimulation(countSimulation) {}

FingerprintArguments::~FingerprintArguments() {}

AtomEnvironmentGenerator::~AtomEnvironmentGenerator() {}

AtomEnvironment::~AtomEnvironment() {}

FingerprintGenerator::FingerprintGenerator(
    AtomEnvironmentGenerator *atomEnvironmentGenerator,
    FingerprintArguments *fingerprintArguments,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator) {
  this->atomEnvironmentGenerator = atomEnvironmentGenerator;
  this->fingerprintArguments = fingerprintArguments;
  this->atomInvariantsGenerator = atomInvariantsGenerator;
  this->bondInvariantsGenerator = bondInvariantsGenerator;
}

void FingerprintGenerator::cleanUpResources() {
  delete atomEnvironmentGenerator;
  delete fingerprintArguments;
  atomEnvironmentGenerator = 0;
  fingerprintArguments = 0;
}

SparseIntVect<boost::uint32_t> *FingerprintGenerator::getFingerprint(
    const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms,
    const std::vector<boost::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput) const {
  SparseIntVect<boost::uint32_t> *res = new SparseIntVect<boost::uint32_t>(
      this->fingerprintArguments->getResultSize());

  std::vector<AtomEnvironment *> atomEnvironments =
      this->atomEnvironmentGenerator->getEnvironments(
          mol, this->fingerprintArguments, fromAtoms, ignoreAtoms, confId,
          additionalOutput);

  std::vector<boost::uint32_t> *atomInvariantsAsPointer = 0;
  std::vector<boost::uint32_t> atomInvariants;
  if (this->atomInvariantsGenerator) {
    atomInvariants = this->atomInvariantsGenerator->getAtomInvariants(mol);
    atomInvariantsAsPointer = &atomInvariants;
  }

  std::vector<boost::uint32_t> *bondInvariantsAsPointer = 0;
  std::vector<boost::uint32_t> bondInvariants;
  if (this->atomInvariantsGenerator) {
    bondInvariants = this->bondInvariantsGenerator->getBondInvariants(mol);
    bondInvariantsAsPointer = &bondInvariants;
  }

  for (std::vector<AtomEnvironment *>::iterator it = atomEnvironments.begin();
       it != atomEnvironments.end(); it++) {
    boost::uint32_t bitId =
        (*it)->getBitId(this->fingerprintArguments, atomInvariantsAsPointer,
                        bondInvariantsAsPointer);

    if (this->fingerprintArguments->countSimulation) {
      res->setVal(bitId, res->getVal(bitId) + 1);
    } else {
      res->setVal(bitId, 1);
    }
  }

  this->atomEnvironmentGenerator->cleanUpEnvironments(atomEnvironments);

  for (std::vector<AtomEnvironment *>::iterator it = atomEnvironments.begin();
       it != atomEnvironments.end(); it++) {
    delete (*it);
  }

  return res;
}

SparseBitVect *FingerprintGenerator::getFingerprintAsBitVect(
    const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms,
    const std::vector<boost::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput) const {
  return 0;
}

SparseIntVect<boost::uint32_t> *FingerprintGenerator::getFoldedFingerprint(
    const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms,
    const std::vector<boost::uint32_t> *ignoreAtoms,
    const AdditionalOutput *additionalOutput) const {
  return 0;
}

ExplicitBitVect *FingerprintGenerator::getFoldedFingerprintAsBitVect(
    const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms,
    const std::vector<boost::uint32_t> *ignoreAtoms,
    const AdditionalOutput *additionalOutput) const {
  return 0;
}

}  // namespace RDKit