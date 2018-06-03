
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

  std::vector<boost::uint32_t> atomInvariants =
      this->atomInvariantsGenerator->getAtomInvariants(mol);
  std::vector<boost::uint32_t> bondInvariants =
      this->bondInvariantsGenerator->getBondInvariants(mol);

  for (std::vector<AtomEnvironment *>::iterator it = atomEnvironments.begin();
       it != atomEnvironments.end(); it++) {
    boost::uint32_t bitId = (*it)->getBitId(this->fingerprintArguments,
                                            &atomInvariants, &bondInvariants);

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