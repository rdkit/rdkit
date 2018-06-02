
#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseBitVect.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>

namespace RDKit {

FingerprintArguments::FingerprintArguments(const bool countSimulation)
    : countSimulation(countSimulation) {}

template <class T1, class T2, class T3>
FingerprintGenerator<T1, T2, T3>::FingerprintGenerator(
    T1 atomEnvironmentGenerator, T2 fingerprintArguments,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator) {
  this->atomEnvironmentGenerator = atomEnvironmentGenerator;
  this->fingerprintArguments = fingerprintArguments;
  this->asCommonArguments = &this->fingerprintArguments;
  this->atomInvariantsGenerator = atomInvariantsGenerator;
  this->bondInvariantsGenerator = bondInvariantsGenerator;
}

template <class T1, class T2, class T3>
SparseIntVect<boost::uint32_t>
    *FingerprintGenerator<T1, T2, T3>::getFingerprint(
        const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms,
        const std::vector<boost::uint32_t> *ignoreAtoms, const int confId,
        const AdditionalOutput *additionalOutput) const {
  SparseIntVect<boost::uint32_t> *res = new SparseIntVect<boost::uint32_t>(
      this->asCommonArguments->getResultSize());

  std::vector<T3> atomEnvironments =
      this->atomEnvironmentGenerator.getEnvironments(
          mol, this->fingerprintArguments, fromAtoms, ignoreAtoms, confId,
          additionalOutput);

  std::vector<boost::uint32_t> *atomInvariants =
      this->atomInvariants->getAtomInvariants(mol);
  std::vector<boost::uint32_t> *bondInvariants =
      this->atomInvariants->getBondInvariants(mol);

  for (unsigned int i = 0; i != atomEnvironments.size(); i++) {
    boost::uint32_t bitId = atomEnvironments[i].getBitId(
        this->fingerprintArguments, atomInvariants, bondInvariants);

    if (this->asCommonArguments->countSimulation) {
      res->setVal(bitId, res->getVal(bitId) + 1);
    } else {
      res->setVal(bitId, 1);
    }
  }

  delete atomInvariants;
  delete bondInvariants;

  this->atomEnvironmentGeerator.cleanUpEnvironments(atomEnvironments);

  return res;
}

template <class T1, class T2, class T3>
SparseBitVect *FingerprintGenerator<T1, T2, T3>::getFingerprintAsBitVect(
    const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms,
    const std::vector<boost::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput) const {
  return 0;
}

template <class T1, class T2, class T3>
SparseIntVect<boost::uint32_t>
    *FingerprintGenerator<T1, T2, T3>::getFoldedFingerprint(
        const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms,
        const std::vector<boost::uint32_t> *ignoreAtoms,
        const AdditionalOutput *additionalOutput) const {
  return 0;
}

template <class T1, class T2, class T3>
ExplicitBitVect *
FingerprintGenerator<T1, T2, T3>::getFoldedFingerprintAsBitVect(
    const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms,
    const std::vector<boost::uint32_t> *ignoreAtoms,
    const AdditionalOutput *additionalOutput) const {
  return 0;
}

}  // namespace RDKit