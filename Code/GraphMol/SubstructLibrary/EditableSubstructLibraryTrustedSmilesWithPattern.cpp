//
// Created by gareth on 2/10/19.
//

#include "EditableSubstructLibraryTrustedSmilesWithPattern.h"

namespace RDKit {

int EditableSubstructLibraryTrustedSmilesWithPattern::addSmiles(
    const std::string &smi, const std::string &fp) {
  const std::vector<std::string> smiles{smi};
  const std::vector<std::string> fps{fp};
  return addSmiles(smiles, fps);
}

int EditableSubstructLibraryTrustedSmilesWithPattern::addSmiles(
    const std::vector<std::string> smilesVector,
    const std::vector<std::string> fpVector) {
  boost::unique_lock<boost::shared_mutex> lock(access_mutex);

  auto start = ids.size();
  for (const auto &smiStr : smilesVector) {
    auto spacePos = smiStr.find(" ");
    if (spacePos == std::string::npos) {
      throw std::invalid_argument("No structure name in smiles " + smiStr);
    }
    auto smi = smiStr.substr(0, spacePos);
    auto name = smiStr.substr(spacePos+1);
    checkIdNotPresent(name);
    ids.insert(idInfo(ids.size(), name));
    molHolder->addSmiles(smi);
    assert(molHolder->size() == ids.size());
  }
  auto nAdded = ids.size() - start;
  assert(nAdded == smilesVector.size());
  addFingerprints(fpVector);
  return static_cast<int>(nAdded);
}

 ExplicitBitVect *EditableSubstructLibraryTrustedSmilesWithPattern::makeFingerprint(const std::string &smiles) const {
    auto mol = SmilesToMol(smiles, 0, false);
    mol->updatePropertyCache();
    return fpHolder->makeFingerprint(*mol);
  }

  std::string EditableSubstructLibraryTrustedSmilesWithPattern::makeStringFingerprint(const std::string &smiles) const {
    auto fp = makeFingerprint(smiles);
    auto stringFp = fp->toString();
    delete fp;
    return stringFp;
  }
  
  }  // namespace RDKit