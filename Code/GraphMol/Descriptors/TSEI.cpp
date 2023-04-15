//
// Created by Jason Biggs on 4/15/23.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Descriptors/MolDescriptors.h>

namespace RDKit {
namespace Descriptors {

std::vector<double> calcTSEI(const ROMol& mol, bool force) {
  std::vector<double> res;
  if (!force && mol.getPropIfPresent("_TSEI", res)) {
    return res;
  }
  RDKit::MolOps::getDistanceMat(mol, false, false, force);
  boost::shared_array<int> pathMat;
  mol.getProp(RDKit::common_properties::DistanceMatrix_Paths, pathMat);

  auto nAtoms = mol.getNumAtoms();
  auto table = RDKit::PeriodicTable::getTable();

  for (auto a1Idx = 0; a1Idx < nAtoms; ++a1Idx) {
    auto val = 0.0;
    auto iR = table->getRcovalent(mol.getAtomWithIdx(a1Idx)->getAtomicNum());
    for (auto j = 0; j < nAtoms; ++j) {
      if (pathMat[a1Idx * nAtoms + j] == -1) {
        continue;
      }
      auto distance = 0.0;
      auto jR = table->getRcovalent(mol.getAtomWithIdx(j)->getAtomicNum());
      auto k = j;
      while (k > -1) {
        // for every step on the path add twice the covalent radius
        distance +=
            2 * (table->getRcovalent(mol.getAtomWithIdx(k)->getAtomicNum()));
        k = pathMat[a1Idx * nAtoms + k];
      }
      distance = distance - iR - jR;  // correct for the endpoints
      val += pow(jR, 3) / pow(distance, 3);
    }
    res.push_back(8.0 * val);
  }

  mol.setProp("_TSEI", res, true);

  return res;
}

}  // end of namespace Descriptors
}  // end of namespace RDKit