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

  for (auto idx1 = 0; idx1 < nAtoms; ++idx1) {
    auto val = 0.0;
    auto iR = table->getRcovalent(mol.getAtomWithIdx(idx1)->getAtomicNum());
    for (auto idx2 = 0; idx2 < nAtoms; ++idx2) {
      if (pathMat[idx1 * nAtoms + idx2] == -1) {
        continue;
      }
      auto distance = 0.0;
      auto jR = table->getRcovalent(mol.getAtomWithIdx(idx2)->getAtomicNum());
      auto k = idx2;
      while (k > -1) {
        // for every step on the path add twice the covalent radius
        distance +=
            2 * (table->getRcovalent(mol.getAtomWithIdx(k)->getAtomicNum()));
        k = pathMat[idx1 * nAtoms + k];
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