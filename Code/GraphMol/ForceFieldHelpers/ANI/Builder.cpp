#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <ForceField/ForceField.h>
#include <ForceField/ANI/AtomicContrib.h>
#include <GraphMol/Descriptors/AtomicEnvironmentVector.h>

#include "Builder.h"

namespace RDKit {
namespace ANI {
using namespace ForceFields::ANI;
namespace Tools {
void addANIContribs(const ROMol &mol, ForceFields::ForceField *field,
                    std::string modelType, unsigned int numLayers,
                    unsigned int ensembleSize, int confId) {
  PRECONDITION(field, "bad ForceField");
  auto conf = mol.getConformer(confId);
  auto numAtoms = mol.getNumAtoms();

  auto speciesVec = RDKit::Descriptors::ANI::GenerateSpeciesVector(mol);

  for (unsigned int i = 0; i < numAtoms; i++) {
    auto atom = conf.getAtomPos(i);
    ANIAtomContrib *ac;
    ac = new ForceFields::ANI::ANIAtomContrib(field, speciesVec(i), i,
                                              speciesVec, numAtoms, numLayers,
                                              ensembleSize, modelType);
    field->contribs().push_back(ForceFields::ContribPtr(ac));
  }
}
}  // namespace Tools

ForceFields::ForceField *constructForceField(ROMol &mol, std::string modelType,
                                             unsigned int ensembleSize,
                                             int confId) {
  auto *res = new ForceFields::ForceField();
  Conformer &conf = mol.getConformer(confId);
  for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    res->positions().push_back(&conf.getAtomPos(i));
  }
  // res->numPoints() = mol.getNumAtoms();
  // res->dimension() = 3;
  Tools::addANIContribs(mol, res, modelType, 0, ensembleSize, confId);
  return res;
}

}  // namespace ANI
}  // namespace RDKit