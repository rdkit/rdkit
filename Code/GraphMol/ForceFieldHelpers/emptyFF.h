#include <RDGeneral/export.h>

namespace ForceFields {
class ForceField;
}
namespace RDKit {
class ROMol;
RDKIT_FORCEFIELDHELPERS_EXPORT ForceFields::ForceField *
constructEmptyForceField(ROMol &mol, int confID = -1) {
  auto *res = new ForceFields::ForceField();
  auto &conf = mol.getConformer(confID);
  for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    res->positions().push_back(&(conf.getAtomPos(i)));
  }
  return res;
}
}  // namespace RDKit