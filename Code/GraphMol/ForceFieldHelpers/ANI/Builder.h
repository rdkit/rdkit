#include <RDGeneral/export.h>
#ifndef RD_ANIBUILDER_H_
#define RD_ANIBUILDER_H_

namespace ForceFields {
class ForceField;
}

namespace RDKit {
class ROMol;
namespace ANI {
RDKIT_FORCEFIELDHELPERS_EXPORT ForceFields::ForceField *constructForceField(
    ROMol &mol, std::string modelType, unsigned int ensembleSize,
    int confId = -1);
namespace Tools {
RDKIT_FORCEFIELDHELPERS_EXPORT void addANIContribs(
    const ROMol &mol, ForceFields::ForceField *field, std::string modelType,
    unsigned int numLayers, unsigned int ensembleSize, int confId = -1);
}  // namespace Tools
}  // namespace ANI
}  // namespace RDKit

#endif