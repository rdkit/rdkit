//
//  Copyright (C) 2020 Manan Goel
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RD_ANIBUILDER_H_
#define RD_ANIBUILDER_H_

namespace ForceFields {
class ForceField;
}

namespace RDKit {
class ROMol;
namespace ANI {

//! Builds and returns an ANI force field for a molecule
/*!
  \param mol              the molecule to use
  \param modelType        type of model used
  \param ensembleSize     number of models in the ensemble
  \param confId           Conformer ID

  \return the new force field
*/
RDKIT_FORCEFIELDHELPERS_EXPORT ForceFields::ForceField *constructForceField(
    ROMol &mol, std::string modelType, unsigned int ensembleSize,
    int confId = -1);
namespace Tools {
//! Adds ANI style contribution to the ForceField
/*!
  \param mol              the molecule to use
  \param field            Pointer to the ForceField object
  \param modelType        type of model used
  \param numLayers        Number of layers in the Neural Network
  \param ensembleSize     number of models in the ensemble
  \param confId           Conformer ID
*/
RDKIT_FORCEFIELDHELPERS_EXPORT void addANIContribs(
    const ROMol &mol, ForceFields::ForceField *field, std::string modelType,
    unsigned int numLayers, unsigned int ensembleSize, int confId = -1);
}  // namespace Tools
}  // namespace ANI
}  // namespace RDKit

#endif