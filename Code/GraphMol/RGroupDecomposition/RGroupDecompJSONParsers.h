#pragma once

#include "RGroupDecompParams.h"

namespace RDKit {

RDKIT_RGROUPDECOMPOSITION_EXPORT void
updateRGroupDecompositionParametersFromJSON(
    RGroupDecompositionParameters &params, const std::string &details_json);
RDKIT_RGROUPDECOMPOSITION_EXPORT void
updateRGroupDecompositionParametersFromJSON(
    RGroupDecompositionParameters &params, const char *details_json);

}  // end namespace RDKit
