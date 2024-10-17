#pragma once

#include "SmilesWrite.h"

namespace RDKit {

RDKIT_SMILESPARSE_EXPORT void updateSmilesWriteParamsFromJSON(
    SmilesWriteParams &params, const std::string &details_json);
RDKIT_SMILESPARSE_EXPORT void updateSmilesWriteParamsFromJSON(
    SmilesWriteParams &params, const char *details_json);
RDKIT_SMILESPARSE_EXPORT void updateCXSmilesFieldsFromJSON(
    std::uint32_t &cxSmilesFields, unsigned int &restoreBondDirs,
    const std::string &details_json);
RDKIT_SMILESPARSE_EXPORT void updateCXSmilesFieldsFromJSON(
    std::uint32_t &cxSmilesFields, unsigned int &restoreBondDirs,
    const char *details_json);

}  // end namespace RDKit
