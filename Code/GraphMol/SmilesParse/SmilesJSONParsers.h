//
//  Copyright (C) 2024 Novartis Biomedical Research and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

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
