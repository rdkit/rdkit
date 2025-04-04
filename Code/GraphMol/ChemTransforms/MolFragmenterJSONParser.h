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

#include "MolFragmenter.h"

namespace RDKit {

//! \brief Parse MolzipParams from JSON.
/*!  The passed MolzipParams instance is updated from
 *   the JSON-parsed content.
 *
 * @param params - molzip parameters
 * @param details_json - JSON string
 */
RDKIT_CHEMTRANSFORMS_EXPORT void parseMolzipParametersJSON(
    MolzipParams &params, const char *details_json);

}  // end namespace RDKit
