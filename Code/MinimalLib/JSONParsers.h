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

#include <GraphMol/FileParsers/PNGParser.h>
#include <GraphMol/MolOps.h>

namespace RDKit {
namespace MinimalLib {

void updatePropertyPickleOptionsFromJSON(unsigned int &propFlags,
                                         const char *details_json);

void updateSanitizeFlagsFromJSON(unsigned int &sanitizeFlags,
                                 const char *details_json);

void updateRemoveHsParametersFromJSON(MolOps::RemoveHsParameters &ps,
                                      bool &sanitize, const char *details_json);

}  // end namespace MinimalLib
}  // end namespace RDKit
