//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#include "StructChecker.h"

namespace RDKit {
namespace StructureCheck {
// used in unit test
RDKIT_STRUCTCHECKER_EXPORT bool StringToAugmentedAtom(const char *str,
                                                      AugmentedAtom &aa);
}  // namespace StructureCheck
}  // namespace RDKit
