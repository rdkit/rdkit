//
//  Copyright (C) 2004-2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef __RD_FFPARAMS_H__
#define __RD_FFPARAMS_H__

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace ForceFields {
constexpr double DEG2RAD = M_PI / 180.0;
constexpr double RAD2DEG = 180 / M_PI;
}  // namespace ForceFields
#endif