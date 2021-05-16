//
// Copyright (c) 2010-2018 greg Landrum
//
//   @@ All Rights Reserved  @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// inspired by:
// https://github.com/openbabel/openbabel/blob/master/src/config.h.cmake
#include <RDGeneral/export.h>
/* Version check macro
   Can be used like #if (RDKIT_VERSION >= RDKIT_VERSION_CHECK(2018, 3, 1)) */
#define RDKIT_VERSION_CHECK(year, month, rev) \
  ((year * 1000) + (month * 10) + (rev))

/* RDKIT_VERSION is (year*1000) + (month*10) + (rev) */
#define RDKIT_VERSION \
  RDKIT_VERSION_CHECK(@RDKit_Year @, @RDKit_intMonth @, @RDKit_Revision @)

namespace RDKit {
RDKIT_RDGENERAL_EXPORT extern const char* rdkitVersion;
RDKIT_RDGENERAL_EXPORT extern const char* boostVersion;
RDKIT_RDGENERAL_EXPORT extern const char* rdkitBuild;
}  // namespace RDKit
