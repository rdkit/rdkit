//
//  Copyright (C) 2021 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include <RDGeneral/export.h>
#ifdef RDKIT_RDKITCFFI_BUILD
#define RDKIT_RDKITCFFI_EXPORT RDKIT_EXPORT_API
#else
#define RDKIT_RDKITCFFI_EXPORT RDKIT_IMPORT_API
#endif

#ifdef __cplusplus
extern "C" {
#endif
RDKIT_RDKITCFFI_EXPORT char *get_smiles(const char *pkl, size_t len);
RDKIT_RDKITCFFI_EXPORT char *get_mol(const char *input, size_t *len);
RDKIT_RDKITCFFI_EXPORT char *version();
#ifdef __cplusplus
}
#endif
