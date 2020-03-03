//
//  Copyright (C) 2020 Brian P. Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDK_SSSFACTORY
#define RDK_SSSFACTORY
#include <RDGeneral/export.h>
#include "SubstructLibrary.h"
#include <GraphMol/FileParsers/MolSupplier.h>

namespace RDKit {
//! Create pattern fingerprints for the given substructure library
//!  The substructure library must not already have fingepints
/*
     \param sslib The substructure library (without pattern fingerprints)
     \param numThreads the number of threads to use, -1 for all threads [default 1]
*/
RDKIT_SUBSTRUCTLIBRARY_EXPORT void addPatterns(SubstructLibrary &sslib, int numThreads=1);  
}
#endif
