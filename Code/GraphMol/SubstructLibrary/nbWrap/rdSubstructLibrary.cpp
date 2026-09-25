//  Copyright (c) 2017-2026, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>

namespace nb = nanobind;

namespace RDKit {
void wrap_substructlibrary(nb::module_ &m);
}

NB_MODULE(rdSubstructLibrary, m) { RDKit::wrap_substructlibrary(m); }
