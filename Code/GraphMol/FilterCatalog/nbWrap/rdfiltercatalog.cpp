//
//  Copyright (C) 2015-2026 Novartis Institutes for BioMedical Research Inc.
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
void wrap_filtercat(nb::module_ &m);
}

NB_MODULE(rdfiltercatalog, m) {
  m.doc() = "Module containing FilterCatalog functionality for filtering "
            "molecules based on structural patterns.";
  RDKit::wrap_filtercat(m);
}
