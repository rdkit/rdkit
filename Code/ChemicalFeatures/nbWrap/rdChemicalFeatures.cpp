//
//  Copyright (C) 2003-2026 Rational Discovery LLC and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>

namespace nb = nanobind;

void wrap_freefeat(nb::module_ &m);

NB_MODULE(rdChemicalFeatures, m) {
  m.doc() = R"DOC(Module containing free chemical feature functionality
These are features that are not associated with molecules. They are
typically derived from pharmacophores and site-maps.
)DOC";

  wrap_freefeat(m);
}
