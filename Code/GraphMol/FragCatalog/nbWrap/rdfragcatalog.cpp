//
//  Copyright (C) 2003-2026 Rational Discovery LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>

namespace nb = nanobind;

void wrap_fragcat(nb::module_ &m);
void wrap_fragparams(nb::module_ &m);
void wrap_fragcatgen(nb::module_ &m);
void wrap_fragFPgen(nb::module_ &m);

NB_MODULE(rdfragcatalog, m) {
  m.doc() = R"DOC(Module containing FragCatalog, FragCatParams, FragCatGenerator,
and FragFPGenerator classes for fragment-based catalog generation and
fingerprinting.)DOC";
  wrap_fragcat(m);
  wrap_fragparams(m);
  wrap_fragcatgen(m);
  wrap_fragFPgen(m);
}
