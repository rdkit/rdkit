//
//  Copyright (C) 2003-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace nb::literals;

void wrap_maxminpick(nb::module_ &m);
void wrap_leaderpick(nb::module_ &m);
void wrap_HierarchCP(nb::module_ &m);

NB_MODULE(rdSimDivPickers, m) {
  m.doc() = "Module containing the diversity and similarity pickers";

  wrap_maxminpick(m);
  wrap_leaderpick(m);
  wrap_HierarchCP(m);
}
