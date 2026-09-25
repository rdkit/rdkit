//
//  Copyright (C) 2013-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>

// ours
#include <RDGeneral/Exceptions.h>
#include <GraphMol/SanitException.h>

namespace nb = nanobind;

namespace RDKit {
void wrap_queries(nb::module_ &m);
}

NB_MODULE(rdqueries, m) {
  m.doc() =
      R"DOC(Module containing RDKit functionality for querying molecules.)DOC";

  RDKit::wrap_queries(m);
}
