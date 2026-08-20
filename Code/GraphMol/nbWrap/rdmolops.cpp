//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>

#include <RDGeneral/types.h>

#include <RDBoost/import_array.h>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/SanitException.h>
#ifdef RDK_BUILD_THREADSAFE_SSS
#include <mutex>
#endif

namespace nb = nanobind;
using namespace RDKit;

#ifdef RDK_BUILD_THREADSAFE_SSS
static std::once_flag s_rdmolops_numpy_init_flag;
#endif

void rdkit_rdmolops_ensure_numpy() {
#ifdef RDK_BUILD_THREADSAFE_SSS
  std::call_once(s_rdmolops_numpy_init_flag, rdkit_import_array);
#else
  static bool initialized = false;
  if (!initialized) {
    initialized = true;
    rdkit_import_array();
  }
#endif
}

void wrap_molops(nb::module_ &m);
void wrap_chiralityops(nb::module_ &m);

NB_MODULE(rdmolops, m) {
  m.doc() =
      R"DOC(Module containing RDKit functionality for manipulating molecules.)DOC";

  // ******************************
  // Functions from MolOps
  //****************************
  wrap_molops(m);
  wrap_chiralityops(m);
}
