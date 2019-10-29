//
//  Copyright (C) 2003-2013 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "rdmolops.h"
#include <RDBoost/python.h>

#include <RDGeneral/types.h>

#include <RDBoost/Wrap.h>
#include <RDBoost/import_array.h>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/SanitException.h>

namespace python = boost::python;
using namespace RDKit;

void rdSanitExceptionTranslator(RDKit::MolSanitizeException const& x) {
  std::ostringstream ss;
  ss << "Sanitization error: " << x.message();
  PyErr_SetString(PyExc_ValueError, ss.str().c_str());
}

void wrap_molops();

BOOST_PYTHON_MODULE(rdmolops) {
  python::scope().attr("__doc__") =
      "Module containing RDKit functionality for manipulating molecules.";
  rdkit_import_array();

  // ******************************
  // Functions from MolOps
  //****************************
  wrap_molops();
}
