//
//  Copyright (C) 2003-2019 Greg Landrum and Rational Discovery LLC
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

void wrap_molops();
void wrap_chiralityops();

BOOST_PYTHON_MODULE(rdmolops) {
  python::scope().attr("__doc__") =
      "Module containing RDKit functionality for manipulating molecules.";
  rdkit_import_array();

  // ******************************
  // Functions from MolOps
  //****************************
  wrap_molops();
  wrap_chiralityops();
}
