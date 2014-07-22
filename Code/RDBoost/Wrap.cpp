// $Id$
//
// Copyright (c) 2003-2006 greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

//
// Generic Wrapper utility functionality
//
#define RDKIT_WRAP_SOURCE
#include "Wrap.h"
#include "pyint_api.h"
#include "Exceptions.h"
#include <sstream>
#include <iostream>

// A helper function for dealing with errors. Throw a Python IndexError
RDKIT_WRAP_DECL void throw_index_error(int key)
{
  PyErr_SetObject(PyExc_IndexError, PyInt_FromLong(key));
  python::throw_error_already_set();
}

// A helper function for dealing with errors. Throw a Python ValueError
RDKIT_WRAP_DECL void throw_value_error(const std::string err)
{
  PyErr_SetString(PyExc_ValueError, err.c_str());
  python::throw_error_already_set();
}

RDKIT_WRAP_DECL void translate_index_error(IndexErrorException const&e){
  throw_index_error(e.index());
}

RDKIT_WRAP_DECL void translate_value_error(ValueErrorException const&e){
  throw_value_error(e.message());
}

