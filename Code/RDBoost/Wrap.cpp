// $Id: Wrap.cpp 4958 2006-02-17 23:48:22Z glandrum $
//
// Copyright (c) 2003-2006 greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

//
// Generic Wrapper utility functionality
//
#include "Wrap.h"
#include "Exceptions.h"
#include <sstream>
#include <iostream>

// A helper function for dealing with errors. Throw a Python IndexError
void throw_index_error(int key)
{
  PyErr_SetObject(PyExc_IndexError, PyInt_FromLong(key));
  python::throw_error_already_set();
}

// A helper function for dealing with errors. Throw a Python ValueError
void throw_value_error(const std::string err)
{
  PyErr_SetString(PyExc_ValueError, err.c_str());
  python::throw_error_already_set();
}

void translate_index_error(IndexErrorException const&e){
  throw_index_error(e.index());
}

void translate_value_error(ValueErrorException const&e){
  throw_value_error(e.message());
}
