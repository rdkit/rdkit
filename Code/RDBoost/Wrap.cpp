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
#include "Wrap.h"
#include "pyint_api.h"
#include <sstream>
#include <iostream>

// A helper function for dealing with errors. Throw a Python IndexError
void throw_index_error(int key) {
  PyErr_SetObject(PyExc_IndexError, PyInt_FromLong(key));
  python::throw_error_already_set();
}

// A helper function for dealing with errors. Throw a Python ValueError
void throw_value_error(const std::string err) {
  PyErr_SetString(PyExc_ValueError, err.c_str());
  python::throw_error_already_set();
}

// A helper function for dealing with errors. Throw a Python KeyError
void throw_key_error(const std::string key) {
  PyErr_SetString(PyExc_KeyError, key.c_str());
  python::throw_error_already_set();
}

void translate_index_error(IndexErrorException const& e) {
  throw_index_error(e.index());
}

void translate_value_error(ValueErrorException const& e) {
  throw_value_error(e.what());
}

void translate_key_error(KeyErrorException const& e) {
  throw_key_error(e.key());
}

#ifdef INVARIANT_EXCEPTION_METHOD
// A helper function for dealing with errors. Throw a Python RuntimeError
void throw_runtime_error(const std::string err) {
  PyErr_SetString(PyExc_RuntimeError, err.c_str());
  python::throw_error_already_set();
}

void translate_invariant_error(Invar::Invariant const& e) {
  throw_runtime_error(e.toUserString());
}

boost::dynamic_bitset<> pythonObjectToDynBitset(const python::object &obj,
                                                   boost::dynamic_bitset<>::size_type maxV) {
  boost::dynamic_bitset<> res(maxV);
  if (obj) {
    python::stl_input_iterator<boost::dynamic_bitset<>::size_type> beg(obj), end;
    while (beg != end) {
      auto v = *beg;
      if (v >= maxV) {
        throw_value_error("list element larger than allowed value");
      }
      res.set(v);
      ++beg;
    }
  }
  return res;
}


#endif
