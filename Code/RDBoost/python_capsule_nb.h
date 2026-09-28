//
// Copyright (C)  2026 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  The idea of inserting a capsule into the builtins dictionary
//  was taken from nanobind, as well as the mechanism to capture
//  the builtins dict object.
//
#pragma once

#include <Python.h>

#include <nanobind/nanobind.h>

namespace nb = nanobind;

inline nb::dict getBuiltinsDict() {
// PyEval_GetBuiltins only gets the builtins for the current thread,
// which may not be the main thread.
#if defined(PYPY_VERSION)
  PyObject *dict = PyEval_GetBuiltins();
#else
  PyObject *dict = PyInterpreterState_GetDict(PyInterpreterState_Get());
#endif
  if (!dict) {
    throw std::runtime_error("Failed to get Python builtins dictionary");
  }

  return nb::borrow<nb::dict>(dict);
}

inline void installCapsule(nb::capsule capsule, const char *name) {
  nb::dict builtins = getBuiltinsDict();
  builtins[name] = capsule;
}

inline void *getCapsulePtr(const char *name) {
  nb::dict builtins = getBuiltinsDict();
  if (!builtins.contains(name)) {
    return nullptr;
  }

  auto capsule = nb::borrow<nb::capsule>(builtins[name]);
  if (!capsule) {
    throw std::runtime_error(std::string("Capsule ") + name +
                             " exists, but we could not retrieve it.");
  }
  return capsule.data();
}