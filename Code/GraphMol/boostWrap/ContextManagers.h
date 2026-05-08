//
//  Copyright (C) 2021 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_WRAP_CONTEXTMGR_H
#define RD_WRAP_CONTEXTMGR_H
//! Template functions for supporting python context managers

#include <RDBoost/python.h>
namespace python = boost::python;

namespace RDKit {
template <typename T>
T *MolIOEnter(T *self) {
  return self;
}

template <typename T>
bool MolIOExit(T *self, python::object exc_type, python::object exc_val,
               python::object traceback) {
  RDUNUSED_PARAM(exc_type);
  RDUNUSED_PARAM(exc_val);
  RDUNUSED_PARAM(traceback);
  self->close();
  return false;
}
}  // namespace RDKit
#endif
