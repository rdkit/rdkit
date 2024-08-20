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
bool MolIOExit(T *self, [[maybe_unused]] python::object exc_type,
               [[maybe_unused]] python::object exc_val,
               [[maybe_unused]] python::object traceback) {
  self->close();
  return false;
}
}  // namespace RDKit
#endif
