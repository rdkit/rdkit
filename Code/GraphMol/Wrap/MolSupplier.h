//
//  Copyright (C) 2005-2019 Greg Landrumm and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_WRAP_MOLSUPPLIER_H
#define RD_WRAP_MOLSUPPLIER_H
//! Template functions for wrapping suppliers as python iterators.

#include <RDBoost/python.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/FileParseException.h>

namespace RDKit {
// Note that this returns a pointer to the supplier itself, so be careful
// that it doesn't get deleted by python!
template <typename T>
T *MolSupplIter(T *suppl) {
  suppl->reset();
  return suppl;
}

template <typename T>
ROMol *MolForwardSupplNext(T *suppl) {
  ROMol *res = nullptr;
  if (!suppl->atEnd()) {
    try {
      res = suppl->next();
    } catch (const FileParseException &) {
      throw;
    } catch (...) {
      res = nullptr;
    }
  } else {
    PyErr_SetString(PyExc_StopIteration, "End of supplier hit");
    throw boost::python::error_already_set();
  }
  if (suppl->atEnd() && suppl->getEOFHitOnRead()) {
    PyErr_SetString(PyExc_StopIteration, "End of supplier hit");
    throw boost::python::error_already_set();
  }
  return res;
}

template <typename T>
ROMol *MolSupplNext(T *suppl) {
  ROMol *res = nullptr;
  if (!suppl->atEnd()) {
    try {
      res = suppl->next();
    } catch (const FileParseException &) {
      throw;
    } catch (...) {
      res = nullptr;
    }
  } else {
    PyErr_SetString(PyExc_StopIteration, "End of supplier hit");
    throw boost::python::error_already_set();
  }

  return res;
}  // namespace RDKit

template <typename T>
ROMol *MolSupplGetItem(T *suppl, int idx) {
  ROMol *res = nullptr;
  if (idx < 0) {
    idx = suppl->length() + idx;
  }
  try {
    res = (*suppl)[idx];
  } catch (...) {
    // it's kind of doofy that we can't just catch the FileParseException
    // that is thrown when we run off the end of the supplier, but it seems
    // that the cross-shared library exception handling thing is just too
    // hard for me as of boost 1.35.0 (also 1.34.1) and g++  4.1.x on linux
    // this approach works as well:
    if (suppl->atEnd()) {
      PyErr_SetString(PyExc_IndexError, "invalid index");
      throw boost::python::error_already_set();
    } else {
      res = nullptr;
    }
  }
  return res;
}
}  // namespace RDKit
#endif
