//
//  Copyright (C) 2005-2008 Greg Landrumm and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_WRAP_MOLSUPPLIER_H_
#define _RD_WRAP_MOLSUPPLIER_H_
//! Template functions for wrapping suppliers as python iterators.

#include <boost/python.hpp>
#include <GraphMol/RDKitBase.h>

namespace RDKit {
  // Note that this returns a pointer to the supplier itself, so be careful
  // that it doesn't get deleted by python!
  template<typename T>
  T *MolSupplIter(T *suppl){
    suppl->reset();
    return suppl;
  }

  template<typename T>
  ROMol *MolSupplNext(T *suppl){
    if (!suppl->atEnd()) {
      return suppl->next();
    }
    else {
      PyErr_SetString(PyExc_StopIteration,"End of supplier hit");
      throw boost::python::error_already_set();
    }
    return 0;
  }
  template<typename T>
  ROMol *MolSupplGetItem(T *suppl,int idx){
    ROMol *res = (*suppl)[idx];
    return res;
  }

}
#endif
