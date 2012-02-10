//
//  Copyright (C) 2005-2008 Greg Landrumm and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_WRAP_MOLSUPPLIER_H_
#define _RD_WRAP_MOLSUPPLIER_H_
//! Template functions for wrapping suppliers as python iterators.

#include <boost/python.hpp>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/FileParseException.h>

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
    ROMol *res=0;
    if (!suppl->atEnd()) {
      try {
        res=suppl->next();
      } catch(...){
        res=0;
      }
    }
    // FIX: there is an edge case here that we ought to catch:
    //    suppliers where the last molecule has a chemistry problem
    //    With the current behavior, those empty molecules will not
    //    show up in the list
    if(suppl->atEnd() && !res){
      PyErr_SetString(PyExc_StopIteration,"End of supplier hit");
      throw boost::python::error_already_set();
    }
    return res;
  }

  template<typename T>
  ROMol *MolSupplNextAcceptNullLastMolecule(T *suppl){
    ROMol *res=0;
    if (!suppl->atEnd()) {
      try {
        res=suppl->next();
      } catch(...){
        res=0;
      }
    } else {
      PyErr_SetString(PyExc_StopIteration,"End of supplier hit");
      throw boost::python::error_already_set();
    }
    return res;
  }

  template<typename T>
  ROMol *MolSupplGetItem(T *suppl,int idx){
    ROMol *res = 0;
    if(idx<0){
      idx = suppl->length()+idx;
    }
    try{
      res=(*suppl)[idx];
    } catch (...) {
      // it's kind of doofy that we can't just catch the FileParseException
      // that is thrown when we run off the end of the supplier, but it seems
      // that the cross-shared library exception handling thing is just too
      // hard for me as of boost 1.35.0 (also 1.34.1) and g++  4.1.x on linux
      // this approach works as well:
      if(suppl->atEnd()){
        PyErr_SetString(PyExc_IndexError,"invalid index");
        throw boost::python::error_already_set();
      } else {
        res=0;
      }
    }
    return res;
  }

}
#endif
