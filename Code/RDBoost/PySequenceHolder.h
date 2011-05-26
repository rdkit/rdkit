//
// Copyright (c) 2003-2006 greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef _RD_PYSEQUENCEHOLDER_H_
#define _RD_PYSEQUENCEHOLDER_H_

//
// Defines a class to hold sequences passed in from Python
//
#include "Wrap.h"
#include <RDGeneral/Invariant.h>

namespace python = boost::python;

//! \brief Class to hold sequences (lists, tuples, arrays, etc.)
//!         passed from Python -> C++
//!
//!  PySequenceHolder is templated on the type of the contained object.
//!
//!  The class is \em lazy: elements are not evaluated until requested
//!     within the C++ code.
//!
template <typename T>
class PySequenceHolder {
public:
  PySequenceHolder(python::object seq) {
    d_seq = seq;
  };

  // --------------------------------------------------
  //! \brief Returns the size of the contained sequence.
  //!
  //! NOTE: the sequence must have a \c __len__ attribute, otherwise
  //!       a \c ValueError will be raised.
  unsigned int size() const {
    unsigned int res=0;
    try {
      res = python::extract<int>(d_seq.attr("__len__")());
    } catch (...) {
      throw_value_error("sequence does not support length query");
    }
    return res;
  };

  // --------------------------------------------------
  //! \brief Returns an element of the sequence
  //!
  //! ARGUMENTS:
  //!   - which: an integer specifying which element should be returned. 
  //!  
  //! NOTES:
  //!   - if the sequence is not \a which elements long, we raise an
  //!     \c IndexError
  //!   - if the element cannot be converted to type \c T, we raise a
  //!     \c ValueError
  T operator[](unsigned int which) const {
    if(which > size()){
      throw_index_error(which);
    }

    try{
      T res = python::extract<T>(d_seq[which]);
      return res;
    } catch (...) {
      throw_value_error("cannot extract desired type from sequence");
    }

    POSTCONDITION(0,"cannot reach this point");
    return static_cast<T>(0);
  };
private:
  python::object d_seq;
};


#endif
