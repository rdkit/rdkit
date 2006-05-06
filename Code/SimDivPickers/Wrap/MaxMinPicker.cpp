// $Id$
//
//  Copyright (C) 2003-2006  Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#define NO_IMPORT_ARRAY

#define PY_ARRAY_UNIQUE_SYMBOL rdpicker_array_API
#include <boost/python.hpp>
#include <RDBoost/Wrap.h>
#include "Numeric/arrayobject.h"


#include <SimDivPickers/DistPicker.h>
#include <SimDivPickers/MaxMinPicker.h>
#include <iostream>


namespace python = boost::python;
namespace RDPickers {
  
  // REVIEW: the poolSize can be pulled from the numeric array
  RDKit::INT_VECT MaxMinPicks(MaxMinPicker *picker, 
                              python::numeric::array &distMat,
                              int poolSize, 
                              int pickSize) {
    // REVIEW: check pickSize < poolSize, otherwise throw_value_error()
    PyArrayObject *copy;
    copy = (PyArrayObject *)PyArray_ContiguousFromObject(distMat.ptr(), 
							 PyArray_DOUBLE, 1,1);
    double *dMat = (double *)copy->data;
    
    RDKit::INT_VECT res=picker->pick(dMat, poolSize, pickSize);
    Py_DECREF(copy);
    return res;
  }
                        
  struct MaxMin_wrap {
    static void wrap() {
      python::class_<MaxMinPicker>("MaxMinPicker", 
                                   "A class for diversity picking of items using the MaxMin Algorithm\n")
        .def("Pick", MaxMinPicks,
             "Pick a subset of items from a pool of items using the MaxMin Algorithm\n"
             "Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604 \n"
             "ARGUMENTS:\n\n"
             "  - distMat: 1D distance matrix (only the lower triangle elements)\n"
             "  - poolSize: number of items in the pool\n"
             "  - pickSize: number of items to pick from the pool\n")
        ;
    };
  };
}

void wrap_maxminpick() {
  RDPickers::MaxMin_wrap::wrap();
}

  
