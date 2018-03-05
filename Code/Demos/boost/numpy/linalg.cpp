//
//  Copyright (C) 2003 Rational Discovery LLC
//

#include <boost/python.hpp>
#include <RDBoost/boost_numpy.h>
#define PY_ARRAY_UNIQUE_SYMBOL RD_array_API
#include <numpy/arrayobject.h>

namespace python = boost::python;

double GetFirstElement(NumpyArrayType &x) {
  PyArrayObject *ptr = (PyArrayObject *)x.ptr();
  void *data = PyArray_DATA(ptr);
  double res = 0.0;

  switch (ptr->descr->type_num) {
    case PyArray_DOUBLE:
      res = ((double *)data)[0];
      break;
    case PyArray_FLOAT:
      res = (double)((float *)data)[0];
      break;
    case PyArray_LONG:
      res = (double)((long *)data)[0];
      break;
    case PyArray_INT:
      res = (double)((int *)data)[0];
      break;
  }
  return res;
}

BOOST_PYTHON_MODULE(linalg) { python::def("GetFirstElement", GetFirstElement); }
