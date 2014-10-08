// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rdinfotheory_array_API
#include <RDBoost/Wrap.h>
#include "numpy/oldnumeric.h"
#include <RDBoost/import_array.h>
#include <ML/InfoTheory/InfoBitRanker.h>
#include <ML/InfoTheory/InfoGainFuncs.h>

namespace python = boost::python;
using namespace RDInfoTheory;

namespace RDInfoTheory {
  double infoEntropy(python::object resArr) {
    PyObject *matObj = resArr.ptr();
    if (!PyArray_Check(matObj)) {
      throw_value_error("Expecting a Numeric array object");
    }
    PyArrayObject *copy;
    copy = (PyArrayObject *)PyArray_ContiguousFromObject(matObj, 
                                                         ((PyArrayObject *)matObj)->descr->type_num,
                                                         1,1);
    double res=0.0;
    // we are expecting a 1 dimensional array
    long int ncols = (long int)((PyArrayObject *)matObj)->dimensions[0];
    CHECK_INVARIANT(ncols > 0, "");
    if (((PyArrayObject *)matObj)->descr->type_num == PyArray_DOUBLE) {
      double *data = (double *)copy->data;
      res = InfoEntropy(data, ncols);
    } else if (((PyArrayObject *)matObj)->descr->type_num == PyArray_FLOAT) {
      float *data = (float *)copy->data;
      res = InfoEntropy(data, ncols);
    } else if (((PyArrayObject *)matObj)->descr->type_num == PyArray_INT) {
      int *data = (int *)copy->data;
      res = InfoEntropy(data, ncols);
    } else if (((PyArrayObject *)matObj)->descr->type_num == PyArray_LONG) {
      long int *data = (long int *)copy->data;
      res = InfoEntropy(data, ncols);
    }
    Py_DECREF(copy);
    return res;
  }
   
  double infoGain(python::object resArr) {
    PyObject *matObj = resArr.ptr();
    if (!PyArray_Check(matObj)) {
      throw_value_error("Expecting a Numeric array object");
    }
    PyArrayObject *copy;
    copy = (PyArrayObject *)PyArray_ContiguousFromObject(matObj, 
                                                         ((PyArrayObject *)matObj)->descr->type_num,
                                                         2,2);
    long int rows = (long int)((PyArrayObject *)matObj)->dimensions[0];
    long int cols = (long int)((PyArrayObject *)matObj)->dimensions[1];
    double res=0.0;
    if (((PyArrayObject *)matObj)->descr->type_num == PyArray_DOUBLE) {
      double *data = (double *)copy->data;
      res = InfoEntropyGain(data, rows, cols);
    } else if (((PyArrayObject *)matObj)->descr->type_num == PyArray_FLOAT) {
      float *data = (float *)copy->data;
      res = InfoEntropyGain(data, rows, cols);
    } else if (((PyArrayObject *)matObj)->descr->type_num == PyArray_INT) {
      int *data = (int *)copy->data;
      res = InfoEntropyGain(data, rows, cols);
    } else if (((PyArrayObject *)matObj)->descr->type_num == PyArray_LONG) {
      long int *data = (long int *)copy->data;
      res = InfoEntropyGain(data, rows, cols);
    } else {
      throw_value_error("Numeric array object of type int or long or float or double");
    }
    Py_DECREF(copy);
    return res;
  }

  double chiSquare(python::object resArr) {
    PyObject *matObj = resArr.ptr();
    if (!PyArray_Check(matObj)) {
      throw_value_error("Expecting a Numeric array object");
    }
    PyArrayObject *copy;
    copy = (PyArrayObject *)PyArray_ContiguousFromObject(matObj, 
                                                         ((PyArrayObject *)matObj)->descr->type_num,
                                                         2,2);
    long int rows = (long int)((PyArrayObject *)matObj)->dimensions[0];
    long int cols = (long int)((PyArrayObject *)matObj)->dimensions[1];
    double res=0.0;
    if (((PyArrayObject *)matObj)->descr->type_num == PyArray_DOUBLE) {
      double *data = (double *)copy->data;
      res = ChiSquare(data, rows, cols);
    } else if (((PyArrayObject *)matObj)->descr->type_num == PyArray_FLOAT) {
      float *data = (float *)copy->data;
      res = ChiSquare(data, rows, cols);
    } else if (((PyArrayObject *)matObj)->descr->type_num == PyArray_INT) {
      int *data = (int *)copy->data;
      res = ChiSquare(data, rows, cols);
    } else if (((PyArrayObject *)matObj)->descr->type_num == PyArray_LONG) {
      long int *data = (long int *)copy->data;
      res = ChiSquare(data, rows, cols);
    } else {
      throw_value_error("Numeric array object of type int or long or float or double");
    }
    Py_DECREF(copy);
    return res;
  }
}

void wrap_ranker();
void wrap_corrmatgen();

BOOST_PYTHON_MODULE(rdInfoTheory)
{
  python::scope().attr("__doc__") =
    "Module containing bunch of functions for information metrics and a ranker to rank bits"
    ;
  
  rdkit_import_array();
  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);

  wrap_ranker();
  wrap_corrmatgen();

  std::string docString="calculates the informational entropy of the values in an array\n\n\
  ARGUMENTS:\n\
    \n\
    - resMat: pointer to a long int array containing the data\n\
    - dim: long int containing the length of the _tPtr_ array.\n\n\
  RETURNS:\n\n\
    a double\n";
  python::def("InfoEntropy", RDInfoTheory::infoEntropy,
              docString.c_str());

  docString="Calculates the information gain for a variable\n\n\
   ARGUMENTS:\n\n\
     - varMat: a Numeric Array object\n\
       varMat is a Numeric array with the number of possible occurances\n\
         of each result for reach possible value of the given variable.\n\n\
       So, for a variable which adopts 4 possible values and a result which\n\
         has 3 possible values, varMat would be 4x3\n\n\
   RETURNS:\n\n\
     - a Python float object\n\n\
   NOTES\n\n\
     - this is a dropin replacement for _PyInfoGain()_ in entropy.py\n";
  python::def("InfoGain", RDInfoTheory::infoGain,
              docString.c_str());


  docString="Calculates the chi squared value for a variable\n\n\
   ARGUMENTS:\n\n\
     - varMat: a Numeric Array object\n\
       varMat is a Numeric array with the number of possible occurances\n\
         of each result for reach possible value of the given variable.\n\n\
       So, for a variable which adopts 4 possible values and a result which\n\
         has 3 possible values, varMat would be 4x3\n\n\
   RETURNS:\n\n\
     - a Python float object\n";
  python::def("ChiSquare", RDInfoTheory::chiSquare,
              docString.c_str());

}

