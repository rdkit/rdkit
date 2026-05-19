//
//  Copyright (C) 2003-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rdinfotheory_array_API
#include <nanobind/nanobind.h>
#include <RDBoost/import_array.h>
#include <ML/InfoTheory/InfoBitRanker.h>
#include <ML/InfoTheory/InfoGainFuncs.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDInfoTheory;

namespace RDInfoTheory {
double infoEntropy(nb::object resArr) {
  PyObject *matObj = resArr.ptr();
  if (!PyArray_Check(matObj)) {
    throw nb::value_error("Expecting a Numeric array object");
  }
  PyArrayObject *copy;
  copy = (PyArrayObject *)PyArray_ContiguousFromObject(
      matObj, PyArray_DESCR((PyArrayObject *)matObj)->type_num, 1, 1);
  double res = 0.0;
  // we are expecting a 1 dimensional array
  auto ncols = (long int)PyArray_DIM((PyArrayObject *)matObj, 0);
  CHECK_INVARIANT(ncols > 0, "");
  if (PyArray_DESCR((PyArrayObject *)matObj)->type_num == NPY_DOUBLE) {
    auto *data = (double *)PyArray_DATA(copy);
    res = InfoEntropy(data, ncols);
  } else if (PyArray_DESCR((PyArrayObject *)matObj)->type_num == NPY_FLOAT) {
    auto *data = (float *)PyArray_DATA(copy);
    res = InfoEntropy(data, ncols);
  } else if (PyArray_DESCR((PyArrayObject *)matObj)->type_num == NPY_INT) {
    int *data = (int *)PyArray_DATA(copy);
    res = InfoEntropy(data, ncols);
  } else if (PyArray_DESCR((PyArrayObject *)matObj)->type_num == NPY_LONG) {
    auto *data = (long int *)PyArray_DATA(copy);
    res = InfoEntropy(data, ncols);
  }
  Py_DECREF(copy);
  return res;
}

double infoGain(nb::object resArr) {
  PyObject *matObj = resArr.ptr();
  if (!PyArray_Check(matObj)) {
    throw nb::value_error("Expecting a Numeric array object");
  }
  PyArrayObject *copy;
  copy = (PyArrayObject *)PyArray_ContiguousFromObject(
      matObj, PyArray_DESCR((PyArrayObject *)matObj)->type_num, 2, 2);
  auto rows = (long int)PyArray_DIM((PyArrayObject *)matObj, 0);
  auto cols = (long int)PyArray_DIM((PyArrayObject *)matObj, 1);
  double res = 0.0;
  if (PyArray_DESCR((PyArrayObject *)matObj)->type_num == NPY_DOUBLE) {
    auto *data = (double *)PyArray_DATA(copy);
    res = InfoEntropyGain(data, rows, cols);
  } else if (PyArray_DESCR((PyArrayObject *)matObj)->type_num == NPY_FLOAT) {
    auto *data = (float *)PyArray_DATA(copy);
    res = InfoEntropyGain(data, rows, cols);
  } else if (PyArray_DESCR((PyArrayObject *)matObj)->type_num == NPY_INT) {
    int *data = (int *)PyArray_DATA(copy);
    res = InfoEntropyGain(data, rows, cols);
  } else if (PyArray_DESCR((PyArrayObject *)matObj)->type_num == NPY_LONG) {
    auto *data = (long int *)PyArray_DATA(copy);
    res = InfoEntropyGain(data, rows, cols);
  } else {
    throw nb::value_error(
        "Numeric array object of type int or long or float or double");
  }
  Py_DECREF(copy);
  return res;
}

double chiSquare(nb::object resArr) {
  PyObject *matObj = resArr.ptr();
  if (!PyArray_Check(matObj)) {
    throw nb::value_error("Expecting a Numeric array object");
  }
  PyArrayObject *copy;
  copy = (PyArrayObject *)PyArray_ContiguousFromObject(
      matObj, PyArray_DESCR((PyArrayObject *)matObj)->type_num, 2, 2);
  auto rows = (long int)PyArray_DIM((PyArrayObject *)matObj, 0);
  auto cols = (long int)PyArray_DIM((PyArrayObject *)matObj, 1);
  double res = 0.0;
  if (PyArray_DESCR((PyArrayObject *)matObj)->type_num == NPY_DOUBLE) {
    auto *data = (double *)PyArray_DATA(copy);
    res = ChiSquare(data, rows, cols);
  } else if (PyArray_DESCR((PyArrayObject *)matObj)->type_num == NPY_FLOAT) {
    auto *data = (float *)PyArray_DATA(copy);
    res = ChiSquare(data, rows, cols);
  } else if (PyArray_DESCR((PyArrayObject *)matObj)->type_num == NPY_INT) {
    int *data = (int *)PyArray_DATA(copy);
    res = ChiSquare(data, rows, cols);
  } else if (PyArray_DESCR((PyArrayObject *)matObj)->type_num == NPY_LONG) {
    auto *data = (long int *)PyArray_DATA(copy);
    res = ChiSquare(data, rows, cols);
  } else {
    throw nb::value_error(
        "Numeric array object of type int or long or float or double");
  }
  Py_DECREF(copy);
  return res;
}
}  // namespace RDInfoTheory

void wrap_ranker(nb::module_ &m);
void wrap_corrmatgen(nb::module_ &m);

NB_MODULE(rdInfoTheory, m) {
  m.doc() =
      "Module containing bunch of functions for information metrics and a "
      "ranker to rank bits";

  rdkit_import_array();

  wrap_ranker(m);
  wrap_corrmatgen(m);

  m.def("InfoEntropy", RDInfoTheory::infoEntropy,
        R"DOC(calculates the informational entropy of the values in an array

ARGUMENTS:

  - resMat: pointer to a long int array containing the data
  - dim: long int containing the length of the _tPtr_ array.

RETURNS:

  a double
)DOC",
        "resArr"_a);

  m.def("InfoGain", RDInfoTheory::infoGain,
        R"DOC(Calculates the information gain for a variable

ARGUMENTS:

  - varMat: a Numeric Array object
    varMat is a Numeric array with the number of possible occurrences
      of each result for reach possible value of the given variable.

    So, for a variable which adopts 4 possible values and a result which
      has 3 possible values, varMat would be 4x3

RETURNS:

  - a Python float object

NOTES

  - this is a dropin replacement for _PyInfoGain()_ in entropy.py
)DOC",
        "resArr"_a);

  m.def("ChiSquare", RDInfoTheory::chiSquare,
        R"DOC(Calculates the chi squared value for a variable

ARGUMENTS:

  - varMat: a Numeric Array object
    varMat is a Numeric array with the number of possible occurrences
      of each result for reach possible value of the given variable.

    So, for a variable which adopts 4 possible values and a result which
      has 3 possible values, varMat would be 4x3

RETURNS:

  - a Python float object
)DOC",
        "resArr"_a);
}
