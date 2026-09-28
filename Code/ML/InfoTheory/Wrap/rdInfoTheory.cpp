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
#include <RDBoost/import_array.h>
#include <ML/InfoTheory/InfoBitRanker.h>
#include <ML/InfoTheory/InfoGainFuncs.h>

namespace python = boost::python;
using namespace RDInfoTheory;

namespace RDInfoTheory {
namespace {
// InfoEntropy(), InfoEntropyGain() and ChiSquare() are instantiated for these
// dtypes, so arrays using them are read directly. Anything else is converted to
// doubles. Note that NPY_LONG is only 32 bits wide on Windows, where numpy's
// default integer is NPY_LONGLONG (github #8421).
int nativeOrDoubleType(const python::object &resArr) {
  if (!PyArray_Check(resArr.ptr())) {
    throw_value_error("Expecting a Numeric array object");
  }
  int typeNum = PyArray_DESCR((PyArrayObject *)resArr.ptr())->type_num;
  switch (typeNum) {
    case NPY_DOUBLE:
    case NPY_FLOAT:
    case NPY_INT:
    case NPY_LONG:
    case NPY_LONGLONG:
      return typeNum;
    default:
      return NPY_DOUBLE;
  }
}

PyArrayObject *contiguousArray(const python::object &resArr, int typeNum,
                               int minDim, int maxDim) {
  auto *copy = (PyArrayObject *)PyArray_ContiguousFromObject(
      resArr.ptr(), typeNum, minDim, maxDim);
  if (!copy) {
    throw_value_error("could not convert argument to a numeric array");
  }
  return copy;
}

// Calls func with the array's data cast to the pointer type matching typeNum.
template <typename Func>
double withTypedData(PyArrayObject *arr, int typeNum, Func func) {
  void *data = PyArray_DATA(arr);
  switch (typeNum) {
    case NPY_FLOAT:
      return func((float *)data);
    case NPY_INT:
      return func((int *)data);
    case NPY_LONG:
      return func((long int *)data);
    case NPY_LONGLONG:
      return func((long long *)data);
    default:
      return func((double *)data);
  }
}
}  // namespace

double infoEntropy(python::object resArr) {
  int typeNum = nativeOrDoubleType(resArr);
  PyArrayObject *copy = contiguousArray(resArr, typeNum, 1, 1);
  // we are expecting a 1 dimensional array
  auto ncols = (long int)PyArray_DIM(copy, 0);
  CHECK_INVARIANT(ncols > 0, "");
  double res = withTypedData(
      copy, typeNum, [ncols](auto *data) { return InfoEntropy(data, ncols); });
  Py_DECREF(copy);
  return res;
}

double infoGain(python::object resArr) {
  int typeNum = nativeOrDoubleType(resArr);
  PyArrayObject *copy = contiguousArray(resArr, typeNum, 2, 2);
  auto rows = (long int)PyArray_DIM(copy, 0);
  auto cols = (long int)PyArray_DIM(copy, 1);
  double res = withTypedData(copy, typeNum, [rows, cols](auto *data) {
    return InfoEntropyGain(data, rows, cols);
  });
  Py_DECREF(copy);
  return res;
}

double chiSquare(python::object resArr) {
  int typeNum = nativeOrDoubleType(resArr);
  PyArrayObject *copy = contiguousArray(resArr, typeNum, 2, 2);
  auto rows = (long int)PyArray_DIM(copy, 0);
  auto cols = (long int)PyArray_DIM(copy, 1);
  double res = withTypedData(copy, typeNum, [rows, cols](auto *data) {
    return ChiSquare(data, rows, cols);
  });
  Py_DECREF(copy);
  return res;
}
}  // namespace RDInfoTheory

void wrap_ranker();
void wrap_corrmatgen();

BOOST_PYTHON_MODULE(rdInfoTheory) {
  python::scope().attr("__doc__") =
      "Module containing bunch of functions for information metrics and a "
      "ranker to rank bits";

  rdkit_import_array();

  wrap_ranker();
  wrap_corrmatgen();

  std::string docString =
      "calculates the informational entropy of the values in an array\n\n\
  ARGUMENTS:\n\
    \n\
    - resMat: pointer to a long int array containing the data\n\
    - dim: long int containing the length of the _tPtr_ array.\n\n\
  RETURNS:\n\n\
    a double\n";
  python::def("InfoEntropy", RDInfoTheory::infoEntropy, docString.c_str(),
              python::args("resArr"));

  docString =
      "Calculates the information gain for a variable\n\n\
   ARGUMENTS:\n\n\
     - varMat: a Numeric Array object\n\
       varMat is a Numeric array with the number of possible occurrences\n\
         of each result for reach possible value of the given variable.\n\n\
       So, for a variable which adopts 4 possible values and a result which\n\
         has 3 possible values, varMat would be 4x3\n\n\
   RETURNS:\n\n\
     - a Python float object\n\n\
   NOTES\n\n\
     - this is a dropin replacement for _PyInfoGain()_ in entropy.py\n";
  python::def("InfoGain", RDInfoTheory::infoGain, docString.c_str(),
              python::args("resArr"));

  docString =
      "Calculates the chi squared value for a variable\n\n\
   ARGUMENTS:\n\n\
     - varMat: a Numeric Array object\n\
       varMat is a Numeric array with the number of possible occurrences\n\
         of each result for reach possible value of the given variable.\n\n\
       So, for a variable which adopts 4 possible values and a result which\n\
         has 3 possible values, varMat would be 4x3\n\n\
   RETURNS:\n\n\
     - a Python float object\n";
  python::def("ChiSquare", RDInfoTheory::chiSquare, docString.c_str(),
              python::args("resArr"));
}
