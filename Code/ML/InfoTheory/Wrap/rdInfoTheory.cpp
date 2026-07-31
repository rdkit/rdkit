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
#include <RDBoost/import_array.h>
#include <ML/InfoTheory/InfoBitRanker.h>
#include <ML/InfoTheory/InfoGainFuncs.h>

namespace python = boost::python;
using namespace RDInfoTheory;

namespace RDInfoTheory {
namespace {
// Dispatching on the input's own type_num used to leave uncommon dtypes
// (NPY_LONGLONG, i.e. numpy's int64 on Windows, but also int8/uint*/float16)
// unhandled, which silently produced a result of 0.0 (github #8421).
PyArrayObject *contiguousDoubles(const python::object &resArr, int minDim,
                                 int maxDim) {
  if (!PyArray_Check(resArr.ptr())) {
    throw_value_error("Expecting a Numeric array object");
  }
  auto *copy = (PyArrayObject *)PyArray_ContiguousFromObject(
      resArr.ptr(), NPY_DOUBLE, minDim, maxDim);
  if (!copy) {
    throw_value_error("could not convert argument to an array of doubles");
  }
  return copy;
}
}  // namespace

double infoEntropy(python::object resArr) {
  PyArrayObject *copy = contiguousDoubles(resArr, 1, 1);
  auto ncols = (long int)PyArray_DIM(copy, 0);
  CHECK_INVARIANT(ncols > 0, "");
  double res = InfoEntropy((double *)PyArray_DATA(copy), ncols);
  Py_DECREF(copy);
  return res;
}

double infoGain(python::object resArr) {
  PyArrayObject *copy = contiguousDoubles(resArr, 2, 2);
  auto rows = (long int)PyArray_DIM(copy, 0);
  auto cols = (long int)PyArray_DIM(copy, 1);
  double res = InfoEntropyGain((double *)PyArray_DATA(copy), rows, cols);
  Py_DECREF(copy);
  return res;
}

double chiSquare(python::object resArr) {
  PyArrayObject *copy = contiguousDoubles(resArr, 2, 2);
  auto rows = (long int)PyArray_DIM(copy, 0);
  auto cols = (long int)PyArray_DIM(copy, 1);
  double res = ChiSquare((double *)PyArray_DATA(copy), rows, cols);
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
