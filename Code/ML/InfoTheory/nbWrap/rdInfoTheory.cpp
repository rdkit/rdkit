//
//  Copyright (C) 2003-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <ML/InfoTheory/InfoBitRanker.h>
#include <ML/InfoTheory/InfoGainFuncs.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDInfoTheory {

template <class T>
double infoEntropyHelper(nb::ndarray<nb::numpy, nb::ndim<1>> arr) {
  auto ncols = (long int)arr.shape(0);
  CHECK_INVARIANT(ncols > 0, "");
  auto *data = reinterpret_cast<T *>(arr.data());
  return InfoEntropy(data, ncols);
}

double infoEntropy(nb::ndarray<nb::numpy, nb::ndim<1>> resArr) {
  if (resArr.dtype() == nb::dtype<double>()) {
    return infoEntropyHelper<double>(resArr);
  } else if (resArr.dtype() == nb::dtype<float>()) {
    return infoEntropyHelper<float>(resArr);
  } else if (resArr.dtype() == nb::dtype<int>()) {
    return infoEntropyHelper<int>(resArr);
  } else if (resArr.dtype() == nb::dtype<long int>()) {
    return infoEntropyHelper<long int>(resArr);
  } else {
    throw nb::value_error(
        "Expecting a Numeric array object of type int, long, float, or double");
  }
}

template <class T>
double infoGainHelper(nb::ndarray<nb::numpy, nb::ndim<2>> arr) {
  auto rows = arr.shape(0);
  auto cols = arr.shape(1);
  auto *data = reinterpret_cast<T *>(arr.data());
  return InfoEntropyGain(data, rows, cols);
}

double infoGain(nb::ndarray<nb::numpy, nb::ndim<2>> resArr) {
  if (resArr.dtype() == nb::dtype<double>()) {
    return infoGainHelper<double>(resArr);
  } else if (resArr.dtype() == nb::dtype<float>()) {
    return infoGainHelper<float>(resArr);
  } else if (resArr.dtype() == nb::dtype<int>()) {
    return infoGainHelper<int>(resArr);
  } else if (resArr.dtype() == nb::dtype<long int>()) {
    return infoGainHelper<long int>(resArr);
  } else {
    throw nb::value_error(
        "Numeric array object of type int or long or float or double");
  }
}

template <class T>
double chiSquareHelper(nb::ndarray<nb::numpy, nb::ndim<2>> arr) {
  auto rows = arr.shape(0);
  auto cols = arr.shape(1);
  auto *data = reinterpret_cast<T *>(arr.data());
  return ChiSquare(data, rows, cols);
}

double chiSquare(nb::ndarray<nb::numpy, nb::ndim<2>> resArr) {
  if (resArr.dtype() == nb::dtype<double>()) {
    return chiSquareHelper<double>(resArr);
  } else if (resArr.dtype() == nb::dtype<float>()) {
    return chiSquareHelper<float>(resArr);
  } else if (resArr.dtype() == nb::dtype<int>()) {
    return chiSquareHelper<int>(resArr);
  } else if (resArr.dtype() == nb::dtype<long int>()) {
    return chiSquareHelper<long int>(resArr);
  } else {
    throw nb::value_error(
        "Numeric array object of type int or long or float or double");
  }
}

}  // namespace RDInfoTheory

void wrap_ranker(nb::module_ &m);
void wrap_corrmatgen(nb::module_ &m);

NB_MODULE(rdInfoTheory, m) {
  m.doc() =
      "Module containing bunch of functions for information metrics and a "
      "ranker to rank bits";

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
