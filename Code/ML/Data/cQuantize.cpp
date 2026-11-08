//
// Copyright 2003-2008 Rational Discovery LLC and Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <cstring>

#include <RDBoost/Wrap.h>
#include <RDBoost/import_array.h>

namespace python = boost::python;
#include "Quantize.h"

/***********************************************

   Recursively finds the best quantization boundaries

   **Arguments**

     - vals: a 1D Numeric array with the values of the variables,
       this should be sorted

     - cuts: a list with the indices of the quantization bounds
       (indices are into _starts_ )

     - which: an integer indicating which bound is being adjusted here
       (and index into _cuts_ )

     - starts: a list of potential starting points for quantization bounds

     - results: a 1D Numeric array of integer result codes

     - nPossibleRes: an integer with the number of possible result codes

   **Returns**

     - a 2-tuple containing:

       1) the best information gain found so far

       2) a list of the quantization bound indices ( _cuts_ for the best case)

   **Notes**

    - this is not even remotely efficient, which is why a C replacement
      was written

    - this is a drop-in replacement for *ML.Data.Quantize._PyRecurseBounds*

 ***********************************************/
static python::tuple cQuantize_RecurseOnBounds(python::object vals,
                                               python::list pyCuts, int which,
                                               python::list pyStarts,
                                               python::object results,
                                               int nPossibleRes) {
  PyArrayObject *contigVals, *contigResults;
  long int *cuts, *starts;

  /*
    -------

    Setup code

    -------
  */
  contigVals = reinterpret_cast<PyArrayObject *>(
      PyArray_ContiguousFromObject(vals.ptr(), NPY_DOUBLE, 1, 1));
  if (!contigVals) {
    throw_value_error("could not convert value argument");
  }

  contigResults = reinterpret_cast<PyArrayObject *>(
      PyArray_ContiguousFromObject(results.ptr(), NPY_LONG, 1, 1));
  if (!contigResults) {
    throw_value_error("could not convert results argument");
  }

  python::ssize_t nCuts = python::len(pyCuts);
  cuts = (long int *)calloc(nCuts, sizeof(long int));
  CHECK_INVARIANT(cuts, "failed to allocate memory");
  for (python::ssize_t i = 0; i < nCuts; i++) {
    python::object elem = pyCuts[i];
    cuts[i] = python::extract<long int>(elem);
  }

  python::ssize_t nStarts = python::len(pyStarts);
  starts = (long int *)calloc(nStarts, sizeof(long int));
  CHECK_INVARIANT(starts, "failed to allocate memory");
  for (python::ssize_t i = 0; i < nStarts; i++) {
    python::object elem = pyStarts[i];
    starts[i] = python::extract<long int>(elem);
  }

  // do the real work
  double gain = RecurseHelper(
      (double *)PyArray_DATA(contigVals), PyArray_DIM(contigVals, 0), cuts,
      nCuts, which, starts, nStarts, (long int *)PyArray_DATA(contigResults),
      nPossibleRes);

  /*
    -------

    Construct the return value

    -------
  */
  python::list cutObj;
  for (python::ssize_t i = 0; i < nCuts; i++) {
    cutObj.append(cuts[i]);
  }
  free(cuts);
  free(starts);
  return python::make_tuple(gain, cutObj);
}

static python::list cQuantize_FindStartPoints(python::object values,
                                              python::object results,
                                              int nData) {
  python::list startPts;

  if (nData < 2) {
    return startPts;
  }

  auto *contigVals = reinterpret_cast<PyArrayObject *>(
      PyArray_ContiguousFromObject(values.ptr(), NPY_DOUBLE, 1, 1));
  if (!contigVals) {
    throw_value_error("could not convert value argument");
  }

  auto *vals = (double *)PyArray_DATA(contigVals);

  auto *contigResults = reinterpret_cast<PyArrayObject *>(
      PyArray_ContiguousFromObject(results.ptr(), NPY_LONG, 1, 1));
  if (!contigResults) {
    throw_value_error("could not convert results argument");
  }

  long *res = (long *)PyArray_DATA(contigResults);

  bool firstBlock = true;
  long lastBlockAct = -2, blockAct = res[0];
  int lastDiv = -1;
  double tol = 1e-8;

  int i = 1;
  while (i < nData) {
    while (i < nData && vals[i] - vals[i - 1] <= tol) {
      if (res[i] != blockAct) {
        blockAct = -1;
      }
      ++i;
    }
    if (firstBlock) {
      firstBlock = false;
      lastBlockAct = blockAct;
      lastDiv = i;
    } else {
      if (blockAct == -1 || lastBlockAct == -1 || blockAct != lastBlockAct) {
        startPts.append(lastDiv);
        lastDiv = i;
        lastBlockAct = blockAct;
      } else {
        lastDiv = i;
      }
    }
    if (i < nData) {
      blockAct = res[i];
    }
    ++i;
  }

  // catch the case that the last point also sets a bin:
  if (blockAct != lastBlockAct) {
    startPts.append(lastDiv);
  }

  return startPts;
}

BOOST_PYTHON_MODULE(cQuantize) {
  rdkit_import_array();

  python::def("_RecurseOnBounds", cQuantize_RecurseOnBounds,
              (python::arg("vals"), python::arg("pyCuts"), python::arg("which"),
               python::arg("pyStarts"), python::arg("results"),
               python::arg("nPossibleRes")),
              "TODO: provide docstring");
  python::def(
      "_FindStartPoints", cQuantize_FindStartPoints,
      (python::arg("values"), python::arg("results"), python::arg("nData")),
      "TODO: provide docstring");
}
