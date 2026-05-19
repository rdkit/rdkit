//
//  Copyright (C) 2003-2026 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rdmetric_array_API
#include <nanobind/nanobind.h>

#include <numpy/ndarraytypes.h>
#include <numpy/npy_common.h>
#include <RDBoost/import_array.h>

#include <RDGeneral/types.h>

#include <DataManip/MetricMatrixCalc/MetricMatrixCalc.h>
#include <DataManip/MetricMatrixCalc/MetricFuncs.h>
#include <DataStructs/BitVects.h>

using namespace RDDataManip;

namespace nb = nanobind;
using namespace nb::literals;

namespace RDDataManip {

//! Minimal sequence holder backed by a nanobind object.
//! Provides size() and operator[] that lazily extract elements as type T.
template <typename T>
class NbSequenceHolder {
 public:
  NbSequenceHolder(nb::object seq) : d_seq(std::move(seq)) {}

  unsigned int size() const {
    return static_cast<unsigned int>(nb::len(d_seq));
  }

  T operator[](unsigned int which) const {
    return nb::cast<T>(d_seq[which]);
  }

 private:
  nb::object d_seq;
};

nb::object getEuclideanDistMat(nb::object descripMat) {
  // Bit of a pain involved here, we accept three types of PyObjects here
  // 1. A Numeric Array
  //     - first find what 'type' of entry we have (float, double and int is all
  //     we recognize for now)
  //     - then point to contiguous piece of memory from the array that contains
  //     the data with a type*
  //     - then make a new type** pointer so that double index into this
  //     contiguous memory will work
  //       and then pass it along to the distance calculator
  // 2. A list of Numeric Vector (or 1D arrays)
  //     - in this case wrap descripMat with a NbSequenceHolder<type*> where
  //     type is the
  //       type of entry in vector (accepted types are int, double and float
  //     - Then pass the NbSequenceHolder to the metric calculator
  // 3. A list (or tuple) of lists (or tuple)
  //     - In this case other than wrapping descripMat with a NbSequenceHolder
  //       each of the individual list in there are also wrapped by a
  //       NbSequenceHolder
  //     - so the distance calculator is passed in a
  //     "NbSequenceHolder<NbSequenceHolder<double>>"
  //     - FIX: not that we always convert entry values to double here, even if
  //     we passed
  //       in a list of list of ints (or floats). Given that lists can be
  //       heterogeneous, I do not
  //       know how to ask a list what type of entries if contains.
  //
  //  OK my brain is going to explode now

  // first deal with situation where we have a Numeric Array
  PyObject *descMatObj = descripMat.ptr();
  PyArrayObject *distResRaw = nullptr;
  if (PyArray_Check(descMatObj)) {
    // get the dimensions of the array
    int nrows = PyArray_DIM((PyArrayObject *)descMatObj, 0);
    int ncols = PyArray_DIM((PyArrayObject *)descMatObj, 1);
    int i;
    CHECK_INVARIANT((nrows > 0) && (ncols > 0), "");

    npy_intp dMatLen = nrows * (nrows - 1) / 2;

    // now that we have the dimensions declare the distance matrix which is
    // always a 1D double array
    distResRaw = (PyArrayObject *)PyArray_SimpleNew(1, &dMatLen, NPY_DOUBLE);

    // grab a pointer to the data in the array so that we can directly put
    // values in there and avoid copying
    auto *dMat = (double *)PyArray_DATA(distResRaw);

    PyArrayObject *copyRaw = (PyArrayObject *)PyArray_ContiguousFromObject(
        descMatObj, PyArray_DESCR((PyArrayObject *)descMatObj)->type_num, 2, 2);

    // if we have a double array
    if (PyArray_DESCR((PyArrayObject *)descMatObj)->type_num == NPY_DOUBLE) {
      auto *desc = (double *)PyArray_DATA((PyArrayObject *)descMatObj);

      // here is the 2D array trick so that when the distance calculator
      // asks for desc2D[i] we basically get the ith row as double*
      std::unique_ptr<double *[]> desc2D(new double *[nrows]);
      for (i = 0; i < nrows; i++) {
        desc2D[i] = desc;
        desc += ncols;
      }
      MetricMatrixCalc<double **, double *> mmCalc;
      mmCalc.setMetricFunc(&EuclideanDistanceMetric<double *, double *>);
      mmCalc.calcMetricMatrix(desc2D.get(), nrows, ncols, dMat);

      Py_XDECREF((PyObject *)copyRaw);
      // we got the distance matrix, we are happy so return
      return nb::steal<nb::object>(PyArray_Return(distResRaw));
    }

    // if we have a float array
    else if (PyArray_DESCR((PyArrayObject *)descMatObj)->type_num == NPY_FLOAT) {
      auto *desc = (float *)PyArray_DATA(copyRaw);
      std::unique_ptr<float *[]> desc2D(new float *[nrows]);
      for (i = 0; i < nrows; i++) {
        desc2D[i] = desc;
        desc += ncols;
      }
      MetricMatrixCalc<float **, float *> mmCalc;
      mmCalc.setMetricFunc(&EuclideanDistanceMetric<float *, float *>);
      mmCalc.calcMetricMatrix(desc2D.get(), nrows, ncols, dMat);
      Py_XDECREF((PyObject *)copyRaw);
      return nb::steal<nb::object>(PyArray_Return(distResRaw));
    }

    // if we have an integer array
    else if (PyArray_DESCR((PyArrayObject *)descMatObj)->type_num == NPY_INT) {
      int *desc = (int *)PyArray_DATA(copyRaw);
      std::unique_ptr<int *[]> desc2D(new int *[nrows]);
      for (i = 0; i < nrows; i++) {
        desc2D[i] = desc;
        desc += ncols;
      }
      MetricMatrixCalc<int **, int *> mmCalc;
      mmCalc.setMetricFunc(&EuclideanDistanceMetric<int *, int *>);
      mmCalc.calcMetricMatrix(desc2D.get(), nrows, ncols, dMat);
      Py_XDECREF((PyObject *)copyRaw);
      return nb::steal<nb::object>(PyArray_Return(distResRaw));
    } else {
      Py_XDECREF((PyObject *)copyRaw);
      Py_XDECREF((PyObject *)distResRaw);
      // unrecognized type for the matrix, throw up
      throw nb::value_error(
          "The array has to be of type int, float, or double for "
          "GetEuclideanDistMat");
    }
  }  // done with an array input
  else {
    // we probably have a list or a tuple

    unsigned int ncols = 0;
    unsigned int nrows = static_cast<unsigned int>(nb::len(descripMat));
    CHECK_INVARIANT(nrows > 0, "Empty list passed in");

    npy_intp dMatLen = nrows * (nrows - 1) / 2;
    distResRaw = (PyArrayObject *)PyArray_SimpleNew(1, &dMatLen, NPY_DOUBLE);
    auto *dMat = (double *)PyArray_DATA(distResRaw);

    // assume that we have a list of list of values (that can be extracted to
    // double)
    std::vector<NbSequenceHolder<double>> dData;
    dData.reserve(nrows);
    for (unsigned int i = 0; i < nrows; i++) {
      NbSequenceHolder<double> row(descripMat[i]);
      if (i == 0) {
        ncols = row.size();
      } else if (row.size() != ncols) {
        Py_XDECREF((PyObject *)distResRaw);
        throw nb::value_error("All subsequences must be the same length");
      }
      dData.push_back(std::move(row));
    }

    MetricMatrixCalc<std::vector<NbSequenceHolder<double>>,
                     NbSequenceHolder<double>>
        mmCalc;
    mmCalc.setMetricFunc(&EuclideanDistanceMetric<NbSequenceHolder<double>,
                                                  NbSequenceHolder<double>>);
    mmCalc.calcMetricMatrix(dData, nrows, ncols, dMat);
  }
  return nb::steal<nb::object>(PyArray_Return(distResRaw));
}

nb::object getTanimotoDistMat(nb::object bitVectList) {
  // we will assume here that we have a either a list of ExplicitBitVectors or
  // SparseBitVects
  unsigned int nrows = static_cast<unsigned int>(nb::len(bitVectList));
  CHECK_INVARIANT(nrows > 1, "");

  // First check what type of vector we have
  nb::object v1 = bitVectList[0];
  bool ebvWorks = nb::isinstance<ExplicitBitVect>(v1);
  bool sbvWorks = !ebvWorks && nb::isinstance<SparseBitVect>(v1);
  if (!ebvWorks && !sbvWorks) {
    throw nb::value_error(
        "GetTanimotoDistMat can only take a sequence of ExplicitBitVects or "
        "SparseBitvects");
  }

  npy_intp dMatLen = nrows * (nrows - 1) / 2;
  auto *simRes = (PyArrayObject *)PyArray_SimpleNew(1, &dMatLen, NPY_DOUBLE);
  auto *sMat = (double *)PyArray_DATA(simRes);

  if (ebvWorks) {
    NbSequenceHolder<ExplicitBitVect> dData(bitVectList);
    MetricMatrixCalc<NbSequenceHolder<ExplicitBitVect>, ExplicitBitVect> mmCalc;
    mmCalc.setMetricFunc(
        &TanimotoDistanceMetric<ExplicitBitVect, ExplicitBitVect>);
    mmCalc.calcMetricMatrix(dData, nrows, 0, sMat);
  } else {
    NbSequenceHolder<SparseBitVect> dData(bitVectList);
    MetricMatrixCalc<NbSequenceHolder<SparseBitVect>, SparseBitVect> mmCalc;
    mmCalc.setMetricFunc(&TanimotoDistanceMetric<SparseBitVect, SparseBitVect>);
    mmCalc.calcMetricMatrix(dData, nrows, 0, sMat);
  }
  return nb::steal<nb::object>(PyArray_Return(simRes));
}

nb::object getTanimotoSimMat(nb::object bitVectList) {
  // we will assume here that we have a either a list of ExplicitBitVectors or
  // SparseBitVects
  unsigned int nrows = static_cast<unsigned int>(nb::len(bitVectList));
  CHECK_INVARIANT(nrows > 1, "");

  // First check what type of vector we have
  nb::object v1 = bitVectList[0];
  bool ebvWorks = nb::isinstance<ExplicitBitVect>(v1);
  bool sbvWorks = !ebvWorks && nb::isinstance<SparseBitVect>(v1);
  if (!ebvWorks && !sbvWorks) {
    throw nb::value_error(
        "GetTanimotoSimMat can only take a sequence of ExplicitBitVects or "
        "SparseBitvects");
  }

  npy_intp dMatLen = nrows * (nrows - 1) / 2;
  auto *simRes = (PyArrayObject *)PyArray_SimpleNew(1, &dMatLen, NPY_DOUBLE);
  auto *sMat = (double *)PyArray_DATA(simRes);

  if (ebvWorks) {
    NbSequenceHolder<ExplicitBitVect> dData(bitVectList);
    MetricMatrixCalc<NbSequenceHolder<ExplicitBitVect>, ExplicitBitVect> mmCalc;
    mmCalc.setMetricFunc(
        &TanimotoSimilarityMetric<ExplicitBitVect, ExplicitBitVect>);
    mmCalc.calcMetricMatrix(dData, nrows, 0, sMat);
  } else {
    NbSequenceHolder<SparseBitVect> dData(bitVectList);
    MetricMatrixCalc<NbSequenceHolder<SparseBitVect>, SparseBitVect> mmCalc;
    mmCalc.setMetricFunc(
        &TanimotoSimilarityMetric<SparseBitVect, SparseBitVect>);
    mmCalc.calcMetricMatrix(dData, nrows, 0, sMat);
  }
  return nb::steal<nb::object>(PyArray_Return(simRes));
}
}  // namespace RDDataManip

NB_MODULE(rdMetricMatrixCalc, m) {
  rdkit_import_array();

  m.doc() = R"DOC(Module containing the calculator for metric matrix calculation,
e.g. similarity and distance matrices)DOC";

  m.def("GetEuclideanDistMat", RDDataManip::getEuclideanDistMat,
        R"DOC(Compute the distance matrix from a descriptor matrix using the Euclidean distance metric

ARGUMENTS:

  descripMat - A python object of any one of the following types
               1. A numeric array of dimensions n by m where n is the number of items in the data set
                  and m is the number of descriptors
               2. A list of Numeric Vectors (or 1D arrays), each entry in the list corresponds
                  to descriptor vector for one item
               3. A list (or tuple) of lists (or tuples) of values, where the values can be extracted to
                  double.

RETURNS:
  A numeric one-dimensional array containing the lower triangle elements of the symmetric distance matrix
)DOC",
        "descripMat"_a);

  m.def("GetTanimotoDistMat", RDDataManip::getTanimotoDistMat,
        R"DOC(Compute the distance matrix from a list of BitVects using the Tanimoto distance metric

ARGUMENTS:

  bitVectList - a list of bit vectors. Currently this works only for a list of explicit bit vectors,
                needs to be expanded to support a list of SparseBitVects

RETURNS:
  A numeric 1 dimensional array containing the lower triangle elements of the
  symmetric distance matrix
)DOC",
        "bitVectList"_a);

  m.def("GetTanimotoSimMat", RDDataManip::getTanimotoSimMat,
        R"DOC(Compute the similarity matrix from a list of BitVects

ARGUMENTS:

  bitVectList - a list of bit vectors. Currently this works only for a list of explicit bit vectors,
                needs to be expanded to support a list of SparseBitVects

RETURNS:
  A numeric 1 dimensional array containing the lower triangle elements of the symmetric similarity matrix
)DOC",
        "bitVectList"_a);
}
