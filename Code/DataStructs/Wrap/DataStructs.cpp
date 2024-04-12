// $Id$
//
// Copyright (c) 2001-2006 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rddatastructs_array_API

#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>
#include <DataStructs/BitVects.h>
#include <DataStructs/DiscreteValueVect.h>
#include <DataStructs/SparseIntVect.h>

#include <RDBoost/boost_numpy.h>
#include <numpy/npy_common.h>
#include <RDBoost/import_array.h>
#include <RDBoost/pyint_api.h>

namespace python = boost::python;

void wrap_SBV();
void wrap_EBV();
void wrap_BitOps();
void wrap_Utils();
void wrap_discreteValVect();
void wrap_sparseIntVect();
void wrap_FPB();

template <typename T>
void convertToNumpyArray(const T &v, python::object destArray) {
  if (!PyArray_Check(destArray.ptr())) {
    throw_value_error("Expecting a Numeric array object");
  }
  auto *destP = (PyArrayObject *)destArray.ptr();
  npy_intp ndims[1];
  ndims[0] = v.size();
  PyArray_Dims dims;
  dims.ptr = ndims;
  dims.len = 1;
  PyArray_Resize(destP, &dims, 0, NPY_ANYORDER);
  for (unsigned int i = 0; i < v.size(); ++i) {
    PyObject *iItem = PyInt_FromLong(v[i]);
    PyArray_SETITEM(destP, static_cast<char *>(PyArray_GETPTR1(destP, i)),
                    iItem);
    Py_DECREF(iItem);
  }
}

BOOST_PYTHON_MODULE(cDataStructs) {
  rdkit_import_array();
  python::scope().attr("__doc__") =
      "Module containing an assortment of functionality for basic data "
      "structures.\n"
      "\n"
      "At the moment the data structures defined are:\n"
      "  Bit Vector classes (for storing signatures, fingerprints and the "
      "like:\n"
      "    - ExplicitBitVect: class for relatively small (10s of thousands of "
      "bits) or\n"
      "                       dense bit vectors.\n"
      "    - SparseBitVect:   class for large, sparse bit vectors\n"
      "  DiscreteValueVect:   class for storing vectors of integers\n"
      "  SparseIntVect:       class for storing sparse vectors of integers\n";

  wrap_Utils();
  wrap_SBV();
  wrap_EBV();
  wrap_BitOps();
  wrap_discreteValVect();
  wrap_sparseIntVect();
  wrap_FPB();

  python::def(
      "ConvertToNumpyArray",
      (void (*)(const ExplicitBitVect &, python::object))convertToNumpyArray,
      (python::arg("bv"), python::arg("destArray")));
  python::def("ConvertToNumpyArray",
              (void (*)(const RDKit::DiscreteValueVect &,
                        python::object))convertToNumpyArray,
              (python::arg("bv"), python::arg("destArray")));
  python::def("ConvertToNumpyArray",
              (void (*)(const RDKit::SparseIntVect<std::int32_t> &,
                        python::object))convertToNumpyArray,
              (python::arg("bv"), python::arg("destArray")));
  python::def("ConvertToNumpyArray",
              (void (*)(const RDKit::SparseIntVect<boost::int64_t> &,
                        python::object))convertToNumpyArray,
              (python::arg("bv"), python::arg("destArray")));
  python::def("ConvertToNumpyArray",
              (void (*)(const RDKit::SparseIntVect<std::uint32_t> &,
                        python::object))convertToNumpyArray,
              (python::arg("bv"), python::arg("destArray")));
  python::def("ConvertToNumpyArray",
              (void (*)(const RDKit::SparseIntVect<boost::uint64_t> &,
                        python::object))convertToNumpyArray,
              (python::arg("bv"), python::arg("destArray")));
}
