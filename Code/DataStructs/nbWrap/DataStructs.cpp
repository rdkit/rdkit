//
// Copyright (C) 2026 greg Landrum
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>

#include <DataStructs/DiscreteValueVect.h>
#include <DataStructs/RealValueVect.h>
#include <DataStructs/BitVects.h>
#include <DataStructs/SparseIntVect.h>
#include <numpy/npy_common.h>
#include <RDBoost/import_array.h>
// #include <RDBoost/pyint_api.h>

namespace nb = nanobind;
using namespace nb::literals;

void wrap_discreteValVect(nb::module_ &m);
void wrap_SBV(nb::module_ &m);
void wrap_EBV(nb::module_ &m);
void wrap_BitOps(nb::module_ &m);
void wrap_Utils(nb::module_ &m);
void wrap_realValVect(nb::module_ &m);
void wrap_sparseIntVect(nb::module_ &m);
void wrap_FPB(nb::module_ &m);

#if 1
namespace {
template <typename T, typename U>
void converter(const T &v, nb::object destArray, U func) {
  if (!PyArray_Check(destArray.ptr())) {
    throw nb::value_error("Expecting a Numeric array object");
  }
  auto *destP = (PyArrayObject *)destArray.ptr();
  npy_intp ndims[1];
  ndims[0] = v.size();
  PyArray_Dims dims;
  dims.ptr = ndims;
  dims.len = 1;
  PyArray_Resize(destP, &dims, 0, NPY_ANYORDER);
  for (unsigned int i = 0; i < v.size(); ++i) {
    PyObject *iItem = func(v[i]);
    PyArray_SETITEM(destP, static_cast<char *>(PyArray_GETPTR1(destP, i)),
                    iItem);
    Py_DECREF(iItem);
  }
}
}  // namespace

template <typename T>
void convertToIntNumpyArray(const T &v, nb::object destArray) {
  converter(v, destArray, PyLong_FromLong);
}
template <typename T>
void convertToDoubleNumpyArray(const T &v, nb::object destArray) {
  converter(v, destArray, PyFloat_FromDouble);
}
#endif
NB_MODULE(cDataStructs, m) {
  rdkit_import_array();
  m.doc() =
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

  wrap_Utils(m);
  wrap_SBV(m);
  wrap_EBV(m);
  wrap_BitOps(m);
  wrap_discreteValVect(m);
  wrap_realValVect(m);
  wrap_sparseIntVect(m);
  wrap_FPB(m);
  m.def("ConvertToNumpyArray",
        (void (*)(const RDKit::DiscreteValueVect &,
                  nb::object))convertToIntNumpyArray,
        "bv"_a, "destArray"_a);
  m.def("ConvertToNumpyArray",
        (void (*)(const RDKit::RealValueVect &,
                  nb::object))convertToDoubleNumpyArray,
        "rvv"_a, "destArray"_a);
  m.def("ConvertToNumpyArray",
        (void (*)(const ExplicitBitVect &, nb::object))convertToIntNumpyArray,
        "bv"_a, "destArray"_a);
  m.def("ConvertToNumpyArray",
        (void (*)(const RDKit::SparseIntVect<std::int32_t> &,
                  nb::object))convertToIntNumpyArray,
        "bv"_a, "destArray"_a);
  m.def("ConvertToNumpyArray",
        (void (*)(const RDKit::SparseIntVect<boost::int64_t> &,
                  nb::object))convertToIntNumpyArray,
        "bv"_a, "destArray"_a);
  m.def("ConvertToNumpyArray",
        (void (*)(const RDKit::SparseIntVect<std::uint32_t> &,
                  nb::object))convertToIntNumpyArray,
        "bv"_a, "destArray"_a);
  m.def("ConvertToNumpyArray",
        (void (*)(const RDKit::SparseIntVect<boost::uint64_t> &,
                  nb::object))convertToIntNumpyArray,
        "bv"_a, "destArray"_a);
}
