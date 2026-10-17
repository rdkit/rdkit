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

namespace {
template <typename T>
void convertToNumpyArray(const T &v, nb::object destArray) {
  if (!PyArray_Check(destArray.ptr())) {
    throw nb::value_error("Expecting a Numeric array object");
  }
  auto *destP = reinterpret_cast<PyArrayObject *>(destArray.ptr());
  npy_intp size = static_cast<npy_intp>(v.size());
  PyArray_Dims dims{&size, 1};
  PyObject *resized = PyArray_Resize(destP, &dims, 0, NPY_ANYORDER);
  if (!resized) {
    throw nb::python_error();
  }
  Py_DECREF(resized);
  for (unsigned int i = 0; i < v.size(); ++i) {
    nb::object item = nb::cast(+v[i]);
    if (PyArray_SETITEM(destP,
                       static_cast<char *>(PyArray_GETPTR1(destP, i)),
                       item.ptr()) < 0) {
      throw nb::python_error();
    }
  }
}
}  // namespace

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
        &convertToNumpyArray<RDKit::DiscreteValueVect>, "bv"_a, "destArray"_a);
  m.def("ConvertToNumpyArray", &convertToNumpyArray<RDKit::RealValueVect>,
        "rvv"_a, "destArray"_a);
  m.def("ConvertToNumpyArray", &convertToNumpyArray<ExplicitBitVect>,
        "bv"_a, "destArray"_a);
  m.def("ConvertToNumpyArray",
        &convertToNumpyArray<RDKit::SparseIntVect<std::int32_t>>,
        "bv"_a, "destArray"_a);
  m.def("ConvertToNumpyArray",
        &convertToNumpyArray<RDKit::SparseIntVect<boost::int64_t>>,
        "bv"_a, "destArray"_a);
  m.def("ConvertToNumpyArray",
        &convertToNumpyArray<RDKit::SparseIntVect<std::uint32_t>>,
        "bv"_a, "destArray"_a);
  m.def("ConvertToNumpyArray",
        &convertToNumpyArray<RDKit::SparseIntVect<boost::uint64_t>>,
        "bv"_a, "destArray"_a);
}
