//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rdeht_array_API
#include <RDBoost/python.h>
#include <boost/python/list.hpp>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/RDKitBase.h>
#include "EHTTools.h"

#include <RDBoost/boost_numpy.h>
#include <RDBoost/PySequenceHolder.h>
#include <RDBoost/Wrap.h>
#include <RDBoost/import_array.h>

namespace python = boost::python;

namespace RDKit {

namespace {

// from: https://stackoverflow.com/a/35960614
/// @brief Transfer ownership to a Python object.  If the transfer fails,
///        then object will be destroyed and an exception is thrown.
template <typename T>
boost::python::object transfer_to_python(T *t) {
  // Transfer ownership to a smart pointer, allowing for proper cleanup
  // incase Boost.Python throws.
  std::unique_ptr<T> ptr(t);

  // Use the manage_new_object generator to transfer ownership to Python.
  namespace python = boost::python;
  typename python::manage_new_object::apply<T *>::type converter;

  // Transfer ownership to the Python handler and release ownership
  // from C++.
  python::handle<> handle(converter(*ptr));
  ptr.release();

  return python::object(handle);
}

PyObject *getMatrixProp(const double *mat, unsigned int sz) {
  if (!mat) throw_value_error("matrix has not be initialized");
  npy_intp dims[2];
  dims[0] = sz;
  dims[1] = sz;

  PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);

  memcpy(PyArray_DATA(res), static_cast<const void *>(mat),
         sz * sz * sizeof(double));

  return PyArray_Return(res);
}
PyObject *getSymmMatrixProp(const double *mat, unsigned int sz) {
  if (!mat) throw_value_error("matrix has not be initialized");
  npy_intp dims[1];
  dims[0] = sz * (sz + 1) / 2;

  PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);

  memcpy(PyArray_DATA(res), static_cast<const void *>(mat),
         dims[0] * sizeof(double));

  return PyArray_Return(res);
}
PyObject *getVectorProp(const double *mat, unsigned int sz) {
  if (!mat) throw_value_error("vector has not be initialized");
  npy_intp dims[1];
  dims[0] = sz;

  PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);

  memcpy(PyArray_DATA(res), static_cast<const void *>(mat),
         sz * sizeof(double));

  return PyArray_Return(res);
}
PyObject *getChargeMatrix(EHTTools::EHTResults &self) {
  return getMatrixProp(self.reducedChargeMatrix.get(), self.numAtoms);
}

PyObject *getOPMatrix(EHTTools::EHTResults &self) {
  return getSymmMatrixProp(self.reducedOverlapPopulationMatrix.get(),
                           self.numAtoms);
}

PyObject *getCharges(EHTTools::EHTResults &self) {
  return getVectorProp(self.atomicCharges.get(), self.numAtoms);
}
python::tuple runCalc(const RDKit::ROMol &mol, int confId) {
  auto eRes = new RDKit::EHTTools::EHTResults();
  bool ok = RDKit::EHTTools::runMol(mol, *eRes, confId);
  return python::make_tuple(ok, transfer_to_python(eRes));
}

}  // end of anonymous namespace

struct EHT_wrapper {
  static void wrap() {
    std::string docString = "";

    python::class_<RDKit::EHTTools::EHTResults, boost::noncopyable>(
        "EHTResults", docString.c_str(), python::no_init)
        .def_readonly("fermiEnergy", &RDKit::EHTTools::EHTResults::fermiEnergy)
        .def_readonly("totalEnergy", &RDKit::EHTTools::EHTResults::totalEnergy)
        .def("GetReducedChargeMatrix", getChargeMatrix,
             "returns the reduced charge matrix")
        .def("GetReducedOverlapPopulationMatrix", getOPMatrix,
             "returns the reduced overlap population matrix")
        .def("GetAtomicCharges", getCharges,
             "returns the calculated atomic charges");

    docString =
        R"DOC(Runs an extended Hueckel calculation for a molecule.
The molecule should have at least one conformation

ARGUMENTS:
   - mol: molecule to use
   - confId: (optional) conformation to use

RETURNS: a 2-tuple:
   - a boolean indicating whether or not the calculation succeeded
   - an EHTResults object with the results
)DOC";
    python::def("RunMol", runCalc,
                (python::arg("mol"), python::arg("confId") = -1),
                docString.c_str());
  }
};

}  // end of namespace RDKit
BOOST_PYTHON_MODULE(rdEHTTools) {
  python::scope().attr("__doc__") =
      R"DOC(Module containing interface to the YAeHMOP extended Hueckel library.
Please note that this interface should still be considered experimental and may
change from one release to the next.)DOC";
  rdkit_import_array();

  RDKit::EHT_wrapper::wrap();
}