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
#include <YAeHMOP/EHTTools.h>

#include <RDBoost/boost_numpy.h>
#include <RDBoost/PySequenceHolder.h>
#include <RDBoost/Wrap.h>
#include <RDBoost/import_array.h>

namespace python = boost::python;

namespace RDKit {

namespace {

PyObject *getMatrixProp(ROMol &mol, std::string propName) {
  if (!mol.hasProp(propName)) {
    throw ValueErrorException("eHT calculation not yet run");
  }

  int nats = mol.getNumAtoms();
  npy_intp dims[2];
  dims[0] = nats;
  dims[1] = nats;
  double *mat = mol.getProp<boost::shared_array<double>>(propName).get();
  CHECK_INVARIANT(mat, "invalid matrix");

  PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);

  memcpy(PyArray_DATA(res), static_cast<void *>(mat),
         nats * nats * sizeof(double));

  return PyArray_Return(res);
}
PyObject *getVectorProp(ROMol &mol, std::string propName) {}
PyObject *getChargeMatrix(ROMol &mol) {
  return getMatrixProp(mol, EHTTools::_EHTChargeMatrix);
}

PyObject *getOPMatrix(ROMol &mol) {
  std::string propName = EHTTools::_EHTMullikenOP;
  if (!mol.hasProp(propName)) {
    throw ValueErrorException("eHT calculation not yet run");
  }

  int nats = mol.getNumAtoms();
  npy_intp dims[1];
  dims[0] = nats * (nats + 1) / 2;
  double *mat = mol.getProp<boost::shared_array<double>>(propName).get();
  CHECK_INVARIANT(mat, "invalid matrix");

  PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);

  memcpy(PyArray_DATA(res), static_cast<void *>(mat), dims[0] * sizeof(double));

  return PyArray_Return(res);
}

}  // end of anonymous namespace

struct EHT_wrapper {
  static void wrap() {
    std::string docString = "";

    docString =
        "Runs the extended Hueckel calculation.\n"
        "ARGUMENTS:\n"
        "   - mol: molecule to use\n"
        "   - confId: (optional) conformation to use\n"
        "\n";
    python::def("RunMol", RDKit::EHTTools::runMol,
                (python::arg("mol"), python::arg("confId") = -1),
                docString.c_str());
    docString =
        "Returns the extended Hueckel charge matrix for a molecule.\n"
        "Note that you have to call RunMol() first."
        "ARGUMENTS:\n"
        "   - mol: molecule to use\n"
        "\n";
    python::def("GetChargeMatrix", getChargeMatrix, (python::arg("mol")),
                docString.c_str());
    docString =
        "Returns the extended Hueckel overlap population matrix for a "
        "molecule.\n"
        "Note that you have to call RunMol() first."
        "ARGUMENTS:\n"
        "   - mol: molecule to use\n"
        "\n";
    python::def("GetOverlapPopulationMatrix", getOPMatrix, (python::arg("mol")),
                docString.c_str());
  }
};

}  // end of namespace RDKit
BOOST_PYTHON_MODULE(rdEHTTools) {
  python::scope().attr("__doc__") =
      "Module containing interface to the YAeHMOP extended Hueckel library.";
  rdkit_import_array();

  RDKit::EHT_wrapper::wrap();
}