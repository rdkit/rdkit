// $Id$
//
//  Copyright (C) 2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define PY_ARRAY_UNIQUE_SYMBOL rdreducedgraphs_array_API
#include <RDBoost/python.h>
#include <RDBoost/boost_numpy.h>

#include <RDBoost/Wrap.h>
#include <GraphMol/GraphMol.h>
#include <RDBoost/import_array.h>

#include <GraphMol/ReducedGraphs/ReducedGraphs.h>
#include <Numerics/Vector.h>

#include <vector>

namespace python = boost::python;

namespace {
RDKit::ROMol *GenerateMolExtendedReducedGraphHelper(const RDKit::ROMol &mol,
                                                    python::object atomTypes) {
  if (atomTypes) {
    throw_value_error("specification of atom types not yet supported");
  }
  RDKit::ROMol *res =
      RDKit::ReducedGraphs::generateMolExtendedReducedGraph(mol);
  return res;
}
PyObject *GenerateErGFingerprintForReducedGraphHelper(const RDKit::ROMol &mol,
                                                      python::object atomTypes,
                                                      double fuzzIncrement,
                                                      int minPath,
                                                      int maxPath) {
  if (atomTypes) {
    throw_value_error("specification of atom types not yet supported");
  }
  RDNumeric::DoubleVector *dv =
      RDKit::ReducedGraphs::generateErGFingerprintForReducedGraph(
          mol, nullptr, fuzzIncrement, minPath, maxPath);
  npy_intp dim = dv->size();
  auto *res = (PyArrayObject *)PyArray_SimpleNew(1, &dim, NPY_DOUBLE);
  memcpy(static_cast<void *>(PyArray_DATA(res)),
         static_cast<void *>(dv->getData()), dv->size() * sizeof(double));
  delete dv;
  return PyArray_Return(res);
}
PyObject *GetErGFingerprintHelper(const RDKit::ROMol &mol,
                                  python::object atomTypes,
                                  double fuzzIncrement, int minPath,
                                  int maxPath) {
  if (atomTypes) {
    throw_value_error("specification of atom types not yet supported");
  }
  RDNumeric::DoubleVector *dv = RDKit::ReducedGraphs::getErGFingerprint(
      mol, nullptr, fuzzIncrement, minPath, maxPath);
  npy_intp dim = dv->size();
  auto *res = (PyArrayObject *)PyArray_SimpleNew(1, &dim, NPY_DOUBLE);
  memcpy(static_cast<void *>(PyArray_DATA(res)),
         static_cast<void *>(dv->getData()), dv->size() * sizeof(double));
  delete dv;
  return PyArray_Return(res);
}
}  // namespace

BOOST_PYTHON_MODULE(rdReducedGraphs) {
  python::scope().attr("__doc__") =
      "Module containing functions to generate and work with reduced graphs";

  rdkit_import_array();

  std::string docString = "";

  docString = "Returns the reduced graph for a molecule";
  python::def(
      "GenerateMolExtendedReducedGraph", GenerateMolExtendedReducedGraphHelper,
      (python::arg("mol"), python::arg("atomTypes") = 0), docString.c_str(),
      python::return_value_policy<python::manage_new_object>());

  docString = "Returns the ErG fingerprint vector for a reduced graph";
  python::def("GenerateErGFingerprintForReducedGraph",
              GenerateErGFingerprintForReducedGraphHelper,
              (python::arg("mol"), python::arg("atomTypes") = 0,
               python::arg("fuzzIncrement") = 0.3, python::arg("minPath") = 1,
               python::arg("maxPath") = 15),
              docString.c_str());
  docString = "Returns the ErG fingerprint vector for a molecule";
  python::def("GetErGFingerprint", GetErGFingerprintHelper,
              (python::arg("mol"), python::arg("atomTypes") = 0,
               python::arg("fuzzIncrement") = 0.3, python::arg("minPath") = 1,
               python::arg("maxPath") = 15),
              docString.c_str());
}
