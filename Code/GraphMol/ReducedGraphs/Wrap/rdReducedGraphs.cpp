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
#include <RDBoost/Wrap.h>
#include <GraphMol/GraphMol.h>
#include <numpy/arrayobject.h>
#include <boost/foreach.hpp>

#include <GraphMol/ReducedGraphs/ReducedGraphs.h>
#include <Numerics/Vector.h>

#include <vector>

namespace python = boost::python;

namespace {
  RDKit::ROMol *GenerateMolExtendedReducedGraphHelper(const RDKit::ROMol &mol,
                                               python::object atomTypes){
    RDKit::ROMol *res=RDKit::ReducedGraphs::generateMolExtendedReducedGraph(mol);
    return res;
  }
  PyArrayObject *GenerateErGFingerprintForReducedGraphHelper(const RDKit::ROMol &mol,
                                                        python::object atomTypes,
                                                        double fuzzIncrement,
                                                        int minPath,
                                                        int maxPath){
    npy_intp dims[2];
    dims[0] = 4;
    dims[1] = 4;
    PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    return res;
  }
  PyArrayObject *GetErGFingerprintHelper(const RDKit::ROMol &mol,
                                         python::object atomTypes,
                                         double fuzzIncrement,
                                         int minPath,
                                         int maxPath){
    npy_intp dims[2];
    dims[0] = 4;
    dims[1] = 4;
    PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    return res;
  }


}

BOOST_PYTHON_MODULE(rdReducedGraphs) {
  python::scope().attr("__doc__") =
    "Module containing functions to generate and work with reduced graphs"
    ;

  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);

  std::string docString = "";

  docString="Returns the reduced graph for a molecule";
  python::def("GenerateMolExtendedReducedGraph", GenerateMolExtendedReducedGraphHelper,
              (python::arg("mol"),
               python::arg("atomTypes")=0
               ),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString="Returns the ErG fingerprint vector for a reduced graph";
  python::def("GenerateErGFingerprintForReducedGraph",
	      GenerateErGFingerprintForReducedGraphHelper,
	      (python::arg("mol"),
               python::arg("atomTypes")=0,
               python::arg("fuzzIncrement")=0.3,
               python::arg("minPath")=1,
               python::arg("maxPath")=15
               ),
              docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());
  docString="Returns the ErG fingerprint vector for a molecule";
  python::def("GetErGFingerprint",
	      GetErGFingerprintHelper,
	      (python::arg("mol"),
               python::arg("atomTypes")=0,
               python::arg("fuzzIncrement")=0.3,
               python::arg("minPath")=1,
               python::arg("maxPath")=15
               ),
              docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());

  
}
