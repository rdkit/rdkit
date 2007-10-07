// $Id$
//
//  Copyright (C) 2007 Greg Landrum
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/Atom.h>
#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap.h>

#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/AtomPairs.h>

#include <vector>

namespace python = boost::python;

namespace {
  std::vector<unsigned int> atomPairTypes(RDKit::Descriptors::AtomPairs::atomNumberTypes,
					  RDKit::Descriptors::AtomPairs::atomNumberTypes+sizeof(RDKit::Descriptors::AtomPairs::atomNumberTypes)/sizeof(unsigned int));
}

BOOST_PYTHON_MODULE(rdMolDescriptors) {
  python::scope().attr("__doc__") =
    "Module containing functions to compute molecular descriptors"
    ;
  std::string docString = "";

  python::class_<python::object>("AtomPairsParameters")
    .def_readonly("version",
		  RDKit::Descriptors::AtomPairs::atomPairsVersion)
    .def_readonly("numTypeBits",
		  RDKit::Descriptors::AtomPairs::numTypeBits)
    .def_readonly("numPiBits",
		  RDKit::Descriptors::AtomPairs::numPiBits)
    .def_readonly("numBranchBits",
		  RDKit::Descriptors::AtomPairs::numBranchBits)
    .def_readonly("codeSize",
		  RDKit::Descriptors::AtomPairs::codeSize)
    .def_readonly("atomTypes",
		  atomPairTypes)
    .def_readonly("numPathBits",
		  RDKit::Descriptors::AtomPairs::numPathBits)
    .def_readonly("numAtomPairFingerprintBits",
		  RDKit::Descriptors::AtomPairs::numAtomPairFingerprintBits)
    ;
  python::def("GetAtomPairAtomCode", RDKit::Descriptors::AtomPairs::getAtomCode,
              (python::arg("atom"), python::arg("branchSubtract")=0),
              docString.c_str());
  python::def("GetAtomPairFingerprint",
	      RDKit::Descriptors::AtomPairs::getAtomPairFingerprint,
	      python::arg("mol"),
              docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());

}
