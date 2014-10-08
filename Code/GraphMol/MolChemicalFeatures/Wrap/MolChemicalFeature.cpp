// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define NO_IMPORT_ARRAY
#include <boost/python.hpp>

#include <RDGeneral/types.h>
#include <RDBoost/pyint_api.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>

namespace python = boost::python;

namespace RDKit {
  PyObject *getFeatAtomIds(const MolChemicalFeature &feat) {
    const MolChemicalFeature::AtomPtrContainer &atoms = feat.getAtoms();
    PyObject *res = PyTuple_New(atoms.size());
    MolChemicalFeature::AtomPtrContainer_CI aci;
    int idx=0;
    for (aci = atoms.begin(); aci != atoms.end(); ++aci) {
      PyTuple_SetItem(res,idx,PyInt_FromLong((*aci)->getIdx()));
      idx++;
    }
    return res;
  }

  std::string featClassDoc="Class to represent a chemical feature.\n\
    These chemical features may or may not have been derived from molecule object;\n\
    i.e. it is possible to have a chemical feature that was created just from its type\n\
    and location.\n";
  struct feat_wrapper {
    static void wrap() {
      python::class_<MolChemicalFeature,FeatSPtr>("MolChemicalFeature", featClassDoc.c_str(),
						  python::no_init)
						  
        .def("GetId", &MolChemicalFeature::getId,
	         "Returns the identifier of the feature\n")
        .def("GetFamily", &MolChemicalFeature::getFamily,
             "Get the family to which the feature belongs; donor, acceptor, etc.",
             python::return_value_policy<python::copy_const_reference>())
        .def("GetType",  &MolChemicalFeature::getType,
             "Get the specific type for the feature",
             python::return_value_policy<python::copy_const_reference>())
        .def("GetPos",
	     (RDGeom::Point3D (MolChemicalFeature::*)(int) const)&MolChemicalFeature::getPos,
             (python::arg("self"),python::arg("confId")=-1),
             "Get the location of the chemical feature")
        .def("GetAtomIds", getFeatAtomIds,
             "Get the IDs of the atoms that participate in the feature")
        .def("GetMol", &MolChemicalFeature::getMol,
             "Get the molecule used to derive the features",
             python::return_value_policy<python::reference_existing_object>())
        .def("GetFactory", &MolChemicalFeature::getFactory,
             "Get the factory used to generate this feature",
             python::return_value_policy<python::reference_existing_object>())
	.def("ClearCache",&MolChemicalFeature::clearCache,
	     "Clears the cache used to store position information.")
        ;
    };
  };

}
void wrap_MolChemicalFeat() {
  RDKit::feat_wrapper::wrap();
}

    
