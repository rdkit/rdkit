// $Id: rdForceFields.cpp 4967 2006-02-18 00:24:33Z glandrum $
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <boost/python.hpp>
#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap.h>


#include <ForceField/ForceField.h>
#include <ForceField/Wrap/PyForceField.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>

namespace python = boost::python;

namespace RDKit {
  int UFFOptimizeMolecule(ROMol &mol, int maxIters=200,
			  double vdwThresh=10.0, int confId=-1 ){
    ForceFields::ForceField *ff=UFF::constructForceField(mol,vdwThresh, confId);
    ff->initialize();
    int res=ff->minimize(maxIters);
    delete ff;
    return res;
  }
  ForceFields::PyForceField *UFFGetMoleculeForceField(ROMol &mol,
						    double vdwThresh=10.0,
						    int confId=-1 ){
    ForceFields::ForceField *ff=UFF::constructForceField(mol,vdwThresh, confId);
    ForceFields::PyForceField *res=new ForceFields::PyForceField(ff);
    res->initialize();
    return res;
  }
}

BOOST_PYTHON_MODULE(rdForceFieldHelpers) {
  python::scope().attr("__doc__") =
    "Module containing functions to handle molecular force fields (currently UFF)"
    ;

  std::string docString = "Use UFF to optimize a molecule's structure\n\n\
 \n\
 ARGUMENTS:\n\n\
    - mol : the molecule of interrest\n\
    - maxIters : the maximum number of iterations (defaults to 100)\n\
    - vdwThresh : used to exclude long-range van der Waals interactions\n\
                  (defaults to 10.0)\n\
    - confId : indicates which conformer to optimize\n\
\n";
  python::def("UFFOptimizeMolecule", RDKit::UFFOptimizeMolecule,
	      (python::arg("self"),python::arg("maxIters")=200,
	       python::arg("vdwThresh")=10.0,python::arg("confId")=-1),
	      docString.c_str());

 docString = "returns a UFF force field for a molecule\n\n\
 \n\
 ARGUMENTS:\n\n\
    - mol : the molecule of interrest\n\
    - vdwThresh : used to exclude long-range van der Waals interactions\n\
                  (defaults to 10.0)\n\
    - confId : indicates which conformer to optimize\n\
\n";
  python::def("UFFGetMoleculeForceField", RDKit::UFFGetMoleculeForceField,
	      (python::arg("mol"),python::arg("vdwThresh")=10.0,python::arg("confId")=-1),
	      python::return_value_policy<python::manage_new_object>(),
	      docString.c_str());

}
