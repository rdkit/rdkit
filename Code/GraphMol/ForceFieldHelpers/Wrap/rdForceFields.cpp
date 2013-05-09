// $Id$
//
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <boost/python.hpp>
#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap.h>


#include <ForceField/ForceField.h>
#include <ForceField/Wrap/PyForceField.h>
#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>

namespace python = boost::python;

namespace RDKit {
  int UFFOptimizeMolecule(ROMol &mol, int maxIters=200,
			  double vdwThresh=10.0, int confId=-1,
                          bool ignoreInterfragInteractions=true ){
    ForceFields::ForceField *ff=UFF::constructForceField(mol,vdwThresh, confId,
                                                         ignoreInterfragInteractions);
    ff->initialize();
    int res=ff->minimize(maxIters);
    delete ff;
    return res;
  }

  ForceFields::PyForceField *UFFGetMoleculeForceField(ROMol &mol,
                                                      double vdwThresh=10.0,
                                                      int confId=-1,
                                                      bool ignoreInterfragInteractions=true ){

    ForceFields::ForceField *ff=UFF::constructForceField(mol,vdwThresh, confId,
                                                         ignoreInterfragInteractions);
    ForceFields::PyForceField *res=new ForceFields::PyForceField(ff);
    res->initialize();
    return res;
  }

  bool UFFHasAllMoleculeParams(const ROMol &mol){
    UFF::AtomicParamVect types;
    bool foundAll;
    boost::tie(types,foundAll)=UFF::getAtomTypes(mol);
    return foundAll;
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
    - maxIters : the maximum number of iterations (defaults to 200)\n\
    - vdwThresh : used to exclude long-range van der Waals interactions\n\
                  (defaults to 10.0)\n\
    - confId : indicates which conformer to optimize\n\
    - ignoreInterfragInteractions : if true, nonbonded terms between \n\
                  fragments will not be added to the forcefield.\n\
\n\
 RETURNS: 0 if the optimization converged, 1 if more iterations are required.\n\
\n";
  python::def("UFFOptimizeMolecule", RDKit::UFFOptimizeMolecule,
	      (python::arg("self"),python::arg("maxIters")=200,
	       python::arg("vdwThresh")=10.0,python::arg("confId")=-1,
               python::arg("ignoreInterfragInteractions")=true),
	      docString.c_str());

 docString = "returns a UFF force field for a molecule\n\n\
 \n\
 ARGUMENTS:\n\n\
    - mol : the molecule of interrest\n\
    - vdwThresh : used to exclude long-range van der Waals interactions\n\
                  (defaults to 10.0)\n\
    - confId : indicates which conformer to optimize\n\
    - ignoreInterfragInteractions : if true, nonbonded terms between \n\
                  fragments will not be added to the forcefield.\n\
\n";
  python::def("UFFGetMoleculeForceField", RDKit::UFFGetMoleculeForceField,
	      (python::arg("mol"),python::arg("vdwThresh")=10.0,python::arg("confId")=-1,
               python::arg("ignoreInterfragInteractions")=true),
	      python::return_value_policy<python::manage_new_object>(),
	      docString.c_str());

 docString = "checks if UFF parameters are available for all of a molecule's atoms\n\n\
 \n\
 ARGUMENTS:\n\n\
    - mol : the molecule of interrest\n\
\n";
  python::def("UFFHasAllMoleculeParams", RDKit::UFFHasAllMoleculeParams,
	      (python::arg("mol")),
	      docString.c_str());

}
