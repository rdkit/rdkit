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

#include <boost/python.hpp>
#include <RDBoost/Wrap.h>
#include <RDBoost/Exceptions.h>
#include <ForceField/ForceField.h>
#include <ForceField/UFF/DistanceConstraint.h>
#include <boost/shared_ptr.hpp>
#include "PyForceField.h"
using namespace ForceFields;
namespace python = boost::python;

void ForceFieldAddDistanceConstraint(PyForceField *self,unsigned int idx1,unsigned int idx2,
				     double minLen,double maxLen,double forceConstant){
  UFF::DistanceConstraintContrib *constraint;
  constraint = new UFF::DistanceConstraintContrib(self->field.get(),idx1,idx2,minLen,maxLen,
						  forceConstant);
  self->field->contribs().push_back(ForceFields::ContribPtr(constraint));
}

PyObject *ForceFieldGetExtraPointLoc(PyForceField *self,unsigned int idx){
  if(idx >= self->extraPoints.size()){
    throw IndexErrorException(idx);
  }
  PyObject *res = PyTuple_New(3);
  PyTuple_SetItem(res,0,PyFloat_FromDouble(self->extraPoints[idx]->x));
  PyTuple_SetItem(res,1,PyFloat_FromDouble(self->extraPoints[idx]->y));
  PyTuple_SetItem(res,2,PyFloat_FromDouble(self->extraPoints[idx]->z));
  return res;
}

BOOST_PYTHON_MODULE(rdForceField) {
  python::scope().attr("__doc__") =
    "Exposes the ForceField class"
    ;
  
  std::string docString;

  python::class_<PyForceField>("ForceField","A force field",python::no_init)
    .def("CalcEnergy",(double (PyForceField::*)() const)&PyForceField::calcEnergy,
	 "Returns the energy of the current arrangement")
    .def("Minimize",&PyForceField::minimize,(python::arg("maxIts")=200,
					   python::arg("forceTol")=1e-4,
					   python::arg("energyTol")=1e-6),
	 "Runs some minimization iterations.\n\n  Returns 0 if the minimization succeeded.")
    .def("AddDistanceConstraint",ForceFieldAddDistanceConstraint,
	 (python::arg("self"),python::arg("idx1"),python::arg("idx2"),
	  python::arg("minLen"),python::arg("maxLen"),
	  python::arg("forceConstant")),
	 "Adds a distance constraint to the force field.")
    .def("Initialize",&PyForceField::initialize,
	 "initializes the force field (call this before minimizing)")
    .def("AddExtraPoint",&PyForceField::addExtraPoint,
	 (python::arg("self"),python::arg("x"),python::arg("y"),python::arg("z"),
	  python::arg("fixed")=true),
	 "Adds an extra point, this can be useful for adding constraints.")
    .def("GetExtraPointPos",ForceFieldGetExtraPointLoc,
	 (python::arg("self"),python::arg("idx")),
	 "returns the location of an extra point as a tuple")
    ;
  python::class_<PyMMFFMolProperties>("MMFFMolProperties",
    "MMFF molecular properties", python::no_init)
    .def("GetMMFFAtomType", &PyMMFFMolProperties::getMMFFAtomType,
    (python::arg("self"), python::arg("idx")),
    "Retrieves MMFF atom type for atom with index idx")
    .def("GetMMFFFormalCharge", &PyMMFFMolProperties::getMMFFFormalCharge,
    (python::arg("self"), python::arg("idx")),
    "Retrieves MMFF formal charge for atom with index idx")
    .def("GetMMFFPartialCharge", &PyMMFFMolProperties::getMMFFPartialCharge,
    (python::arg("self"), python::arg("idx")),
    "Retrieves MMFF partial charge for atom with index idx")
    .def("SetMMFFDielectricModel", &PyMMFFMolProperties::setMMFFDielectricModel,
    (python::arg("self"), python::arg("dielModel") = 1),
    "Sets the DielModel MMFF property (1: constant; 2: distance-dependent; "
    "defaults to constant)")
    .def("SetMMFFDielectricConstant", &PyMMFFMolProperties::setMMFFDielectricConstant,
    (python::arg("self"), python::arg("dielConst") = 1.0),
    "Sets the DielConst MMFF property (defaults to 1.0)")
    .def("SetMMFFBondTerm", &PyMMFFMolProperties::setMMFFBondTerm,
    (python::arg("self"), python::arg("state") = true),
    "Sets the bond term to be included in the MMFF equation (defaults to True)")
    .def("SetMMFFAngleTerm", &PyMMFFMolProperties::setMMFFAngleTerm,
    (python::arg("self"), python::arg("state") = true),
    "Sets the angle term to be included in the MMFF equation (defaults to True)")
    .def("SetMMFFStretchBendTerm", &PyMMFFMolProperties::setMMFFStretchBendTerm,
    (python::arg("self"), python::arg("state") = true),
    "Sets the stretch-bend term to be included in the MMFF equation (defaults to True)")
    .def("SetMMFFOopTerm", &PyMMFFMolProperties::setMMFFOopTerm,
    (python::arg("self"), python::arg("state") = true),
    "Sets the out-of-plane bend term to be included in the MMFF equation (defaults to True)")
    .def("SetMMFFTorsionTerm", &PyMMFFMolProperties::setMMFFTorsionTerm,
    (python::arg("self"), python::arg("state") = true),
    "Sets the torsional term to be included in the MMFF equation (defaults to True)")
    .def("SetMMFFVdWTerm", &PyMMFFMolProperties::setMMFFVdWTerm,
    (python::arg("self"), python::arg("state") = true),
    "Sets the Van der Waals term to be included in the MMFF equation (defaults to True)")
    .def("SetMMFFEleTerm", &PyMMFFMolProperties::setMMFFEleTerm,
    (python::arg("self"), python::arg("state") = true),
    "Sets the electrostatic term to be included in the MMFF equation (defaults to True)")
    .def("SetMMFFVariant", &PyMMFFMolProperties::setMMFFVariant,
    (python::arg("self"), python::arg("mmffVariant") = "MMFF94"),
    "Sets the MMFF variant to be used (\"MMFF94\" or \"MMFF94s\"; defaults to \"MMFF94\")")
    .def("SetMMFFVerbosity", &PyMMFFMolProperties::setMMFFVerbosity,
    (python::arg("self"), python::arg("verbosity") = 0),
    "Sets the MMFF verbosity (0: none; 1: low; 2: high; defaults to 0)")
    ;

}
