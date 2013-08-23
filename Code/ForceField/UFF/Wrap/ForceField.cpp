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

}
