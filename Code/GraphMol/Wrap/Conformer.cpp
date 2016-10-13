// $Id$
//
//  Copyright (C) 2004-2008 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#include <RDBoost/python.h>
#include <string>
#include "rdchem.h"
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <Geometry/point.h>
#include <GraphMol/Conformer.h>
#include <RDBoost/PySequenceHolder.h>
//#include "seqs.hpp"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
namespace python = boost::python;

namespace RDKit {
RDGeom::Point3D GetAtomPos(const Conformer *conf, unsigned int aid) {
  RDGeom::Point3D res = conf->getAtomPos(aid);
  return res;
}

PyObject* GetPos(const Conformer *conf) {
    const RDGeom::POINT3D_VECT &pos = conf->getPositions();

    // define a 2D array with the following size
    npy_intp dims[2];
    dims[0] = pos.size();
    dims[1] = 3;

    // initialize the array
    PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);

    // represent the array as a 1D/flat array of doubles
    double *resData = reinterpret_cast<double *>(PyArray_DATA(res));

    // manually insert the data, 3 corresponsd to the x, y and z dimensions
    for (unsigned int i = 0; i < pos.size(); ++i) {
        resData[3*i+0] = pos[i].x;
        resData[3*i+1] = pos[i].y;
        resData[3*i+2] = pos[i].z;
    }
    return PyArray_Return(res);
}

void SetAtomPos(Conformer *conf, unsigned int aid, python::object loc) {
  // const std::vector<double> &loc) {
  int dim = python::extract<int>(loc.attr("__len__")());
  CHECK_INVARIANT(dim == 3, "");
  PySequenceHolder<double> pdata(loc);
  RDGeom::Point3D pt(pdata[0], pdata[1], pdata[2]);
  conf->setAtomPos(aid, pt);
}

std::string confClassDoc =
    "The class to store 2D or 3D conformation of a molecule\n";

struct conformer_wrapper {
  static void wrap() {
    python::class_<Conformer, CONFORMER_SPTR>("Conformer", confClassDoc.c_str(),
                                              python::init<>())
        .def(python::init<unsigned int>(
            "Constructor with the number of atoms specified"))
        .def(python::init<const Conformer &>())

        .def("GetNumAtoms", &Conformer::getNumAtoms,
             "Get the number of atoms in the conformer\n")

        .def("GetOwningMol", &Conformer::getOwningMol,
             "Get the owning molecule\n",
             python::return_value_policy<python::reference_existing_object>())

        .def("GetId", &Conformer::getId, "Get the ID of the conformer")
        .def("SetId", &Conformer::setId, "Set the ID of the conformer\n")

        .def("GetAtomPosition", GetAtomPos, "Get the posistion of an atom\n")
	.def("GetPositions", GetPos, "Get positions of all the atoms\n")
        .def("SetAtomPosition", SetAtomPos,
             "Set the position of the specified atom\n")
        .def("SetAtomPosition", (void (Conformer::*)(unsigned int, const RDGeom::Point3D&)) &
                               Conformer::setAtomPos,
             "Set the position of the specified atom\n")

        .def("Set3D", &Conformer::set3D, "Set the 3D flag of the conformer\n")
        .def("Is3D", &Conformer::is3D,
             "returns the 3D flag of the conformer\n");
  };
};
}

void wrap_conformer() { RDKit::conformer_wrapper::wrap(); }
