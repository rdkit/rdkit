//
//  Copyright (C) 2004-2019 Greg Ladrum and Rational Discovery LLC
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
#include "props.hpp"

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <Geometry/point.h>
#include <GraphMol/Conformer.h>
#include <RDBoost/PySequenceHolder.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/python/numpy.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace python = boost::python;
namespace np = boost::python::numpy;

namespace RDKit {
RDGeom::Point3D GetAtomPos(const Conformer *conf, unsigned int aid) {
  RDGeom::Point3D res = conf->getAtomPos(aid);
  return res;
}

PyObject *GetPos(const Conformer *conf) {
  const RDGeom::POINT3D_VECT &pos = conf->getPositions();

  // define a 2D array with the following size
  npy_intp dims[2];
  dims[0] = pos.size();
  dims[1] = 3;

  // initialize the array
  auto *res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);

  // represent the array as a 1D/flat array of doubles
  auto *resData = reinterpret_cast<double *>(PyArray_DATA(res));

  // manually insert the data, 3 corresponds to the x, y and z dimensions
  for (unsigned int i = 0; i < pos.size(); ++i) {
    resData[3 * i + 0] = pos[i].x;
    resData[3 * i + 1] = pos[i].y;
    resData[3 * i + 2] = pos[i].z;
  }
  return PyArray_Return(res);
}

void SetPos(Conformer *conf, np::ndarray const & array) {
  if (array.get_dtype() != np::dtype::get_builtin<double>()) {
    PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
    python::throw_error_already_set();
  }

  if(array.get_nd() != 2) {
    PyErr_SetString(PyExc_TypeError, "Input array shape must be of rank 2");
    python::throw_error_already_set();
  }

  if(array.shape(0) != conf->getNumAtoms()) {
    PyErr_SetString(PyExc_ValueError, "Position array shape doesn't equal the number of atoms in the conformer");
    python::throw_error_already_set();
  }

  if(array.shape(1) < 2 || array.shape(1) > 3) {
    PyErr_SetString(PyExc_ValueError, "Position array point dimension must be 2 or 3 (2d or 3d)");
    python::throw_error_already_set();
  }

  const auto *data = reinterpret_cast<double*>(array.get_data());
  RDGeom::POINT3D_VECT &pos = conf->getPositions();
  if(array.shape(1) == 2) {
    for(size_t i=0; i<conf->getNumAtoms(); ++i) {
      pos[i].x = data[i*2];
      pos[i].y = data[i*2+1];
      pos[i].z = 0.0;
    }
  } else {
    for(size_t i=0; i<conf->getNumAtoms(); ++i) {
      pos[i].x = data[i*3];
      pos[i].y = data[i*3+1];
      pos[i].z = data[i*3+2];
    }
  }
  
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
    python::class_<Conformer, CONFORMER_SPTR>(
        "Conformer", confClassDoc.c_str(), python::init<>(python::args("self")))
        .def(python::init<unsigned int>(
            python::args("self", "numAtoms"),
            "Constructor with the number of atoms specified"))
        .def(python::init<const Conformer &>(python::args("self", "other")))

        .def("GetNumAtoms", &Conformer::getNumAtoms, python::args("self"),
             "Get the number of atoms in the conformer\n")

        .def("HasOwningMol", &Conformer::hasOwningMol, python::args("self"),
             "Returns whether or not this instance belongs to a molecule.\n")
        .def("GetOwningMol", &Conformer::getOwningMol,
             "Get the owning molecule\n",
             python::return_value_policy<python::reference_existing_object>(),
             python::args("self"))

        .def("GetId", &Conformer::getId, python::args("self"),
             "Get the ID of the conformer")
        .def("SetId", &Conformer::setId, python::args("self", "id"),
             "Set the ID of the conformer\n")

        .def("GetAtomPosition", GetAtomPos, python::args("self", "aid"),
             "Get the posistion of an atom\n")
        .def("GetPositions", GetPos, python::args("self"),
             "Get positions of all the atoms\n")
        .def("SetPositions", SetPos,
	     (python::args("self"), python::args("positions")),
             "Set positions of all the atoms given a 2D or 3D numpy array of type double\n")
        .def("SetAtomPosition", SetAtomPos, python::args("self", "aid", "loc"),
             "Set the position of the specified atom\n")
        .def("SetAtomPosition",
             (void(Conformer::*)(unsigned int, const RDGeom::Point3D &)) &
                 Conformer::setAtomPos,
             python::args("self", "atomId", "position"),
             "Set the position of the specified atom\n")

        .def("Set3D", &Conformer::set3D, python::args("self", "v"),
             "Set the 3D flag of the conformer\n")
        .def("Is3D", &Conformer::is3D, python::args("self"),
             "returns the 3D flag of the conformer\n")

        // properties
        .def("SetProp", MolSetProp<Conformer, std::string>,
             (python::arg("self"), python::arg("key"), python::arg("val"),
              python::arg("computed") = false),
             "Sets a molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value (a string).\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to False.\n\n")
        .def("SetDoubleProp", MolSetProp<Conformer, double>,
             (python::arg("self"), python::arg("key"), python::arg("val"),
              python::arg("computed") = false),
             "Sets a double valued molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value as a double.\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to 0.\n\n")
        .def("SetIntProp", MolSetProp<Conformer, int>,
             (python::arg("self"), python::arg("key"), python::arg("val"),
              python::arg("computed") = false),
             "Sets an integer valued molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (an unsigned "
             "number).\n"
             "    - value: the property value as an integer.\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to False.\n\n")
        .def("SetUnsignedProp", MolSetProp<Conformer, unsigned int>,
             (python::arg("self"), python::arg("key"), python::arg("val"),
              python::arg("computed") = false),
             "Sets an unsigned integer valued molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value as an unsigned integer.\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to False.\n\n")
        .def("SetBoolProp", MolSetProp<Conformer, bool>,
             (python::arg("self"), python::arg("key"), python::arg("val"),
              python::arg("computed") = false),
             "Sets a boolean valued molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value as a bool.\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to False.\n\n")
        .def("HasProp", MolHasProp<Conformer>, python::args("self", "key"),
             "Queries a conformer to see if a particular property has been "
             "assigned.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to check for (a string).\n")
        .def(
            "GetProp", GetPyProp<Conformer>,
            (python::arg("self"), python::arg("key"),
             python::arg("autoConvert") = false),
            "Returns the value of the property.\n\n"
            "  ARGUMENTS:\n"
            "    - key: the name of the property to return (a string).\n\n"
            "    - autoConvert: if True attempt to convert the property into a python object\n\n"
            "  RETURNS: a string\n\n"
            "  NOTE:\n"
            "    - If the property has not been set, a KeyError exception "
            "will be raised.\n")
        .def("GetDoubleProp", GetProp<Conformer, double>,
             python::args("self", "key"),
             "Returns the double value of the property if possible.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: a double\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")
        .def("GetIntProp", GetProp<Conformer, int>, python::args("self", "key"),
             "Returns the integer value of the property if possible.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: an integer\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")
        .def("GetUnsignedProp", GetProp<Conformer, unsigned int>,
             python::args("self", "key"),
             "Returns the unsigned int value of the property if possible.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: an unsigned integer\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")
        .def("GetBoolProp", GetProp<Conformer, bool>,
             python::args("self", "key"),
             "Returns the Bool value of the property if possible.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: a bool\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")
        .def("ClearProp", MolClearProp<Conformer>, python::args("self", "key"),
             "Removes a property from the conformer.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to clear (a string).\n")

        .def("ClearComputedProps", MolClearComputedProps<Conformer>,
             python::args("self"),
             "Removes all computed properties from the conformer.\n\n")
        .def("GetPropNames", &Conformer::getPropList,
             (python::arg("self"), python::arg("includePrivate") = false,
              python::arg("includeComputed") = false),
             "Returns a tuple with all property names for this conformer.\n\n"
             "  ARGUMENTS:\n"
             "    - includePrivate: (optional) toggles inclusion of private "
             "properties in the result set.\n"
             "                      Defaults to 0.\n"
             "    - includeComputed: (optional) toggles inclusion of computed "
             "properties in the result set.\n"
             "                      Defaults to 0.\n\n"
             "  RETURNS: a tuple of strings\n")

        .def("GetPropsAsDict", GetPropsAsDict<Conformer>,
             (python::arg("self"), python::arg("includePrivate") = false,
              python::arg("includeComputed") = false,
              python::arg("autoConvertStrings") = true),
             "Returns a dictionary populated with the conformer's properties.\n"
             " n.b. Some properties are not able to be converted to python "
             "types.\n\n"
             "  ARGUMENTS:\n"
             "    - includePrivate: (optional) toggles inclusion of private "
             "properties in the result set.\n"
             "                      Defaults to False.\n"
             "    - includeComputed: (optional) toggles inclusion of computed "
             "properties in the result set.\n"
             "                      Defaults to False.\n\n"
             "  RETURNS: a dictionary\n");
  };
};
}  // namespace RDKit

void wrap_conformer() { RDKit::conformer_wrapper::wrap(); }
