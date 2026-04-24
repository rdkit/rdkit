//
//  Copyright (C) 2026 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <string>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/ndarray.h>

#include "props.hpp"

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <Geometry/point.h>
#include <GraphMol/Conformer.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {

nb::ndarray<nb::numpy, double, nb::shape<-1, 3>> GetPos(const Conformer *conf) {
  const RDGeom::POINT3D_VECT &pos = conf->getPositions();

  auto np = nb::module_::import_("numpy");
  auto arrObj =
      np.attr("empty")(nb::make_tuple(pos.size(), 3), "dtype"_a = "float64");
  auto arr =
      nb::cast<nb::ndarray<nb::numpy, double, nb::shape<-1, 3>, nb::c_contig>>(
          arrObj);
  auto *resData = arr.data();

  // Fill an array of shape (numAtoms,3) with XYZ coordinates.
  for (size_t i = 0; i < pos.size(); ++i) {
    resData[3 * i + 0] = pos[i].x;
    resData[3 * i + 1] = pos[i].y;
    resData[3 * i + 2] = pos[i].z;
  }
  return nb::cast<nb::ndarray<nb::numpy, double, nb::shape<-1, 3>>>(arrObj);
}

void SetPos(Conformer *conf,
            const nb::ndarray<nb::numpy, const double, nb::shape<-1, -1>,
                              nb::c_contig> &array) {
  if ((size_t)array.shape(0) != conf->getNumAtoms()) {
    throw ValueErrorException(
        "Position array shape doesn't equal the number of atoms in the conformer");
  }

  if (array.shape(1) < 2 || array.shape(1) > 3) {
    throw ValueErrorException(
        "Position array point dimension must be 2 or 3 (2d or 3d)");
  }

  RDGeom::POINT3D_VECT &pos = conf->getPositions();
  const double *data = array.data();
  const size_t width = (size_t)array.shape(1);

  if (array.shape(1) == 2) {
    for (size_t i = 0; i < conf->getNumAtoms(); ++i) {
      pos[i].x = data[i * width];
      pos[i].y = data[i * width + 1];
      pos[i].z = 0.0;
    }
  } else {
    for (size_t i = 0; i < conf->getNumAtoms(); ++i) {
      pos[i].x = data[i * width];
      pos[i].y = data[i * width + 1];
      pos[i].z = data[i * width + 2];
    }
  }
}

void SetAtomPos(Conformer *conf, unsigned int aid, nb::object loc) {
  try {
    std::vector<double> coords;
    for (auto item : loc) {
      coords.push_back(nb::cast<double>(item));
    }
    if (coords.size() != 3) {
      throw ValueErrorException("Atom position must have 3 coordinates");
    }
    RDGeom::Point3D pt(coords[0], coords[1], coords[2]);
    conf->setAtomPos(aid, pt);
  } catch (const nb::cast_error &) {
    throw ValueErrorException(
        "Could not convert position coordinates to doubles");
  }
}

std::string confClassDoc =
    "The class to store 2D or 3D conformation of a molecule\n";

struct conformer_wrapper {
  static void wrap(nb::module_ &m) {
    nb::class_<Conformer>(m, "Conformer")
        .def(nb::init<>(), "Constructor with the number of atoms specified")
        .def(nb::init<unsigned int>(), "numAtoms"_a,
             "Constructor with the number of atoms specified")
        .def(nb::init<const Conformer &>(), "other"_a)

        .def("GetNumAtoms", &Conformer::getNumAtoms,
             "Get the number of atoms in the conformer\n")

        .def("HasOwningMol", &Conformer::hasOwningMol,
             "Returns whether or not this instance belongs to a molecule.\n")
        .def("GetOwningMol", &Conformer::getOwningMol,
             "Get the owning molecule\n", nb::rv_policy::reference_internal)

        .def("GetId", &Conformer::getId, "Get the ID of the conformer")
        .def("SetId", &Conformer::setId, "id"_a,
             "Set the ID of the conformer\n")

        .def(
            "GetAtomPosition",
            [](const Conformer *conf, unsigned int aid) {
              return conf->getAtomPos(aid);
            },
            "aid"_a, "Get the position of an atom\n")
        .def("GetPositions", GetPos, "Get positions of all the atoms\n")
        .def(
            "SetPositions", SetPos, "positions"_a,
            "Set positions of all the atoms given a 2D or 3D numpy array of type double\n")
        .def("SetAtomPosition", SetAtomPos, "aid"_a, "loc"_a,
             "Set the position of the specified atom\n")
        .def("SetAtomPosition",
             (void (Conformer::*)(
                 unsigned int, const RDGeom::Point3D &))&Conformer::setAtomPos,
             "atomId"_a, "position"_a,
             "Set the position of the specified atom\n")

        .def("Set3D", &Conformer::set3D, "v"_a,
             "Set the 3D flag of the conformer\n")
        .def("Is3D", &Conformer::is3D, "returns the 3D flag of the conformer\n")

        // properties
        .def("SetProp", MolSetProp<Conformer, std::string>, "key"_a, "val"_a,
             "computed"_a = false,
             "Sets a molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value (a string).\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to False.\n\n")
        .def("SetDoubleProp", MolSetProp<Conformer, double>, "key"_a, "val"_a,
             "computed"_a = false,
             "Sets a double valued molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value as a double.\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to 0.\n\n")
        .def("SetIntProp", MolSetProp<Conformer, int>, "key"_a, "val"_a,
             "computed"_a = false,
             "Sets an integer valued molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (an unsigned "
             "number).\n"
             "    - value: the property value as an integer.\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to False.\n\n")
        .def("SetUnsignedProp", MolSetProp<Conformer, unsigned int>, "key"_a,
             "val"_a, "computed"_a = false,
             "Sets an unsigned integer valued molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value as an unsigned integer.\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to False.\n\n")
        .def("SetBoolProp", MolSetProp<Conformer, bool>, "key"_a, "val"_a,
             "computed"_a = false,
             "Sets a boolean valued molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value as a bool.\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to False.\n\n")
        .def("HasProp", MolHasProp<Conformer>, "key"_a,
             "Queries a conformer to see if a particular property has been "
             "assigned.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to check for (a string).\n")
        .def(
            "GetProp", GetPyProp<Conformer>, "key"_a, "autoConvert"_a = false,
            "Returns the value of the property.\n\n"
            "  ARGUMENTS:\n"
            "    - key: the name of the property to return (a string).\n\n"
            "    - autoConvert: if True attempt to convert the property into a python object\n\n"
            "  RETURNS: a string\n\n"
            "  NOTE:\n"
            "    - If the property has not been set, a KeyError exception "
            "will be raised.\n")
        .def("GetDoubleProp", GetProp<Conformer, double>, "key"_a,
             "Returns the double value of the property if possible.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: a double\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")
        .def("GetIntProp", GetProp<Conformer, int>, "key"_a,
             "Returns the integer value of the property if possible.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: an integer\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")
        .def("GetUnsignedProp", GetProp<Conformer, unsigned int>, "key"_a,
             "Returns the unsigned int value of the property if possible.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: an unsigned integer\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")
        .def("GetBoolProp", GetProp<Conformer, bool>, "key"_a,
             "Returns the Bool value of the property if possible.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: a bool\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")
        .def("ClearProp", MolClearProp<Conformer>, "key"_a,
             "Removes a property from the conformer.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to clear (a string).\n")

        .def("ClearComputedProps", MolClearComputedProps<Conformer>,
             "Removes all computed properties from the conformer.\n\n")
        .def("GetPropNames", &Conformer::getPropList,
             "includePrivate"_a = false, "includeComputed"_a = false,
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
             "includePrivate"_a = false, "includeComputed"_a = false,
             "autoConvertStrings"_a = true, getPropsAsDictDocString.c_str())
        .doc() = confClassDoc.c_str();
  };
};
}  // namespace RDKit

void wrap_conformer(nb::module_ &m) { RDKit::conformer_wrapper::wrap(m); }
