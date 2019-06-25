//
//  Copyright (C) 2018 Dan Nealschneider
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#include <boost/python.hpp>
#include <string>

// ours
#include <RDBoost/Wrap.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/StereoGroup.h>

namespace python = boost::python;

using boost::python::converter::rvalue_from_python_stage1_data;

namespace RDKit {
namespace {
std::string stereoGroupClassDoc =
    "A collection of atoms with a defined stereochemical relationship.\n\n"
    "Used to help represent a sample with unknown stereochemistry, or that "
    "is a mix\nof diastereomers.\n";

static StereoGroup stereoGroupFactory(StereoGroupType grouptype,
                                      python::object atoms) {
  std::vector<Atom*> cppAtoms;
  pythonObjectToVect<Atom*>(atoms, cppAtoms);
  return {grouptype, std::move(cppAtoms)};
}

// Is obj convertable to a vector?
template <typename T>
void* vector_convertable(PyObject* obj) {
  return PyList_Check(obj) ? obj : nullptr;
}

// Convert obj to a vector.
template <typename T>
void vector_convert(PyObject* obj, rvalue_from_python_stage1_data* data) {
  // See
  // https://www.boost.org/doc/libs/1_61_0/libs/python/doc/html/faq/how_can_i_automatically_convert_.html

  // Boost uses a borrowed reference here.
  boost::python::object borrowed_obj(
      boost::python::handle<>(boost::python::borrowed(obj)));
  python::stl_input_iterator<T> beg(borrowed_obj), end;

  // Allocate the C++ type into the converter's memory block
  void* storage =
      ((boost::python::converter::rvalue_from_python_storage<std::vector<T>>*)
           data)
          ->storage.bytes;
  new (storage) std::vector<T>(beg, end);
  data->convertible = storage;
}

// Wrap vector<T> for use in Python
//
// C++ to Python: Python gets a wrapped C++ vector<T> using the
// boost vector_indexing_suite.
//
// Python to C++: Either a wrapped C++ vector or a native Python list is
// allowed.
template <typename T>
void register_vector2way(const char* name) {
  RegisterVectorConverter<T>(name);

  boost::python::converter::registry::push_back(
      &vector_convertable<T>, &vector_convert<T>,
      boost::python::type_id<std::vector<T>>());
}

StereoGroup* createAndReturnStereoGroup(ROMol* mol, python::object atomIds,
                                        StereoGroupType typ) {
  PRECONDITION(mol, "no molecule");
  std::vector<Atom*> cppAtoms;
  python::stl_input_iterator<unsigned int> beg(atomIds), end;
  while (beg != end) {
    unsigned int v = *beg;
    if (v >= mol->getNumAtoms())
      throw_value_error("atom index exceeds mol.GetNumAtoms()");
    cppAtoms.push_back(mol->getAtomWithIdx(v));
    ++beg;
  }
  StereoGroup* sg = new StereoGroup(typ, cppAtoms);
  return sg;
}

}  // namespace
struct stereogroup_wrap {
  static void wrap() {
    python::enum_<RDKit::StereoGroupType>("StereoGroupType")
        .value("STEREO_ABSOLUTE", RDKit::StereoGroupType::STEREO_ABSOLUTE)
        .value("STEREO_OR", RDKit::StereoGroupType::STEREO_OR)
        .value("STEREO_AND", RDKit::StereoGroupType::STEREO_AND)
        .export_values();

    register_vector2way<Atom*>("Atom_vect");

    python::class_<StereoGroup, boost::shared_ptr<StereoGroup>>(
        "StereoGroup", stereoGroupClassDoc.c_str(), python::no_init)
        .def(
            python::init<RDKit::StereoGroupType, const ROMol::ATOM_PTR_VECT&>())
        .def("GetGroupType", &StereoGroup::getGroupType,
             "Returns the StereoGroupType.\n")
        .def("GetAtoms", &StereoGroup::getAtoms,
             "Access the atoms in the StereoGroup.\n",
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1>>());
    python::def("CreateStereoGroup", &createAndReturnStereoGroup,
                "creates a StereoGroup associated with a molecule from a list "
                "of atom Ids",
                (python::arg("mol"), python::arg("atomIds"),
                 python::arg("stereoGroupType")),
                python::return_value_policy<
                    python::manage_new_object,
                    python::with_custodian_and_ward_postcall<0, 1>>());
  }
};
}  // namespace RDKit

void wrap_stereogroup() { RDKit::stereogroup_wrap::wrap(); }
