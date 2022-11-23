//
//  Copyright (C) 2018-2025 Dan Nealschneider and other RDKit contributors
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

namespace RDKit {

namespace {
std::string stereoGroupClassDoc =
    "A collection of atoms with a defined stereochemical relationship.\n\n"
    "Used to help represent a sample with unknown stereochemistry, or that "
    "is a mix\nof diastereomers.\n";

StereoGroup *createStereoGroup(StereoGroupType typ, ROMol &mol,
                               python::object atomIds, python::object bondIds, unsigned readId) {
  std::vector<Atom *> cppAtoms;
  std::vector<Bond *> cppBonds;
  python::stl_input_iterator<unsigned int> beg(atomIds), end;
  while (beg != end) {
    unsigned int v = *beg;
    if (v >= mol.getNumAtoms()) {
      throw_value_error("atom index exceeds mol.GetNumAtoms()");
    }
    cppAtoms.push_back(mol.getAtomWithIdx(v));
    ++beg;
  }
  python::stl_input_iterator<unsigned int> bbeg(bondIds), bend;
  while (bbeg != bend) {
    unsigned int v = *bbeg;
    if (v >= mol.getNumBonds()) {
      throw_value_error("bond index exceeds mol.GetNumBonds()");
    }
    cppBonds.push_back(mol.getBondWithIdx(v));
    ++bbeg;
  }
  if (cppAtoms.empty() && cppBonds.empty()) {
    throw_value_error("New StereoGroup must contain at least one atom or bond.");
  }
  auto *sg = new StereoGroup(typ, cppAtoms, cppBonds, readId);
  return sg;
}

python::object getAtomsHelper(StereoGroup &sg) {
  python::list res;
  for (auto at : sg.getAtoms()) {
    res.append(boost::ref(*at));
  }
  return python::tuple(res);
}
python::object getBondsHelper(StereoGroup &sg) {
  python::list res;
  for (auto bnd : sg.getBonds()) {
    res.append(boost::ref(*bnd));
  }
  return python::tuple(res);
}
}  // namespace

struct stereogroup_wrap {
  static void wrap() {
    python::enum_<RDKit::StereoGroupType>("StereoGroupType")
        .value("STEREO_ABSOLUTE", RDKit::StereoGroupType::STEREO_ABSOLUTE)
        .value("STEREO_OR", RDKit::StereoGroupType::STEREO_OR)
        .value("STEREO_AND", RDKit::StereoGroupType::STEREO_AND)
        .export_values();

    python::class_<StereoGroup, std::shared_ptr<StereoGroup>>(
        "StereoGroup", stereoGroupClassDoc.c_str(), python::no_init)
        .def("GetGroupType", &StereoGroup::getGroupType, python::args("self"),
             "Returns the StereoGroupType.\n")
        .def("GetAtoms", getAtomsHelper, python::args("self"),
             "access the atoms in the StereoGroup.\n")
        .def("GetBonds", getBondsHelper, python::args("self"),
             "access the bonds in the StereoGroup.\n")
        .def("GetReadId", &StereoGroup::getReadId, python::args("self"),
             "return the StereoGroup's original ID.\n"
             "Note that the ID only makes sense for AND/OR groups.\n")
        .def("GetWriteId", &StereoGroup::getWriteId, python::args("self"),
             "return the StereoGroup's ID that will be exported.\n"
             "Note that the ID only makes sense for AND/OR groups.\n")
        .def("SetWriteId", &StereoGroup::setWriteId, python::args("self", "id"),
             "return the StereoGroup's ID that will be exported.\n"
             "Note that the ID only makes sense for AND/OR groups.\n");

    python::def("CreateStereoGroup", &createStereoGroup,
                "creates a StereoGroup associated with a molecule from a list "
                "of atom Ids",
                (python::arg("stereoGroupType"), python::arg("mol"),
                 python::arg("atomIds") = boost::python::list(), python::arg("bondIds") = boost::python::list(),
                 python::arg("readId") = 0),
                python::return_value_policy<
                    python::manage_new_object,
                    python::with_custodian_and_ward_postcall<0, 2>>());

    python::def(
        "ForwardStereoGroupIds", &RDKit::forwardStereoGroupIds,
        python::args("mol"),
        "Forward the original Stereo Group IDs when exporting the Mol.");
  }
};
}  // namespace RDKit

void wrap_stereogroup() { RDKit::stereogroup_wrap::wrap(); }
