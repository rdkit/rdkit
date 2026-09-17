//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>
#include <string>

// ours
#include <RDBoost/Wrap_nb.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/StereoGroup.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {

namespace {
std::string stereoGroupClassDoc =
    R"DOC(A collection of atoms with a defined stereochemical relationship.

Used to help represent a sample with unknown stereochemistry, or that is a mix
of diastereomers.
)DOC";

StereoGroup *createStereoGroup(StereoGroupType typ, ROMol &mol,
                               nb::object atomIds, nb::object bondIds,
                               unsigned readId) {
  std::vector<Atom *> cppAtoms;
  std::vector<Bond *> cppBonds;
  if (!atomIds.is_none()) {
    auto atomIdVec =
        pythonObjectToVect<unsigned int>(atomIds, mol.getNumAtoms());
    if (atomIdVec) {
      for (const auto v : *atomIdVec) {
        cppAtoms.push_back(mol.getAtomWithIdx(v));
      }
    }
  }
  if (!bondIds.is_none()) {
    auto bondIdVec =
        pythonObjectToVect<unsigned int>(bondIds, mol.getNumBonds());
    if (bondIdVec) {
      for (const auto v : *bondIdVec) {
        cppBonds.push_back(mol.getBondWithIdx(v));
      }
    }
  }
  if (cppAtoms.empty() && cppBonds.empty()) {
    throw ValueErrorException(
        "New StereoGroup must contain at least one atom or bond.");
  }
  auto *sg = new StereoGroup(typ, cppAtoms, cppBonds, readId);
  return sg;
}

nb::tuple getAtomsHelper(StereoGroup &sg) {
  nb::list res;
  for (auto at : sg.getAtoms()) {
    res.append(at);
  }
  return nb::tuple(res);
}
nb::tuple getBondsHelper(StereoGroup &sg) {
  nb::list res;
  for (auto bnd : sg.getBonds()) {
    res.append(bnd);
  }
  return nb::tuple(res);
}
}  // namespace

struct stereogroup_wrap {
  static void wrap(nb::module_ &m) {
    nb::enum_<RDKit::StereoGroupType>(m, "StereoGroupType")
        .value("STEREO_ABSOLUTE", RDKit::StereoGroupType::STEREO_ABSOLUTE)
        .value("STEREO_OR", RDKit::StereoGroupType::STEREO_OR)
        .value("STEREO_AND", RDKit::StereoGroupType::STEREO_AND)
        .export_values();

    nb::class_<StereoGroup>(m, "StereoGroup", stereoGroupClassDoc.c_str())
        .def("GetGroupType", &StereoGroup::getGroupType,
             R"DOC(Returns the StereoGroupType.)DOC")
        .def("GetAtoms", getAtomsHelper,
             R"DOC(access the atoms in the StereoGroup.)DOC")
        .def("GetBonds", getBondsHelper,
             R"DOC(access the bonds in the StereoGroup.)DOC")
        .def("GetReadId", &StereoGroup::getReadId,
             R"DOC(return the StereoGroup's original ID.
Note that the ID only makes sense for AND/OR groups.)DOC")
        .def("GetWriteId", &StereoGroup::getWriteId,
             R"DOC(return the StereoGroup's ID that will be exported.
Note that the ID only makes sense for AND/OR groups.)DOC")
        .def("SetWriteId", &StereoGroup::setWriteId, "id"_a,
             R"DOC(return the StereoGroup's ID that will be exported.
Note that the ID only makes sense for AND/OR groups.)DOC");

    m.def(
        "CreateStereoGroup", &createStereoGroup, "stereoGroupType"_a, "mol"_a,
        "atomIds"_a = nb::list(), "bondIds"_a = nb::list(), "readId"_a = 0,
        R"DOC(creates a StereoGroup associated with a molecule from a list of atom Ids)DOC",
        nb::rv_policy::take_ownership, nb::keep_alive<0, 2>());

    m.def(
        "ForwardStereoGroupIds", &RDKit::forwardStereoGroupIds, "mol"_a,
        R"DOC(Forward the original Stereo Group IDs when exporting the Mol.)DOC");
  }
};
}  // namespace RDKit

void wrap_stereogroup(nb::module_ &m) { RDKit::stereogroup_wrap::wrap(m); }
