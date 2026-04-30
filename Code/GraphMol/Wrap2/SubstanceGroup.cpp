//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <string>
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>

// ours
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SubstanceGroup.h>
#include <RDBoost/Wrap_nb.h>
#include "props.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {

namespace {

SubstanceGroup *getMolSubstanceGroupWithIdx(ROMol &mol, unsigned int idx) {
  auto &sgs = getSubstanceGroups(mol);
  if (idx >= sgs.size()) {
    auto msg = "SubstanceGroup index out of range: " + std::to_string(idx);
    throw nb::index_error(msg.c_str());
  }
  return &(sgs[idx]);
}

std::vector<SubstanceGroup> getMolSubstanceGroups(ROMol &mol) {
  return getSubstanceGroups(mol);
}
void clearMolSubstanceGroups(ROMol &mol) {
  std::vector<SubstanceGroup> &sgs = getSubstanceGroups(mol);
  sgs.clear();
}

SubstanceGroup *createMolSubstanceGroup(ROMol &mol, std::string type) {
  SubstanceGroup sg(&mol, type);
  addSubstanceGroup(mol, sg);
  return &(getSubstanceGroups(mol).back());
}

SubstanceGroup *createMolDataSubstanceGroup(ROMol &mol, std::string fieldName,
                                            std::string value) {
  SubstanceGroup sg(&mol, "DAT");
  sg.setProp("FIELDNAME", fieldName);
  STR_VECT dataFields{value};
  sg.setProp("DATAFIELDS", dataFields);
  addSubstanceGroup(mol, sg);
  return &(getSubstanceGroups(mol).back());
}

SubstanceGroup *addMolSubstanceGroup(ROMol &mol, const SubstanceGroup &sgroup) {
  addSubstanceGroup(mol, sgroup);
  return &(getSubstanceGroups(mol).back());
}

void addBracketHelper(SubstanceGroup &self, const nb::object &pts) {
  unsigned int sz = static_cast<unsigned int>(nb::len(pts));
  if (sz != 2 && sz != 3) {
    throw ValueErrorException("pts object have a length of 2 or 3");
  }

  SubstanceGroup::Bracket bkt;
  auto ptVec = pythonObjectToVect<RDGeom::Point3D>(pts);
  if (!ptVec || ptVec->size() != sz) {
    throw ValueErrorException(
        "could not interpret pts as a sequence of 3D points");
  }
  for (unsigned int i = 0; i < sz; ++i) {
    bkt[i] = (*ptVec)[i];
  }
  self.addBracket(bkt);
}

nb::tuple getCStatesHelper(const SubstanceGroup &self) {
  nb::list res;
  for (const auto &cs : self.getCStates()) {
    res.append(cs);
  }
  return nb::tuple(res);
}

nb::tuple getBracketsHelper(const SubstanceGroup &self) {
  nb::list res;
  for (const auto &brk : self.getBrackets()) {
    res.append(nb::make_tuple(brk[0], brk[1], brk[2]));
  }
  return nb::tuple(res);
}

nb::tuple getAttachPointsHelper(const SubstanceGroup &self) {
  nb::list res;
  for (const auto &ap : self.getAttachPoints()) {
    res.append(ap);
  }
  return nb::tuple(res);
}

void SetAtomsHelper(SubstanceGroup &self, const nb::object &iterable) {
  std::vector<unsigned int> atoms;
  pythonObjectToVect(iterable, atoms);
  self.setAtoms(atoms);
}

void SetParentAtomsHelper(SubstanceGroup &self, const nb::object &iterable) {
  std::vector<unsigned int> patoms;
  pythonObjectToVect(iterable, patoms);
  self.setParentAtoms(patoms);
}

void SetBondsHelper(SubstanceGroup &self, const nb::object &iterable) {
  std::vector<unsigned int> bonds;
  pythonObjectToVect(iterable, bonds);
  self.setBonds(bonds);
}

}  // namespace

std::string sGroupClassDoc =
    R"DOC(A collection of atoms and bonds with associated properties
)DOC";

struct sgroup_wrap {
  static void wrap(nb::module_ &m) {
    nb::class_<SubstanceGroup::CState>(m, "SubstanceGroupCState",
                                       R"DOC(CSTATE for a SubstanceGroup)DOC")
        .def(nb::init<>())
        .def_ro("bondIdx", &SubstanceGroup::CState::bondIdx)
        .def_ro("vector", &SubstanceGroup::CState::vector);

    nb::class_<SubstanceGroup::AttachPoint>(
        m, "SubstanceGroupAttach", R"DOC(AttachPoint for a SubstanceGroup)DOC")
        .def(nb::init<>())
        .def_ro("aIdx", &SubstanceGroup::AttachPoint::aIdx,
                R"DOC(attachment index)DOC")
        .def_ro("lvIdx", &SubstanceGroup::AttachPoint::lvIdx,
                R"DOC(leaving atom or index (0 for implied))DOC")
        .def_ro("id", &SubstanceGroup::AttachPoint::id,
                R"DOC(attachment id)DOC");

    nb::class_<SubstanceGroup>(m, "SubstanceGroup", sGroupClassDoc.c_str())
        .def("GetOwningMol", &SubstanceGroup::getOwningMol,
             nb::rv_policy::reference_internal,
             R"DOC(returns the molecule owning this SubstanceGroup)DOC")
        .def(
            "GetIndexInMol", &SubstanceGroup::getIndexInMol,
            R"DOC(returns the index of this SubstanceGroup in the owning molecule's list.)DOC")
        .def(
            "GetAtoms", &SubstanceGroup::getAtoms,
            nb::rv_policy::reference_internal,
            R"DOC(returns a list of the indices of the atoms in this SubstanceGroup)DOC")
        .def(
            "GetParentAtoms", &SubstanceGroup::getParentAtoms,
            nb::rv_policy::reference_internal,
            R"DOC(returns a list of the indices of the parent atoms in this SubstanceGroup)DOC")
        .def(
            "GetBonds", &SubstanceGroup::getBonds,
            nb::rv_policy::reference_internal,
            R"DOC(returns a list of the indices of the bonds in this SubstanceGroup)DOC")
        .def(
            "SetAtoms", SetAtomsHelper, "iterable"_a,
            R"DOC(Set the list of the indices of the atoms in this SubstanceGroup.
Note that this does not update properties, CStates or Attachment Points.)DOC")
        .def(
            "SetParentAtoms", SetParentAtomsHelper, "iterable"_a,
            R"DOC(Set the list of the indices of the parent atoms in this SubstanceGroup.
Note that this does not update properties, CStates or Attachment Points.)DOC")
        .def(
            "SetBonds", SetBondsHelper, "iterable"_a,
            R"DOC(Set the list of the indices of the bonds in this SubstanceGroup.
Note that this does not update properties, CStates or Attachment Points.)DOC")
        .def("AddAtomWithIdx", &SubstanceGroup::addAtomWithIdx, "idx"_a)
        .def("AddBondWithIdx", &SubstanceGroup::addBondWithIdx, "idx"_a)
        .def("AddParentAtomWithIdx", &SubstanceGroup::addParentAtomWithIdx,
             "idx"_a)
        .def("AddAtomWithBookmark", &SubstanceGroup::addAtomWithBookmark,
             "mark"_a)
        .def("AddParentAtomWithBookmark",
             &SubstanceGroup::addParentAtomWithBookmark, "mark"_a)
        .def("AddCState", &SubstanceGroup::addCState, "bondIdx"_a, "vector"_a)
        .def("GetCStates", getCStatesHelper)
        .def("AddBondWithBookmark", &SubstanceGroup::addBondWithBookmark,
             "mark"_a)
        .def("AddAttachPoint", &SubstanceGroup::addAttachPoint, "aIdx"_a,
             "lvIdx"_a, "idStr"_a)
        .def("GetAttachPoints", getAttachPointsHelper)
        .def("AddBracket", addBracketHelper, "pts"_a)
        .def("GetBrackets", getBracketsHelper)
        .def("ClearBrackets", &SubstanceGroup::clearBrackets,
             R"DOC(Clear bracket definitions.)DOC")
        .def("ClearCStates", &SubstanceGroup::clearCStates,
             R"DOC(Clear CSTATE entries.)DOC")
        .def("ClearAttachPoints", &SubstanceGroup::clearAttachPoints,
             R"DOC(Clear attachment points.)DOC")

        .def("SetProp",
             (void (RDProps::*)(const std::string_view, std::string, bool)
                  const) &
                 SubstanceGroup::setProp<std::string>,
             "key"_a, "val"_a, "computed"_a = false,
             R"DOC(sets the value of a particular property)DOC")
        .def("SetDoubleProp",
             (void (RDProps::*)(const std::string_view, double, bool) const) &
                 SubstanceGroup::setProp<double>,
             "key"_a, "val"_a, "computed"_a = false,
             R"DOC(sets the value of a particular property)DOC")
        .def("SetIntProp",
             (void (RDProps::*)(const std::string_view, int, bool) const) &
                 SubstanceGroup::setProp<int>,
             "key"_a, "val"_a, "computed"_a = false,
             R"DOC(sets the value of a particular property)DOC")
        .def("SetUnsignedProp",
             (void (RDProps::*)(const std::string_view, unsigned int, bool)
                  const) &
                 SubstanceGroup::setProp<unsigned int>,
             "key"_a, "val"_a, "computed"_a = false,
             R"DOC(sets the value of a particular property)DOC")
        .def("SetBoolProp",
             (void (RDProps::*)(const std::string_view, bool, bool) const) &
                 SubstanceGroup::setProp<bool>,
             "key"_a, "val"_a, "computed"_a = false,
             R"DOC(sets the value of a particular property)DOC")
        .def("HasProp",
             (bool (RDProps::*)(const std::string_view) const) &
                 SubstanceGroup::hasProp,
             "key"_a,
             R"DOC(returns whether or not a particular property exists)DOC")
        .def("GetProp", GetPyProp<SubstanceGroup>, "key"_a,
             "autoConvert"_a = false,
             R"DOC(Returns the value of the property.

  ARGUMENTS:
    - key: the name of the property to return (a string).

    - autoConvert: if True attempt to convert the property into a python object

  RETURNS: a string

  NOTE:
    - If the property has not been set, a KeyError exception will be raised.
)DOC")
        .def("GetIntProp",
             (int (RDProps::*)(const std::string_view) const) &
                 SubstanceGroup::getProp<int>,
             "key"_a, R"DOC(returns the value of a particular property)DOC")
        .def("GetUnsignedProp",
             (unsigned int (RDProps::*)(const std::string_view) const) &
                 SubstanceGroup::getProp<unsigned int>,
             "key"_a, R"DOC(returns the value of a particular property)DOC")
        .def("GetDoubleProp",
             (double (RDProps::*)(const std::string_view) const) &
                 SubstanceGroup::getProp<double>,
             "key"_a, R"DOC(returns the value of a particular property)DOC")
        .def("GetBoolProp",
             (bool (RDProps::*)(const std::string_view) const) &
                 SubstanceGroup::getProp<bool>,
             "key"_a, R"DOC(returns the value of a particular property)DOC")
        .def("GetUnsignedVectProp",
             (std::vector<unsigned int> (RDProps::*)(const std::string_view)
                  const) &
                 SubstanceGroup::getProp<std::vector<unsigned int>>,
             "key"_a, R"DOC(returns the value of a particular property)DOC")
        .def("GetStringVectProp",
             (std::vector<std::string> (RDProps::*)(const std::string_view)
                  const) &
                 SubstanceGroup::getProp<std::vector<std::string>>,
             "key"_a, R"DOC(returns the value of a particular property)DOC")
        .def("GetPropNames", &SubstanceGroup::getPropList,
             "includePrivate"_a = false, "includeComputed"_a = false,
             R"DOC(Returns a list of the properties set on the SubstanceGroup.

)DOC")
        .def(
            "GetPropsAsDict", GetPropsAsDict<SubstanceGroup>,
            "includePrivate"_a = true, "includeComputed"_a = true,
            "autoConvertStrings"_a = true,
            R"DOC(Returns a dictionary of the properties set on the SubstanceGroup.
 n.b. some properties cannot be converted to python types.
)DOC")
        .def("ClearProp",
             (void (RDProps::*)(const std::string_view) const) &
                 SubstanceGroup::clearProp,
             "key"_a,
             R"DOC(Removes a particular property (does nothing if not set).

)DOC");

    m.def("GetMolSubstanceGroups", &getMolSubstanceGroups, "mol"_a,
          R"DOC(returns a copy of the molecule's SubstanceGroups (if any))DOC",
          nb::keep_alive<0, 1>());
    m.def("GetMolSubstanceGroupWithIdx", &getMolSubstanceGroupWithIdx, "mol"_a,
          "idx"_a,
          R"DOC(returns a particular SubstanceGroup from the molecule)DOC",
          nb::rv_policy::reference_internal, nb::keep_alive<0, 1>());
    m.def("ClearMolSubstanceGroups", &clearMolSubstanceGroups, "mol"_a,
          R"DOC(removes all SubstanceGroups from a molecule (if any))DOC");
    m.def(
        "CreateMolSubstanceGroup", &createMolSubstanceGroup, "mol"_a, "type"_a,
        R"DOC(creates a new SubstanceGroup associated with a molecule, returns the new SubstanceGroup)DOC",
        nb::rv_policy::reference, nb::keep_alive<0, 1>());
    m.def(
        "CreateMolDataSubstanceGroup", &createMolDataSubstanceGroup, "mol"_a,
        "fieldName"_a, "value"_a,
        R"DOC(creates a new DATA SubstanceGroup associated with a molecule, returns the new SubstanceGroup)DOC",
        nb::rv_policy::reference, nb::keep_alive<0, 1>());
    m.def(
        "AddMolSubstanceGroup", &addMolSubstanceGroup, "mol"_a, "sgroup"_a,
        R"DOC(adds a copy of a SubstanceGroup to a molecule, returns the new SubstanceGroup)DOC",
        nb::rv_policy::reference, nb::keep_alive<0, 1>());
  }
};
}  // namespace RDKit

void wrap_sgroup(nb::module_ &m) { RDKit::sgroup_wrap::wrap(m); }
