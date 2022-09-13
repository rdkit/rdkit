//
//  Copyright (C) 2019 Greg Landrum
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
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SubstanceGroup.h>
#include <RDBoost/Wrap.h>
#include "props.hpp"

namespace python = boost::python;

namespace RDKit {

namespace {

SubstanceGroup *getMolSubstanceGroupWithIdx(ROMol &mol, unsigned int idx) {
  auto &sgs = getSubstanceGroups(mol);
  if (idx >= sgs.size()) {
    throw_index_error(idx);
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

void addBracketHelper(SubstanceGroup &self, python::object pts) {
  unsigned int sz = python::extract<unsigned int>(pts.attr("__len__")());
  if (sz != 2 && sz != 3) {
    throw_value_error("pts object have a length of 2 or 3");
  }

  SubstanceGroup::Bracket bkt;
  python::stl_input_iterator<RDGeom::Point3D> beg(pts);
  for (unsigned int i = 0; i < sz; ++i) {
    bkt[i] = *beg;
    ++beg;
  }
  self.addBracket(bkt);
}

python::tuple getCStatesHelper(const SubstanceGroup &self) {
  python::list res;
  for (const auto &cs : self.getCStates()) {
    res.append(cs);
  }
  return python::tuple(res);
}

python::tuple getBracketsHelper(const SubstanceGroup &self) {
  python::list res;
  for (const auto &brk : self.getBrackets()) {
    res.append(python::make_tuple(brk[0], brk[1], brk[2]));
  }
  return python::tuple(res);
}

python::tuple getAttachPointsHelper(const SubstanceGroup &self) {
  python::list res;
  for (const auto &ap : self.getAttachPoints()) {
    res.append(ap);
  }
  return python::tuple(res);
}

void SetAtomsHelper(SubstanceGroup &self, const python::object &iterable) {
  std::vector<unsigned int> atoms;
  pythonObjectToVect(iterable, atoms);
  self.setAtoms(atoms);
}

void SetParentAtomsHelper(SubstanceGroup &self,
                          const python::object &iterable) {
  std::vector<unsigned int> patoms;
  pythonObjectToVect(iterable, patoms);
  self.setParentAtoms(patoms);
}

void SetBondsHelper(SubstanceGroup &self, const python::object &iterable) {
  std::vector<unsigned int> bonds;
  pythonObjectToVect(iterable, bonds);
  self.setBonds(bonds);
}

}  // namespace

std::string sGroupClassDoc =
    "A collection of atoms and bonds with associated properties\n";

struct sgroup_wrap {
  static void wrap() {
    RegisterVectorConverter<SubstanceGroup>("SubstanceGroup_VECT");

    python::class_<SubstanceGroup::CState,
                   boost::shared_ptr<SubstanceGroup::CState>>(
        "SubstanceGroupCState", "CSTATE for a SubstanceGroup", python::init<>())
        .def_readonly("bondIdx", &SubstanceGroup::CState::bondIdx)
        .def_readonly("vector", &SubstanceGroup::CState::vector);

    python::class_<SubstanceGroup::AttachPoint,
                   boost::shared_ptr<SubstanceGroup::AttachPoint>>(
        "SubstanceGroupAttach", "AttachPoint for a SubstanceGroup",
        python::init<>())
        .def_readonly("aIdx", &SubstanceGroup::AttachPoint::aIdx,
                      "attachment index")
        .def_readonly("lvIdx", &SubstanceGroup::AttachPoint::lvIdx,
                      "leaving atom or index (0 for implied)")
        .def_readonly("id", &SubstanceGroup::AttachPoint::id, "attachment id");

    python::class_<SubstanceGroup, boost::shared_ptr<SubstanceGroup>>(
        "SubstanceGroup", sGroupClassDoc.c_str(), python::no_init)
        .def("GetOwningMol", &SubstanceGroup::getOwningMol,
             "returns the molecule owning this SubstanceGroup",
             python::return_internal_reference<>())
        .def("GetIndexInMol", &SubstanceGroup::getIndexInMol,
             "returns the index of this SubstanceGroup in the owning "
             "molecule's list.")
        .def("GetAtoms", &SubstanceGroup::getAtoms,
             "returns a list of the indices of the atoms in this "
             "SubstanceGroup",
             python::return_value_policy<python::copy_const_reference>())
        .def("GetParentAtoms", &SubstanceGroup::getParentAtoms,
             "returns a list of the indices of the parent atoms in this "
             "SubstanceGroup",
             python::return_value_policy<python::copy_const_reference>())
        .def("GetBonds", &SubstanceGroup::getBonds,
             "returns a list of the indices of the bonds in this "
             "SubstanceGroup",
             python::return_value_policy<python::copy_const_reference>())
        .def("SetAtoms", SetAtomsHelper,
             "Set the list of the indices of the atoms in this "
             "SubstanceGroup.\nNote that this does not update "
             "properties, CStates or Attachment Points.")
        .def("SetParentAtoms", SetParentAtomsHelper,
             "Set the list of the indices of the parent atoms in this "
             "SubstanceGroup.\nNote that this does not update "
             "properties, CStates or Attachment Points.")
        .def("SetBonds", SetBondsHelper,
             "Set the list of the indices of the bonds in this "
             "SubstanceGroup.\nNote that this does not update "
             "properties, CStates or Attachment Points.")
        .def("AddAtomWithIdx", &SubstanceGroup::addAtomWithIdx)
        .def("AddBondWithIdx", &SubstanceGroup::addBondWithIdx)
        .def("AddParentAtomWithIdx", &SubstanceGroup::addParentAtomWithIdx)
        .def("AddAtomWithBookmark", &SubstanceGroup::addAtomWithBookmark)
        .def("AddParentAtomWithBookmark",
             &SubstanceGroup::addParentAtomWithBookmark)
        .def("AddCState", &SubstanceGroup::addCState)
        .def("GetCStates", getCStatesHelper)
        .def("AddBondWithBookmark", &SubstanceGroup::addBondWithBookmark)
        .def("AddAttachPoint", &SubstanceGroup::addAttachPoint)
        .def("GetAttachPoints", getAttachPointsHelper)
        .def("AddBracket", addBracketHelper)
        .def("GetBrackets", getBracketsHelper)
        .def("ClearBrackets", &SubstanceGroup::clearBrackets)
        .def("ClearCStates", &SubstanceGroup::clearCStates)
        .def("ClearAttachPoints", &SubstanceGroup::clearAttachPoints)

        .def("SetProp",
             (void (RDProps::*)(const std::string &, std::string, bool) const) &
                 SubstanceGroup::setProp<std::string>,
             (python::arg("self"), python::arg("key"), python::arg("val"),
              python::arg("computed") = false),
             "sets the value of a particular property")
        .def("SetDoubleProp",
             (void (RDProps::*)(const std::string &, double, bool) const) &
                 SubstanceGroup::setProp<double>,
             (python::arg("self"), python::arg("key"), python::arg("val"),
              python::arg("computed") = false),
             "sets the value of a particular property")
        .def("SetIntProp",
             (void (RDProps::*)(const std::string &, int, bool) const) &
                 SubstanceGroup::setProp<int>,
             (python::arg("self"), python::arg("key"), python::arg("val"),
              python::arg("computed") = false),
             "sets the value of a particular property")
        .def(
            "SetUnsignedProp",
            (void (RDProps::*)(const std::string &, unsigned int, bool) const) &
                SubstanceGroup::setProp<unsigned int>,
            (python::arg("self"), python::arg("key"), python::arg("val"),
             python::arg("computed") = false),
            "sets the value of a particular property")
        .def("SetBoolProp",
             (void (RDProps::*)(const std::string &, bool, bool) const) &
                 SubstanceGroup::setProp<bool>,
             (python::arg("self"), python::arg("key"), python::arg("val"),
              python::arg("computed") = false),
             "sets the value of a particular property")
        .def("HasProp",
             (bool (RDProps::*)(const std::string &) const) &
                 SubstanceGroup::hasProp,
             "returns whether or not a particular property exists")
        .def("GetProp",
             (std::string(RDProps::*)(const std::string &) const) &
                 SubstanceGroup::getProp<std::string>,
             "returns the value of a particular property")
        .def("GetIntProp",
             (int (RDProps::*)(const std::string &) const) &
                 SubstanceGroup::getProp<int>,
             "returns the value of a particular property")
        .def("GetUnsignedProp",
             (unsigned int (RDProps::*)(const std::string &) const) &
                 SubstanceGroup::getProp<unsigned int>,
             "returns the value of a particular property")
        .def("GetDoubleProp",
             (double (RDProps::*)(const std::string &) const) &
                 SubstanceGroup::getProp<double>,
             "returns the value of a particular property")
        .def("GetBoolProp",
             (bool (RDProps::*)(const std::string &) const) &
                 SubstanceGroup::getProp<bool>,
             "returns the value of a particular property")
        .def(
            "GetUnsignedVectProp",
            (std::vector<unsigned int>(RDProps::*)(const std::string &) const) &
                SubstanceGroup::getProp<std::vector<unsigned int>>,
            "returns the value of a particular property")
        .def("GetStringVectProp",
             (std::vector<std::string>(RDProps::*)(const std::string &) const) &
                 SubstanceGroup::getProp<std::vector<std::string>>,
             "returns the value of a particular property")
        .def("GetPropNames", &SubstanceGroup::getPropList,
             (python::arg("self"), python::arg("includePrivate") = false,
              python::arg("includeComputed") = false),
             "Returns a list of the properties set on the "
             "SubstanceGroup.\n\n")
        .def("GetPropsAsDict", GetPropsAsDict<SubstanceGroup>,
             (python::arg("self"), python::arg("includePrivate") = true,
              python::arg("includeComputed") = true),
             "Returns a dictionary of the properties set on the "
             "SubstanceGroup.\n"
             " n.b. some properties cannot be converted to python types.\n")
        .def("ClearProp",
             (void (RDProps::*)(const std::string &) const) &
                 SubstanceGroup::clearProp,
             "Removes a particular property (does nothing if not set).\n\n");

    python::def("GetMolSubstanceGroups", &getMolSubstanceGroups,
                "returns a copy of the molecule's SubstanceGroups (if any)",
                python::with_custodian_and_ward_postcall<0, 1>());
    python::def("GetMolSubstanceGroupWithIdx", &getMolSubstanceGroupWithIdx,
                "returns a particular SubstanceGroup from the molecule",
                python::return_internal_reference<
                    1, python::with_custodian_and_ward_postcall<0, 1>>());
    python::def("ClearMolSubstanceGroups", &clearMolSubstanceGroups,
                "removes all SubstanceGroups from a molecule (if any)");
    python::def("CreateMolSubstanceGroup", &createMolSubstanceGroup,
                (python::arg("mol"), python::arg("type")),
                "creates a new SubstanceGroup associated with a molecule, "
                "returns the new SubstanceGroup",
                python::return_value_policy<
                    python::reference_existing_object,
                    python::with_custodian_and_ward_postcall<0, 1>>());
    python::def(
        "CreateMolDataSubstanceGroup", &createMolDataSubstanceGroup,
        (python::arg("mol"), python::arg("fieldName"), python::arg("value")),
        "creates a new DATA SubstanceGroup associated with a molecule, "
        "returns the new SubstanceGroup",
        python::return_value_policy<
            python::reference_existing_object,
            python::with_custodian_and_ward_postcall<0, 1>>());
    python::def("AddMolSubstanceGroup", &addMolSubstanceGroup,
                (python::arg("mol"), python::arg("sgroup")),
                "adds a copy of a SubstanceGroup to a molecule, returns the "
                "new SubstanceGroup",
                python::return_value_policy<
                    python::reference_existing_object,
                    python::with_custodian_and_ward_postcall<0, 1>>());
  }
};
}  // namespace RDKit

void wrap_sgroup() { RDKit::sgroup_wrap::wrap(); }
