//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/Wrap.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/MolStandardize/Metal.h>

namespace python = boost::python;

namespace {
RDKit::ROMol *disconnect(RDKit::MolStandardize::MetalDisconnector &self,
                         RDKit::ROMol &mol) {
  RDKit::ROMol *nm = self.disconnect(mol);
  return nm;
}

std::string getMetalNofHelper(RDKit::MolStandardize::MetalDisconnector &self) {
  return RDKit::MolToSmarts(*(self.getMetalNof()));
}

std::string getMetalNonHelper(RDKit::MolStandardize::MetalDisconnector &self) {
  return RDKit::MolToSmarts(*(self.getMetalNon()));
}

void setMetalNonHelper(RDKit::MolStandardize::MetalDisconnector &self,
                       const RDKit::ROMol &mol) {
  self.setMetalNon(mol);
}

void setMetalNofHelper(RDKit::MolStandardize::MetalDisconnector &self,
                       const RDKit::ROMol &mol) {
  self.setMetalNof(mol);
}

}  // namespace

struct metal_wrapper {
  static void wrap() {
    python::scope().attr("__doc__") =
        "Module containing functions for molecular standardization";

    std::string docString = "Metal Disconnector Options";
    python::class_<RDKit::MolStandardize::MetalDisconnectorOptions,
                   boost::noncopyable>("MetalDisconnectorOptions",
                                       docString.c_str(), python::init<>())
        .def_readwrite(
            "splitGrignards",
            &RDKit::MolStandardize::MetalDisconnectorOptions::splitGrignards,
            "Whether to split Grignard-type complexes. Default false.")
        .def_readwrite(
            "splitAromaticC",
            &RDKit::MolStandardize::MetalDisconnectorOptions::splitAromaticC,
            "Whether to split metal-aromatic C bonds.  Default false.")
        .def_readwrite(
            "adjustCharges",
            &RDKit::MolStandardize::MetalDisconnectorOptions::adjustCharges,
            "Whether to adjust charges on ligand atoms.  Default true.")
        .def_readwrite("removeHapticDummies",
                       &RDKit::MolStandardize::MetalDisconnectorOptions::
                           removeHapticDummies,
                       "Whether to remove the dummy atoms representing haptic"
                       " bonds.  Such dummies are bonded to the metal with a"
                       " bond that has the MolFileBondEndPts prop set."
                       "  Default false.");

    docString =
        "a class to disconnect metals that are defined as covalently bonded to"
        " non-metals";
    python::class_<RDKit::MolStandardize::MetalDisconnector,
                   boost::noncopyable>(
        "MetalDisconnector", docString.c_str(),
        python::init<python::optional<
            const RDKit::MolStandardize::MetalDisconnectorOptions &>>(
            (python::arg("options") = python::object())))
        .add_property("MetalNof", &getMetalNofHelper,
                      "SMARTS defining the metals to disconnect if attached to "
                      "Nitrogen, Oxygen or Fluorine")
        .add_property(
            "MetalNon", &getMetalNonHelper,
            "SMARTS defining the metals to disconnect other inorganic elements")
        .def(
            "SetMetalNon", &setMetalNonHelper,
            (python::arg("self"), python::arg("mol")),
            "Set the query molecule defining the metals to disconnect from other"
            " inorganic elements.")
        .def(
            "SetMetalNof", &setMetalNofHelper,
            (python::arg("self"), python::arg("mol")),
            "Set the query molecule defining the metals to disconnect if attached"
            " to Nitrogen, Oxygen or Fluorine.")
        .def("Disconnect", &disconnect,
             (python::arg("self"), python::arg("mol")), docString.c_str(),
             python::return_value_policy<python::manage_new_object>());
  }
};

void wrap_metal() { metal_wrapper::wrap(); }
