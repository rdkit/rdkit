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

class MetalDisconnectorWrap {
 public:
  MetalDisconnectorWrap(python::object options = python::object()) {
    if (options.is_none()) {
      md_.reset(new RDKit::MolStandardize::MetalDisconnector());
    } else {
      RDKit::MolStandardize::MetalDisconnectorOptions md_opts;
      md_opts.splitGrignards =
          python::extract<bool>(options.attr("splitGrignards"));
      md_opts.splitAromaticC =
          python::extract<bool>(options.attr("splitAromaticC"));
      md_opts.adjustCharges =
          python::extract<bool>(options.attr("adjustCharges"));
      md_opts.removeHapticDummies =
          python::extract<bool>(options.attr("removeHapticDummies"));
      md_.reset(new RDKit::MolStandardize::MetalDisconnector(md_opts));
    }
  }

  RDKit::ROMol *getMetalNof() { return md_->getMetalNof(); }
  RDKit::ROMol *getMetalNon() { return md_->getMetalNon(); }
  void setMetalNof(const RDKit::ROMol &mol) { md_->setMetalNof(mol); }
  void setMetalNon(const RDKit::ROMol &mol) { md_->setMetalNon(mol); }

  RDKit::ROMol *disconnect(const RDKit::ROMol &mol) {
    return md_->disconnect(mol);
  }
  void disconnectInPlace(RDKit::ROMol &mol) {
    return md_->disconnectInPlace(static_cast<RDKit::RWMol &>(mol));
  }

 private:
  std::unique_ptr<RDKit::MolStandardize::MetalDisconnector> md_;
};

std::string getMetalNofHelper(MetalDisconnectorWrap &self) {
  return RDKit::MolToSmarts(*(self.getMetalNof()));
}

std::string getMetalNonHelper(MetalDisconnectorWrap &self) {
  return RDKit::MolToSmarts(*(self.getMetalNon()));
}

void setMetalNonHelper(MetalDisconnectorWrap &self, const RDKit::ROMol &mol) {
  self.setMetalNon(mol);
}

void setMetalNofHelper(MetalDisconnectorWrap &self, const RDKit::ROMol &mol) {
  self.setMetalNof(mol);
}
}  // namespace

struct metal_wrapper {
  static void wrap() {
    python::scope().attr("__doc__") =
        "Module containing functions for molecular standardization";

    std::string docString = "Metal Disconnector Options";
    python::class_<RDKit::MolStandardize::MetalDisconnectorOptions>(
        "MetalDisconnectorOptions", docString.c_str(),
        python::init<>(python::args("self")))
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
    python::class_<MetalDisconnectorWrap, boost::noncopyable>(
        "MetalDisconnector", docString.c_str(),
        python::init<python::optional<python::object>>(
            (python::arg("self"), python::arg("options") = python::object())))
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
        .def("Disconnect", &MetalDisconnectorWrap::disconnect,
             (python::arg("self"), python::arg("mol")),
             "performs the disconnection",
             python::return_value_policy<python::manage_new_object>())
        .def("DisconnectInPlace", &MetalDisconnectorWrap::disconnectInPlace,
             (python::arg("self"), python::arg("mol")),
             "performs the disconnection, modifies the input molecule");
  }
};

void wrap_metal() { metal_wrapper::wrap(); }
