//
//  Copyright (C) 2018-2026 Susan H. Leung and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>

#include <optional>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/MolStandardize/Metal.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

class MetalDisconnectorWrap {
 public:
  MetalDisconnectorWrap(
      std::optional<RDKit::MolStandardize::MetalDisconnectorOptions *> options =
          std::nullopt) {
    if (!options.has_value() || !options.value()) {
      md_.reset(new RDKit::MolStandardize::MetalDisconnector());
    } else {
      md_.reset(new RDKit::MolStandardize::MetalDisconnector(*options.value()));
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

void wrap_metal(nb::module_ &m) {
  nb::class_<RDKit::MolStandardize::MetalDisconnectorOptions>(
      m, "MetalDisconnectorOptions", "Metal Disconnector Options")
      .def(nb::init<>())
      .def_rw("splitGrignards",
              &RDKit::MolStandardize::MetalDisconnectorOptions::splitGrignards,
              "Whether to split Grignard-type complexes. Default false.")
      .def_rw("splitAromaticC",
              &RDKit::MolStandardize::MetalDisconnectorOptions::splitAromaticC,
              "Whether to split metal-aromatic C bonds.  Default false.")
      .def_rw("adjustCharges",
              &RDKit::MolStandardize::MetalDisconnectorOptions::adjustCharges,
              "Whether to adjust charges on ligand atoms.  Default true.")
      .def_rw(
          "removeHapticDummies",
          &RDKit::MolStandardize::MetalDisconnectorOptions::removeHapticDummies,
          "Whether to remove the dummy atoms representing haptic"
          " bonds.  Such dummies are bonded to the metal with a"
          " bond that has the MolFileBondEndPts prop set."
          "  Default false.");

  nb::class_<MetalDisconnectorWrap>(
      m, "MetalDisconnector",
      "a class to disconnect metals that are defined as covalently bonded to"
      " non-metals")
      .def(nb::init<std::optional<
               RDKit::MolStandardize::MetalDisconnectorOptions *>>(),
           "options"_a = nb::none())
      .def_prop_ro("MetalNof", &getMetalNofHelper,
                   "SMARTS defining the metals to disconnect if attached to "
                   "Nitrogen, Oxygen or Fluorine")
      .def_prop_ro(
          "MetalNon", &getMetalNonHelper,
          "SMARTS defining the metals to disconnect other inorganic elements")
      .def("SetMetalNon", &setMetalNonHelper, "mol"_a,
           "Set the query molecule defining the metals to disconnect from other"
           " inorganic elements.")
      .def("SetMetalNof", &setMetalNofHelper, "mol"_a,
           "Set the query molecule defining the metals to disconnect if "
           "attached to Nitrogen, Oxygen or Fluorine.")
      .def("Disconnect", &MetalDisconnectorWrap::disconnect, "mol"_a,
           "performs the disconnection", nb::rv_policy::take_ownership)
      .def("DisconnectInPlace", &MetalDisconnectorWrap::disconnectInPlace,
           "mol"_a, "performs the disconnection, modifies the input molecule");
}
