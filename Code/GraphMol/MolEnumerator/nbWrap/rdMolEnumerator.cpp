//
//  Copyright (C) 2020-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <GraphMol/MolEnumerator/MolEnumerator.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {

enum class EnumeratorTypes { LinkNode, PositionVariation, RepeatUnit };

std::shared_ptr<MolEnumerator::MolEnumeratorOp> opFromName(
    EnumeratorTypes typ) {
  switch (typ) {
    case EnumeratorTypes::LinkNode:
      return std::make_shared<MolEnumerator::LinkNodeOp>();
    case EnumeratorTypes::PositionVariation:
      return std::make_shared<MolEnumerator::PositionVariationOp>();
    case EnumeratorTypes::RepeatUnit:
      return std::make_shared<MolEnumerator::RepeatUnitOp>();
    default:
      throw std::invalid_argument("unrecognized operator type");
  }
}

}  // namespace

NB_MODULE(rdMolEnumerator, m) {
  m.doc() = "Module containing classes and functions for enumerating molecules";

  nb::enum_<EnumeratorTypes>(m, "EnumeratorType")
      .value("LinkNode", EnumeratorTypes::LinkNode)
      .value("PositionVariation", EnumeratorTypes::PositionVariation)
      .value("RepeatUnit", EnumeratorTypes::RepeatUnit);

  nb::class_<MolEnumerator::MolEnumeratorParams>(m, "MolEnumeratorParams",
                                                  "Molecular enumerator parameters")
      .def(nb::init<>())
      .def("__init__",
           [](MolEnumerator::MolEnumeratorParams *self, EnumeratorTypes typ) {
             new (self) MolEnumerator::MolEnumeratorParams();
             self->dp_operation = opFromName(typ);
           },
           "typ"_a)
      .def_rw("sanitize", &MolEnumerator::MolEnumeratorParams::sanitize,
              "sanitize molecules after enumeration")
      .def_rw("maxToEnumerate",
              &MolEnumerator::MolEnumeratorParams::maxToEnumerate,
              "maximum number of molecules to enumerate")
      .def_rw("doRandom", &MolEnumerator::MolEnumeratorParams::doRandom,
              "do random enumeration (not yet implemented)")
      .def_rw("randomSeed", &MolEnumerator::MolEnumeratorParams::randomSeed,
              "seed for the random enumeration (not yet implemented)")
      .def(
          "SetEnumerationOperator",
          [](MolEnumerator::MolEnumeratorParams *self, EnumeratorTypes typ) {
            self->dp_operation = opFromName(typ);
          },
          "typ"_a, "set the operator to be used for enumeration");

  m.def(
      "Enumerate",
      [](const ROMol &mol, unsigned int maxPerOperation) {
        return MolEnumerator::enumerate(mol, maxPerOperation);
      },
      "mol"_a, "maxPerOperation"_a = 0u,
      R"DOC(do an enumeration and return a MolBundle.
If maxPerOperation is >0 that will be used as the maximum number of molecules which
can be returned by any given operation.
Limitations:
  - the current implementation does not support molecules which include both
    SRUs and LINKNODEs
  - Overlapping SRUs, i.e. where one monomer is contained within another, are
    not supported)DOC");

  m.def(
      "Enumerate",
      [](const ROMol &mol, const MolEnumerator::MolEnumeratorParams &ps) {
        return MolEnumerator::enumerate(mol, ps);
      },
      "mol"_a, "enumParams"_a,
      R"DOC(do an enumeration for the supplied parameter type and return a MolBundle
Limitations:
  - Overlapping SRUs, i.e. where one monomer is contained within another, are
    not supported)DOC");
}
