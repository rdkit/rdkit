//
//  Copyright (C) 2020-2021 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/python.h>
#include <GraphMol/MolEnumerator/MolEnumerator.h>

namespace python = boost::python;
using namespace RDKit;
namespace {

enum class EnumeratorTypes { LinkNode, PositionVariation, RepeatUnit };

std::shared_ptr<MolEnumerator::MolEnumeratorOp> opFromName(
    EnumeratorTypes typ) {
  std::shared_ptr<MolEnumerator::MolEnumeratorOp> res;
  switch (typ) {
    case EnumeratorTypes::LinkNode:
      res.reset(new MolEnumerator::LinkNodeOp());
      break;
    case EnumeratorTypes::PositionVariation:
      res.reset(new MolEnumerator::PositionVariationOp());
      break;
    case EnumeratorTypes::RepeatUnit:
      res.reset(new MolEnumerator::RepeatUnitOp());
      break;
    default:
      throw ValueErrorException("unrecognized operator type");
  }
  return res;
}
MolEnumerator::MolEnumeratorParams *createParamsFromName(EnumeratorTypes typ) {
  auto res = new MolEnumerator::MolEnumeratorParams();
  res->dp_operation = opFromName(typ);
  return res;
}
void setEnumerationHelper(MolEnumerator::MolEnumeratorParams *self,
                          EnumeratorTypes typ) {
  self->dp_operation = opFromName(typ);
}
MolBundle *enumerateHelper1(const ROMol &mol, unsigned int maxPerOperation) {
  auto res = MolEnumerator::enumerate(mol, maxPerOperation);
  return new MolBundle(res);
}
MolBundle *enumerateHelper2(const ROMol &mol,
                            const MolEnumerator::MolEnumeratorParams &ps) {
  auto res = MolEnumerator::enumerate(mol, ps);
  return new MolBundle(res);
}
}  // namespace
BOOST_PYTHON_MODULE(rdMolEnumerator) {
  python::scope().attr("__doc__") =
      "Module containing classes and functions for enumerating molecules";
  python::enum_<EnumeratorTypes>("EnumeratorType")
      .value("LinkNode", EnumeratorTypes::LinkNode)
      .value("PositionVariation", EnumeratorTypes::PositionVariation)
      .value("RepeatUnit", EnumeratorTypes::RepeatUnit);

  python::class_<MolEnumerator::MolEnumeratorParams>(
      "MolEnumeratorParams", "Molecular enumerator parameters",
      python::init<>())
      .def("__init__", python::make_constructor(createParamsFromName))
      .def_readwrite("sanitize", &MolEnumerator::MolEnumeratorParams::sanitize,
                     "sanitize molecules after enumeration")
      .def_readwrite("maxToEnumerate",
                     &MolEnumerator::MolEnumeratorParams::maxToEnumerate,
                     "maximum number of molecules to enumerate")
      .def_readwrite("doRandom", &MolEnumerator::MolEnumeratorParams::doRandom,
                     "do random enumeration (not yet implemented")
      .def_readwrite("randomSeed",
                     &MolEnumerator::MolEnumeratorParams::randomSeed,
                     "seed for the random enumeration (not yet implemented")
      .def("SetEnumerationOperator", &setEnumerationHelper,
           "set the operator to be used for enumeration");
  python::def("Enumerate", &enumerateHelper1,
              (python::arg("mol"), python::arg("maxPerOperation") = 0),
              python::return_value_policy<python::manage_new_object>(),
              R"DOC(do an enumeration and return a MolBundle.
  If maxPerOperation is >0 that will be used as the maximum number of molecules which 
    can be returned by any given operation.
Limitations:
  - the current implementation does not support molecules which include both
    SRUs and LINKNODEs
  - Overlapping SRUs, i.e. where one monomer is contained within another, are
    not supported)DOC");
  python::def(
      "Enumerate", &enumerateHelper2,
      (python::arg("mol"), python::arg("enumParams")),
      python::return_value_policy<python::manage_new_object>(),
      R"DOC(do an enumeration for the supplied parameter type and return a MolBundle
Limitations:
  - Overlapping SRUs, i.e. where one monomer is contained within another, are
    not supported)DOC");
}
