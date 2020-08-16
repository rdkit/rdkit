//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
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

enum class EnumeratorTypes { LinkNode, PositionVariation };

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
MolBundle *enumerateHelper(const ROMol &mol,
                           const MolEnumerator::MolEnumeratorParams &ps) {
  auto res = MolEnumerator::enumerate(mol, ps);
  return new MolBundle(res);
}
}  // namespace
BOOST_PYTHON_MODULE(rdMolEnumerator) {
  python::scope().attr("__doc__") =
      "Module containing classes and functions for enumerating query molecules";
  python::enum_<EnumeratorTypes>("EnumeratorType")
      .value("LinkNode", EnumeratorTypes::LinkNode)
      .value("PositionVariation", EnumeratorTypes::PositionVariation);

  python::class_<MolEnumerator::MolEnumeratorParams>(
      "MolEnumeratorParams", "Molecular enumerator parameters",
      python::init<>())
      .def("__init__", python::make_constructor(createParamsFromName))
      .def_readwrite("sanitize", &MolEnumerator::MolEnumeratorParams::sanitize,
                     "sanitize molecules after enumeration")
      .def_readwrite("maxToEnumerate",
                     &MolEnumerator::MolEnumeratorParams::maxToEnumerate,
                     "maximum number of molecules to enuemrate")
      .def_readwrite("doRandom", &MolEnumerator::MolEnumeratorParams::doRandom,
                     "do random enumeration (not yet implemented")
      .def_readwrite("randomSeed",
                     &MolEnumerator::MolEnumeratorParams::randomSeed,
                     "seed for the random enumeration (not yet implemented")
      .def("SetEnumerationOperator", &setEnumerationHelper,
           "set the operator to be used for enumeration");
  python::def("Enumerate", &enumerateHelper,
              python::return_value_policy<python::manage_new_object>(),
              "do an enumeration and return a MolBundle");
}
