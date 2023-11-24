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

enum class EnumeratorTypes {
  LinkNode,
  PositionVariation,
  RepeatUnit,
  StereoIsomer
};

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
    case EnumeratorTypes::StereoIsomer:
      res = MolEnumerator::StereoIsomerOp::createOp();
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
MolBundle *enumerateHelper1(const ROMol &mol, unsigned int maxPerOperation,
                            bool enumerateStereo) {
  auto res = MolEnumerator::enumerate(mol, maxPerOperation, enumerateStereo);
  return new MolBundle(res);
}
MolBundle *enumerateHelper2(const ROMol &mol,
                            const MolEnumerator::MolEnumeratorParams &ps) {
  auto res = MolEnumerator::enumerate(mol, ps);
  return new MolBundle(res);
}

MolBundle *enumerate_stereoisomers_helper(
    const ROMol &mol, const MolEnumerator::StereoEnumerationOptions &options) {
  return new MolBundle(MolEnumerator::enumerate_stereoisomers(mol, options));
}

}  // namespace
BOOST_PYTHON_MODULE(rdMolEnumerator) {
  python::scope().attr("__doc__") =
      "Module containing classes and functions for enumerating molecules";
  python::enum_<EnumeratorTypes>("EnumeratorType")
      .value("LinkNode", EnumeratorTypes::LinkNode)
      .value("PositionVariation", EnumeratorTypes::PositionVariation)
      .value("RepeatUnit", EnumeratorTypes::RepeatUnit)
      .value("StereoIsomer", EnumeratorTypes::StereoIsomer);

  python::class_<MolEnumerator::MolEnumeratorParams>(
      "MolEnumeratorParams", "Molecular enumerator parameters",
      python::init<>(python::args("self")))
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
           python::args("self", "typ"),
           "set the operator to be used for enumeration");
  python::def("Enumerate", &enumerateHelper1,
              (python::arg("mol"), python::arg("maxPerOperation") = 0,
               python::arg("enumerateStereo") = false),
              python::return_value_policy<python::manage_new_object>(),
              R"DOC(do an enumeration and return a MolBundle.
  If maxPerOperation is >0 that will be used as the maximum number of molecules which
    can be returned by any given operation.
  If enumerateStereo == true, the procedure will enumerate stereoisomers in addition
    to the other enumerator types.
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

  python::class_<MolEnumerator::StereoEnumerationOptions>(
      "StereoEnumerationOptions", "Stereoisomer enumeration parameters",
      python::init<>())
      .def("__init__", python::make_constructor(createParamsFromName))
      .def_readwrite(
          "tryEmbedding",
          &MolEnumerator::StereoEnumerationOptions::tryEmbedding,
          R"(if set the process attempts to generate a standard RDKit distance geometry
                conformation for the stereisomer. If this fails, we assume that the stereoisomer is
                non-physical and don't return it. NOTE that this is computationally expensive and is
                just a heuristic that could result in stereoisomers being lost.
              )")
      .def_readwrite(
          "onlyUnassigned",
          &MolEnumerator::StereoEnumerationOptions::onlyUnassigned,
          R"(if set (the default), stereocenters which have specified stereochemistry
                will not be perturbed unless they are part of a relative stereo
                group.)")
      .def_readwrite(
          "onlyStereoGroups",
          &MolEnumerator::StereoEnumerationOptions::onlyStereoGroups,
          R"(Only find stereoisomers that differ at the StereoGroups associated with the molecule.)")
      .def_readwrite(
          "unique", &MolEnumerator::StereoEnumerationOptions::unique,
          R"(If set, removes duplicate stereoisomers from the result.)")
      .def_readwrite("maxIsomers",
                     &MolEnumerator::StereoEnumerationOptions::maxIsomers,
                     R"(the maximum number of isomers to yield, if the
                number of possible isomers is greater than maxIsomers, a
                random subset will be yielded. If 0, all isomers are
                yielded. Since every additional stereo center doubles the
                number of results (and execution time) it's important to
                keep an eye on this.)");

  python::def(
      "GetStereoisomerCount",
      (unsigned int (*)(const ROMol &,
                        const MolEnumerator::StereoEnumerationOptions))
          MolEnumerator::get_stereoisomer_count,
      (python::arg("mol"),
       python::arg("options") = MolEnumerator::StereoEnumerationOptions()),
      R"DOC(get the total number of non-unique stereoisomers that can be be generated from mol.)DOC");

  python::def(
      "EnumerateStereoisomers", &enumerate_stereoisomers_helper,
      (python::arg("mol"),
       python::arg("options") = MolEnumerator::StereoEnumerationOptions()),
      python::return_value_policy<python::manage_new_object>(),
      R"DOC(enumerate stereoisomers and return a MolBundle)DOC");
}
