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
#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap.h>

#include <GraphMol/Abbreviations/Abbreviations.h>

namespace python = boost::python;
using namespace RDKit;

namespace {

ROMol *condenseMolAbbreviationsHelper(const ROMol *mol,
                                      python::object pyabbrevs,
                                      double maxCoverage, bool sanitize) {
  RWMol *res = new RWMol(*mol);
  auto abbrevs =
      pythonObjectToVect<Abbreviations::AbbreviationDefinition>(pyabbrevs);
  if (abbrevs) {
    Abbreviations::condenseMolAbbreviations(*res, *abbrevs, maxCoverage,
                                            sanitize);
  }
  return rdcast<ROMol *>(res);
}

ROMol *condenseAbbreviationSGroupHelper(const ROMol *mol) {
  RWMol *res = new RWMol(*mol);
  Abbreviations::condenseAbbreviationSubstanceGroups(*res);
  return rdcast<ROMol *>(res);
}

ROMol *labelMolAbbreviationsHelper(const ROMol *mol, python::object pyabbrevs,
                                   double maxCoverage) {
  RWMol *res = new RWMol(*mol);
  auto abbrevs =
      pythonObjectToVect<Abbreviations::AbbreviationDefinition>(pyabbrevs);
  if (abbrevs) {
    Abbreviations::labelMolAbbreviations(*res, *abbrevs, maxCoverage);
  }
  return rdcast<ROMol *>(res);
}
}  // namespace

BOOST_PYTHON_MODULE(rdAbbreviations) {
  python::scope().attr("__doc__") =
      "Module containing functions for working with molecular abbreviations";
  // RegisterVectorConverter<Abbreviations::AbbreviationMatch>();
  RegisterVectorConverter<Abbreviations::AbbreviationDefinition>();

  python::class_<Abbreviations::AbbreviationDefinition>(
      "AbbreviationDefinition", "Abbreviation Definition", python::init<>())
      .def_readwrite("label", &Abbreviations::AbbreviationDefinition::label,
                     "the label")
      .def_readwrite(
          "displayLabel", &Abbreviations::AbbreviationDefinition::displayLabel,
          "the label in a drawing when the bond comes from the right")
      .def_readwrite("displayLabelW",
                     &Abbreviations::AbbreviationDefinition::displayLabelW,
                     "the label in a drawing when the bond comes from the west")
      .def_readwrite(
          "mol", &Abbreviations::AbbreviationDefinition::mol,
          "the query molecule (should have a dummy as the first atom)");

  python::def("GetDefaultAbbreviations",
              &Abbreviations::Utils::getDefaultAbbreviations,
              "returns a list of the default abbreviation definitions");
  python::def("GetDefaultLinkers", &Abbreviations::Utils::getDefaultLinkers,
              "returns a list of the default linker definitions");
  python::def("ParseAbbreviations", &Abbreviations::Utils::parseAbbreviations,
              (python::arg("text"), python::arg("removeExtraDummies") = false,
               python::arg("allowConnectionToDummies") = false),
              "returns a set of abbreviation definitions from a string");
  python::def("ParseLinkers", &Abbreviations::Utils::parseLinkers,
              (python::arg("text")),
              "returns a set of linker definitions from a string");
  python::def(
      "CondenseMolAbbreviations", &condenseMolAbbreviationsHelper,
      (python::arg("mol"), python::arg("abbrevs"),
       python::arg("maxCoverage") = 0.4, python::arg("sanitize") = true),
      python::return_value_policy<python::manage_new_object>(),
      "Finds and replaces abbreviations in a molecule. The result is not "
      "sanitized.");
  python::def("LabelMolAbbreviations", &labelMolAbbreviationsHelper,
              (python::arg("mol"), python::arg("abbrevs"),
               python::arg("maxCoverage") = 0.4),
              python::return_value_policy<python::manage_new_object>(),
              "Finds abbreviations and adds to them to a molecule as \"SUP\" "
              "SubstanceGroups");
  python::def(
      "CondenseAbbreviationSubstanceGroups", &condenseAbbreviationSGroupHelper,
      (python::arg("mol")),
      python::return_value_policy<python::manage_new_object>(),
      "Finds and replaces abbreviation (i.e. \"SUP\") substance groups in a "
      "molecule. The result is not sanitized.");
}
