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
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap_nb.h>

#include <GraphMol/Abbreviations/Abbreviations.h>

using namespace RDKit;

NB_MODULE(rdAbbreviations, m) {
  m.doc() =
      "Module containing functions for working with molecular abbreviations";
  nb::class_<Abbreviations::AbbreviationDefinition>(m, "AbbreviationDefinition",
                                                    "Abbreviation Definition")
      .def(nb::init<>())
      .def_rw("label", &Abbreviations::AbbreviationDefinition::label,
              "the label")
      .def_rw("displayLabel",
              &Abbreviations::AbbreviationDefinition::displayLabel,
              "the label in a drawing when the bond comes from the right")
      .def_rw("displayLabelW",
              &Abbreviations::AbbreviationDefinition::displayLabelW,
              "the label in a drawing when the bond comes from the west")
      .def_rw(
          "mol", &Abbreviations::AbbreviationDefinition::mol,
          "the query molecule (should have a dummy as the first atom if "
          "includesXBonds is true)")
      .def_rw("includesXBonds",
              &Abbreviations::AbbreviationDefinition::includesXBonds,
              "whether or not the abbreviation definition includes "
              "bonds to non-abbreviation atoms");

  m.def("GetDefaultAbbreviations",
        &Abbreviations::Utils::getDefaultAbbreviations,
        "returns a list of the default abbreviation definitions");
  m.def("GetDefaultLinkers", &Abbreviations::Utils::getDefaultLinkers,
        "returns a list of the default linker definitions");
  m.def(
      "ParseAbbreviations", &Abbreviations::Utils::parseAbbreviations,
      nb::arg("text"), nb::arg("removeExtraDummies") = false,
      nb::arg("allowConnectionToDummies") = false,
      R"DOC(Returns a set of abbreviation definitions from a string.
Format of the text data: A series of lines, each of which contains:

* label
* SMARTS
* displayLabel
* displayLabelW

Where 'label' is the label used for the abbreviation,
'SMARTS' is the SMARTS definition of the abbreviation,
'displayLabel' is used in drawings to render the abbreviations and
'displayLabelW' is the display label if a bond comes in from the right.
The 'displayLabel' and 'displayLabelW' fields are optional.
Use dummies in the SMARTS to indicate attachment points. The assumption
is that the first atom is a dummy (one will be added if this is not
true) and that the second atom is the surrogate for the rest of
the group.)DOC");
  m.def("ParseLinkers", &Abbreviations::Utils::parseLinkers, nb::arg("text"),
        "Returns a set of linker definitions from a string. Equivalent to "
        "calling ParseAbbreviations(text, True True).");
  m.def(
      "CondenseMolAbbreviations",
      [](const ROMol *mol,
         std::vector<Abbreviations::AbbreviationDefinition> abbrevs,
         double maxCoverage, bool sanitize) {
        RWMol *res = new RWMol(*mol);
        Abbreviations::condenseMolAbbreviations(*res, abbrevs, maxCoverage,
                                                sanitize);
        return rdcast<ROMol *>(res);
      },
      nb::arg("mol"), nb::arg("abbrevs"), nb::arg("maxCoverage") = 0.4,
      nb::arg("sanitize") = true, nb::rv_policy::take_ownership,
      "Finds and replaces abbreviations in a molecule. The result is not "
      "sanitized.");
  m.def(
      "LabelMolAbbreviations",
      [](const ROMol *mol,
         std::vector<Abbreviations::AbbreviationDefinition> abbrevs,
         double maxCoverage) {
        RWMol *res = new RWMol(*mol);
        Abbreviations::labelMolAbbreviations(*res, abbrevs, maxCoverage);
        return rdcast<ROMol *>(res);
      },
      nb::arg("mol"), nb::arg("abbrevs"), nb::arg("maxCoverage") = 0.4,
      nb::rv_policy::take_ownership,
      "Finds abbreviations and adds to them to a molecule as \"SUP\" "
      "SubstanceGroups.");
  m.def(
      "CondenseAbbreviationSubstanceGroups",
      [](const ROMol *mol) {
        RWMol *res = new RWMol(*mol);
        Abbreviations::condenseAbbreviationSubstanceGroups(*res);
        return rdcast<ROMol *>(res);
      },
      nb::arg("mol"), nb::rv_policy::take_ownership,
      "Finds and replaces abbreviation (i.e. \"SUP\") substance groups in a "
      "molecule. The result is not sanitized.");
}
