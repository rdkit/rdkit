//
//  Copyright (C) 2018-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <RDBoost/boost_shared_ptr.h>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolInterchange/MolInterchange.h>

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(rdMolInterchange, m) {
  m.doc() = R"DOC(Module containing functions for interchange of molecules.
Note that this should be considered beta and that the format
and API will very likely change in future releases.)DOC";

  nb::class_<RDKit::MolInterchange::JSONParseParameters>(
      m, "JSONParseParameters", "Parameters controlling the JSON parser")
      .def(nb::init<>())
      .def_rw("setAromaticBonds",
              &RDKit::MolInterchange::JSONParseParameters::setAromaticBonds,
              "set bond types to aromatic for bonds flagged aromatic")
      .def_rw("strictValenceCheck",
              &RDKit::MolInterchange::JSONParseParameters::strictValenceCheck,
              "be strict when checking atom valences")
      .def_rw("parseConformers",
              &RDKit::MolInterchange::JSONParseParameters::parseConformers,
              "parse conformers in the JSON")
      .def_rw("parseProperties",
              &RDKit::MolInterchange::JSONParseParameters::parseProperties,
              "parse molecular properties in the JSON")
      .def_rw("useHCounts",
              &RDKit::MolInterchange::JSONParseParameters::useHCounts,
              R"DOC(use atomic H counts from the JSON. You may want to set
this to False when parsing queries.)DOC");

  nb::class_<RDKit::MolInterchange::JSONWriteParameters>(
      m, "JSONWriteParameters", "Parameters controlling the JSON writer")
      .def(nb::init<>())
      .def_rw("useRDKitExtensions",
              &RDKit::MolInterchange::JSONWriteParameters::useRDKitExtensions,
              "use RDKit extensions to the commonchem format");

  m.def(
      "MolToJSON",
      [](const RDKit::ROMol &mol, nb::object pyparams) {
        RDKit::MolInterchange::JSONWriteParameters params =
            RDKit::MolInterchange::defaultJSONWriteParameters;
        if (!pyparams.is_none()) {
          params =
              nb::cast<RDKit::MolInterchange::JSONWriteParameters>(pyparams);
        }
        return RDKit::MolInterchange::MolToJSONData(mol, params);
      },
      "mol"_a, "params"_a = nb::none(),
      R"DOC(Convert a single molecule to JSON

ARGUMENTS:
  - mol: the molecule to work with
RETURNS:
  a string)DOC");

  m.def(
      "MolsToJSON",
      [](nb::object mols_obj, nb::object pyparams) {
        std::vector<const RDKit::ROMol *> mols;
        for (nb::handle h : nb::iter(mols_obj)) {
          mols.push_back(nb::cast<const RDKit::ROMol *>(h));
        }
        RDKit::MolInterchange::JSONWriteParameters params =
            RDKit::MolInterchange::defaultJSONWriteParameters;
        if (!pyparams.is_none()) {
          params =
              nb::cast<RDKit::MolInterchange::JSONWriteParameters>(pyparams);
        }
        return RDKit::MolInterchange::MolsToJSONData(mols, params);
      },
      "mols"_a, "params"_a = nb::none(),
      R"DOC(Convert a set of molecules to JSON

ARGUMENTS:
  - mols: the molecules to work with
RETURNS:
  a string)DOC");

  m.def(
      "JSONToMols",
      [](const std::string &jsonBlock, nb::object pyparams) {
        RDKit::MolInterchange::JSONParseParameters params;
        if (!pyparams.is_none()) {
          params =
              nb::cast<RDKit::MolInterchange::JSONParseParameters>(pyparams);
        } else {
          params = RDKit::MolInterchange::defaultJSONParseParameters;
        }
        auto mols = RDKit::MolInterchange::JSONDataToMols(jsonBlock, params);
        nb::list result;
        for (auto &mol : mols) {
          result.append(mol);
        }
        return nb::tuple(result);
      },
      "jsonBlock"_a, "params"_a = nb::none(),
      R"DOC(Convert JSON to a tuple of molecules

ARGUMENTS:
  - jsonBlock: the molecule to work with
  - params: (optional) JSONParseParameters controlling the JSON parsing
RETURNS:
  a tuple of Mols)DOC");
}
