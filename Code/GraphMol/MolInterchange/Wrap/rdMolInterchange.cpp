//
//  Copyright (C) 2018 Greg Landrum
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

#include <GraphMol/MolInterchange/MolInterchange.h>

namespace python = boost::python;

namespace {
python::tuple JSONToMols(const std::string &jsonBlock,
                         python::object pyparams) {
  RDKit::MolInterchange::JSONParseParameters params;
  if (pyparams) {
    params =
        python::extract<RDKit::MolInterchange::JSONParseParameters>(pyparams);
  } else {
    params = RDKit::MolInterchange::defaultJSONParseParameters;
  }
  auto mols = RDKit::MolInterchange::JSONDataToMols(jsonBlock, params);
  python::list result;
  for (auto &mol : mols) {
    result.append(mol);
  }
  return python::tuple(result);
}

std::string MolsToJSON(const python::object &mols) {
  auto pymols = pythonObjectToVect<const RDKit::ROMol *>(mols);
  if (!pymols) {
    return "";
  }
  return RDKit::MolInterchange::MolsToJSONData(*pymols);
}
}  // namespace

BOOST_PYTHON_MODULE(rdMolInterchange) {
  python::scope().attr("__doc__") =
      "Module containing functions for interchange of molecules.\n"
      "Note that this should be considered beta and that the format\n"
      "  and API will very likely change in future releases.";

  python::class_<RDKit::MolInterchange::JSONParseParameters,
                 boost::noncopyable>("JSONParseParameters",
                                     "Parameters controlling the JSON parser")
      .def_readwrite(
          "setAromaticBonds",
          &RDKit::MolInterchange::JSONParseParameters::setAromaticBonds,
          "set bond types to aromatic for bonds flagged aromatic")
      .def_readwrite(
          "strictValenceCheck",
          &RDKit::MolInterchange::JSONParseParameters::strictValenceCheck,
          "be strict when checking atom valences")
      .def_readwrite(
          "parseConformers",
          &RDKit::MolInterchange::JSONParseParameters::parseConformers,
          "parse conformers in the JSON")
      .def_readwrite(
          "parseProperties",
          &RDKit::MolInterchange::JSONParseParameters::parseProperties,
          "parse molecular properties in the JSON")
      .def_readwrite("useHCounts",
                     &RDKit::MolInterchange::JSONParseParameters::useHCounts,
                     "use atomic H counts from the JSON. You may want to set "
                     "this to False when parsing queries.");

  std::string docString;
  docString =
      "Convert a single molecule to JSON\n\
\n\
    ARGUMENTS:\n\
      - mol: the molecule to work with\n\
    RETURNS:\n\
      a string\n";
  python::def("MolToJSON",
              (std::string(*)(
                  const RDKit::ROMol &))RDKit::MolInterchange::MolToJSONData,
              (python::arg("mol")), docString.c_str());
  docString =
      "Convert a set of molecules to JSON\n\
\n\
    ARGUMENTS:\n\
      - mols: the molecules to work with\n\
    RETURNS:\n\
      a string\n";
  python::def("MolsToJSON", MolsToJSON, (python::arg("mols")),
              docString.c_str());
  docString =
      "Convert JSON to a tuple of molecules\n\
\n\
    ARGUMENTS:\n\
      - jsonBlock: the molecule to work with\n\
      - params: (optional) JSONParseParameters controlling the JSON parsing\n\
    RETURNS:\n\
      a tuple of Mols\n";
  python::def(
      "JSONToMols", JSONToMols,
      (python::arg("jsonBlock"), python::arg("params") = python::object()),
      docString.c_str());
}
