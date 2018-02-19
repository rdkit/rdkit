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

std::string MolsToJSON(const python::object &mols,
                       const std::string &datasetName) {
  auto pymols = pythonObjectToVect<const RDKit::ROMol *>(mols);
  return RDKit::MolInterchange::MolsToJSONData(*pymols, datasetName.c_str());
}
}

BOOST_PYTHON_MODULE(rdMolInterchange) {
  python::scope().attr("__doc__") =
      "Module containing functions for interchange of molecules";

  python::class_<RDKit::MolInterchange::JSONParseParameters,
                 boost::noncopyable>("JSONParseParameters",
                                     "Paramters controlling the JSON parser")
      .def_readwrite(
          "setAromaticBonds",
          &RDKit::MolInterchange::JSONParseParameters::setAromaticBonds,
          "set bond types to aromatic for bonds flagged aromatic")
      .def_readwrite(
          "strictValenceCheck",
          &RDKit::MolInterchange::JSONParseParameters::strictValenceCheck,
          "be strict when checking atom valences");

  std::string docString;
  docString =
      "Convert a single molecule to JSON\n\
\n\
    ARGUMENTS:\n\
      - mol: the molecule to work with\n\
      - datasetName: (optional) the name of the molecule set\n\
    RETURNS:\n\
      a string\n";
  python::def("MolToJSON", (std::string(*)(const RDKit::ROMol &, const char *))
                               RDKit::MolInterchange::MolToJSONData,
              (python::arg("mol"), python::arg("datasetName") = "rdkit mol"),
              docString.c_str());
  docString =
      "Convert a set of molecules to JSON\n\
\n\
    ARGUMENTS:\n\
      - mols: the molecules to work with\n\
      - datasetName: (optional) the name of the molecule set\n\
    RETURNS:\n\
      a string\n";
  python::def("MolsToJSON", MolsToJSON,
              (python::arg("mols"), python::arg("datasetName") = "rdkit mols"),
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
