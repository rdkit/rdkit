//
//  Copyright (c) 2025, Glysade Inc.
//  and other RDKit contributors
//
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutues for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include <ChemDraw/chemdraw.h>
#include <ChemDraw/chemdrawreaction.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/Wrap/props.hpp>
#include <RDBoost/python.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionRunner.h>
#include <GraphMol/ChemReactions/PreprocessRxn.h>
#include <GraphMol/ChemReactions/SanitizeRxn.h>
#include <GraphMol/MarvinParse/MarvinParser.h>
#include <GraphMol/Depictor/DepictUtils.h>
#include <GraphMol/FilterCatalog/FunctionalGroupHierarchy.h>

#include <RDBoost/Wrap.h>

#include <RDGeneral/Exceptions.h>
#include <GraphMol/SanitException.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/ChemReactions/ReactionFingerprints.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>

namespace python = boost::python;
using namespace RDKit;
namespace {
std::string pyObjectToString(python::object input) {
  python::extract<std::string> ex(input);
  if (ex.check()) {
    return ex();
  }
  std::wstring ws = python::extract<std::wstring>(input);
  return std::string(ws.begin(), ws.end());
}

python::object MolsFromChemDrawBlockHelper(const std::string &filename, bool sanitize,
                                 bool removeHs) {
  std::vector<std::unique_ptr<RWMol>> mols;
  try {
    mols = RDKit::v2::MolsFromChemDrawBlock(filename, {sanitize, removeHs});
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw python::error_already_set();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  python::list res;
  for (auto &mol : mols) {
    // take ownership of the data from the unique_ptr
    ROMOL_SPTR sptr(static_cast<ROMol *>(mol.release()));
    res.append(sptr);
  }
  return python::tuple(res);
}

python::tuple MolsFromChemDrawFileHelper(python::object cdxml, bool sanitize,
                            bool removeHs) {
  auto mols = RDKit::v2::MolsFromChemDrawFile(pyObjectToString(cdxml), {sanitize, removeHs});
  python::list res;
  for (auto &mol : mols) {
    // take ownership of the data from the unique_ptr
    ROMOL_SPTR sptr(static_cast<ROMol *>(mol.release()));
    res.append(sptr);
  }
  return python::tuple(res);
}

python::object ReactionsFromChemDrawFileHelper(const char *filename, bool sanitize,
                                      bool removeHs) {
  std::vector<std::unique_ptr<ChemicalReaction>> rxns;
  try {
    rxns = RDKit::v2::ChemDrawFileToChemicalReactions(filename, sanitize, removeHs);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw python::error_already_set();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  python::list res;
  for (auto &rxn : rxns) {
    // take ownership of the data from the unique_ptr
    res.append(std::shared_ptr<ChemicalReaction>(rxn.release()));
  }
  return python::tuple(res);
}

python::object ReactionsFromChemDrawBlockHelper(python::object imolBlock, bool sanitize,
                                       bool removeHs) {
  std::istringstream inStream(pyObjectToString(imolBlock));
  std::vector<std::unique_ptr<ChemicalReaction>> rxns;
  try {
    rxns = RDKit::v2::ChemDrawDataStreamToChemicalReactions(inStream, sanitize, removeHs);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  python::list res;
  for (auto &rxn : rxns) {
    // take ownership of the data from the unique_ptr
    res.append(std::shared_ptr<ChemicalReaction>(rxn.release()));
  }
  return python::tuple(res);
}
}

BOOST_PYTHON_MODULE(rdChemDraw) {
  python::scope().attr("__doc__") =
    "Module containing classes and functions for working with ChemDraw files.";

  // Molecule Interface
  std::string docString =
      R"DOC(Extract all molecules from a ChemDraw file.

     Note that the ChemDraw format is large and complex, the RDKit doesn't support
     full functionality, just the base ones required for molecule and
     reaction parsing.

     ARGUMENTS:

       - filename: the chemdraw filename (.cdx/.cdxml)

       - sanitize: if True, sanitize the molecules [default True]

       - removeHs: if True, convert explicit Hs into implicit Hs. [default True]

     RETURNS:
       a tuple of parsed Mol objects.)DOC";

  python::def("MolsFromChemDrawFile", MolsFromChemDrawFileHelper,
              (python::arg("filename"), python::arg("sanitize") = true,
               python::arg("removeHs") = true),
              docString.c_str());

  docString =
      R"DOC(Extract all molecules from a ChemDraw file.

     Note that the ChemDraw format is large and complex, the RDKit doesn't support
     full functionality, just the base ones required for molecule and
     reaction parsing.

     ARGUMENTS:

       - block: the CDX/CDXML block

       - sanitize: if True, sanitize the molecules [default True]

       - removeHs: if True, convert explicit Hs into implicit Hs. [default True]

     RETURNS:
       a tuple of parsed Mol objects.)DOC";

  python::def("MolsFromChemDrawBlock", MolsFromChemDrawBlockHelper,
              (python::arg("block"), python::arg("sanitize") = true,
               python::arg("removeHs") = true),
              docString.c_str());
  
  docString =
      R"DOC(Extract all reactions from a ChemDraw file.

     Note that the ChemDraw format is large and complex, the RDKit doesn't support
     full functionality, just the base ones required for molecule and
     reaction parsing.

     ARGUMENTS:

       - filename: the chemdraw filename (.cdx/.cdxml)

       - sanitize: if True, sanitize the molecules [default True]

       - removeHs: if True, convert explicit Hs into implicit Hs. [default True]

     RETURNS:
       an iterator of parsed ChemicalReaction objects.)DOC";

  // Reaction Interface
  python::def("ReactionsFromChemDrawFile", ReactionsFromChemDrawFileHelper,
              (python::arg("filename"), python::arg("sanitize") = false,
               python::arg("removeHs") = false),
              "construct a tuple of ChemicalReactions from a ChemDraw rxn file");

  python::def(
      "ReactionsFromChemDrawBlock", ReactionsFromChemDrawBlockHelper,
      (python::arg("rxnblock"), python::arg("sanitize") = false,
       python::arg("removeHs") = false),
      "construct a tuple of ChemicalReactions from a string in ChemDraw format");


  python::enum_<v2::CDXFormat>("CDXFormat")
    .value("CDX", v2::CDXFormat::CDX)
    .value("CDXML", v2::CDXFormat::CDXML);

    python::def(
		"MolToChemDrawBlock", v2::MolToChemDrawBlock,
		(python::arg("mol"), python::arg("format")=v2::CDXFormat::CDXML),
		"Convert a molecule into a chemdraw string using the specified format");
}
